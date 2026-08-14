#!/usr/bin/env python3
"""Rebuild per-subject motion summary from DPARSF RealignParameter outputs.

Changes vs prior version:
  - Power-style neighbor-expanded scrubbing (1-back, 2-forward; Power et al. 2014)
  - Contiguous-run accounting (drops residual segments shorter than MIN_SEGMENT_LEN)
  - Minutes of usable data remaining (requires TR_SECONDS)
  - Late-vs-early FD metrics: mean of first/last thirds, ratio, linear slope
  - Revised exclusion rule:  mean FD > 0.5 mm  OR  < 5 min clean data remaining
  - FD threshold raised from 0.2 to 0.3 mm
  - Plots now shade scrubbed frames
  - Legacy exclusion columns kept side-by-side for sanity comparison
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ----------------------------- CONFIG -----------------------------
REALIGN_DIR = Path("/home/martin/data_imaging/CovidProject/FullStudy/fMRI/DPARSF/RealignParameter")
OUT_CSV     = Path("/home/martin/data_imaging/CovidProject/FullStudy/fMRI/DPARSF/motion_summary.csv")
OUT_EXCLUDED_CSV = Path("/home/martin/data_imaging/CovidProject/FullStudy/fMRI/DPARSF/excluded_subjects_by_criteria.csv")
PLOT_DIR    = Path("/home/martin/data_imaging/CovidProject/FullStudy/fMRI/DPARSF/FD_power_plots")

TR_SECONDS       = 3.0    # *** SET TO YOUR SEQUENCE TR (s) ***  placeholder value
FD_THRESH        = 0.3    # per-volume censoring threshold (mm) — revised up from 0.2
SCRUB_BACK       = 1      # frames censored BEFORE each high-FD frame (Power 2014)
SCRUB_FORWARD    = 2      # frames censored AFTER  each high-FD frame
MIN_SEGMENT_LEN  = 5      # min contiguous clean run to count toward usable data
MIN_CLEAN_MIN    = 5.0    # min minutes of usable data after scrubbing for inclusion
MEAN_FD_LIMIT    = 0.5    # subject-level mean FD ceiling
LATE_RATIO_FLAG  = 2.0    # flag late-deterioration if last-third mean > this × first-third
# ------------------------------------------------------------------


def load_fd(subject_dir: Path, metric: str) -> np.ndarray | None:
    matches = list(subject_dir.glob(f"FD_{metric}_*.txt"))
    if not matches:
        return None
    return np.loadtxt(matches[0])


def power_scrub_mask(fd: np.ndarray, threshold: float,
                     back: int = 1, forward: int = 2) -> np.ndarray:
    """Boolean mask (True = KEEP) after Power-style neighbor-expanded scrubbing.

    For every frame i where FD[i] > threshold, also censors frames
    [i-back, ..., i+forward].
    """
    n = len(fd)
    bad = fd > threshold
    keep = ~bad.copy()
    for i in np.where(bad)[0]:
        lo = max(0, i - back)
        hi = min(n, i + forward + 1)
        keep[lo:hi] = False
    return keep


def run_lengths(mask: np.ndarray) -> list[int]:
    """Lengths of contiguous True runs in a boolean array."""
    if mask.size == 0:
        return []
    padded = np.concatenate(([False], mask, [False]))
    diffs = np.diff(padded.astype(int))
    starts = np.where(diffs == 1)[0]
    ends   = np.where(diffs == -1)[0]
    return (ends - starts).tolist()


def late_early_metrics(fd: np.ndarray) -> dict:
    """Quantify progressive deterioration across the run."""
    n = len(fd)
    third = max(1, n // 3)
    first = fd[:third]
    last  = fd[-third:]
    slope, _ = np.polyfit(np.arange(n), fd, 1)
    ratio = float(np.mean(last) / np.mean(first)) if np.mean(first) > 0 else np.nan
    return {
        "mean_FD_first_third": float(np.mean(first)),
        "mean_FD_last_third":  float(np.mean(last)),
        "FD_late_early_ratio": ratio,
        "FD_slope_per_frame":  float(slope),
    }


def summarize_subject(subject_dir: Path):
    """Return (row_dict, fd_power_full, keep_mask) or (None, None, None)."""
    sid = subject_dir.name
    fd_power = load_fd(subject_dir, "Power")
    if fd_power is None:
        return None, None, None

    fd_jenk    = load_fd(subject_dir, "Jenkinson")
    fd_vandijk = load_fd(subject_dir, "VanDijk")

    # Drop leading 0 (volume 1 has no predecessor)
    fd_p = fd_power[1:]
    fd_j = fd_jenk[1:]    if fd_jenk    is not None else None
    fd_v = fd_vandijk[1:] if fd_vandijk is not None else None
    n = len(fd_p)

    over_thresh = fd_p > FD_THRESH

    keep_mask = power_scrub_mask(fd_p, FD_THRESH, SCRUB_BACK, SCRUB_FORWARD)
    runs = run_lengths(keep_mask)
    n_clean        = int(keep_mask.sum())
    n_clean_valid  = int(sum(r for r in runs if r >= MIN_SEGMENT_LEN))
    longest_run    = max(runs) if runs else 0

    le = late_early_metrics(fd_p)

    row = {
        "SubjectID":            sid,
        "nVols":                len(fd_power),
        "mean_FD_Power":        float(np.mean(fd_p)),
        "max_FD_Power":         float(np.max(fd_p)),
        "median_FD_Power":      float(np.median(fd_p)),
        "mean_FD_Jenk":         float(np.mean(fd_j)) if fd_j is not None else np.nan,
        "mean_FD_VanDijk":      float(np.mean(fd_v)) if fd_v is not None else np.nan,
        # Per-volume threshold counts (raw, before neighbor expansion)
        "n_FD_gt_thresh":       int(over_thresh.sum()),
        "pct_FD_gt_thresh":     float(over_thresh.mean() * 100),
        # Scrubbing results (with neighbor expansion)
        "n_vols_scrubbed":      n - n_clean,
        "n_vols_clean":         n_clean,
        "n_vols_clean_in_valid_segments": n_clean_valid,
        "longest_clean_run":    longest_run,
        "minutes_clean":        n_clean * TR_SECONDS / 60.0,
        "minutes_clean_valid":  n_clean_valid * TR_SECONDS / 60.0,
        # Late vs early
        **le,
    }
    return row, fd_power, keep_mask


def save_fd_power_plot(subject_id: str, fd_power: np.ndarray,
                       keep_mask: np.ndarray) -> None:
    """Per-frame FD trace with scrubbed frames shaded."""
    fd_p = fd_power[1:]
    n = len(fd_p)
    frames = np.arange(2, n + 2)   # 1-indexed; volume 1 was dropped

    fig, ax = plt.subplots(figsize=(10, 3.5))
    ax.plot(frames, fd_p, color="tab:blue", linewidth=1.2)
    ax.axhline(FD_THRESH, color="tab:red", linestyle="--", linewidth=1.0,
               label=f"FD = {FD_THRESH} mm")

    scrubbed = ~keep_mask
    n_scrub = int(scrubbed.sum())
    if n_scrub > 0:
        idx = np.where(scrubbed)[0]
        groups = np.split(idx, np.where(np.diff(idx) != 1)[0] + 1)
        for g in groups:
            if len(g) == 0:
                continue
            ax.axvspan(frames[g[0]] - 0.5, frames[g[-1]] + 0.5,
                       color="tab:red", alpha=0.15, linewidth=0)

    ax.set_xlabel("Frame")
    ax.set_ylabel("FD Power (mm)")
    ax.set_title(f"{subject_id}  —  scrubbed {n_scrub}/{n} frames "
                 f"(thr {FD_THRESH} mm, ±{SCRUB_BACK}/{SCRUB_FORWARD})")
    ax.legend(loc="upper right", fontsize=8)
    plt.tight_layout()
    plt.savefig(PLOT_DIR / f"{subject_id}_FD_power.png", dpi=150)
    plt.close()


def main():
    subject_dirs = sorted(d for d in REALIGN_DIR.iterdir() if d.is_dir())
    print(f"Found {len(subject_dirs)} subject directories")
    print(f"TR = {TR_SECONDS} s  |  FD_thresh = {FD_THRESH} mm  "
          f"|  Scrub ±{SCRUB_BACK}/{SCRUB_FORWARD}  "
          f"|  Min clean = {MIN_CLEAN_MIN} min in segments ≥ {MIN_SEGMENT_LEN} frames")
    PLOT_DIR.mkdir(parents=True, exist_ok=True)

    rows, missing = [], []
    for sdir in subject_dirs:
        row, fd_full, keep_mask = summarize_subject(sdir)
        if row is None:
            missing.append(sdir.name)
            continue
        save_fd_power_plot(sdir.name, fd_full, keep_mask)
        rows.append(row)

    df = pd.DataFrame(rows).sort_values("SubjectID").reset_index(drop=True)

    # ---- Revised exclusion logic ----
    df["exclude_high_mean_fd"]      = df["mean_FD_Power"]       > MEAN_FD_LIMIT
    df["exclude_insufficient_data"] = df["minutes_clean_valid"] < MIN_CLEAN_MIN
    df["exclude_any"] = df["exclude_high_mean_fd"] | df["exclude_insufficient_data"]
    # Late-deterioration: flag only, not for exclusion
    df["flag_late_deterioration"] = (
        (df["FD_late_early_ratio"] > LATE_RATIO_FLAG) &
        (df["mean_FD_last_third"]  > FD_THRESH)
    )
    # Legacy flags for direct comparison with the old criteria
    df["exclude_legacy_strict"]  = (df["mean_FD_Power"] > 0.25) | (df["pct_FD_gt_thresh"] > 20)
    df["exclude_legacy_liberal"] = (df["mean_FD_Power"] > 0.5)  | (df["pct_FD_gt_thresh"] > 50)

    df.to_csv(OUT_CSV, index=False)
    print(f"\nWrote {OUT_CSV}  ({len(df)} subjects)")

    # Excluded-subject lists by requested criteria.
    # Mapping used here:
    #   strict  -> legacy strict thresholding
    #   legacy  -> legacy liberal thresholding
    #   revised -> current revised exclusion rule
    excluded_rows = []
    criteria_map = {
        "strict": "exclude_legacy_strict",
        "legacy": "exclude_legacy_liberal",
        "revised": "exclude_any",
    }
    for criterion_name, col in criteria_map.items():
        excluded_ids = df.loc[df[col], "SubjectID"].astype(str).tolist()
        for sid in excluded_ids:
            excluded_rows.append({"criterion": criterion_name, "SubjectID": sid})

    excluded_df = pd.DataFrame(excluded_rows).sort_values(["criterion", "SubjectID"]).reset_index(drop=True)
    excluded_df.to_csv(OUT_EXCLUDED_CSV, index=False)
    print(f"Wrote {OUT_EXCLUDED_CSV}  ({len(excluded_df)} rows)")

    if missing:
        print(f"\nWARNING: no FD_Power file for {len(missing)} subjects:")
        for m in missing:
            print(f"  - {m}")

    # ---- Console summary ----
    print("\nDistribution summary:")
    print(df[["mean_FD_Power", "max_FD_Power", "pct_FD_gt_thresh",
              "minutes_clean_valid", "FD_late_early_ratio"]].describe().round(3))

    print("\nExclusion comparison:")
    n = len(df)
    print(f"  Legacy strict  (meanFD>0.25 or >20% at FD>0.2): {df['exclude_legacy_strict'].sum():3d}/{n}")
    print(f"  Legacy liberal (meanFD>0.5  or >50% at FD>0.2): {df['exclude_legacy_liberal'].sum():3d}/{n}")
    print(f"  REVISED        (meanFD>{MEAN_FD_LIMIT} or <{MIN_CLEAN_MIN} min clean): "
          f"{df['exclude_any'].sum():3d}/{n}")
    print(f"      ↳ excluded by high mean FD:    {df['exclude_high_mean_fd'].sum():3d}")
    print(f"      ↳ excluded by insufficient data: {df['exclude_insufficient_data'].sum():3d}")
    print(f"  Late-deterioration FLAG (kept):  {df['flag_late_deterioration'].sum():3d}")


if __name__ == "__main__":
    main()