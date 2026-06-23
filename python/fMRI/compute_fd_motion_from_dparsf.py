#!/usr/bin/env python3
"""Rebuild per-subject motion summary from DPARSF RealignParameter outputs."""

from pathlib import Path
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

REALIGN_DIR = Path("/home/martin/data_imaging/CovidProject/FullStudy/fMRI/DPARSF/RealignParameter")
OUT_CSV     = Path("/home/martin/data_imaging/CovidProject/FullStudy/fMRI/DPARSF/motion_summary.csv")
PLOT_DIR    = Path("/home/martin/data_imaging/CovidProject/FullStudy/fMRI/DPARSF/FD_power_plots")
FD_THRESH   = 0.2   # Power threshold (mm)


def load_fd(subject_dir: Path, metric: str) -> np.ndarray | None:
    """Load FD vector for a subject; returns None if file missing."""
    matches = list(subject_dir.glob(f"FD_{metric}_*.txt"))
    if not matches:
        return None
    return np.loadtxt(matches[0])


def summarize_subject(subject_dir: Path) -> dict | None:
    sid = subject_dir.name
    fd_power = load_fd(subject_dir, "Power")
    if fd_power is None:
        return None

    fd_jenk   = load_fd(subject_dir, "Jenkinson")
    fd_vandijk = load_fd(subject_dir, "VanDijk")

    # Drop leading 0 (volume 1 has no predecessor)
    fd_p = fd_power[1:]
    fd_j = fd_jenk[1:] if fd_jenk is not None else None
    fd_v = fd_vandijk[1:] if fd_vandijk is not None else None

    over_thresh = fd_p > FD_THRESH

    return {
        "SubjectID":        sid,
        "nVols":            len(fd_power),
        "mean_FD_Power":    float(np.mean(fd_p)),
        "max_FD_Power":     float(np.max(fd_p)),
        "median_FD_Power":  float(np.median(fd_p)),
        "mean_FD_Jenk":     float(np.mean(fd_j)) if fd_j is not None else np.nan,
        "mean_FD_VanDijk":  float(np.mean(fd_v)) if fd_v is not None else np.nan,
        "n_FD_gt_thresh":   int(over_thresh.sum()),
        "pct_FD_gt_thresh": float(over_thresh.mean() * 100),
        "n_vols_surviving": int((~over_thresh).sum()),
    }


def save_fd_power_plot(subject_id: str, fd_power: np.ndarray) -> None:
    """Save per-frame FD_power trace for one subject."""
    frames = np.arange(1, len(fd_power) + 1)

    plt.figure(figsize=(10, 3.5))
    plt.plot(frames, fd_power, color="tab:blue", linewidth=1.2)
    plt.axhline(FD_THRESH, color="tab:red", linestyle="--", linewidth=1.0)
    plt.xlabel("Frame")
    plt.ylabel("FD Power (mm)")
    plt.title(f"{subject_id} FD Power per frame")
    plt.tight_layout()
    plt.savefig(PLOT_DIR / f"{subject_id}_FD_power.png", dpi=150)
    plt.close()


def main():
    subject_dirs = sorted(d for d in REALIGN_DIR.iterdir() if d.is_dir())
    print(f"Found {len(subject_dirs)} subject directories")
    PLOT_DIR.mkdir(parents=True, exist_ok=True)

    rows = []
    missing = []
    for sdir in subject_dirs:
        fd_power = load_fd(sdir, "Power")
        if fd_power is not None:
            save_fd_power_plot(sdir.name, fd_power)

        row = summarize_subject(sdir)
        if row is None:
            missing.append(sdir.name)
        else:
            rows.append(row)

    df = pd.DataFrame(rows).sort_values("SubjectID").reset_index(drop=True)

    # Flag-style columns for quick filtering
    df["exclude_strict"]  = (df["mean_FD_Power"] > 0.25) | (df["pct_FD_gt_thresh"] > 20)
    df["exclude_liberal"] = (df["mean_FD_Power"] > 0.5)  | (df["pct_FD_gt_thresh"] > 50)

    df.to_csv(OUT_CSV, index=False)
    print(f"Wrote {OUT_CSV}  ({len(df)} subjects)")

    if missing:
        print(f"\nWARNING: no FD_Power file for {len(missing)} subjects:")
        for m in missing:
            print(f"  - {m}")

    # Quick console summary
    print("\nMotion summary:")
    print(df[["mean_FD_Power", "max_FD_Power", "pct_FD_gt_thresh"]].describe().round(3))
    print(f"\nFlagged strict  (meanFD>0.25 or >20% censored): {df['exclude_strict'].sum()}")
    print(f"Flagged liberal (meanFD>0.5  or >50% censored): {df['exclude_liberal'].sum()}")


if __name__ == "__main__":
    main()