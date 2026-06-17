"""
match_preprocessed_batches_to_adnimerge.py

Read all preprocessed batch CSVs from a folder, extract ADNI subject/site IDs
from Subject fields, and match them to ADNIMERGE rows using SUBJID.

Input batch CSVs must contain:
- Subject field: "SubjectID" or "Subject" (format like 003_S_4644)
- Image field: "ImageID" or "Image Data ID"

ADNIMERGE file must contain:
- SUBJID
- SITEID

Output:
A single merged CSV that keeps one row per input row x matching ADNIMERGE row,
with identifier columns from input plus all ADNIMERGE columns.

Note:
Input site IDs from XXX_S_YYYY and ADNIMERGE SITEID are both preserved, but
matching is done by SUBJID to avoid site coding mismatches between datasets.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Optional

import pandas as pd

DEFAULT_BATCH_FOLDER = Path(
    "/home/martin/data_imaging/ADNIdata/RevisionPaperfMRI_2026_04/AllCsvs/Batches"
)
DEFAULT_ADNIMERGE = Path(
    "/home/martin/data_imaging/ADNIdata/StudyInfo_2026/AdniMergeDefault_2026_06_11.csv"
)
DEFAULT_SLICE_TIMING_CSV = Path(
    "/home/martin/data_imaging/ADNIdata/RevisionPaperfMRI_2026_04/AllCsvs/fmri_slice_timing_with_imageid.csv"
)
DEFAULT_OUTPUT = "/home/martin/data_imaging/ADNIdata/RevisionPaperfMRI_2026_04/AllCsvs/preprocessed_batches_with_adnimerge.csv"

SUBJECT_COL_CANDIDATES = ["SubjectID", "Subject"]
IMAGE_COL_CANDIDATES = ["ImageID", "Image Data ID"]

# Accept site as 1-3 digits to be robust; normalize for matching later.
SUBJECT_PATTERN = re.compile(r"^\s*(\d{1,3})_S_(\d+)\s*$", re.IGNORECASE)


def normalize_numeric_id(value: object) -> str:
    """Normalize IDs for matching (remove spaces, .0, leading zeros)."""
    if pd.isna(value):
        return ""

    text = str(value).strip()
    if text == "":
        return ""

    # Keep digits only to tolerate formats like "4644.0" or "RID=4644".
    digits = "".join(ch for ch in text if ch.isdigit())
    if digits == "":
        return ""

    return str(int(digits))


def normalize_site_id(value: object) -> str:
    """Normalize SITEID for matching without relying on zero padding."""
    return normalize_numeric_id(value)


def normalize_image_id(value: object) -> str:
    """Normalize ImageID values for robust matching across CSVs."""
    if pd.isna(value):
        return ""

    text = str(value).strip()
    if text == "":
        return ""

    # Convert pure numeric-like strings (including float notation) to integer-like IDs.
    try:
        return str(int(float(text)))
    except ValueError:
        pass

    match = re.match(r"^[iI](\d+)$", text)
    if match:
        return match.group(1)

    # Fallback: keep only digits so IDs like "I12345" or "ID=12345" match.
    digits = "".join(ch for ch in text if ch.isdigit())
    return str(int(digits)) if digits else text


def find_first_existing_column(df: pd.DataFrame, candidates: list[str]) -> Optional[str]:
    for col in candidates:
        if col in df.columns:
            return col
    return None


def extract_subject_parts(subject_value: object) -> tuple[str, str, bool]:
    """
    Return (site_raw, subjid_raw, valid_format) from values like "003_S_4644".
    """
    if pd.isna(subject_value):
        return "", "", False

    text = str(subject_value).strip()
    match = SUBJECT_PATTERN.match(text)
    if not match:
        return "", "", False

    site_raw, subjid_raw = match.group(1), match.group(2)
    return site_raw, subjid_raw, True


def read_batch_csvs(batch_folder: Path) -> pd.DataFrame:
    csv_files = sorted(batch_folder.glob("*.csv"))
    csv_files = [
        f
        for f in csv_files
        if not f.name.startswith("preprocessed_batches_with_adnimerge")
    ]
    if not csv_files:
        raise FileNotFoundError(f"No CSV files found in: {batch_folder}")

    extracted_rows: list[pd.DataFrame] = []

    for csv_file in csv_files:
        batch_df = pd.read_csv(csv_file)
        batch_df.columns = batch_df.columns.str.strip()

        subject_col = find_first_existing_column(batch_df, SUBJECT_COL_CANDIDATES)
        image_col = find_first_existing_column(batch_df, IMAGE_COL_CANDIDATES)

        if subject_col is None:
            print(
                f"[SKIP] {csv_file.name}: missing subject/image columns. "
                f"Found subject={subject_col}, image={image_col}"
            )
            continue

        temp = pd.DataFrame(index=batch_df.index)
        temp["SourceCSV"] = csv_file.name
        temp["SubjectRaw"] = batch_df[subject_col]
        if image_col is None:
            print(f"[WARN] {csv_file.name}: image column not found, writing empty ImageID")
            temp["ImageID"] = ""
        else:
            temp["ImageID"] = batch_df[image_col].map(normalize_image_id)

        parsed = temp["SubjectRaw"].map(extract_subject_parts)
        temp["InputSiteRaw"] = parsed.map(lambda x: x[0])
        temp["InputSubjRaw"] = parsed.map(lambda x: x[1])
        temp["SubjectFormatValid"] = parsed.map(lambda x: x[2])

        temp["InputSiteID"] = temp["InputSiteRaw"].map(normalize_site_id)
        temp["SUBJID_MATCH"] = temp["InputSubjRaw"].map(normalize_numeric_id)

        extracted_rows.append(temp)

    if not extracted_rows:
        raise ValueError(
            "No usable CSVs found. Ensure files contain SubjectID/Subject and ImageID/Image Data ID."
        )

    return pd.concat(extracted_rows, ignore_index=True)


def read_and_prepare_adnimerge(adnimerge_csv: Path) -> pd.DataFrame:
    adni = pd.read_csv(adnimerge_csv, low_memory=False)
    adni.columns = adni.columns.str.strip()

    required = ["SUBJID", "SITEID"]
    missing = [col for col in required if col not in adni.columns]
    if missing:
        raise ValueError(f"ADNIMERGE missing required columns: {missing}")

    adni["SUBJID_MATCH"] = adni["SUBJID"].map(normalize_numeric_id)
    adni["SITEID_MATCH"] = adni["SITEID"].map(normalize_site_id)
    return adni


def read_and_prepare_slice_timing(slice_timing_csv: Path) -> pd.DataFrame:
    slice_df = pd.read_csv(slice_timing_csv, low_memory=False)
    slice_df.columns = slice_df.columns.str.strip()

    if "IMAGEUID_fmri" not in slice_df.columns:
        raise ValueError("Slice timing CSV missing required column: IMAGEUID_fmri")

    slice_df["IMAGEID_MATCH"] = slice_df["IMAGEUID_fmri"].map(normalize_image_id)
    return slice_df


def build_output(
    batch_folder: Path,
    adnimerge_csv: Path,
    slice_timing_csv: Path,
    output_csv: Path,
) -> None:
    batch_ids = read_batch_csvs(batch_folder)
    adni = read_and_prepare_adnimerge(adnimerge_csv)
    slice_df = read_and_prepare_slice_timing(slice_timing_csv)

    batch_ids["IMAGEID_MATCH"] = batch_ids["ImageID"].map(normalize_image_id)

    merged = batch_ids.merge(
        adni,
        on=["SUBJID_MATCH"],
        how="left",
        suffixes=("", "_ADNIMERGE"),
    )

    merged = merged.merge(
        slice_df,
        on=["IMAGEID_MATCH"],
        how="left",
        suffixes=("", "_SLICE"),
    )

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(output_csv, index=False)

    total = len(batch_ids)
    valid_subject = int(batch_ids["SubjectFormatValid"].sum())
    adni_keys = set(adni["SUBJID_MATCH"]) - {""}
    slice_keys = set(slice_df["IMAGEID_MATCH"]) - {""}
    matched = int(batch_ids["SUBJID_MATCH"].isin(adni_keys).sum())
    unmatched = total - matched
    matched_slice = int(batch_ids["IMAGEID_MATCH"].isin(slice_keys).sum())

    print(f"Input folder: {batch_folder}")
    print(f"ADNIMERGE:   {adnimerge_csv}")
    print(f"Slice timing CSV: {slice_timing_csv}")
    print(f"Output CSV:  {output_csv}")
    print("-")
    print(f"Rows from input CSVs: {total}")
    print(f"Rows with valid subject format XXX_S_YYYY: {valid_subject}")
    print(f"Rows matched to ADNIMERGE by subjid: {matched}")
    print(f"Rows unmatched: {unmatched}")
    print(f"Rows matched to slice timing by image id: {matched_slice}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Merge preprocessed batch CSVs with ADNIMERGE using subject IDs in format "
            "XXX_S_YYYY and ADNIMERGE fields SITEID/SUBJID."
        )
    )
    parser.add_argument(
        "--batch-folder",
        type=Path,
        default=DEFAULT_BATCH_FOLDER,
        help="Folder containing batch CSV files.",
    )
    parser.add_argument(
        "--adnimerge-csv",
        type=Path,
        default=DEFAULT_ADNIMERGE,
        help="ADNIMERGE CSV path.",
    )
    parser.add_argument(
        "--slice-timing-csv",
        type=Path,
        default=DEFAULT_SLICE_TIMING_CSV,
        help="Slice timing CSV path with IMAGEUID_fmri.",
    )
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=DEFAULT_OUTPUT,
        help="Output merged CSV path.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    build_output(
        args.batch_folder,
        args.adnimerge_csv,
        args.slice_timing_csv,
        args.output_csv,
    )


if __name__ == "__main__":
    main()
