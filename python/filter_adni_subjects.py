"""
filter_adni_subjects.py

Reads the new ADNI subject search CSV and filters out subjects that are already
present in the previously processed batch CSV (same Subject ID AND same Visit).

If a subject ID is already in the batch CSV but for a *different* visit, that
row is kept and appended at the end with a flag column `already_in_batch=True`.

Also extracts the filtered subjects from a zip file containing ADNI DICOM data,
separating them into:
    - New/ADNI/<SubjectID>/...          (truly new subjects)
    - AdditionalTimePoints/ADNI/<SubjectID>/...  (same subject, different visit)

Inputs:
    new_csv   : /home/martin/data_imaging/ADNIdata/AD_fMRI_ida_search_4_28_2026.csv
    batch_csv : /home/martin/data_imaging/ADNIdata/fMRI_first_batch_to_process_2023_04_27/fMRI_first_batch_to_process.csv
    zip_file  : /home/martin/data_imaging/ADNIdata/AD fMRI.zip

Output:
    filtered CSV saved next to the new_csv file.
    extracted data in /home/martin/data_imaging/ADNIdata/RevisionPaperfMRI_AD/
"""

import zipfile
import pandas as pd
from pathlib import Path

# ── Paths ─────────────────────────────────────────────────────────────────────
NEW_CSV   = Path("/home/martin/data_imaging/ADNIdata/AD_fMRI_ida_search_4_28_2026.csv")
BATCH_CSV = Path("/home/martin/data_imaging/ADNIdata/fMRI_ADNI2_ADNI3_initial_visit/SubjectsDataAndTests.csv")
OUT_CSV   = NEW_CSV.parent / "AD_fMRI_filtered_4_28_2026.csv"
ZIP_FILE  = Path("/home/martin/data_imaging/ADNIdata/AD fMRI.zip")
EXTRACT_ROOT = Path("/home/martin/data_imaging/ADNIdata/RevisionPaperfMRI_AD")

# ── Load ──────────────────────────────────────────────────────────────────────
new_df   = pd.read_csv(NEW_CSV)
batch_df = pd.read_csv(BATCH_CSV)

# Normalise column names (strip whitespace, unify to snake_case keys we need)
new_df.columns   = new_df.columns.str.strip()
batch_df.columns = batch_df.columns.str.strip()

# Identify the Subject ID and Visit columns in each file
# new_csv   → "Subject ID",  "Visit"
# batch_csv → "SubjectID",   "Visit"
NEW_SUBJ_COL  = "Subject ID"
NEW_VISIT_COL = "Visit"
BAT_SUBJ_COL  = "SubjectID"
BAT_VISIT_COL = "Visit"

# Build sets from the batch CSV for fast lookup
batch_exact   = set(zip(batch_df[BAT_SUBJ_COL], batch_df[BAT_VISIT_COL]))   # (subj, visit)
batch_subjects = set(batch_df[BAT_SUBJ_COL])                                 # subject only

# ── Classify each row in the new CSV ─────────────────────────────────────────
def classify(row):
    key = (row[NEW_SUBJ_COL], row[NEW_VISIT_COL])
    if key in batch_exact:
        return "exact_match"          # same subject AND same visit → exclude
    elif row[NEW_SUBJ_COL] in batch_subjects:
        return "different_visit"      # same subject, different visit → keep + flag
    else:
        return "new"                  # completely new subject → keep

new_df["_status"] = new_df.apply(classify, axis=1)

# ── Build output ──────────────────────────────────────────────────────────────
truly_new   = new_df[new_df["_status"] == "new"].copy()
diff_visit  = new_df[new_df["_status"] == "different_visit"].copy()
excluded    = new_df[new_df["_status"] == "exact_match"].copy()

# Add flag column only to the "different visit" rows, False for new rows
truly_new["already_in_batch"]  = False
diff_visit["already_in_batch"] = True

filtered_df = pd.concat([truly_new, diff_visit], ignore_index=True)
filtered_df = filtered_df.drop(columns=["_status"])

# ── Save ──────────────────────────────────────────────────────────────────────
filtered_df.to_csv(OUT_CSV, index=False)

# ── Summary ───────────────────────────────────────────────────────────────────
print(f"New CSV total rows         : {len(new_df)}")
print(f"Batch CSV total rows       : {len(batch_df)}")
print(f"Excluded (exact matches)   : {len(excluded)}")
print(f"Kept – truly new subjects  : {len(truly_new)}")
print(f"Kept – same subj, diff visit (flagged): {len(diff_visit)}")
print(f"Output rows                : {len(filtered_df)}")
print(f"\nOutput saved to: {OUT_CSV}")

if not diff_visit.empty:
    print("\nSubjects present in batch with a different visit:")
    print(diff_visit[[NEW_SUBJ_COL, NEW_VISIT_COL]].to_string(index=False))

# ── Extract from zip ──────────────────────────────────────────────────────────
print(f"\n{'─'*60}")
print(f"Extracting from: {ZIP_FILE}")

new_subjects       = set(truly_new[NEW_SUBJ_COL])
flagged_subjects   = set(diff_visit[NEW_SUBJ_COL])

dest_new        = EXTRACT_ROOT / "New"
dest_additional = EXTRACT_ROOT / "AdditionalTimePoints"
dest_new.mkdir(parents=True, exist_ok=True)
dest_additional.mkdir(parents=True, exist_ok=True)

extracted_new        = set()
extracted_additional = set()
skipped              = []

with zipfile.ZipFile(ZIP_FILE, "r") as zf:
    all_entries = zf.namelist()
    for entry in all_entries:
        # Expected structure: ADNI/<SubjectID>/...
        parts = Path(entry).parts
        if len(parts) < 2 or parts[0] != "ADNI":
            continue
        subject_id = parts[1]

        if subject_id in new_subjects:
            dest = dest_new
            extracted_new.add(subject_id)
        elif subject_id in flagged_subjects:
            dest = dest_additional
            extracted_additional.add(subject_id)
        else:
            skipped.append(subject_id)
            continue

        target = dest / entry
        if entry.endswith("/"):
            target.mkdir(parents=True, exist_ok=True)
        else:
            target.parent.mkdir(parents=True, exist_ok=True)
            with zf.open(entry) as src, open(target, "wb") as dst:
                dst.write(src.read())

print(f"\nExtraction complete.")
print(f"  New subjects extracted          : {len(extracted_new)}")
print(f"  AdditionalTimePoints extracted  : {len(extracted_additional)}")

# Subjects in filtered list but NOT found in zip
missing_new  = new_subjects - extracted_new
missing_add  = flagged_subjects - extracted_additional
if missing_new:
    print(f"\n  WARNING – new subjects not found in zip ({len(missing_new)}):")
    for s in sorted(missing_new):
        print(f"    {s}")
if missing_add:
    print(f"\n  WARNING – flagged subjects not found in zip ({len(missing_add)}):")
    for s in sorted(missing_add):
        print(f"    {s}")

print(f"\n  Output → {EXTRACT_ROOT}")
