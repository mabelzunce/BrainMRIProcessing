# ============================================================================
#  Build a per-IMAGEUID fMRI metadata CSV from ADNIMERGE2
#  Combines slice timing / order / scanner info, keyed by IMAGEUID
# ============================================================================

library(ADNIMERGE2)
library(dplyr)

# ----------------------------------------------------------------------------
# 1. Anchor: every fMRI scan that has slice timing recorded
# ----------------------------------------------------------------------------
# MAYOADIRL_MRI_FMRI_NFQ is the most complete fMRI table (~3000 rows),
# but lacks IMAGEUID — we'll recover it from two bridge tables below.

anchor <- MAYOADIRL_MRI_FMRI_NFQ %>%
  transmute(
    RID, VISCODE, VISCODE2,
    SCANDATE = as.Date(SCANDATE),
    SERIESNUMBER,
    MANUFACTURER, MANUFACTURERSMODELNAME,
    REPETITIONTIME, ECHOTIME,
    SLICETIMING_NFQ = SLICETIMING,
    NFQ, OVERALLQC
  )

# ----------------------------------------------------------------------------
# 2. Bridge 1: MAYOADIRL_MRI_QUALITY_ADNI3 → IMAGEUID for ADNI3 scans
# ----------------------------------------------------------------------------
# LONI_IMAGE is the IMAGEUID. Filter to fMRI series so structural scans don't
# get joined by accident.

qc_bridge <- MAYOADIRL_MRI_QUALITY_ADNI3 %>%
  filter(
    grepl("fMRI|BOLD|resting|rsfMRI|fcMRI", SERIES_DESCRIPTION, ignore.case = TRUE),
    !is.na(LONI_IMAGE)
  ) %>%
  transmute(
    RID,
    SCANDATE = as.Date(SERIES_DATE),
    IMAGEUID_qc = as.character(LONI_IMAGE),
    SERIES_DESCRIPTION, SERIES_QUALITY
  ) %>%
  distinct(RID, SCANDATE, .keep_all = TRUE)

# ----------------------------------------------------------------------------
# 3. Bridge 2: MAYOADIRL_MRI_FMRI → IMAGEUID for ADNI1 / ADNI-GO / ADNI2
# ----------------------------------------------------------------------------
# This table only has ~900 rows, but covers the older phases that ADNI3 QC
# misses. It also carries Mayo's SLICEORD field.

fmri_bridge <- MAYOADIRL_MRI_FMRI %>%
  transmute(
    RID,
    SCANDATE = as.Date(SCANDATE),
    IMAGEUID_fmri = as.character(IMAGEUID),
    SLICEORD_mayo = SLICEORD,
    MAGSTRENG,
    MEANTSNR
  ) %>%
  distinct(RID, SCANDATE, .keep_all = TRUE)

# ----------------------------------------------------------------------------
# 4. Merge anchor + bridges, then coalesce IMAGEUID from the two sources
# ----------------------------------------------------------------------------

merged <- anchor %>%
  left_join(qc_bridge,   by = c("RID", "SCANDATE")) %>%
  left_join(fmri_bridge, by = c("RID", "SCANDATE")) %>%
  mutate(IMAGEUID = coalesce(IMAGEUID_qc, IMAGEUID_fmri))

# ----------------------------------------------------------------------------
# 5. Slice-order columns
# ----------------------------------------------------------------------------
# SLICEORD_mayo (from MAYOADIRL_MRI_FMRI) — Mayo's reported order, may be empty
#                or, at some sites, disagrees with the actual timings.
# SLICEORD_FROM_TIMES — derived from SLICETIMING_NFQ. Times are indexed by
#                slice number (times[i] = time slice i was acquired), so
#                order(times) gives the true acquisition order.
# SLICEORD_FINAL — prefer the timing-derived order (ground truth from
#                physical timestamps); fall back to Mayo's where times missing.

slice_order_from_times <- function(slicetiming_str) {
  if (is.na(slicetiming_str) || !nzchar(slicetiming_str)) return(NA_character_)
  times <- suppressWarnings(as.numeric(strsplit(slicetiming_str, "_")[[1]]))
  if (length(times) == 0 || any(is.na(times))) return(NA_character_)
  # Break ties by slice index so multiband packs list in slice-number order
  ord <- order(times, seq_along(times))
  paste(ord, collapse = "_")
}

merged <- merged %>%
  rowwise() %>%
  mutate(
    SLICEORD_NORM = if (is.na(SLICEORD_mayo) || !nzchar(SLICEORD_mayo))
                      NA_character_
                    else
                      gsub("\\s+", "_", trimws(SLICEORD_mayo)),
    SLICEORD_FROM_TIMES = slice_order_from_times(SLICETIMING_NFQ)
  ) %>%
  ungroup() %>%
  mutate(
    SLICEORD_AGREES = case_when(
      is.na(SLICEORD_NORM) | is.na(SLICEORD_FROM_TIMES) ~ NA,
      SLICEORD_NORM == SLICEORD_FROM_TIMES              ~ TRUE,
      TRUE                                              ~ FALSE
    ),
    SLICEORD_FINAL = coalesce(SLICEORD_FROM_TIMES, SLICEORD_NORM)
  )

# ----------------------------------------------------------------------------
# 6. Tidy column order, write CSV
# ----------------------------------------------------------------------------

final <- merged %>%
  select(
    RID, VISCODE, VISCODE2, SCANDATE, IMAGEUID,
    MANUFACTURER, MANUFACTURERSMODELNAME,
    MAGSTRENG, REPETITIONTIME, ECHOTIME,
    SLICEORD_FINAL, SLICEORD_FROM_TIMES, SLICEORD_NORM, SLICEORD_AGREES,
    SLICETIMING_NFQ,
    NFQ, OVERALLQC, SERIES_QUALITY, MEANTSNR,
    SERIES_DESCRIPTION,
    IMAGEUID_qc, IMAGEUID_fmri   # keep for auditing; drop later if desired
  )

write.csv(final, "fmri_slice_timing_with_imageid.csv", row.names = FALSE)

# ----------------------------------------------------------------------------
# 7. Sanity report
# ----------------------------------------------------------------------------

cat("\n========== COVERAGE REPORT ==========\n")
cat("Total rows                     :", nrow(final), "\n")
cat("Have IMAGEUID                  :", sum(!is.na(final$IMAGEUID)), "\n")
cat("Have slice timing              :", sum(!is.na(final$SLICETIMING_NFQ)), "\n")
cat("Have slice order (any source)  :", sum(!is.na(final$SLICEORD_FINAL)), "\n")
cat("Have slice order from times    :", sum(!is.na(final$SLICEORD_FROM_TIMES)), "\n")
cat("Have Mayo's SLICEORD           :", sum(!is.na(final$SLICEORD_NORM)), "\n\n")

cat("========== SLICEORD AGREEMENT ==========\n")
print(table(final$SLICEORD_AGREES, useNA = "ifany"))

cat("\n========== DISAGREEMENTS BY MANUFACTURER ==========\n")
print(
  final %>%
    filter(SLICEORD_AGREES == FALSE) %>%
    count(MANUFACTURER, MANUFACTURERSMODELNAME) %>%
    arrange(desc(n))
)

cat("\n========== DUPLICATE IMAGEUIDS (should be 0) ==========\n")
print(
  final %>%
    filter(!is.na(IMAGEUID)) %>%
    count(IMAGEUID) %>%
    filter(n > 1) %>%
    nrow()
)

cat("\nWrote: fmri_slice_timing_with_imageid.csv\n")