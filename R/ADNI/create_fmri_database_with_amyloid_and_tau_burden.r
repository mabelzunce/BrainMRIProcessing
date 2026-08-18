# ============================================================================
#  build_fmri_registry.R
#  ----------------------------------------------------------------------------
#  Build a per-IMAGEUID registry of every fMRI (BOLD) scan in ADNI, enriched
#  with scanner metadata and Mayo-derived slice timing/order where available.
#
#  Output: all_fmri_images.csv  — one row per fMRI IMAGEUID
#
#  Data sources (all from the ADNIMERGE2 R package):
#    - MRIQC                       canonical fMRI image list + scanner metadata
#    - MAYOADIRL_MRI_FMRI_NFQ      Mayo's slice timing (Philips-populated)
#    - MAYOADIRL_MRI_QUALITY_ADNI3 IMAGEUID bridge for ADNI3+
#    - MAYOADIRL_MRI_FMRI          IMAGEUID bridge for ADNI1/GO/2 + SLICEORD
#    - PTDEMOG                     site lookup (RID → SITEID)
# ============================================================================

library(ADNIMERGE2)
library(dplyr)


# ============================================================================
# PART A — Build Mayo-derived slice-timing table (keyed by IMAGEUID)
# ============================================================================

# Per ADNI data dictionary (MAYOADIRL_MRI_FMRI_NFQ.SLICETIMING):
#   "...slices are axial and always counted from inferior to superior.
#    This is invariant to the way the data is stored in DICOM."
# Therefore SLICETIMING_NFQ[i] = absolute acquisition time of slice i
# (slice 1 = most inferior). To get acquisition order:
#   order(times)  →  vector of slice numbers in the order they were acquired.
#
# For Philips, neither SLICEORD nor SLICETIMING_NFQ comes from DICOM — both
# are inferred by Mayo from protocol knowledge. At some sites these two
# fields disagree (see SLICEORD_AGREES). When they do, we prefer
# SLICETIMING_NFQ because its semantics are documented unambiguously and
# its values produce physically-clean inter-slice intervals.

slice_order_from_times <- function(slicetiming_str) {
  if (is.na(slicetiming_str) || !nzchar(slicetiming_str)) return(NA_character_)
  times <- suppressWarnings(as.numeric(strsplit(slicetiming_str, "_")[[1]]))
  if (length(times) == 0 || any(is.na(times))) return(NA_character_)
  # Break ties by slice index so multiband packs list in slice-number order
  ord <- order(times, seq_along(times))
  paste(ord, collapse = "_")
}

# Anchor: every fMRI scan that has slice timing recorded by Mayo
mayo_anchor <- MAYOADIRL_MRI_FMRI_NFQ %>%
  transmute(
    RID, VISCODE, VISCODE2,
    SCANDATE = as.Date(SCANDATE),
    SERIESNUMBER,
    MANUFACTURER, MANUFACTURERSMODELNAME,
    REPETITIONTIME, ECHOTIME,
    SLICETIMING_NFQ = SLICETIMING,
    NFQ, OVERALLQC
  )

# Bridge 1: IMAGEUID for ADNI3+ scans via the QC table
qc_bridge <- MAYOADIRL_MRI_QUALITY_ADNI3 %>%
  filter(
    grepl("fMRI|BOLD|resting|rsfMRI|fcMRI", SERIES_DESCRIPTION, ignore.case = TRUE),
    !is.na(LONI_IMAGE)
  ) %>%
  transmute(
    RID,
    SCANDATE = as.Date(SERIES_DATE),
    IMAGEUID_qc = as.character(LONI_IMAGE),
    SERIES_DESCRIPTION,
    SERIES_QUALITY
  ) %>%
  distinct(RID, SCANDATE, .keep_all = TRUE)

# Bridge 2: IMAGEUID + Mayo's SLICEORD for ADNI1/GO/2
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

# Merge anchor + bridges; coalesce IMAGEUID from the two sources
mayo_merged <- mayo_anchor %>%
  left_join(qc_bridge,   by = c("RID", "SCANDATE")) %>%
  left_join(fmri_bridge, by = c("RID", "SCANDATE")) %>%
  mutate(IMAGEUID = coalesce(IMAGEUID_qc, IMAGEUID_fmri))

# Compute slice-order columns
mayo_merged <- mayo_merged %>%
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
    # Prefer timing-derived order; for disagreements, leave NA so the user
    # decides per-case (see Part C diagnostics).
    SLICEORD_FINAL = case_when(
      SLICEORD_AGREES %in% TRUE  ~ SLICEORD_FROM_TIMES,
      is.na(SLICEORD_NORM)       ~ SLICEORD_FROM_TIMES,
      is.na(SLICEORD_FROM_TIMES) ~ SLICEORD_NORM,
      TRUE                       ~ NA_character_      # conflict → unresolved
    )
  )


# ============================================================================
# PART B — Build the canonical fMRI image list from MRIQC
# ============================================================================

# MRIQC is the most complete source: ~4000+ BOLD series across ADNI1/GO/2/3/4,
# regardless of whether Mayo processed them downstream.

all_fmri <- MRIQC %>%
  filter(SeriesType == "BOLD") %>%
  transmute(
    IMAGEUID = as.character(image_id),
    PTID     = ParticipantID,
    RID      = as.integer(sub(".*_S_", "", ParticipantID)),
    VISCODE2,
    SCANDATE = as.Date(StudyDate),
    SeriesDescription,
    is_multiband = grepl("\\bMB\\b|multiband", SeriesDescription, ignore.case = TRUE),
    is_extended  = grepl("Extended", SeriesDescription, ignore.case = TRUE),
    ScannerManufacturer, ScannerModel, SoftwareVersion,
    MagneticFieldStrength, Acceleration,
    AcquisitionType, AcquisitionPlane,
    NumberVolumes, SlicesPerVolume, SliceThickness,
    ReceiveCoilName,
    StudyInstanceUID, SeriesInstanceUID,
    LONIStudy, LONISeries
  ) %>%
  arrange(RID, SCANDATE, IMAGEUID)


# ============================================================================
# PART C — Enrich the MRIQC list with Mayo metadata, by IMAGEUID
# ============================================================================

mayo_for_join <- mayo_merged %>%
  filter(!is.na(IMAGEUID)) %>%
  select(IMAGEUID,
         REPETITIONTIME_mayo = REPETITIONTIME,
         ECHOTIME_mayo       = ECHOTIME,
         SLICEORD_FINAL, SLICEORD_FROM_TIMES,
         SLICEORD_mayo, SLICEORD_NORM, SLICEORD_AGREES,
         SLICETIMING_NFQ,
         NFQ, OVERALLQC, SERIES_QUALITY,
         MEANTSNR) %>%
  distinct(IMAGEUID, .keep_all = TRUE)

# Add site information from PTDEMOG
site_lookup <- PTDEMOG %>%
  select(RID, SITEID) %>%
  distinct(RID, .keep_all = TRUE)

all_fmri_enriched <- all_fmri %>%
  left_join(mayo_for_join, by = "IMAGEUID") %>%
  left_join(site_lookup,   by = "RID")


# ============================================================================
# PART D — Write outputs
# ============================================================================

write.csv(all_fmri_enriched,
          "all_fmri_images.csv",
          row.names = FALSE)

# Slim version (drop audit columns)
all_fmri_enriched %>%
  select(-SLICEORD_NORM, -SLICEORD_AGREES) %>%
  write.csv("all_fmri_images_slim.csv", row.names = FALSE)


# ============================================================================
# PART E — Sanity reports
# ============================================================================

cat("\n========== OVERALL COVERAGE ==========\n")
cat("Total BOLD series:        ", nrow(all_fmri_enriched), "\n")
cat("Unique subjects:          ", n_distinct(all_fmri_enriched$RID), "\n")
cat("Multiband series:         ", sum(all_fmri_enriched$is_multiband), "\n")
cat("Extended (long) protocol: ", sum(all_fmri_enriched$is_extended), "\n\n")

cat("========== BY SCANNER MANUFACTURER ==========\n")
print(all_fmri_enriched %>% count(ScannerManufacturer) %>% arrange(desc(n)))

cat("\n========== MAYO METADATA COVERAGE ==========\n")
print(
  all_fmri_enriched %>%
    summarise(
      total              = n(),
      with_slice_timing  = sum(!is.na(SLICETIMING_NFQ)),
      with_slice_order   = sum(!is.na(SLICEORD_FINAL)),
      with_mayo_qc       = sum(!is.na(OVERALLQC)),
      with_neither       = sum(is.na(SLICETIMING_NFQ) & is.na(SLICEORD_FINAL))
    )
)

cat("\n========== SLICETIMING COVERAGE BY MANUFACTURER ==========\n")
print(
  all_fmri_enriched %>%
    mutate(has_timing = !is.na(SLICETIMING_NFQ)) %>%
    count(ScannerManufacturer, has_timing) %>%
    arrange(ScannerManufacturer, desc(has_timing))
)

cat("\n========== SLICEORD AGREEMENT (Mayo vs timing-derived) ==========\n")
print(table(all_fmri_enriched$SLICEORD_AGREES, useNA = "ifany"))

cat("\n========== DISAGREEMENTS BY SITE + MODEL ==========\n")
disagreeing_models <- all_fmri_enriched %>%
  filter(SLICEORD_AGREES == FALSE) %>%
  pull(ScannerModel) %>%
  unique()

if (length(disagreeing_models) > 0) {
  print(
    all_fmri_enriched %>%
      filter(ScannerModel %in% disagreeing_models,
             !is.na(SLICEORD_FROM_TIMES), !is.na(SLICEORD_NORM)) %>%
      group_by(ScannerModel, SITEID) %>%
      summarise(
        n             = n(),
        n_agreeing    = sum(SLICEORD_AGREES, na.rm = TRUE),
        n_disagreeing = sum(!SLICEORD_AGREES, na.rm = TRUE),
        pct_agree     = round(100 * n_agreeing / n, 1),
        .groups = "drop"
      ) %>%
      arrange(ScannerModel, desc(n_disagreeing)),
    n = 100
  )
} else {
  cat("(No disagreements found)\n")
}

cat("\n========== DUPLICATE IMAGEUIDS (should be 0) ==========\n")
n_dup <- all_fmri_enriched %>%
  filter(!is.na(IMAGEUID)) %>%
  count(IMAGEUID) %>%
  filter(n > 1) %>%
  nrow()
cat("Number of duplicated IMAGEUIDs:", n_dup, "\n")

cat("\nDone. Wrote:\n")
cat("  all_fmri_images.csv      (full, with audit columns)\n")
cat("  all_fmri_images_slim.csv (without audit columns)\n")


# ============================================================================
# PART F — Add amyloid PET status (A+/A−) matched to nearest fMRI scan
# ============================================================================

amy <- UCBERKELEY_AMY_6MM %>%
  transmute(
    RID,
    AMY_SCANDATE = as.Date(SCANDATE),
    AMYLOID_STATUS,               # 0 = A-, 1 = A+
    AMYLOID_STATUS_COMPOSITE_REF,
    CENTILOIDS,
    SUMMARY_SUVR,
    AMY_TRACER = TRACER,
    AMY_QC = qc_flag
  ) %>%
  filter(!is.na(AMYLOID_STATUS))

# For each fMRI scan, find nearest amyloid PET (within 24 months)
fmri_x_amy <- all_fmri_enriched %>%
  select(IMAGEUID, RID, SCANDATE) %>%
  inner_join(amy, by = "RID", relationship = "many-to-many") %>%
  mutate(days_diff = as.numeric(AMY_SCANDATE - SCANDATE)) %>%
  filter(abs(days_diff) <= 730) %>%
  group_by(IMAGEUID) %>%
  slice_min(abs(days_diff), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(IMAGEUID,
         AMY_SCANDATE, days_diff_amy = days_diff,
         AMYLOID_STATUS, AMYLOID_STATUS_COMPOSITE_REF,
         CENTILOIDS, SUMMARY_SUVR, AMY_TRACER, AMY_QC)

all_fmri_enriched <- all_fmri_enriched %>%
  left_join(fmri_x_amy, by = "IMAGEUID")

# Coverage report
cat("\n========== AMYLOID PET COVERAGE ==========\n")
cat("fMRI scans with amyloid within 2 years:",
    sum(!is.na(all_fmri_enriched$AMYLOID_STATUS)), "/",
    nrow(all_fmri_enriched), "\n\n")
print(
  all_fmri_enriched %>%
    filter(!is.na(AMYLOID_STATUS)) %>%
    count(AMYLOID_STATUS) %>%
    mutate(pct = round(100 * n / sum(n), 1))
)


# ============================================================================
# PART G — Add tau PET status (T+/T−) matched to nearest fMRI scan
# ============================================================================

# UCBERKELEY_TAU_6MM: FTP (and later MK-6240) with inferior cerebellum reference.
# No pre-computed positivity flag — we apply a threshold to META_TEMPORAL_SUVR.
# Default cutpoint 1.32 (Jack et al. 2017). Change TAU_CUTPOINT to explore.

TAU_CUTPOINT <- 1.32

tau <- UCBERKELEY_TAU_6MM %>%
  transmute(
    RID,
    TAU_SCANDATE = as.Date(SCANDATE),
    META_TEMPORAL_SUVR,
    CTX_ENTORHINAL_SUVR,       # Braak I proxy — early tau
    INFERIORCEREBELLUM_SUVR,   # reference (should be ~1)
    TAU_TRACER = TRACER,
    TAU_QC = qc_flag
  ) %>%
  filter(!is.na(META_TEMPORAL_SUVR)) %>%
  # Only FTP for now — MK-6240 SUVR values are on a different scale and
  # would need a different cutpoint (do not mix without harmonization).
  filter(TAU_TRACER == "FTP") %>%
  mutate(TAU_STATUS = as.integer(META_TEMPORAL_SUVR > TAU_CUTPOINT))

fmri_x_tau <- all_fmri_enriched %>%
  select(IMAGEUID, RID, SCANDATE) %>%
  inner_join(tau, by = "RID", relationship = "many-to-many") %>%
  mutate(days_diff = as.numeric(TAU_SCANDATE - SCANDATE)) %>%
  filter(abs(days_diff) <= 730) %>%
  group_by(IMAGEUID) %>%
  slice_min(abs(days_diff), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(IMAGEUID,
         TAU_SCANDATE, days_diff_tau = days_diff,
         TAU_STATUS,
         META_TEMPORAL_SUVR, CTX_ENTORHINAL_SUVR,
         TAU_TRACER, TAU_QC)

all_fmri_enriched <- all_fmri_enriched %>%
  left_join(fmri_x_tau, by = "IMAGEUID")

# Coverage report
cat("\n========== TAU PET COVERAGE (FTP, cutpoint =", TAU_CUTPOINT, ") ==========\n")
cat("fMRI scans with tau within 2 years:",
    sum(!is.na(all_fmri_enriched$TAU_STATUS)), "/",
    nrow(all_fmri_enriched), "\n\n")
print(
  all_fmri_enriched %>%
    filter(!is.na(TAU_STATUS)) %>%
    count(TAU_STATUS) %>%
    mutate(pct = round(100 * n / sum(n), 1))
)

all_fmri_enriched <- all_fmri_enriched %>%
  mutate(
    AT_STATUS = case_when(
      is.na(AMYLOID_STATUS) | is.na(TAU_STATUS)   ~ NA_character_,
      AMYLOID_STATUS == 0 & TAU_STATUS == 0       ~ "A-T-",
      AMYLOID_STATUS == 1 & TAU_STATUS == 0       ~ "A+T-",
      AMYLOID_STATUS == 0 & TAU_STATUS == 1       ~ "A-T+",   # suspected non-AD tauopathy
      AMYLOID_STATUS == 1 & TAU_STATUS == 1       ~ "A+T+"
    )
  )

cat("\n========== AT COMBINED STATUS ==========\n")
print(
  all_fmri_enriched %>%
    filter(!is.na(AT_STATUS)) %>%
    count(AT_STATUS) %>%
    mutate(pct = round(100 * n / sum(n), 1))
)

# Write enriched fMRI database to CSV
write.csv(all_fmri_enriched, 
          file = "fmri_database_enriched.csv", 
          row.names = FALSE)
cat("\nEnriched fMRI database written to: fmri_database_enriched.csv\n")