# 1. Dependencies (one-time)
#install.packages("Hmisc")
#install.packages(c("dplyr", "tibble", "assertr", "tidyselect",
#                   "purrr", "tidyr", "readr", "bookdown", "Hmisc"))

# 2. Install the package from the file you downloaded
#install.packages("/home/martin/data_imaging/ADNIdata/StudyInfo_2026/ADNIMERGE2.tar.gz", repos = NULL, type = "source")

# 3. Load it and export the merged table
library(ADNIMERGE2)
library(dplyr)
# Confirm package is actually attached
search()                               # should include "package:ADNIMERGE2"
"ADNIMERGE2" %in% loadedNamespaces()

# Try alternative ways to list datasets
data(package = "ADNIMERGE2")$results[, "Item"]

all_objs <- ls("package:ADNIMERGE2")
# Filter to likely imaging tables
all_objs[grepl("MRI|IMAGE|FS|UCSF|FBB|AV45|FDG|PET|TAU", all_objs, ignore.case = TRUE)]

write.csv(ADSL, file = "ADSL.csv", row.names = FALSE)

# Imaging metadata (one row per scan)
#mri <- MRILIST %>%
#  select(RID, VISCODE, EXAMDATE, IMAGEUID, SCANDATE, FIELD_STRENGTH, SEQUENCE)
# adjust column names to what names(MRILIST) actually shows

# Longitudinal diagnosis at that visit
#dx <- DXSUM %>% select(RID, VISCODE, DIAGNOSIS)

# Subject-level demographics from the CSV you already have
#adsl <- ADSL %>%
#  select(SUBJID, AGE, SEX, EDUC, RACE, ETHNIC, APOE, ORIGPROT) %>%
#  mutate(RID = as.integer(SUBJID))

#merged <- mri %>%
#  left_join(dx,   by = c("RID", "VISCODE")) %>%
#  left_join(adsl, by = "RID")

#write.csv(merged, "ADNIMERGE_by_image.csv", row.names = FALSE)

slice_timing %>%
  count(ScannerManufacturer,
        has_slice_timing = !is.na(SliceTiming_MRINFQ) | !is.na(SLICETIMING_NFQ))
# 1. Primary anchor — one row per fMRI IMAGEUID
fmri <- MAYOADIRL_MRI_FMRI %>%
  select(RID, VISCODE, VISCODE2, SCANDATE, IMAGEUID,
         MAGSTRENG, SERDESC, SLICEORD, PHASEDIR,
         MEANTSNR, MEDTSNR, SDTSNR, STATUS)

# 2. ADNI3 QC table — filter to fMRI series only, key on LONI_IMAGE
qc <- MAYOADIRL_MRI_QUALITY_ADNI3 %>%
  filter(grepl("fMRI|BOLD|resting|rsfMRI", SERIES_DESCRIPTION, ignore.case = TRUE)) %>%
  transmute(RID,
            IMAGEUID = as.character(LONI_IMAGE),
            SERIES_DESCRIPTION, SERIES_QUALITY, SERIES_SELECTED,
            FIELD_STRENGTH, SERIES_DATE)

# 3. MRINFQ — rich scanner + SliceTiming, key on LONIImage
mrinfq <- MRINFQ %>%
  transmute(RID, VISCODE,
            IMAGEUID = as.character(LONIImage),
            ScannerManufacturer, ScannerModel, SoftwareVersion,
            MagneticFieldStrength, Acceleration,
            AcquisitionType, AcquisitionPlane,
            RepetitionTime, EchoTime,
            NumberVolumes, SlicesPerVolume, SliceThickness,
            SliceTiming_MRINFQ = SliceTiming)

# 4. Earlier-phase NFQ — no IMAGEUID, join by RID + SCANDATE later
fmri_nfq <- MAYOADIRL_MRI_FMRI_NFQ %>%
  select(RID, VISCODE, SCANDATE, SERIESNUMBER,
         MANUFACTURER, MANUFACTURERSMODELNAME,
         REPETITIONTIME, ECHOTIME,
         SLICETIMING_NFQ = SLICETIMING) %>%
  # collapse duplicates if multiple series per session
  distinct(RID, SCANDATE, .keep_all = TRUE)

# Make IMAGEUID types match before joining
fmri  <- fmri  %>% mutate(IMAGEUID = as.character(IMAGEUID))

# Merge
slice_timing <- fmri %>%
  left_join(qc,       by = c("RID", "IMAGEUID")) %>%
  left_join(mrinfq,   by = c("RID", "VISCODE", "IMAGEUID")) %>%
  left_join(fmri_nfq, by = c("RID", "VISCODE", "SCANDATE"))

cat("rows:", nrow(slice_timing),
    "  unique IMAGEUIDs:", length(unique(slice_timing$IMAGEUID)), "\n")

write.csv(slice_timing, "fmri_slice_timing.csv", row.names = FALSE)