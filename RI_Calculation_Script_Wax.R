####### Script to calculate RI #######

# Installation of R package to calculate linear retention index (RI).

# Installation of "MetaboCoreUtils" package
#install.packages("remotes")
#remotes::install_github("rformassspectrometry/MetaboCoreUtils")

# Loading "MetaboCoreUtils" library
library("MetaboCoreUtils")

# Extra libraries
library("readxl")
library("writexl")

# Experimental RI for features that match with NIST libraries
## Loadding the retention time (RT) of each n-alkane
alkane_rt_wx <- c(3.870, 5.412, 7.318, 9.469, 11.664, 13.870, 16.031, 18.122,
                  20.139, 22.076, 23.937, 25.728, 27.624, 29.995, 33.071, 36.990,
                  40.020, 42.278, 44.163, 45.958, 48.044, 50.540)
## Loadding the RI of each n-alkane
alkane_ri_wx <- c(1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100,
                  2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000, 3100,
                  3200, 3300)
## n-alkanes DataFrame
alkane_data_wx <- data.frame(rtime = alkane_rt_wx, rindex = alkane_ri_wx)

## Experimental Retention Index (RI) calculation
################################# MS-DIAL ######################################
##Feature list deconvoluted with MS-DIAL
msdial3to38min_wx <- read_excel("../X_hylaeae_metabolomics/Data/MSDIAL_GCMS_polar_feat_list_3-38min.xlsx", sheet = 2)
msdial3to38min_wx <- data.frame(Alignment_ID = msdial3to38min_wx$`Alignment ID`,
                                RT = msdial3to38min_wx$`Average Rt(min)`,
                                MSDIAL_RI = msdial3to38min_wx$`Average RI`,
                                Exp_EI = msdial3to38min_wx$`EI spectrum`,
                                RI_tag = msdial3to38min_wx$Comment)
msdial38to50min_wx <- read_excel("../X_hylaeae_metabolomics/Data/MSDIAL_GCMS_polar_feat_list_38-50min.xlsx", sheet = 2)
msdial38to50min_wx <- data.frame(Alignment_ID = msdial38to50min_wx$`Alignment ID`,
                                 RT = msdial38to50min_wx$`Average Rt(min)`,
                                 MSDIAL_RI = msdial38to50min_wx$`Average RI`,
                                 Exp_EI = msdial38to50min_wx$`EI spectrum`,
                                 RI_tag = msdial38to50min_wx$Comment)
##Merge files
msdial_rtwx <- bind_rows(
  msdial3to38min_wx %>% mutate(Segment = "3-38min"),
  msdial38to50min_wx %>% mutate(Segment = "38-50min"))
# Adding a global ID
msdial_rtwx <- msdial_rtwx %>%
  mutate(Feature_ID = row_number())
##RI calculation
msdial_riwx <- indexRtime(msdial_rtwx$RT, alkane_data_wx)
##Adding the RI to the MS-DIAL feature list
msdial_rtwx$R_RI <- msdial_riwx
## Moving columns
msdial_rtwx <- msdial_rtwx[,c(7, 6, 1:2, 8, 3:5)]
##Exporting the MS-DIAL feature list with RI
write_xlsx(msdial_rtwx, "../X_hylaeae_metabolomics/Result/MSDIAL_GCMS_polar_feat_list_with_RI.xlsx")

################################# eRah ######################################
##Feature list deconvoluted with eRah
erah3to38min_wx <- read.csv("../X_hylaeae_metabolomics/Data/eRah_GCMS_polar_feat_list_3to38min.csv",
                            sep = ",", check.names = FALSE)
erah3to38min_wx <- data.frame(Alignment_ID = erah3to38min_wx$AlignID,
                              RT = erah3to38min_wx$tmean,
                              Exp_EI = erah3to38min_wx$Spectra)
erah3to38min_wx$RI_tag <- NA_character_
erah38to50min_wx <- read.csv("../X_hylaeae_metabolomics/Data/eRah_GCMS_polar_feat_list_38to50min.csv",
                             sep = ",", check.names = FALSE)
erah38to50min_wx <- data.frame(Alignment_ID = erah38to50min_wx$AlignID,
                               RT = erah38to50min_wx$tmean,
                               Exp_EI = erah38to50min_wx$Spectra)
erah38to50min_wx$RI_tag <- NA_character_
##Merge files
erah_rtwx <- bind_rows(
  erah3to38min_wx %>% mutate(Segment = "3-38min"),
  erah38to50min_wx %>% mutate(Segment = "38-50min"))
# Adding a global ID
erah_rtwx <- erah_rtwx %>%
  mutate(Feature_ID = row_number())
##RI calculation
erah_riwx <- indexRtime(erah_rtwx$RT, alkane_data_wx)
##Adding the RI to the eRah feature list
erah_rtwx$R_RI <- erah_riwx
## Moving columns
erah_rtwx <- erah_rtwx[,c(6, 5, 1:2, 7, 3:4)]
##Exporting the eRah feature list with RI
write_xlsx(erah_rtwx, "../X_hylaeae_metabolomics/Result/erah_GCMS_polar_feat_list_with_RI_EI.xlsx")

################################# MZmine ######################################
##Feature list deconvoluted with MZmine
mzmine3to38min_wx <- read.csv("../X_hylaeae_metabolomics/Data/MZmine_GCMS_polar_feat_list_3-38min.csv",
                              sep = ",", check.names = FALSE)
mzmine3to38min_wx <- data.frame(Alignment_ID = mzmine3to38min_wx$`row ID`,
                                RT = mzmine3to38min_wx$`row retention time`)
mzmine3to38min_wx$RI_tag <- NA
mzmine38to50min_wx <- read.csv("../X_hylaeae_metabolomics/Data/MZmine_GCMS_polar_feat_list_38-50min.csv",
                               sep = ",", check.names = FALSE)
mzmine38to50min_wx <- data.frame(Alignment_ID = mzmine38to50min_wx$`row ID`,
                                 RT = mzmine38to50min_wx$`row retention time`)
mzmine38to50min_wx$RI_tag <- NA
##Merge files
mzmine_rtwx <- bind_rows(
  mzmine3to38min_wx %>% mutate(Segment = "3-38min"),
  mzmine38to50min_wx %>% mutate(Segment = "38-50min"))
# Adding a global ID
mzmine_rtwx <- mzmine_rtwx %>%
  mutate(Feature_ID = row_number())
##RI calculation
mzmine_riwx <- indexRtime(mzmine_rtwx$RT, alkane_data_wx)
##Adding the RI to the MZmine feature list
mzmine_rtwx$R_RI <- mzmine_riwx
## Moving columns
mzmine_rtwx <- mzmine_rtwx[,c(5, 4, 1:2, 6, 3)]
##Exporting the MZmine feature list with RI
write_xlsx(mzmine_rtwx, "../X_hylaeae_metabolomics/Result/MZmine_GCMS_polar_feat_list_with_RI.xlsx")

################################# MSHub ######################################
##Feature list deconvoluted with MSHub
mshub_wx <- read.csv("../X_hylaeae_metabolomics/Data/MSHub_GCMS_polar_feat_list.csv",
                     sep = ",", check.names = FALSE)
mshub_wx <- data.frame(Alignment_ID = mshub_wx$`row ID`,
                       RT = mshub_wx$`row retention time`)
mshub_wx$RI_tag <- NA_character_
# Adding a global ID
mshub_rtwx <- mshub_wx %>%
  mutate(Feature_ID = row_number())
##RI calculation
mshub_riwx <- indexRtime(mshub_rtwx$RT, alkane_data_wx)
##Adding the RI to the MSHub feature list
mshub_rtwx$R_RI <- mshub_riwx
## Moving columns
mshub_rtwx <- mshub_rtwx[,c(4, 1:2, 5, 3)]
##Exporting the MSHub feature list with RI
write_xlsx(mshub_rtwx, "../X_hylaeae_metabolomics/Result/MSHub_GCMS_polar_feat_list_with_RI.xlsx")

