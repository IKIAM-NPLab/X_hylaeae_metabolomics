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



























################################# MZmine #######################################
##Feature list deconvoluted with MZmine (3 to 43 min dataset)
mzmine_rt <- read_excel("../B_grandiflora_metabolomics/Data/B_gradiflora_MZmine_Feat_list_3to43min.xlsx", sheet = 1)
##RI calculation (3 to 43 min dataset)
mzmine_ri <- indexRtime(mzmine_rt$`row retention time`, alkane_data)
##Adding the RI to the MZmine feature list (3 to 43 min dataset)
mzmine_rt$RI <- mzmine_ri 
## Moving the RI column close to the RT column (3 to 43 min dataset)
mzmine_rt <- mzmine_rt[,c(1:2, 47, 3:46)]
##Exporting the MZmine feature list with RI
write_xlsx(mzmine_rt, "../B_grandiflora_metabolomics/Result/B_gradiflora_MZmine_Feat_list_3to43min_RI.xlsx")

## Experimental Retention Index (RI) calculation
##Feature list deconvoluted with MZmine (43 to 53 min dataset)
mzmine_rt43to53 <- read_excel("../B_grandiflora_metabolomics/Data/B_gradiflora_MZmine_Feat_list_43to53min.xlsx", sheet = 1)
##RI calculation (43 to 53 min dataset)
mzmine_ri43to53 <- indexRtime(mzmine_rt43to53$`row retention time`, alkane_data)
##Adding the RI to the MZmine feature list (43 to 53 min dataset)
mzmine_rt43to53$RI <- mzmine_ri43to53 
## Moving the RI column close to the RT column (43 to 53 min dataset)
mzmine_rt43to53 <- mzmine_rt43to53[,c(1:3, 48, 4:47)]
##Exporting the MZmine feature list with RI
write_xlsx(mzmine_rt43to53, "../B_grandiflora_metabolomics/Result/B_gradiflora_MZmine_Feat_list_43to53min_RI.xlsx")

## Experimental Retention Index (RI) calculation
##Feature list deconvoluted with MZmine (53 to 58 min dataset)
mzmine_rt53to58 <- read_excel("E:/B_grandiflora/B_grandiflora_metabolomics/Data/B_gradiflora_MZmine_Feat_list_53to58min.xlsx", sheet = 1)
##RI calculation (53 to 58 min dataset)
mzmine_ri53to58 <- indexRtime(mzmine_rt53to58$`row retention time`, alkane_data)
##Adding the RI to the MZmine feature list (53 to 58 min dataset)
mzmine_rt53to58$RI <- mzmine_ri53to58 
## Moving the RI column close to the RT column (53 to 58 min dataset)
mzmine_rt53to58 <- mzmine_rt53to58[,c(1:3, 48, 4:47)]
##Exporting the MZmine feature list with RI
write_xlsx(mzmine_rt53to58, "E:/B_grandiflora/B_grandiflora_metabolomics/Result/B_gradiflora_MZmine_Feat_list_53to58min_RI.xlsx")


################################# eRah #########################################
## Experimental Retention Index (RI) calculation
##Feature list deconvoluted with eRah (3 to 43 min dataset)
erah_rt3to43 <- read_excel("E:/B_grandiflora/B_grandiflora_metabolomics/Data/B_gradiflora_eRah_Feat_list_3to43min.xlsx", sheet = 1)
##RI calculation (3 to 43 min dataset)
erah_ri3to43 <- indexRtime(erah_rt3to43$tmean, alkane_data)
##Adding the RI to the eRah feature list (3 to 43 min dataset)
erah_rt3to43$RI <- erah_ri3to43
## Moving the RI column close to the RT column (3 to 43 min dataset)
erah_rt3to43 <- erah_rt3to43[,c(1:4, 41, 5:40)]
##Exporting the eRah feature list with RI
write_xlsx(erah_rt3to43, "E:/B_grandiflora/B_grandiflora_metabolomics/Result/B_gradiflora_eRah_Feat_list_3to43min_RI.xlsx")

## Experimental Retention Index (RI) calculation
##Feature list deconvoluted with eRah (43 to 53 min dataset)
erah_rt43to53 <- read_excel("E:/B_grandiflora/B_grandiflora_metabolomics/Data/B_gradiflora_eRah_Feat_list_43to53.xlsx", sheet = 1)
##RI calculation (43 to 53 min dataset)
erah_ri43to53 <- indexRtime(erah_rt43to53$tmean, alkane_data)
##Adding the RI to the eRah feature list (43 to 53 min dataset)
erah_rt43to53$RI <- erah_ri43to53
## Moving the RI column close to the RT column (43 to 53 min dataset)
erah_rt43to53 <- erah_rt43to53[,c(1:4, 40, 5:39)]
##Exporting the eRah feature list with RI
write_xlsx(erah_rt43to53, "E:/B_grandiflora/B_grandiflora_metabolomics/Result/B_gradiflora_eRah_Feat_list_43to53min_RI.xlsx")

## Experimental Retention Index (RI) calculation
##Feature list deconvoluted with eRah (53 to 58 min dataset)
erah_rt53to58 <- read_excel("E:/B_grandiflora/B_grandiflora_metabolomics/Data/B_gradiflora_eRah_Feat_list_53to58.xlsx", sheet = 1)
##RI calculation (53 to 58 min dataset)
erah_ri53to58 <- indexRtime(erah_rt53to58$tmean, alkane_data)
##Adding the RI to the eRah feature list (53 to 58 min dataset)
erah_rt53to58$RI <- erah_ri53to58
## Moving the RI column close to the RT column (53 to 58 min dataset)
erah_rt53to58 <- erah_rt53to58[,c(1:4, 39, 5:38)]
##Exporting the eRah feature list with RI
write_xlsx(erah_rt53to58, "E:/B_grandiflora/B_grandiflora_metabolomics/Result/B_gradiflora_eRah_Feat_list_53to58min_RI.xlsx")

################################# MSHub ########################################
## Experimental Retention Index (RI) calculation
##Feature list deconvoluted with MShub
mshub_rt <- read_excel("E:/B_grandiflora/B_grandiflora_metabolomics/Data/B_gradiflora_MSHub_Feat_list.xlsx", sheet = 1)
##RI calculation
mshub_ri <- indexRtime(mshub_rt$`row retention time`, alkane_data)
##Adding the RI to the MSHub feature list
mshub_rt$RI <- mshub_ri
## Moving the RI column close to the RT column
mshub_rt <- mshub_rt[,c(1:4, 39, 5:38)]
##Exporting the eRah feature list with RI
write_xlsx(mshub_rt, "E:/B_grandiflora/B_grandiflora_metabolomics/Result/B_gradiflora_MSHub_Feat_list_RI.xlsx")

