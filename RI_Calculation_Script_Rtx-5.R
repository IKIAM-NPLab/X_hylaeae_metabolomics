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
alkane_rt <- c(3.448, 5.066, 7.205, 9.636, 12.181, 14.707, 17.168, 19.502,
               21.744, 23.880, 25.915, 28.150, 31.073, 35.044, 39.108, 41.831,
               43.964, 45.754, 47.323, 48.739, 50.040, 51.267, 52.414, 53.565,
               54.862)
## Loadding the RI of each n-alkane
alkane_ri <- c(900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800,
               1900, 2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800,
               2900, 3000, 3100, 3200, 3300)
## n-alkanes DataFrame
alkane_data <- data.frame(rtime = alkane_rt, rindex = alkane_ri)


## Experimental Retention Index (RI) calculation
################################# MS-DIAL ######################################
##Feature list deconvoluted with MS-DIAL
msdial3to49min <- read_excel("../X_hylaeae_metabolomics/Data/MSDIAL_feat_list_3-49min.xlsx", sheet = 2)
msdial3to49min <- data.frame(Alignment_ID = msdial3to49min$`Alignment ID`,
                             RT = msdial3to49min$`Average Rt(min)`,
                             MSDIAL_RI = msdial3to49min$`Average RI`,
                             Exp_EI = msdial3to49min$`EI spectrum`,
                             RI_tag = msdial3to49min$Comment)
msdial49to54min <- read_excel("../X_hylaeae_metabolomics/Data/MSDIAL_feat_list_49-54min.xlsx", sheet = 2)
msdial49to54min <- data.frame(Alignment_ID = msdial49to54min$`Alignment ID`,
                             RT = msdial49to54min$`Average Rt(min)`,
                             MSDIAL_RI = msdial49to54min$`Average RI`,
                             Exp_EI = msdial49to54min$`EI spectrum`,
                             RI_tag = msdial49to54min$Comment)
msdial54to56min <- read_excel("../X_hylaeae_metabolomics/Data/MSDIAL_feat_list_54-56min.xlsx", sheet = 2)
msdial54to56min <- data.frame(Alignment_ID = msdial54to56min$`Alignment ID`,
                              RT = msdial54to56min$`Average Rt(min)`,
                              MSDIAL_RI = msdial54to56min$`Average RI`,
                              Exp_EI = msdial54to56min$`EI spectrum`,
                              RI_tag = msdial54to56min$Comment)
##Merge files
msdial_rt <- bind_rows(
  msdial3to49min %>% mutate(Segment = "3-49min"),
  msdial49to54min %>% mutate(Segment = "49-54min"),
  msdial54to56min %>% mutate(Segment = "54-56min"))
# Adding a global ID
msdial_rt <- msdial_rt %>%
  mutate(Feature_ID = row_number())
##RI calculation
msdial_ri <- indexRtime(msdial_rt$RT, alkane_data)
##Adding the RI to the MS-DIAL feature list
msdial_rt$R_RI <- msdial_ri
## Moving columns
msdial_rt <- msdial_rt[,c(7, 6, 1:2, 8, 3:5)]
##Exporting the MS-DIAL feature list with RI
write_xlsx(msdial_rt, "../X_hylaeae_metabolomics/Result/MSDIAL_feat_list_with_RI_EI.xlsx")



























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

