
# this file is to clean the CFPS database in order to generate a dataset with continuous Y and multiple natural IVs

library(tidyverse)
library(haven)
library(plyr)

prefix <- "/Users/junhui/Documents/Rutgers/data/"

# read data
cfps_famconf <- read_dta(paste0(prefix, "CFPS/cfps2018famconf_202008.dta"))
cfps_famecon <- read_dta(paste0(prefix, "CFPS/cfps2018famecon_202101.dta"))
cfps_person <- read_dta(paste0(prefix, "CFPS/cfps2018person_202012.dta"))


# if a variable has more than 30000 missing data, delete it
idx <- c()
for (i in 1:ncol(cfps_person)) {
  tmp <- sum(is.na(cfps_person[,i])) + sum(cfps_person[,i]<0, na.rm = T)
  if(tmp > 30000){
    idx <- c(idx, i)
  }
}
cfps_person_na_clean <- cfps_person[, -idx]
unlist_label_person <- lapply(cfps_person_na_clean, function(x) attr(x, "label"))
label1 <- ldply (unlist_label_person, data.frame)

idx <- c()
for (i in 1:ncol(cfps_famconf)) {
  tmp <- sum(is.na(cfps_famconf[,i])) +  sum(cfps_famconf[,i]<0, na.rm = T)
  if(tmp > 55000){
    idx <- c(idx, i)
  }
}
cfps_famconf_na_clean <- cfps_famconf[, -idx]
unlist_label_famconf <- lapply(cfps_famconf_na_clean, function(x) attr(x, "label"))
label2 <- ldply (unlist_label_famconf, data.frame)

idx <- c()
for (i in 1:ncol(cfps_famecon)) {
  tmp <- sum(is.na(cfps_famecon[,i])) + sum(cfps_famecon[,i]<0, na.rm = T)
  if(tmp > 10000){
    idx <- c(idx, i)
  }
}
cfps_famecon_na_clean <- cfps_famecon[, -idx]
unlist_label_famconf <- lapply(cfps_famecon_na_clean, function(x) attr(x, "label"))
label3 <- ldply (unlist_label_famconf, data.frame)

# select needed variables
cfps_person_na_clean <- cfps_person_na_clean %>% 
  select(pid, fid18, provcd18, countyid18, cid18, urban18, ibirthy, gender, fml_count, qa301, qc1,kw2y, 
         # x and y
         cfps2018edu, cfps2018eduy_im, egc1011, income, 
         # parents' variables as IV, individual's variables as IV or confounders representing ability
         qu701, qu301, wv102, wv105, qq1101, qq1102, qm510, qn4001,
         qbb002, qbb003, qbb004, qbb005, 
         whichwordlist, whichmathlist, wordtest18_sc2, mathtest18_sc2,
         # confounding covariates
         qga1, egc1052y, egc1053y, egc201, qg2, qka202) 

cfps_famconf_na_clean <- cfps_famconf_na_clean %>% 
  select(pid, pid_a_f, pid_a_m, tb4_a18_f, hukou_a18_f, tb4_a18_m, hukou_a18_m, tb4_a18_s, hukou_a18_s)

cfps_famecon_na_clean <- cfps_famecon_na_clean %>% 
  select(fid18, resp1pid, fp510)

# merge the datasets
cfps <- left_join(cfps_person_na_clean, cfps_famconf_na_clean, by = "pid")
cfps <- left_join(cfps, cfps_famecon_na_clean, by = "fid18")


# add parents' test results
# father
tmp_f <- cfps %>% select(pid, qu701, qu301, wv102, wv105, qq1101, qq1102, qm510, qn4001,
                         qbb002, qbb003, qbb004, qbb005, 
                         whichwordlist, whichmathlist, wordtest18_sc2, mathtest18_sc2) %>% 
  dplyr::rename(qu701_f = qu701, qu301_f = qu301, wv102_f = wv102, wv105_f = wv105, 
                qq1101_f = qq1101, qq1102_f = qq1102, qm510_f = qm510, qn4001_f = qn4001,
                qbb002_f = qbb002, qbb003_f = qbb003, qbb004_f = qbb004, qbb005_f = qbb005, 
                whichwordlist_f = whichwordlist, whichmathlist_f = whichmathlist, 
                wordtest18_sc2_f = wordtest18_sc2, mathtest18_sc2_f = mathtest18_sc2) 
cfps <- left_join(cfps, tmp_f, by = c("pid_a_f" = "pid"))

# mother
tmp_m <- cfps %>% select(pid, qu701, qu301, wv102, wv105, qq1101, qq1102, qm510, qn4001,
                         qbb002, qbb003, qbb004, qbb005, 
                         whichwordlist, whichmathlist, wordtest18_sc2, mathtest18_sc2) %>% 
  dplyr::rename(qu701_m = qu701, qu301_m = qu301, wv102_m = wv102, wv105_m = wv105, 
                qq1101_m = qq1101, qq1102_m = qq1102, qm510_m = qm510, qn4001_m = qn4001,
                qbb002_m = qbb002, qbb003_m = qbb003, qbb004_m = qbb004, qbb005_m = qbb005, 
                whichwordlist_m = whichwordlist, whichmathlist_m = whichmathlist, 
                wordtest18_sc2_m = wordtest18_sc2, mathtest18_sc2_m = mathtest18_sc2) 
cfps <- left_join(cfps, tmp_m, by = c("pid_a_m" = "pid"))


# clean the variables related to financial test
cfps <- cfps %>% mutate( # answer correct is 1, otherwise 0
  qbb002 = ifelse(qbb002 == 2, 1, 0),
  qbb002_m = ifelse(qbb002_m == 2, 1, 0),
  qbb002_f = ifelse(qbb002_f == 2, 1, 0),
  qbb003 = ifelse(qbb003 == 3, 1, 0),
  qbb003_m = ifelse(qbb003_m == 3, 1, 0),
  qbb003_f = ifelse(qbb003_f == 3, 1, 0),
  qbb004 = ifelse(qbb004 == 5, 1, 0),
  qbb004_m = ifelse(qbb004_m == 5, 1, 0),
  qbb004_f = ifelse(qbb004_f == 5, 1, 0))

# to show a list of column labels
unlist_label_df <- lapply(cfps, function(x) attr(x, "label"))
label <- ldply (unlist_label_df, data.frame)


saveRDS(cfps, paste0(prefix,"cleaned/2022/CFPS_education_return.rds"))
write_csv(label, paste0(prefix,"cleaned/2022/CFPS_education_return_codebook.csv"))