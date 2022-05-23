# load libraries
library(tidyverse)

# read data
df <- readRDS('CFPS_education_return.rds')

dim(df)
# [1] 37354    80

#### continuous Y: individual income --------------------------------------------
income <- df$income


#### endogenous variable X: the number of years educated ------------------------
education <- df$cfps2018eduy_im


#### IVs ------------------------------------------------------------------------
### From Literature
## 1. family background

# father's education
father_edu <- df$tb4_a18_f
# father's hukou
father_hukou <- df$hukou_a18_f
# mother's education
mother_edu <- df$tb4_a18_m
# mother's hukou
mother_hukou <- df$hukou_a18_m
# spouse's education
spouse_edu <- df$tb4_a18_s
# family size
fml_size <- df$fml_count


## 2. aggregate IV: the average years of schooling by age, gender, and region 
var <- c("ibirthy", "gender", "provcd18")
df <- df %>% dplyr::group_by(across(all_of(var))) %>% dplyr::mutate(group_edu = mean(cfps2018eduy_im, na.rm =T)) %>% ungroup
group_edu = df$group_edu


## 3. inter-generational ability endowments (the variables selected below into this category are based on my personal judgment)
# parents' frequency of using the internet to study
internet_study_f <- df$qu701_f
internet_study_m <- df$qu701_m

# parents' evaluation on the importance of study when using the internet
# Likert scale with 1 not important and 5 very important
internet_study_importance_f <- df$qu301_f
internet_study_importance_m <- df$qu301_m

# parents' answers on to what extent do you agree the statement
# Likert scale with 1 strongly disagree and 5 strongly agree
# fair competition leads to harmonious relationships
statement_fair_comp_f <- df$wv102_f
statement_fair_comp_m <- df$wv102_m
# talent and cleverness pay off
statement_talent_f <- df$wv105_f
statement_talent_m <- df$wv105_m
# hoping children have achievements
children_achievements_f <- df$qm510_f
children_achievements_m <- df$qm510_m

# parents' reading books or not
reading_f <- df$qq1101_f
reading_m <- df$qq1101_m

# the number of books parents have read
read_num_f <- df$qq1102_f
read_num_m <- df$qq1102_m

# parents' answer correctness on financial knowledge test 
saving_interest_f <- df$qbb002_f
saving_interest_m <- df$qbb002_m
inflation_f <- df$qbb003_f
inflation_m <- df$qbb003_m
stock_risk_f <- df$qbb004_f
stock_risk_m <- df$qbb004_m

# parents' willingness to use credit to meet their consumption needs
# Likert scale with 1 strongly disagree and 5 strongly agree
credit_consumption_f <- df$qbb005_f
credit_consumption_m <- df$qbb005_m

# parents' test scores on words and math (the number of answering correct, in total 34 questions)
wordtest_score_f <- df$wordtest18_sc2_f
wordtest_score_m <- df$wordtest18_sc2_m
mathtest_score_f <- df$mathtest18_sc2_f
mathtest_score_m <- df$mathtest18_sc2_m



### IVs we can select from based on our discretion
# individual's motivation for study
# frequency of using the internet to study
internet_study <- df$qu701

# the importance of study when using the internet
# Likert scale with 1 not important and 5 very important
internet_study_importance <- df$qu301

# individual's willingness to compete through ability/education
# to what extent do you agree the statement
# Likert scale with 1 strongly disagree and 5 strongly agree
# fair competition leads to harmonious relationships
statement_fair_comp <- df$wv102
# talent and cleverness pay off
statement_talent <- df$wv105

# if reading books
reading <- df$qq1101

# the number of books having read
read_num <- df$qq1102

# expected number of children
num_children <- df$qka202

# education expenditure in the last 12 months of one family member
edu_exp <- df$fp510





#### covariates -----------------------------------------------------------------
# urban or rural (0: rural, 1: urban)
urban <- df$urban18

# hukou (1: agriculture hukou, 3: non-agriculture hukou, 5: no hukou, 79: non Chinese citizen)
hukou <- df$qa301

# gender
gender <- df$gender

# full-time job expereince (0: don't have, 1: have)
full_time <- df$qga1

# the length of years working for the most recent job 
df <- df %>% mutate(egc1052y = ifelse(egc1052y < 0, NA, egc1052y),
                    egc1053y = ifelse(egc1053y < 0, NA, egc1053y)) %>% 
  # if the start year missing: NA; the start year is non-missing and the job continues: 2018 - the start year;    
  # the start year is non-missing and the job ends: the end year - the start year
  mutate(work_years = ifelse(is.na(egc1052y), NA, ifelse(egc1053y < 0, (2018 - egc1052y), (egc1053y - egc1052y))))
work_years = df$work_years

# the number of other jobs
num_other_jobs <- df$egc201

# the type of employer in the most recent job
employer_type <- df$qg2

# industry classification code
industry <- df$qg302code

## individual's ability
# answer correctness on financial knowledge test 
saving_interest <- df$qbb002
inflation <- df$qbb003
stock_risk <- df$qbb004

# willingness to use credit to meet their consumption needs
# Likert scale with 1 strongly disagree and 5 strongly agree
credit_consumption <- df$qbb005

# test scores on words and math (the number of answering correct, in total 34 questions)
wordtest_score <- df$wordtest18_sc2
mathtest_score <- df$mathtest18_sc2

##########################################################
source("src/main.R")
source("src/helpers.R")
source("src/TSHT-ldim.R")
source("src/invalidIV.R")
source("https://raw.githubusercontent.com/xlbristol/CIIV/main/R/CIIV_Functions.R")

library(MASS)
library(intervals)
Y_cand = income
D_cand = education
Z_cand = cbind(father_edu, father_hukou, mother_edu, mother_hukou, spouse_edu, fml_size,
               group_edu, internet_study_f, internet_study_m, internet_study_importance_f,
               internet_study_importance_m, statement_fair_comp_f, statement_fair_comp_m,
               statement_talent_f, statement_talent_m, children_achievements_f, children_achievements_m,
               reading_f, reading_m, read_num_f, read_num_m, saving_interest_f, saving_interest_m,
               inflation_f, inflation_m, stock_risk_f, stock_risk_m, credit_consumption_f, credit_consumption_m,
               wordtest_score_f, wordtest_score_m, mathtest_score_f, mathtest_score_m,
               internet_study, internet_study_importance, statement_fair_comp, statement_talent,
               reading, read_num, num_children, edu_exp)
X_cand = cbind(urban, hukou, gender, full_time, work_years, num_other_jobs, employer_type, industry,
               saving_interest, inflation, stock_risk, credit_consumption, wordtest_score, mathtest_score)

colSums(is.na(Z_cand))
## based on the number of missings, re-adjust Z_cand


# test on 4/26
set.seed(2022)
Z_cand1 = cbind(father_edu, mother_edu, spouse_edu, group_edu, fml_size,
                statement_fair_comp, statement_talent,reading,
                edu_exp)
X_cand1 = cbind(urban, hukou, gender)

DF = cbind(Y_cand, D_cand, Z_cand1, X_cand1)
colnames(DF) = c("Y","D","father_edu", "mother_edu","spouse_edu", "group_edu", "fml_size",
                 "statement_fair_comp", "statement_talent","reading",
                 "edu_exp","urban","hukou","gender")
DF = DF[!rowSums(is.na(DF)),]
dim(DF)
DF = DF[DF[,"Y"]>0,]
dim(DF)

apply(DF, MARGIN=2, FUN=function(X) sum(X>=0))
sum(apply(DF, MARGIN=1, FUN=function(X) all(X>=0)))

# select non-negative rows
temp1 = apply(DF, MARGIN=1, FUN=function(X) all(X>=0))
DF = DF[temp1,]
# select <79
DF = DF[DF[,"hukou"]<5,]
dim(DF)
DF[DF[,"hukou"] == 1,"hukou"]=0
DF[DF[,"hukou"] == 3,"hukou"]=1
# education gather
DF[DF[,"father_edu"] >= 5, "father_edu"] = 5
DF[DF[,"mother_edu"] >= 5, "mother_edu"] = 5
DF[DF[,"spouse_edu"] >= 5, "spouse_edu"] = 5
DF[,"edu_exp"] = log(1+DF[,"edu_exp"])
dim(DF)
Y = log(DF[,1]); D = DF[,2]; 
Z = DF[,c("father_edu", "mother_edu", "spouse_edu" ,"group_edu", "fml_size", 
          "statement_fair_comp", "statement_talent","reading",
          "edu_exp")]; 
X = DF[,c("urban","hukou","gender")]

CI.mat = matrix(NA, nrow=8, ncol=2)
rownames(CI.mat) = c("OLS", "TSLS", "TSHT", "CIIV", "Search", "Sample", "Union-1", "Union-2")

## OLS
CI.mat[1,] = confint(lm(Y~D+Z+X))[2,]
## TSLS
CI.mat[2,] = confint(ivreg(Y~D+X|Z+X))[2,]
## TSHT
out = TSHT(matrix(Y,ncol=1), matrix(D,ncol=1), Z, X)
CI.mat[3,] = out$ci
## CIIV
CIIV_out = CIIV(Y,D,Z,X,robust=TRUE, firststage=TRUE)
CIIV_out$Valid_Instruments
CI.mat[4,] = CIIV_out$ci_CIM
#CIIV_out$ci_CIM_GMM

## Searching
test1 = SearchingSampling(Y, D, Z, X, Sampling=FALSE); 
CI.mat[5,] = test1$CI
test1$VHat; test1$TSHT.out$SHat
## Sampling
test2 = SearchingSampling(Y, D, Z, X, Sampling=TRUE);
CI.mat[6,] = test2$CI

## Union-1
method=c("AR","CLR","TSLS","SarganTSLS","SarganCLR")
s_bar = ncol(Z)-1
resultOut.temp = tryCatch_E(invalidIVCI(Y, D, Z, X, U = s_bar, alpha = 0.05, alpha.overid = 0.05/2, intercept=TRUE,convexHull=FALSE,method=method),
                            ret.obj = list(NULL))
if (is.null(resultOut.temp$error)) {
  resultOut <- resultOut.temp$value
} else {
  next
}
CI.mat[7,] = resultOut$SarganTSLS
## Union-2
s_bar = ceiling(ncol(Z)/2)
resultOut.temp = tryCatch_E(invalidIVCI(Y, D, Z, X, U = s_bar, alpha = 0.05, alpha.overid = 0.05/2, intercept=TRUE,convexHull=FALSE,method=method),
                            ret.obj = list(NULL))
if (is.null(resultOut.temp$error)) {
  resultOut <- resultOut.temp$value
} else {
  next
}
resultOut
CI.mat[8,] = resultOut$SarganTSLS
CI.mat 
## note:
## that the sampling method has a different CI in paper, we pick the one in the file
## illustrate_sampling_RealData.R to comply with sampling figure.

## calculate concentration parameter ##
W = cbind(Z, X)
W = cbind(W, 1)
pz = ncol(Z)
qrW = qr(W)
ITT_D = qr.coef(qrW, D)[1:pz]
resid_D = as.vector(qr.resid(qrW, D))
SigmaSqD = mean(resid_D^2)

qrX = qr(cbind(X,1))
resid_Z = qr.resid(qrX, Z)
numerator1 = sum((resid_Z %*% ITT_D)^2)
param1 = numerator1 / SigmaSqD # 2906.36

Z_valid = Z[,test1$VHat]
X_valid = cbind(Z[,-test1$VHat], X, 1)
qrX_valid = qr(X_valid)
resid_Z_valid = qr.resid(qrX_valid, Z_valid)
numerator2 = sum((resid_Z_valid %*% ITT_D[test1$VHat])^2)
param2 = numerator2 / SigmaSqD # 2850.57, VHat is (1,2,3,4,8)

CIIV_VHat = which(colnames(Z) %in% CIIV_out$Valid_Instruments)
Z_valid_CIIV = Z[,CIIV_VHat]
X_valid_CIIV = cbind(Z[,-CIIV_VHat], X, 1)
qrX_valid_CIIV = qr(X_valid_CIIV)
resid_Z_valid_CIIV = qr.resid(qrX_valid_CIIV, Z_valid_CIIV)
numerator3 = sum((resid_Z_valid_CIIV %*% ITT_D[CIIV_VHat])^2)
param3 = numerator3 / SigmaSqD # 2850.751 VHat is 




################# Previously Tested, Do not use anymore ##############
if(case==1){
  Z_cand1 = cbind(father_edu, mother_edu, spouse_edu, fml_size, group_edu,
                  internet_study, statement_fair_comp, statement_talent,
                  reading, num_children, edu_exp)
  X_cand1 = cbind(urban, hukou, gender)
  
  
  DF = cbind(Y_cand, D_cand, Z_cand1, X_cand1)
  DF = DF[!rowSums(is.na(DF)),]
  dim(DF)
  DF = DF[DF[,1]>0,]
  dim(DF)
  
  apply(DF, MARGIN=2, FUN=function(X) sum(X>=0))
  sum(apply(DF, MARGIN=1, FUN=function(X) all(X>=0))) # 2956
  
  # select non-negative rows
  temp1 = apply(DF, MARGIN=1, FUN=function(X) all(X>=0))
  DF = DF[temp1,]
  # select <79
  temp2 = apply(DF[,c(8:12, 15)], MARGIN=1, FUN=function(X) all(X<79))
  DF = DF[temp2,]
  dim(DF)
  
  Y = log(DF[,1]); D = DF[,2]; Z = DF[,3:13]; X = DF[,14:16]
  
  test1 = SearchingSampling(Y, D, Z, X, Sampling=FALSE); 
  test1$CI
  test1$VHat; test1$TSHT.out$SHat
  test2 = SearchingSampling(Y, D, Z, X, Sampling=TRUE); 
  test2$CI
}
if(case==2){
  Z_cand1 = cbind(father_edu, mother_edu, spouse_edu, fml_size, group_edu,
                  internet_study, statement_fair_comp, statement_talent,
                  reading, num_children, edu_exp)
  X_cand1 = cbind(urban, hukou, gender)
  
  
  DF = cbind(Y_cand, D_cand, Z_cand1, X_cand1)
  DF = DF[!rowSums(is.na(DF)),]
  dim(DF)
  DF = DF[DF[,1]>0,]
  DF = DF[DF[,13]>0,]
  dim(DF)
  
  apply(DF, MARGIN=2, FUN=function(X) sum(X>=0))
  sum(apply(DF, MARGIN=1, FUN=function(X) all(X>=0))) # 2956
  
  # select non-negative rows
  temp1 = apply(DF, MARGIN=1, FUN=function(X) all(X>=0))
  DF = DF[temp1,]
  # select <79
  temp2 = apply(DF[,c(8:12, 15)], MARGIN=1, FUN=function(X) all(X<79))
  DF = DF[temp2,]
  dim(DF)
  
  DF[,13] = log(DF[,13])
  Y = log(DF[,1]); D = DF[,2]; Z = DF[,3:13]; X = DF[,14:16]
  
  test1 = SearchingSampling(Y, D, Z, X, Sampling=FALSE); 
  test1$CI
  test1$VHat; test1$TSHT.out$SHat
  test2 = SearchingSampling(Y, D, Z, X, Sampling=TRUE); 
  test2$CI
}
if(case==3){
  Z_cand1 = cbind(father_edu, mother_edu, spouse_edu, fml_size, group_edu,
                  statement_fair_comp, statement_talent,
                  reading, num_children, edu_exp)
  X_cand1 = cbind(urban, hukou, gender)
  
  
  DF = cbind(Y_cand, D_cand, Z_cand1, X_cand1)
  DF = DF[!rowSums(is.na(DF)),]
  dim(DF)
  DF = DF[DF[,1]>0,]
  DF = DF[DF[,12]>0,]
  dim(DF)
  
  apply(DF, MARGIN=2, FUN=function(X) sum(X>=0))
  sum(apply(DF, MARGIN=1, FUN=function(X) all(X>=0))) # 2956
  
  # select non-negative rows
  temp1 = apply(DF, MARGIN=1, FUN=function(X) all(X>=0))
  DF = DF[temp1,]
  # select <79
  temp2 = apply(DF[,c(8:11, 14)], MARGIN=1, FUN=function(X) all(X<79))
  DF = DF[temp2,]
  dim(DF)
  
  DF[,12] = log(DF[,12])
  Y = log(DF[,1]); D = DF[,2]; Z = DF[,3:12]; X = DF[,13:15]
  
  test1 = SearchingSampling(Y, D, Z, X, Sampling=FALSE); 
  test1$CI
  test1$VHat; test1$TSHT.out$SHat
  test2 = SearchingSampling(Y, D, Z, X, Sampling=TRUE); 
  test2$CI
}
if(case==4){
  Z_cand1 = cbind(father_edu, mother_edu, spouse_edu, fml_size, group_edu,
                  statement_fair_comp, statement_talent,
                  reading, edu_exp)
  X_cand1 = cbind(urban, hukou, gender)
  
  
  DF = cbind(Y_cand, D_cand, Z_cand1, X_cand1)
  DF = DF[!rowSums(is.na(DF)),]
  dim(DF)
  DF = DF[DF[,1]>0,]
  DF = DF[DF[,11]>0,]
  dim(DF)
  
  apply(DF, MARGIN=2, FUN=function(X) sum(X>=0))
  sum(apply(DF, MARGIN=1, FUN=function(X) all(X>=0))) # 2956
  
  # select non-negative rows
  temp1 = apply(DF, MARGIN=1, FUN=function(X) all(X>=0))
  DF = DF[temp1,]
  # select <79
  temp2 = apply(DF[,c(8:10, 14)], MARGIN=1, FUN=function(X) all(X<79))
  DF = DF[temp2,]
  dim(DF)
  
  DF[,11] = log(DF[,11])
  Y = log(DF[,1]); D = DF[,2]; Z = DF[,3:11]; X = DF[,12:14]
  
  test1 = SearchingSampling(Y, D, Z, X, Sampling=FALSE); 
  test1$CI
  test1$VHat; test1$TSHT.out$SHat
  test2 = SearchingSampling(Y, D, Z, X, Sampling=TRUE); 
  test2$CI
}
if(case==5){
  Z_cand1 = cbind(father_edu, mother_edu, spouse_edu, num_children, group_edu,
                  statement_fair_comp, statement_talent,
                  reading, edu_exp)
  X_cand1 = cbind(urban, hukou, gender)
  
  
  DF = cbind(Y_cand, D_cand, Z_cand1, X_cand1)
  DF = DF[!rowSums(is.na(DF)),]
  dim(DF)
  DF = DF[DF[,1]>0,]
  DF = DF[DF[,11]>0,]
  dim(DF)
  
  apply(DF, MARGIN=2, FUN=function(X) sum(X>=0))
  sum(apply(DF, MARGIN=1, FUN=function(X) all(X>=0))) # 2956
  
  # select non-negative rows
  temp1 = apply(DF, MARGIN=1, FUN=function(X) all(X>=0))
  DF = DF[temp1,]
  # select <79
  temp2 = apply(DF[,c(8:10, 14)], MARGIN=1, FUN=function(X) all(X<79))
  DF = DF[temp2,]
  dim(DF)
  
  DF[,11] = log(DF[,11])
  Y = log(DF[,1]); D = DF[,2]; Z = DF[,3:11]; X = DF[,12:14]
  
  test1 = SearchingSampling(Y, D, Z, X, Sampling=FALSE); 
  test1$CI
  test1$VHat; test1$TSHT.out$SHat
  test2 = SearchingSampling(Y, D, Z, X, Sampling=TRUE); 
  test2$CI
}
if(case==6){
  
  set.seed(0)
  Z_cand1 = cbind(father_edu, mother_edu, spouse_edu, group_edu,
                  statement_fair_comp, statement_talent,
                  reading, edu_exp)
  X_cand1 = cbind(urban, hukou, gender)
  
  DF = cbind(Y_cand, D_cand, Z_cand1, X_cand1)
  DF = DF[!rowSums(is.na(DF)),]
  dim(DF)
  DF = DF[DF[,1]>0,]
  #DF = DF[DF[,10]>0,]
  dim(DF)
  
  apply(DF, MARGIN=2, FUN=function(X) sum(X>=0))
  sum(apply(DF, MARGIN=1, FUN=function(X) all(X>=0))) # 3775
  
  # select non-negative rows
  temp1 = apply(DF, MARGIN=1, FUN=function(X) all(X>=0))
  DF = DF[temp1,]
  # select <79
  temp2 = apply(DF[,c(7:9, 12)], MARGIN=1, FUN=function(X) all(X<79))
  DF = DF[temp2,]
  dim(DF)
  # select < 5
  temp3 = DF[,12] < 5
  DF = DF[temp3,]
  DF[DF[,12] == 1,12]=0
  DF[DF[,12] == 3,12]=1
  # education gather 5,6,7
  DF[DF[,3] >= 5, 3] = 5
  DF[DF[,4] >= 5, 4] = 5
  DF[DF[,5] >= 5, 5] = 5
  DF[,10] = log(1+DF[,10])
  dim(DF) # 3771
  Y = log(DF[,1]); D = DF[,2]; Z = DF[,3:10]; X = DF[,11:13]
  
  CI.mat = matrix(NA, nrow=8, ncol=2)
  rownames(CI.mat) = c("OLS", "TSLS", "TSHT", "CIIV", "Search", "Sample", "Union-1", "Union-2")
  
  ## OLS
  CI.mat[1,] = confint(lm(Y~D+Z+X))[2,]
  ## TSLS
  CI.mat[2,] = confint(ivreg(Y~D+X|Z+X))[2,]
  ## TSHT
  out = TSHT(matrix(Y,ncol=1), matrix(D,ncol=1), Z, X)
  CI.mat[3,] = out$ci
  ## CIIV
  CIIV_out = CIIV(Y,D,Z,X,robust=TRUE)
  CIIV_out$Valid_Instruments
  CI.mat[4,] = CIIV_out$ci_CIM
  #CIIV_out$ci_CIM_GMM
  
  ## Searching
  test1 = SearchingSampling(Y, D, Z, X, Sampling=FALSE); 
  CI.mat[5,] = test1$CI
  test1$VHat; test1$TSHT.out$SHat
  ## Sampling
  test2 = SearchingSampling(Y, D, Z, X, Sampling=TRUE);
  CI.mat[6,] = test2$CI
  
  ## Union-1
  method=c("AR","CLR","TSLS","SarganTSLS","SarganCLR")
  s_bar = ncol(Z)-1
  resultOut.temp = tryCatch_E(invalidIVCI(Y, D, Z, X, U = s_bar, alpha = 0.05, alpha.overid = 0.05/2, intercept=TRUE,convexHull=FALSE,method=method),
                              ret.obj = list(NULL))
  if (is.null(resultOut.temp$error)) {
    resultOut <- resultOut.temp$value
  } else {
    next
  }
  CI.mat[7,] = resultOut$SarganTSLS
  ## Union-2
  s_bar = ceiling(ncol(Z)/2)
  resultOut.temp = tryCatch_E(invalidIVCI(Y, D, Z, X, U = s_bar, alpha = 0.05, alpha.overid = 0.05/2, intercept=TRUE,convexHull=FALSE,method=method),
                              ret.obj = list(NULL))
  if (is.null(resultOut.temp$error)) {
    resultOut <- resultOut.temp$value
  } else {
    next
  }
  resultOut
  CI.mat[8,] = resultOut$SarganTSLS
  
}
if(case==6){
  
  set.seed(0)
  Z_cand1 = cbind(father_edu, mother_edu, spouse_edu, group_edu,
                  statement_fair_comp, statement_talent,
                  reading, edu_exp, num_children)
  X_cand1 = cbind(urban, hukou, gender)
  
  DF = cbind(Y_cand, D_cand, Z_cand1, X_cand1)
  DF = DF[!rowSums(is.na(DF)),]
  dim(DF)
  DF = DF[DF[,1]>0,]
  #DF = DF[DF[,10]>0,]
  dim(DF)
  
  apply(DF, MARGIN=2, FUN=function(X) sum(X>=0))
  sum(apply(DF, MARGIN=1, FUN=function(X) all(X>=0))) # 3764
  
  # select non-negative rows
  temp1 = apply(DF, MARGIN=1, FUN=function(X) all(X>=0))
  DF = DF[temp1,]
  # select <79
  temp2 = apply(DF[,c(7:9, 13)], MARGIN=1, FUN=function(X) all(X<79))
  DF = DF[temp2,]
  dim(DF)
  # select < 5
  temp3 = DF[,13] < 5
  DF = DF[temp3,]
  DF[DF[,13] == 1,13]=0
  DF[DF[,13] == 3,13]=1
  # education gather 5,6,7
  DF[DF[,3] >= 5, 3] = 5
  DF[DF[,4] >= 5, 4] = 5
  DF[DF[,5] >= 5, 5] = 5
  DF[,10] = log(1+DF[,10])
  dim(DF) # 3760
  Y = log(DF[,1]); D = DF[,2]; Z = DF[,3:11]; X = DF[,12:14]
  
  CI.mat = matrix(NA, nrow=8, ncol=2)
  rownames(CI.mat) = c("OLS", "TSLS", "TSHT", "CIIV", "Search", "Sample", "Union-1", "Union-2")
  
  ## OLS
  CI.mat[1,] = confint(lm(Y~D+Z+X))[2,]
  ## TSLS
  CI.mat[2,] = confint(ivreg(Y~D+X|Z+X))[2,]
  ## TSHT
  out = TSHT(matrix(Y,ncol=1), matrix(D,ncol=1), Z, X)
  CI.mat[3,] = out$ci
  ## CIIV
  CIIV_out = CIIV(Y,D,Z,X,robust=TRUE)
  CIIV_out$Valid_Instruments
  CI.mat[4,] = CIIV_out$ci_CIM
  #CIIV_out$ci_CIM_GMM
  
  ## Searching
  test1 = SearchingSampling(Y, D, Z, X, Sampling=FALSE); 
  CI.mat[5,] = test1$CI
  test1$VHat; test1$TSHT.out$SHat
  ## Sampling
  test2 = SearchingSampling(Y, D, Z, X, Sampling=TRUE);
  CI.mat[6,] = test2$CI
  
  ## Union-1
  method=c("AR","CLR","TSLS","SarganTSLS","SarganCLR")
  s_bar = ncol(Z)-1
  resultOut.temp = tryCatch_E(invalidIVCI(Y, D, Z, X, U = s_bar, alpha = 0.05, alpha.overid = 0.05/2, intercept=TRUE,convexHull=FALSE,method=method),
                              ret.obj = list(NULL))
  if (is.null(resultOut.temp$error)) {
    resultOut <- resultOut.temp$value
  } else {
    next
  }
  CI.mat[7,] = resultOut$SarganTSLS
  ## Union-2
  s_bar = ceiling(ncol(Z)/2)
  resultOut.temp = tryCatch_E(invalidIVCI(Y, D, Z, X, U = s_bar, alpha = 0.05, alpha.overid = 0.05/2, intercept=TRUE,convexHull=FALSE,method=method),
                              ret.obj = list(NULL))
  if (is.null(resultOut.temp$error)) {
    resultOut <- resultOut.temp$value
  } else {
    next
  }
  resultOut
  CI.mat[8,] = resultOut$SarganTSLS
  
}
if(case==6){
  
  set.seed(0)
  Z_cand1 = cbind(father_edu, mother_edu, spouse_edu, group_edu,
                  statement_fair_comp, reading,
                  num_children, edu_exp)
  X_cand1 = cbind(urban, hukou, gender)
  
  DF = cbind(Y_cand, D_cand, Z_cand1, X_cand1)
  DF = DF[!rowSums(is.na(DF)),]
  dim(DF)
  DF = DF[DF[,1]>0,]
  #DF = DF[DF[,10]>0,]
  dim(DF)
  
  apply(DF, MARGIN=2, FUN=function(X) sum(X>=0))
  sum(apply(DF, MARGIN=1, FUN=function(X) all(X>=0))) # 3775
  
  # select non-negative rows
  temp1 = apply(DF, MARGIN=1, FUN=function(X) all(X>=0))
  DF = DF[temp1,]
  # select <79
  temp2 = apply(DF[,c(7:9, 12)], MARGIN=1, FUN=function(X) all(X<79))
  DF = DF[temp2,]
  dim(DF)
  # select < 5
  temp3 = DF[,12] < 5
  DF = DF[temp3,]
  DF[DF[,12] == 1,12]=0
  DF[DF[,12] == 3,12]=1
  # education gather 5,6,7
  DF[DF[,3] >= 5, 3] = 5
  DF[DF[,4] >= 5, 4] = 5
  DF[DF[,5] >= 5, 5] = 5
  DF[,10] = log(1+DF[,10])
  dim(DF) # 3771
  Y = log(DF[,1]); D = DF[,2]; Z = DF[,3:10]; X = DF[,11:13]
  
  CI.mat = matrix(NA, nrow=8, ncol=2)
  rownames(CI.mat) = c("OLS", "TSLS", "TSHT", "CIIV", "Search", "Sample", "Union-1", "Union-2")
  
  ## OLS
  CI.mat[1,] = confint(lm(Y~D+Z+X))[2,]
  ## TSLS
  CI.mat[2,] = confint(ivreg(Y~D+X|Z+X))[2,]
  ## TSHT
  out = TSHT(matrix(Y,ncol=1), matrix(D,ncol=1), Z, X)
  CI.mat[3,] = out$ci
  ## CIIV
  CIIV_out = CIIV(Y,D,Z,X,robust=TRUE)
  CIIV_out$Valid_Instruments
  CI.mat[4,] = CIIV_out$ci_CIM
  #CIIV_out$ci_CIM_GMM
  
  ## Searching
  test1 = SearchingSampling(Y, D, Z, X, Sampling=FALSE); 
  CI.mat[5,] = test1$CI
  test1$VHat; test1$TSHT.out$SHat
  ## Sampling
  test2 = SearchingSampling(Y, D, Z, X, Sampling=TRUE);
  CI.mat[6,] = test2$CI
  
  ## Union-1
  method=c("AR","CLR","TSLS","SarganTSLS","SarganCLR")
  s_bar = ncol(Z)-1
  resultOut.temp = tryCatch_E(invalidIVCI(Y, D, Z, X, U = s_bar, alpha = 0.05, alpha.overid = 0.05/2, intercept=TRUE,convexHull=FALSE,method=method),
                              ret.obj = list(NULL))
  if (is.null(resultOut.temp$error)) {
    resultOut <- resultOut.temp$value
  } else {
    next
  }
  CI.mat[7,] = resultOut$SarganTSLS
  ## Union-2
  s_bar = ceiling(ncol(Z)/2)
  resultOut.temp = tryCatch_E(invalidIVCI(Y, D, Z, X, U = s_bar, alpha = 0.05, alpha.overid = 0.05/2, intercept=TRUE,convexHull=FALSE,method=method),
                              ret.obj = list(NULL))
  if (is.null(resultOut.temp$error)) {
    resultOut <- resultOut.temp$value
  } else {
    next
  }
  resultOut
  CI.mat[8,] = resultOut$SarganTSLS
  
}