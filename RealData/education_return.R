# load libraries
library(tidyverse)

# read data
df <- readRDS('C:/Users/Zhenyu Wang/Dropbox/Zhenyu/Instrumental Variable/RealData-Econ/2-CFPS (Continuous_outcome_multiple_IV)/CFPS_education_return.rds')

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








