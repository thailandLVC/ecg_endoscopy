################################################################################
## Descriptive analysis of patient baseline characteristics
################################################################################


## -----------------------------------------------------------------------------
## SECTION 1: SETUP AND DATA PREPARATION
## -----------------------------------------------------------------------------

# Load required libraries
library(tidyverse)
library(comorbidity)
library(tableone)

# Load cohort data from script 1
cohort <- read.csv("~/cohort.csv")

# Subset unique observations per PID to analyse frequency of demographic characteristics
cohortUnq <- cohort %>%
  select(pid, is_expPatient, age_new, gender_new, rgn1_af) %>%
  # Categorize patients into age groups
  mutate(age_group = case_when(age_new >= 0 & age_new < 10 ~ '01-10',
                               age_new >= 10 & age_new < 20 ~ '10-20',
                               age_new >= 20 & age_new < 30 ~ '20-30',
                               age_new >= 30 & age_new < 40 ~ '30-40',
                               age_new >= 40 & age_new < 50 ~ '40-50',
                               age_new >= 50 & age_new < 60 ~ '50-60',
                               age_new >= 60 & age_new < 70 ~ '60-70',
                               age_new >= 70 & age_new < 80 ~ '70-80',
                               age_new >= 80 ~ '80+')) %>%
  distinct(pid, .keep_all = TRUE)

## -----------------------------------------------------------------------------
## SECTION 2: CALCULATION OF CHARLSON COMORBIDITY INDEX 
## -----------------------------------------------------------------------------

# Reshape dataframe for comorbidity package
dx <- cohort %>%
  filter(case_when(is_expPatient == 1 ~ postEKG == 0, is_expPatient == 0 ~ postEndo == 0)) %>%
  select(pid, pdx, sdx_adm) %>%
  unite("dx", c(pdx, sdx_adm), na.rm = TRUE, sep = " ") %>%
  mutate(dx = strsplit(dx, " ")) %>%
  unnest(dx) %>%
  distinct()

# Categorize comorbid conditions into Charlson domains using ICD-10 diagnoses
charlson <- comorbidity(x = dx, id = "pid", code = "dx", map = "charlson_icd10_quan", assign0 = FALSE) 

# Calculate the Charlson Comorbidity Index 
cci <- score(charlson, weights = NULL, assign0 = FALSE)

# Collect all scores in a new dataframe
comorbidities <- data.frame(pid = unique(dx$pid), cci = cci)

# Identify patients who have a risk factor for cardiovascular disease: hypertension (HBP), hypercholesterolemia (HCL), or diabetes (dm)
riskf <- dx %>%
  mutate(hbp = ifelse(str_detect(dx, "I10|I11|I12|I13|I14|I15|I16"),1,0),
         hcl = ifelse(str_detect(dx, "E780"),1,0),
         dm = ifelse(str_detect(dx, "E11"),1,0)) %>%
  select(-dx) %>%
  group_by(pid) %>%
  summarise(hbp = max(hbp), hcl = max(hcl), dm = max(dm))

# Count number of patients with each risk factor
riskf %>%
  summarise(across(hbp:dm, sum))

# Join comorbidities data and group comorbidity scores
cohortUnq <- cohortUnq %>%
  left_join(comorbidities, by = "pid") %>%
  mutate_at(c("cci","eci"), ~replace_na(.,0)) %>%
  # Group patients with CCI equal to 0, 1, 2, or 3+
  mutate(cci_group = case_when(cci < 3 ~ as.character(cci), cci >= 3 ~ "3+")) 

# Join risk factor data
cohortUnq <- cohortUnq %>%
  left_join(riskf, by = "pid") %>%
  mutate_at(c("hbp","hcl","dm"), ~replace_na(.,0))

#Subset exposed and comparison groups
cohortUnqE <- cohortUnq[cohortUnq$is_expPatient == 1,]
cohortUnqC <- cohortUnq[cohortUnq$is_expPatient == 0,]


## -----------------------------------------------------------------------------
## SECTION 3: ANALYSIS OF PATIENT CHARACTERISTICS 
## -----------------------------------------------------------------------------

# Convert categorical variables into factor
cohortUnq <- cohortUnq %>%
  mutate(gender_new = factor(gender_new, levels = c(1,2)),
         cci_group = factor(cci_group, levels = c("0","1","2","3+"), ordered = TRUE),
         rgn1_af = factor(rgn1_af, levels = seq(1,13,1)),
         hbp = factor(hbp, levels = c(0,1)),
         hcl = factor(hcl, levels = c(0,1)),
         dm = factor(dm, levels = c(0,1))) %>%       
  mutate(age_group = factor(age_group, levels = c("01-10","10-20","20-30","30-40",
                                                  "40-50","50-60","60-70","70-80","80+"), ordered = TRUE)) 

# Calculate means/standard deviations of continuous variables and counts/percentages of categorical variables
tab_base <- CreateTableOne(vars = c("age_new", "age_group", "gender_new", 
                                    "rgn1_af", "cci_group", "hbp", "hcl", "dm", "obs"), 
                           strata = "is_expPatient",
                           data = cohortUnq, test = FALSE)
print(tab_base, smd = TRUE, catDigits = 2)


# Significance tests
t.test(cohortUnq$is_expPatient, cohortUnq$age_new)
chisq.test(cohortUnq$is_expPatient, cohortUnq$gender_new, correct = FALSE)
chisq.test(cohortUnq$is_expPatient, cohortUnq$cci_group, correct = FALSE)
chisq.test(cohortUnq$is_expPatient, cohortUnq$rgn1_af, correct = FALSE)
chisq.test(cohortUnq$is_expPatient, cohortUnq$hbp, correct = FALSE)
chisq.test(cohortUnq$is_expPatient, cohortUnq$hcl, correct = FALSE)
chisq.test(cohortUnq$is_expPatient, cohortUnq$dm, correct = FALSE)
chisq.test(cohortUnq$is_expPatient, cohortUnq$obs, correct = FALSE)


## -----------------------------------------------------------------------------
## SECTION 4: AGE GROUP DISTRIBUTION
## -----------------------------------------------------------------------------

# Create age distribution data frame for males
ageDist.m <- cohortUnq %>%
  filter(gender_new == 1) %>%
  # Count number in each age group by cohort (exposed or comparison)
  group_by(is_expPatient, age_group) %>%
  summarise(n = n()) %>%
  # Calculate percentage. Set percentage as negative for comparison group to prep for age pyramid
  mutate(perc = case_when(is_expPatient == 1 ~ 100*n/sum(n), is_expPatient == 0 ~ -100*n/sum(n))) %>%
  # Convert age group and cohort to factors for plotting
  mutate(age_group = as_factor(age_group), is_expPatient = as_factor(is_expPatient)) %>%
  mutate(gender = "Male")

# Create age distribution data frame for females
ageDist.f <- cohortUnq %>%
  filter(gender_new == 2) %>%
  # Count number in each age group by cohort (exposed or comparison)
  group_by(is_expPatient, age_group) %>%
  summarise(n = n()) %>%
  # Calculate percentage. Set percentage as negative for comparison group to prep for age pyramid
  mutate(perc = case_when(is_expPatient == 1 ~ 100*n/sum(n), is_expPatient == 0 ~ -100*n/sum(n))) %>%
  # Convert age group and cohort to factors for plotting
  mutate(age_group = as_factor(age_group), is_expPatient = as_factor(is_expPatient)) %>%
  mutate(gender = "Female")

ageDist <- rbind(ageDist.m,ageDist.f) %>%
  mutate(shift = case_when(is_expPatient == 0 ~ -sd(perc)/3, TRUE ~ sd(perc)/3))


## -----------------------------------------------------------------------------
## SECTION 5: LENGTH OF STAY OF INDEX ENDOSCOPY VISIT
## -----------------------------------------------------------------------------

# Calculate summary statistics of length of stay of endoscopy visit
cohort %>%
  filter(dateadm_new == date_indexEndo) %>%
  group_by(is_expPatient) %>%
  summarise(mean = mean(los_new), sd = sd(los_new),
            median = median(los_new), iqr = IQR(los_new), 
            min = min(los_new), max = max(los_new))

# Calculate number and proportion of patients with each length of stay by receipt of ECG
los <- cohort %>%
  filter(dateadm_new == date_indexEndo & endoscopy == 1) %>%
  distinct(pid, .keep_all = TRUE) %>%
  group_by(is_expPatient, los_new) %>%
  summarise(n = n()) %>%
  mutate(perc = 100*n/sum(n))


## -----------------------------------------------------------------------------
## SECTION 6: REIMBURESEMENT AMOUNT OF ENDOSCOPY VISIT
## -----------------------------------------------------------------------------

# Calculate summary statistics of reimbursement amount in USD
cohort %>%
  filter(dateadm_new == date_indexEndo & endoscopy == 1) %>%
  distinct(pid, .keep_all = TRUE) %>%
  group_by(is_expPatient) %>%
  mutate(act_amt = coalesce(act_amt,0)) %>%
  summarise(mean = mean(act_amt_usd, na.rm = TRUE), sd = sd(act_amt_usd, na.rm = TRUE), 
            median = median(act_amt_usd, na.rm = TRUE), iqr = IQR(act_amt_usd, na.rm = TRUE), 
            min = min(act_amt_usd, na.rm = TRUE), max = max(act_amt_usd, na.rm = TRUE),
            perc0 = 100*sum(act_amt_usd == 0)/n())


## -----------------------------------------------------------------------------
## SECTION 7: SAME-DAY ECG
## -----------------------------------------------------------------------------

# Calculate the percentage of ECGs that occurred during the same admission as the endscopy
cohort %>%
  filter(dateadm_new == date_indexEndo & endoscopy == 1 & is_expPatient == 1) %>%
  distinct(pid, .keep_all = TRUE) %>%
  mutate(samedayEKG = case_when(date_indexEndo == date_preEKG ~ 1, TRUE ~ 0)) %>%
  group_by(samedayEKG) %>%
  summarise(n = n()) %>%
  mutate(perc = 100*n/sum(n))
