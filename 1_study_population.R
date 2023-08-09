################################################################################
## Selection of study population and analysis of the prevalence of low-risk
## endoscopy and low-value ECGs
################################################################################

## -----------------------------------------------------------------------------
## SECTION 1: SETUP AND PREPARATION OF MASTER DATASET  
## -----------------------------------------------------------------------------

# Load required packages
library(tidyverse)

# Import data of all admissions of patients who have received a GI endoscopy, chest x-ray, and/or a electrocardiogram
allAdm <- read.csv("~/gi_endoscopy_cxr_ECG_icd9.csv")
allAdm[allAdm == ""] <- NA

# Count number of unique PIDs
n_distinct(allAdm$pid)

# Collate all procedures and all diagnoses for each admission
# Identify GI endoscopy recipients. Classify as therapeutic or diagnostic endoscopy.
allAdm <- allAdm %>% 
  group_by(pid) %>%
  arrange(dateadm_new, .by_group = TRUE) %>%
  # Clean procedure variables
  mutate_at(vars(starts_with('proc')), ~gsub("\\+.*","",.)) %>%
  # Merge all procedure columns
  unite("proc_adm", c(proc1:proc21), na.rm = TRUE, sep = " ") %>%
  # Merge all secondary diagnosis columns
  unite("sdx_adm", c(sdx1:sdx21), na.rm = TRUE, sep = " ") %>% 
  # Create an indicator variable for any endoscopic procedure
  mutate(endoscopy = ifelse(str_detect(proc_adm,"4221|4222|4223|4224|4231|4232|4233|4281|
                                       4291|4292|4311|4341|4342|4349|4411|4412|4413|4414|
                                       4422|4440|4441|4442|4443|4444|4491|4493|4494|4511|
                                       4512|4513|4514|4516|4521|4522|4523|4524|4525|4530|
                                       4531|4532|4533|4534|4541|4542|4543|4549|4685|4686|
                                       4687|4821|4822|4823|4824|4831|4832|4833|4834|4835|
                                       4836|4921|4922|4923|4931|4939|5110|5111|5164|5181|
                                       5182|5183|5184|5185|5186|5187|5188|5198|5199|5213|
                                       5214|5221|5222|5292|5293|5294|5297|5298|9622|9623|
                                       9802|9803|9804|9805"),1,0), .after = proc_adm) 


## -----------------------------------------------------------------------------
## SECTION 2: IDENTIFICATION OF GI ENDOSCOPY PATIENTS  
## -----------------------------------------------------------------------------

# Subset all admissions (i.e., complete hospitalisation history) of GI endoscopy recipients 
endo <- allAdm %>% 
  # Keep variables of interest
  select(c(1:20)) %>% 
  # Keep all rows of patients who received at least 1 endoscopy
  group_by(pid) %>%
  arrange(dateadm_new, .by_group = TRUE) %>%
  filter(any(endoscopy == 1)) %>%
  # Remove any patients with admissions before 2016 (told to be unreliable) 
  filter(!any(yearAdm < 2016)) %>%
  # Format admission and discharge dates as dates
  mutate(dateadm_new = as.POSIXct(dateadm_new, format = "%d%b%Y"), datedsc = as.POSIXct(datedsc_new, format = "%d%b%Y")) %>%
  mutate(yr_month = substr(dateadm_new,1,7), .after = dateadm_new) %>%
  # Determine fiscal year
  mutate(fiscalyr = ifelse(month(dateadm_new)>9, year(dateadm_new)+1, year(dateadm_new)), .after = yearAdm)

# Count number of GI endoscopy admissions
sum(endo$endoscopy)

# Count number of unique GI endoscopy recipients
n_distinct(endo$pid)


## -----------------------------------------------------------------------------
## SECTION 3: IDENTIFICATION OF ADMISSIONS THAT INCLUDE ECG, OTHER CARDIAC TEST,
## AND/OR CARDIOVASCULAR DIAGNOSIS
## -----------------------------------------------------------------------------

# Create a new index for PID for quicker processing time
pid_index <- data.frame(pid = unique(endo$pid)) %>%
  mutate(pid_new = row_number())

# Identify admissions that include ECG and cardiovascular diagnosis (an exclusion criteria)
endo <- endo %>%
  # Recode transaction ID and join new PID
  ungroup() %>%
  mutate(tran_id_new = row_number(), .after = tran_id) %>%
  left_join(pid_index, by = "pid") %>%
  relocate(pid_new, .after = pid)  %>%
  select(-c(pid, tran_id)) %>%
  rename(pid = pid_new, tran_id = tran_id_new) %>%
  # Create an indicator variable for ECG
  mutate(ECG = ifelse(str_detect(proc_adm, "8951|8952"),1,0), .after = endoscopy) %>%
  # Create an indicator variable for cardiovascular diagnosis based on primary and secondary diagnoses
  mutate(cv_dx = ifelse(
    str_detect(pdx, "I05|I06|I07|I08|I09|I11|I13|I20|I21|I22|I23|I24|I25|I33|I34|I35|I36|I37|I38|
               I39|I44|I45|I47|I48|I49|I50|R071|R072|R078|R079|R060|R001|I95|R224|R000|R002|R42|
               R55|I469|Z950") |
    str_detect(sdx_adm,"I05|I06|I07|I08|I09|I11|I13|I20|I21|I22|I23|I24|I25|I33|I34|I35|I36|I37|
               I38|I39|I44|I45|I47|I48|I49|I50|R071|R072|R078|R079|R060|R001|I95|R224|R000|R002|
               R42|R55|I469|Z950"),1,0), .after = sdx_adm) %>%
  # Create an indicator variable for other cardiac tests
  mutate(cv_pretest = ifelse(str_detect(proc_adm, 
                                        "8941|8942|8943|8944|3728"),1,0), .after = ECG) 


## -----------------------------------------------------------------------------
## SECTION 4: IDENTIFICATION OF LOW-RISK ENDOSCOPY RECIPIENTS
## -----------------------------------------------------------------------------

# Identify patients who underwent diagnostic gastroscopy and colonoscopy (endoscopy with low cardiac risk) 
endo <- endo %>%
  # Create an indicator variable for diagnostic gastroscopy and colonoscopy
  mutate(endoscopy_dx_gc = ifelse(str_detect(proc_adm, "4311|4411|4412|4413|4414|4521|4522|4523|
                                             4524|4525"),1,0), .after = endoscopy_dx)


## -----------------------------------------------------------------------------
## SECTION 5: EXTRACTION OF DATES OF INDEX ENDOSCOPY, INDEX ECG, AND ANY PRIOR
## CARDIOVASCULAR DISEASE DIAGNOSES
## -----------------------------------------------------------------------------

# Extract all GI endoscopy dates 
datesEndo <- endo %>% 
  mutate(date_endo = case_when(endoscopy == 1 ~ dateadm_new)) %>%
  select(pid, date_endo, fiscalyr, endoscopy_dx_gc) %>%
  filter(!is.na(date_endo))

# Extract all ECG dates
datesECG <- endo %>% 
  mutate(date_ECG = case_when(ECG == 1 ~ dateadm_new)) %>%
  select(pid, date_ECG) %>%
  filter(!is.na(date_ECG))

# Extract all dates of other cardiac tests
datesOtherPre <- endo %>% 
  mutate(date_otherPre = case_when(cv_pretest == 1 ~ dateadm_new)) %>%
  select(pid, date_otherPre) %>%
  filter(!is.na(date_otherPre))

# Extract all dates of cardiovascular disease diagnoses
datesCv <- endo %>% 
  mutate(date_cv = case_when(cv_dx == 1 ~ dateadm_new)) %>%
  group_by(pid) %>%
  mutate(first_cv = min(date_cv)) %>%
  filter(!is.na(first_cv)) %>%
  select(pid, first_cv) %>%
  distinct()

# Pair GI dates with ECG and other preprocedural cardiac test dates such that all combinations are kept
dates <- merge(datesEndo, datesECG, by = 'pid', all = TRUE)
dates <- merge(dates, datesOtherPre, by = 'pid', all = TRUE)

# Join cardiovascular diagnosis dates
dates <- left_join(dates, datesCv, by = 'pid')


## -----------------------------------------------------------------------------
## SECTION 6: IDENTIFICATION OF ECG AND COMPARISON GROUPS USING DATE DATA
## -----------------------------------------------------------------------------

# Identify patients who received a ECG within 90 days from GI endoscopy
# Identify patients who received a cardiovascular diagnosis before ECG
dates <- dates %>%
  group_by(pid) %>%
  # Calculate days between ECG and endoscopy
  mutate(timeDiff_ECG_endo = as.numeric(date_endo - date_ECG)) %>%
  # Create an indicator variable that equals 1 when the ECG is pre-procedural (within 90 days of endoscopy)
  mutate(is_preECG = case_when(between(timeDiff_ECG_endo,0,90) ~ 1, TRUE ~ 0)) %>%
  # Create an indicator variable that equals 1 when the patient has received at least 1 pre-procedural ECG
  mutate(is_preECG_pid = max(is_preECG), .after = is_preECG) %>%
  
  # Calculate days between ECG and endoscopy
  mutate(timeDiff_otherPre_endo = as.numeric(date_endo - date_otherPre)) %>%
  # Create an indicator variable that equals 1 when the other cardiac test is pre-procedural (within 90 days of endoscopy)
  mutate(is_otherPre = case_when(between(timeDiff_otherPre_endo,0,90) ~ 1, TRUE ~ 0)) %>%
  
  # Create an indicator variable for a cardiovascular diagnosis before ECG
  mutate(is_cv = case_when(date_ECG > first_cv ~ 1, TRUE ~ 0)) %>% 
  # Create an indicator variable for patients that fall within the study period
  mutate(study_period = case_when(fiscalyr == 2018 | fiscalyr == 2019 ~ 1, TRUE ~ 0))

# Subset ECG patients excluded due to cardiovascular diagnosis
ECG_cv_pid <- dates %>%
  group_by(pid) %>%
  filter(is_preECG == 1 & is_cv == 1 & study_period == 1) %>%
  distinct(pid) %>%
  pull()

# Subset ECG patients excluded due to endoscopic procedure with moderate-to-high cardiac risk
ECG_procRisk_pid <- dates %>%
  group_by(pid) %>%
  filter(!pid %in% ECG_cv_pid) %>%
  filter(is_preECG == 1 & is_cv == 0 & study_period == 1) %>%
  filter(!any(endoscopy_dx_gc == 1)) %>%
  distinct(pid) %>%
  pull()

# Subset index endoscopy dates of exposed group
indexE <- dates %>%
  # Keep unique patients with pre-procedural ECG before GI endoscopy and no CV diagnosis
  filter(!pid %in% c(ECG_cv_pid,ECG_procRisk_pid)) %>%
  filter(is_preECG == 1 & is_cv == 0 & endoscopy_dx_gc == 1 & study_period == 1) %>%
  group_by(pid) %>%
  # Keep unique patients
  mutate(count = row_number(pid)) %>%
  filter(count == 1) %>%
  rename("date_indexEndo" = "date_endo", "date_preECG" = "date_ECG") %>%
  mutate(is_expPatient = 1) %>%
  select(pid, date_indexEndo, date_preECG, is_expPatient)

# Extract PIDs 
pidE <- indexE %>%
  distinct(pid) %>%
  pull()

# Count PIDs in exposed group
n_ECG <- length(pidE)

# Subset non-ECG patients excluded due to cardiovascular diagnosis
noECG_cv_pid <- dates %>%
  filter(!pid %in% c(pidE,ECG_cv_pid,ECG_procRisk_pid)) %>%
  filter(study_period == 1 & is_cv == 1) %>%
  distinct(pid) %>%
  pull()

# Subset non-ECG patients excluded due to receipt of another cardiac test 
noECG_otherPre_pid <- dates %>%
  filter(!pid %in% c(pidE,ECG_cv_pid,ECG_procRisk_pid,noECG_cv_pid)) %>%
  filter(study_period == 1 & is_otherPre == 1) %>%
  distinct(pid) %>%
  pull()

# Subset non-ECG patients excluded due to endoscopic procedure with moderate-to-high cardiac risk
noECG_procRisk_pid <- dates %>%
  group_by(pid) %>%
  filter(!pid %in% c(pidE,ECG_cv_pid,ECG_procRisk_pid,noECG_cv_pid,noECG_otherPre_pid)) %>%
  filter(is_cv == 0 & study_period == 1) %>%
  filter(!any(endoscopy_dx_gc == 1)) %>%
  distinct(pid) %>%
  pull()

# Subset index endoscopy dates of comparison group
indexC <- dates %>%
  filter(!pid %in% c(pidE,ECG_cv_pid,ECG_procRisk_pid,noECG_cv_pid,noECG_otherPre_pid,noECG_procRisk_pid)) %>%
  filter(study_period == 1 & endoscopy_dx_gc == 1) %>%
  group_by(pid) %>%
  # Keep unique patients
  mutate(count = row_number(pid)) %>%
  filter(count == 1) %>%
  rename(date_indexEndo = date_endo) %>%
  mutate(is_expPatient = 0) %>%
  select(pid, date_indexEndo, is_expPatient)

# Extract PIDs 
pidC <- indexC %>%
  distinct(pid) %>%
  pull()

# Join index endoscopy dates of exposed and comparison groups
dates_indexEndo <- rbind(indexE, indexC)

# Count annual number of GI endoscopy admissions, admissions of patients without prior cardiovascular diagnosis, and admissions with pre-procedural ECG
endoAdm_yr <- dates %>%
  filter(endoscopy_dx_gc == 1 & is_cv == 0) %>%
  group_by(fiscalyr) %>%
  tally() %>%
  rename(all_adm = n)

preECGAdm_yr <- dates %>%
  filter(endoscopy_dx_gc == 1 & is_cv == 0 & is_preECG == 1) %>%
  group_by(fiscalyr) %>%
  tally() %>%
  rename(preECG_adm = n)

# Merge all admissions count data
endoAdm_yr <- endoAdm_yr %>%
  left_join(preECGAdm_yr, by = "fiscalyr") %>%
  mutate(perc_preECG = preECG_adm/all_adm)


## -----------------------------------------------------------------------------
## SECTION 7: INTEGRATION OF GROUP CLASSIFICATION INTO DATASET WITH BASELINE 
## PATIENT CHARACTERISTICS
## -----------------------------------------------------------------------------

# Mid-year exchange rate from Thai baht to US dollars in 2019
xr <- 31.05

# Subset study cohort (exposed and comparison group). Identify admissions timing relative to index ECG for exposed group and index endoscopy for comparison group.
cohort <- endo %>%
  inner_join(dates_indexEndo, by = 'pid')  %>% 
  group_by(pid) %>%
  arrange(dateadm_new, .by_group = TRUE) %>%
  mutate(pid_count = row_number(pid), .after = pid) %>% 
  # Create indicator variables for whether the admission was the same day as or after the index ECG and index endoscopy
  mutate(postECG = case_when(dateadm_new >= date_preECG ~ 1, TRUE ~ 0)) %>% 
  mutate(postEndo = case_when(dateadm_new >= date_indexEndo ~ 1, TRUE ~ 0)) %>%
  # Convert reimbursed amounts into USD 
  mutate(act_amt_usd = act_amt/xr)

# Count number of unique patients in exposed group
n_distinct(cohort$pid[cohort$is_expPatient == "1"])

# Count number of unique patients in comparison group
n_distinct(cohort$pid[cohort$is_expPatient == "0"])

# Count number of unique patients in exposed group per year
cohort %>% 
  filter(is_expPatient == 1) %>%
  group_by(fiscalyr) %>%
  summarise(unique_pid = n_distinct(pid))

# Count number of excluded patients due to unspecified health region
cohort %>%
  group_by(pid) %>%
  filter(any(rgn1_af == 14)) %>%
  distinct(pid, .keep_all = TRUE) %>%
  group_by(is_expPatient) %>%
  summarise(n = n())

# Exclude patients in unspecified health region
cohort <- cohort %>%
  group_by(pid) %>%
  filter(!any(rgn1_af == 14)) 

## -----------------------------------------------------------------------------
## SECTION 8: PREVALENCE OF PREPROCEDURAL ECG BY ENDOSCOPIC PROCEDURE
## -----------------------------------------------------------------------------

# Count the number of each type of GI endoscopy
procTotals <- cohort %>%
  mutate(proc_adm = strsplit(proc_adm, " ")) %>%
  unnest(proc_adm) %>%
  mutate(endoscopy_dx_gc = ifelse(str_detect(proc_adm,"4311|4411|4412|4413|4414|4521|4522|4523|4524|4525"),1,0), .after = proc_adm) %>%
  filter(endoscopy_dx_gc == 1) %>%
  select(proc_adm) %>%
  group_by(proc_adm) %>%
  count(proc_adm)

# Count number of pre-procedural ECG per procedure type
ECGTotals <- cohort %>%
  mutate(proc_adm = strsplit(proc_adm, " ")) %>%
  unnest(proc_adm) %>% 
  mutate(endoscopy = ifelse(str_detect(proc_adm,"4311|4411|4412|4413|4414|4521|4522|4523|4524|4525"),1,0), .after = proc_adm) %>%
  filter(endoscopy == 1 & is_expPatient == 1) %>%
  select(proc_adm) %>%
  group_by(proc_adm) %>%
  count(proc_adm) %>%
  rename(preECG = n) 

prev_by_proc <- left_join(procTotals, ECGTotals, by = "proc_adm") %>%
  mutate(preECG = ifelse(is.na(preECG), 0, preECG))
prev_by_proc$percentage <- prev_by_proc$preECG/prev_by_proc$n

## -----------------------------------------------------------------------------
## SECTION 9: MONTHLY PREVALENCE OF GI ENDOSCOPY
## -----------------------------------------------------------------------------
# Count monthly GI endoscopy admissions of patients without prior cardiovascular diagnosis
endo_monthly <- dates %>%
  filter(endoscopy_dx_gc == 1 & is_cv == 0) %>%
  mutate(yr_month = substr(date_endo, 1,7)) %>%
  group_by(yr_month) %>%
  count() %>%
  rename(total = n) 

# Count monthly GI endoscopy admissions with pre-procedural ECG
ECG_monthly <- dates %>%
  filter(endoscopy_dx_gc == 1 & is_preECG == 1 & is_cv == 0) %>%
  mutate(yr_month = substr(date_ECG, 1,7)) %>%
  group_by(yr_month) %>%
  count() %>%
  rename(preECG = n)

# Join all monthly admissions data
prevM <- endo_monthly %>%
  left_join(ECG_monthly, by = "yr_month") %>%
  pivot_longer(names_to = "adm_type", cols = !1) %>%
  mutate(yr = substr(yr_month,1,4)) %>%
  mutate(month = as.numeric(substr(yr_month,6,7))) %>%
  mutate(adm_type = factor(adm_type , levels = c("total","preECG")))


## -----------------------------------------------------------------------------
## SECTION 10: ANNUAL PREVALENCE OF GI ENDOSCOPY
## -----------------------------------------------------------------------------
# Count number of GI endoscopy admissions with pre-procedural ECG per fiscal year
ECG_yearly <- indexE %>%
  mutate(fiscalyr = ifelse(month(date_indexEndo)>9, year(date_indexEndo)+1, year(date_indexEndo))) %>%
  group_by(fiscalyr) %>%
  count() %>%
  rename(preECG = n)

# Count number of GI endoscopy admissions with pre-procedural ECG per fiscal year
noECG_yearly <- indexC %>%
  mutate(fiscalyr = ifelse(month(date_indexEndo)>9, year(date_indexEndo)+1, year(date_indexEndo))) %>%
  group_by(fiscalyr) %>%
  count() %>%
  rename(noECG = n)

# Join ECG and non-ECG prevalence counts
prevFY <- noECG_yearly  %>%
  left_join(ECG_yearly, by = "fiscalyr") %>%
  pivot_longer(names_to = "group", cols = !1)  %>%
  group_by(fiscalyr) %>%
  mutate(label = scales::percent(value/sum(value), accuracy = 0.1L),fiscalyr,group) 
