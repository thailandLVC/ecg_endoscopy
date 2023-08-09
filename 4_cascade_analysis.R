################################################################################
## Identification of downstream cardiac tests and procedures, cascade tests, 
## nosocomial infection, and procedure-related complications
################################################################################

## -----------------------------------------------------------------------------
## SECTION 1: SETUP AND DATA PREPARATION
## -----------------------------------------------------------------------------

# Load required libaries
library(tidyverse)

# Import cohort data from script 1
cohort <- read.csv("~/cohort.csv")

# Subset admissions within cascade period: 1 year of index EKG (exposed group) or index endoscopy (comparison group)
postAdm <- cohort %>%
  # First subset admissions after index EKG  or index endoscopy for quicker processing
  filter(postEKG == 1 | postEndo == 1) %>%
  # Calculate when visit occurred relative to index EKG and index endoscopy
  mutate(daysPostEKG = as.numeric(as.Date(dateadm_new) - as.Date(date_preEKG)),
         daysPostEndo = as.numeric(as.Date(dateadm_new) - as.Date(date_indexEndo))) %>%
  # Create indicator variables for 30-day readmissions and 14-day readmissions
  mutate(readm_30d = case_when(between(daysPostEndo,1,30) ~ 1, TRUE ~ 0),
         readm_14d = case_when(between(daysPostEndo,1,14) ~ 1, TRUE ~ 0)) %>%
  # Create an indicator variable for visits within the cascade period
  mutate(cascade_period = case_when(daysPostEndo >= 0 & daysPostEndo <= 90 | daysPostEKG >= 0 & daysPostEndo <= 90 ~ 1,
                                    TRUE ~ 0)) %>%
  # Combine primary and secondary diagnoses into one column
  unite("dx", c(pdx, sdx_adm), na.rm = TRUE, sep = " ")


# Calculate mean of pre-procedural period length
postAdm %>%
  # Keep if EKG was completed before endoscopy (i.e., not same day)
  filter(daysPostEKG == 0 & daysPostEndo <= 0) %>%
  ungroup() %>%
  summarise(median = median(daysPostEndo ))


## -----------------------------------------------------------------------------
## SECTION 2: IDENTIFICATION AND CHARACTERISATION OF KEY EVENTS DOWNSTREAM OF ECG
## -----------------------------------------------------------------------------

# Identify cardiac care cascade events
cascadeCounts <- postAdm %>%
  filter(daysPostEKG > 0 | daysPostEndo > 0) %>%
  select(tran_id, pid, pid_count, is_expPatient, proc_adm, dx, dateadm_new, daysPostEKG, daysPostEndo) %>%
  mutate(cascade_tx = str_count(proc_adm, "^350|^351|^352|^355|^356|
                                           ^357|^358|^359|^360|^361|^362|^363|
                                           ^369|^377|^378|3964"),
         cascade_test = str_count(proc_adm, "8951|8952|3728|3720|3726|
                                             3727|3721|3722|3723|3824|^885|8892|9205|
                                             8941|8942|8943|8944|8945|8946|8947|8948|
                                             8872|8749|8950|8953|8955|8956|8957|3729"),
         cascade_dx = str_count(dx, "I05|I06|I07|I08|I09|I11|I13|I20|I21|I22|I23|I24|
                           I25|I33|I34|I35|I36|I37|I38|I39|I44|I45|I47|I48|I49|I50|R071|R072|
                           R078|R079|R060|R001|I95|R224|R000|R002|R42|R55|I469|Z950")) %>%
  group_by(pid,is_expPatient) %>%
  summarise(n_cascade_tx = sum(cascade_tx),
            n_cascade_test = sum(cascade_test),
            n_cascade_dx = sum(cascade_dx))


# Identify potential procedure-related complications
proc.cxCounts <- postAdm %>%
  filter(readm_14d == 1)  %>%
  select(tran_id, pid, pid_count, is_expPatient, proc_adm, dx, dateadm_new, daysPostEKG, daysPostEndo) %>%
  mutate(perforation = str_count(dx, "K223|K251|K252|K255|K256|K261|K262|K265|K266|K271|
                                 K272|K275|K276|K281|K282|K285|K286|K570|K572|K574|K575|
                                 K578|K630|K631|K632|P780|S363|S364|S365|S366|S367|S369"),
         haemorrhage = str_count(dx, "D699|K250|K252|K254|K256|K260|K262|K264|K266|K270|
                                 K272|K274|K276|K280|K282|K284|K286|K290|K625|K661|K922"),
         transfusion = str_count(proc_adm, "^990| 990")) %>%
  group_by(pid,is_expPatient) %>%
  summarise(n_perforation = sum(perforation),
            n_haemorrhage = sum(haemorrhage),
            n_transfusion = sum(transfusion))

# Identify incidences of hospital-acquired infections and merge with procedure-related complication data
cxCounts <- postAdm %>%
  filter(daysPostEndo == 0 | readm_30d == 1) %>%
  select(tran_id, pid, pid_count, is_expPatient, dx, dateadm_new, daysPostEKG, daysPostEndo) %>%
  mutate(infection = str_count(dx, "A047|A048|A084|A090|A099|A419|T814|
                               B599|J100|J101|J111|J118|J121|J123|J129|J139|J149|J150|
                               J151|J152|J153|J154|J155|J156|J157|J158|J159|J168|J172|
                               J178|J180|J181|J188|J189|N109|N129|N151|N300|N308|N309|
                               N340|N341|N342|N390|O862|T835|J985|K630|K650|K658|K659|
                               K750|L033|L089|L730|L033|L089|L730|L732|M000|M001|M002|
                               M008|M009|M861|O860|O861|O864|O868|A021|A400|A401|A402|
                               A403|A408|A409|A410|A411|A412|A413|A414|A415|A418|A419|
                               B007|B377|O859")) %>%
  group_by(pid,is_expPatient) %>%
  summarise(n_infection = sum(infection)) %>%
  full_join(proc.cxCounts, by = "pid") %>%
  mutate_if(is.numeric, coalesce, 0)


# Calculate the number of infections that occurred during the index endoscopy admission versus a 30-day readmission
n_infectionTime <- postAdm %>%
  filter(daysPostEndo == 0 | readm_30d == 1) %>%
  select(tran_id, pid, pid_count, is_expPatient, dx, dateadm_new, daysPostEKG, daysPostEndo) %>%
  mutate(dx = strsplit(dx, " ")) %>%
  unnest(dx) %>% 
  mutate(infection = case_when(str_detect(dx, "A047|A048|A084|A090|A099|A419|T814|
                               B599|J100|J101|J111|J118|J121|J123|J129|J139|J149|J150|
                               J151|J152|J153|J154|J155|J156|J157|J158|J159|J168|J172|
                               J178|J180|J181|J188|J189|N109|N129|N151|N300|N308|N309|
                               N340|N341|N342|N390|O862|T835|J985|K630|K650|K658|K659|
                               K750|L033|L089|L730|L033|L089|L730|L732|M000|M001|M002|
                               M008|M009|M861|O860|O861|O864|O868|A021|A400|A401|A402|
                               A403|A408|A409|A410|A411|A412|A413|A414|A415|A418|A419|
                               B007|B377|O859") ~ 1, TRUE ~ 0)) %>%
  filter(infection == 1) %>%
  select(daysPostEndo,is_expPatient) %>%
  group_by(daysPostEndo,is_expPatient) %>%
  count(daysPostEndo)


# Count number of cardiac procedures per patient
n_cascade_tx <- as.data.frame(table(cascadeCounts$n_cascade_tx, cascadeCounts$is_expPatient)) %>%
  rename(n_cascade_tx = Var1, is_expPatient = Var2)

# Count number of cardiac tests per patient
n_cascade_test <- as.data.frame(table(cascadeCounts$n_cascade_test, cascadeCounts$is_expPatient)) %>%
  rename(n_cascade_tx = Var1, is_expPatient = Var2)

# Count number of cascade tests (cardiac tests not followed by treatment) per patient
n_cascade_test_notx <- as.data.frame(table(cascadeCounts$n_cascade_test[cascadeCounts$n_cascade_tx == 0], 
                                           cascadeCounts$is_expPatient[cascadeCounts$n_cascade_tx == 0])) %>%
  rename(n_cascade_tx = Var1, is_expPatient = Var2)


## -----------------------------------------------------------------------------
## SECTION 3: IDENTIFICATION OF CASCADE TEST SEQUENCES FOR DECISION TREE MODEL
## -----------------------------------------------------------------------------

# Identify and group tests into one of five categories: A = electrocardiogram, 
# B = electrophysiology or catheterisation, C = imaging, D = exercise stress test or monitoring, E = xray

cascade <-  postAdm %>%
  filter(cascade_period == 1) %>%
  group_by(pid) %>%
  select(tran_id, pid, pid_count, is_expPatient, proc_adm, act_amt, dateadm_new, daysPostEKG, daysPostEndo) %>%
  mutate(proc_adm = strsplit(proc_adm, " ")) %>%
  unnest(proc_adm) %>%
  mutate(event = case_when(str_detect(proc_adm, "8951|8952") ~ "A",
                           str_detect(proc_adm, "3720|3721|3722|3723|3726|3727|3729") ~ "B",
                           str_detect(proc_adm, "3824|^885|8892|9205") ~ "C",
                           str_detect(proc_adm, "8941|8942|8943|8944|8950|8953|8955|8956|8957") ~ "D",
                           proc_adm == "8749" ~ "E",
                           str_detect(proc_adm, "^350|^351|^352|^361|^355|^356|^357|^358|^359|
                                      ^360|^362|^363|^369|^377|^378|3964") ~ "X")) %>%
  filter(!is.na(event))


# Subset cascade sequences of ECG group
cascade.E <- cascade %>%
  filter(is_expPatient == 1) %>%
  group_by(pid) %>%
  arrange(dateadm_new, event, .by_group = TRUE) %>%
  mutate(pid_count = row_number()) %>%
  mutate(daysPostEKG_range = case_when(daysPostEKG >= 0 & daysPostEKG < 30 ~ "01",
                                       daysPostEKG >= 30 & daysPostEKG < 60 ~ "02",
                                       daysPostEKG >= 60 & daysPostEKG < 90 ~ "03",
                                       daysPostEKG >= 90 & daysPostEKG < 120 ~ "04",
                                       daysPostEKG >= 120 & daysPostEKG < 150 ~ "05",
                                       daysPostEKG >= 150 & daysPostEKG < 180 ~ "06",
                                       daysPostEKG >= 180 & daysPostEKG < 210 ~ "07",
                                       daysPostEKG >= 210 & daysPostEKG < 240 ~ "08",
                                       daysPostEKG >= 240 & daysPostEKG < 270 ~ "09",
                                       daysPostEKG >= 270 & daysPostEKG < 300 ~ "10",
                                       daysPostEKG >= 300 & daysPostEKG < 330 ~ "11",
                                       daysPostEKG >= 330 & daysPostEKG < 365 ~ "12"))


# Restructure data so that each row contains one cascade sequence
cascade.E.wide <- cascade.E %>%
  group_by(pid) %>%
  mutate(cascade_seq = str_trim(paste(event, collapse = "-"))) %>%
  select(pid, cascade_seq) %>%
  distinct(pid, .keep_all = TRUE)

# Count frequency of each sequence for ECG group
patternFreq.E <- cascade.E.wide %>%
  group_by(cascade_seq) %>%
  tally()

# Subset cascade sequences of comparison group
cascade.C <- cascade %>%
  filter(is_expPatient == 0) %>%
  group_by(pid) %>%
  arrange(dateadm_new, event, .by_group = TRUE) %>%
  mutate(pid_count = row_number()) %>%
  mutate(daysPostEKG_range = case_when(daysPostEKG >= 0 & daysPostEKG < 30 ~ "01",
                                       daysPostEKG >= 30 & daysPostEKG < 60 ~ "02",
                                       daysPostEKG >= 60 & daysPostEKG < 90 ~ "03",
                                       daysPostEKG >= 90 & daysPostEKG < 120 ~ "04",
                                       daysPostEKG >= 120 & daysPostEKG < 150 ~ "05",
                                       daysPostEKG >= 150 & daysPostEKG < 180 ~ "06",
                                       daysPostEKG >= 180 & daysPostEKG < 210 ~ "07",
                                       daysPostEKG >= 210 & daysPostEKG < 240 ~ "08",
                                       daysPostEKG >= 240 & daysPostEKG < 270 ~ "09",
                                       daysPostEKG >= 270 & daysPostEKG < 300 ~ "10",
                                       daysPostEKG >= 300 & daysPostEKG < 330 ~ "11",
                                       daysPostEKG >= 330 & daysPostEKG < 365 ~ "12"))

# Restructure data so that each row contains one cascade sequence
cascade.C.wide <- cascade.C %>%
  group_by(pid) %>%
  mutate(cascade_seq = str_trim(paste(event, collapse = "-"))) %>%
  select(pid, cascade_seq) %>%
  distinct(pid, .keep_all = TRUE)

# Count frequency of each sequence for comparison group
patternFreq.C <- cascade.C.wide %>%
  group_by(cascade_seq) %>%
  tally()


