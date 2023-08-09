################################################################################
## Regression model
################################################################################

## -----------------------------------------------------------------------------
## SECTION 1: SETUP
## -----------------------------------------------------------------------------

# Load required libraries
library(sandwich)
library(survey)
library(tidyverse)
library(tableone)
library(car)

# Import cohort (version with unique PID) from script 1 and cascade count data from script 4
cohortUnq <- read.csv("~/cohortUnq.csv")
cascadeCounts <- read.csv("~/cascade_counts.csv")


## -----------------------------------------------------------------------------
## SECTION 2: PROPENSITY SCORE MODEL USING INVERSE PROBABILITY TREATMENT WEIGHTING
## (final model which includes clinically significant variables)
## -----------------------------------------------------------------------------

# Covariates
vars.ps <- c("age_new", "gender_new", "rgn1_af", "cci_group","hbp","hcl")

# Join cohort and cascade count data
cohortUnq <- cohortUnq %>%
  left_join(cascadeCounts, by = "pid") 

# Create binary variables
cohortUnq <- cohortUnq %>%
  mutate(cascade_tx = ifelse(n_cascade_tx > 0, 1, 0),
         cascade_test = ifelse(n_cascade_test > 0, 1, 0),
         cascade_dx = ifelse(n_cascade_dx > 0, 1, 0),
         infection = ifelse(n_infection > 0, 1, 0),
         perforation = ifelse(n_perforation > 0, 1, 0),
         haemorrhage = ifelse(n_haemorrhage > 0, 1, 0),
         transfusion = ifelse(n_transfusion > 0, 1, 0),
         cascade_any = ifelse(n_cascade_any > 0, 1, 0)) %>%

# Reformat variables as factors   
cohortUnq <- cohortUnq %>%
  mutate_if(is.numeric, coalesce, 0) %>%
  mutate_at(c("is_expPatient","gender_new","rgn1_af", "age_group", 
              "cci_group", "hbp", "hcl", "dm",
              "cascade_tx", "cascade_test", "cascade_dx", 
              "infection", "perforation", "haemorrhage", "transfusion", "cascade_any"), as.factor) %>% 
  # Set Bangkok as reference level for health region
  mutate(rgn1_af = relevel(rgn1_af, ref = 13))


# Run regression
model.ps <- glm(is_expPatient ~ age_new + gender_new + rgn1_af + cci_group + hbp + hcl, 
                family = binomial(link="logit"), data = cohortUnq)
summary(model.ps)

model.ps2 <- glm(is_expPatient ~ age_new + gender_new + rgn1_af + cci_group, 
                 family = binomial(link="logit"), data = cohortUnq)
summary(model.ps2)

# Calculate odds ratio from coefficients
exp(coefficients(model.ps))

# Calculate 95% confidence interval for odds ratio
exp(confint(model.ps))

# Calculate log likelihood
logLik(model.ps)

# Calculate propensity score for each patient
ps <- predict(model.ps, type = "response")

# Test for influential observations based on Cook's distance
plot(model.ps, which = 4, id.n = 3)

# Test for multicollinearity
car::vif(model.ps)

# Create weights 
iptw <- ifelse(cohortUnq$is_expPatient == 1, 1/(ps), 1/(1-ps))

# Apply weights to data
cohortUnq.wt <- svydesign(ids = ~1, data = cohortUnq, weights = ~iptw)

# Create weighted table to check balance
tab.wt <- svyCreateTableOne(vars = vars.ps, strata = "is_expPatient", data = cohortUnq.wt, test = FALSE)
print(tab.wt, smd = TRUE)

# Truncate weights at 99 percentile
max_wt <- quantile(iptw, 0.99)
iptw.trunc <- replace(iptw, iptw > max_wt, max_wt)

# Standardised residuals plot
augment(model.ps) %>%
  mutate(index = 1:n()) %>%
  ggplot(aes(index, .std.resid)) +
  geom_point(aes(color = is_expPatient), alpha = 0.6) +
  theme_bw() +
  labs(x = "Index", y = "Standardized residuals") +
  scale_color_manual(values = c("darkgrey","#4480B0"), labels = c("EKG","Comparison"), name = "Group") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10))


## -----------------------------------------------------------------------------
## SECTION 3: PROPENSITY SCORE MODEL USING INVERSE PROBABILITY TREATMENT WEIGHTING
## (model achieved by backward elimination)
## -----------------------------------------------------------------------------

# Covariates 
vars.ps2 <- c("age_new", "gender_new", "rgn1_af", "cci_group")

# Run regression
model.ps2 <- glm(is_expPatient ~ age_new + gender_new + rgn1_af + cci_group, 
                 family = binomial(link="logit"), data = cohortUnq)
summary(model.ps2)

# Calculate odds ratio from coefficients
exp(coefficients(model.ps2))

# Calculate 95% confidence interval for odds ratio
exp(confint(model.ps2))

# Calculate log likelihood
logLik(model.ps2)

# Calculate propensity score for each patient
ps2 <- predict(model.ps2, type = "response")

# Test for multicollinearity
car::vif(model.ps2)

# Create weights 
iptw2 <- ifelse(cohortUnq$is_expPatient == 1, 1/(ps2), 1/(1-ps2))

# Apply weights to data
cohortUnq.wt2 <- svydesign(ids = ~1, data = cohortUnq, weights = ~iptw2)

# Create weighted table to check balance
tab.wt2 <- svyCreateTableOne(vars = vars.ps2, strata = "is_expPatient", data = cohortUnq.wt2, test = FALSE)
print(tab.wt2, smd = TRUE)

# Truncate weights at 99 percentile
max_wt2 <- quantile(iptw2, 0.99)
iptw.trunc2 <- replace(iptw2, iptw2 > max_wt2, max_wt2)

# Standardised residuals plot
augment(model.ps2) %>%
  mutate(index = 1:n()) %>%
  ggplot(aes(index, .std.resid*0.85)) +
  geom_point(aes(color = is_expPatient), alpha = 0.6) +
  theme_bw() +
  labs(x = "Index", y = "Standardized residuals") +
  scale_color_manual(values = c("darkgrey","#4480B0"), labels = c("EKG","Comparison"), name = "Group") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10))


## -----------------------------------------------------------------------------
## SECTION 4: CALCULATION OF RELATIVE RISKS 
## -----------------------------------------------------------------------------

# Regression model for relative risk of a cascade test given receipt of ECG
regRR.test <- glm(cascade_test ~ is_expPatient, weights = iptw.trunc, family = quasibinomial(link = log), data = cohortUnq)
betaRR.test <- coef(regRR.test)[2]
seRR.test <- sqrt(diag(vcovHC(regRR.test, type = "HC0")))
# Relative risk with 95% CI
rr.test <- exp(betaRR.test)
lclRR.test <- exp(betaRR.test - 1.96*seRR.test[2])
uclRR.test <- exp(betaRR.test + 1.96*seRR.test[2])
c(lclRR.test, rr.test, uclRR.test)

# Regression model for relative risk of a cascade test (test + no treatment) given receipt of ECG
cohortUnq.noCTx <- cohortUnq[cohortUnq$cascade_tx == 0,]
iptw.trunc.noCTx <- iptw.trunc[which(cohortUnq$cascade_tx == 0)]
regRR.test.noCTx <- glm(cascade_test ~ is_expPatient, 
                        weights = iptw.trunc.noCTx, family = quasibinomial(link = log), data = cohortUnq.noCTx )
betaRR.test.noCTx <- coef(regRR.test.noCTx)[2]
seRR.test.noCTx <- sqrt(diag(vcovHC(regRR.test.noCTx, type = "HC0")))
# Relative risk with 95% CI
rr.test.noCTx <- exp(betaRR.test.noCTx)
lclRR.test.noCTx <- exp(betaRR.test.noCTx - 1.96*seRR.test.noCTx[2])
uclRR.test.noCTx <- exp(betaRR.test.noCTx + 1.96*seRR.test.noCTx[2])
c(lclRR.test.noCTx, rr.test.noCTx, uclRR.test.noCTx)

# Regression model for relative risk of a cascade treatment
regRR.tx <- glm(cascade_tx ~ is_expPatient, weights = iptw.trunc, family = quasibinomial(link = log), data = cohortUnq)
betaRR.tx <- coef(regRR.tx)[2]
seRR.tx <- sqrt(diag(vcovHC(regRR.tx, type = "HC0")))
rr.tx <- exp(betaRR.tx)
# Relative risk with 95% CI
lclRR.tx <- exp(betaRR.tx - 1.96*seRR.tx[2])
uclRR.tx <- exp(betaRR.tx + 1.96*seRR.tx[2])
c(lclRR.tx, rr.tx, uclRR.tx)


# Regression model for relative risk of a cascade treatment given cascade test
cohortUnq.cTest <- cohortUnq[cohortUnq$cascade_test == 1,]
iptw.trunc.cTest <- iptw.trunc[which(cohortUnq$cascade_test == 1)]
regRR.tx.givenTest <- glm(cascade_tx ~ is_expPatient, weights = iptw.trunc.cTest, family = quasibinomial(link = log), data = cohortUnq.cTest)
betaRR.tx.givenTest <- coef(regRR.tx.givenTest)[2]
seRR.tx.givenTest <- sqrt(diag(vcovHC(regRR.tx.givenTest, type = "HC0")))
# Relative risk with 95% CI
rr.tx.givenTest <- exp(betaRR.tx.givenTest)
lclRR.tx.givenTest <- exp(betaRR.tx.givenTest - 1.96*seRR.tx.givenTest[2])
uclRR.tx.givenTest <- exp(betaRR.tx.givenTest + 1.96*seRR.tx.givenTest[2])
c(lclRR.tx.givenTest, rr.tx.givenTest, uclRR.tx.givenTest)


# Regression model for relative risk of a hospital-acquired infection
regRR.inf <- glm(infection ~ is_expPatient, weights = iptw.trunc, family = quasibinomial(link = log), data = cohortUnq)
betaRR.inf <- coef(regRR.inf)[2]
seRR.inf <- sqrt(diag(vcovHC(regRR.inf, type = "HC0")))
# Relative risk with 95% CI
rr.inf <- exp(betaRR.inf)
lclRR.inf <- exp(betaRR.inf - 1.96*seRR.inf[2])
uclRR.inf <- exp(betaRR.inf + 1.96*seRR.inf[2])
c(lclRR.inf, rr.inf, uclRR.inf)


# Regression model for relative risk of a perforation
regRR.perf <- glm(perforation ~ is_expPatient, weights = iptw.trunc, family = quasibinomial(link = log), data = cohortUnq)
betaRR.perf <- coef(regRR.perf)[2]
seRR.perf <- sqrt(diag(vcovHC(regRR.perf, type = "HC0")))
# Relative risk with 95% CI
rr.perf <- exp(betaRR.perf)
lclRR.perf <- exp(betaRR.perf - 1.96*seRR.perf[2])
uclRR.perf <- exp(betaRR.perf + 1.96*seRR.perf[2])
c(lclRR.perf, rr.perf, uclRR.perf)


# Regression model for relative risk of a haemorrhage
regRR.hh <- glm(haemorrhage ~ is_expPatient, weights = iptw.trunc, family = quasibinomial(link = log), data = cohortUnq)
betaRR.hh <- coef(regRR.hh)[2]
seRR.hh <- sqrt(diag(vcovHC(regRR.hh, type = "HC0")))
# Relative risk with 95% CI
rr.hh <- exp(betaRR.hh)
lclRR.hh <- exp(betaRR.hh - 1.96*seRR.hh[2])
uclRR.hh <- exp(betaRR.hh + 1.96*seRR.hh[2])
c(lclRR.hh, rr.hh, uclRR.hh)


# Regression model for relative risk of a transfusion
regRR.trans <- glm(transfusion ~ is_expPatient, weights = iptw.trunc, family = quasibinomial(link = log), data = cohortUnq)
betaRR.trans <- coef(regRR.trans)[2]
seRR.trans <- sqrt(diag(vcovHC(regRR.trans, type = "HC0")))
rr.trans <- exp(betaRR.trans)
# Relative risk with 95% CI
lclRR.trans <- exp(betaRR.trans - 1.96*seRR.trans[2])
uclRR.trans <- exp(betaRR.trans + 1.96*seRR.trans[2])
c(lclRR.trans, rr.trans, uclRR.trans)


