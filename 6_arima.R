################################################################################
## ARIMA model
################################################################################

## -----------------------------------------------------------------------------
## SECTION 1: SETUP
## -----------------------------------------------------------------------------
# Load required libraries
library(astsa)
library(forecast)
library(zoo)
library(plyr)
library(tidyverse)
library(lubridate)
library(showtext)

# Customise theme
font_add_google(name = "Roboto", family = "Roboto")
showtext_auto()

theme_mgh <- function(){
  font <- "Roboto"
  base_size <- 11
  theme_classic() + 
    theme(panel.grid.major.y = element_line(color = "gray80", linetype = "dotted"),
          axis.text = element_text(base_size, family = font),
          axis.ticks.y = element_blank(),
          axis.line.x = element_line(),
          axis.line.y = element_blank(),
          axis.title = element_text(base_size, family = font),
          axis.title.x = element_text(margin = margin(t = 8, r = 0, b = 0, l = 0)),
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
          legend.title = element_blank(),
          legend.background = element_blank())
}

# Import monthly prevalence data
raw_df <- read.csv("~/timeseries.csv")


## -----------------------------------------------------------------------------
## SECTION 2: DATA PREPARATION
## -----------------------------------------------------------------------------
# Format data for time series analysis
raw_df <- raw_df %>%
  # Format time as date
  mutate(time = as.Date(time)) %>%
  # Filter for dates to exclude 2021 (due to incomplete data) and COVID-19 lockdown period
  filter(time < "2020-03-01") 

# Divide data into a training period and test period. Convert both into time series objects.
ts_train <- ts(raw_df$prev, start = c(2016,1), end = c(2019,2), frequency = 12)
ts_test <- ts(raw_df$prev, start = c(2019,3), end = c(2020,02), frequency = 12)

# View training data
plot(ts_train, type='l')

## -----------------------------------------------------------------------------
## SECTION 3: MODEL IDENTIFICATION
## -----------------------------------------------------------------------------

# Plot ACF and PACF of undifferenced data
acf(ts_train, lag.max = 50)
pacf(ts_train, lag.max = 50)

# Determine the differencing order
ndiffs(ts_train, test = "adf")

# Make time series stationary by taking a seasonal difference then first difference
ts_train_sty <- diff(diff(ts_train,12))

# Plot ACF and PACF of differenced data
acf(ts_train_sty, lag.max = 50)
pacf(ts_train_sty, lag.max = 50)

# Identify seasonal ARIMA model parameters (p,q,P,Q) using an algorithm that iteratively searches for the best fitting model. Set order of first differencing (d) = 1 and order of seasonal differencing (D) = 1
arima1 <- auto.arima(ts_train, seasonal = TRUE, max.d = 1, max.D = 1, stepwise = FALSE, trace = TRUE)

# Check residuals for autocorrelation
checkresiduals(arima1)
Box.test(arima1$residuals, type = "Ljung-Box")

# Return parameters and confidence intervals
summary(arima1)
confint(arima1)

# Store model results in a dataframe
fitted_df <- data.frame(time = seq.Date(from = as.Date("2016-01-01"), length.out = 38, by = "month"),
                        prev = raw_df$prev[1:38],
                        fitted = as.numeric(arima1$fitted))

# Compare model output to input data
ggplot(fitted_df[c(13:38),], aes(x = time)) +
  geom_line(aes(y = prev, col = "prev")) +
  geom_line(aes(y = fitted, colour = "fitted"), lty = 2) +
  scale_x_date(limit = c(as.Date("2017-01-01"), as.Date("2019-02-01")),
               breaks = seq.Date(as.Date("2017-01-01"), as.Date("2019-02-01"), by = "month"), 
               date_labels = "%b %Y", expand = c(0,0)) +
  scale_y_continuous(limits = c(0,4000), expand = c(0,0)) + 
  scale_color_manual("", breaks = c("prev", "fitted"),
                     values = c("prev" = "dimgray", "fitted" = "steelblue"), 
                     labels = c("Observed", "Fitted")) + 
  theme_mgh() +
  theme(axis.text.x = element_text(angle = 90), text = element_text(size = 12), axis.title = element_text(size = 10),
        legend.position = c(0.85,0.15), plot.margin = unit(c(0.5,0.5,0,0), "cm")) +
  labs(x = "Time", y = "Number of low-risk \n endoscopy patients") 


## -----------------------------------------------------------------------------
## SECTION 4: MODEL VALIDATION
## -----------------------------------------------------------------------------

# Assess model using test data
pred_testdata <- forecast(arima1, h = 12, level = 95)

# Store validation results in a dataframe
pred_testdata_df <- data.frame(time = seq.Date(from = as.Date("2019-03-01"), length.out = 12, by = "month"),
                               pred_mean = as.numeric(pred_testdata$mean), 
                               pred_lower = as.numeric(pred_testdata$lower), 
                               pred_upper = as.numeric(pred_testdata$upper))

# Bind validation and test results in one dataframe
out1 <- rbind.fill(fitted_df[c(1:38),c("time","fitted")], pred_testdata_df)

# Plot
ggplot() +
  geom_ribbon(data = pred_testdata_df, aes(x = time, ymin = pred_lower, ymax = pred_upper), 
              fill = "steelblue", alpha = 0.2) +
  geom_line(data = pred_testdata_df, aes(x = time, y = pred_mean), col = "steelblue", lty = 2, size = 0.5) +
  geom_line(data = raw_df[39:50,], aes(x = time, y = prev), col = "dimgray", lty = 1, size = 0.5) +
  scale_x_date(breaks = seq.Date(from = as.Date("2019-03-01"), to = as.Date("2020-02-01"), by = "month"), 
               date_labels = "%b %Y", expand = c(0,0)) +
  scale_y_continuous(limits = c(0,8000), breaks = seq(0,8000,2000), 
                     labels = seq(0,8000,2000), expand = c(0,0)) +
  theme_mgh() +
  theme(axis.text.x = element_text(angle = 90), text = element_text(size = 12), axis.title = element_text(size = 10),
        plot.margin = unit(c(0.5,0.5,0,0), "cm")) +
  labs(x = "Time", y = "Number of low-risk \n endoscopy patients") 


## -----------------------------------------------------------------------------
## SECTION 5: PROJECTIONS
## -----------------------------------------------------------------------------

# Project prevalence until 2028
pred <- forecast(arima1, h = 106, level = 75)

# Store projection results in a dataframe
pred_df <- data.frame(time = seq.Date(from = as.Date("2020-03-01"), length.out = 106, by = "month"),
                      pred_mean = as.numeric(pred$mean), 
                      pred_lower = as.numeric(pred$lower), 
                      pred_upper = as.numeric(pred$upper))


# Append projection results to original time series data
out2 <- rbind.fill(fitted_df[c(1:38),c("time","fitted")], pred_testdata_df, pred_df)

# Plot time series 
ggplot() +
  geom_rect(aes(xmin = as.Date("2020-03-01"), xmax = as.Date("2020-06-30"), ymin = -Inf, ymax = Inf), 
            fill = "#f2f2f0", alpha = 0.1) +
  geom_ribbon(data = out2[c(39:156),], aes(x = time, ymin = pred_lower, ymax = pred_upper), 
              fill = "#1567b6", alpha = 0.2) +
  geom_line(data = out2[c(39:156),], aes(x = time, y = pred_mean), col = "#2B3956", lty = 2, size = 0.6) +
  geom_line(data = out2[c(1:38),], aes(x = time, y = fitted), col = "#2B3956", lty = 1, size = 0.75) +
  scale_x_date(limits = c(as.Date("2016-01-01"), as.Date("2029-01-01")), 
               breaks = seq.Date(from = as.Date("2016-01-01"), to = as.Date("2029-01-01"), 
                                 by = "year"), date_labels = "%Y", expand = c(0,0)) +
  scale_y_continuous(limits = c(0,1.6e4), breaks = seq(0,1.4e4,2e3), 
                     labels = seq(0,1.4e4,2e3), expand = c(0,0)) +
  theme_mgh() +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), axis.title = element_text(size = 10)) +
  labs(x = "Year", y = "Number of low-risk endoscopy patients") 


# Calculate projected annual prevalence from 2024 to 2028
predYr <- out2 %>%
  filter(time >= "2023-10-01" & time < "2028-10-01") %>%
  select(-fitted)  %>%
  mutate(fiscalyr = ifelse(month(time)>9, year(time)+1, year(time))) %>%
  group_by(fiscalyr) %>%
  summarise(total = round(sum(pred_mean),0), min = round(sum(pred_lower),0), max = round(sum(pred_upper),0)) %>%
  mutate(lvc = round(total*0.062,0))
