################################################################################
## Sensitivity analysis for costs of decision tree model
## (probabilities and prevalence held at baseline)
################################################################################


## -----------------------------------------------------------------------------
## SECTION 1: SETUP
## -----------------------------------------------------------------------------

# Load required libraries
library(data.tree)
library(DiagrammeR)
library(stringi)
library(tidyverse)
library(rsvg)
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


# Import cascade data for exposed and comparison groups
cascadeE <- read.csv("~/cascade_seq_ecg.csv")
cascadeC <- read.csv("~/cascade_seq_comp.csv")

# Import cost data
costs <- read.csv("~/cascade_costs.csv")

# Import projected prevalence from ARIMA model
proj <- read.csv("~/projected_prev.csv")


## -----------------------------------------------------------------------------
## SECTION 2: DATA PREPARATION
## -----------------------------------------------------------------------------

# Count number of each group that was admitted during the cascade period
n_EKG <- sum(cascadeE$n)
n_noEKG <- sum(cascadeC$n)

# Add 'No EKG' as the first event of cascade sequence for comparison group
cascadeC <- cascadeC %>%
  mutate(cascade_seq = paste0("N-",cascade_seq)) %>%
  # Add row for patients who did receive an EKG and were not readmitted during the cascade period
  add_row(cascade_seq = "N", n = 55443 - n_noEKG)

# Merge cascade datasets for exposed and comparison groups
cascade <- rbind(cascadeE, cascadeC)

# Remove events after treatment and sum counts for unique cascade sequences
cascade <- cascade %>%
  mutate_at("cascade_seq", ~gsub("F","E",.)) %>%
  mutate_at("cascade_seq", ~gsub("X.*","X",.)) %>%
  group_by(cascade_seq) %>%
  summarise(n = sum(n))

# Get maximum length of cascade sequence
max_length <- cascade %>%
  mutate(cascade_seq = gsub("-", "", cascade_seq)) %>%
  summarise(max = max(nchar(cascade_seq))) %>%
  pull()

# Make column names of new dataset
new_cols <- c(paste0("event",as.character(1:max_length)))

# Pivot dataset to wide format and separate cascade sequence into different columns
cascadeWide <- cascade %>%
  separate(cascade_seq, new_cols, fill = "right")

# Set seed
set.seed(6)



## -----------------------------------------------------------------------------
## SECTION 3: PROBABILITIES
## -----------------------------------------------------------------------------

# Node 1: EKG vs. no EKG
node1 <- cascadeWide %>%
  # Count EKG versus no EKG
  group_by(event1) %>% 
  summarise(count = sum(n)) %>%
  # Calculate probability
  ungroup() %>%
  mutate(prob = prop.table(count)) %>%
  rename(pathString = event1)


# Node 2
node2 <- cascadeWide %>%
  # Count each 2-event sequence
  group_by(event1, event2) %>% 
  summarise(count = sum(n)) %>%
  # Calculate probability
  ungroup() %>%
  group_by(event1) %>%
  mutate(prob = prop.table(count)) %>%
  # If cascade discontinued, name 2nd event as "0"
  mutate_at("event2", ~replace_na(.,"0")) %>%
  unite("pathString", c(event1:event2), na.rm = TRUE, sep = "/")


# Node 3
node3 <- cascadeWide %>%
  # Exclude if no 2nd event occurred or the 2nd event was a treatment
  filter(!is.na(event2) & event2 != "X") %>%
  # Count each 3-event sequence
  group_by_at(vars(c(event1:event3))) %>% 
  summarise(count = sum(n)) %>%
  # Calculate probability
  ungroup() %>%
  group_by_at(vars(c(event1:event2))) %>%
  mutate(prob = prop.table(count)) %>%
  # If cascade discontinued, name 3rd event as "0"
  mutate_at("event3", ~replace_na(.,"0")) %>%
  unite("pathString", c(event1:event3), na.rm = TRUE, sep = "/")


# Node 4
node4 <- cascadeWide %>%
  # Exclude if no 3rd event occurred or the 3rd event was a treatment
  filter(!is.na(event3) & event3 != "X") %>%
  # Count each 4-event sequence
  group_by_at(vars(c(event1:event4))) %>% 
  summarise(count = sum(n)) %>%
  # Calculate probability
  ungroup() %>%
  group_by_at(vars(c(event1:event3))) %>%
  mutate(prob = prop.table(count)) %>%
  # If cascade discontinued, name 4th event as "0"
  mutate_at("event4", ~replace_na(.,"0")) %>%
  unite("pathString", c(event1:event4), na.rm = TRUE, sep = "/")


# Node 5
node5 <- cascadeWide %>%
  # Exclude if no 4th event occurred or the 4th event was a treatment
  filter(!is.na(event4) & event4 != "X") %>%
  # Count each 5-event sequence
  group_by_at(vars(c(event1:event5))) %>% 
  summarise(count = sum(n)) %>%
  # Calculate probability
  ungroup() %>%
  group_by_at(vars(c(event1:event4))) %>%
  mutate(prob = prop.table(count)) %>%
  # If cascade discontinued, name 5th event as "0"
  mutate_at("event5", ~replace_na(.,"0")) %>%
  unite("pathString", c(event1:event5), na.rm = TRUE, sep = "/")


# Node 6
node6 <- cascadeWide %>%
  # Exclude if no 5th event occurred or the 5th event was a treatment
  filter(!is.na(event5) & event5 != "X") %>%
  # Count each 6-event sequence
  group_by_at(vars(c(event1:event6))) %>% 
  summarise(count = sum(n)) %>%
  # Calculate probability
  ungroup() %>%
  group_by_at(vars(c(event1:event5))) %>%
  mutate(prob = prop.table(count)) %>%
  # If cascade discontinued, name 6th event as "0"
  mutate_at("event6", ~replace_na(.,"0")) %>%
  unite("pathString", c(event1:event6), na.rm = TRUE, sep = "/")


# Node 7
node7 <- cascadeWide %>%
  # Exclude if no 6th event occurred or the 6th event was a treatment
  filter(!is.na(event6) & event6 != "X") %>%
  # Count each 7-event sequence
  group_by_at(vars(c(event1:event7))) %>%  
  summarise(count = sum(n)) %>%
  # Calculate probability
  ungroup() %>%
  group_by_at(vars(c(event1:event6))) %>%
  mutate(prob = prop.table(count)) %>%
  # If cascade discontinued, name 7th event as "0"
  mutate_at("event7", ~replace_na(.,"0")) %>%
  unite("pathString", c(event1:event7), na.rm = TRUE, sep = "/")


# Node 8
node8 <- cascadeWide %>%
  # Exclude if no 7th event occurred or the 7th event was a treatment
  filter(!is.na(event7) & event7 != "X") %>%
  # Count each 8-event sequence
  group_by_at(vars(c(event1:event8))) %>%
  summarise(count = sum(n)) %>%
  # Calculate probability
  ungroup() %>%
  group_by_at(vars(c(event1:event7))) %>%
  mutate(prob = prop.table(count)) %>%
  # If cascade discontinued, name 8th event as "0"
  mutate_at("event8", ~replace_na(.,"0")) %>%
  unite("pathString", c(event1:event8), na.rm = TRUE, sep = "/")


# Node 9
node9 <- cascadeWide %>%
  # Exclude if no 8th event occurred or the 8th event was a treatment
  filter(!is.na(event8) & event8 != "X") %>%
  # Count each 9-event sequence
  group_by_at(vars(c(event1:event9))) %>%
  summarise(count = sum(n)) %>%
  # Calculate probability
  ungroup() %>%
  group_by_at(vars(c(event1:event8))) %>%
  mutate(prob = prop.table(count)) %>%
  # If cascade discontinued, name 9th event as "0"
  mutate_at("event9", ~replace_na(.,"0")) %>%
  unite("pathString", c(event1:event9), na.rm = TRUE, sep = "/")


# Node 10
node10 <- cascadeWide %>%
  # Exclude if no 9th event occurred or the 9th event was a treatment
  filter(!is.na(event9) & event9 != "X") %>%
  # Count each 10-event sequence
  group_by_at(vars(c(event1:event10))) %>% 
  summarise(count = sum(n)) %>%
  # Calculate probability
  ungroup() %>%
  group_by_at(vars(c(event1:event9))) %>%
  mutate(prob = prop.table(count)) %>%
  # If cascade discontinued, name 10th event as "0"
  mutate_at("event10", ~replace_na(.,"0")) %>%
  unite("pathString", c(event1:event10), na.rm = TRUE, sep = "/")

prob_data <- rbind(node1,node2,node3,node4,node5,node6,node7,node8,node9,node10) 

prob_data <- prob_data %>%
  add_row(pathString = paste0(node10$pathString[1],"/0"), prob = 1)

# Replace index letters with procedure categories. First replace with numbers to avoid double replacement.
prob_data$pathString <- stri_replace_all_regex(prob_data$pathString,
                                               pattern=c('A', 'B', 'C', 'D', 'E', 'X', 'N'),
                                               replacement=c('1', '2', '3', '4', '5', '6', '7'),
                                               vectorize=FALSE)
prob_data$pathString <- stri_replace_all_regex(prob_data$pathString,
                                               pattern=c('1', '2', '3', '4', '5', '6', '7', '0'),
                                               replacement=c('EKG', 'EP or cath', 'Imaging', 
                                                             'EST or monitor', 'XR', 'Tx', 'No EKG' , 'No tx'),
                                               vectorize=FALSE)

# Create root of tree 
prob_data <- prob_data %>%
  mutate(pathString = paste0("Endoscopy/", pathString))

## -----------------------------------------------------------------------------
## SECTION 4: CONVERT DATAFRAME TO TREE STRUCTURE
## -----------------------------------------------------------------------------
# Add terminal node of 'no treatment' for all paths that do not end in treatment
tree_df <- cascade %>%
  mutate(cascade_seq = ifelse(grepl("X$", cascade_seq), 
                              cascade_seq,
                              paste0(cascade_seq,"-0")))

# Replace index letters with procedure categories. First replace with numbers to avoid double replacement.
tree_df$cascade_seq <- stri_replace_all_regex(tree_df$cascade_seq,
                                              pattern=c('A', 'B', 'C', 'D', 'E', 'X', 'N'),
                                              replacement=c('1', '2', '3', '4', '5', '6', '7'),
                                              vectorize=FALSE)
tree_df$cascade_seq <- stri_replace_all_regex(tree_df$cascade_seq,
                                              pattern=c('1', '2', '3', '4', '5', '6', '7', '0'),
                                              replacement=c('EKG', 'EP or cath', 'Imaging', 
                                                            'EST or monitor', 'XR', 'Tx', 'No EKG' , 'No tx'),
                                              vectorize=FALSE)

# Create root of tree 
tree_df <- tree_df %>%
  mutate(cascade_seq = paste0("Endoscopy-", cascade_seq))

# Convert dataframe into tree 
tree <- as.Node(tree_df, pathName = "cascade_seq", pathDelimiter = "-")

print(tree,"pathString", "prob","cum_prob_1", "level")




## -----------------------------------------------------------------------------
## SECTION 5: ATTACH PROBABILITIES
## -----------------------------------------------------------------------------
# Attach probabilities as an attribute to each node
tree$Do(function(node) node$prob <- prob_data$prob[match(node$pathString, prob_data$pathString)])

# Set probabilities of level 1 nodes to one to make ECG vs. no ECG a decision node rather than a chance node
tree$'EKG'$prob <- 1
tree$'No EKG'$prob <- 1


## -----------------------------------------------------------------------------
## SECTION 6: ATTACH CUMULATIVE PROBABILITIES
## -----------------------------------------------------------------------------
tree$Do(function(node) node$cum_prob_1 <- node$parent$prob * node$prob, 
        traversal = "pre-order", filterFun = function(node) node$level == 3)
tree$Do(function(node) node$cum_prob <- node$cum_prob_1, 
        traversal = "pre-order", filterFun = function(node) node$level == 3 & node$count == 0)

tree$Do(function(node) node$cum_prob_2 <- node$parent$cum_prob_1 * node$prob, 
        traversal = "pre-order", filterFun = function(node) node$level == 4)
tree$Do(function(node) node$cum_prob <- node$cum_prob_2, 
        traversal = "pre-order", filterFun = function(node) node$level == 4 & node$count == 0)

tree$Do(function(node) node$cum_prob_3 <- node$parent$cum_prob_2 * node$prob, 
        traversal = "pre-order", filterFun = function(node) node$level == 5)
tree$Do(function(node) node$cum_prob <- node$cum_prob_3, 
        traversal = "pre-order", filterFun = function(node) node$level == 5 & node$count == 0)

tree$Do(function(node) node$cum_prob_4 <- node$parent$cum_prob_3 * node$prob, 
        traversal = "pre-order", filterFun = function(node) node$level == 6)
tree$Do(function(node) node$cum_prob <- node$cum_prob_4, 
        traversal = "pre-order", filterFun = function(node) node$level == 6 & node$count == 0)

tree$Do(function(node) node$cum_prob_5 <- node$parent$cum_prob_4 * node$prob, 
        traversal = "pre-order", filterFun = function(node) node$level == 7)
tree$Do(function(node) node$cum_prob <- node$cum_prob_5, 
        traversal = "pre-order", filterFun = function(node) node$level == 7 & node$count == 0)

tree$Do(function(node) node$cum_prob_6 <- node$parent$cum_prob_5 * node$prob, 
        traversal = "pre-order", filterFun = function(node) node$level == 8)
tree$Do(function(node) node$cum_prob <- node$cum_prob_6, 
        traversal = "pre-order", filterFun = function(node) node$level == 8 & node$count == 0)

tree$Do(function(node) node$cum_prob_7 <- node$parent$cum_prob_6 * node$prob, 
        traversal = "pre-order", filterFun = function(node) node$level == 9)
tree$Do(function(node) node$cum_prob <- node$cum_prob_7, 
        traversal = "pre-order", filterFun = function(node) node$level == 9 & node$count == 0)

tree$Do(function(node) node$cum_prob_8 <- node$parent$cum_prob_7 * node$prob, 
        traversal = "pre-order", filterFun = function(node) node$level == 10)
tree$Do(function(node) node$cum_prob <- node$cum_prob_8, 
        traversal = "pre-order", filterFun = function(node) node$level == 10 & node$count == 0)

tree$Do(function(node) node$cum_prob_9 <- node$parent$cum_prob_8 * node$prob, 
        traversal = "pre-order", filterFun = function(node) node$level == 11)
tree$Do(function(node) node$cum_prob <- node$cum_prob_9, 
        traversal = "pre-order", filterFun = function(node) node$level == 11 & node$count == 0)

tree$Do(function(node) node$cum_prob_10 <- node$parent$cum_prob_9 * node$prob, 
        traversal = "pre-order", filterFun = function(node) node$level == 12)
tree$Do(function(node) node$cum_prob <- node$cum_prob_10, 
        traversal = "pre-order", filterFun = function(node) node$level == 12 & node$count == 0)


## -----------------------------------------------------------------------------
## SECTION 7: ATTACH COSTS
## -----------------------------------------------------------------------------

# Calculate gamma based on mean and SD of each test category
gamma_parms <- function(mean, sd){
  beta <- sd^2/mean
  alpha <- mean/beta
  return(list(alpha = alpha, beta = beta))
}

sd <- apply(costs[,-1], 1, sd, na.rm = TRUE)  
costs <- cbind(costs, sd = sqrt(sd)) 

costs[,c(12,13)] <- gamma_parms(costs$mean, costs$sd)

# Sample from gamma distribution
ecg <- rgamma(n = 1000, shape = costs$alpha[1], scale = costs$beta[1])
ep_cath <- rgamma(n = 1000, shape = costs$alpha[2], scale = costs$beta[2])
imaging <- rgamma(n = 1000, shape = costs$alpha[3], scale = costs$beta[3])
est_monitor <- rgamma(n = 1000, shape = costs$alpha[4], scale = costs$beta[4])
xr <- rgamma(n = 1000, shape = costs$alpha[5], scale = costs$beta[5])



## -----------------------------------------------------------------------------
## SECTION 6: CUMULATIVE COSTS
## -----------------------------------------------------------------------------

# Calculate gamma based on mean and SD of each test category
gamma_parms <- function(mean, sd){
  beta <- sd^2/mean
  alpha <- mean/beta
  return(list(alpha = alpha, beta = beta))
}

sd <- apply(costs[,-1], 1, sd, na.rm = TRUE)  
costs <- cbind(costs, sd = sqrt(sd)) 

costs[,c(12,13)] <- gamma_parms(costs$mean, costs$sd)

# Sample from gamma distribution
ecg <- rgamma(n = 1000, shape = costs$alpha[1], scale = costs$beta[1])
ep_cath <- rgamma(n = 1000, shape = costs$alpha[2], scale = costs$beta[2])
imaging <- rgamma(n = 1000, shape = costs$alpha[3], scale = costs$beta[3])
est_monitor <- rgamma(n = 1000, shape = costs$alpha[4], scale = costs$beta[4])
xr <- rgamma(n = 1000, shape = costs$alpha[5], scale = costs$beta[5])


# Attach costs from 1000 simulations
sim_cumCost <- ToDataFrameTree(tree,"pathString")[2]

for (i in 1:1000){
  
  attach_costs <- function(node){
    if (node$name == "ecg"){node$cost <-  ecg[i]}
    if (node$name == "EP or cath"){node$cost <- ep_cath[i]}
    if (node$name == "Imaging"){node$cost <- imaging[i]}
    if (node$name == "EST or monitor"){node$cost <- est_monitor[i]}
    if (node$name == "XR"){node$cost <- xr[i]}
    if (node$name == "No tx" | node$name == "Tx"){node$cost <- 0}
  }
  
  
  tree$Do(attach_costs)
  tree$'No ecg'$cost <- 0
  
  tree$Do(function(node) node$cum_cost_1 <- node$parent$cost + node$cost, 
          traversal = "pre-order", filterFun = function(node) node$level == 3)
  tree$Do(function(node) if (node$name != "Tx") node$cum_cost <- node$cum_cost_1 else node$cum_cost <- 0, 
          traversal = "pre-order", filterFun = function(node) node$level == 3 & node$count == 0)
  
  tree$Do(function(node) node$cum_cost_2 <- node$parent$cum_cost_1 + node$cost, 
          traversal = "pre-order", filterFun = function(node) node$level == 4)  
  tree$Do(function(node) if (node$name != "Tx") node$cum_cost <- node$cum_cost_2 else node$cum_cost <- 0, 
          traversal = "pre-order", filterFun = function(node) node$level == 4 & node$count == 0)
  
  tree$Do(function(node) node$cum_cost_3 <- node$parent$cum_cost_2 + node$cost, 
          traversal = "pre-order", filterFun = function(node) node$level == 5)  
  tree$Do(function(node) if (node$name != "Tx") node$cum_cost <- node$cum_cost_3 else node$cum_cost <- 0, 
          traversal = "pre-order", filterFun = function(node) node$level == 5 & node$count == 0)
  
  tree$Do(function(node) node$cum_cost_4 <- node$parent$cum_cost_3 + node$cost, 
          traversal = "pre-order", filterFun = function(node) node$level == 6)  
  tree$Do(function(node) if (node$name != "Tx") node$cum_cost <- node$cum_cost_4 else node$cum_cost <- 0, 
          traversal = "pre-order", filterFun = function(node) node$level == 6 & node$count == 0)
  
  tree$Do(function(node) node$cum_cost_5 <- node$parent$cum_cost_4 + node$cost, 
          traversal = "pre-order", filterFun = function(node) node$level == 7)  
  tree$Do(function(node) if (node$name != "Tx") node$cum_cost <- node$cum_cost_5 else node$cum_cost <- 0, 
          traversal = "pre-order", filterFun = function(node) node$level == 7 & node$count == 0)
  
  tree$Do(function(node) node$cum_cost_6 <- node$parent$cum_cost_5 + node$cost, 
          traversal = "pre-order", filterFun = function(node) node$level == 8)  
  tree$Do(function(node) if (node$name != "Tx") node$cum_cost <- node$cum_cost_6 else node$cum_cost <- 0, 
          traversal = "pre-order", filterFun = function(node) node$level == 8 & node$count == 0)
  
  tree$Do(function(node) node$cum_cost_7 <- node$parent$cum_cost_6 + node$cost, 
          traversal = "pre-order", filterFun = function(node) node$level == 9)  
  tree$Do(function(node) if (node$name != "Tx") node$cum_cost <- node$cum_cost_7 else node$cum_cost <- 0, 
          traversal = "pre-order", filterFun = function(node) node$level == 9 & node$count == 0)
  
  tree$Do(function(node) node$cum_cost_8 <- node$parent$cum_cost_7 + node$cost, 
          traversal = "pre-order", filterFun = function(node) node$level == 10)  
  tree$Do(function(node) if (node$name != "Tx") node$cum_cost <- node$cum_cost_8 else node$cum_cost <- 0, 
          traversal = "pre-order", filterFun = function(node) node$level == 10 & node$count == 0)
  
  tree$Do(function(node) node$cum_cost_9 <- node$parent$cum_cost_8 + node$cost, 
          traversal = "pre-order", filterFun = function(node) node$level == 11)  
  tree$Do(function(node) if (node$name != "Tx") node$cum_cost <- node$cum_cost_9 else node$cum_cost <- 0, 
          traversal = "pre-order", filterFun = function(node) node$level == 11 & node$count == 0)
  
  tree$Do(function(node) node$cum_cost_10 <- node$parent$cum_cost_9 + node$cost, 
          traversal = "pre-order", filterFun = function(node) node$level == 12)  
  tree$Do(function(node) if (node$name != "Tx") node$cum_cost <- node$cum_cost_10 else node$cum_cost <- 0, 
          traversal = "pre-order", filterFun = function(node) node$level == 12 & node$count == 0)
  
  tmp_df <- ToDataFrameTree(tree,"pathString","cum_cost")
  #  tmp_df <- tmp_df[-1,]
  names(tmp_df)[3] <- paste0("V",i)
  tmp_df <- tmp_df[match(tmp_df$pathString, sim_cumCost$pathString),]
  sim_cumCost <- cbind(sim_cumCost, tmp_df[3])
}


## -----------------------------------------------------------------------------
## SECTION 7: ESTIMATE EXPECTED COSTS AT TERMINAL NODES
## -----------------------------------------------------------------------------

sim_cumCost <- na.omit(sim_cumCost)
sim_cumProb <- na.omit(sim_cumProb)

# Set node type as terminal or chance
tree$Do(function(node) if(node$count == 0) node$type <- "terminal" else node$type <- "chance")

# Calculate expected cost given cumulative probabilities and cumulative costs at terminal nodes
expected_cost <- ToDataFrameTree(tree,"pathString")[2]

for (i in 1:1000){
  tree$Do(function(node) node$cum_prob <- sim_cumProb[match(node$pathString, sim_cumProb$pathString), i+1])
  tree$Do(function(node) node$cum_cost <- sim_cumCost[match(node$pathString, sim_cumCost$pathString), i+1])
  
  n_lvc <- 59121*node1_sim[1,i+1]
  
  estimate_payoff <- function(node){
    node$expected_cost <- node$cum_prob * node$cum_cost * n_lvc
  }
  
  tree$Do(estimate_payoff, filterFun = function(node) node$type == "terminal")
  
  tmp_df <- ToDataFrameTree(tree,"pathString","expected_cost")
  
  names(tmp_df)[3] <- paste0("V",i)
  tmp_df <- tmp_df[match(tmp_df$pathString, expected_cost$pathString),]
  expected_cost <- cbind(expected_cost, tmp_df[3])
  
}

# Export results
psa_results <- expected_cost %>%  
  mutate(pathString = str_remove(pathString , "Endoscopy/")) %>%
  filter(!is.na(V1)) %>%
  mutate(group = sub("\\/.*", "",pathString)) %>%  
  group_by(group) %>% 
  summarise(across(where(is.numeric), sum))

psa_results <- aggregate(.~group, psa_results, sum)

# Summary statistics for cascade costs from ECG group 
ecg_mean <- mean(as.numeric(psa_results[1,-1]))
ecg_ci <- quantile(as.numeric(psa_results[1,-1]), probs = c(0.025,0.975))
ecg_median <- median(as.numeric(psa_results[1,-1]))
ecg_min <- min(as.numeric(psa_results[1,-1]))
ecg_max <- max(as.numeric(psa_results[1,-1]))

# Summary statistics for cascade costs from comparison group
comp_mean <- mean(as.numeric(psa_results[2,-1]))
comp_ci <- quantile(as.numeric(psa_results[2,-1]), probs = c(0.025,0.975))
comp_median <- median(as.numeric(psa_results[2,-1]))
comp_min <- min(as.numeric(psa_results[2,-1]))
comp_max <- max(as.numeric(psa_results[2,-1]))

# Summary statistics for excess costs
excess_cost <- as.data.frame(diff(as.matrix(psa_results[c(2:1),-1])))
excess_cost_mean <- mean(as.numeric(excess_cost[1,]))
excess_cost_ci <- quantile(as.numeric(excess_cost), probs = c(0.025,0.975))
excess_cost_median <- median(as.numeric(excess_cost))
excess_cost_min <- min(as.numeric(excess_cost))
excess_cost_max <- max(as.numeric(excess_cost))


## -----------------------------------------------------------------------------
## SECTION 8: PROJECTED COSTS
## -----------------------------------------------------------------------------

# Calculate +/- 25% of projected prevalence
proj <- proj %>%
  mutate(lower = 0.75 * total, upper = 1.25 * total)

# Sample from uniform distribution
n_2024 <- runif(n = 1000, min = proj$lower[1], max = proj$upper[1])
n_2025 <- runif(n = 1000, min = proj$lower[2], max = proj$upper[2])
n_2026 <- runif(n = 1000, min = proj$lower[3], max = proj$upper[3])
n_2027 <- runif(n = 1000, min = proj$lower[4], max = proj$upper[4])
n_2028 <- runif(n = 1000, min = proj$lower[5], max = proj$upper[5])

projected_cost_2024 <- ToDataFrameTree(tree,"pathString")[2]
projected_cost_2025 <- ToDataFrameTree(tree,"pathString")[2]
projected_cost_2026 <- ToDataFrameTree(tree,"pathString")[2]
projected_cost_2027 <- ToDataFrameTree(tree,"pathString")[2]
projected_cost_2028 <- ToDataFrameTree(tree,"pathString")[2]

# Project costs
for (i in 1:1000){
  tree$Do(function(node) node$cum_prob <- sim_cumProb[match(node$pathString, sim_cumProb$pathString), i+1])
  tree$Do(function(node) node$cum_cost <- sim_cumCost[match(node$pathString, sim_cumCost$pathString), i+1])
  
  n_lvc_2024 <- n_2024[i]*node1_sim[1,i+1]
  n_lvc_2025 <- n_2025[i]*node1_sim[1,i+1]
  n_lvc_2026 <- n_2026[i]*node1_sim[1,i+1]
  n_lvc_2027 <- n_2027[i]*node1_sim[1,i+1]
  n_lvc_2028 <- n_2028[i]*node1_sim[1,i+1]
  
  estimate_payoff <- function(node){
    node$expected_cost_2024 <- node$cum_prob * node$cum_cost * n_lvc_2024
    node$expected_cost_2025 <- node$cum_prob * node$cum_cost * n_lvc_2025
    node$expected_cost_2026 <- node$cum_prob * node$cum_cost * n_lvc_2026
    node$expected_cost_2027 <- node$cum_prob * node$cum_cost * n_lvc_2027
    node$expected_cost_2028 <- node$cum_prob * node$cum_cost * n_lvc_2028
  }
  
  tree$Do(estimate_payoff, filterFun = function(node) node$type == "terminal")
  
  tmp_df_2024 <- ToDataFrameTree(tree,"pathString","expected_cost_2024")
  names(tmp_df_2024)[3] <- paste0("V",i)
  tmp_df_2024 <- tmp_df_2024[match(tmp_df_2024$pathString, projected_cost_2024$pathString),]
  projected_cost_2024 <- cbind(projected_cost_2024, tmp_df_2024[3])
  
  tmp_df_2025 <- ToDataFrameTree(tree,"pathString","expected_cost_2025")
  names(tmp_df_2025)[3] <- paste0("V",i)
  tmp_df_2025 <- tmp_df_2025[match(tmp_df_2025$pathString, projected_cost_2025$pathString),]
  projected_cost_2025 <- cbind(projected_cost_2025, tmp_df_2025[3])
  
  tmp_df_2026 <- ToDataFrameTree(tree,"pathString","expected_cost_2026")
  names(tmp_df_2026)[3] <- paste0("V",i)
  tmp_df_2026 <- tmp_df_2026[match(tmp_df_2026$pathString, projected_cost_2026$pathString),]
  projected_cost_2026 <- cbind(projected_cost_2026, tmp_df_2026[3])
  
  tmp_df_2027 <- ToDataFrameTree(tree,"pathString","expected_cost_2027")
  names(tmp_df_2027)[3] <- paste0("V",i)
  tmp_df_2027 <- tmp_df_2027[match(tmp_df_2027$pathString, projected_cost_2027$pathString),]
  projected_cost_2027 <- cbind(projected_cost_2027, tmp_df_2027[3])
  
  tmp_df_2028 <- ToDataFrameTree(tree,"pathString","expected_cost_2028")
  names(tmp_df_2028)[3] <- paste0("V",i)
  tmp_df_2028 <- tmp_df_2028[match(tmp_df_2028$pathString, projected_cost_2028$pathString),]
  projected_cost_2028 <- cbind(projected_cost_2028, tmp_df_2028[3])
}

# Projected cost for 2024 ####
results_2024 <- projected_cost_2024 %>%  
  mutate(pathString = str_remove(pathString , "Endoscopy/")) %>%
  filter(!is.na(V1)) %>%
  mutate(group = sub("\\/.*", "",pathString)) %>%  
  group_by(group) %>% 
  summarise(across(where(is.numeric), sum))

results_2024 <- aggregate(.~group, results_2024, sum)

ecg_mean_2024 <- mean(as.numeric(results_2024[1,-1]))
ecg_ci_2024 <- quantile(as.numeric(results_2024[1,-1]), probs = c(0.025,0.975))
ecg_median_2024 <- median(as.numeric(results_2024[1,-1]))
ecg_min_2024 <- min(as.numeric(results_2024[1,-1]))
ecg_max_2024 <- max(as.numeric(results_2024[1,-1]))

comp_mean_2024 <- mean(as.numeric(results_2024[2,-1]))
comp_ci_2024 <- quantile(as.numeric(results_2024[2,-1]), probs = c(0.025,0.975))
comp_median_2024 <- median(as.numeric(results_2024[2,-1]))
comp_min_2024 <- min(as.numeric(results_2024[2,-1]))
comp_max_2024 <- max(as.numeric(results_2024[2,-1]))

excess_cost_2024 <- as.data.frame(diff(as.matrix(results_2024[c(2:1),-1])))
excess_cost_mean_2024 <- mean(as.numeric(excess_cost_2024[1,]))
excess_cost_ci_2024 <- quantile(as.numeric(excess_cost_2024), probs = c(0.025,0.975))
excess_cost_median_2024 <- median(as.numeric(excess_cost_2024))
excess_cost_min_2024 <- min(as.numeric(excess_cost_2024))
excess_cost_max_2024 <- max(as.numeric(excess_cost_2024))

# Projected cost for 2025 ####
results_2025 <- projected_cost_2025 %>%  
  mutate(pathString = str_remove(pathString , "Endoscopy/")) %>%
  filter(!is.na(V1)) %>%
  mutate(group = sub("\\/.*", "",pathString)) %>%  
  group_by(group) %>% 
  summarise(across(where(is.numeric), sum))

results_2025 <- aggregate(.~group, results_2025, sum)


ecg_mean_2025 <- mean(as.numeric(results_2025[1,-1]))
ecg_ci_2025 <- quantile(as.numeric(results_2025[1,-1]), probs = c(0.025,0.975))
ecg_median_2025 <- median(as.numeric(results_2025[1,-1]))
ecg_min_2025 <- min(as.numeric(results_2025[1,-1]))
ecg_max_2025 <- max(as.numeric(results_2025[1,-1]))

comp_mean_2025 <- mean(as.numeric(results_2025[2,-1]))
comp_ci_2025 <- quantile(as.numeric(results_2025[2,-1]), probs = c(0.025,0.975))
comp_median_2025 <- median(as.numeric(results_2025[2,-1]))
comp_min_2025 <- min(as.numeric(results_2025[2,-1]))
comp_max_2025 <- max(as.numeric(results_2025[2,-1]))

excess_cost_2025 <- as.data.frame(diff(as.matrix(results_2025[c(2:1),-1])))
excess_cost_mean_2025 <- mean(as.numeric(excess_cost_2025[1,]))
excess_cost_ci_2025 <- quantile(as.numeric(excess_cost_2025), probs = c(0.025,0.975))
excess_cost_median_2025 <- median(as.numeric(excess_cost_2025))
excess_cost_min_2025 <- min(as.numeric(excess_cost_2025))
excess_cost_max_2025 <- max(as.numeric(excess_cost_2025))


# Projected cost for 2026 ####
results_2026 <- projected_cost_2026 %>%  
  mutate(pathString = str_remove(pathString , "Endoscopy/")) %>%
  filter(!is.na(V1)) %>%
  mutate(group = sub("\\/.*", "",pathString)) %>%  
  group_by(group) %>% 
  summarise(across(where(is.numeric), sum))

results_2026 <- aggregate(.~group, results_2026, sum)

ecg_mean_2026 <- mean(as.numeric(results_2026[1,-1]))
ecg_ci_2026 <- quantile(as.numeric(results_2026[1,-1]), probs = c(0.025,0.975))
ecg_median_2026 <- median(as.numeric(results_2026[1,-1]))
ecg_min_2026 <- min(as.numeric(results_2026[1,-1]))
ecg_max_2026 <- max(as.numeric(results_2026[1,-1]))

comp_mean_2026 <- mean(as.numeric(results_2026[2,-1]))
comp_ci_2026 <- quantile(as.numeric(results_2026[2,-1]), probs = c(0.025,0.975))
comp_median_2026 <- median(as.numeric(results_2026[2,-1]))
comp_min_2026 <- min(as.numeric(results_2026[2,-1]))
comp_max_2026 <- max(as.numeric(results_2026[2,-1]))

excess_cost_2026 <- as.data.frame(diff(as.matrix(results_2026[c(2:1),-1])))
excess_cost_mean_2026 <- mean(as.numeric(excess_cost_2026[1,]))
excess_cost_ci_2026 <- quantile(as.numeric(excess_cost_2026), probs = c(0.025,0.975))
excess_cost_median_2026 <- median(as.numeric(excess_cost_2026))
excess_cost_min_2026 <- min(as.numeric(excess_cost_2026))
excess_cost_max_2026 <- max(as.numeric(excess_cost_2026))


# Projected cost for 2027 ####
results_2027 <- projected_cost_2027 %>%  
  mutate(pathString = str_remove(pathString , "Endoscopy/")) %>%
  filter(!is.na(V1)) %>%
  mutate(group = sub("\\/.*", "",pathString)) %>%  
  group_by(group) %>% 
  summarise(across(where(is.numeric), sum))

results_2027 <- aggregate(.~group, results_2027, sum)

ecg_mean_2027 <- mean(as.numeric(results_2027[1,-1]))
ecg_ci_2027 <- quantile(as.numeric(results_2027[1,-1]), probs = c(0.025,0.975))
ecg_median_2027 <- median(as.numeric(results_2027[1,-1]))
ecg_min_2027 <- min(as.numeric(results_2027[1,-1]))
ecg_max_2027 <- max(as.numeric(results_2027[1,-1]))

comp_mean_2027 <- mean(as.numeric(results_2027[2,-1]))
comp_ci_2027 <- quantile(as.numeric(results_2027[2,-1]), probs = c(0.025,0.975))
comp_median_2027 <- median(as.numeric(results_2027[2,-1]))
comp_min_2027 <- min(as.numeric(results_2027[2,-1]))
comp_max_2027 <- max(as.numeric(results_2027[2,-1]))

excess_cost_2027 <- as.data.frame(diff(as.matrix(results_2027[c(2:1),-1])))
excess_cost_mean_2027 <- mean(as.numeric(excess_cost_2027[1,]))
excess_cost_ci_2027 <- quantile(as.numeric(excess_cost_2027), probs = c(0.025,0.975))
excess_cost_median_2027 <- median(as.numeric(excess_cost_2027))
excess_cost_min_2027 <- min(as.numeric(excess_cost_2027))
excess_cost_max_2027 <- max(as.numeric(excess_cost_2027))


# Projected cost for 2028 ####
results_2028 <- projected_cost_2028 %>%  
  mutate(pathString = str_remove(pathString , "Endoscopy/")) %>%
  filter(!is.na(V1)) %>%
  mutate(group = sub("\\/.*", "",pathString)) %>%  
  group_by(group) %>% 
  summarise(across(where(is.numeric), sum))

results_2028 <- aggregate(.~group, results_2028, sum)

ecg_mean_2028 <- mean(as.numeric(results_2028[1,-1]))
ecg_ci_2028 <- quantile(as.numeric(results_2028[1,-1]), probs = c(0.025,0.975))
ecg_median_2028 <- median(as.numeric(results_2028[1,-1]))
ecg_min_2028 <- min(as.numeric(results_2028[1,-1]))
ecg_max_2028 <- max(as.numeric(results_2028[1,-1]))

comp_mean_2028 <- mean(as.numeric(results_2028[2,-1]))
comp_ci_2028 <- quantile(as.numeric(results_2028[2,-1]), probs = c(0.025,0.975))
comp_median_2028 <- median(as.numeric(results_2028[2,-1]))
comp_min_2028 <- min(as.numeric(results_2028[2,-1]))
comp_max_2028 <- max(as.numeric(results_2028[2,-1]))

excess_cost_2028 <- as.data.frame(diff(as.matrix(results_2028[c(2:1),-1])))
excess_cost_mean_2028 <- mean(as.numeric(excess_cost_2028[1,]))
excess_cost_ci_2028 <- quantile(as.numeric(excess_cost_2028), probs = c(0.025,0.975))
excess_cost_median_2028 <- median(as.numeric(excess_cost_2028))
excess_cost_min_2028 <- min(as.numeric(excess_cost_2028))
excess_cost_max_2028 <- max(as.numeric(excess_cost_2028))


# Collate costs from 2024-2028 ####
excess_cost_2024_2028 <- data.frame(year = seq(2024,2028,1),
                                    median = c(excess_cost_median_2024,
                                               excess_cost_median_2025,
                                               excess_cost_median_2026,
                                               excess_cost_median_2027,
                                               excess_cost_median_2028),
                                    min = c(excess_cost_min_2024,
                                            excess_cost_min_2025,
                                            excess_cost_min_2026,
                                            excess_cost_min_2027,
                                            excess_cost_min_2028),
                                    max = c(excess_cost_max_2024,
                                            excess_cost_max_2025,
                                            excess_cost_max_2026,
                                            excess_cost_max_2027,
                                            excess_cost_max_2028))


# Convert to dataframe and format for plotting
excess_cost_2024_2028 <- excess_cost_2024_2028 %>%
  mutate(time = year - 2022,
         discount_rate = 0.03,
         pv_median = median / (1+discount_rate)^time,
         pv_min = min / (1+discount_rate)^time,
         pv_max = max / (1+discount_rate)^time)

# Sum projected excess cost over 5 years 
excess_cost_2024_2028 %>%
  summarise(mean = sum(pv_median), lower = sum(pv_min), upper = sum(pv_max))

ggplot(excess_cost_2024_2028) +
  geom_ribbon(aes(x = year, ymin = min, ymax = max), 
              fill = "#E8464E", alpha = 0.2) +
  geom_line(aes(x = year, y = median), col = "#E8464E", lty = 1, linewidth = 0.5) +
  geom_point(aes(x = year, y = median), col = "#E8464E", size = 1) +
  scale_x_continuous(expand = c(0.01,0.01)) + 
  scale_y_continuous(limits = c(0,5e5), breaks = seq(0,5e5,1e5),
                     expand = c(0,0), labels = scales::comma) +
  labs(x = "Year", y = "Projected excess cost, USD") +
  theme_mgh() +
  theme(axis.text.y = element_text(margin = margin(t = 0, r = 8, b = 0, l = 0)),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), axis.title = element_text(size = 10))

