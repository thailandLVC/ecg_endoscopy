################################################################################
## Decision tree model for baseline estimates
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
cascadeE <- read.csv("~/cascade_seq_ECG.csv")
cascadeC <- read.csv("~/cascade_seq_comp.csv")

# Import cost data
costs <- read.csv("~/cascade_costs.csv")

# Import projected prevalence from ARIMA model
proj <- read.csv("~/projected_prev.csv")

## -----------------------------------------------------------------------------
## SECTION 2: DATA PREPARATION
## -----------------------------------------------------------------------------
# Count number of each group that was admitted during the cascade period
n_ECG <- sum(cascadeE$n)
n_noECG <- sum(cascadeC$n)

# Add 'No ECG' as the first event of cascade sequence for comparison group
cascadeC <- cascadeC %>%
  mutate(cascade_seq = paste0("N-",cascade_seq)) %>%
  # Add row for patients who did receive an ECG and were not readmitted during the cascade period
  add_row(cascade_seq = "N", n = 55443 - n_noECG)

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


## -----------------------------------------------------------------------------
## SECTION 3: CALCULATE NODE PROBABILITIES
## -----------------------------------------------------------------------------

# Node 1: ECG vs. no ECG
node1 <- cascadeWide %>%
  # Count ECG versus no ECG
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
                                               replacement=c('ECG', 'EP or cath', 'Imaging', 
                                                             'EST or monitor', 'XR', 'Tx', 'No ECG' , 'No tx'),
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
                                              replacement=c('ECG', 'EP or cath', 'Imaging', 
                                                            'EST or monitor', 'XR', 'Tx', 'No ECG' , 'No tx'),
                                              vectorize=FALSE)

# Create root of tree 
tree_df <- tree_df %>%
  mutate(cascade_seq = paste0("Endoscopy-", cascade_seq))

# Convert dataframe into tree 
tree <- as.Node(tree_df, pathName = "cascade_seq", pathDelimiter = "-")


## -----------------------------------------------------------------------------
## SECTION 5: POPULATE TREE WITH PROBABILITIES
## -----------------------------------------------------------------------------

# Attach probabilities as an attribute to each node
tree$Do(function(node) node$prob <- prob_data$prob[match(node$pathString, prob_data$pathString)])

# Set probabilities of level 1 nodes to one to make ECG vs. no ECG a decision node rather than a chance node
tree$'ECG'$prob <- 1
tree$'No ECG'$prob <- 1


## -----------------------------------------------------------------------------
## SECTION 6: ATTACH COSTS
## -----------------------------------------------------------------------------

# Create function for extracting mean cost of each test category
attach_costs <- function(node){
  if (node$name == "ECG"){node$cost_mean <- costs$mean[costs$test == "ECG"]}
  if (node$name == "EP or cath"){node$cost_mean <- costs$mean[costs$test == "EP or cath"]}
  if (node$name == "Imaging"){node$cost_mean <- costs$mean[costs$test == "Imaging"]}
  if (node$name == "EST or monitor"){node$cost_mean <- costs$mean[costs$test == "EST or monitor"]}
  if (node$name == "XR"){node$cost_mean <- costs$mean[costs$test == "XR"]}
  if (node$name == "No tx" | node$name == "Tx"){node$cost_mean <- 0}
}

# Attach costs
tree$Do(attach_costs)

# Set cost of no ECG to 0 
tree$'No ECG'$cost_mean <- 0


## -----------------------------------------------------------------------------
## SECTION 7: CALCULATE CUMULATIVE PROBABILITIES AT TERMINAL NODES
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
## SECTION 8: CALCULATE CUMULATIVE COSTS AT TERMINAL NODES
## -----------------------------------------------------------------------------

tree$Do(function(node) node$cum_cost_1 <- node$parent$cost_mean + node$cost_mean, 
        traversal = "pre-order", filterFun = function(node) node$level == 3)
tree$Do(function(node) if (node$name != "Tx") node$cum_cost <- node$cum_cost_1 else node$cum_cost <- 0, 
        traversal = "pre-order", filterFun = function(node) node$level == 3 & node$count == 0)

tree$Do(function(node) node$cum_cost_2 <- node$parent$cum_cost_1 + node$cost_mean, 
        traversal = "pre-order", filterFun = function(node) node$level == 4)  
tree$Do(function(node) if (node$name != "Tx") node$cum_cost <- node$cum_cost_2 else node$cum_cost <- 0, 
        traversal = "pre-order", filterFun = function(node) node$level == 4 & node$count == 0)

tree$Do(function(node) node$cum_cost_3 <- node$parent$cum_cost_2 + node$cost_mean, 
        traversal = "pre-order", filterFun = function(node) node$level == 5)  
tree$Do(function(node) if (node$name != "Tx") node$cum_cost <- node$cum_cost_3 else node$cum_cost <- 0, 
        traversal = "pre-order", filterFun = function(node) node$level == 5 & node$count == 0)

tree$Do(function(node) node$cum_cost_4 <- node$parent$cum_cost_3 + node$cost_mean, 
        traversal = "pre-order", filterFun = function(node) node$level == 6)  
tree$Do(function(node) if (node$name != "Tx") node$cum_cost <- node$cum_cost_4 else node$cum_cost <- 0, 
        traversal = "pre-order", filterFun = function(node) node$level == 6 & node$count == 0)

tree$Do(function(node) node$cum_cost_5 <- node$parent$cum_cost_4 + node$cost_mean, 
        traversal = "pre-order", filterFun = function(node) node$level == 7)  
tree$Do(function(node) if (node$name != "Tx") node$cum_cost <- node$cum_cost_5 else node$cum_cost <- 0, 
        traversal = "pre-order", filterFun = function(node) node$level == 7 & node$count == 0)

tree$Do(function(node) node$cum_cost_6 <- node$parent$cum_cost_5 + node$cost_mean, 
        traversal = "pre-order", filterFun = function(node) node$level == 8)  
tree$Do(function(node) if (node$name != "Tx") node$cum_cost <- node$cum_cost_6 else node$cum_cost <- 0, 
        traversal = "pre-order", filterFun = function(node) node$level == 8 & node$count == 0)

tree$Do(function(node) node$cum_cost_7 <- node$parent$cum_cost_6 + node$cost_mean, 
        traversal = "pre-order", filterFun = function(node) node$level == 9)  
tree$Do(function(node) if (node$name != "Tx") node$cum_cost <- node$cum_cost_7 else node$cum_cost <- 0, 
        traversal = "pre-order", filterFun = function(node) node$level == 9 & node$count == 0)

tree$Do(function(node) node$cum_cost_8 <- node$parent$cum_cost_7 + node$cost_mean, 
        traversal = "pre-order", filterFun = function(node) node$level == 10)  
tree$Do(function(node) if (node$name != "Tx") node$cum_cost <- node$cum_cost_8 else node$cum_cost <- 0, 
        traversal = "pre-order", filterFun = function(node) node$level == 10 & node$count == 0)

tree$Do(function(node) node$cum_cost_9 <- node$parent$cum_cost_8 + node$cost_mean, 
        traversal = "pre-order", filterFun = function(node) node$level == 11)  
tree$Do(function(node) if (node$name != "Tx") node$cum_cost <- node$cum_cost_9 else node$cum_cost <- 0, 
        traversal = "pre-order", filterFun = function(node) node$level == 11 & node$count == 0)

tree$Do(function(node) node$cum_cost_10 <- node$parent$cum_cost_9 + node$cost_mean, 
        traversal = "pre-order", filterFun = function(node) node$level == 12)  
tree$Do(function(node) if (node$name != "Tx") node$cum_cost <- node$cum_cost_10 else node$cum_cost <- 0, 
        traversal = "pre-order", filterFun = function(node) node$level == 12 & node$count == 0)


## -----------------------------------------------------------------------------
## SECTION 9: PLOT TREE
## -----------------------------------------------------------------------------

# Set node type as terminal or chance
tree$Do(function(node) if(node$count == 0) node$type <- "terminal" else node$type <- "chance")

tree$'ECG'$prob <- NA
tree$'No ECG'$prob <- NA

get_edge_label <- function(node) { 
  if(!is.na(node$prob)){
    edge_label <- paste0("p = ", round(as.numeric(node$prob),3))
    return(edge_label)
  }
}

get_node_label <- function(node) {
  if(node$type == "terminal" & node$name != "Tx"){
    node_label <- paste0("No tx \n$", as.numeric(node$cum_cost), " (",
                         round(100*as.numeric(node$cum_prob), 4), "%)")
  }
}


get_font_color <- function(node) 
  switch(node$name, `Endoscopy` = "black", `ECG` = "black", `No ECG` = "black", `EP or cath` = "black", 
         `Imaging` = "black", `EST or monitor` = "black", `XR` = "black", `Tx` =  "black", `No tx` = "#A4001E")

get_node_color <- function(node) 
  switch(node$name, `ECG` = "white", `EP or cath` = "white", `Imaging` = "white", 
         `EST or monitor` = "white", `XR` = "white", `Tx` =  "white", `No tx` = "#A4001E")

SetEdgeStyle(tree, fontname = 'helvetica', label = get_edge_label, fontsize = 80)
SetNodeStyle(tree, fontname = 'helvetica', label = get_node_label, color = get_node_color, penwidth = 5,
             shape = "box", fontsize = 100, fontcolor = get_font_color)

SetGraphStyle(tree, rankdir = "LR")
plot(tree)


## -----------------------------------------------------------------------------
## SECTION 10: ESTIMATE COST OF CASCADES, FY2018-19
## -----------------------------------------------------------------------------

# Size of ECG group
n_lvc <- 59121*node1$prob[1]

# Multiply cumulative costs by cumulative probabilities
estimate_payoff <- function(node){
  node$expected_cost <- node$cum_prob * node$cum_cost * n_lvc
  node$expected_minCost <- node$cum_prob * node$cum_minCost * n_lvc
  node$expected_maxCost <- node$cum_prob * node$cum_maxCost * n_lvc
}

tree$Do(estimate_payoff, filterFun = function(node) node$type == "terminal")

tree_results_df <- ToDataFrameTree(tree,"pathString","expected_cost","expected_minCost","expected_maxCost")

# Calculate excess cost
tree_results_df <- tree_results_df %>%  
  select(-levelName) %>%
  mutate(pathString = str_remove(pathString , "Endoscopy/")) %>%
  filter(!is.na(expected_cost) & pathString != "Endoscopy") %>%
  mutate(group = sub("\\/.*", "",pathString)) %>%  
  group_by(group) %>% 
  summarise(total_cost = sum(expected_cost))

net_costs <- diff(as.matrix(tree_results_df[c(2,1),-1]))



## -----------------------------------------------------------------------------
## SECTION 11: PROJECT COSTS FOR FY2024-28
## -----------------------------------------------------------------------------

# Size of projected low-value ECG group
n_lvc_2024 <- proj$total[1]*node1$prob[1]
n_lvc_2025 <- proj$total[2]*node1$prob[1]
n_lvc_2026 <- proj$total[3]*node1$prob[1]
n_lvc_2027 <- proj$total[4]*node1$prob[1]
n_lvc_2028 <- proj$total[5]*node1$prob[1]

# Multiply cumulative costs by cumulative probabilities
project_payoff <- function(node){
  node$expected_cost_2024 <- node$cum_prob * node$cum_cost * n_lvc_2024
  node$expected_cost_2025 <- node$cum_prob * node$cum_cost * n_lvc_2025
  node$expected_cost_2026 <- node$cum_prob * node$cum_cost * n_lvc_2026
  node$expected_cost_2027 <- node$cum_prob * node$cum_cost * n_lvc_2027
  node$expected_cost_2028 <- node$cum_prob * node$cum_cost * n_lvc_2028
}

tree$Do(project_payoff, filterFun = function(node) node$type == "terminal")

# Store projection results in a dataframe
tree_proj_df <- ToDataFrameTree(tree,"pathString",
                                "expected_cost_2024","expected_minCost_2024","expected_maxCost_2024",
                                "expected_cost_2025","expected_minCost_2025","expected_maxCost_2025",
                                "expected_cost_2026","expected_minCost_2026","expected_maxCost_2026",
                                "expected_cost_2027","expected_minCost_2027","expected_maxCost_2027",
                                "expected_cost_2028","expected_minCost_2028","expected_maxCost_2028")

# Sum costs by year
tree_proj_df <- tree_proj_df %>%  
  select(-levelName) %>%
  mutate(pathString = str_remove(pathString , "Endoscopy/")) %>%
  filter(!is.na(expected_cost_2024) & pathString != "Endoscopy") %>%
  mutate(group = sub("\\/.*", "",pathString)) %>%  
  group_by(group) %>% 
  summarise(total_cost_2024 = sum(expected_cost_2024), 
            total_cost_2025 = sum(expected_cost_2025), 
            total_cost_2026 = sum(expected_cost_2026), 
            total_cost_2027 = sum(expected_cost_2027), 
            total_cost_2028 = sum(expected_cost_2028))

# Store projected net costs
net_projcosts <- diff(as.matrix(tree_proj_df[c(2,1),-1]))

# Convert to dataframe and format for plotting
net_projcosts_df <- as.data.frame(net_projcosts) %>%
  pivot_longer(cols = c(1:5), names_to = "year", values_to = "cost") %>%
  mutate(year = as.numeric(str_sub(year,-4,-1))) %>%
  mutate(time = year - 2022,
         discount_rate = 0.03,
         pv = cost / (1+discount_rate)^time)

# Sum projected costs across five years
net_projcosts_df %>%
  summarise(base = sum(pv))

# Plot
ggplot(net_projcosts_df) +
  geom_line(aes(x = year, y = cost), col = "#E8464E", lty = 1, linewidth = 0.75) +
  geom_point(aes(x = year, y = cost), col = "#E8464E") +
  scale_x_continuous(expand = c(0.005,0.005)) + 
  scale_y_continuous(limits = c(0,4e5), breaks = seq(0,4e5,1e5),
                     expand = c(0,0), labels = scales::comma) +
  labs(x = "Year", y = "Projected excess cost, USD") +
  theme_mgh() +
  theme(axis.text.y = element_text(margin = margin(t = 0, r = 8, b = 0, l = 0)),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))


