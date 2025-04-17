############################
### logit coverage rates ###
############################
library(dplyr)
library(tidyr)
library(readr)

# Read the CI data
ci_data <- read_csv("all_ci_long_logit.csv")

# True values
true_values <- c(
  INNE = 6.525974,
  IEIN = 6.282517,
  INNT = 6.372880,
  DNNE = 3.081796,
  DEIN = 3.066874,
  DNNT = 3.072529,
  NNE  = 2.093277,
  EIN  = 2.060850,
  NNT  = 2.073056
)

# Step 1: Reshape to long and calculate coverage
ci_long <- ci_data %>%
  pivot_longer(cols = -c(n, iteration),
               names_to = "IndexBound",
               values_to = "Value") %>%
  mutate(
    Index = sub("_[LU]$", "", IndexBound),
    Bound = ifelse(grepl("_L$", IndexBound), "Lower", "Upper")
  ) %>%
  select(-IndexBound) %>%
  pivot_wider(names_from = Bound, values_from = Value) %>%
  filter(!(is.infinite(Lower) & is.infinite(Upper))) %>%  # <-- Exclude full-Inf intervals
  mutate(True = true_values[Index],
         Covered = (True >= Lower & True <= Upper))

# Step 2: Compute coverage by index and sample size
coverage_by_n <- ci_long %>%
  group_by(Index, n) %>%
  summarise(Coverage = mean(Covered, na.rm = TRUE), .groups = "drop")

# Step 3: Reshape to wide format
coverage_wide <- coverage_by_n %>%
  pivot_wider(names_from = n, values_from = Coverage)

# Step 4: Reorder by true_values vector
coverage_wide <- coverage_wide %>%
  mutate(Index = factor(Index, levels = names(true_values))) %>%
  arrange(Index)

# Step 5: Print coverage rates table
print(as.data.frame(coverage_wide), row.names = FALSE)

#######################
#### reality check ####
#######################
##### CIs #####
# Step 6: Count number of full-Inf CI rows per sample size and index
inf_count_by_index <- ci_data %>%
  pivot_longer(cols = -c(n, iteration),
               names_to = "IndexBound",
               values_to = "Value") %>%
  mutate(
    Index = sub("_[LU]$", "", IndexBound),
    Bound = ifelse(grepl("_L$", IndexBound), "Lower", "Upper")
  ) %>%
  select(-IndexBound) %>%
  pivot_wider(names_from = Bound, values_from = Value) %>%
  filter(is.infinite(Lower) & is.infinite(Upper)) %>%
  group_by(n, Index) %>%
  summarise(Num_Full_Inf = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = n, values_from = Num_Full_Inf, values_fill = 0)

cat("\nNumber of observations with both CI limits equal to Inf, by index and sample size:\n")
print(inf_count_by_index, row.names = FALSE)

##### POINT ESTIMATORS #####
# Read the data
est_data <- read_csv("all_est_long_logit.csv")

# Define desired order
desired_order <- c("DEIN", "DNNE", "DNNT", "EIN", "IEIN", "INNT", "NNE", "NNT")

# Count negative or infinite estimates
bad_estimates <- est_data %>%
  pivot_longer(cols = -c(n, iteration),
               names_to = "Index",
               values_to = "Estimate") %>%
  filter(is.infinite(Estimate) | Estimate < 1) %>%
  group_by(Index, n) %>%
  summarise(Num_Bad = n(), .groups = "drop") %>%
  filter(Index %in% desired_order)

# Reshape to wide format and reorder
bad_estimates_wide <- bad_estimates %>%
  pivot_wider(names_from = n, values_from = Num_Bad, values_fill = 0) %>%
  slice(match(desired_order, Index))

# Display
cat("Number of negative or infinite point estimates by sample size and index:\n")
print(as.data.frame(bad_estimates_wide), row.names = FALSE)
colSums(bad_estimates_wide[,-1])/(100 * 9)