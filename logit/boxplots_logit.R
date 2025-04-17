#####################
### boxplot logit ###
#####################
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(grid)

# 1. Read the CSV
all_est_long <- read.csv("all_est_long_logit.csv")

# 2. Reshape to long format
est_long_long <- all_est_long %>%
  pivot_longer(cols = c(INNE, IEIN, INNT, DNNE, DEIN, DNNT, NNE, EIN, NNT),
               names_to = "Index",
               values_to = "Estimate")

# 3. Define true values (before applying truncation logic)
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
true_values <- round(true_values, 2)

# 4. Truncate based on condition: if < 1 → Inf; if > 3 × true_value → NA
est_long_long <- est_long_long %>%
  rowwise() %>%
  mutate(Estimate = ifelse(Estimate < 1, Inf,
                           ifelse(Estimate > 3 * true_values[Index], NA, Estimate))) %>%
  ungroup()

# 5. Ensure sample size is treated as factor
est_long_long$n <- factor(est_long_long$n, levels = sort(unique(est_long_long$n)))

# 6. Create plots in 3x3 layout
indices_order <- list(
  c("INNE", "IEIN", "INNT"),
  c("DNNE", "DEIN", "DNNT"),
  c("NNE",  "EIN",  "NNT")
)

plots <- lapply(unlist(indices_order), function(index) {
  ggplot(est_long_long %>% filter(Index == index),
         aes(x = n, y = Estimate)) +
    geom_boxplot(fill = "grey80") +
    geom_hline(yintercept = true_values[index], 
               color = "red", linetype = "dashed", linewidth = 1) +
    labs(title = index, x = "Sample size", y = "Estimator") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(size = 12))
})

# 7. Two-line title
title_grob <- textGrob(
  "Point Estimators of the Indirect, Direct, and Marginal Indices \n in Double Logit Model as a Function of the Sample Size",
  gp = gpar(fontsize = 12), just = "center"
)

# 8. Save as png
png("logit_estimator_boxplots.png", width = 1600, height = 1200, res = 200)
grid.arrange(grobs = plots, nrow = 3, ncol = 3, top = title_grob)
dev.off()

# 9. Also show on screen
grid.arrange(grobs = plots, nrow = 3, ncol = 3, top = title_grob)
