#####################
### boxplot logit ###
#####################
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(grid)

# 1. Read the CSV
all_est_long <- read.csv("all_est_long_probit.csv")

# 2. Reshape to long format
est_long_long <- all_est_long %>%
  pivot_longer(cols = c(INNE, IEIN, INNT, DNNE, DEIN, DNNT, NNE, EIN, NNT),
               names_to = "Index",
               values_to = "Estimate")

# 3. Truncate values: if < 1, set to 1; if > 30, set to Inf
est_long_long <- est_long_long %>%
  mutate(Estimate = ifelse(Estimate < 1, 1,
                           ifelse(Estimate > 10, Inf, Estimate)))

# 4. Define true values
true_values <- c(
  INNE = 4.493874,
  IEIN = 4.176306,
  INNT = 4.291575,
  DNNE = 2.062010,
  DEIN = 2.056743,
  DNNT = 2.058742,
  NNE  = 1.413450,
  EIN  = 1.378072,
  NNT  = 1.391308
)
true_values <- round(true_values, 2)

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
    theme_minimal(base_size = 12)
})

# 7. Two-line title
title_grob <- textGrob(
  "Point Estimators of the Indirect, Direct, and Marginal Indices \n in Double Probit Model as a Function of the Sample Size",
  gp = gpar(fontsize = 12), just = "center"
)

# 8. Save as JPEG
png("probit_estimator_boxplots.png", width = 1600, height = 1200, res = 200)
grid.arrange(grobs = plots, nrow = 3, ncol = 3, top = title_grob)
dev.off()

# 10. Also show on screen
grid.arrange(grobs = plots, nrow = 3, ncol = 3, top = title_grob)
