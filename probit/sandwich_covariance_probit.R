######################################
### pnorm function - NNE, EIN, NNT ###
### ESTIMATION #######################
######################################
library(pracma)
library(nleqslv)

N <- c(200, 400, 800, 1600)
#N <- c(800)
K <- 100

all_ci  <- list()
all_est <- list()

set.seed(123)

for (n in 1:length(N)) {
  mat_est <- data.frame(matrix(NA, nrow = K, ncol = 9))
  
  mat_ci <- data.frame(matrix(NA, nrow = K, ncol = 9 * 2))
  
  colnames(mat_ci) <- c("INNE_L", "INNE_U",
                        "IEIN_L", "IEIN_U",
                        "INNT_L", "INNT_U",
                        "DNNE_L", "DNNE_U",
                        "DEIN_L", "DEIN_U",
                        "DNNT_L", "DNNT_U",
                        "NNE_L",  "NNE_U",
                        "EIN_L",  "EIN_U",
                        "NNT_L",  "NNT_U")
  
  colnames(mat_est) <- c("INNE", 
                         "IEIN",
                         "INNT",
                         "DNNE",
                         "DEIN",
                         "DNNT",
                         "NNE", 
                         "EIN", 
                         "NNT")
  
  for (k in 1:K) {
    print(c(n, k))
    
    L <- rnorm(N[n], 0.5, 0.1)
    
    A <- rbinom(N[n], 1, plogis(2 - 3 * L ))
    
    # CONDITIONAL MEDIATOR-EXPOSURE MODEL
    M <- rbinom(N[n], 1, pnorm(-1 + 3 * A - 2 * L))
    
    # CONDITIONAL OUTCOME-MEDIATOR MODEL
    Y <- rbinom(N[n], 1, pnorm(-1 + 1.5 * A + 1.5 * M - 2 * L))
    
    ########################################################   
    # estimating conditional outcome-mediator model
  
    y = Y
    a = A 
    m = M
    l = L
    
    sol = nleqslv(x = rep(1, 32), fn = qvec2_sum)$x
    # 
    #print(sol)
    # true_ind <- sol[22:30]
    
    mat_est[k,] <- sol[24:32]
    
    ########################  
    ##### CIs SANDWICH #####
    ########################
    
    ##### THE A BREAD MATRIX #####
    bread_mat <- matrix(0, 32, 32)
    
    for (b in 1:length(Y)) {
      y <-  Y[b]
      a <-  A[b]
      m <-  M[b]
      l <-  L[b]
      
      jac_mat <- -jacobian(f = qvec2, x0 = sol)
      bread_mat <- bread_mat + jac_mat
    }
    
    bread_mat <- 1/length(Y) * bread_mat
    #bread_mat
    
    if (!is.matrix(bread_mat) || any(!is.finite(bread_mat)) || nrow(bread_mat) != ncol(bread_mat)) {
      warning("Bread matrix invalid: not square or contains NA/NaN/Inf")
      next
    }
    
    if (rcond(bread_mat) < .Machine$double.eps) {
      warning(paste("Iteration skipped due to singular bread matrix"))
      next
    }
    
    inv_a <- solve(bread_mat)
    
    ##### THE B MEAT MATRIX #####
    meat_mat <- matrix(0, 32, 32)
    
    for (b in 1:length(Y)) {
      y <-  Y[b]
      a <-  A[b]
      m <-  M[b]
      l <-  L[b]
      
      
      out_prod <-  outer(qvec2(x = sol), 
                         qvec2(x = sol) )
      
      meat_mat <- meat_mat + out_prod
    }
    
    
    meat_mat <- 1/length(Y) * meat_mat
    # meat_mat
    
    ##### THE SANDWICH inv(A) %*% B %*% t(inv(A)) MATRIX ##### 
    
    sand_mat <- 1 / length(Y) * inv_a %*% meat_mat %*% t(inv_a) 
    
    mat_ci[k,] = c(
      if (sol[24] >= 1) c(max(sol[24] - 1.96 * sqrt(sand_mat[24, 24]), 1), sol[24] + 1.96 * sqrt(sand_mat[24, 24])) else c(Inf, Inf),
      if (sol[25] >= 1) c(max(sol[25] - 1.96 * sqrt(sand_mat[25, 25]), 1), sol[25] + 1.96 * sqrt(sand_mat[25, 25])) else c(Inf, Inf),
      if (sol[26] >= 1) c(max(sol[26] - 1.96 * sqrt(sand_mat[26, 26]), 1), sol[26] + 1.96 * sqrt(sand_mat[26, 26])) else c(Inf, Inf),
      if (sol[27] >= 1) c(max(sol[27] - 1.96 * sqrt(sand_mat[27, 27]), 1), sol[27] + 1.96 * sqrt(sand_mat[27, 27])) else c(Inf, Inf),
      if (sol[28] >= 1) c(max(sol[28] - 1.96 * sqrt(sand_mat[28, 28]), 1), sol[28] + 1.96 * sqrt(sand_mat[28, 28])) else c(Inf, Inf),
      if (sol[29] >= 1) c(max(sol[29] - 1.96 * sqrt(sand_mat[29, 29]), 1), sol[29] + 1.96 * sqrt(sand_mat[29, 29])) else c(Inf, Inf),
      if (sol[30] >= 1) c(max(sol[30] - 1.96 * sqrt(sand_mat[30, 30]), 1), sol[30] + 1.96 * sqrt(sand_mat[30, 30])) else c(Inf, Inf),
      if (sol[31] >= 1) c(max(sol[31] - 1.96 * sqrt(sand_mat[31, 31]), 1), sol[31] + 1.96 * sqrt(sand_mat[31, 31])) else c(Inf, Inf),
      if (sol[32] >= 1) c(max(sol[32] - 1.96 * sqrt(sand_mat[32, 32]), 1), sol[32] + 1.96 * sqrt(sand_mat[32, 32])) else c(Inf, Inf)
    )
    
    #print(i)
  }
  all_ci[[as.character(N[n])]]  <- mat_ci
  all_est[[as.character(N[n])]] <- mat_est
}

# Named vector of true values (not estimating equations - MC, sigma = 0.1, N = 10^7)
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


# # Compute coverage using a for loop
# index <- names(true_values)
# true_value <- true_values
# coverage <- numeric(length(true_values))
# 
# for (i in seq_along(true_values)) {
#   name <- index[i]
#   coverage[i] <- mean(true_value[i] >= mat_ci[, paste0(name, "_L")] &
#                         true_value[i] <= mat_ci[, paste0(name, "_U")], na.rm = T)
# }
# 
# # Combine into a data frame
# coverage_df <- data.frame(index, true_value, coverage)
# 
# # View result
# print(coverage_df)
# 
# # number of NAs
# na_rows <- sum(apply(mat_ci, 1, anyNA))
# cat("Number of rows with NA values:", na_rows, "\n")


coverage_df <- data.frame(Index = names(true_values))

for (n in names(all_ci)) {
  mat <- all_ci[[n]]
  coverage <- sapply(names(true_values), function(name) {
    mean(true_values[name] >= mat[, paste0(name, "_L")] & true_values[name] <= mat[, paste0(name, "_U")], na.rm = T)
  })
  coverage_df[[paste0("N", n)]] <- coverage
}

# View
print(coverage_df)

######################
### csv of all CIs ###
######################

ci_long <- do.call(rbind, lapply(names(all_ci), function(n) {
  df <- all_ci[[n]]
  df$n <- as.integer(n)
  df$iteration <- 1:nrow(df)
  df
}))

# Move `n` to first column
ci_long <- ci_long[, c("n", "iteration", setdiff(names(ci_long), c("n", "iteration")))]

# Save
write.csv(ci_long, "all_ci_long_probit.csv", row.names = FALSE)

#######################
### csv of all ESTs ###
#######################

est_long <- do.call(rbind, lapply(names(all_est), function(n) {
  df <- all_est[[n]]
  df$n <- as.integer(n)
  df$iteration <- 1:nrow(df)
  df
}))

# Reorder columns
est_long <- est_long[, c("n", "iteration", setdiff(names(est_long), c("n", "iteration")))]

# Save
write.csv(est_long, "all_est_long_probit.csv", row.names = FALSE)


# library(xtable)
# xtable(all_ci)
