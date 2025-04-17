####################
### TRUE INDICES ###
####################
####################
### PROBIT MODEL ###
####################

set.seed(123)
N = 10 ^ 7
    L <- rnorm(N, 0.5, 0.1)
    
    A <- rbinom(N, 1, plogis(2  - 3 * L ))
    
    # CONDITIONAL MEDIATOR-EXPOSURE MODEL
    M <- rbinom(N, 1, pnorm(-1 + 3 * A - 2 * L))
    
    # CONDITIONAL OUTCOME-MEDIATOR MODEL
    Y <- rbinom(N, 1, pnorm(-1 + 1.5 * A + 1.5 * M - 2 * L))
    
    # estimating conditional outcome-mediator model
    #drct_b <- coef(glm(Y ~ A + M + L, family = "binomial"))
    
    # estimating conditional mediator-exposure model
    #mdtr_b <- coef(glm(M ~ A + L, family = "binomial"))
    
    # marginal crude pb
    
    # conditional mediator-exposure model as a function of L
    pimL <- function(L) { 
      p_out <-  pnorm(-1 + 3 * 1 - 2 * L) - 
                pnorm(-1 + 3 * 0 - 2 * L)
      return(p_out)
    }
    
    # conditional outcome-mediator model for A=0 as a function of L
    pioML <- function(L) { 
      p_out <-  pnorm(-1 + 1.5 * 0 + 1.5 * 1 - 2 * L ) - 
        pnorm(-1 + 1.5 * 0 + 1.5 * 0 - 2 * L)
      return(p_out)
    }
    
    # EST IEIN
    p_i1 <- mean(pimL(L[A == 1])) * mean(pioML(L[A == 1]))
    
    # IEIN[k, n] <- 1 / p_i1
    IEIN <- 1 / p_i1
    
    # EST INNE
    p_i0 <- mean(pimL(L[A == 0])) * mean(pioML(L[A == 0]))
    # INNE[k, n] <- 1 / p_i0 
    INNE <- 1 / p_i0
    
    # EST INNT g(E[p(A)])
    p_i <- p_i0 * mean(A == 0) + p_i1 * mean(A == 1)
    # INNT[k, n] <- 1 / p_i
    INNT <- 1 / p_i
    
    # TRUE INNT direct calculation (mediation formula)
    # 1/( mean(pimL(L)) * mean(pioL(L)) )
    
    #############################
    # outcome model for M=0 as a function of L
    pioAM0L <- function(L) { 
      p_out <-  pnorm(-1 + 1.5 * 1 + 1.5 * 0 - 2 * L) -
        pnorm(-1 + 1.5 * 0 + 1.5 * 0 - 2 * L)
      return(p_out)
    }
    
    # outcome model for M=1 as a function of L
    pioAM1L <- function(L) { 
      p_out <-  pnorm(-1 + 1.5 * 1 + 1.5 * 1 - 2 * L) - 
        pnorm(-1 + 1.5 * 0 + 1.5 * 1 - 2 * L)
      return(p_out)
    }
    
    # TRUE DEIN
    p_d1 <- ( mean(pioAM0L(L[A == 1])) * ( 1 - mean(pnorm( -1 + 3 * 1 - 2 * L[A == 1])) ) +
                mean(pioAM1L(L[A == 1])) * ( mean(pnorm(-1 + 3 * 1 - 2 * L[A == 1] )) ) )       
    # DEIN[k, n] <- 1 / p_d1
    DEIN <- 1 / p_d1
    
    # TRUE DNNE
    p_d0 <-  mean(pioAM0L(L[A == 0])) * ( 1 - mean(pnorm( -1 + 3 * 1 - 2 * L[A == 0] )) ) +
                mean(pioAM1L(L[A == 0])) * mean(pnorm(-1 + 3 * 1 - 2 * L[A == 0]))          
    # DNNE[k, n] <- 1 / p_d0
    DNNE <- 1 / p_d0
    
    # TRUE DNNT marginalized over A
    p_d <- ( p_d0 * mean(A == 0) + p_d1 * mean(A == 1) )
    # DNNT[k, n] <- 1 / p_d
    DNNT <- 1 / p_d
    
    # TRUE INNT - direct calculations
    # 1/( mean(pioAM0L(L)) * ( 1 - mean(pnorm(-1 + 3 * 1 - 2 * L )) ) +
    #       mean(pioAM1L(L)) * ( mean(pnorm(-1 + 3 * 1 - 2 * L )) ) )
    
    ### MARGINAL EIN, NNE, and NNT ###
    # EIN
    p_b1 <- p_i1 + p_d1
    # EIN[k, n] <- 1 / p_b1
    EIN <- 1 / p_b1
    
    # NNE
    p_b0 <- p_i0 + p_d0
    # NNE[k, n] <- 1 / p_b0
    NNE <- 1 / p_b0
    
    # NNT
    p_b <- p_i + p_d
    # NNT[k, n] <- 1 / p_b
    NNT <- 1 / p_b
    
# Create the data frame of true index values
true_indices_df <- data.frame(
  INNE = INNE,
  IEIN = IEIN,
  INNT = INNT,
  DNNE = DNNE,
  DEIN = DEIN,
  DNNT = DNNT,
  NNE  = NNE,
  EIN  = EIN,
  NNT  = NNT
)

# View it
print(true_indices_df)

########################################
### solving system of est. equations ###
########################################
#################################
##### SANDWICH MATRIX LOGIT #####
#################################
library(pracma)
library(nleqslv)

set.seed(123)
N = 10 ^ 4

L <- rnorm(N, 0.5, 0.1)

A <- rbinom(N, 1, plogis(2  - 3 * L ))

# CONDITIONAL MEDIATOR-EXPOSURE MODEL
M <- rbinom(N, 1, pnorm(-1 + 3 * A - 2 * L))

# CONDITIONAL OUTCOME-MEDIATOR MODEL
Y <- rbinom(N, 1, pnorm(-1 + 1.5 * A + 1.5 * M - 2 * L))

y = Y
a = A 
m = M
l = L

qvec_sum <- function(
    b0, ba, bm, bl,
    gam0, gama, gaml,
    E_M1_1, E_M1_0, 
    E_M0_1, E_M0_0,
    E_I01I00_0, E_I01I00_1,
    E_I10I00_0, E_I10I00_1, 
    E_I11I01_0, E_I11I01_1,
    pi0, pi1, pi,
    pd0, pd1, pd,
    INNE, IEIN, INNT,
    DNNE, DEIN, DNNT,
    NNE, EIN, NNT
) {
  out <- c(
    # score outcome model (probit)
    sum((y - pnorm(b0 + ba*a + bm*m + bl*l)) / dnorm(b0 + ba*a + bm*m + bl*l)),
    sum((y - pnorm(b0 + ba*a + bm*m + bl*l)) / dnorm(b0 + ba*a + bm*m + bl*l) * a),
    sum((y - pnorm(b0 + ba*a + bm*m + bl*l)) / dnorm(b0 + ba*a + bm*m + bl*l) * m),
    sum((y - pnorm(b0 + ba*a + bm*m + bl*l)) / dnorm(b0 + ba*a + bm*m + bl*l) * l),
    
    # score mediator model (probit)
    sum((m - pnorm(gam0 + gama*a + gaml*l)) / dnorm(gam0 + gama*a + gaml*l)),
    sum((m - pnorm(gam0 + gama*a + gaml*l)) / dnorm(gam0 + gama*a + gaml*l) * a),
    sum((m - pnorm(gam0 + gama*a + gaml*l)) / dnorm(gam0 + gama*a + gaml*l) * l),
    
    # Mediator potential outcomes
    sum((pnorm(gam0 + gama*1 + gaml*l) - E_M1_1) * a),
    sum((pnorm(gam0 + gama*1 + gaml*l) - E_M1_0) * (1-a)),
    sum((pnorm(gam0 + gama*0 + gaml*l) - E_M0_1) * a),
    sum((pnorm(gam0 + gama*0 + gaml*l) - E_M0_0) * (1-a)),
    
    # mediator effect on the outcome 
    sum((pnorm(b0 + ba*0 + bm*1 + bl*l) - pnorm(b0 + ba*0 + bm*0 + bl*l) - E_I01I00_1) * a),
    sum((pnorm(b0 + ba*0 + bm*1 + bl*l) - pnorm(b0 + ba*0 + bm*0 + bl*l) - E_I01I00_0) * (1 - a)),
    
    # exposure effect on the outcome
    sum((pnorm(b0 + ba*1 + bm*0 + bl*l) - pnorm(b0 + ba*0 + bm*0 + bl*l) - E_I10I00_1) * a),
    sum((pnorm(b0 + ba*1 + bm*0 + bl*l) - pnorm(b0 + ba*0 + bm*0 + bl*l) - E_I10I00_0) * (1 - a)),
    sum((pnorm(b0 + ba*1 + bm*1 + bl*l) - pnorm(b0 + ba*0 + bm*1 + bl*l) - E_I11I01_1) * a),
    sum((pnorm(b0 + ba*1 + bm*1 + bl*l) - pnorm(b0 + ba*0 + bm*1 + bl*l) - E_I11I01_0) * (1 - a)),    
    # Indirect effects
    sum(((E_M1_1 - E_M0_1) * E_I01I00_1  - pi1) * a),
    sum(((E_M1_0 - E_M0_0) * E_I01I00_0  - pi0) * ( 1 - a)),
    sum( pi0 * (1-a) + pi1 * a - pi ), 
    # Direct effects
    sum((E_I10I00_1 * (1 - E_M1_1) + E_I11I01_1 *  E_M1_1 - pd1) * a),
    sum((E_I10I00_0 * (1 - E_M1_0) + E_I11I01_0 *  E_M1_0 - pd0) * (1-a)),
    sum( pd0 * (1-a) + pd1 * a - pd ), 
    # Indirect indices    
    sum(1/pi0         - INNE),
    sum(1/pi1         - IEIN),
    sum(1/pi          - INNT),
    # Direct indices    
    sum(1/pd0         - DNNE),
    sum(1/pd1         - DEIN),
    sum(1/pd          - DNNT),
    # marginal indices
    sum(1/(pi0 + pd0)                    - NNE),
    sum(1/(pi1 + pd1)                    - EIN),
    sum(1/(pi  + pd)                     - NNT)
  )
  return(out)
}


qvec2_sum <- function(x) {
  out <- qvec_sum(#y = 1, a = 1, z = 1, 
    b0 = x[1], ba = x[2], bm = x[3], bl = x[4],
    gam0 = x[5], gama = x[6], gaml = x[7],
    E_M1_1 = x[8], E_M1_0 = x[9], 
    E_M0_1 = x[10], E_M0_0 = x[11],
    E_I01I00_0 = x[12], E_I01I00_1 = x[13],
    E_I10I00_0 = x[14], E_I10I00_1 = x[15], 
    E_I11I01_0 = x[16], E_I11I01_1 = x[17],
    pi0 = x[18], pi1 = x[19], pi = x[20], 
    pd0 = x[21], pd1 = x[22], pd = x[23],
    INNE = x[24], IEIN = x[25], INNT = x[26],
    DNNE = x[27], DEIN = x[28], DNNT = x[29],
    NNE = x[30], EIN = x[31], NNT = x[32]
  )
  return(out)
}

#####################
### reality check ###
#####################

# ### derivative
# inv_a <- solve(-jacobian(f = qvec2, x0 = rep(1, 30)))
# 
# ### outer product
# outer(qvec2(x = rep(1, 30)), 
#       qvec2(x = rep(1, 30)) )
# 
### solving the system of linear equations
sol = nleqslv(x = c(rep(1, 32)), fn = qvec2_sum)$x
# 

true_ind <- sol[24:32]
true_ind
