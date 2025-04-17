##################################
##### SANDWICH MATRIX PROBIT #####
##################################
library(pracma)
library(nleqslv)

y = 1
a = 1 
m = 1
l = 1

qvec <- function(
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
    ((y - pnorm(b0 + ba*a + bm*m + bl*l)) / dnorm(b0 + ba*a + bm*m + bl*l)),
    ((y - pnorm(b0 + ba*a + bm*m + bl*l)) / dnorm(b0 + ba*a + bm*m + bl*l) * a),
    ((y - pnorm(b0 + ba*a + bm*m + bl*l)) / dnorm(b0 + ba*a + bm*m + bl*l) * m), 
    ((y - pnorm(b0 + ba*a + bm*m + bl*l)) / dnorm(b0 + ba*a + bm*m + bl*l) * l),
    
    # score mediator model (probit)
    ((m - pnorm(gam0 + gama*a + gaml*l)) / dnorm(gam0 + gama*a + gaml*l)),
    ((m - pnorm(gam0 + gama*a + gaml*l)) / dnorm(gam0 + gama*a + gaml*l) * a), 
    ((m - pnorm(gam0 + gama*a + gaml*l)) / dnorm(gam0 + gama*a + gaml*l) * l),
    
    # Mediator potential outcomes
    ((pnorm(gam0 + gama*1 + gaml*l) - E_M1_1) * a),
    ((pnorm(gam0 + gama*1 + gaml*l) - E_M1_0) * (1-a)),
    ((pnorm(gam0 + gama*0 + gaml*l) - E_M0_1) * a),
    ((pnorm(gam0 + gama*0 + gaml*l) - E_M0_0) * (1-a)),
    # mediator effect on the outcome 
    ((pnorm(b0 + ba*0 + bm*1 + bl*l) - pnorm(b0 + ba*0 + bm*0 + bl*l) - E_I01I00_1) * a),
    ((pnorm(b0 + ba*0 + bm*1 + bl*l) - pnorm(b0 + ba*0 + bm*0 + bl*l) - E_I01I00_0) * (1 - a)),
    # exposure effect on the outcome
    ((pnorm(b0 + ba*1 + bm*0 + bl*l) - pnorm(b0 + ba*0 + bm*0 + bl*l) - E_I10I00_1) * a),
    ((pnorm(b0 + ba*1 + bm*0 + bl*l) - pnorm(b0 + ba*0 + bm*0 + bl*l) - E_I10I00_0) * (1 - a)),
    ((pnorm(b0 + ba*1 + bm*1 + bl*l) - pnorm(b0 + ba*0 + bm*1 + bl*l) - E_I11I01_1) * a),
    ((pnorm(b0 + ba*1 + bm*1 + bl*l) - pnorm(b0 + ba*0 + bm*1 + bl*l) - E_I11I01_0) * (1 - a)),    
    # Indirect effects
    (((E_M1_1 - E_M0_1) * E_I01I00_1  - pi1) * a),
    (((E_M1_0 - E_M0_0) * E_I01I00_0  - pi0) * ( 1 - a)),
    ( pi0 * (1-a) + pi1 * a - pi ), 
    # Direct effects
    ((E_I10I00_1 * (1 - E_M1_1) + E_I11I01_1 *  E_M1_1 - pd1) * a),
    ((E_I10I00_0 * (1 - E_M1_0) + E_I11I01_0 *  E_M1_0 - pd0) * (1-a)),
    ( pd0 * (1-a) + pd1 * a - pd ), 
    # Indirect indices    
    (1/pi0         - INNE),
    (1/pi1         - IEIN),
    (1/pi          - INNT),
    # Direct indices    
    (1/pd0         - DNNE),
    (1/pd1         - DEIN),
    (1/pd          - DNNT),
    # marginal indices
    (1/(pi0 + pd0)                    - NNE),
    (1/(pi1 + pd1)                    - EIN),
    (1/(pi  + pd)                     - NNT)
  )
  return(out)
}


qvec2 <- function(x) {
  out <- qvec(#y = 1, a = 1, z = 1, 
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
#nleqslv(x = rep(1, 30), fn = qvec2)$x
# 
