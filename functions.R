# Custom functions for
# R Script for Power Simulation for interaction effect in 2 x 2 between subjects ANOVA Design 
# Authors: Salomé Li Keintzel, Frauke Langhein and Paula G. Rosendahl (2023)


#### functions ----

# calculate residual sd (sde) and betas (effect coding) for specific mean 
# structure and effect size f
make_parameters_for_f <- function(mu, # means vector (cells 00,01,10,11) 
                                  n_rel_matr, # 3x3 matrix of relative n per cell + marginal
                                  f = 0.4 # Cohens f
                                  ) { 
  
  # means vector to 3x3 means matrix with marginal means 
  mu_matr_small <- matrix(mu, nrow = 2, byrow = T, dimnames = list(c("C_b = 0", "C_b = 1"), 
                                                                   c("C_a = 0", "C_a = 1")))
  C_b_means <- rowMeans(mu_matr_small)
  C_a_means <- colMeans(mu_matr_small)
  mu_matr <- rbind(cbind(mu_matr_small, C_b_means), C_a_means = c(C_a_means, mean(mu_matr_small)))
  
  # calculate betas (based on (unweiged) effect coding)
  b0 = mean(mu_matr[1:2, 1:2]) # unweighed mean 
  b1 = mu_matr["C_a_means", "C_a = 1"] - b0 # difference C_a = 1 from unweighed mean
  b2 = mu_matr["C_b = 1", "C_b_means"] - b0 # difference C_b = 1 from unweighed mean
  b3 = mu_matr["C_b = 1", "C_a = 1"] - b0 - b1 - b2 # for interaction term
  
  # Formula for retrieving sd_res for given f (see Cohens_d_sde.md for more info)
  sde = b3/f
  
  return(list(mu = mu,
              mu_matr = mu_matr,
              n_rel_matr = n_rel_matr,
              f = f,
              betas = c(b0, b1, b2, b3),
              sde = sde))
}

# Create data frame for simulation
sampling_dicho <- function(n_rel_matr, # relative n per cell & marginal
                           n_sample, # total n 
                           meane = 0, # mean of residuals set to 0
                           sde, # sd of residuals 
                           betas # vector of beta weights for effect coding
) {  
  
  n_tot_matr <- n_rel_matr * n_sample # total n per cell & marginal
  
  # Matrix with efefct coding for cell A=0, B=0
  # 2 columns (factor A and B) and ngroup["C_b = 1","C_a = 1"] rows (number of 
  # participant in cell)
  A0B0 <- cbind(rep(-1, times = n_tot_matr["C_b = 0","C_a = 0"]),       
                rep(-1, times = n_tot_matr["C_b = 0","C_a = 0"])) 
  
  # Matrix with efefct coding for cell A=1, B=0
  A1B0 <- cbind(rep(1, times = n_tot_matr["C_b = 0","C_a = 1"]),    
                rep(-1, times = n_tot_matr["C_b = 0","C_a = 1"]))
  
  # Matrix with efefct coding for cell A=0, B=1
  A0B1 <- cbind(rep(-1, times = n_tot_matr["C_b = 1", "C_a = 0"]),
                rep(1, times = n_tot_matr["C_b = 1", "C_a = 0"]))            
  
  # Matrix with efefct coding for cell A=1, B=1
  A1B1 <- cbind(rep(1, times = n_tot_matr["C_b = 1","C_a = 1"]),  
                rep(1, times = n_tot_matr["C_b = 1","C_a = 1"]))     
  
  # rbind those matrices together
  xmat <- rbind(A0B0, A1B0, A0B1, A1B1)
  colnames(xmat) <- c("A", "B")
  
  # add column for interaction term
  AxB <- xmat[,"A"] * xmat[,"B"]
  xmat <- cbind(xmat, AxB)  
  
  # define all parameters and objects
  e <- rnorm(n = n_sample, mean = meane, sd = sde)   # vector with residuals
  u <- rep(1, n_sample)                              # 1 vector for intercept
  X <- cbind(u, xmat)                                # design matrix for betas
  y <- X %*% betas + e                               # calculate dv 
  id <- 1:n_sample                                   # add id column (for afex)
  
  # bind final df together
  oursample <- data.frame(id,xmat, y)      
  return(oursample)
}

# calculates anova via afex
calculate_anova <- function(Data, Modell) {
  suppressMessages( # suppresses contrast and factor setting message in console
    anova_result <- aov_car(Modell, Data, 
                            type=3, # type 3 anova
                            anova_table = list(es = "pes"))) # return partial eta squared
  
  coeff<-anova_result[["lm"]][["coefficients"]] # save empirical betas
  pvalue<-anova_result[["Anova"]][["Pr(>F)"]] # save p values 
  pes <-anova_result[["anova_table"]][["pes"]] # save emiricap partial eta squared
  all<-cbind(coeff = coeff, pvalue=pvalue[1:4], pes = c(NA, pes)) # bind to one object
  return(all)
}
