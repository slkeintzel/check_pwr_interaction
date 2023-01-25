# Power calculation by tools
# R Script for Power Simulation for interaction effect in 2 x 2 between subjects ANOVA Design 
# Authors: Salomé Li Keintzel, Frauke Langhein and Paula G. Rosendahl (2023)


# G*Power 3.1.9.7 ####
# Selected options:
# ANOVA fixed effects, special, main effects and interaction
# Post hoc: Compute achieved power – given alpha, sample size and effect size
# Effect size f = 0.4, alpha = 0.05, total sample size n_vektor, numerator 
# df = 1, number of groups = 4

power_GPower <- c(0.3905143, # n = 20
                  0.6921406,
                  0.8610110,
                  0.9420455,
                  0.9772343,
                  0.9914667,
                  0.9969209,
                  0.9989237,
                  0.9996338,
                  0.9998782,
                  0.9999603,
                  0.9999873,
                  0.9999960 # n = 260
)

# quick plot
plot(x = n_vektor, y = power_GPower, 
     type = "l",
     xlab = "total sample size",
     ylab = "Power Calculated",
     main  = "G*Power 3.1.9.7",
     col = "darkorchid4"
)


# Webpower ####
# Selected options:
# Two-Way, three-Way, and k-Way ANOVA
# Number of groups = 4, numerator df = 1, Effect size (f) = 0.4, 
# significance level = 0.05, total sample size n_vektor

power_Webpower <- c(0.3905, # n = 20
                    0.6921,
                    0.8610,
                    0.9420,
                    0.9772,
                    0.9915,
                    0.9969,
                    0.9989,
                    0.9996,
                    0.9999,
                    1.0000,
                    1.0000,
                    1.0000 # n = 260
)

plot(x = n_vektor, y = power_Webpower, 
     type = "l",
     xlab = "total sample size",
     ylab = "Power Calculated",
     main  = "Webpower",
     col = "darkorchid4"
)

# Superpower ####
# library(Superpower) # ‘0.2.0’

#1) ordinal Interaction
# get means vector for this interaction type
mu_ordinal <- mu_list[[1]]

# get sde for this interacton type and fixed f
parameters_ordinal <- make_parameters_for_f(mu = mu_ordinal, n_rel_matr = balanced_n_matrix_rel, f = f)
sd_ordinal <- parameters_ordinal[["sde"]]

# powercurve vector (takes some time)
power_Superpower_ordinal <- sapply(n_vektor, function(current_total_n) { # loop for sample sizes
  
  design_ordinal <- ANOVA_design(design = "2b*2b",      # 2x2 between-subjects Design
                                 n = current_total_n/4, # needs sample size per cell as input
                                 mu = mu_ordinal,       # means per cell
                                 sd = sd_ordinal)       # residual sd = sde (that is, sd in each cell)
  
  power_result <- ANOVA_power(
    design_ordinal,
    alpha_level = 0.05,
    correction = "none", # unnecessary because between design
    nsims = runs, # same as in custom simulation 
    seed = NULL, # no seed as in custom simulation
    verbose = F,
    emm = F
  )
  
  return(power_result[["main_results"]][["power"]][3]) # return power
})

# save
save(list = c("power_Superpower_ordinal"), file = paste0("./", Sys.Date(), "power_Superpower_ordinal.Rdata"))

#2) hybrid interaction
# get means vector for this interaction type
mu_hybrid <- mu_list[[2]]

# get sde for this interacton type and fixed f
parameters_hybrid <- make_parameters_for_f(mu = mu_hybrid, n_rel_matr = balanced_n_matrix_rel, f = f)
sd_hybrid <- parameters_hybrid[["sde"]]

# powercurve vector
power_Superpower_hybrid <- sapply(n_vektor, function(current_total_n) { # current_total_n = total n
  
  design_hybrid <- ANOVA_design(design = "2b*2b",        # 2x2 between-subjects Design
                                n = current_total_n/4,   # needs sample size per cell as input
                                mu = mu_hybrid,          # means per cell
                                sd = sd_hybrid)          # residual sd = sde (that is, sd in each cell)
  
  power_result <- ANOVA_power(
    design_hybrid,
    alpha_level = 0.05,
    correction = "none", # unnecessary because between design
    nsims = runs, # same as in custom simulation 
    seed = NULL, # no seed as in custom simulation
    verbose = F,
    emm = F
  )
  
  return(power_result[["main_results"]][["power"]][3])
})

# save
save(list = c("power_Superpower_hybrid"), file = paste0("./", Sys.Date(), "power_Superpower_hybrid.Rdata"))

#2) Disordinal interaction
# get means vector for this interaction type
mu_disordinal <- mu_list[[3]]

# get sde for this interacton type and fixed f
parameters_disordinal <- make_parameters_for_f(mu = mu_disordinal, n_rel_matr = balanced_n_matrix_rel, f = f)
sd_disordinal <- parameters_disordinal[["sde"]]

# powercurve vector
power_Superpower_disordinal <- sapply(n_vektor, function(current_total_n) { # current_total_n = total n
  
  design_disordinal <- ANOVA_design(design = "2b*2b",      # 2x2 between-subjects Design
                                    n = current_total_n/4, # needs sample size per cell as input
                                    mu = mu_disordinal,    # means per cell
                                    sd = sd_disordinal)    # residual sd = sde (that is, sd in each cell)
  
  
  power_result <- ANOVA_power(
    design_disordinal,
    alpha_level = 0.05,
    correction = "none", # unnecessary because between design
    nsims = runs, # same as in custom simulation 
    seed = NULL, # no seed as in custom simulation
    verbose = F,
    emm = F
  )
  
  return(power_result[["main_results"]][["power"]][3])
})

# save
save(list = c("power_Superpower_disordinal"), file = paste0("./", Sys.Date(), "power_Superpower_disordinal.Rdata"))

# quick plots Superpower

plot(x = n_vektor, y = power_Superpower_ordinal/100, # superpower returns in %, get to same metric as other powercurves
     type = "l",
     xlab = "total sample size",
     ylab = "Power Calculated",
     main  = "Superpower",
     col = "blue4"
)
lines(x = n_vektor, y = power_Superpower_hybrid/100, 
      type = "l",
      xlab = "total sample size",
      ylab = "Power Calculated",
      col = "chartreuse3"
)
lines(x = n_vektor, y = power_Superpower_disordinal/100, 
      type = "l",
      xlab = "total sample size",
      ylab = "Power Calculated",
      col = "firebrick"
)
legend("bottomright",
       pch = 16,
       legend = c("ordinal", "hybrid", "disordinal"),
       col = c("blue4", "chartreuse3", "firebrick")
) 
