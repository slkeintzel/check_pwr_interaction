# R Script for Power Simulation for interaction effect in 2 x 2 between subjects ANOVA Design 
# Authors: Salomé Li Keintzel, Frauke Langhein and Paula G. Rosendahl (2023)

# R and packages used:
# R Version 4.2.2
library(afex) # 1.2.0
library(parallel) # 4.2.1
library(Superpower) # ‘0.2.0’

# load custom functions
source("functions.R")

#### Simulation design parameters ----
# vector of (total) sample sizes to test power for
n_vektor <- seq(from = 20, to = 260, by = 20)

# vector of means c(00, 01, 10, 11) for our 3 interaction designs
mu_list <- list(
  mu_ordinal = c(0, 0, 0, 20), # ordinal
  mu_hybrid = c(0, -20, 0, 20), # hybrid
  mu_disordinal = c(20, -20, -20, 20) # disordinal 
)

# load different balancing matrices
source("matrices.R")

# model syntax for afex
model_gleichung <- as.formula("y ~ A * B + Error(id|1)")

# cohens for interaction effect to power for
f = 0.4

# mean of residuals set to 0 
meane <- 0 

#### sanity check: ---- 
# should return correct betas and effect size for balanced design
parameters <- make_parameters_for_f(mu = mu_list[[1]], n_rel_matr = balanced_n_matrix_rel, f = f)

# check input:
parameters["mu"] 
parameters["betas"] 
parameters["sde"]

# create sample
testsample <- sampling_dicho(n_rel_matr = balanced_n_matrix_rel, 
                              n_sample = 20000, 
                              meane = 0,
                              sde = parameters[["sde"]], 
                              betas = parameters[["betas"]] 
) 

# check output mean structure: 
aggregate(y ~ A + B, testsample, mean) 

# anova (coeff are calculated using different reversed contrast coding than for betas, this reverses some signs)
calculate_anova(Data = testsample, Modell = model_gleichung)


#### Simulation Loop ----

# number of simulated data frames per simulation condition
runs <- 20000

# for running parallelized:
cl <- makeForkCluster(nnodes = getOption("mc.cores", 10L))
clusterExport(cl = cl, varlist = ls())
# clusterEvalQ(cl,library(Iso)) # check if correct packages are loaded in cl 

print(paste("started", Sys.time()))

result_array = sapply(rel_matrices, FUN=function(n_rel_matr){ # loop for balancing matrices
  
  simulationen = sapply(mu_list, FUN=function(mu_vektor) { # loop for interaction designs
    
    # calculate betas and sde for fixed f and varying interaction design and balancing matrix
    parameters <- make_parameters_for_f(mu = mu_vektor, n_rel_matr = n_rel_matr, f = f)
    betas <- parameters[["betas"]]
    sde <- parameters[["sde"]]
  
    sapply(n_vektor, FUN=function(n_sample){ # loop for sample sizes
    
      # parallelized: 
      sample_matrix = parSapply(cl = cl, 1:runs, FUN=function(isample){ # loop for runs
      # for running without parallelization:
      # sample_matrix = sapply(1:runs, FUN=function(isample){ # loop for runs
        
        # create sample df  
        current_sample <- sampling_dicho(n_rel_matr = n_rel_matr, n_sample = n_sample,
                                       meane = 0, 
                                       sde = sde, 
                                       betas = betas)
        # calculate ANOVA
        results<-calculate_anova(Modell=model_gleichung, Data=current_sample)
        results
        
      }, simplify= "array")
    }, simplify="array")
  },simplify="array")
}, simplify="array")

stopCluster(cl = cl)

print(paste("done", Sys.time()))

# check dimensions of array
dimnames(result_array)
dim(result_array)

# name columns of result_array
dimnames(result_array)[3:4] <- list(paste0("Run ", 1:runs),paste0("Samplesize = ", n_vektor))

#### save result_array ----
save(list = c("result_array"), file = paste0("./", Sys.Date(), "simulationen.Rdata"), compress = "xz", compression_level = 9)

# recode pvalue to code significant (1) or not (0)
result_array[,"pvalue",,,,]<-result_array[,"pvalue",,,,] < 0.05

# mean array per interaction design + balanced power curve for interaction effect
res_arr_ordinal <- apply(result_array, MARGIN = c(1,2,4,5,6), FUN = mean)[,,,"mu_ordinal",]
power_ordinal_balanced <- c(res_arr_ordinal["A1:B1","pvalue",,"balanced_n_matrix_rel"])

res_arr_hybrid <- apply(result_array, MARGIN = c(1,2,4,5,6), FUN = mean)[,,,"mu_hybrid",]
power_hybrid_balanced <- c(res_arr_hybrid["A1:B1","pvalue",,"balanced_n_matrix_rel"])

res_arr_disordinal <- apply(result_array, MARGIN = c(1,2,4,5,6), FUN = mean)[,,,"mu_disordinal",]
power_disordinal_balanced <- c(res_arr_disordinal["A1:B1","pvalue",,"balanced_n_matrix_rel"])

# load results from packages (might take a while due to superpower simulating power)
source("tools.R")

#### plots for a quick overview ----
## Powercurves ####

### Ordinal ----
# balanced
title = "Ordinal design, balanced"
jpeg(paste0("./plots/", title, ".jpg"), width = 700, height = 700, quality = 100)

plot(x = n_vektor, y = power_ordinal_balanced, 
     type = "l",
     xlab = "total sample size",
     ylab = "Simulated Power",
     main = title,
     ylim = c(0.3, 1),
     col = "blue")
# add g*power
lines(x = n_vektor, y = power_GPower, 
      type = "l",
      col = "chartreuse3")
# add superpower 
lines(x = n_vektor, y = power_Superpower_ordinal*0.01, 
      type = "l",
      col = "firebrick")
# add webpower 
lines(x = n_vektor, y = power_Webpower, 
      type = "l",
      col = "orange")
# legend
legend("bottomright",
       pch = 16,
       legend = c("Simulation", "G*Power", "Superpower", "Webpower"),
       col = c("blue", "chartreuse3", "firebrick", "orange")
) 
dev.off()

# factors 0.6/0.4
title = "Ordinal design, factors 0.4 vs 0.6"
jpeg(paste0("./plots/", title, ".jpg"), width = 700, height = 700, quality = 100)

plot(x = n_vektor, y = res_arr_ordinal["A1:B1","pvalue",,"a0_matrix_rel"], 
     type = "l",
     xlab = "total sample size",
     ylab = "Simulated Power",
     main = title,
     ylim = c(0.3, 1),
     col = "blue")
lines(x = n_vektor, y = res_arr_ordinal["A1:B1","pvalue",,"a1_matrix_rel"], 
      type = "l",
      col = "lightblue")
lines(x = n_vektor, y = res_arr_ordinal["A1:B1","pvalue",,"b0_matrix_rel"], 
      type = "l",
      col = "turquoise")
lines(x = n_vektor, y = res_arr_ordinal["A1:B1","pvalue",,"b1_matrix_rel"], 
      type = "l",
      col = "darkblue")
# add g*power
lines(x = n_vektor, y = power_GPower, 
      type = "l",
      col = "chartreuse3")
# add superpower 
lines(x = n_vektor, y = power_Superpower_ordinal*0.01, 
      type = "l",
      col = "firebrick")
# add webpower 
lines(x = n_vektor, y = power_Webpower, 
      type = "l",
      col = "orange")
# legend
legend("bottomright",
       pch = 16,
       legend = c("simulation a0 0.6", "simulation a1 0.6", 
                  "simulation b0 0.6", "simulation b1 0.6", 
                  "G*Power", "Superpower", "Webpower"),
       col = c("blue", "lightblue", "turquoise", "darkblue", 
               "chartreuse3", "firebrick", "orange")
) 
dev.off()

# one cell 0.3/0.1 minor
title = "Ordinal design, one cell 0.3 vs 0.1 minor" 
jpeg(paste0("./plots/", title, ".jpg"), width = 700, height = 700, quality = 100)

plot(x = n_vektor, y = res_arr_ordinal["A1:B1","pvalue",,"a0b0_minor_matrix_rel"], 
     type = "l",
     xlab = "total sample size",
     ylab = "Simulated Power",
     main = title,
     ylim = c(0.3, 1),
     col = "blue")
lines(x = n_vektor, y = res_arr_ordinal["A1:B1","pvalue",,"a1b0_minor_matrix_rel"], 
      type = "l",
      col = "lightblue")
lines(x = n_vektor, y = res_arr_ordinal["A1:B1","pvalue",,"a0b1_minor_matrix_rel"], 
      type = "l",
      col = "turquoise")
lines(x = n_vektor, y = res_arr_ordinal["A1:B1","pvalue",,"a1b1_minor_matrix_rel"], 
      type = "l",
      col = "darkblue")
# add g*power
lines(x = n_vektor, y = power_GPower, 
      type = "l",
      col = "chartreuse3")
# add superpower 
lines(x = n_vektor, y = power_Superpower_ordinal*0.01, 
      type = "l",
      col = "firebrick")
# add webpower 
lines(x = n_vektor, y = power_Webpower, 
      type = "l",
      col = "orange")
# legend
legend("bottomright",
       pch = 16,
       legend = c("simulation a0b0 0.1", "simulation a1b0 0.1", 
                  "simulation a0b1 0.1", "simulation a1b1 0.1", 
                  "G*Power", "Superpower", "Webpower"),
       col = c("blue", "lightblue", "turquoise", "darkblue", 
                      "chartreuse3", "firebrick", "orange")
)  
dev.off()

# one cell 0.4/0.2 major
title = "Ordinal design, one cell 0.4 vs 0.2 major"
jpeg(paste0("./plots/", title, ".jpg"), width = 700, height = 700, quality = 100)

plot(x = n_vektor, y = res_arr_ordinal["A1:B1","pvalue",,"a0b0_major_matrix_rel"], 
     type = "l",
     xlab = "total sample size",
     ylab = "Simulated Power",
     main = title,
     ylim = c(0.3, 1),
     col = "blue")
lines(x = n_vektor, y = res_arr_ordinal["A1:B1","pvalue",,"a1b0_major_matrix_rel"], 
      type = "l",
      col = "lightblue")
lines(x = n_vektor, y = res_arr_ordinal["A1:B1","pvalue",,"a0b1_major_matrix_rel"], 
      type = "l",
      col = "turquoise")
lines(x = n_vektor, y = res_arr_ordinal["A1:B1","pvalue",,"a1b1_major_matrix_rel"], 
      type = "l",
      col = "darkblue")
# add g*power
lines(x = n_vektor, y = power_GPower, 
      type = "l",
      col = "chartreuse3")
# add superpower 
lines(x = n_vektor, y = power_Superpower_ordinal*0.01, 
      type = "l",
      col = "firebrick")
# add webpower 
lines(x = n_vektor, y = power_Webpower, 
      type = "l",
      col = "orange")
# legend
legend("bottomright",
       pch = 16,
       legend = c("simulation a0b0 0.4", "simulation a1b0 0.4", 
                  "simulation a0b1 0.4", "simulation a1b1 0.4", 
                  "G*Power", "Superpower", "Webpower"),
       col = c("blue", "lightblue", "turquoise", "darkblue", 
                      "chartreuse3", "firebrick", "orange")
)  
dev.off()

### Hybrid ----
# balanced
title = "Hybrid design, balanced"
jpeg(paste0("./plots/", title, ".jpg"), width = 700, height = 700, quality = 100)

plot(x = n_vektor, y = power_hybrid_balanced, 
     type = "l",
     xlab = "total sample size",
     ylab = "Simulated Power",
     main = title,
     ylim = c(0.3, 1),
     col = "blue")
# add g*power
lines(x = n_vektor, y = power_GPower, 
      type = "l",
      col = "chartreuse3")
# add superpower 
lines(x = n_vektor, y = power_Superpower_hybrid*0.01, 
      type = "l",
      col = "firebrick")
# add webpower 
lines(x = n_vektor, y = power_Webpower, 
      type = "l",
      col = "orange")
# legend
legend("bottomright",
       pch = 16,
       legend = c("simulation", "G*Power", "Superpower", "Webpower"),
       col = c("blue", "chartreuse3", "firebrick", "orange")
) 
dev.off()

# factors 0.6/0.4
title = "Hybrid design, factors 0.4 vs 0.6"
jpeg(paste0("./plots/", title, ".jpg"), width = 700, height = 700, quality = 100)

plot(x = n_vektor, y = res_arr_hybrid["A1:B1","pvalue",,"a0_matrix_rel"], 
     type = "l",
     xlab = "total sample size",
     ylab = "Simulated Power",
     main = title,
     ylim = c(0.3, 1),
     col = "blue")
lines(x = n_vektor, y = res_arr_hybrid["A1:B1","pvalue",,"a1_matrix_rel"], 
      type = "l",
      col = "lightblue")
lines(x = n_vektor, y = res_arr_hybrid["A1:B1","pvalue",,"b0_matrix_rel"], 
      type = "l",
      col = "turquoise")
lines(x = n_vektor, y = res_arr_hybrid["A1:B1","pvalue",,"b1_matrix_rel"], 
      type = "l",
      col = "darkblue")
# add g*power
lines(x = n_vektor, y = power_GPower, 
      type = "l",
      col = "chartreuse3")
# add superpower 
lines(x = n_vektor, y = power_Superpower_hybrid*0.01, 
      type = "l",
      col = "firebrick")
# add webpower 
lines(x = n_vektor, y = power_Webpower, 
      type = "l",
      col = "orange")
# legend
legend("bottomright",
       pch = 16,
       legend = c("simulation a0 0.6", "simulation a1 0.6", 
                  "simulation b0 0.6", "simulation b1 0.6", 
                  "G*Power", "Superpower", "Webpower"),
       col = c("blue", "lightblue", "turquoise", "darkblue", 
                      "chartreuse3", "firebrick", "orange")
) 
dev.off()

# one cell 0.3/0.1 minor
title = "Hybrid design, one cell 0.3 vs 0.1 minor"
jpeg(paste0("./plots/", title, ".jpg"), width = 700, height = 700, quality = 100)

plot(x = n_vektor, y = res_arr_hybrid["A1:B1","pvalue",,"a0b0_minor_matrix_rel"], 
     type = "l",
     xlab = "total sample size",
     ylab = "Simulated Power",
     main = title,
     ylim = c(0.3, 1),
     col = "blue")
lines(x = n_vektor, y = res_arr_hybrid["A1:B1","pvalue",,"a1b0_minor_matrix_rel"], 
      type = "l",
      col = "lightblue")
lines(x = n_vektor, y = res_arr_hybrid["A1:B1","pvalue",,"a0b1_minor_matrix_rel"], 
      type = "l",
      col = "turquoise")
lines(x = n_vektor, y = res_arr_hybrid["A1:B1","pvalue",,"a1b1_minor_matrix_rel"], 
      type = "l",
      col = "darkblue")
# add g*power
lines(x = n_vektor, y = power_GPower, 
      type = "l",
      col = "chartreuse3")
# add superpower 
lines(x = n_vektor, y = power_Superpower_hybrid*0.01, 
      type = "l",
      col = "firebrick")
# add webpower 
lines(x = n_vektor, y = power_Webpower, 
      type = "l",
      col = "orange")
# legend
legend("bottomright",
       pch = 16,
       legend = c("simulation a0b0 0.1", "simulation a1b0 0.1", 
                  "simulation a0b1 0.1", "simulation a1b1 0.1", 
                  "G*Power", "Superpower", "Webpower"),
       col = c("blue", "lightblue", "turquoise", "darkblue", 
                      "chartreuse3", "firebrick", "orange")
)  
dev.off()

# one cell 0.4/0.2 major
title = "Hybrid design, one cell 0.4 vs 0.2 major"
jpeg(paste0("./plots/", title, ".jpg"), width = 700, height = 700, quality = 100)

plot(x = n_vektor, y = res_arr_hybrid["A1:B1","pvalue",,"a0b0_major_matrix_rel"], 
     type = "l",
     xlab = "total sample size",
     ylab = "Simulated Power",
     main = title,
     ylim = c(0.3, 1),
     col = "blue")
lines(x = n_vektor, y = res_arr_hybrid["A1:B1","pvalue",,"a1b0_major_matrix_rel"], 
      type = "l",
      col = "lightblue")
lines(x = n_vektor, y = res_arr_hybrid["A1:B1","pvalue",,"a0b1_major_matrix_rel"], 
      type = "l",
      col = "turquoise")
lines(x = n_vektor, y = res_arr_hybrid["A1:B1","pvalue",,"a1b1_major_matrix_rel"], 
      type = "l",
      col = "darkblue")
# add g*power
lines(x = n_vektor, y = power_GPower, 
      type = "l",
      col = "chartreuse3")
# add superpower 
lines(x = n_vektor, y = power_Superpower_hybrid*0.01, 
      type = "l",
      col = "firebrick")
# add webpower 
lines(x = n_vektor, y = power_Webpower, 
      type = "l",
      col = "orange")
# legend
legend("bottomright",
       pch = 16,
       legend = c("simulation a0b0 0.4", "simulation a1b0 0.4", 
                  "simulation a0b1 0.4", "simulation a1b1 0.4", 
                  "G*Power", "Superpower", "Webpower"),
       col = c("blue", "lightblue", "turquoise", "darkblue", 
                      "chartreuse3", "firebrick", "orange")
)  
dev.off()

### Disordinal ----
# balanced
title = "Disordinal design, balanced"
jpeg(paste0("./plots/", title, ".jpg"), width = 700, height = 700, quality = 100)

plot(x = n_vektor, y = power_disordinal_balanced, 
     type = "l",
     xlab = "total sample size",
     ylab = "Simulated Power",
     main = title,
     ylim = c(0.3, 1),
     col = "blue")
# add g*power
lines(x = n_vektor, y = power_GPower, 
      type = "l",
      col = "chartreuse3")
# add superpower 
lines(x = n_vektor, y = power_Superpower_disordinal*0.01, 
      type = "l",
      col = "firebrick")
# add webpower 
lines(x = n_vektor, y = power_Webpower, 
      type = "l",
      col = "orange")
# legend
legend("bottomright",
       pch = 16,
       legend = c("simulation", "G*Power", "Superpower", "Webpower"),
       col = c("blue", "chartreuse3", "firebrick", "orange")
) 
dev.off()

# factors 0.6/0.4
title = "Disordinal design, factors 0.4 vs 0.6"
jpeg(paste0("./plots/", title, ".jpg"), width = 700, height = 700, quality = 100)

plot(x = n_vektor, y = res_arr_disordinal["A1:B1","pvalue",,"a0_matrix_rel"], 
     type = "l",
     xlab = "total sample size",
     ylab = "Simulated Power",
     main = title,
     ylim = c(0.3, 1),
     col = "blue")
lines(x = n_vektor, y = res_arr_disordinal["A1:B1","pvalue",,"a1_matrix_rel"], 
      type = "l",
      col = "lightblue")
lines(x = n_vektor, y = res_arr_disordinal["A1:B1","pvalue",,"b0_matrix_rel"], 
      type = "l",
      col = "turquoise")
lines(x = n_vektor, y = res_arr_disordinal["A1:B1","pvalue",,"b1_matrix_rel"], 
      type = "l",
      col = "darkblue")
# add g*power
lines(x = n_vektor, y = power_GPower, 
      type = "l",
      col = "chartreuse3")
# add superpower 
lines(x = n_vektor, y = power_Superpower_disordinal*0.01, 
      type = "l",
      col = "firebrick")
# add webpower 
lines(x = n_vektor, y = power_Webpower, 
      type = "l",
      col = "orange")
# legend
legend("bottomright",
       pch = 16,
       legend = c("simulation a0 0.6", "simulation a1 0.6", 
                  "simulation b0 0.6", "simulation b1 0.6", 
                  "G*Power", "Superpower", "Webpower"),
       col = c("blue", "lightblue", "turquoise", "darkblue", 
                      "chartreuse3", "firebrick", "orange")
)
dev.off()

# one cell 0.3/0.1 minor
title = "Disordinal design, one cell 0.3 vs 0.1 minor"
jpeg(paste0("./plots/", title, ".jpg"), width = 700, height = 700, quality = 100)

plot(x = n_vektor, y = res_arr_disordinal["A1:B1","pvalue",,"a0b0_minor_matrix_rel"], 
     type = "l",
     xlab = "total sample size",
     ylab = "Simulated Power",
     main = title,
     ylim = c(0.3, 1),
     col = "blue")
lines(x = n_vektor, y = res_arr_disordinal["A1:B1","pvalue",,"a1b0_minor_matrix_rel"], 
      type = "l",
      col = "lightblue")
lines(x = n_vektor, y = res_arr_disordinal["A1:B1","pvalue",,"a0b1_minor_matrix_rel"], 
      type = "l",
      col = "turquoise")
lines(x = n_vektor, y = res_arr_disordinal["A1:B1","pvalue",,"a1b1_minor_matrix_rel"], 
      type = "l",
      col = "darkblue")
# add g*power
lines(x = n_vektor, y = power_GPower, 
      type = "l",
      col = "chartreuse3")
# add superpower 
lines(x = n_vektor, y = power_Superpower_disordinal*0.01, 
      type = "l",
      col = "firebrick")
# add webpower 
lines(x = n_vektor, y = power_Webpower, 
      type = "l",
      col = "orange")
# legend
legend("bottomright",
       pch = 16,
       legend = c("simulation a0b0 0.1", "simulation a1b0 0.1", 
                  "simulation a0b1 0.1", "simulation a1b1 0.1", 
                  "G*Power", "Superpower", "Webpower"),
       col = c("blue", "lightblue", "turquoise", "darkblue", 
                      "chartreuse3", "firebrick", "orange")
)  
dev.off()

# one cell 0.4/0.2 major
title = "Disordinal design, one cell 0.4 vs 0.2 major"
jpeg(paste0("./plots/", title, ".jpg"), width = 700, height = 700, quality = 100)

plot(x = n_vektor, y = res_arr_disordinal["A1:B1","pvalue",,"a0b0_major_matrix_rel"], 
     type = "l",
     xlab = "total sample size",
     ylab = "Simulated Power",
     main = title,
     ylim = c(0.3, 1),
     col = "blue")
lines(x = n_vektor, y = res_arr_disordinal["A1:B1","pvalue",,"a1b0_major_matrix_rel"], 
      type = "l",
      col = "lightblue")
lines(x = n_vektor, y = res_arr_disordinal["A1:B1","pvalue",,"a0b1_major_matrix_rel"], 
      type = "l",
      col = "turquoise")
lines(x = n_vektor, y = res_arr_disordinal["A1:B1","pvalue",,"a1b1_major_matrix_rel"], 
      type = "l",
      col = "darkblue")
# add g*power
lines(x = n_vektor, y = power_GPower, 
      type = "l",
      col = "chartreuse3")
# add superpower 
lines(x = n_vektor, y = power_Superpower_disordinal*0.01, 
      type = "l",
      col = "firebrick")
# add webpower 
lines(x = n_vektor, y = power_Webpower, 
      type = "l",
      col = "orange")
# legend
legend("bottomright",
       pch = 16,
       legend = c("simulation a0b0 0.4", "simulation a1b0 0.4", 
                  "simulation a0b1 0.4", "simulation a1b1 0.4", 
                  "G*Power", "Superpower", "Webpower"),
       col = c("blue", "lightblue", "turquoise", "darkblue", 
                      "chartreuse3", "firebrick", "orange")
)  
dev.off()

## Effect size visualization per balancing ----
# plot f = sqrt(etasq/(1-etasq))
title = "Ordinal design, empirical f per balancing condition"
jpeg(paste0("./plots/", title, ".jpg"), width = 700, height = 700, quality = 100)

# Ordinal
plot(x = n_vektor, 
     y = rep(0.4, length(n_vektor)), 
     type = "l",
     xlab = "total sample size",
     ylab = "Empirical f",
     main = title,
     ylim = c(0.2, 0.6),
     col = "black")
# balanced
lines(x = n_vektor, 
      y = sqrt(res_arr_ordinal["A1:B1","pes",,"balanced_n_matrix_rel"]/(1-res_arr_disordinal["A1:B1","pes",,"balanced_n_matrix_rel"])), 
      type = "l",
      col = "blue")
# unbalanced main effect
lines(x = n_vektor, 
      y = sqrt(res_arr_ordinal["A1:B1","pes",,"a0_matrix_rel"]/(1-res_arr_disordinal["A1:B1","pes",,"a0_matrix_rel"])), 
      type = "l",
      col = "red")
lines(x = n_vektor, 
      y = sqrt(res_arr_ordinal["A1:B1","pes",,"a1_matrix_rel"]/(1-res_arr_disordinal["A1:B1","pes",,"a1_matrix_rel"])), 
      type = "l",
      col = "firebrick")
lines(x = n_vektor, 
      y = sqrt(res_arr_ordinal["A1:B1","pes",,"b0_matrix_rel"]/(1-res_arr_disordinal["A1:B1","pes",,"b0_matrix_rel"])), 
      type = "l",
      col = "orange")
lines(x = n_vektor, 
      y = sqrt(res_arr_ordinal["A1:B1","pes",,"b1_matrix_rel"]/(1-res_arr_disordinal["A1:B1","pes",,"b1_matrix_rel"])), 
      type = "l",
      col = "darkred")
# unbalanced major cell
lines(x = n_vektor, 
      y = sqrt(res_arr_ordinal["A1:B1","pes",,"a0b0_major_matrix_rel"]/(1-res_arr_disordinal["A1:B1","pes",,"a0b0_major_matrix_rel"])), 
      type = "l",
      col = "yellow")
lines(x = n_vektor, 
      y = sqrt(res_arr_ordinal["A1:B1","pes",,"a1b0_major_matrix_rel"]/(1-res_arr_disordinal["A1:B1","pes",,"a1b0_major_matrix_rel"])), 
      type = "l",
      col = "green")
lines(x = n_vektor, 
      y = sqrt(res_arr_ordinal["A1:B1","pes",,"a0b1_major_matrix_rel"]/(1-res_arr_disordinal["A1:B1","pes",,"a0b1_major_matrix_rel"])), 
      type = "l",
      col = "lightgreen")
lines(x = n_vektor, 
      y = sqrt(res_arr_ordinal["A1:B1","pes",,"a1b1_major_matrix_rel"]/(1-res_arr_disordinal["A1:B1","pes",,"a1b1_major_matrix_rel"])), 
      type = "l",
      col = "darkgreen")
# unbalanced minor cell
lines(x = n_vektor, 
      y = sqrt(res_arr_ordinal["A1:B1","pes",,"a0b0_minor_matrix_rel"]/(1-res_arr_disordinal["A1:B1","pes",,"a0b0_minor_matrix_rel"])), 
      type = "l",
      col = "purple")
lines(x = n_vektor, 
      y = sqrt(res_arr_ordinal["A1:B1","pes",,"a1b0_minor_matrix_rel"]/(1-res_arr_disordinal["A1:B1","pes",,"a1b0_minor_matrix_rel"])), 
      type = "l",
      col = "pink")
lines(x = n_vektor, 
      y = sqrt(res_arr_ordinal["A1:B1","pes",,"a0b1_minor_matrix_rel"]/(1-res_arr_disordinal["A1:B1","pes",,"a0b1_minor_matrix_rel"])), 
      type = "l",
      col = "darkmagenta")
lines(x = n_vektor, 
      y = sqrt(res_arr_ordinal["A1:B1","pes",,"a1b1_minor_matrix_rel"]/(1-res_arr_disordinal["A1:B1","pes",,"a1b1_minor_matrix_rel"])), 
      type = "l",
      col = "darkorchid")

# legend
legend("topright",
       pch = 16,
       legend = c("balanced", 
                  "a0_matrix_rel", "a1_matrix_rel", "b0_matrix_rel", "b1_matrix_rel", 
                  "a0b0_major_matrix_rel", "a1b0_major_matrix_rel", "a0b1_major_matrix_rel", "a1b1_major_matrix_rel",
                  "a0b0_minor_matrix_rel", "a1b0_minor_matrix_rel", "a0b1_minor_matrix_rel", "a1b1_minor_matrix_rel"
                  ),
       col = c("blue", 
               "red", "firebrick", "orange", "darkred",
               "yellow", "green", "lightgreen", "darkgreen",
               "purple", "pink", "darkmagenta", "darkorchid"),
       y.intersp = 0.6
)  
dev.off()

# Hybrid
title = "hybrid design, empirical f per balancing condition"
jpeg(paste0("./plots/", title, ".jpg"), width = 700, height = 700, quality = 100)

# hybrid
plot(x = n_vektor, 
     y = rep(0.4, length(n_vektor)), 
     type = "l",
     xlab = "total sample size",
     ylab = "Empirical f",
     main = title,
     ylim = c(0.2, 0.6),
     col = "black")
# balanced
lines(x = n_vektor, 
      y = sqrt(res_arr_hybrid["A1:B1","pes",,"balanced_n_matrix_rel"]/(1-res_arr_hybrid["A1:B1","pes",,"balanced_n_matrix_rel"])), 
      type = "l",
      col = "blue")
# unbalanced main effect
lines(x = n_vektor, 
      y = sqrt(res_arr_hybrid["A1:B1","pes",,"a0_matrix_rel"]/(1-res_arr_hybrid["A1:B1","pes",,"a0_matrix_rel"])), 
      type = "l",
      col = "red")
lines(x = n_vektor, 
      y = sqrt(res_arr_hybrid["A1:B1","pes",,"a1_matrix_rel"]/(1-res_arr_hybrid["A1:B1","pes",,"a1_matrix_rel"])), 
      type = "l",
      col = "firebrick")
lines(x = n_vektor, 
      y = sqrt(res_arr_hybrid["A1:B1","pes",,"b0_matrix_rel"]/(1-res_arr_hybrid["A1:B1","pes",,"b0_matrix_rel"])), 
      type = "l",
      col = "orange")
lines(x = n_vektor, 
      y = sqrt(res_arr_hybrid["A1:B1","pes",,"b1_matrix_rel"]/(1-res_arr_hybrid["A1:B1","pes",,"b1_matrix_rel"])), 
      type = "l",
      col = "darkred")
# unbalanced major cell
lines(x = n_vektor, 
      y = sqrt(res_arr_hybrid["A1:B1","pes",,"a0b0_major_matrix_rel"]/(1-res_arr_hybrid["A1:B1","pes",,"a0b0_major_matrix_rel"])), 
      type = "l",
      col = "yellow")
lines(x = n_vektor, 
      y = sqrt(res_arr_hybrid["A1:B1","pes",,"a1b0_major_matrix_rel"]/(1-res_arr_hybrid["A1:B1","pes",,"a1b0_major_matrix_rel"])), 
      type = "l",
      col = "green")
lines(x = n_vektor, 
      y = sqrt(res_arr_hybrid["A1:B1","pes",,"a0b1_major_matrix_rel"]/(1-res_arr_hybrid["A1:B1","pes",,"a0b1_major_matrix_rel"])), 
      type = "l",
      col = "lightgreen")
lines(x = n_vektor, 
      y = sqrt(res_arr_hybrid["A1:B1","pes",,"a1b1_major_matrix_rel"]/(1-res_arr_hybrid["A1:B1","pes",,"a1b1_major_matrix_rel"])), 
      type = "l",
      col = "darkgreen")
# unbalanced minor cell
lines(x = n_vektor, 
      y = sqrt(res_arr_hybrid["A1:B1","pes",,"a0b0_minor_matrix_rel"]/(1-res_arr_hybrid["A1:B1","pes",,"a0b0_minor_matrix_rel"])), 
      type = "l",
      col = "purple")
lines(x = n_vektor, 
      y = sqrt(res_arr_hybrid["A1:B1","pes",,"a1b0_minor_matrix_rel"]/(1-res_arr_hybrid["A1:B1","pes",,"a1b0_minor_matrix_rel"])), 
      type = "l",
      col = "pink")
lines(x = n_vektor, 
      y = sqrt(res_arr_hybrid["A1:B1","pes",,"a0b1_minor_matrix_rel"]/(1-res_arr_hybrid["A1:B1","pes",,"a0b1_minor_matrix_rel"])), 
      type = "l",
      col = "darkmagenta")
lines(x = n_vektor, 
      y = sqrt(res_arr_hybrid["A1:B1","pes",,"a1b1_minor_matrix_rel"]/(1-res_arr_hybrid["A1:B1","pes",,"a1b1_minor_matrix_rel"])), 
      type = "l",
      col = "darkorchid")


# legend
legend("topright",
       pch = 16,
       legend = c("balanced", 
                  "a0_matrix_rel", "a1_matrix_rel", "b0_matrix_rel", "b1_matrix_rel", 
                  "a0b0_major_matrix_rel", "a1b0_major_matrix_rel", "a0b1_major_matrix_rel", "a1b1_major_matrix_rel",
                  "a0b0_minor_matrix_rel", "a1b0_minor_matrix_rel", "a0b1_minor_matrix_rel", "a1b1_minor_matrix_rel"
       ),
       col = c("blue", 
               "red", "firebrick", "orange", "darkred",
               "yellow", "green", "lightgreen", "darkgreen",
               "purple", "pink", "darkmagenta", "darkorchid"),
       y.intersp = 0.6
)  

dev.off()

# disordinal
title = "Disordinal design, empirical f per balancing condition"
jpeg(paste0("./plots/", title, ".jpg"), width = 700, height = 700, quality = 100)

# disordinal
plot(x = n_vektor, 
     y = rep(0.4, length(n_vektor)), 
     type = "l",
     xlab = "total sample size",
     ylab = "Empirical f",
     main = title,
     ylim = c(0.2, 0.6),
     col = "black")
# balanced
lines(x = n_vektor, 
      y = sqrt(res_arr_disordinal["A1:B1","pes",,"balanced_n_matrix_rel"]/(1-res_arr_disordinal["A1:B1","pes",,"balanced_n_matrix_rel"])), 
      type = "l",
      col = "blue")
# unbalanced main effect
lines(x = n_vektor, 
      y = sqrt(res_arr_disordinal["A1:B1","pes",,"a0_matrix_rel"]/(1-res_arr_disordinal["A1:B1","pes",,"a0_matrix_rel"])), 
      type = "l",
      col = "red")
lines(x = n_vektor, 
      y = sqrt(res_arr_disordinal["A1:B1","pes",,"a1_matrix_rel"]/(1-res_arr_disordinal["A1:B1","pes",,"a1_matrix_rel"])), 
      type = "l",
      col = "firebrick")
lines(x = n_vektor, 
      y = sqrt(res_arr_disordinal["A1:B1","pes",,"b0_matrix_rel"]/(1-res_arr_disordinal["A1:B1","pes",,"b0_matrix_rel"])), 
      type = "l",
      col = "orange")
lines(x = n_vektor, 
      y = sqrt(res_arr_disordinal["A1:B1","pes",,"b1_matrix_rel"]/(1-res_arr_disordinal["A1:B1","pes",,"b1_matrix_rel"])), 
      type = "l",
      col = "darkred")
# unbalanced major cell
lines(x = n_vektor, 
      y = sqrt(res_arr_disordinal["A1:B1","pes",,"a0b0_major_matrix_rel"]/(1-res_arr_disordinal["A1:B1","pes",,"a0b0_major_matrix_rel"])), 
      type = "l",
      col = "yellow")
lines(x = n_vektor, 
      y = sqrt(res_arr_disordinal["A1:B1","pes",,"a1b0_major_matrix_rel"]/(1-res_arr_disordinal["A1:B1","pes",,"a1b0_major_matrix_rel"])), 
      type = "l",
      col = "green")
lines(x = n_vektor, 
      y = sqrt(res_arr_disordinal["A1:B1","pes",,"a0b1_major_matrix_rel"]/(1-res_arr_disordinal["A1:B1","pes",,"a0b1_major_matrix_rel"])), 
      type = "l",
      col = "lightgreen")
lines(x = n_vektor, 
      y = sqrt(res_arr_disordinal["A1:B1","pes",,"a1b1_major_matrix_rel"]/(1-res_arr_disordinal["A1:B1","pes",,"a1b1_major_matrix_rel"])), 
      type = "l",
      col = "darkgreen")
# unbalanced minor cell
lines(x = n_vektor, 
      y = sqrt(res_arr_disordinal["A1:B1","pes",,"a0b0_minor_matrix_rel"]/(1-res_arr_disordinal["A1:B1","pes",,"a0b0_minor_matrix_rel"])), 
      type = "l",
      col = "purple")
lines(x = n_vektor, 
      y = sqrt(res_arr_disordinal["A1:B1","pes",,"a1b0_minor_matrix_rel"]/(1-res_arr_disordinal["A1:B1","pes",,"a1b0_minor_matrix_rel"])), 
      type = "l",
      col = "pink")
lines(x = n_vektor, 
      y = sqrt(res_arr_disordinal["A1:B1","pes",,"a0b1_minor_matrix_rel"]/(1-res_arr_disordinal["A1:B1","pes",,"a0b1_minor_matrix_rel"])), 
      type = "l",
      col = "darkmagenta")
lines(x = n_vektor, 
      y = sqrt(res_arr_disordinal["A1:B1","pes",,"a1b1_minor_matrix_rel"]/(1-res_arr_disordinal["A1:B1","pes",,"a1b1_minor_matrix_rel"])), 
      type = "l",
      col = "darkorchid")

# legend
legend("topright",
       pch = 16,
       legend = c("balanced", 
                  "a0_matrix_rel", "a1_matrix_rel", "b0_matrix_rel", "b1_matrix_rel", 
                  "a0b0_major_matrix_rel", "a1b0_major_matrix_rel", "a0b1_major_matrix_rel", "a1b1_major_matrix_rel",
                  "a0b0_minor_matrix_rel", "a1b0_minor_matrix_rel", "a0b1_minor_matrix_rel", "a1b1_minor_matrix_rel"
       ),
       col = c("blue", 
               "red", "firebrick", "orange", "darkred",
               "yellow", "green", "lightgreen", "darkgreen",
               "purple", "pink", "darkmagenta", "darkorchid"),
       y.intersp = 0.6
)  

dev.off()
