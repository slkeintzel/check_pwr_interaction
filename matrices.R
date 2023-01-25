# Unbalanced Matrices for 
# R Script for Power Simulation for interaction effect in 2 x 2 between subjects ANOVA Design 
# Authors: Salomé Li Keintzel, Frauke Langhein and Paula G. Rosendahl (2023)

# balanced relative n matrix
balanced_n_matrix_rel = matrix(c(0.25,0.25,0.5,0.25,0.25,0.5,0.5,0.5,1), 
                               nrow = 3, byrow = T, 
                               dimnames = list(c("C_b = 0", "C_b = 1", "tot"), 
                                               c("C_a = 0", "C_a = 1", "tot")))


#### main effect for A unequal

#A0 larger (0.6 vs. 0.4)
a0_matrix_rel <- matrix(c(0.3,0.2,0.5,
                          0.3,0.2,0.5,
                          0.6,0.4,1), 
                        nrow = 3, byrow = T, 
                        dimnames = list(c("C_b = 0", "C_b = 1", "tot"), 
                                        c("C_a = 0", "C_a = 1", "tot")))

#A1 larger
a1_matrix_rel <- matrix(c(0.2,0.3,0.5,
                          0.2,0.3,0.5,
                          0.4,0.6,1), 
                        nrow = 3, byrow = T, 
                        dimnames = list(c("C_b = 0", "C_b = 1", "tot"), 
                                        c("C_a = 0", "C_a = 1", "tot")))

#### main effect for B unequal 

#B0 larger (0.6 vs. 0.4)
b0_matrix_rel <- matrix(c(0.3,0.3,0.6,
                          0.2,0.2,0.4,
                          0.5,0.5,1),  
                        nrow = 3, byrow = T, 
                        dimnames = list(c("C_b = 0", "C_b = 1", "tot"), 
                                        c("C_a = 0", "C_a = 1", "tot")))

#B1 larger
b1_matrix_rel <- matrix(c(0.2,0.2,0.4, 
                          0.3,0.3,0.6,
                          0.5,0.5,1), 
                        nrow = 3, byrow = T, 
                        dimnames = list(c("C_b = 0", "C_b = 1", "tot"), 
                                        c("C_a = 0", "C_a = 1", "tot")))

####Interaction A:B - each cell larger once

#A0B0 larger (0.4 vs. 0.2)
a0b0_major_matrix_rel <- matrix(c(0.4,0.2,0.6,
                                  0.2,0.2,0.4,
                                  0.6,0.4,1), 
                                nrow = 3, byrow = T, 
                                dimnames = list(c("C_b = 0", "C_b = 1", "tot"), 
                                                c("C_a = 0", "C_a = 1", "tot")))

#A1B0 larger (0.4 vs. 0.2)
a1b0_major_matrix_rel <- matrix(c(0.2,0.4,0.6,
                                  0.2,0.2,0.4,
                                  0.4,0.6,1),
                                nrow = 3, byrow = T, 
                                dimnames = list(c("C_b = 0", "C_b = 1", "tot"), 
                                                c("C_a = 0", "C_a = 1", "tot")))

#A0B1 larger (0.4 vs. 0.2)
a0b1_major_matrix_rel <- matrix(c(0.2,0.2,0.4,
                                  0.4,0.2,0.6,
                                  0.6,0.4,1),
                                nrow = 3, byrow = T, 
                                dimnames = list(c("C_b = 0", "C_b = 1", "tot"), 
                                                c("C_a = 0", "C_a = 1", "tot")))

#A1B1 larger (0.4 vs. 0.2)
a1b1_major_matrix_rel <- matrix(c(0.2,0.2,0.4,
                                  0.2,0.4,0.6,
                                  0.4,0.6,1),
                                nrow = 3, byrow = T, 
                                dimnames = list(c("C_b = 0", "C_b = 1", "tot"), 
                                                c("C_a = 0", "C_a = 1", "tot")))

####Interaction A:B - each cell smaller once

#A0B0 smaller (0.1 vs. 0.3)
a0b0_minor_matrix_rel <- matrix(c(0.1,0.3,0.4,0.3,0.3,0.6,0.4,0.6,1), 
                                nrow = 3, byrow = T, 
                                dimnames = list(c("C_b = 0", "C_b = 1", "tot"), 
                                                c("C_a = 0", "C_a = 1", "tot")))

#A1B0 smaller (0.1 vs. 0.3)
a1b0_minor_matrix_rel <- matrix(c(0.3,0.1,0.4,0.3,0.3,0.6,0.6,0.4,1), 
                                nrow = 3, byrow = T, 
                                dimnames = list(c("C_b = 0", "C_b = 1", "tot"), 
                                                c("C_a = 0", "C_a = 1", "tot")))

#A0B1 smaller (0.1 vs. 0.3)
a0b1_minor_matrix_rel <- matrix(c(0.3,0.3,0.6,0.1,0.3,0.4,0.4,0.6,1), 
                                nrow = 3, byrow = T, 
                                dimnames = list(c("C_b = 0", "C_b = 1", "tot"), 
                                                c("C_a = 0", "C_a = 1", "tot")))

#A1B1 smaller (0.1 vs. 0.3)
a1b1_minor_matrix_rel <- matrix(c(0.3,0.3,0.6,0.3,0.1,0.4,0.6,0.4,1), 
                                nrow = 3, byrow = T, 
                                dimnames = list(c("C_b = 0", "C_b = 1", "tot"), 
                                                c("C_a = 0", "C_a = 1", "tot")))

# list all matrices
rel_matrices <- list(balanced_n_matrix_rel = balanced_n_matrix_rel,
                     a0_matrix_rel = a0_matrix_rel, 
                     a1_matrix_rel = a1_matrix_rel, 
                     b0_matrix_rel = b0_matrix_rel, 
                     b1_matrix_rel = b1_matrix_rel,
                     a0b0_major_matrix_rel = a0b0_major_matrix_rel, 
                     a1b0_major_matrix_rel = a1b0_major_matrix_rel, 
                     a0b1_major_matrix_rel = a0b1_major_matrix_rel, 
                     a1b1_major_matrix_rel = a1b1_major_matrix_rel,
                     a0b0_minor_matrix_rel = a0b0_minor_matrix_rel, 
                     a1b0_minor_matrix_rel = a1b0_minor_matrix_rel,
                     a0b1_minor_matrix_rel = a0b1_minor_matrix_rel, 
                     a1b1_minor_matrix_rel = a1b1_minor_matrix_rel)

