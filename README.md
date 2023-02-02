# check_pwr_interaction

Project for ARMS course - Psychology - University of Kassel (2023)
Authors: Salomé Li Keintzel, Frauke Langhein and Paula G. Rosendahl 

Check correctness of power tools G*Power, Webpower and Superpower via custom simulation of 2x2 between-subjects ANOVA interaction effect

simulation conditions: 
- different interaction designs (ordinal, hybrid, disordinal)
- different balancing of n (see matrices.R)
- different sample sizes (for plotting power curve: 20:260) 

# Files
`simulation.R`: main R Script with custom simulation and plots:  
`functions.R`: custom functions used  
`matrices.R`: Script defining different balancing matrices  
`tools.R`: results from tested tools  
`figures_for_poster.R`: code for figures on poster  
`Cohens_f_sde.md`: Explanation of how Cohens f can be calculated from mean structure and residual standard deviation (and vice versa) for our simulation  

# Data
`... power_Superpower... .Rdata`: Power curve results for Superpower (generated in tools.R)

# Folders
`/plots`: empty folder to save plots into

