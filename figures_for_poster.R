# Figures for Poster

# Methods figure - plot for visualizing designs ----

# Ordinal
#jpeg(paste0("./plots/", "Interaction Designs", ".jpg"), width = 1200, height = 300, quality = 100)

par(mar=c(5,5,4,4), mfrow = c(1,4))

ylim = c(-22,22)

# residual sd 12.5 ordinal 25 hybrid, 50 disordinal
res_sd_ordinal <- make_parameters_for_f(mu = mu_list[[1]], n_rel_matr = balanced_n_matrix_rel, f = f)[["sde"]]
res_sd_hybrid <- make_parameters_for_f(mu = mu_list[[2]], n_rel_matr = balanced_n_matrix_rel, f = f)[["sde"]]
res_sd_disordinal <- make_parameters_for_f(mu = mu_list[[3]], n_rel_matr = balanced_n_matrix_rel, f = f)[["sde"]]


title = "Ordinal, res sd = 12.5"

plot(x = c(-0.04,1.04),
     y = c(0,0),
     type = "l",
     lwd = 1,
     xlim = c(0,1),
     ylim = ylim,
     main = title,
     cex.main = 2, 
     xlab = "Factor B",
     ylab = "Y",
     cex.lab = 2,
     xaxt = "none", # delete x axis for custom labels
     yaxt = "none" # delete y axis for custom labels
     )
axis(1, at = c(0,1), 
     labels = c("B0", "B1"),
     cex.axis = 2)
axis(2, at = seq(-20, 20, by = 10), 
     labels = c("-20", "-10", "0", "10", "20"),
     cex.axis = 2)

points(x = c(0,1),
       y = c(0,20),
       type = "b",
       pch = 2, 
       cex = 4,
       lwd = 2,
       lty = 1
       )
points(x = c(0,1),
       y = c(0,0),
       type = "b",
       pch = 4, 
       cex = 4,
       lwd = 2,
       lty = 2
)
# pontentially add sd to plot
# arrows(x0 = 0, 
#        y0 = 0 - res_sd,
#        y1 = 0 + res_sd, 
#        angle = 90,
#        code = 3, # arrow shaft on both sides
#        length = 0.1, # length of arrow shaft
#        lwd = 2
#        )
# arrows(x0 = 0, 
#        y0 = 0 - res_sd,
#        y1 = 0 + res_sd, 
#        angle = 90,
#        code = 3, # arrow shaft on both sides
#        length = 0.1, # length of arrow shaft
#        lwd = 2
# )
# arrows(x0 = 1, 
#        y0 = 0 - res_sd,
#        y1 = 0 + res_sd, 
#        angle = 90,
#        code = 3, # arrow shaft on both sides
#        length = 0.1, # length of arrow shaft
#        lwd = 2
# )
# arrows(x0 = 1, 
#        y0 = 20 - res_sd,
#        y1 = 20 + res_sd, 
#        angle = 90,
#        code = 3, # arrow shaft on both sides
#        length = 0.1, # length of arrow shaft
#        lwd = 2
# )

# par(mar=c(5, 5, 4, 8), xpd = TRUE)

# Hybrid
title = "Hybrid, res sd = 25"

plot(x = c(-0.04,1.04),
     y = c(0,0),
     type = "l",
     lwd = 1,
     xlim = c(0,1),
     ylim = ylim,
     main = title,
     cex.main = 2, 
     xlab = "Factor B",
     ylab = "Y",
     cex.lab = 2,
     xaxt = "none", # delete x axis for custom labels
     yaxt = "none" # delete y axis for custom labels
)
axis(1, at = c(0,1), 
     labels = c("B0", "B1"),
     cex.axis = 2)
axis(2, at = seq(-20, 20, by = 10), 
     labels = c("-20", "-10", "0", "10", "20"),
     cex.axis = 2)

points(x = c(0,1),
       y = c(-20,20),
       type = "b",
       pch = 2, 
       cex = 4,
       lwd = 2,
       lty = 1
)
points(x = c(0,1),
       y = c(0,0),
       type = "b",
       pch = 4, 
       cex = 4,
       lwd = 2,
       lty = 2
)


# Disordinal
# par(mar=c(5, 5, 4, 8), xpd = TRUE)

title = "Disordinal, res sd = 50"

plot(x = c(-0.04,1.04),
     y = c(0,0),
     type = "l",
     lwd = 1,
     xlim = c(0,1),
     ylim = ylim,
     main = title,
     cex.main = 2, 
     xlab = "Factor B",
     ylab = "Y",
     cex.lab = 2,
     xaxt = "none", # delete x axis for custom labels
     yaxt = "none" # delete y axis for custom labels
)
axis(1, at = c(0,1), 
     labels = c("B0", "B1"),
     cex.axis = 2)
axis(2, at = seq(-20, 20, by = 10), 
     labels = c("-20", "-10", "0", "10", "20"),
     cex.axis = 2)

points(x = c(0,1),
       y = c(-20,20),
       type = "b",
       pch = 2, 
       cex = 4,
       lwd = 2,
       lty = 1
)
points(x = c(0,1),
       y = c(20,-20),
       type = "b",
       pch = 4, 
       cex = 4,
       lwd = 2,
       lty = 2
)

legend(legend = c("A1", "A2"),
       pch = c(2,4),
       x = 1.1, y = 5,
       cex = 2, 
       xpd = NA)


dev.off()

## Results figure: Power Curve ----

par(mar = c(5,5,3,4))

title = "Power curves Ordinal Design, Selected Conditions"
# in order of appearance
packages_and_conditions = c("G*Power (balanced)", 
                            "Superpower (balanced)", 
                            "Webpower (balanced)", 
                            "Custom Balanced", 
                            "Custom Unbalanced Factor Levels for A", 
                            "Custom Unbalanced Cells, A1B1 Major", 
                            "Custom Unbalanced Cells, A1B1 Minor")

colors_vec = c("#84160E", "#D71367", "orange",
               "#288D23", "#1C4587", "#00C9C4", "#7C2CBE")


#jpeg(paste0("./plots/", title, ".jpg"), width = 919, height = 579, quality = 100)

plot(x = n_vektor, y = res_arr_ordinal["A1:B1","pvalue",,"balanced_n_matrix_rel"], 
     type = "l",
     xlab = "total sample size",
     ylab = "Power in %",
     # main = title,
     cex.main = 1.5,
     ylim = c(0.3, 1),
     xlim = c(18, 200),
     cex.lab = 1.5,
     col = colors_vec[4],
     lwd = 3,
     xaxt = "none", # delete x axis for custom labels
     yaxt = "none" # delete y axis for custom labels
     )
axis(1, at = seq(20, 200, by = 20), 
     labels = c("20", "40", "60", "80", "100", "120", "140", "160", "180", "200"),
     cex.axis = 1.5)
axis(2, at = seq(0.3, 1, by = 0.1), 
     labels = c("30", "40", "50", "60", "70", "80", "90", "100"),
     cex.axis = 1.5)

# add unbalanced main effect: a1_matrix_rel
lines(x = n_vektor, y = res_arr_ordinal["A1:B1","pvalue",,"a1_matrix_rel"], 
      type = "l",
      lwd = 3,
      col = colors_vec[5])
# a1b1_major_matrix_rel
lines(x = n_vektor, y = res_arr_ordinal["A1:B1","pvalue",,"a1b1_major_matrix_rel"], 
      type = "l",
      lwd = 3,
      col = colors_vec[6])
# a1b1_minor_matrix_rel
lines(x = n_vektor, y = res_arr_ordinal["A1:B1","pvalue",,"a1b1_minor_matrix_rel"], 
      type = "l",
      lwd = 3,
      col = colors_vec[7])

# add g*power
lines(x = n_vektor, y = power_GPower, 
      type = "l",
      lwd = 3,
      col = colors_vec[1])
# add superpower 
lines(x = n_vektor, y = power_Superpower_ordinal*0.01, 
      type = "l",
      lwd = 3,
      col = colors_vec[2])
# add webpower 
lines(x = n_vektor, y = power_Webpower, 
      type = "l",
      lwd = 3,
      col = colors_vec[3])

# legend
legend(legend = packages_and_conditions,
       col = colors_vec,
       x = 70, y = 0.7, # position
       lty = 1, # linetype in legend
       lwd = 3, # thickness of line in legend
       pch = 15, # symbol in legend
       y.intersp = 1, # vertical space between lines
       cex = 1.5,
       bty = "n" # no box around legend
) 

#dev.off()

# Results figure: Cohens f 

title = "Cohens f Ordinal Design, Selected Conditions"
# in order of appearance
conditions = c("Custom Balanced", 
               "Custom Unbalanced Factor Levels for A", 
               "Custom Unbalanced Cells, A1B1 Major", 
               "Custom Unbalanced Cells, A1B1 Minor")
colors_vec = c("#288D23", "#1C4587", "#00C9C4", "#7C2CBE")

# jpeg(paste0("./plots/", title, ".jpg"), width = 870, height = 560, quality = 100)

plot(x = n_vektor, 
     y = rep(0.4, length(n_vektor)), # simulated f 
     type = "l",
     lwd = 2,
     xlab = "total sample size",
     ylab = "Empirical f",
     cex.lab = 1.5,
     #main = title,
     cex.main = 1.5,
     ylim = c(0.3, 0.55),
     xlim = c(20,200),
     col = "black",
     xaxt = "none", # delete x axis for custom labels
     yaxt = "none" # delete y axis for custom labels
     )
axis(1, at = seq(20, 200, by = 20), 
     labels = c("20", "40", "60", "80", "100", "120", "140", "160", "180", "200"),
     cex.axis = 1.5)
axis(2, at = seq(0.3, 0.55, by = 0.05), 
     labels = c(0.30,0.35,0.40,0.45,0.50,0.55),
     cex.axis = 1.5)

# add balanced design
lines(x = n_vektor, 
      y = sqrt(res_arr_ordinal["A1:B1","pes",,"balanced_n_matrix_rel"]/(1-res_arr_ordinal["A1:B1","pes",,"balanced_n_matrix_rel"])), 
      type = "l",
      lwd = 3,
      col = colors_vec[1])

# add unbalanced main effect: a1_matrix_rel
lines(x = n_vektor, 
      y = sqrt(res_arr_ordinal["A1:B1","pes",,"a1_matrix_rel"]/(1-res_arr_ordinal["A1:B1","pes",,"a1_matrix_rel"])), 
      type = "l",
      lwd = 3,
      col = colors_vec[2])
# a1b1_major_matrix_rel
lines(x = n_vektor, 
      y = sqrt(res_arr_ordinal["A1:B1","pes",,"a1b1_major_matrix_rel"]/(1-res_arr_ordinal["A1:B1","pes",,"a1b1_major_matrix_rel"])), 
      type = "l",
      lwd = 3,
      col = colors_vec[3])
# a1b1_minor_matrix_rel
lines(x = n_vektor, 
      y = sqrt(res_arr_ordinal["A1:B1","pes",,"a1b1_minor_matrix_rel"]/(1-res_arr_ordinal["A1:B1","pes",,"a1b1_minor_matrix_rel"])), 
      type = "l",
      lwd = 3,
      col = colors_vec[4])

# legend
legend(legend = conditions,
       col = colors_vec,
       x = 80, y = 0.55, # position
       lty = 1, # linetype in legend
       lwd = 3, # thickness of line in legend
       pch = 15, # symbol in legend
       y.intersp = 1, # vertical space between lines
       cex = 1.5,
       bty = "n" # no box around legend
) 

# dev.off()







