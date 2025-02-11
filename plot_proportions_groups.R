rm(list = ls())

load("PLOTS/colors.Rdata")
load("RESULTS/check_proportions_groups.Rdata")

# Plot the prevalence by risk group
y_end_plot <- 2024
pdf(file = ("PLOTS/precure/FigS7.pdf"), width = size_plot, height = size_plot * rescale*0.82)
par(mar = c(4,4,.5,6), xpd=T)
par(mgp = c(3,1,0))
step <- 0.2
y0 <- 0
y1 <- 1
plot(0, xlim = c(0, (y_end_plot-years_calibration[first_year_fit_idx])*dt), ylim = c(0, y1 + (y1 - y0) * 0.), col = "white", xlab = "Year", ylab = "Proportion of cases", xaxt = "n", yaxt = "n", bty = "n", cex.lab=1.25)
axis(1, at = seq(0, (y_end_plot-years_calibration[first_year_fit_idx])*dt , dt), labels = NA, col = "gray72", pos=0)
axis(2, at = (seq(0, 1, step)), cex.axis=1.25, col = "gray72", las=2, pos=0)
text(x = seq(1, (y_end_plot-years_calibration[first_year_fit_idx]+1)*dt , dt), y = par("usr")[3] - 0.125 * diff(par("usr")[3:4]),
     labels = seq(years_calibration[first_year_fit_idx], y_end_plot, 1),
     srt = 45, adj = c(0.5, 0.5), xpd = TRUE, cex = 1.25)
for (i in 1:((y_end_plot-years_calibration[first_year_fit_idx])*dt)) {
  vec <- 0
  rect(xleft =  i - 1., xright = i + 0, y0, mean_prop_inf1[i] , col = adjustcolor(colbarprop[1], alpha.f = 0.8), border = NA)
  vec <- vec + mean_prop_inf1[i]
  rect(xleft = i - 1, xright = i + 0, y0+vec, mean_prop_inf2[i]+vec , col = adjustcolor(colbarprop[2], alpha.f = 0.8), border = NA)
  vec <- vec + mean_prop_inf2[i]
  rect(xleft = i - 1, xright = i + 0, y0+vec, mean_prop_inf3[i]+vec , col = adjustcolor(colbarprop[3], alpha.f = 0.8), border = NA)
  vec <- vec + mean_prop_inf3[i]
  rect(xleft = i - 1, xright = i + 0, y0+vec, mean_prop_inf4[i]+vec , col = adjustcolor(colbarprop[4], alpha.f = 0.8), border = NA)
}
legend(x = (y_end_plot-years_calibration[first_year_fit_idx])*dt + 1, y = y1, legend = c("1, low", "2, mid-low", "3, mid-high", "4, high"), col = c(colbarprop[1], colbarprop[2], colbarprop[3], colbarprop[4]), pch = rep(15,4), title = expression(bold("Risk group")), pt.cex = rep(2,4), cex = 1, bty = "n")
dev.off()
