rm(list = ls())

# This file is used to process the fitted weibull distributions, which were producessed with Mathematica. 
library(readxl)
library(graphics)
library(scales)
library(viridis)
library(RColorBrewer)

load("PLOTS/colors.Rdata")

generate_plot_cdf <- function(datafit, scenario) {
  sum_datafit <- sum(datafit)
  histfit <- table(datafit)
  sum_histfit <- sum(histfit)
  histfit <- histfit / sum_histfit
  
  custom_breaks <- seq(0, max(datafit), 1)
  
  x_histfit <- as.integer(names(histfit))
  xx <- 0:max(x_histfit)
  yy <- rep(0, max(x_histfit) + 1)
  yy[x_histfit + 1] <- histfit
  
  weibull <- read_excel(paste0("DATA/SEXUAL_BEHAVIOR/weibull_total_",scenario,".xlsx"))
  shape <- as.numeric(names(weibull))
  scale <- weibull[[1]]
  
  cumhistfit <- cumsum(yy)
  xcont <- seq(0, max(datafit), 0.001)
  
  ints <- read_excel(paste0("DATA/SEXUAL_BEHAVIOR/intervalspartnerchangerates_total_",scenario,".xlsx"))
  ints <- as.numeric(names(ints))
  ints <- ints[-c(1,length(ints))]
  
  size_plot=6
  rescale=0.6666
  if(scenario == 0){
    panel <- "A"
  }else if(scenario == 1){
    panel <- "B"
  }else if(scenario == 2){
    panel <- "C"
  }
  pdf(file = paste0("PLOTS/precure/FigS2",panel,".pdf"), width = size_plot, height = size_plot * rescale*0.82)
  par(mar = c(4,4,.5,.5))
  par(mgp = c(2,1,0))
  plot(0, xlim = c(0, 40), ylim = c(0, 1), col = "white", xlab = "Partner change rate (1/six months)", ylab = "Cumulative density function", xaxt = "n", yaxt = "n", bty = 'n', cex.lab=1.25)
  axis(1, at = seq(0, 50, 5), cex.axis=1.25, col = "gray72", las=1, pos=0)
  axis(2, at = seq(0, 1, 0.2), cex.axis=1.25, col = "gray72", las=2, pos=0)
  abline(v = ints, col = adjustcolor(colrom3, alpha.f = 1), lty = 5, lwd=1.)
  abline(h = seq(0, 1, 0.2), col = "gray72", lty = 3)
  lines(xcont, pweibull(xcont, shape = shape, scale = scale), col = colnew2, lwd = 2)
  points(custom_breaks, cumhistfit, col = colrom1, pch = 1, cex = 1.25, lwd=2)
  legend(x = "bottomright", legend = c("Data", paste0("Weibull fit: shape = ",round(shape,2),", scale = ", round(scale,2) )), col = c(colrom1, colnew2), pch = c(1, NA), lty = c(NA, 1), lwd=c(2,2), cex = 1, bty = "n")
  dev.off()
}

load("DATA/SEXUAL_BEHAVIOR/outdfs_total.Rdata")

generate_plot_cdf(sdf$outdf_total_0,0)
generate_plot_cdf(sdf$outdf_total_1,1)
generate_plot_cdf(sdf$outdf_total_2,2)
