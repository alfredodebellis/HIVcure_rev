rm(list = ls())

library(RColorBrewer)
library("viridis") 
library(ggplot2)
library(patchwork)
library(ppcor)
library(patchwork)

colors_inf <- c(rgb(194, 170, 202, maxColorValue = 255), # Purple
                rgb(224, 172, 180, maxColorValue = 255), # Pink
                rgb(162, 0, 37, maxColorValue = 255), # Red
                rgb(220, 100, 50, maxColorValue = 255)) # Orange
colors_inf <- colors_inf[c(1,3,4,2)]

colrom2 <- rgb(241/255, 189/255, 64/255)
colrom1 <- rgb(140/255, 31/255, 47/255)
colrom3 <- "#11315D"
colbar2 <- rev(brewer.pal(11, "RdBu"))
colnew0 <- colbar2[10]
colnew1 <- colbar2[8]
colnew2 <- colbar2[3]
colnew3 <- colbar2[4]
colbarprop <- rev(brewer.pal(4,"RdBu"))
colcurr <-"#F4A582"
  
colors1 <- brewer.pal(n = 9, name = 'Blues')[c(4,6,9)]
colors2 <- rev(brewer.pal(n = 9, name = 'Greens')[c(4,6,9)])
colors2 <- c(colcurr,colors2)
colors2_reduct <- colors2[-1]

colcol <- rev(brewer.pal(11, "RdBu"))
colcol_reb <- (brewer.pal(11, "RdPu"))
colcol_new <- rev(brewer.pal(11, "PuOr"))

size_plot=6
rescale=0.6666

infty <- 10000000

save.image("PLOTS/colors.RData")
