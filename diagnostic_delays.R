rm(list = ls())

library(readxl)
library(dplyr)
# This code estimates the diagnostic delays distributions for acute and chronic stages

# Number of cases
Z <- 1000

# Number of sampled set of diagnostics delays
n_samp_sub <- 10000

# Generate random sets from uniform priors
mint1 <- 1/12
mint2 <- 1/12
maxt1 <- 2
maxt2 <- 1
tau1 <- runif(n_samp_sub, min = mint1, max = maxt1)
tau2 <- runif(n_samp_sub, min = mint2, max = maxt2)
tau3 <- rep(12,n_samp_sub)
taus <- cbind(tau1, tau2, tau3)

# Importing diagnoses proportion data from 2015 to 2022
dat_prop <- read_excel("DATA/EPIDEMIOLOGICAL/HIV_epidemiological_data.xlsx", sheet = "SHM", col_names = F)
dat_prop <- as.data.frame(dat_prop)
dat_prop <- t(dat_prop)
colnames(dat_prop) <- dat_prop[1,]
idxdiagprop1 <- which(colnames(dat_prop) == "diagnosed in early stage I (6 month before negative), proportion")
idxdiagprop2 <- which(colnames(dat_prop) == "diagnosed in early stage II (12 month before negative), proportion")
dat_prop <- dat_prop[-1,]
dat_prop <- as.data.frame(dat_prop)
dat_prop <- dat_prop[,c(1,idxdiagprop1,idxdiagprop2)]
dat_prop <- dat_prop %>%
  mutate(across(everything(), ~ as.numeric(.) %>% round(3)))
dat_prop <- dat_prop[nrow(dat_prop):1, ]
lastpropfit <- which(dat_prop$Year == 2019)

# Importing data, from 2015 to 2019 (pre-COVID-19), on the average proportion of diagnoses before 6 months, between 6 and 12 months and more than a year
prop_6m <- mean(dat_prop$`diagnosed in early stage I (6 month before negative), proportion`[1:lastpropfit])
prop_12m <- mean(dat_prop$`diagnosed in early stage II (12 month before negative), proportion`[1:lastpropfit]) - prop_6m
prop_long <- 1 - prop_12m - prop_6m
hist_data <- c(prop_6m,prop_12m,prop_long)

progression <- c(1/0.142, 1/8.439, 1/1.184, 1/1.316) #see Table 1

ddnew <- c()
dd <- c()
dd_tmp <- rep(NA, n_samp_sub*Z)

count_na <- rep(0, n_samp_sub)

accp_submodel <- c()
accp_submodel_dd <- c()
diff_all <- rep(NA,3)

# Threshold of error for acceptance of parameters
thresh_sub <- 0.1

dd_95 <- rep(0,n_samp_sub)

count <- 1
count1 <- 1
count2 <- 1
count3 <- 1
# Stochastic model of progression of the virus/diagnoses
for(i in 1:n_samp_sub){
  print(i)
  dd_tmp <- rep(NA,Z)
  dd_tmp2 <- rep(NA,Z)
  # dd_tmp3 <- rep(NA,Z)
  for (j in 1:Z) {
    tp1 <- rexp(1, rate = progression[1])
    tt1 <- rexp(1, rate = tau1[i])
    if(tp1 < tt1){
      tp2 <- rexp(1, rate = progression[2])
      tt2 <- rexp(1, rate = tau2[i])
      if(tp2 < tt2){
        tp3 <- rexp(1, rate = progression[3])
        tt3 <- rexp(1, rate = tau3[i])
        if(tp3 < tt3){
          tp4 <- rexp(1, rate = progression[4])
          tt4 <- rexp(1, rate = tau3[i])
          if(tp4 < tt4){
            count_na[i] <- count_na[i] + 1
            next
          }else{
            u <- tp1 + tp2 + tp3 + tt4
            dd_tmp[j] <- u
            count <- count + 1
            count3 <- count3 + 1
          }
        }else{
          u <- tp1 + tp2 + tt3
          dd_tmp[j] <- u
          count <- count + 1
          count3 <- count3 + 1
        }
      }else{
        u <- tp1 + tt2
        dd_tmp[j] <- u
        count <- count + 1
        count2 <- count2 + 1
      }
    }else{
      u <- tt1
      dd_tmp[j] <- u
      count <- count + 1
      count1 <- count1 + 1
    }
  }
  dd_tmp2 <- dd_tmp
  dd_tmp <- dd_tmp[which(!is.na(dd_tmp))]
  dd_95[i] <- mean(dd_tmp)
  
  hist_man <- c(length(which(dd_tmp<=1/2)),
                length(which(dd_tmp>1/2 & dd_tmp<=1)),
                length(which(dd_tmp>1))
  )
  hist_man <- hist_man / sum(hist_man)
  diff_all <- abs(hist_man - hist_data) / hist_data
  if(all(diff_all < thresh_sub)){
    accp_submodel <- rbind(accp_submodel, taus[i,])
    accp_submodel_dd <- rbind(accp_submodel_dd, dd_95[i])
    dd <- c(dd,dd_tmp)
    ddnew <- rbind(ddnew,dd_tmp2)
  }
}

mean(accp_submodel_dd) * 12
quantile(accp_submodel_dd,0.025) * 12
quantile(accp_submodel_dd,0.975) * 12

# Select only sets for which tau1 > tau2 and recalculate the number of accepted sets
accp_submodel <- as.data.frame(accp_submodel)
accp_submodel <- accp_submodel[which(accp_submodel$tau1 > accp_submodel$tau2),]
len_post_sub <- nrow(accp_submodel)

# Diagnostic delays of each time interval
ddfirst <- rep(NA,len_post_sub)
ddsecond <- rep(NA,len_post_sub)
ddthird <- rep(NA,len_post_sub)
for (i in 1:len_post_sub) {
  ddfirst[i] <- length(which(!is.na(ddnew[i,]) & ddnew[i,]<=1/2)) / (Z-count_na[i])
  ddsecond[i] <- length(which(!is.na(ddnew[i,]) & ddnew[i,]>1/2 & ddnew[i,]<=1)) / (Z-count_na[i])
  ddthird[i] <- length(which(!is.na(ddnew[i,]) & ddnew[i,]>1)) / (Z-count_na[i])
}

load("PLOTS/colors.Rdata")

# Plot diagnostic delay for each time interval
size_plot <- 6
rescale <- 0.66666
addresc <- 0.8
pdf(file = "PLOTS/precure/FigS3.pdf", width = size_plot, height = size_plot * rescale)
par(mar = c(6.5,4.5,.5,.5))
par(mgp = c(3,1,0))
vec <- dd
hist_man <- c(length(which(dd<=1/2)),
              length(which(dd>1/2 & dd<=1)),
              length(which(dd>1))
)
hist_man <- hist_man / sum(hist_man)
step2 <- 0.2
y1 <- max(hist_man)
plot(0, xlim = c(1-0.5, 3+0.5), ylim = c(0, y1+step2), col = "white", xlab = "", ylab = "Frequency", xaxt = "n", yaxt = "n", bty = "n")
for (i in 1:3) {
  rect(xleft = i - 0.4, xright = i , 0, hist_data[i] , col = adjustcolor(colnew1,1), border = NA)
}
rect(xleft = 1, xright = 1 + 0.4 , 0, mean(ddfirst) , col = adjustcolor(colnew3,1), border = NA)
rect(xleft = 2, xright = 2 + 0.4 , 0, mean(ddsecond) , col = adjustcolor(colnew3,1), border = NA)
rect(xleft = 3, xright = 3 + 0.4 , 0, mean(ddthird) , col = adjustcolor(colnew3,1), border = NA)
segments(x0=1+0.2, x1=1+0.2, y0=quantile(ddfirst,0.025), y1=quantile(ddfirst,0.975), col = adjustcolor(colrom3,alpha.f = 1), lwd=2.)
segments(x0=2+0.2, x1=2+0.2, y0=quantile(ddsecond,0.025), y1=quantile(ddsecond,0.975), col = adjustcolor(colrom3,alpha.f = 1), lwd=2.)
segments(x0=3+0.2, x1=3+0.2, y0=quantile(ddthird,0.025), y1=quantile(ddthird,0.975), col = adjustcolor(colrom3,alpha.f = 1), lwd=2.)
axis(1, at = 1:3, labels = NA, col="white", las=2, cex.axis=1.25)
axis(2, at = seq(0,y1+step2,step2), cex.axis=1, las=2, col="gray72")
text(x = 1:3, y = par("usr")[3] - 0.15 * diff(par("usr")[3:4]),
     labels = c("0-6 months", "6-12 months", "> 1 year"),
     srt = 45, adj = c(0.5, 0.5), xpd = TRUE, cex = 1)
text(x = 3/2 + 1/2, y = par("usr")[3] - (0.4)* diff(par("usr")[3:4]),
     labels = "Average (over 2017-2022) proportion of diagnoses",
     srt = 0, adj = c(0.5, 0.5), xpd = TRUE, cex = 1)
legend("topleft", c("Data", "Stochastic submodel"), pch = rep(15,3), col = c(colnew1, colnew3), pt.cex = rep(2,3), cex = 1, bty = "n" )
dev.off()

# Save accepted sets of diagnosis rates
outres <- accp_submodel[,1:2]
save(outres, file = ("DATA/EPIDEMIOLOGICAL/taus12.RData"))

# Plot diagnosis rates distribution reconstructed
pdf(file = "PLOTS/precure/FigS4.pdf", width = 2*size_plot, height = size_plot * rescale)
layout.matrix<-matrix(c(1,2), nrow = 1, ncol = 2)
layout(mat = layout.matrix)
par(mar = c(4.5*(2)-2-1,4.5*(2)-2-1,3,0.5))
layout(mat = layout.matrix, 
       heights = c(1,1),widths = c(1,1))
par(mgp = c(4, 1, 0))
vec <- outres[,1]
step <- 0.1
y0 <- 0
y1 <- maxt1
bins <- seq(y0,y1+step,step)
hist(vec, breaks=bins, col = colnew3, border = "white", xlab="", ylab="", xaxt="n", yaxt="n", bty='n', main = NULL)
hist_data <- hist(vec, breaks = bins, plot = FALSE)
counts <- hist_data$counts
max_counts <- max(counts)
norm_counts <- counts / sum(counts)
max_freq <- max(norm_counts)
bins <- seq(y0,y1+step*4,step*4)
axis(1, at = bins, col="gray72", las=1, cex.axis=1.25, pos=0)
axis(2, at = seq(0,max_counts+step*0.3*sum(counts),step*0.3*sum(counts)), label=seq(0,max_freq+step*0.3,step*0.3), cex.axis=1.25, las=2, col="gray72", pos=0)
mtext(expression(tau[1]), side = 1, line = 3.5, cex = 1.75)
mtext("Frequency", side = 2, line = 3.5, cex = 1.75)
abline(v=mean(vec), col=colrom3, lwd=2)
abline(v=quantile(vec,0.025), col=colrom3, lty=2, lwd=2)
abline(v=quantile(vec,0.975), col=colrom3, lty=2, lwd=2)
mtext("a", side = 3, line = 1.5, adj = 0.05, cex = 1.5, font = 2)
vec <- outres[,2]
step <- 0.05
y0 <- 0
y1 <- maxt2
bins <- seq(y0,y1+step,step)
hist(vec, breaks=bins, col = colnew3, border = "white", xlab="", ylab="", xaxt="n", yaxt="n", bty='n', main = NULL)
hist_data <- hist(vec, breaks = bins, plot = FALSE)
counts <- hist_data$counts
max_counts <- max(counts)
norm_counts <- counts / sum(counts)
max_freq <- max(norm_counts)
bins <- seq(y0,y1+step*4,step*4)
axis(1, at = bins, col="gray72", las=1, cex.axis=1.25, pos=0)
axis(2, at = seq(0,max_counts+step*2*sum(counts),step*2*sum(counts)), label=seq(0,max_freq+step*2,step*2), cex.axis=1.25, las=2, col="gray72", pos=0)
mtext(expression(tau[2]), side = 1, line = 3.5, cex = 1.75)
mtext("Frequency", side = 2, line = 3.5, cex = 1.75)
abline(v=mean(vec), col=colrom3, lwd=2)
abline(v=quantile(vec,0.025), col=colrom3, lty=2, lwd=2)
abline(v=quantile(vec,0.975), col=colrom3, lty=2, lwd=2)
mtext("b", side = 3, line = 1.5, adj = 0.05, cex = 1.5, font = 2)
dev.off()

