rm(list = ls())

library(ppcor)
library(RColorBrewer)
library("viridis") 
library(ggplot2)

load("RESULTS/fit_threshold15.RData")
load("PLOTS/colors.Rdata")

# bests <- 100
bests <- nrow(accp_list)

mean(counter_fit[order(counter_fit)[1:bests]])
which(counter_fit[order(counter_fit)]>0.3)[1]

accp_list <- accp_list[order(counter_fit)[1:bests],]

counter_incidence <- counter_incidence[order(counter_fit)[1:bests],]
counter_undiag <- counter_undiag[order(counter_fit)[1:bests],]
counter_newdiag <- counter_newdiag[order(counter_fit)[1:bests],]
counter_prevalence <- counter_prevalence[order(counter_fit)[1:bests],]
counter_art <- counter_art[order(counter_fit)[1:bests],]
counter_prep <- counter_prep[order(counter_fit)[1:bests],]
counter_pop <- counter_pop[order(counter_fit)[1:bests],]
counter_diag <- counter_diag[order(counter_fit)[1:bests],]
counter_impdiag <- counter_impdiag[order(counter_fit)[1:bests],]
counter_impundiag <- counter_impundiag[order(counter_fit)[1:bests],]
counter_newtreat <- counter_newtreat[order(counter_fit)[1:bests],]

# Define indexes and boundaries (years) of the plots
y_start_plot <- years_fit[1]
y_end_plot <- y_end
idx_end_plot <- y_end_plot - years_calibration[1] + 1
idxs_years_plot <- seq(dt*(first_year_fit_idx-3)+1,dt*(last_year_fit_idx),by=dt)
idxs_years_plot_2022 <- idxs_years_plot[length(idxs_years_plot)] + dt
yf_idxs_plus <- seq(dt*(first_year_fit_idx-1)+1,dt*(last_year_fit_idx+1),by=dt)

# Rescale of the incidences
perinc <- parms["NN"] / 1e5

# Calculate means and CrIs of the observables
mean_undiag <- colMeans(counter_undiag)
ci_undiag <- matrix(0,2,length(mean_undiag))
for (i in 1:length(mean_undiag)) {
  ci_undiag[1,i] <- quantile(counter_undiag[,i], 0.975)
  ci_undiag[2,i] <- quantile(counter_undiag[,i], 0.025)
}
mean_incidence <- colMeans(counter_incidence)
ci_incidence <- matrix(0,2,length(mean_incidence))
for (i in 1:length(mean_incidence)) {
  ci_incidence[1,i] <- quantile(counter_incidence[,i], 0.975)
  ci_incidence[2,i] <- quantile(counter_incidence[,i], 0.025)
}
counter_art2 <- counter_art / counter_prevalence
mean_art2 <- colMeans(counter_art2)
mean_art <- colMeans(counter_art)
ci_art <- matrix(0,2,length(mean_art))
ci_art2 <- matrix(0,2,length(mean_art2))
for (i in 1:length(mean_art)) {
  ci_art[1,i] <- quantile(counter_art[,i], 0.975)
  ci_art[2,i] <- quantile(counter_art[,i], 0.025)
  ci_art2[1,i] <- quantile(counter_art2[,i], 0.975)
  ci_art2[2,i] <- quantile(counter_art2[,i], 0.025)
}
mean_prev <- colMeans(counter_prevalence)
ci_prev <- matrix(0,2,length(mean_prev))
for (i in 1:length(mean_prev)) {
  ci_prev[1,i] <- quantile(counter_prevalence[,i], 0.975)
  ci_prev[2,i] <- quantile(counter_prevalence[,i], 0.025)
}
mean_prep <- colMeans(counter_prep)
ci_prep <- matrix(0,2,length(mean_prep))
for (i in 1:length(mean_prep)) {
  ci_prep[1,i] <- quantile(counter_prep[,i], 0.975)
  ci_prep[2,i] <- quantile(counter_prep[,i], 0.025)
}
mean_impdiag <- colMeans(counter_impdiag)
median_impdiag <- apply(counter_impdiag, 2, median)
ci_impdiag <- matrix(0,2,length(mean_impdiag))
for (i in 1:length(mean_impdiag)) {
  ci_impdiag[1,i] <- quantile(counter_impdiag[,i], 0.975)
  ci_impdiag[2,i] <- quantile(counter_impdiag[,i], 0.025)
}
mean_newdiag <- colMeans(counter_newdiag)
ci_newdiag <- matrix(0,2,length(mean_newdiag))
for (i in 1:length(mean_newdiag)) {
  ci_newdiag[1,i] <- quantile(counter_newdiag[,i], 0.975)
  ci_newdiag[2,i] <- quantile(counter_newdiag[,i], 0.025)
}

#Plot distributions of accepted parameters
rescale2 <- 0.75
n_columns <- 3
n_rows <- 3
cex_scale <- n_columns
pdf(file = "PLOTS/precure/FigS5.pdf", width = n_columns*size_plot, height = n_rows*size_plot * rescale2)
layout.matrix<-matrix(c(1,4,7,2,5,8,3,6,9), nrow = 3, ncol = 3)
layout(mat = layout.matrix)
par(mar = c(5.5*(cex_scale)-cex_scale-1,5.5*(cex_scale)-cex_scale-1,6,1.5))
layout(mat = layout.matrix, 
       heights = c(1,1),widths = c(1,1))
par(mgp = c(4, 1, 0))
vec <- accp_list$lambda
bin1 <- 6
bin2 <- 10
step <- 0.0075
y0 <- parms_fit_ranges[1,1]
y0 <- 0
y1 <- parms_fit_ranges[1,2]
bins <- seq(y0,y1+step,step)
hist(vec, breaks = bins, col = colnew3, border = "white", xlab="", ylab="", xaxt="n", yaxt="n", bty='n', main = NULL)
hist_data <- hist(vec, breaks = bins, plot = FALSE)
counts <- hist_data$counts
max_counts <- max(counts)
norm_counts <- counts / sum(counts)
max_freq <- max(norm_counts)
bins <- seq(y0,y1+step*bin1,step*bin1)
bins <- seq(0,y1,0.05)
axis(1, at = bins, col="gray72", las=2, cex.axis=cex_scale, pos=0)
axis(2, at = seq(0,max_counts+step*bin2*sum(counts),100), label=round(seq(0,max_counts+step*bin2*sum(counts),100)/sum(counts),2), cex.axis=cex_scale, las=2, col="gray72", pos=0)
mtext(expression(lambda), side = 1, line = 10, cex = 2.5)
mtext("Frequency", side = 2, line = 8.5, cex = 2.5)
abline(v=mean(vec), col=colrom3, lwd=2)
abline(v=quantile(vec,0.025), col=colrom3, lty=2, lwd=2)
abline(v=quantile(vec,0.975), col=colrom3, lty=2, lwd=2)
mtext("a", side = 3, line = 3, adj = 0.05, cex = cex_scale-1, font = 2)
vec <- accp_list$omega
bin1 <- 5
bin2 <- 0.5
step <- 0.1
y0 <- parms_fit_ranges[2,1]
y1 <- parms_fit_ranges[2,2]
bins <- seq(y0,y1+step,step)
hist(vec, breaks = bins, col = colnew3, border = "white", xlab="", ylab="", xaxt="n", yaxt="n", bty='n', main = NULL)
hist_data <- hist(vec, breaks = bins, plot = FALSE)
counts <- hist_data$counts
max_counts <- max(counts)
norm_counts <- counts / sum(counts)
max_freq <- max(norm_counts)
bins <- seq(y0,y1+step*bin1,step*bin1)
bins <- seq(y0,y1+step*bin1,step*bin1)
axis(1, at = bins, col="gray72", las=2, cex.axis=cex_scale, pos=0)
axis(2, at = seq(0,max_counts+step*bin2*sum(counts),step*bin2*sum(counts)), label=seq(0,max_freq+step*bin2,step*bin2), cex.axis=cex_scale, las=2, col="gray72", pos=0)
mtext(expression(omega), side = 1, line = 10, cex = 2.5)
mtext("Frequency", side = 2, line = 8.5, cex = 2.5)
abline(v=mean(vec), col=colrom3, lwd=2)
abline(v=quantile(vec,0.025), col=colrom3, lty=2, lwd=2)
abline(v=quantile(vec,0.975), col=colrom3, lty=2, lwd=2)
mtext("b", side = 3, line = 3, adj = 0.05, cex = cex_scale-1, font = 2)
vec <- accp_list$eps
bin1 <- 5
bin2 <- 10
step <- 0.005
y0 <- parms_fit_ranges[3,1]
y1 <- parms_fit_ranges[3,2]
y0 <- 0
bins <- seq(y0,y1+step,step)
hist(vec, breaks = bins, col = colnew3, border = "white", xlab="", ylab="", xaxt="n", yaxt="n", bty='n', main = NULL)
hist_data <- hist(vec, breaks = bins, plot = FALSE)
counts <- hist_data$counts
max_counts <- max(counts)
norm_counts <- counts / sum(counts)
max_freq <- max(norm_counts)
bins <- seq(y0,y1+step*bin1,step*bin1)
axis(1, at = bins, col="gray72", las=2, cex.axis=cex_scale, pos=0)
axis(2, at = seq(0,max_counts+step*bin2*sum(counts),step*bin2*sum(counts)), label=seq(0,max_freq+step*bin2,step*bin2), cex.axis=cex_scale, las=2, col="gray72", pos=0)
mtext(expression(epsilon), side = 1, line = 10, cex = 2.5)
mtext("Frequency", side = 2, line = 8.5, cex = 2.5)
abline(v=mean(vec), col=colrom3, lwd=2)
abline(v=quantile(vec,0.025), col=colrom3, lty=2, lwd=2)
abline(v=quantile(vec,0.975), col=colrom3, lty=2, lwd=2)
mtext("c", side = 3, line = 3, adj = 0.05, cex = cex_scale-1, font = 2)
vec <- accp_list$ur
bin1 <- 4
bin2 <- 0.001
step <- 300
y0 <- parms_fit_ranges[4,1]
y1 <- parms_fit_ranges[4,2]
bins <- seq(y0,y1+step,step)
hist(vec, breaks = bins, col = colnew3, border = "white", xlab="", ylab="", xaxt="n", yaxt="n", bty='n', main = NULL)
hist_data <- hist(vec, breaks = bins, plot = FALSE)
counts <- hist_data$counts
max_counts <- max(counts)
norm_counts <- counts / sum(counts)
max_freq <- max(norm_counts)
bins <- seq(y0,y1+step*bin1,step*bin1)
axis(1, at = bins, col="gray72", las=2, cex.axis=cex_scale, pos=0)
axis(2, at = seq(0,max_counts,40), label=round(seq(0,max_counts,40)/sum(counts),2), cex.axis=cex_scale, las=2, col="gray72", pos=0)
mtext("Frequency", side = 2, line = 8.5, cex = 2.5)
mtext(expression(U[0]), side = 1, line = 10, cex = 2.5)
abline(v=mean(vec), col=colrom3, lwd=2)
abline(v=quantile(vec,0.025), col=colrom3, lty=2, lwd=2)
abline(v=quantile(vec,0.975), col=colrom3, lty=2, lwd=2)
mtext("d", side = 3, line = 3, adj = 0.05, cex = cex_scale-1, font = 2)
vec <- accp_list$imp
step <- 0.003
bin1 <- 4
bin2 <- 15
y0 <- parms_fit_ranges[5,1]
y0 <- 0
y1 <- parms_fit_ranges[5,2]
bins <- seq(y0,y1+step,step)
hist(vec, breaks = bins, col = colnew3, border = "white", xlab="", ylab="", xaxt="n", yaxt="n", bty='n', main = NULL)
hist_data <- hist(vec, breaks = bins, plot = FALSE)
counts <- hist_data$counts
max_counts <- max(counts)
norm_counts <- counts / sum(counts)
max_freq <- max(norm_counts)
bins <- seq(y0,y1+step*bin1,step*bin1)
bins <- seq(0,y1,0.01)
axis(1, at = bins, col="gray72", las=2, cex.axis=cex_scale, pos=0)
axis(2, at = seq(0,max_counts+step*bin2*sum(counts),50), label=round(seq(0,max_counts+step*bin2*sum(counts),50)/sum(counts),2), cex.axis=cex_scale, las=2, col="gray72", pos=0)
mtext("Frequency", side = 2, line = 8.5, cex = 2.5)
mtext(expression(M[I] / mu * N[0] ~ (year^{-2})), side = 1, line = 10, cex = 2.5)
abline(v=mean(vec), col=colrom3, lwd=2)
abline(v=quantile(vec,0.025), col=colrom3, lty=2, lwd=2)
abline(v=quantile(vec,0.975), col=colrom3, lty=2, lwd=2)
mtext("e", side = 3, line = 3, adj = 0.05, cex = cex_scale-1, font = 2)
vec <- accp_list$qadj1
bin1 <- 5
bin2 <- 4
step <- 0.05
y0 <- parms_fit_ranges[6,1]
y1 <- parms_fit_ranges[6,2]
y0 <- 0
y1 <- 1
bins <- seq(y0,y1+step,step)
hist(vec, breaks = bins, col = colnew3, border = "white", xlab="", ylab="", xaxt="n", yaxt="n", bty='n', main = NULL)
hist_data <- hist(vec, breaks = bins, plot = FALSE)
counts <- hist_data$counts
max_counts <- max(counts)
norm_counts <- counts / sum(counts)
max_freq <- max(norm_counts)
bins <- seq(y0,y1+step*bin1,step*bin1)
axis(1, at = bins, col="gray72", las=2, cex.axis=cex_scale, pos=0)
axis(2, at = seq(0,max_counts+step*bin2*sum(counts),100), label=round(seq(0,max_counts+step*bin2*sum(counts),100)/sum(counts),2), cex.axis=cex_scale, las=2, col="gray72", pos=0)
mtext("Frequency", side = 2, line = 8.5, cex = 2.5)
mtext(expression(Q[1]), side = 1, line = 10, cex = 2.5)
abline(v=mean(vec), col=colrom3, lwd=2)
abline(v=quantile(vec,0.025), col=colrom3, lty=2, lwd=2)
abline(v=quantile(vec,0.975), col=colrom3, lty=2, lwd=2)
mtext("f", side = 3, line = 3, adj = 0.05, cex = cex_scale-1, font = 2)
vec <- accp_list$qadj2
bin1 <- 5
bin2 <- 2
step <- 0.05
y0 <- parms_fit_ranges[7,1]
y1 <- parms_fit_ranges[7,2]
y0 <- 0
y1 <- 1
bins <- seq(y0,y1+step,step)
hist(vec, breaks = bins, col = colnew3, border = "white", xlab="", ylab="", xaxt="n", yaxt="n", bty='n', main = NULL)
hist_data <- hist(vec, breaks = bins, plot = FALSE)
counts <- hist_data$counts
max_counts <- max(counts)
norm_counts <- counts / sum(counts)
max_freq <- max(norm_counts)
bins <- seq(y0,y1+step*bin1,step*bin1)
axis(1, at = bins, col="gray72", las=2, cex.axis=cex_scale, pos=0)
axis(2, at = seq(0,max_counts+step*bin2*sum(counts),step*bin2*sum(counts)), label=seq(0,max_freq+step*bin2,step*bin2), cex.axis=cex_scale, las=2, col="gray72", pos=0)
mtext("Frequency", side = 2, line = 8.5, cex = 2.5)
mtext(expression(Q[2]), side = 1, line = 10, cex = 2.5)
abline(v=mean(vec), col=colrom3, lwd=2)
abline(v=quantile(vec,0.025), col=colrom3, lty=2, lwd=2)
abline(v=quantile(vec,0.975), col=colrom3, lty=2, lwd=2)
mtext("g", side = 3, line = 3, adj = 0.05, cex = cex_scale-1, font = 2)
vec <- accp_list$qadj3
bin1 <- 5
bin2 <- 2
step <- 0.05
y0 <- parms_fit_ranges[8,1]
y1 <- parms_fit_ranges[8,2]
y0 <- 0
y1 <- 1
bins <- seq(y0,y1+step,step)
hist(vec, breaks = bins, col = colnew3, border = "white", xlab="", ylab="", xaxt="n", yaxt="n", bty='n', main = NULL)
hist_data <- hist(vec, breaks = bins, plot = FALSE)
counts <- hist_data$counts
max_counts <- max(counts)
norm_counts <- counts / sum(counts)
max_freq <- max(norm_counts)
bins <- seq(y0,y1+step*bin1,step*bin1)
axis(1, at = bins, col="gray72", las=2, cex.axis=cex_scale, pos=0)
axis(2, at = seq(0,max_counts+step*bin2*sum(counts),step*bin2*sum(counts)), label=seq(0,max_freq+step*bin2,step*bin2), cex.axis=cex_scale, las=2, col="gray72", pos=0)
mtext("Frequency", side = 2, line = 8.5, cex = 2.5)
mtext(expression(Q[3]), side = 1, line = 10, cex = 2.5)
abline(v=mean(vec), col=colrom3, lwd=2)
abline(v=quantile(vec,0.025), col=colrom3, lty=2, lwd=2)
abline(v=quantile(vec,0.975), col=colrom3, lty=2, lwd=2)
mtext("h", side = 3, line = 3, adj = 0.05, cex = cex_scale-1, font = 2)
vec <- accp_list$qadj4
bin1 <- 5
bin2 <- 2
step <- 0.05
y0 <- parms_fit_ranges[9,1]
y1 <- parms_fit_ranges[9,2]
y0 <- 0
y1 <- 1
bins <- seq(y0,y1+step,step)
hist(vec, breaks = bins, col = colnew3, border = "white", xlab="", ylab="", xaxt="n", yaxt="n", bty='n', main = NULL)
hist_data <- hist(vec, breaks = bins, plot = FALSE)
counts <- hist_data$counts
max_counts <- max(counts)
norm_counts <- counts / sum(counts)
max_freq <- max(norm_counts)
bins <- seq(y0,y1+step*bin1,step*bin1)
axis(1, at = bins, col="gray72", las=2, cex.axis=cex_scale, pos=0)
axis(2, at = seq(0,max_counts+step*bin2*sum(counts),100), label=round(seq(0,max_counts+step*bin2*sum(counts),100)/sum(counts),2), cex.axis=cex_scale, las=2, col="gray72", pos=0)
mtext("Frequency", side = 2, line = 8.5, cex = 2.5)
mtext(expression(Q[4]), side = 1, line = 10, cex = 2.5)
abline(v=mean(vec), col=colrom3, lwd=2)
abline(v=quantile(vec,0.025), col=colrom3, lty=2, lwd=2)
abline(v=quantile(vec,0.975), col=colrom3, lty=2, lwd=2)
mtext("i", side = 3, line = 3, adj = 0.05, cex = cex_scale-1, font = 2)
dev.off()

# Calculate and plot Spearman correlations between accepted parameters
ycor <- which(years_sim == c(y_end_plot))
datacor <- as.data.frame(counter_incidence[,(ycor-1)*dt+1])
for(i in 1:ncol(accp_list)){
  datacor <- cbind(datacor,accp_list[,i])
}
names(datacor)[1] <-"inc"
for(i in 1:ncol(accp_list)){
  names(datacor)[i+1] <- names(accp_list)[i]
}
cor_matrix <- matrix(NA, nrow = ncol(datacor) - 1, ncol = ncol(datacor) - 1)
cor_names <- matrix("", nrow = ncol(datacor) - 1, ncol = ncol(datacor) - 1)
p_value_matrix <- matrix(NA, nrow = ncol(datacor) - 1, ncol = ncol(datacor) - 1)
for (i in 2:ncol(datacor)) {
  for (j in 2:ncol(datacor)) {
    cor_test <- cor.test(datacor[, i], datacor[, j], method = "spearman")
    cor_name <- paste(names(datacor)[i], names(datacor)[j], sep = "_")
    cor_matrix[i - 1, j - 1] <- cor_test$estimate
    p_value_matrix[i - 1, j - 1] <- cor_test$p.value
    cor_names[i - 1, j - 1] <- cor_name
  }
}
cor_df <- expand.grid(x = 1:nrow(cor_matrix), y = 1:ncol(cor_matrix))
cor_df$z <- as.vector(cor_matrix)
lb <- c(expression(lambda),expression(omega),expression(epsilon), expression(U[0]),
        expression(M[I] / beta * N[0] ), expression(Q[1]), expression(Q[2]), expression(Q[3]), expression(Q[4]))
heatmap <- ggplot(data = cor_df, aes(x = x, y = y, fill = z)) +
  geom_tile() +
  scale_fill_gradientn(colours = colbar2, limits = c(-1, 1), name = "") +
  labs(title = "", x = "", y = "") +
  theme_minimal() +
  theme(
    axis.line = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.key.width = unit(1, "cm"),  # Adjust the width of legend keys
    legend.key.height = unit(2, "cm"),  # Adjust the height of legend keys
    legend.text = element_text(size = 12),  # Adjust the size of legend text
    axis.text = element_text(size = 12)  # Increase the size of tick labels
  ) +
  scale_x_continuous(breaks = 1:length(lb), labels = lb) +
  scale_y_continuous(breaks = 1:length(lb), labels = lb) +
  annotate("text", x = max(cor_df$x) + 1, y = nrow(cor_matrix)/2, label = "Spearman correlation", angle = -270, hjust = 0.5, vjust = 0.5, size = 6)
ggsave("PLOTS/precure/FigS6.pdf", plot = heatmap, width = 1.25*size_plot, height = size_plot)

# Plot fitted curves and HIV incidence
rescmtext <- 0.85
rescleg <- 0.95
morew <- 2
moreh <- 2.25
textxspace <- 0.12
titxspace <- 0.13
rtrt <- moreh/2
parmar <-  c(moreh*3,morew*6.3,moreh*1.1,morew*0.3)
parmgp <- c(3.6,1,0)*morew
paroma <- c(1.5*moreh,0,0,0)
resclinew <- 0.75
cexlab <- moreh
rescmtext <- 0.8*morew/2 
rescyear <- 0.69*morew/2
rescpo <- 1.5
pdf(file = ("PLOTS/precure/Fig1.pdf"), width = morew*size_plot, height = moreh*size_plot * rescale)
layout(matrix(c(1,1, 2,2, 3, 3,0,0), 2, 4, byrow = TRUE))
par(mar = parmar)
par(mgp = parmgp)
par(oma = paroma)  # Bottom, left, top, right
y0 <- 0
y1 <- max(ci_newdiag[1,(first_year_fit_idx-1)*dt])
y1 <- 400
plot(0, xlim = c(years_calibration[first_year_fit_idx], y_end_plot), ylim = c(y0, y1 + (y1 - y0) * 0.), col = "white", xlab = "", ylab = expression(atop("New diagnoses", "(1/100,000 person-years)")), 
     xaxt = "n", yaxt = "n", bty = "n", cex.lab=cexlab)
axis(1, at = seq(years_calibration[first_year_fit_idx]-1, y_end_plot+1, 1), labels = NA, col = "gray72", pos=0)
step <- 100
axis(2, at = round(seq(0, y1, step),0), cex.axis = cexlab, col = "gray72", las=2)
text(x = seq(years_calibration[first_year_fit_idx], y_end_plot, 1), y = par("usr")[3] - textxspace * diff(par("usr")[3:4]),
     labels = seq(years_calibration[first_year_fit_idx], y_end_plot, 1),
     srt = 45, adj = c(0.5, 0.5), xpd = TRUE, cex = cexlab)
abline(h=seq(y0,y1,step), lty=3, col="gray72")
lines(seq(years_sim[1],years_sim[length(years_sim)],by=1/dt), mean_newdiag / perinc, col = adjustcolor(colcurr, alpha = 1), lwd = moreh*morew*resclinew)
polygon(c(seq(years_sim[1],years_sim[length(years_sim)],by=1/dt), rev(seq(years_sim[1],years_sim[length(years_sim)],by=1/dt))), c(ci_newdiag[1,] / perinc, rev(ci_newdiag[2,] / perinc)), col = adjustcolor(colcurr, alpha = 0.1), border = NA)
points(c(years_calibration,2022), c(newinc,newinc_2022) / perinc, col = colrom1, pch = 19, lwd=.5, cex = (moreh+morew)/rescpo)
legend("topright", legend = c("Data", "Model"), col = c(colrom1, colcurr), pch = c(19, NA), lty = c(NA, 1), lwd=c(.75,(moreh+morew)/rescpo), cex = cexlab, bty = "n")
mtext("a", side = 3, line = 0.45, adj = 0.0, cex = cexlab*rescmtext, font = 2)
step <- 500
y0 <- 0
yy0 <- floor(y0 / step) * step
y1 <- max(c((max(ci_undiag[1,]))))
y1 <- ci_undiag[1,(idx_end_plot-1)*dt]
y1 <- 2000
plot(0, xlim = c(y_start_plot, y_end_plot), ylim = c(yy0, y1), col = "white", xlab = "", ylab = expression(atop("Number of", "undiagnosed infections")), 
     xaxt = "n", yaxt = "n", bty = "n", cex.lab=cexlab)
axis(1, at = seq(y_start_plot-1, y_end_plot+1, 1), labels = NA, col = "gray72", pos=0)
axis(2, at = round(seq(yy0, y1, step),0), cex.axis=cexlab, col = "gray72", las=2)
text(x = seq(y_start_plot, y_end_plot, 1), y = par("usr")[3] - textxspace * diff(par("usr")[3:4]),
     labels = seq(y_start_plot, y_end_plot, 1),
     srt = 45, adj = c(0.5, 0.5), xpd = TRUE, cex = cexlab)
abline(h=seq(yy0,y1,step), lty=3, col=adjustcolor("gray72", alpha.f = 0.75))
lines(seq(years_sim[1],years_sim[length(years_sim)],by=1/dt), mean_undiag, col = adjustcolor(colcurr, alpha = 1), lwd = moreh*morew*resclinew)
polygon(c(seq(years_sim[1],years_sim[length(years_sim)],by=1/dt), rev(seq(years_sim[1],years_sim[length(years_sim)],by=1/dt))), c(ci_undiag[1,], rev(ci_undiag[2,])), col = adjustcolor(colcurr, alpha = 0.1), border = NA)
points(years_calibration, currently_undiag, col = colrom1, pch = 19, lwd=.5, cex = (moreh+morew)/rescpo)
points(2022, undiag_2022, col = colrom1, pch = 19, lwd=.5, cex = (moreh+morew)/rescpo)
segments(x0=years_calibration, x1=years_calibration, y0=cilow_undiag, y1=cihigh_undiag, col = adjustcolor(colrom1,alpha.f = 1), lwd=(moreh+morew)/rescpo)
segments(x0=2022, x1=2022, y0=lo_undiag_2022, y1=hi_undiag_2022, col = adjustcolor(colrom1,alpha.f = 1), lwd=(moreh+morew)/rescpo)
legend("topright", legend = c("Data", "Model"), col = c(colrom1, colcurr), pch = c(19, NA), lty = c(NA, 1), lwd=c(.75,(moreh+morew)/2), cex = cexlab, bty = "n")
mtext("b", side = 3, line = 0.45, adj = 0.0, cex = cexlab*rescmtext, font = 2)
y0 <- 0
y1 <- max(c((max(ci_incidence[1,])),max(newly_infected_est)))
y1 <- 150
step <- 30
plot(0, xlim = c(y_start_plot, y_end_plot), ylim = c(y0, y1), col = "white", 
     xlab = "", 
     ylab = expression(atop("New infections", "(1/100,000 person-years)")), 
     xaxt = "n", yaxt = "n", bty = "n", cex.lab=cexlab)
axis(1, at = seq(y_start_plot-1, y_end_plot+1, 1), labels = NA, col = "gray72", pos=0)
yy0 <- floor(y0 / 100) * 100
axis(2, at = (seq(yy0, y1, step)), cex.axis=cexlab, col = "gray72", las=2)
text(x = seq(y_start_plot, y_end_plot, 1), y = par("usr")[3] - textxspace * diff(par("usr")[3:4]), 
     labels = seq(y_start_plot, y_end_plot, 1), 
     srt = 45, adj = c(0.5, 0.5), xpd = TRUE, cex = cexlab)
abline(h=seq(yy0,y1+step,step), lty=3, col="gray72")
lines(seq(years_sim[1],years_sim[length(years_sim)],by=1/dt), mean_incidence / perinc, col = adjustcolor(colcurr, alpha = 1), lwd = moreh*morew*resclinew)
polygon(c(seq(years_sim[1],years_sim[length(years_sim)],by=1/dt), rev(seq(years_sim[1],years_sim[length(years_sim)],by=1/dt))), c(ci_incidence[1,] / perinc, rev(ci_incidence[2,] / perinc)), col = adjustcolor(colcurr, alpha = 0.1), border = NA)
mtext("c", side = 3, line = 0.45, adj = 0.0, cex = cexlab*rescmtext, font = 2)
mtext("Year", side=1, outer=TRUE, line=1., at=0.32, cex=cexlab*rescyear, adj = c(0.5, 0.5), font=1)  # Below the D panel
legend("topright", legend = "Model", col = colcurr, lty = 1, lwd=(moreh+morew)/rescpo, cex = cexlab, bty = "n")
dev.off()

# Plot validation curves
pdf(file = ("PLOTS/precure/FigS1.pdf"), width = morew*size_plot, height = moreh*size_plot * rescale)
layout(matrix(c(1,1, 2,2, 0, 3,3,0), 2, 4, byrow = TRUE))
par(mar = parmar)
par(mgp = parmgp)
par(oma = paroma)  # Bottom, left, top, right
y0 <- 0
y1 <- max(ci_prep[1,])
y1 <- 10000
step <- 2500
plot(0, xlim = c(y_start_plot, y_end_plot ), ylim = c(y0, y1), col = "white", xlab = "", ylab = "PrEP users", xaxt = "n", yaxt = "n", bty='n', cex.lab=cexlab)
axis(1, at = seq(y_start_plot-1, y_end_plot+1, 1), labels = NA, col="gray72", pos=0)
axis(2, at = round(seq(0, y1, step),0), cex.axis = cexlab, las=2, col = "gray72")
text(x = seq(y_start_plot, y_end_plot, 1), y = par("usr")[3] - textxspace * diff(par("usr")[3:4]), 
     labels = seq(y_start_plot, y_end_plot, 1), 
     srt = 45, adj = c(0.5, 0.5), xpd = TRUE, cex = cexlab)
abline(h=seq(y0,y1,step), lty=3, col="gray72")
lines(seq(years_sim[1],years_sim[length(years_sim)],by=1/dt), mean_prep, col = adjustcolor(colcurr, alpha = 1), lwd = moreh*morew*resclinew)
polygon(c(seq(years_sim[1],years_sim[length(years_sim)],by=1/dt), rev(seq(years_sim[1],years_sim[length(years_sim)],by=1/dt))), c(ci_prep[1,], rev(ci_prep[2,])), col = adjustcolor(colrom3, alpha = 0.1), border = NA)
points(2019:2021, prep_y[prep_fit_idxs], col = colrom1, pch = 19, lwd=.5, cex = (moreh+morew)/rescpo)
mtext("a", side = 3, line = 0.5, adj = 0.0, cex = cexlab*rescmtext, font = 2)
legend(x = "bottomright", legend = c("Data", "Model"), col = c(colrom1, colcurr), pch = c(19, NA), lty = c(NA, 1), lwd=c(.75,moreh), cex = cexlab*rescleg, bty = "n")
y0 <- 0
y1 <- 400
step <- 100
yy0 <- floor(y0 / step) * step
plot(0, xlim = c(y_start_plot, y_end_plot), ylim = c(y0, y1 + (y1 - y0) * 0.), col = "white", xlab = "", ylab = expression(atop("New individuals on ART", "immigrating to the NL (1/year)")), xaxt = "n", yaxt = "n", bty = "n", cex.lab=cexlab)
axis(1, at = seq(y_start_plot-1, y_end_plot+1, 1), labels = NA, col = "gray72", pos=0)
axis(2, at = round(seq(yy0, y1, step),0), cex.axis = cexlab, col = "gray72", las=2)
text(x = seq(y_start_plot, y_end_plot, 1), y = par("usr")[3] - textxspace * diff(par("usr")[3:4]), 
     labels = seq(y_start_plot, y_end_plot, 1), 
     srt = 45, adj = c(0.5, 0.5), xpd = TRUE, cex = cexlab)
abline(h=seq(yy0,y1,step), lty=3, col="gray72")
lines(seq(years_sim[1],years_sim[length(years_sim)],by=1/dt), mean_impdiag, col = adjustcolor(colcurr, alpha = 1), lwd = moreh*morew*resclinew)
polygon(c(seq(years_sim[1],years_sim[length(years_sim)],by=1/dt), rev(seq(years_sim[1],years_sim[length(years_sim)],by=1/dt))), c(ci_impdiag[1,], rev(ci_impdiag[2,])), col = adjustcolor(colcurr, alpha = 0.1), border = NA)
points(c(years_calibration,2022), c(newimp,newimp_2022), col = colrom1, pch = 19, lwd=.5, cex = (moreh+morew)/rescpo)
mtext("b", side = 3, line = 0.5, adj = 0.0, cex = cexlab*rescmtext, font = 2)
legend("bottomright", legend = c("Data", "Model"), col = c(colrom1, colcurr), pch = c(19, NA), lty = c(NA, 1), lwd=c(.75,moreh), cex = cexlab*rescleg, bty = "n")
y0 <- 0
step <- 0.2
yy0 <- floor(y0 / step) * step
y1 <- 1
plot(0, xlim = c(y_start_plot, y_end_plot), ylim = c(yy0, y1 + (y1 - y0) * 0.), col = "white", xlab = "", ylab = "ART coverage (%)", xaxt = "n", yaxt = "n", bty = 'n', cex.lab=cexlab)
axis(1, at = seq(y_start_plot-1, y_end_plot+1, 1), labels = NA, col = "gray72", pos=0)
axis(2, at = seq(y0, y1, step), labels = round(seq(y0, y1, step)*100,0), cex.axis=cexlab, col = "gray72", las=2)
text(x = seq(y_start_plot, y_end_plot, 1), y = par("usr")[3] - textxspace * diff(par("usr")[3:4]), 
     labels = seq(y_start_plot, y_end_plot, 1), 
     srt = 45, adj = c(0.5, 0.5), xpd = TRUE, cex = cexlab)
abline(h=seq(yy0,y1+step,step), lty=3, col="gray72")
lines(seq(years_sim[1],years_sim[length(years_sim)],by=1/dt), mean_art2, col = adjustcolor(colcurr, alpha = 1), lwd = moreh*morew*resclinew)
polygon(c(seq(years_sim[1],years_sim[length(years_sim)],by=1/dt), rev(seq(years_sim[1],years_sim[length(years_sim)],by=1/dt))), c(ci_art2[1,], rev(ci_art2[2,])), col = adjustcolor(colcurr, alpha = 0.2), border = NA)
points(years_calibration, treat_currently / currently_infected_est, col = colrom1, pch = 19, lwd=.5, cex = (moreh+morew)/rescpo)
points(2022, treat_2022 / prev_2022, col = colrom1, pch = 19, lwd=.5, cex = (moreh+morew)/rescpo)
mtext("c", side = 3, line = 0.5, adj = 0.0, cex = cexlab*rescmtext, font = 2)
mtext("Year", side=1, outer=TRUE, line=1, at=0.565, cex=cexlab*rescyear, adj = c(0.5, 0.5))  # Below the D panel
legend("bottomright", legend = c("Data", "Model"), col = c(colrom1, colcurr), pch = c(19, NA), lty = c(NA, 1), lwd=c(.75,moreh), cex = cexlab*rescleg, bty = "n")

dev.off()
