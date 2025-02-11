rm(list = ls())

load("PLOTS/colors.Rdata")

# Select monitoring and time to rebound of the baseline
gg <- 2
ttr <- 2

eff <- 0.9
ar <- 0.9
load(paste0("RESULTS/SCENARIOS/lambda_increased_SC1.RData"))
years_sim <- years_calibration[1]:y_end
perinc <- parms["NN"] / 1e5

mean_inc_nocure2 <- array(apply(incidence_nocure_time, 2, mean), dim=c(nt))
lo_inc_nocure2 <- array(apply(incidence_nocure_time, 2, quantile, prob=0.025), dim=c(nt))
hi_inc_nocure2 <- array(apply(incidence_nocure_time, 2, quantile, prob=0.975), dim=c(nt))
mean_inc_new_g2_ttr6 <- mean_inc_new[gg,ttr,]
lo_inc_new_g2_ttr6 <- q1_inc_new[gg,ttr,]
hi_inc_new_g2_ttr6 <- q2_inc_new[gg,ttr,]
mean_inc_reb_g2_ttr6 <- mean_inc_reb[gg,ttr,]
lo_inc_reb_g2_ttr6 <- q1_inc_reb[gg,ttr,]
hi_inc_reb_g2_ttr6 <- q2_inc_reb[gg,ttr,]

mean_inc_new_g2_ttr2 <- mean_inc_new[gg,1,]
lo_inc_new_g2_ttr2 <- q1_inc_new[gg,1,]
hi_inc_new_g2_ttr2 <- q2_inc_new[gg,1,]
mean_inc_reb_g2_ttr2 <- mean_inc_reb[gg,1,]
lo_inc_reb_g2_ttr2 <- q1_inc_reb[gg,1,]
hi_inc_reb_g2_ttr2 <- q2_inc_reb[gg,1,]

mean_inc_new_g2_ttr0 <- mean_inc_new[gg,3,]
lo_inc_new_g2_ttr0 <- q1_inc_new[gg,3,]
hi_inc_new_g2_ttr0 <- q2_inc_new[gg,3,]
mean_inc_reb_g2_ttr0 <- mean_inc_reb[gg,3,]
lo_inc_reb_g2_ttr0 <- q1_inc_reb[gg,3,]
hi_inc_reb_g2_ttr0 <- q2_inc_reb[gg,3,]


# Plot parameters
morew <- 2
moreh <- 2.25/2
cexlab <- moreh
textxspace <- 0.1
titxspace <- 0.15
parmar <-  c(5.,6.5,0.5,0.5)
parmgp <- c(3.5,1,0)

leg_sa <- c("no cure", "never", "6 years", "2 years")
cums <- rep(NA,4)

first_year_show <- 2020

pdf(file = ("PLOTS/postcure/FigLIA.pdf"), width = size_plot, height = size_plot * rescale)
par(mar = parmar)
par(mgp = parmgp)
y0 <- 0
y1 <- 150
plot(0, xlim = c(first_year_show, y_end ), ylim = c(y0, y1), col = "white", xlab = "", ylab = expression(atop("New infections", "(1/100,000 person-years)")), xaxt = "n", yaxt = "n", bty = "n", cex.lab=cexlab)
axis(1, at = seq(first_year_show-1, y_end, 1), labels = NA, col = "gray72",pos=0)
step <- 30
axis(2, at = round(seq(y0, y1, step),0), cex.axis=cexlab, col = "gray72", las=2)
text(x = seq(first_year_show, y_end, 2), y = par("usr")[3] - textxspace * diff(par("usr")[3:4]),
     labels = seq(first_year_show, y_end, 2),
     srt = 45, adj = c(0.5, 0.5), xpd = TRUE, cex = cexlab)
text(x = y_end/2 + first_year_show/2, y = par("usr")[3] - (textxspace + titxspace)* diff(par("usr")[3:4]),
     labels = "Year",
     srt = 0, adj = c(0.5, 0.5), xpd = TRUE, cex = 1.25)
abline(h=seq(y0,y1+step,step), lty=3, col="gray72")
lines(seq(years_sim[first_year_fit_idx]+1,years_sim[length(years_sim)],by=1/dt), mean_inc_new_g2_ttr2[(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)] / perinc, col = colors1[1], lwd = 2.5)
polygon(c(seq(years_sim[first_year_fit_idx]+1,years_sim[length(years_sim)],by=1/dt), rev(seq(years_sim[first_year_fit_idx]+1,years_sim[length(years_sim)],by=1/dt))), c(lo_inc_new_g2_ttr2[(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)] / perinc, rev(hi_inc_new_g2_ttr2[(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)]) / perinc), col = adjustcolor(colors1[1], alpha = 0.1), border = NA)
curet <- which(years_sim==y_cure)
cureidx <- (curet-1)*dt + 1
# cums[2] <- sum(mean_inc_new_g2_ttr6[cureidx:nt]) / dt
lines(seq(years_sim[first_year_fit_idx]+1,years_sim[length(years_sim)],by=1/dt), mean_inc_new_g2_ttr6[(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)] / perinc, col = colors1[2], lwd = 2.5)
polygon(c(seq(years_sim[first_year_fit_idx]+1,years_sim[length(years_sim)],by=1/dt), rev(seq(years_sim[first_year_fit_idx]+1,years_sim[length(years_sim)],by=1/dt))), c(lo_inc_new_g2_ttr6[(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)] / perinc, rev(hi_inc_new_g2_ttr6[(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)]) / perinc), col = adjustcolor(colors1[2], alpha = 0.1), border = NA)
curet <- which(years_sim==y_cure)
cureidx <- (curet-1)*dt + 1
lines(seq(years_sim[first_year_fit_idx]+1,years_sim[length(years_sim)],by=1/dt), mean_inc_new_g2_ttr0[(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)] / perinc, col = colors1[3], lwd = 2.5)
polygon(c(seq(years_sim[first_year_fit_idx]+1,years_sim[length(years_sim)],by=1/dt), rev(seq(years_sim[first_year_fit_idx]+1,years_sim[length(years_sim)],by=1/dt))), c(lo_inc_new_g2_ttr0[(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)] / perinc, rev(hi_inc_new_g2_ttr0[(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)]) / perinc), col = adjustcolor(colors1[3], alpha = 0.1), border = NA)
curet <- which(years_sim==y_cure)
cureidx <- (curet-1)*dt + 1
lines(seq(years_sim[first_year_fit_idx]+1,years_sim[length(years_sim)],by=1/dt), mean_inc_nocure2[(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)] / perinc, lty=1, col = colcurr, lwd = 2.5)
polygon(c(seq(years_sim[first_year_fit_idx]+1,years_sim[length(years_sim)],by=1/dt), rev(seq(years_sim[first_year_fit_idx]+1,years_sim[length(years_sim)],by=1/dt))), c(lo_inc_nocure2[(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)] / perinc, rev(hi_inc_nocure2[(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)]) / perinc), col = adjustcolor(colcurr, alpha = 0.1), border = NA)
cums[1] <- sum(mean_inc_nocure2[cureidx:nt]) / dt
legends <- rep("",length(leg_sa))
for(l in 1:length(leg_sa)){
  legends[l] <- leg_sa[l]
}
arrows(y_cure, par("usr")[3] + 0.25 * diff(par("usr")[3:4]), y_cure, 0, col="darkred", lwd=2, length=0.1)
legend("topleft", legend = legends, col = c(colcurr, colors1), lty = rep(1, length(leg_sa)), lwd=rep(2.5, length(leg_sa)), bty = "n", title = "Time to rebound")
dev.off()

cum_new <- cums
cum_new <- 100 * (cum_new - cums[1]) / cums[1]

pdf(file = ("PLOTS/postcure/FigLIC.pdf"), width = size_plot, height = size_plot * rescale)
par(mar = parmar)
par(mgp = parmgp)
y0 <- 0
y1 <- 2000 
plot(0, xlim = c(first_year_show, y_end ), ylim = c(y0, y1), col = "white", xlab = "", ylab = expression(atop("New rebounds", "(1/100,000 person-years)")), xaxt = "n", yaxt = "n", bty = "n", cex.lab=cexlab)
axis(1, at = seq(first_year_show-1, y_end, 1), labels = NA, col = "gray72",pos=0)
step <- 400
axis(2, at = round(seq(y0, y1, step),0), cex.axis=cexlab, col = "gray72", las=2)
text(x = seq(first_year_show, y_end, 2), y = par("usr")[3] - textxspace * diff(par("usr")[3:4]),
     labels = seq(first_year_show, y_end, 2),
     srt = 45, adj = c(0.5, 0.5), xpd = TRUE, cex = cexlab)
text(x = y_end/2 + first_year_show/2, y = par("usr")[3] - (textxspace + titxspace)* diff(par("usr")[3:4]),
     labels = "Year",
     srt = 0, adj = c(0.5, 0.5), xpd = TRUE, cex = 1.25)
abline(h=seq(y0,y1+step,step), lty=3, col="gray72")
lines(seq(years_sim[first_year_fit_idx]+1,years_sim[length(years_sim)],by=1/dt), mean_inc_reb_g2_ttr2[(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)] / perinc, col = colors1[1], lwd = 2.5)
polygon(c(seq(years_sim[first_year_fit_idx]+1,years_sim[length(years_sim)],by=1/dt), rev(seq(years_sim[first_year_fit_idx]+1,years_sim[length(years_sim)],by=1/dt))), c(lo_inc_reb_g2_ttr2[(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)] / perinc, rev(hi_inc_reb_g2_ttr2[(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)]) / perinc), col = adjustcolor(colors1[1], alpha = 0.1), border = NA)
curet <- which(years_sim==y_cure)
cureidx <- (curet-1)*dt + 1
cums[2] <- sum(mean_inc_reb_g2_ttr6[cureidx:nt]) / dt
lines(seq(years_sim[first_year_fit_idx]+1,years_sim[length(years_sim)],by=1/dt), mean_inc_reb_g2_ttr6[(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)] / perinc, col = colors1[2], lwd = 2.5)
polygon(c(seq(years_sim[first_year_fit_idx]+1,years_sim[length(years_sim)],by=1/dt), rev(seq(years_sim[first_year_fit_idx]+1,years_sim[length(years_sim)],by=1/dt))), c(lo_inc_reb_g2_ttr6[(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)] / perinc, rev(hi_inc_reb_g2_ttr6[(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)]) / perinc), col = adjustcolor(colors1[2], alpha = 0.1), border = NA)
lines(seq(years_sim[first_year_fit_idx]+1,years_sim[length(years_sim)],by=1/dt), mean_inc_reb_g2_ttr0[(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)] / perinc, col = colors1[3], lwd = 2.5)
polygon(c(seq(years_sim[first_year_fit_idx]+1,years_sim[length(years_sim)],by=1/dt), rev(seq(years_sim[first_year_fit_idx]+1,years_sim[length(years_sim)],by=1/dt))), c(lo_inc_reb_g2_ttr0[(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)] / perinc, rev(hi_inc_reb_g2_ttr0[(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)]) / perinc), col = adjustcolor(colors1[3], alpha = 0.1), border = NA)
lines(seq(years_sim[first_year_fit_idx]+1,years_sim[length(years_sim)],by=1/dt), rep(0,length(seq(years_sim[first_year_fit_idx]+1,years_sim[length(years_sim)],by=1/dt))), lty=1, col = colcurr, lwd = 2.5)
cums[1] <- 0
arrows(y_cure, par("usr")[3] + 0.25 * diff(par("usr")[3:4]), y_cure, 0, col="darkred", lwd=2, length=0.1)
dev.off()

cum_reb <- cums
