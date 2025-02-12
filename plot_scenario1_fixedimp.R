rm(list = ls())

load("PLOTS/colors.Rdata")

efficacies <- c(0.2,0.9)
# cure_efficacy_option <- "minimum"
cure_efficacy_option <- "optimal"
if(cure_efficacy_option == "optimal"){
  eff <- efficacies[2]
}else{
  eff <- efficacies[1]
}

ar <- 0.9
# ar <- 0.5
# ar <- 0.1

# Infectivities of Ic
hf <- 0 # baseline
# hf <- 1 # sensitivity analysis
if(hf == 0){
  hcfile <- "MAIN"
}else{
  hcfile <- "SA"
}

y_cure <- 2026

load(paste0("RESULTS/SCENARIOS/SC1_eff", round(eff*100,0),"_cu", round(ar*100,0), "_hc", hcfile,"_curestart",y_cure,"_fixedimp.RData"))

# Select monitoring strategy
# g <- 1 # standard of care
g <- 2 # monitoring
# g <- 3 # follow-up

years_sim <- years_calibration[1]:y_end
perinc <- parms["NN"] / 1e5
pop <- parms["NN"]

mean_inc_nocure <- array(apply(incidence_nocure_time, 2, mean), dim=c(nt))
lo_inc_nocure <- array(apply(incidence_nocure_time, 2, quantile, prob=0.025), dim=c(nt))
hi_inc_nocure <- array(apply(incidence_nocure_time, 2, quantile, prob=0.975), dim=c(nt))
mean_prev_nocure <- array(apply(prevalence_nocure_time, 2, mean), dim=c(nt))
lo_prev_nocure <- array(apply(prevalence_nocure_time, 2, quantile, prob=0.025), dim=c(nt))
hi_prev_nocure <- array(apply(prevalence_nocure_time, 2, quantile, prob=0.975), dim=c(nt))

mean_inc_art_nocure <- array(apply(incidence_art_nocure_time, 2, mean), dim=c(nt))
mean_inc_diag_nocure <- array(apply(incidence_diag_nocure_time, 2, mean), dim=c(nt))
mean_inc_undiag_nocure <- array(apply(incidence_undiag_nocure_time, 2, mean), dim=c(nt))

curet <- which(years_sim==y_cure)
cureidx <- (curet-1)*dt + 1

subset_art_given_nocure <- art_given_nocure_time[,cureidx:nt]
cumulative_art_given_nocure <- apply(art_given_nocure_time, 1, sum) / dt
save(cumulative_art_given_nocure, file = ("RESULTS/cumartgivennocure_post_fixedimp.RData"))
subset_art_given <- art_given_time_pool[,,,cureidx:nt]
cumulative_art_given <- array(apply(subset_art_given, c(1, 2, 3), sum), dim=c(n_monitors, n_sa, len_post)) / dt

subset_newly_cured <- newly_cured_time_pool[,,,cureidx:nt]
cumulative_newly_cured <- array(apply(subset_newly_cured, c(1, 2, 3), sum), dim=c(n_monitors, n_sa, len_post)) / dt

cumulative_art_averted <- cumulative_art_given
for (p in 1:len_post) {
  cumulative_art_averted[,,p] <- 100*(-cumulative_art_given[,,p] + cumulative_art_given_nocure[p]) / cumulative_newly_cured[,,p]
}

rescleg <- 0.95
morew <- 2
moreh <- 2.25
rescyear <- 0.69*morew/2

if(hcfile=="MAIN"){
  if(y_cure == 2026){
    pdf(file = ("PLOTS/postcure/Fig2_fixedimp.pdf"), width = morew*size_plot, height = moreh*size_plot * rescale)
  }else{
    pdf(file = (paste0("PLOTS/postcure/FigS11","_fixedimp.pdf")), width = morew*size_plot, height = moreh*size_plot * rescale)
  }
}else{
  pdf(file = ("PLOTS/postcure/FigS9_fixedimp.pdf"), width = morew*size_plot, height = moreh*size_plot * rescale)
}

first_year_show <- 2024

# Figure parameters
textxspace <- 0.12
parmar <-  c(moreh*3,morew*6.3,moreh*1.1,morew*0.3)
parmgp <- c(3.6,1,0)*morew
paroma <- c(1.5*moreh,0,0,0)
cexlab <- moreh
rescmtext <- 0.8*morew/2
resclinew <- 0.75
layout(matrix(c(1,1, 2,2, 3, 3,4,4), 2, 4, byrow = TRUE))
par(mar = parmar)
par(mgp = parmgp)
par(oma = paroma)

# Panel a
y0 <- 0
y1 <- 120 
step <- 30
if(hcfile=="SA"){
  y1 <- y1*2
  step <- step*2
}
plot(0, xlim = c(first_year_show, y_end ), ylim = c(y0, y1), col = "white", xlab = "", ylab = expression(atop("New infections", "(1/100,000 person-years)")), xaxt = "n", yaxt = "n", bty = "n", cex.lab=cexlab)
axis(1, at = seq(first_year_show-1, y_end, 1), labels = NA, col = "gray72",pos=0)
axis(2, at = round(seq(y0, y1, step),0), cex.axis=cexlab, col = "gray72", las=2)
text(x = seq(first_year_show, y_end, 2), y = par("usr")[3] - textxspace * diff(par("usr")[3:4]),
     labels = seq(first_year_show, y_end, 2),
     srt = 45, adj = c(0.5, 0.5), xpd = TRUE, cex = cexlab)
abline(h=seq(y0,y1+step,step), lty=3, col="gray72")
leg_sa <- rep("",n_sa)
for(i in 1:n_sa){
  lines(seq(years_sim[first_year_fit_idx]+1,years_sim[length(years_sim)],by=1/dt), mean_inc_new[g,i,(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)] / perinc, col = colors1[i], lwd = moreh*morew*resclinew)
  polygon(c(seq(years_sim[first_year_fit_idx]+1,years_sim[length(years_sim)],by=1/dt), rev(seq(years_sim[first_year_fit_idx]+1,years_sim[length(years_sim)],by=1/dt))), c(q1_inc_new[g,i,(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)] / perinc, rev(q2_inc_new[g,i,(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)]) / perinc), col = adjustcolor(colors1[i], alpha = 0.1), border = NA)
  leg_sa[i] <- 1 / failures[i]
}
lines(seq(years_sim[first_year_fit_idx]+1,years_sim[length(years_sim)],by=1/dt), mean_inc_nocure[(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)] / perinc, lty=1, col = colcurr, lwd = moreh*morew*resclinew)
polygon(c(seq(years_sim[first_year_fit_idx]+1,years_sim[length(years_sim)],by=1/dt), rev(seq(years_sim[first_year_fit_idx]+1,years_sim[length(years_sim)],by=1/dt))), c(lo_inc_nocure[(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)] / perinc, rev(hi_inc_nocure[(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)]) / perinc), col = adjustcolor(colcurr, alpha = 0.1), border = NA)
mtext("a", side = 3, line = 0.5, adj = 0.0, cex = cexlab*rescmtext, font = 2)
arrows(y_cure, par("usr")[3] + 0.25 * diff(par("usr")[3:4]), y_cure, 0, col="darkred", lwd=3, length=0.1)

# Panel b
y0 <- 0
y1 <- 2000
plot(0, xlim = c(first_year_show, y_end ), ylim = c(y0, y1), col = "white", xlab = "", ylab = expression(atop("New rebounds", "(1/100,000 person-years)")),
     xaxt = "n", yaxt = "n", bty = "n", cex.lab=cexlab)
axis(1, at = seq(first_year_show-1, y_end, 1), labels = NA, col = "gray72",pos=0)
step <- 500
axis(2, at = round(seq(y0, y1, step),0), cex.axis=cexlab, col = "gray72", las=2)
text(x = seq(first_year_show, y_end, 2), y = par("usr")[3] - textxspace * diff(par("usr")[3:4]),
     labels = seq(first_year_show, y_end, 2),
     srt = 45, adj = c(0.5, 0.5), xpd = TRUE, cex = cexlab)
abline(h=seq(y0,y1+step,step), lty=3, col="gray72")
leg_sa <- rep("",n_sa)
for(i in 1:n_sa){
  lines(seq(years_sim[first_year_fit_idx]+1,years_sim[length(years_sim)],by=1/dt), mean_inc_reb[g,i,(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)] / perinc, col = colors1[i], lwd = moreh*morew*resclinew)
  polygon(c(seq(years_sim[first_year_fit_idx]+1,years_sim[length(years_sim)],by=1/dt), rev(seq(years_sim[first_year_fit_idx]+1,years_sim[length(years_sim)],by=1/dt))), c(q1_inc_reb[g,i,(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)] / perinc, rev(q2_inc_reb[g,i,(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)]) / perinc), col = adjustcolor(colors1[i], alpha = 0.1), border = NA)
  leg_sa[i] <- 1 / failures[i]
}
lines(years_sim[first_year_fit_idx:length(years_sim)], rep(0,length(first_year_fit_idx:length(years_sim))), lty=1, col = colcurr, lwd = moreh*morew*resclinew)
mtext("b", side = 3, line = 0.5, adj = 0.0, cex = cexlab*rescmtext, font = 2)
arrows(y_cure, par("usr")[3] + 0.25 * diff(par("usr")[3:4]), y_cure, 0, col="darkred", lwd=3, length=0.1)

# Panel c
step <- 5000
step <- 0.02
y0 <- 0
y1 <- 0.08
plot(0, xlim = c(first_year_show, y_end ), ylim = c(y0, y1), col = "white", xlab = "", ylab = "Prevalence (%)", xaxt = "n", yaxt = "n", bty = "n", cex.lab=cexlab)
axis(1, at = seq(first_year_show-1, y_end, 1), labels = NA, col = "gray72",pos=0)
axis(2, at = round(seq(y0, y1, step),2), labels = 100*round(seq(y0, y1, step),2), cex.axis=cexlab, col = "gray72", las=2)
text(x = seq(first_year_show, y_end, 2), y = par("usr")[3] - textxspace * diff(par("usr")[3:4]),
     labels = seq(first_year_show, y_end, 2),
     srt = 45, adj = c(0.5, 0.5), xpd = TRUE, cex = cexlab)
abline(h=seq(y0,y1+step,step), lty=3, col="gray72")
leg_sa <- rep("",n_sa)
for(i in 1:n_sa){
  lines(seq(years_sim[first_year_fit_idx]+1,years_sim[length(years_sim)],by=1/dt), mean_prev_new[g,i,(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)] / pop + mean_prev_reb[g,i,(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)] / pop, col = colors1[i], lwd = moreh*morew*resclinew)
  polygon(c(seq(years_sim[first_year_fit_idx]+1,years_sim[length(years_sim)],by=1/dt), rev(seq(years_sim[first_year_fit_idx]+1,years_sim[length(years_sim)],by=1/dt))), c((q1_prev_new[g,i,(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)] + q1_prev_reb[g,i,(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)]), rev(q2_prev_new[g,i,(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)]+q2_prev_reb[g,i,(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)]))/pop, col = adjustcolor(colors1[i], alpha = 0.1), border = NA)
  leg_sa[i] <- 1 / failures[i]
}
lines(seq(years_sim[first_year_fit_idx]+1,years_sim[length(years_sim)],by=1/dt), mean_prev_nocure[(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)] / pop, lty=1, col = colcurr, lwd = moreh*morew*resclinew)
polygon(c(seq(years_sim[first_year_fit_idx]+1,years_sim[length(years_sim)],by=1/dt), rev(seq(years_sim[first_year_fit_idx]+1,years_sim[length(years_sim)],by=1/dt))), c(lo_prev_nocure[(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)] / pop, rev(hi_prev_nocure[(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)]) / pop), col = adjustcolor(colcurr, alpha = 0.1), border = NA)
mtext("c", side = 3, line = 0.5, adj = 0.0, cex = cexlab*rescmtext, font = 2)
arrows(y_cure, par("usr")[3] + 0.25 * diff(par("usr")[3:4]), y_cure, 0, col="darkred", lwd=3, length=0.1)

# Panel d
y0 <- 0
y1 <- 100
plot(0, xlim = c(first_year_show, y_end ), ylim = c(y0, y1), col = "white", xlab = "", ylab = "Cure coverage (%)", xaxt = "n", yaxt = "n", bty = "n", cex.lab=cexlab)
axis(1, at = seq(first_year_show-1, y_end, 1), labels = NA, col = "gray72",pos=0)
step <- 25
axis(2, at = seq(y0, y1, step), cex.axis=cexlab, col = "gray72", las=2)
text(x = seq(first_year_show, y_end, 2), y = par("usr")[3] - textxspace * diff(par("usr")[3:4]),
     labels = seq(first_year_show, y_end, 2),
     srt = 45, adj = c(0.5, 0.5), xpd = TRUE, cex = cexlab)
abline(h=seq(y0,y1+step,step), lty=3, col="gray72")
leg_sa <- rep("",n_sa)
for(i in 1:n_sa){
  lines(seq(years_sim[first_year_fit_idx]+1,years_sim[length(years_sim)],by=1/dt),round(100*mean_coverage[g,i,(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)], 1), col = colors1[i], lwd = moreh*morew*resclinew)
  polygon(c(seq(years_sim[first_year_fit_idx]+1,years_sim[length(years_sim)],by=1/dt), rev(seq(years_sim[first_year_fit_idx]+1,years_sim[length(years_sim)],by=1/dt))), c(q1_coverage[g,i,(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)], rev(q2_coverage[g,i,(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)])), col = adjustcolor(colors1[i], alpha = 0.1), border = NA)
  if(i == 1){
    leg_sa[i] <- "never"
  }else{
    leg_sa[i] <- paste0(round(as.integer(1 / failures[i]),0)," years")
  }
}
lines(years_sim[first_year_fit_idx:length(years_sim)], rep(0,length(first_year_fit_idx:length(years_sim))), lty=1, col = colcurr, lwd = moreh*morew*resclinew)
legend("bottomright", c(leg_sa,"no cure"), lty = c(rep(1,n_sa),1), lwd = rep(moreh,n_sa+1), col = c(colors1,colcurr), bty = "n", cex = cexlab*rescleg, title = "Time to rebound")
mtext("Year", side=1, outer=TRUE, line=1., at=0.32, cex=cexlab*rescyear, adj = c(0.5, 0.5))  # Below the C panel
mtext("Year", side=1, outer=TRUE, line=1, at=0.82, cex=cexlab*rescyear, adj = c(0.5, 0.5))  # Below the D panel
mtext("d", side = 3, line = 0.5, adj = 0.0, cex = cexlab*rescmtext, font = 2)
arrows(y_cure, par("usr")[3] + 0.25 * diff(par("usr")[3:4]), y_cure, 0, col="darkred", lwd=3, length=0.1)
dev.off()

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
pdf(file = ("PLOTS/postcure/infectivity_nocure_fixedimp.pdf"), width = morew*size_plot, height = moreh*size_plot * rescale)

p1_mat <- matrix(NA,nrow=3,ncol=nt)
p2_mat <- matrix(NA,nrow=3,ncol=nt)
p3_mat <- matrix(NA,nrow=3,ncol=nt)
p_mat <- matrix(NA,nrow=3,ncol=nt)
pl <- 3

for(i in 1:3){
  p1_mat[i,] <- mean_inc_art_new[g,i,] / perinc
  p2_mat[i,] <- mean_inc_undiag_new[g,i,] / perinc
  p3_mat[i,] <- mean_inc_diag_new[g,i,] / perinc
  p_mat[i,] <- mean_inc_new[g,i,] / perinc
}

layout(matrix(c(1,1, 2,2, 3, 3,4,4), 2, 4, byrow = TRUE))
par(mar = parmar)
par(mgp = parmgp)
par(oma = paroma)  # Bottom, left, top, right

y_start_plot <- 2017
y0 <- 0
y1 <- 1
y1 <- 120
plot(0, xlim = c(y_start_plot, 2024 ), ylim = c(y0, y1), col = "white", xlab = "", ylab = "", 
     xaxt = "n", yaxt = "n", bty = "n", cex.lab=cexlab)
plot(0, xlim = c(y_start_plot, 2024 ), ylim = c(y0, y1), col = "white", xlab = "", ylab = "", 
     xaxt = "n", yaxt = "n", bty = "n", cex.lab=cexlab)
plot(0, xlim = c(y_start_plot, 2024 ), ylim = c(y0, y1), col = "white", xlab = "", ylab = "", 
     xaxt = "n", yaxt = "n", bty = "n", cex.lab=cexlab)
plot(0, xlim = c(y_start_plot, 2024 ), ylim = c(y0, y1), col = "white", xlab = "", ylab = expression(atop("New infections by status", "(100,000 person-years)")), 
     xaxt = "n", yaxt = "n", bty = "n", cex.lab=cexlab)
axis(1, at = seq(y_start_plot-1, 2024, 1), labels = NA, col = "gray72",pos=0)
# step <- 0.2
# axis(2, at = seq(y0, y1, step), labels = 100*seq(y0, y1, step), cex.axis=cexlab, col = "gray72", las=2)
step <- 30
axis(2, at = seq(y0, y1, step), labels = seq(y0, y1, step), cex.axis=cexlab, col = "gray72", las=2)
text(x = seq(y_start_plot, 2024, 1), y = par("usr")[3] - textxspace * diff(par("usr")[3:4]),
     labels = seq(y_start_plot, 2024, 1),
     srt = 45, adj = c(0.5, 0.5), xpd = TRUE, cex = cexlab)
abline(h=seq(y0,y1+step,step), lty=3, col="gray72")
polygon(c(seq(years_sim[1],years_sim[length(years_sim)],by=1/dt), rev(seq(years_sim[1],years_sim[length(years_sim)],by=1/dt))), c(p1_mat[pl,], rev(rep(0,length(p_mat[pl,])))), col = colors_inf[1], border = NA)
polygon(c(seq(years_sim[1],years_sim[length(years_sim)],by=1/dt), rev(seq(years_sim[1],years_sim[length(years_sim)],by=1/dt))), c((p1_mat[pl,]+p2_mat[pl,]), rev(p1_mat[pl,])), col = colors_inf[2], border = NA)
polygon(c(seq(years_sim[1],years_sim[length(years_sim)],by=1/dt), rev(seq(years_sim[1],years_sim[length(years_sim)],by=1/dt))), c((p1_mat[pl,]+p2_mat[pl,]+p3_mat[pl,]), rev((p1_mat[pl,]+p2_mat[pl,]))), col = colors_inf[3], border = NA)
mtext("d", side = 3, line = 0.5, adj = 0.0, cex = cexlab*rescmtext, font = 2)
arrows(y_cure, par("usr")[3] + 0.25 * diff(par("usr")[3:4]), y_cure, 0, col="darkred", lwd=3, length=0.1)
legend("topright", legend = c("Diagnosed", "Undiagnosed", "Treated"), col = "white", pt.bg = rev(colors_inf), pch = rep(22,3), pt.cex = 4, cex = cexlab, bty = 'n')
mtext("Year", side=1, outer=TRUE, line=1, at=0.82, cex=cexlab*rescyear, adj = c(0.5, 0.5), font=1)  # Below the D panel
dev.off()
