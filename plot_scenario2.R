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

# Infectivities of Ic
hf <- 0 # baseline
# hf <- 1 # sensitivity analysis
if(hf == 0){
  hcfile <- "MAIN"
}else{
  hcfile <- "SA"
}

y_cure <- 2026
# y_cure <- 2030

load(paste0("RESULTS/SCENARIOS/SC2_eff", round(eff*100,0), "_hc", hcfile,"_curestart",y_cure,".RData"))

# Select monitoring strategy
# g <- 1 # standard of care
g <- 2 # monitoring
# g <- 3 # follow-up

years_sim <- years_calibration[1]:y_end
perinc <- parms["NN"] / 1e5
pop <- parms["NN"]

curet <- which(years_sim==y_cure)
cureidx <- (curet-1)*dt + 1

load("RESULTS/cumartgivennocure_post.RData")

subset_art_given <- art_given_time_pool[,,,cureidx:nt]
cumulative_art_given <- array(apply(subset_art_given, c(1, 2, 3), sum), dim=c(n_monitors, n_sa, len_post)) / dt

subset_newly_cured <- newly_cured_time_pool[,,,cureidx:nt]
cumulative_newly_cured <- array(apply(subset_newly_cured, c(1, 2, 3), sum), dim=c(n_monitors, n_sa, len_post)) / dt

cumulative_art_averted_percure <- cumulative_art_given
cumulative_art_averted <- cumulative_art_given
for (p in 1:len_post) {
  cumulative_art_averted_percure[,2:4,p] <- (-cumulative_art_given[,2:4,p] + cumulative_art_given_nocure[p]) / cumulative_newly_cured[,2:4,p]
  cumulative_art_averted[,2:4,p] <- (-cumulative_art_given[,2:4,p] + cumulative_art_given_nocure[p]) / cumulative_art_given_nocure[p]
}

cumulative_art_averted <- cumulative_art_averted[,2:4,]
cumulative_art_averted_percure <- cumulative_art_averted_percure[,2:4,]

rescleg <- 0.95
morew <- 2
moreh <- 2.25
rescyear <- 0.69*morew/2

if(hcfile=="MAIN"){
  if(y_cure == 2026){
    pdf(file = ("PLOTS/postcure/Fig4.pdf"), width = morew*size_plot, height = moreh*size_plot * rescale)
  }else{
    pdf(file = (paste0("PLOTS/postcure/FigS12",".pdf")), width = morew*size_plot, height = moreh*size_plot * rescale)
  }
}else{
  pdf(file = ("PLOTS/postcure/FigS10.pdf"), width = morew*size_plot, height = moreh*size_plot * rescale)
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
y1 <- 100
plot(0, xlim = c(first_year_show, y_end ), ylim = c(y0, y1), col = "white", xlab = "", ylab = expression(atop("New infections", "(1/100,000 person-years)")), xaxt = "n", yaxt = "n", bty = "n", cex.lab=cexlab)
axis(1, at = seq(first_year_show-1, y_end, 1), labels = NA, col = "gray72", pos=0)
step <- 20
axis(2, at = seq(y0, y1, step), cex.axis=cexlab, col = "gray72", las=2)
text(x = seq(first_year_show, y_end, 2), y = par("usr")[3] - textxspace * diff(par("usr")[3:4]),
     labels = seq(first_year_show, y_end, 2),
     srt = 45, adj = c(0.5, 0.5), xpd = TRUE, cex = cexlab)
abline(h=seq(y0,y1+step,step), lty=3, col="gray72")
leg_sa <- rep("",n_sa)
for(i in 2:n_sa){
  lines(seq(years_sim[first_year_fit_idx]+1,years_sim[length(years_sim)],by=1/dt), mean_inc_new[g,i,(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)] / perinc, col = colors2[i], lwd = moreh*morew*resclinew)
  polygon(c(seq(years_sim[first_year_fit_idx]+1,years_sim[length(years_sim)],by=1/dt), rev(seq(years_sim[first_year_fit_idx]+1,years_sim[length(years_sim)],by=1/dt))), c(q1_inc_new[g,i,(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)], rev(q2_inc_new[g,i,(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)])) / perinc, col = adjustcolor(colors2[i], alpha = 0.1), border = NA)
  leg_sa[i] <- paste0(round((1 - exp(- alphas))[i] * 100,0),"%")
}
i <- 1
lines(seq(years_sim[first_year_fit_idx]+1,years_sim[length(years_sim)],by=1/dt), mean_inc_new[g,i,(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)] / perinc, col = colors2[i], lwd = moreh*morew*resclinew, lty=1)
polygon(c(seq(years_sim[first_year_fit_idx]+1,years_sim[length(years_sim)],by=1/dt), rev(seq(years_sim[first_year_fit_idx]+1,years_sim[length(years_sim)],by=1/dt))), c(q1_inc_new[g,i,(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)], rev(q2_inc_new[g,i,(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)])) / perinc, col = adjustcolor(colors2[i], alpha = 0.1), border = NA)
mtext("a", side = 3, line = 0.5, adj = 0.0, cex = cexlab*rescmtext, font = 2)
arrows(y_cure, par("usr")[3] + 0.25 * diff(par("usr")[3:4]), y_cure, 0, col="darkred", lwd=3, length=0.1)

# Panel b
y0 <- 0
y1 <- 12
plot(0, xlim = c(first_year_show, y_end ), ylim = c(y0, y1), col = "white", xlab = "", ylab = expression(atop("New re-infections", "(1/100,000 person-years)")), xaxt = "n", yaxt = "n", bty = "n", cex.lab=cexlab)
axis(1, at = seq(first_year_show-1, y_end, 1), labels = NA, col = "gray72", pos=0)
step <- 3
axis(2, at = seq(y0, y1, step), cex.axis=cexlab, col = "gray72", las=2)
text(x = seq(first_year_show, y_end, 2), y = par("usr")[3] - textxspace * diff(par("usr")[3:4]),
     labels = seq(first_year_show, y_end, 2),
     srt = 45, adj = c(0.5, 0.5), xpd = TRUE, cex = cexlab)
abline(h=seq(y0,y1+step,step), lty=3, col="gray72")
leg_sa <- rep("",n_sa)
for(i in 2:n_sa){
  lines(seq(years_sim[first_year_fit_idx]+1,years_sim[length(years_sim)],by=1/dt), mean_inc_reb[g,i,(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)] / perinc, col = colors2[i], lwd = moreh*morew*resclinew)
  polygon(c(seq(years_sim[first_year_fit_idx]+1,years_sim[length(years_sim)],by=1/dt), rev(seq(years_sim[first_year_fit_idx]+1,years_sim[length(years_sim)],by=1/dt))), c(q1_inc_reb[g,i,(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)], rev(q2_inc_reb[g,i,(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)])) / perinc, col = adjustcolor(colors2[i], alpha = 0.1), border = NA)
  leg_sa[i] <- paste0(round((1 - exp(- alphas))[i] * 100,0),"%")
}
i <- 1
lines(seq(years_sim[first_year_fit_idx]+1,years_sim[length(years_sim)],by=1/dt), mean_inc_reb[g,i,(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)] / perinc, col = colors2[i], lwd = moreh*morew*resclinew, lty=1)
polygon(c(seq(years_sim[first_year_fit_idx]+1,years_sim[length(years_sim)],by=1/dt), rev(seq(years_sim[first_year_fit_idx]+1,years_sim[length(years_sim)],by=1/dt))), c(q1_inc_reb[g,i,(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)], rev(q2_inc_reb[g,i,(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)])) / perinc, col = adjustcolor(colors2[i], alpha = 0.1), border = NA)
mtext("b", side = 3, line = 0.5, adj = 0.0, cex = cexlab*rescmtext, font = 2)
arrows(y_cure, par("usr")[3] + 0.25 * diff(par("usr")[3:4]), y_cure, 0, col="darkred", lwd=3, length=0.1)

# Panel c
step <- 0.02
y0 <- 0
y1 <- 0.08
plot(0, xlim = c(first_year_show, y_end ), ylim = c(y0, y1), col = "white", xlab = "", ylab = "Prevalence (%)", xaxt = "n", yaxt = "n", bty = "n", cex.lab=cexlab)
axis(1, at = seq(first_year_show-1, y_end, 1), labels = NA, col = "gray72", pos=0)
axis(2, at = round(seq(y0, y1, step),2), labels = 100*round(seq(y0, y1, step),2), cex.axis=cexlab, col = "gray72", las=2)
text(x = seq(first_year_show, y_end, 2), y = par("usr")[3] - textxspace * diff(par("usr")[3:4]),
     labels = seq(first_year_show, y_end, 2),
     srt = 45, adj = c(0.5, 0.5), xpd = TRUE, cex = cexlab)
abline(h=seq(y0,y1+step,step), lty=3, col="gray72")
leg_sa <- rep("",n_sa)
for(i in 2:n_sa){
  lines(seq(years_sim[first_year_fit_idx]+1,years_sim[length(years_sim)],by=1/dt), mean_prev_new[g,i,(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)] / pop + mean_prev_reb[g,i,(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)] / pop, col = colors2[i], lwd = moreh*morew*resclinew)
  polygon(c(seq(years_sim[first_year_fit_idx]+1,years_sim[length(years_sim)],by=1/dt), rev(seq(years_sim[first_year_fit_idx]+1,years_sim[length(years_sim)],by=1/dt))), c((q1_prev_new[g,i,(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)] + q1_prev_reb[g,i,(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)]), rev(q2_prev_new[g,i,(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)]+q2_prev_reb[g,i,(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)]))/pop, col = adjustcolor(colors2[i], alpha = 0.1), border = NA)
  leg_sa[i] <- paste0(round((1 - exp(- alphas))[i] * 100,0),"%")
}
i <- 1
lines(seq(years_sim[first_year_fit_idx]+1,years_sim[length(years_sim)],by=1/dt), mean_prev_new[g,i,(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)] / pop + mean_prev_reb[g,i,(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)] / pop, col = colors2[i], lwd = moreh*morew*resclinew, lty=1)
polygon(c(seq(years_sim[first_year_fit_idx]+1,years_sim[length(years_sim)],by=1/dt), rev(seq(years_sim[first_year_fit_idx]+1,years_sim[length(years_sim)],by=1/dt))), c((q1_prev_new[g,i,(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)] + q1_prev_reb[g,i,(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)]), rev(q2_prev_new[g,i,(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)]+q2_prev_reb[g,i,(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)]))/pop, col = adjustcolor(colors2[i], alpha = 0.1), border = NA)
mtext("c", side = 3, line = 0.5, adj = 0.0, cex = cexlab*rescmtext, font = 2)
arrows(y_cure, par("usr")[3] + 0.25 * diff(par("usr")[3:4]), y_cure, 0, col="darkred", lwd=3, length=0.1)

# Panel d
y0 <- 0
y1 <- 100
plot(0, xlim = c(first_year_show, y_end ), ylim = c(y0, y1), col = "white", xlab = "", ylab = "Cure coverage (%)", xaxt = "n", yaxt = "n", bty = "n", cex.lab=cexlab)
axis(1, at = seq(first_year_show-1, y_end, 1), labels = NA, col = "gray72", pos=0)
step <- 25
axis(2, at = seq(y0, y1, step), cex.axis=cexlab, col = "gray72", las=2)
text(x = seq(first_year_show, y_end, 2), y = par("usr")[3] - textxspace * diff(par("usr")[3:4]),
     labels = seq(first_year_show, y_end, 2),
     srt = 45, adj = c(0.5, 0.5), xpd = TRUE, cex = cexlab)
abline(h=seq(y0,y1+step,step), lty=3, col="gray72")
leg_sa <- rep("",n_sa)
for(i in 1:n_sa){
  if(i==1){
    leg_sa[i] <- "no cure"
  }else{
    lines(seq(years_sim[first_year_fit_idx]+1,years_sim[length(years_sim)],by=1/dt),round(100*mean_coverage[g,i,(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)], 1), col = colors2[i], lwd = moreh*morew*resclinew)
    polygon(c(seq(years_sim[first_year_fit_idx]+1,years_sim[length(years_sim)],by=1/dt), rev(seq(years_sim[first_year_fit_idx]+1,years_sim[length(years_sim)],by=1/dt))), c(q1_coverage[g,i,(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)], rev(q2_coverage[g,i,(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)])), col = adjustcolor(colors2[i], alpha = 0.1), border = NA)
    leg_sa[i] <- paste0(round((1 - exp(- alphas))[i] * 100,0),"%")
  }
}
i <- 1
lines(seq(years_sim[first_year_fit_idx]+1,years_sim[length(years_sim)],by=1/dt),round(100*mean_coverage[g,i,(dt*first_year_fit_idx+1):(dt*length(years_sim)-dt+1)], 1), col = colors2[i], lwd = moreh*morew*resclinew, lty=1)
legend("topleft", rev(leg_sa), lty = (c(rep(1,n_sa-1),2)), lwd = rep(moreh,n_sa), col = rev(colors2), bty = "n", cex = cexlab*rescleg, title = "Cure uptake (%)")
mtext("d", side = 3, line = 0.5, adj = 0.0, cex = cexlab*rescmtext, font = 2)
arrows(y_cure, par("usr")[3] + 0.25 * diff(par("usr")[3:4]), y_cure, 0, col="darkred", lwd=3, length=0.1)
mtext("Year", side=1, outer=TRUE, line=1., at=0.32, cex=cexlab*rescyear, adj = c(0.5, 0.5))  # Below the C panel
mtext("Year", side=1, outer=TRUE, line=1, at=0.82, cex=cexlab*rescyear, adj = c(0.5, 0.5))  # Below the D panel
dev.off()

if(hcfile == "MAIN" & y_cure == 2026){
  
  pdf(file = ("PLOTS/postcure/infectivity2.pdf"), width = morew*size_plot, height = moreh*size_plot * rescale)
  
  layout(matrix(c(1,1, 2,2, 3, 3,4,4), 2, 4, byrow = TRUE))
  par(mar = parmar * c(1.666, 1.1, 1, 1))
  par(mgp = parmgp)
  par(oma = paroma)
  
  p1_mat <- matrix(NA,nrow=4,ncol=nt)
  p2_naive_mat <- matrix(NA,nrow=4,ncol=nt)
  p2_nonnaive_mat <- matrix(NA,nrow=4,ncol=nt)
  p3_mat <- matrix(NA,nrow=4,ncol=nt)
  p_mat <- matrix(NA,nrow=4,ncol=nt)
  
  for(i in 2:4){
    p1_mat[i,] <- mean_inc_art_new[g,i,] / perinc
    p2_naive_mat[i,] <- mean_inc_undiag_naive_new[g,i,] / perinc
    p2_nonnaive_mat[i,] <- mean_inc_undiag_nonnaive_new[g,i,] / perinc
    p3_mat[i,] <- mean_inc_diag_new[g,i,] / perinc
    p_mat[i,] <- mean_inc_new[g,i,] / perinc
  }
  p1_mat_nocure <- rep(NA,nt)
  p2_mat_nocure <- rep(NA,nt)
  p3_mat_nocure <- rep(NA,nt)
  p_mat_nocure <- rep(NA,nt)
  p1_mat_nocure <- mean_inc_art_new[g,1,] / perinc
  p2_mat_nocure <- mean_inc_undiag_naive_new[g,1,] / perinc
  p3_mat_nocure <- mean_inc_diag_new[g,1,] / perinc
  p_mat_nocure <- mean_inc_new[g,1,] / perinc
  
  y_start_plot <- 2017
  panels <- c(NA, "a", "b", "c")
  titles <- c("Cure uptake: 10%", "Cure uptake: 50%", "Cure uptake: 90%")
  y0 <- 0
  y1 <- 1
  y1 <- 120
  y_end2 <- 2036
  plot(0, xlim = c(y_start_plot, y_end2 ), ylim = c(y0, y1), col = "white", xlab = "", ylab = "Contribution to\nnew infections by status\n(100,000 person-years)", 
       xaxt = "n", yaxt = "n", bty = "n", cex.lab=cexlab)
  axis(1, at = seq(y_start_plot-1, y_end2, 1), labels = NA, col = "gray72",pos=0)
  step <- 30
  axis(2, at = seq(y0, y1, step), labels = seq(y0, y1, step), cex.axis=cexlab, col = "gray72", las=2)
  text(x = seq(y_start_plot, y_end2, 2), y = par("usr")[3] - textxspace * diff(par("usr")[3:4]),
       labels = seq(y_start_plot, y_end2, 2),
       srt = 45, adj = c(0.5, 0.5), xpd = TRUE, cex = cexlab)
  abline(h=seq(y0,y1+step,step), lty=3, col="gray72")
  polygon(c(seq(years_sim[1],years_sim[length(years_sim)],by=1/dt), rev(seq(years_sim[1],years_sim[length(years_sim)],by=1/dt))), c(p1_mat_nocure, rev(rep(0,length(p_mat[1,])))), col = colors_inf[1], border = NA)
  polygon(c(seq(years_sim[1],years_sim[length(years_sim)],by=1/dt), rev(seq(years_sim[1],years_sim[length(years_sim)],by=1/dt))), c((p1_mat_nocure+p2_mat_nocure), rev(p1_mat_nocure)), col = colors_inf[2], border = NA)
  polygon(c(seq(years_sim[1],years_sim[length(years_sim)],by=1/dt), rev(seq(years_sim[1],years_sim[length(years_sim)],by=1/dt))), c((p1_mat_nocure+p2_mat_nocure+p3_mat_nocure), rev((p1_mat_nocure+p2_mat_nocure))), col = colors_inf[4], border = NA)
  mtext("a", side = 3, line = 0.5, adj = 0.0, cex = cexlab*rescmtext, font = 2)
  mtext("Year", side=1, outer=TRUE, line=-36.6, at=0.32, cex=cexlab*rescyear, adj = c(0.5, 0.5), font=1)  # Below the D panel
  mtext("No cure", side=1, outer=TRUE, line=-63.9, at=0.32, cex=cexlab*rescyear, adj = c(0.5, 0.5), font=1)  # Below the D panel
  
  for(pl in 2:4){
    
    plot(0, xlim = c(y_start_plot, y_end2 ), ylim = c(y0, y1), col = "white", xlab = "", ylab = "Contribution to\nnew infections by status\n(100,000 person-years)", 
         xaxt = "n", yaxt = "n", bty = "n", cex.lab=cexlab)
    axis(1, at = seq(y_start_plot-1, y_end2, 1), labels = NA, col = "gray72",pos=0)
    axis(2, at = seq(y0, y1, step), labels = seq(y0, y1, step), cex.axis=cexlab, col = "gray72", las=2)
    text(x = seq(y_start_plot, y_end2, 2), y = par("usr")[3] - textxspace * diff(par("usr")[3:4]),
         labels = seq(y_start_plot, y_end2, 2),
         srt = 45, adj = c(0.5, 0.5), xpd = TRUE, cex = cexlab)
    abline(h=seq(y0,y1+step,step), lty=3, col="gray72")
    polygon(c(seq(years_sim[1],years_sim[length(years_sim)],by=1/dt), rev(seq(years_sim[1],years_sim[length(years_sim)],by=1/dt))), c(p1_mat[pl,], rev(rep(0,length(p_mat[pl,])))), col = colors_inf[1], border = NA)
    polygon(c(seq(years_sim[1],years_sim[length(years_sim)],by=1/dt), rev(seq(years_sim[1],years_sim[length(years_sim)],by=1/dt))), c((p1_mat[pl,]+p2_naive_mat[pl,]), rev(p1_mat[pl,])), col = colors_inf[2], border = NA)
    polygon(c(seq(years_sim[1],years_sim[length(years_sim)],by=1/dt), rev(seq(years_sim[1],years_sim[length(years_sim)],by=1/dt))), c((p1_mat[pl,]+p2_naive_mat[pl,]+p2_nonnaive_mat[pl,]), rev(p1_mat[pl,]+p2_naive_mat[pl,])), col = colors_inf[3], border = NA)
    polygon(c(seq(years_sim[1],years_sim[length(years_sim)],by=1/dt), rev(seq(years_sim[1],years_sim[length(years_sim)],by=1/dt))), c((p1_mat[pl,]+p2_naive_mat[pl,]+p2_nonnaive_mat[pl,]+p3_mat[pl,]), rev((p1_mat[pl,]+p2_naive_mat[pl,]+p2_nonnaive_mat[pl,]))), col = colors_inf[4], border = NA)
    mtext(panels[pl], side = 3, line = 0.5, adj = 0.0, cex = cexlab*rescmtext, font = 2)
    arrows(y_cure, par("usr")[3] + 0.25 * diff(par("usr")[3:4]), y_cure, 0, col="darkred", lwd=3, length=0.1)
    if(pl==4){
      legend("topright", legend = c("Diagnosed", "Undiagnosed re-infected", "Undiagnosed naive", "Treated"), col = "white", pt.bg = rev(colors_inf), pch = rep(22,3), pt.cex = 4, cex = cexlab, bty = 'n')
    }
  }
  
  mtext("Year", side=1, outer=TRUE, line=-36.6, at=0.82, cex=cexlab*rescyear, adj = c(0.5, 0.5), font=1)  # Below the D panel
  mtext("Year", side=1, outer=TRUE, line=-4, at=0.32, cex=cexlab*rescyear, adj = c(0.5, 0.5), font=1)  # Below the D panel
  mtext("Year", side=1, outer=TRUE, line=-4, at=0.82, cex=cexlab*rescyear, adj = c(0.5, 0.5), font=1)  # Below the D panel
  
  mtext(titles[1], side=1, outer=TRUE, line=-63.9, at=0.82, cex=cexlab*rescyear, adj = c(0.5, 0.5), font=1)  # Below the D panel
  mtext(titles[2], side=1, outer=TRUE, line=-31.5, at=0.32, cex=cexlab*rescyear, adj = c(0.5, 0.5), font=1)  # Below the D panel
  mtext(titles[3], side=1, outer=TRUE, line=-31.5, at=0.82, cex=cexlab*rescyear, adj = c(0.5, 0.5), font=1)  # Below the D panel
  
  dev.off()
  
  bp_caa <- rowMeans(cumulative_art_averted[2,,])
  bp_caa_lo <- apply(cumulative_art_averted[2,,], 1, quantile, prob = 0.025)
  bp_caa_hi <- apply(cumulative_art_averted[2,,], 1, quantile, prob = 0.975)
  
  pdf(file = ("PLOTS/postcure/artaverted2.pdf"), width = morew*size_plot, height = moreh*size_plot * rescale)
  layout(matrix(c(1,1, 2,2, 3, 3,4,4), 2, 4, byrow = TRUE))
  par(mar = parmar * c(1.666, 1, 1, 1))
  par(mgp = parmgp)
  par(oma = paroma)
  
  plot(0, xlim = c(0.5, 3.5 ), ylim = c(y0, y1+step), col = "white", xlab = "", ylab = expression(atop("Cumulative ART", "person-years averted (%)")), xaxt = "n", yaxt = "n", bty = "n", cex.lab=cexlab)
  plot(0, xlim = c(0.5, 3.5 ), ylim = c(y0, y1+step), col = "white", xlab = "", ylab = expression(atop("Cumulative ART", "person-years averted (%)")), xaxt = "n", yaxt = "n", bty = "n", cex.lab=cexlab)
  plot(0, xlim = c(0.5, 3.5 ), ylim = c(y0, y1+step), col = "white", xlab = "", ylab = expression(atop("Cumulative ART", "person-years averted (%)")), xaxt = "n", yaxt = "n", bty = "n", cex.lab=cexlab)
  
  # Panel f
  y0 <- 0
  step <- 0.2
  y1 <- max(bp_caa_hi)
  if(hcfile=="SA"){
    y1 <- y1*2
    step <- step*2
  }
  plot(0, xlim = c(0.5, 3.5 ), ylim = c(y0, y1+step), col = "white", xlab = "", ylab = expression(atop("Cumulative ART", "person-years averted (%)")), xaxt = "n", yaxt = "n", bty = "n", cex.lab=cexlab)
  axis(1, at = seq(1, 3, 1), labels = NA, col = "gray72",pos=0)
  axis(2, at = (seq(y0, y1+step, step)), labels = 100*(seq(y0, y1+step, step)), cex.axis=cexlab, col = "gray72", las=2)
  text(x = seq(1, 3, 1), y = par("usr")[3] - textxspace * diff(par("usr")[3:4]),
       labels = c("10", "50", "90"),
       srt = 0, adj = c(0.5, 0.5), xpd = TRUE, cex = cexlab)
  abline(h=seq(y0,y1+step,step), lty=3, col="gray72")
  for(i in 1:(n_sa-1))){
    rect(xleft=i-0.48, ybottom = 0, xright = i+0.48, ytop=bp_caa[i], col = colors2_reduct[i], border = "white")
    individual_points <- cumulative_art_averted[2, i, ]
    points(x = i + runif(length(individual_points), -0.47, 0.47),
           y = individual_points,
           pch = 21, bg = "magenta2", col = "magenta2", cex = 0.225)
    segments(x0 = i, y0 = bp_caa_lo[i], x1 = i, y1 = bp_caa_hi[i], col = adjustcolor("darkred",alpha.f = 1), lwd = 3)
  }
  mtext("b", side = 3, line = 0.5, adj = 0.0, cex = cexlab*rescmtext, font = 2)
  arrows(y_cure, par("usr")[3] + 0.25 * diff(par("usr")[3:4]), y_cure, 0, col="darkred", lwd=3, length=0.1)
  mtext("Cure uptake (%)", side=1, outer=TRUE, line=-4, at=0.82, cex=cexlab*rescyear, adj = c(0.5, 0.5), font=1)  # Below the D panel
  
  dev.off()
  
  pdf(file = ("PLOTS/postcure/artpercure2.pdf"), width = morew*size_plot, height = moreh*size_plot * rescale)
  
  layout(matrix(c(1,1, 2,2, 3, 3,4,4), 2, 4, byrow = TRUE))
  par(mar = parmar * c(1.666, 1, 1, 1))
  par(mgp = parmgp)
  par(oma = paroma)
  
  bp_caa_percure <- rowMeans(cumulative_art_averted_percure[2,,])
  bp_caa_percure_lo <- apply(cumulative_art_averted_percure[2,,], 1, quantile, prob = 0.025)
  bp_caa_percure_hi <- apply(cumulative_art_averted_percure[2,,], 1, quantile, prob = 0.975)
  
  plot(0, xlim = c(0.5, 3.5 ), ylim = c(y0, y1+step), col = "white", xlab = "", ylab = expression(atop("Cumulative ART person-years", " averted per cure given")), xaxt = "n", yaxt = "n", bty = "n", cex.lab=cexlab)
  plot(0, xlim = c(0.5, 3.5 ), ylim = c(y0, y1+step), col = "white", xlab = "", ylab = expression(atop("Cumulative ART person-years", " averted per cure given")), xaxt = "n", yaxt = "n", bty = "n", cex.lab=cexlab)
  plot(0, xlim = c(0.5, 3.5 ), ylim = c(y0, y1+step), col = "white", xlab = "", ylab = expression(atop("Cumulative ART person-years", " averted per cure given")), xaxt = "n", yaxt = "n", bty = "n", cex.lab=cexlab)
  
  # Panel f
  y0 <- 0
  step <- 1
  y1 <- 4
  plot(0, xlim = c(0.5, 3.5 ), ylim = c(y0, y1+step), col = "white", xlab = "", ylab = expression(atop("Cumulative ART person-years", " averted per cure given")), xaxt = "n", yaxt = "n", bty = "n", cex.lab=cexlab)
  axis(1, at = seq(1, 3, 1), labels = NA, col = "gray72",pos=0)
  axis(2, at = (seq(y0, y1+step, step)), cex.axis=cexlab, col = "gray72", las=2)
  text(x = seq(1, 3, 1), y = par("usr")[3] - textxspace * diff(par("usr")[3:4]),
       labels = c("10", "50", "90"),
       srt = 0, adj = c(0.5, 0.5), xpd = TRUE, cex = cexlab)
  abline(h=seq(y0,y1+step,step), lty=3, col="gray72")
  for(i in 1:n_sa){
    rect(xleft=i-0.48, ybottom = 0, xright = i+0.48, ytop=bp_caa_percure[i], col = colors2_reduct[i], border = "white")
    segments(x0 = i, y0 = bp_caa_percure_lo[i], x1 = i, y1 = bp_caa_percure_hi[i], col = adjustcolor("darkred",alpha.f = 1), lwd = 3)
  }
  arrows(y_cure, par("usr")[3] + 0.25 * diff(par("usr")[3:4]), y_cure, 0, col="darkred", lwd=3, length=0.1)
  mtext("Cure uptake (%)", side=1, outer=TRUE, line=-4, at=0.82, cex=cexlab*rescyear, adj = c(0.5, 0.5), font=1)  # Below the D panel
  mtext("d", side = 3, line = 0.5, adj = 0.0, cex = cexlab*rescmtext, font = 2)
  
  dev.off()
  
  # 10%
  p1_mat[2,length(p1_mat[1,])] / p_mat[2,length(p1_mat[1,])]
  p2_naive_mat[2,length(p1_mat[1,])] / p_mat[2,length(p1_mat[1,])]
  p2_nonnaive_mat[2,length(p1_mat[1,])] / p_mat[2,length(p1_mat[1,])]
  p3_mat[2,length(p1_mat[1,])] / p_mat[2,length(p1_mat[1,])]
  # 50%
  p1_mat[3,length(p1_mat[1,])] / p_mat[3,length(p1_mat[1,])]
  p2_naive_mat[3,length(p1_mat[1,])] / p_mat[3,length(p1_mat[1,])]
  p2_nonnaive_mat[3,length(p1_mat[1,])] / p_mat[3,length(p1_mat[1,])]
  p3_mat[3,length(p1_mat[1,])] / p_mat[3,length(p1_mat[1,])]
  # 90%
  p1_mat[4,length(p1_mat[1,])] / p_mat[4,length(p1_mat[1,])]
  p2_naive_mat[4,length(p1_mat[1,])] / p_mat[4,length(p1_mat[1,])]
  p2_nonnaive_mat[4,length(p1_mat[1,])] / p_mat[4,length(p1_mat[1,])]
  p3_mat[4,length(p1_mat[1,])] / p_mat[4,length(p1_mat[1,])]
  
  # no cure
  p1_mat_nocure[length(p1_mat[1,])] / p_mat_nocure[length(p1_mat[1,])]
  p2_mat_nocure[length(p1_mat[1,])] / p_mat_nocure[length(p1_mat[1,])]
  p3_mat_nocure[length(p1_mat[1,])] / p_mat_nocure[length(p1_mat[1,])]
  
}


# no cure
mean_inc_new[g,1,length(mean_inc_new[g,1,])] / perinc
q1_inc_new[g,1,length(mean_inc_new[g,1,])] / perinc
q2_inc_new[g,1,length(mean_inc_new[g,1,])] / perinc

# 10%
mean_inc_new[g,2,length(mean_inc_new[g,1,])] / perinc
q1_inc_new[g,2,length(mean_inc_new[g,1,])] / perinc
q2_inc_new[g,2,length(mean_inc_new[g,1,])] / perinc

# 50%
mean_inc_new[g,3,length(mean_inc_new[g,1,])] / perinc
q1_inc_new[g,3,length(mean_inc_new[g,1,])] / perinc
q2_inc_new[g,3,length(mean_inc_new[g,1,])] / perinc

# 90%
mean_inc_new[g,4,length(mean_inc_new[g,1,])] / perinc
q1_inc_new[g,4,length(mean_inc_new[g,1,])] / perinc
q2_inc_new[g,4,length(mean_inc_new[g,1,])] / perinc
