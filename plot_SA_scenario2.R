rm(list = ls())
library(patchwork)
library("ggplot2")


load("PLOTS/colors.Rdata")

efficaciessa <- c(0.2,0.9)
leneff <- length(efficaciessa)

arssa <- c(0.1, 0.5, 0.9)
lenars <- length(arssa)

n_monitors <- 3
gs <- 1:n_monitors

sa_matrix <- array(0,dim=c(n_monitors,leneff,lenars))
sa_matrix_reb <- array(0,dim=c(n_monitors,leneff,lenars))
sa_matrix_lo <- array(0,dim=c(n_monitors,leneff,lenars))
sa_matrix_reb_lo <- array(0,dim=c(n_monitors,leneff,lenars))
sa_matrix_hi <- array(0,dim=c(n_monitors,leneff,lenars))
sa_matrix_reb_hi <- array(0,dim=c(n_monitors,leneff,lenars))

# titles <- c("Standard of care", "Monitoring", "Follow-up")
titles <- c("28 months", "3 months", "2 weeks")

load("RESULTS/cumincnew.RData")
load("RESULTS/cumincnew_post.RData")

limfix <- 61
limfix_reb <- 25

rescleg <- 1.5
resccolh <- 3
resccolw <- 4
rescmarg <- 2
resclab <- resccolw*1.1

perinc <- 210000 / 1e5

for (iii in 1:leneff) {
  efficacy <- efficaciessa[iii]
  
  load(paste0("RESULTS/SCENARIOS/SC2_eff", round(efficacy*100,0), "_hcMAIN_curestart2026.RData"))
  years_sim <- years_calibration[1]:y_end

  mean_inc_new <- array(apply(incidence_after_new_time_pool, c(1, 2, 4), mean), dim=c(n_monitors, n_sa, nt))
  mean_inc_reb <- array(apply(incidence_after_reb_time_pool, c(1, 2, 4), mean), dim=c(n_monitors, n_sa, nt))
  q1_inc_new <- array(apply(incidence_after_new_time_pool, c(1, 2, 4), function(x) quantile(x, 0.025)), dim=c(n_monitors, n_sa, nt))
  q2_inc_new <- array(apply(incidence_after_new_time_pool, c(1, 2, 4), function(x) quantile(x, 0.975)), dim=c(n_monitors, n_sa, nt))
  q1_inc_reb <- array(apply(incidence_after_reb_time_pool, c(1, 2, 4), function(x) quantile(x, 0.025)), dim=c(n_monitors, n_sa, nt))
  q2_inc_reb <- array(apply(incidence_after_reb_time_pool, c(1, 2, 4), function(x) quantile(x, 0.975)), dim=c(n_monitors, n_sa, nt))
  
  curet <- which(years_sim==y_cure)
  cureidx <- (curet-1)*dt + 1
  subset_incidence_new <- incidence_after_new_time_pool[,,,cureidx:nt]
  cumulative_incidence_new <- array(apply(subset_incidence_new, c(1, 2, 3), sum), dim=c(n_monitors, n_sa, len_post))
  
  cumulative_incidence_diff <- cumulative_incidence_new
  for (p in 1:len_post) {
    cumulative_incidence_diff[,,p] <- 100*(cumulative_incidence_new[,,p] - cumulative_nocure_incidence[p]) / cumulative_nocure_incidence[p]
  }
  
  mean_cum_new <- array(apply(cumulative_incidence_diff, c(1, 2), mean), dim=c(n_monitors, n_sa))
  lo_cum_new <- array(apply(cumulative_incidence_diff, c(1, 2), quantile, prob=0.025), dim=c(n_monitors, n_sa))
  hi_cum_new <- array(apply(cumulative_incidence_diff, c(1, 2), quantile, prob=0.975), dim=c(n_monitors, n_sa))
  time_diff <- 1/dt
  subset_incidence_reb <- incidence_after_reb_time_pool[,,,cureidx:nt]
  n_time_points <- dim(subset_incidence_reb)[4]
  cumulative_incidence_reb <- array(apply(subset_incidence_reb, c(1, 2, 3), sum), dim=c(n_monitors, n_sa, len_post)) / dt

  mean_cum_reb <- array(apply(cumulative_incidence_reb, c(1, 2), mean), dim=c(n_monitors, n_sa))
  lo_cum_reb <- array(apply(cumulative_incidence_reb, c(1, 2), quantile, prob=0.025), dim=c(n_monitors, n_sa))
  hi_cum_reb <- array(apply(cumulative_incidence_reb, c(1, 2), quantile, prob=0.975), dim=c(n_monitors, n_sa))
  for(g in 1:n_monitors){
    for(jjj in 2:n_sa){
      tmpnew <- mean_cum_new[g,jjj]
      changerelnew <- 100 * (tmpnew - mean_cum_nocure) / mean_cum_nocure
      sa_matrix[g,iii,jjj-1] <- mean_cum_new[g,jjj]
      sa_matrix_lo[g,iii,jjj-1] <- lo_cum_new[g,jjj]
      sa_matrix_hi[g,iii,jjj-1] <- hi_cum_new[g,jjj]
      tmpreb <- mean_cum_reb[g,jjj]
      sa_matrix_reb[g,iii,jjj-1] <- tmpreb
      sa_matrix_reb_lo[g,iii,jjj-1] <- lo_cum_reb[g,jjj]
      sa_matrix_reb_hi[g,iii,jjj-1] <- hi_cum_reb[g,jjj]
    }
  }
}

rescrev <- 1.2
sa_matrix_reb <- sa_matrix_reb / perinc
sa_matrix_reb_lo <- sa_matrix_reb_lo / perinc
sa_matrix_reb_hi <- sa_matrix_reb_hi / perinc
ratio_man <- 0.73
plots_heatmap <- list()
plots_heatmap_reb <- list()
for (g in 1:n_monitors) {
  matmat <- sa_matrix[g,,]
  sa_df <- expand.grid(x = 1:nrow(matmat), y = 1:ncol(matmat))
  sa_df$z <- as.vector(matmat)
  matmat_reb <- sa_matrix_reb[g,,]
  sa_df_reb <- expand.grid(x = 1:nrow(matmat_reb), y = 1:ncol(matmat_reb))
  sa_df_reb$z <- as.vector(matmat_reb)
  limax <- max(max(matmat), (-min(matmat)))
  limax_reb <- max(max(matmat_reb), (-min(matmat_reb)))
  lby <- round(100*arssa,0)
  lbx <- round(100*efficaciessa,0)
  heatmap <- ggplot(data = sa_df, aes(x = x, y = y, fill = z)) +
    geom_tile() +
    scale_fill_gradientn(colours = colcol_new, limits = c(-limfix, limfix)
                         , name = "") +
    labs(title = titles[g], x = "Cure efficacy (%)", y = if (g == 1) "Cure uptake (%)" else NULL) +
    theme_minimal() +
    theme(
      axis.line = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.key.width = unit(size_plot/resccolw, "cm"),  
      legend.key.height = unit(size_plot/resccolh, "cm"), 
      legend.text = element_text(size = size_plot*resclab*rescrev),  
      legend.position = if (g == n_monitors) "right" else "none",  
      axis.text = element_text(size = size_plot*resclab*rescrev),  
      axis.title = element_text(size = size_plot*resclab*rescrev),
      axis.title.y = element_text(margin = margin(r = size_plot*rescmarg*rescrev)),
      axis.title.x = element_text(margin = margin(t = size_plot*rescmarg*rescrev)),
      plot.title = element_text(size = size_plot*resclab*rescrev), 
    ) +
    scale_x_continuous(breaks = 1:length(lbx), labels = lbx) +
    scale_y_continuous(breaks = 1:length(lby), labels = lby) +
    coord_fixed(ratio = ratio_man)  # Add fixed aspect ratio 
  heatmap_reb <- ggplot(data = sa_df_reb, aes(x = x, y = y, fill = z)) +
    geom_tile() +
    scale_fill_gradientn(colours = colcol_reb, limits = c(0, limfix_reb)
                         , name = "", na.value = "grey25") +
    labs(title = titles[g], x = "Cure efficacy (%)", y = if (g == 1) "Cure uptake (%)" else NULL) +
    theme_minimal() +
    theme(
      axis.line = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.key.width = unit(size_plot/resccolw, "cm"),
      legend.key.height = unit(size_plot/resccolh, "cm"), 
      legend.text = element_text(size = size_plot*resclab*rescrev),
      legend.position = if (g == n_monitors) "right" else "none", 
      axis.text = element_text(size = size_plot*resclab*rescrev),  
      axis.title = element_text(size = size_plot*resclab*rescrev),
      axis.title.y = element_text(margin = margin(r = size_plot*rescmarg*rescrev)),
      axis.title.x = element_text(margin = margin(t = size_plot*rescmarg*rescrev)),
      plot.title = element_text(size = size_plot*resclab*rescrev), 
    ) +
    scale_x_continuous(breaks = 1:length(lbx), labels = lbx) +
    scale_y_continuous(breaks = 1:length(lby), labels = lby) +
    coord_fixed(ratio = ratio_man)  # Add fixed aspect ratio
  plots_heatmap[[g]] <- heatmap
  plots_heatmap_reb[[g]] <- heatmap_reb
}

t(sa_matrix[1,,3:1])
t(sa_matrix[2,,3:1])
t(sa_matrix[3,,3:1])
min(sa_matrix)
max(sa_matrix)
min(sa_matrix_reb)
max(sa_matrix_reb)

combined_plot_heatmap <- wrap_plots(plots_heatmap, ncol = 3)
combined_plot_heatmap_reb <- wrap_plots(plots_heatmap_reb, ncol = 3)
all_plots <- c(plots_heatmap, plots_heatmap_reb)
combined_plot_tot <- wrap_plots(all_plots, ncol = 3, nrow = 2)
ggsave("PLOTS/postcure/Fig5A.pdf", plot = combined_plot_heatmap, width = size_plot * 3, height = size_plot)
ggsave("PLOTS/postcure/Fig5B.pdf", plot = combined_plot_heatmap_reb, width = size_plot * 3, height = size_plot)

