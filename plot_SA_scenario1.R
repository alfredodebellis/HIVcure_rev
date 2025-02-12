rm(list = ls())

library(patchwork)
library("ggplot2")

load("PLOTS/colors.Rdata")
  
efficaciessa <- c(0.2,0.9)
leneff <- length(efficaciessa)

arssa <- c(0.1, 0.5, 0.9)
lenars <- length(arssa)

failuressa <- c(0, 0.16666666, 0.5)
lenfail <- length(failuressa)

n_monitors <- 3
gs <- 1:n_monitors

sa_matrix <- array(0,dim=c(n_monitors,leneff,lenars))
sa_matrix_lo <- array(0,dim=c(n_monitors,leneff,lenars))
sa_matrix_hi <- array(0,dim=c(n_monitors,leneff,lenars))
sa_matrix_reb <- array(0,dim=c(n_monitors,leneff,lenars))
sa_matrix_reb_lo <- array(0,dim=c(n_monitors,leneff,lenars))
sa_matrix_reb_hi <- array(0,dim=c(n_monitors,leneff,lenars))

rescleg <- 1.4
resccolh <- 3
resccolw <- 4
rescmarg <- 2
resclab <- resccolw

load(paste0("RESULTS/SCENARIOS/SC1_eff90_cu90_hcMAIN_curestart2026.RData"))

years_sim <- years_calibration[1]:y_end
curet <- which(years_sim==y_cure)
cureidx <- (curet-1)*dt + 1

# titles <- c("Standard of care", "Monitoring", "Follow-up")
titles <- c("28 months", "3 months", "2 weeks")
subset_nocure_incidence <- incidence_nocure_time[,cureidx:nt]

perinc <- 210000 / 1e5
cirnc <- array(apply(subset_nocure_incidence, c(1), sum), dim=c(len_post)) / dt
mean_cirnc <- mean(cirnc) / perinc
lo_cirnc <- quantile(cirnc, 0.025) / perinc
hi_cirnc <- quantile(cirnc, 0.975) / perinc

# cumulative_nocure_incidence <- apply(subset_nocure_incidence, 2, sum)
cumulative_nocure_incidence <- apply(subset_nocure_incidence, 1, sum)
mean_cum_nocure <- mean(cumulative_nocure_incidence)
lo_cum_nocure <- quantile(cumulative_nocure_incidence,0.025)
hi_cum_nocure <- quantile(cumulative_nocure_incidence,0.975)
save(mean_cum_nocure, file = ("RESULTS/cumincnew.RData"))
save(cumulative_nocure_incidence, file = ("RESULTS/cumincnew_post.RData"))
limfix <- 61
limfix_app <- 60
limcol1 <- 100
limcol2 <- 200
limfix_reb <- 5000
limcol1_reb <- 10000
limcol2_reb <- 1000000


ratio_man <- 0.73
rescrev <- 1.5

matrices_heatmap <- list()
matrices_heatmap_reb <- list()
matrices_heatmap_lo <- list()
matrices_heatmap_reb_lo <- list()
matrices_heatmap_hi <- list()
matrices_heatmap_reb_hi <- list()
plots_heatmap <- list()
plots_heatmap_reb <- list()
idxplot <- 1
for(kkk in 1:lenfail){
  failure <- failuressa[kkk]
  idxfailure <- which(failuressa==failure)
  for (iii in 1:leneff) {
    efficacy <- efficaciessa[iii]
    for (jjj in 1:lenars){
      ar <- arssa[jjj]
      load(paste0("RESULTS/SCENARIOS/SC1_eff", round(efficacy*100,0),"_cu",round(ar*100,0),"_hcMAIN_curestart2026.RData"))
      mean_inc_new <- array(apply(incidence_after_new_time_pool, c(1, 2, 4), mean), dim=c(n_monitors, n_sa, nt))
      mean_inc_reb <- array(apply(incidence_after_reb_time_pool, c(1, 2, 4), mean), dim=c(n_monitors, n_sa, nt))
      q1_inc_new <- array(apply(incidence_after_new_time_pool, c(1, 2, 4), function(x) quantile(x, 0.025)), dim=c(n_monitors, n_sa, nt))
      q2_inc_new <- array(apply(incidence_after_new_time_pool, c(1, 2, 4), function(x) quantile(x, 0.975)), dim=c(n_monitors, n_sa, nt))
      q1_inc_reb <- array(apply(incidence_after_reb_time_pool, c(1, 2, 4), function(x) quantile(x, 0.025)), dim=c(n_monitors, n_sa, nt))
      q2_inc_reb <- array(apply(incidence_after_reb_time_pool, c(1, 2, 4), function(x) quantile(x, 0.975)), dim=c(n_monitors, n_sa, nt))
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
      cumulative_incidence_reb <- array(apply(subset_incidence_reb, c(1, 2, 3), sum), dim=c(n_monitors, n_sa, len_post)) / dt
      n_time_points <- dim(subset_incidence_reb)[4]
      mean_cum_reb <- array(apply(cumulative_incidence_reb, c(1, 2), mean), dim=c(n_monitors, n_sa))
      lo_cum_reb <- array(apply(cumulative_incidence_reb, c(1, 2), quantile, prob=0.025), dim=c(n_monitors, n_sa))
      hi_cum_reb <- array(apply(cumulative_incidence_reb, c(1, 2), quantile, prob=0.975), dim=c(n_monitors, n_sa))
      for(g in 1:n_monitors){
        tmpnew <- mean_cum_new[g,idxfailure]
        changerelnew <- 100 * (tmpnew - mean_cum_nocure) / mean_cum_nocure
        sa_matrix[g,iii,jjj] <- mean_cum_new[g,idxfailure]
        sa_matrix_lo[g,iii,jjj] <- lo_cum_new[g,idxfailure]
        sa_matrix_hi[g,iii,jjj] <- hi_cum_new[g,idxfailure]
        tmpreb <- mean_cum_reb[g,idxfailure]
        sa_matrix_reb[g,iii,jjj] <- tmpreb
        sa_matrix_reb_lo[g,iii,jjj] <- lo_cum_reb[g,idxfailure]
        sa_matrix_reb_hi[g,iii,jjj] <- hi_cum_reb[g,idxfailure]
      }
    }
  }
  
  sa_matrix_reb <- sa_matrix_reb / perinc
  sa_matrix_reb_lo <- sa_matrix_reb_lo / perinc
  sa_matrix_reb_hi <- sa_matrix_reb_hi / perinc
  matrices_heatmap[[kkk]] <- sa_matrix
  matrices_heatmap_lo[[kkk]] <- sa_matrix_lo
  matrices_heatmap_hi[[kkk]] <- sa_matrix_hi
  matrices_heatmap_reb[[kkk]] <- sa_matrix_reb
  matrices_heatmap_reb_lo[[kkk]] <- sa_matrix_reb_lo
  matrices_heatmap_reb_hi[[kkk]] <- sa_matrix_reb_hi
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
    sa_df_filtered0 <- sa_df[sa_df$z > limfix & sa_df$z < limcol1, ]
    sa_df_filtered1 <- sa_df[sa_df$z > limcol1 & sa_df$z < limcol2, ]
    sa_df_filtered2 <- sa_df[sa_df$z > limcol2, ]
    sa_df_filtered0_reb <- sa_df_reb[sa_df_reb$z > limfix_reb & sa_df_reb$z < limcol1_reb, ]
    sa_df_filtered1_reb <- sa_df_reb[sa_df_reb$z > limcol1_reb & sa_df_reb$z < limcol2_reb, ]
    sa_df_filtered2_reb <- sa_df_reb[sa_df_reb$z > limcol2_reb, ]
    heatmap <- ggplot(data = sa_df, aes(x = x, y = y, fill = z)) +
      geom_tile() +
      scale_fill_gradientn(colours = colcol_new, limits = c(-limfix, limfix)
                           , name = "", na.value = "grey25") +
      labs(title = if(kkk == 1) titles[g] else NULL, x = if (kkk == lenfail) "Cure efficacy (%)" else NULL, y = if (g == 1) "Cure uptake (%)" else NULL) +
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
      coord_fixed(ratio = ratio_man) +
      annotate("text", x = sa_df_filtered0$x, y = sa_df_filtered0$y, label = paste0("> ", limfix_app,"%"), hjust = 0.5, vjust = 0.5, size = size_plot*rescleg*rescrev, colour="white")+
      annotate("text", x = sa_df_filtered1$x, y = sa_df_filtered1$y, label = paste0("> ", limcol1,"%"), hjust = 0.5, vjust = 0.5, size = size_plot*rescleg*rescrev, colour="white")+
      annotate("text", x = sa_df_filtered2$x, y = sa_df_filtered2$y, label = paste0("> ", limcol2,"%"), hjust = 0.5, vjust = 0.5, size = size_plot*rescleg*rescrev, colour="white")

    heatmap_reb <- ggplot(data = sa_df_reb, aes(x = x, y = y, fill = z)) +
      geom_tile() +
      scale_fill_gradientn(colours = colcol_reb, limits = c(0, limfix_reb)
                           , name = "", na.value = "grey25") +
      labs(title = if(kkk == 1) titles[g] else NULL, x = if (kkk == lenfail) "Cure efficacy (%)" else NULL, y = if (g == 1) "Cure uptake (%)" else NULL) +
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
      coord_fixed(ratio = ratio_man) +
      annotate("text", x = sa_df_filtered0_reb$x, y = sa_df_filtered0_reb$y, label = paste0("> ", limfix_reb), hjust = 0.5, vjust = 0.5, size = size_plot*rescleg*rescrev, colour="white")+
      annotate("text", x = sa_df_filtered1_reb$x, y = sa_df_filtered1_reb$y, label = paste0("> ", limcol1_reb), hjust = 0.5, vjust = 0.5, size = size_plot*rescleg*rescrev, colour="white")+
      annotate("text", x = sa_df_filtered2_reb$x, y = sa_df_filtered2_reb$y, label = paste0("> ", limcol2_reb), hjust = 0.5, vjust = 0.5, size = size_plot*rescleg*rescrev, colour="white")
    plots_heatmap[[idxplot]] <- heatmap
    plots_heatmap_reb[[idxplot]] <- heatmap_reb
    idxplot <- idxplot + 1
  }
}

all_plots_new <- c(plots_heatmap)
all_plots_reb <- c(plots_heatmap_reb)
combined_plot_new <- wrap_plots(all_plots_new, ncol = 3, nrow = 3)
combined_plot_reb <- wrap_plots(all_plots_reb, ncol = 3, nrow = 3)
ggsave("PLOTS/postcure/Fig3A.pdf", plot = combined_plot_new, width = size_plot * 3, height = size_plot * 3)
ggsave("PLOTS/postcure/Fig3B.pdf", plot = combined_plot_reb, width = size_plot * 3, height = size_plot * 3)


monitoring_strategies <- c("28 months", "3 months", "2 weeks")
time_to_rebound <- c("never", "6 years", "2 years")
cure_uptake_labels <- c("90%", "50%", "10%")  # Rows
cure_efficacy_labels <- c("20%", "90%")       # Columns

matrices <- list()
for (gg in 1:length(monitoring_strategies)) {
  matrices[[monitoring_strategies[gg]]] <- list()
  
  for (ttr in 1:length(time_to_rebound)) {
    matrices[[monitoring_strategies[gg]]][[time_to_rebound[ttr]]] <- 
      t(matrices_heatmap[[ttr]][gg,,])[rev(1:3),]
  }
}

print(matrices[["2 weeks"]][["2 years"]])  # Example matrix output

cure_uptake_labels <- c("90%", "50%", "10%")  # Rows
cure_efficacy_labels <- c("20%", "90%")       # Columns

results <- list()
for (ttr in time_to_rebound) {
  for (i in 1:3) {  # Rows (Cure Uptake: 90%, 50%, 10%)
    for (j in 1:2) {  # Columns (Cure Efficacy: 20%, 90%)
      found <- FALSE
      for (gg in monitoring_strategies) {
        mat <- matrices[[gg]][[ttr]]
        if (mat[i, j] < 0) {
          results <- append(results, list(data.frame(
            "Time to Rebound" = ttr,
            "Cure Uptake (%)" = cure_uptake_labels[i],
            "Cure Efficacy (%)" = cure_efficacy_labels[j],
            "Minimum Monitoring" = gg
          )))
          found <- TRUE
          break  
        }
      }
      if (!found) {
        results <- append(results, list(data.frame(
          "Time to Rebound" = ttr,
          "Cure Uptake (%)" = cure_uptake_labels[i],
          "Cure Efficacy (%)" = cure_efficacy_labels[j],
          "Minimum Monitoring" = "--"
        )))
      }
    }
  }
}

df_results <- do.call(rbind, results)
print(df_results)
write.csv(df_results, "min_monitoring_strategy.csv", row.names = FALSE)

