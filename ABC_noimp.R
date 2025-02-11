rm(list = ls())

library(deSolve)
library(lhs)
library(readxl)
library(dplyr)
library(tidyr)

# Model function for scenario 0
model <- function(t, x, parms) {
  with(as.list(c(parms, x)), {
    
    #parameters
    beta <- parms["beta"]
    gamma <- parms["gamma"]
    gamma_diag <- parms["gamma_diag"]
    mu <- parms["mu"]
    lambda <- parms["lambda"]
    eps <- parms["eps"]
    eps_p <- parms["eps_p"]
    eps_d <- parms["eps_d"]
    eta <- parms["eta"]
    omega <- parms["omega"]
    Omega <- parms["Omega"]
    NN <- parms["NN"]
    imp <- parms["imp"]
    rho <- sapply(1:numstages, function(k) parms[paste0("rho", k)])
    h <- sapply(1:numstages, function(k) parms[paste0("h", k)])
    tau_p <- sapply(1:numgroups, function(l) parms[paste0("tau_p", l)])
    tau <- sapply(1:numstages, function(k) parms[paste0("tau", k)])
    q <- sapply(1:numgroups, function(l) parms[paste0("q", l)])
    c <- sapply(1:numgroups, function(l) parms[paste0("c", l)])
    kon <- sapply(1:numgroups, function(l) parms[paste0("kon", l)])
    kof <- sapply(1:numgroups, function(l) parms[paste0("kof", l)])
    qadj <- sapply(1:numgroups, function(l) parms[paste0("qadj", l)])
    
    #variables
    S <- sapply(1:numgroups, function(l) x[paste0("S", l)])
    Sp <- sapply(1:numgroups, function(l) x[paste0("Sp", l)])
    Ip <- sapply(1:numgroups, function(l) x[paste0("Ip", l)])
    I <- matrix(sapply(1:(numgroups * numstages), function(i) x[paste0("I", i)]), nrow = numgroups, ncol = numstages)
    A <- sapply(1:numgroups, function(l) x[paste0("A", l)])
    D <- sapply(1:numgroups, function(l) x[paste0("D", l)])
    
    Ng <- rowSums(cbind(S, Sp, rowSums(I), Ip, A, D))
    N <- sum(Ng)
    
    # Calculate the mixing matrix
    M <- matrix(0,numgroups,numgroups)
    denM <- sum(c*Ng)
    for(l in 1:numgroups){
      for (ll in 1:numgroups) {
        M[l,ll] <- omega * c[ll]*Ng[ll] / denM
        if(l==ll){
          M[l,ll] <- M[l,ll] + (1 - omega)
        }
      }
    }
    
    # Calculate the force of infection
    J <- rep(0,numgroups)
    for(l in 1:numgroups){
      for (ll in 1:numgroups) {
        J[l] <- J[l] + M[l, ll] * (eps_p * Ip[ll] / Ng[ll] + (eps * A[ll] + eps_d * D[ll]) / Ng[ll])
        for (k in 1:numstages) {
          J[l] <- J[l] + M[l, ll] * h[k] * (I[ll, k]) / Ng[ll]
        }
      }
      J[l] <- J[l] * (lambda * c[l])
    }
    
    dS <- rep(0,numgroups)
    dSp <- rep(0,numgroups)
    dI <- matrix(0,numgroups,numstages)
    dIp <- rep(0,numgroups)
    dA <- rep(0,numgroups)
    dD <- rep(0,numgroups)
    
    kontime <- rep(0,numgroups)
    
    # dynamics without cure
    if(t<times[dt*(prep_fit_idxs[1]-1)+1]){
      kon[3] <- (m_before3 + kof[3]*Sp[3]) / (S[3])
      kon[4] <- (m_before4 + kof[4]*Sp[4]) / (S[4])
      dS <- beta * NN * q + kof * Sp - (mu + J + kon) * S
      dSp <- kon * S - (mu + Omega * J + kof) * Sp
      dI[,1] <- J * S - (mu + rho[1] + tau[1]) * I[,1] + beta * NN * qadj *(imp * prob_stages[1])
      for (k in 2:numstages) {
        dI[,k] <- rho[k-1]*I[,k-1] - (mu+rho[k]+tau[k])*I[,k] + beta * NN * qadj *(imp * prob_stages[k])
      }
      dIp <- Omega * J * Sp - (mu + tau_p) * Ip
      dA <- eta * D - (mu + gamma) * A + qadj * m_imp * t
      dD <- tau_p * Ip - (eta + mu + gamma_diag) * D
      for (k in 1:numstages) {
        dD <- dD + I[,k]*tau[k]
      }
    }else if(t<t_end_prep & t>=times[dt*(prep_fit_idxs[1]-1)+1]){
      kon[3] <- (m_after3 + kof[3]*Sp[3]) / (S[3])
      kon[4] <- (m_after4 + kof[4]*Sp[4]) / (S[4])
      dS <- beta * NN * q + kof * Sp - (mu + J + kon) * S
      dSp <- kon * S - (mu + Omega * J + kof) * Sp
      dI[,1] <- J * S - (mu + rho[1] + tau[1]) * I[,1] + beta * NN * qadj *(imp * prob_stages[1])
      for (k in 2:numstages) {
        dI[,k] <- rho[k-1]*I[,k-1] - (mu+rho[k]+tau[k])*I[,k] + beta * NN * qadj *(imp * prob_stages[k])
      }
      dIp <- Omega * J * Sp - (mu + tau_p) * Ip
      dA <- eta * D - (mu + gamma) * A + qadj * m_imp * t 
      dD <- tau_p * Ip - (eta + mu + gamma_diag) * D
      for (k in 1:numstages) {
        dD <- dD + I[,k]*tau[k]
      }
    }else if(t>=t_end_prep){
      kon[3] <- (m_after3 + kof[3]*Sp[3]) / (S[3])
      kon[4] <- (m_after4 + kof[4]*Sp[4]) / (S[4])
      kontime[3] <- kon[3] / (1 + exp(-1*(-sum(Sp)+thr_prep)) )
      kontime[4] <- kon[4] / (1 + exp(-1*(-sum(Sp)+thr_prep)) )
      dS <- beta * NN * q + kof * Sp - (mu + J + kontime) * S
      dSp <- kontime * S - (mu + Omega * J + kof) * Sp
      dI[,1] <- J * S - (mu + rho[1] + tau[1]) * I[,1] + beta * NN * qadj *(imp * prob_stages[1])
      for (k in 2:numstages) {
        dI[,k] <- rho[k-1]*I[,k-1] - (mu+rho[k]+tau[k])*I[,k] + beta * NN * qadj *(imp * prob_stages[k])
      }
      dIp <- Omega * J * Sp - (mu + tau_p) * Ip
      dA <- eta * D - (mu + gamma) * A + qadj * m_imp * t
      dD <- tau_p * Ip - (eta + mu + gamma_diag) * D
      for (k in 1:numstages) {
        dD <- dD + I[,k]*tau[k]
      }
    }
    res <- c(dS, dSp, dI, dIp, dA, dD)
    list(res)
  })
}

# Define sex risk groups
numstages <- 4
stages <- 1:numstages

# Define population fractions based on sex risk groups
population_fractions <- c(0.451, 0.353, 0.125, 0.071)

# Define HIV stages
numgroups <- length(population_fractions)
groups <- 1:numgroups

# Importing epidemiological data from 2015 to 2022
dat_epi <- read_excel("DATA/EPIDEMIOLOGICAL/HIV_epidemiological_data.xlsx", sheet = "SHM", col_names = F)
dat_epi <- as.data.frame(dat_epi)
dat_epi <- t(dat_epi)
colnames(dat_epi) <- dat_epi[1,]
dat_epi <- dat_epi[-1,]
dat_epi <- as.data.frame(dat_epi)
idxna <- which(is.na(dat_epi[1,]))[1]
dat_epi <- dat_epi[,-c(idxna:ncol(dat_epi))]
dat_epi <- separate(dat_epi, `estimated number of currently infected`, into = c("mean_prev", "lo_prev"), sep = " \\(")
dat_epi <- separate(dat_epi, lo_prev, into = c("lo_prev", "hi_prev"), sep = " - |\\) ")
dat_epi$hi_prev <- gsub(")", "", dat_epi$hi_prev)
dat_epi <- separate(dat_epi, `estimated number of newly infected`, into = c("mean_inc", "lo_inc"), sep = " \\(")
dat_epi <- separate(dat_epi, lo_inc, into = c("lo_inc", "hi_inc"), sep = " - |\\) ")
dat_epi$hi_inc <- gsub(")", "", dat_epi$hi_inc)
dat_epi <- separate(dat_epi, `estimated number of undiagnosed`, into = c("mean_undiag", "lo_undiag"), sep = " \\(")
dat_epi <- separate(dat_epi, lo_undiag, into = c("lo_undiag", "hi_undiag"), sep = "-|\\) ")
dat_epi$hi_undiag <- gsub(")", "", dat_epi$hi_undiag)
dat_epi <- mutate_all(dat_epi, ~gsub(",", ".", .))
dat_epi <- dat_epi %>%
  mutate(across(everything(), ~ as.numeric(.) %>% round(3)))
dat_epi <- dat_epi[nrow(dat_epi):1, ]

idx_2022 <- which(dat_epi$Year == 2022)

# Diagnosed cases
diag_currently <- dat_epi$`number of diagnosed`[-idx_2022]
diag_2022 <- dat_epi$`number of diagnosed`[idx_2022]

# Treated cases
treat_currently <- dat_epi$`number of treated`[-idx_2022]
treat_2022 <- dat_epi$`number of treated`[idx_2022]

# Estimated undiagnosed cases
currently_undiag <- dat_epi$mean_undiag[-idx_2022]
undiag_2022 <- dat_epi$mean_undiag[idx_2022]
cilow_undiag <- dat_epi$lo_undiag[-idx_2022]
cihigh_undiag <- dat_epi$hi_undiag[-idx_2022]
lo_undiag_2022 <- dat_epi$lo_undiag[idx_2022]
hi_undiag_2022 <- dat_epi$hi_undiag[idx_2022]

# Estimated incidence
newly_infected_est <- dat_epi$mean_inc[-idx_2022]
low_newly_infected_est <- dat_epi$lo_inc[-idx_2022]
high_newly_infected_est <- dat_epi$hi_inc[-idx_2022]
inc_2022 <- dat_epi$mean_inc[idx_2022]
lo_inc_2022 <- dat_epi$lo_inc[idx_2022]
hi_inc_2022 <- dat_epi$hi_inc[idx_2022]

# Estimated prevalence
currently_infected_est <- dat_epi$mean_prev[-idx_2022]
low_currently_infected_est <- dat_epi$lo_prev[-idx_2022]
high_currently_infected_est <- dat_epi$hi_prev[-idx_2022]
prev_2022 <- dat_epi$mean_prev[idx_2022]
lo_prev_2022 <- dat_epi$lo_prev[idx_2022]
hi_prev_2022 <- dat_epi$hi_prev[idx_2022]

# Newly diagnosed cases
newinc <- dat_epi$`newly diagnosed`[-idx_2022]
newinc_2022 <- dat_epi$`newly diagnosed`[idx_2022]

# Newly imported cases already on on ART
newimp <- dat_epi$`newly imported already on treatment`[-idx_2022]
newimp_2022 <- dat_epi$`newly imported already on treatment`[idx_2022]

# Initial number of cases in the compartment D
init_diag <- diag_currently - treat_currently

# Define contact rates data of scenario 0
contact_rates <- read_excel("DATA/SEXUAL_BEHAVIOR/meanpartnerchangerates_total_0.xlsx")
contact_rates <- as.numeric(names(contact_rates))
contact_rates <- unlist(contact_rates)[1:(length(contact_rates)-1)]
contact_rates <- contact_rates * 1.5

# Define fixed parameters (see Table 1)
effectiveness_prep <- 0.14
flux_off_prep <- c(0,0,1/5,1/5)
birth_rate <- 1 / 45
mortality_rate <- 1 / 45
progression <- c(1/0.142, 1/8.439, 1/1.184, 1/1.316)
recovery_rate <- 1 / 61
infectivity <- c(0.62, 0.12, 0.642, 0)
transmissibility_prep <- 0.62 / 2
delay_treat <- 8
diag_rate <- rep(0,numstages)
diag_rate[numstages-1] <- 12
diag_rate[numstages] <- 12
diag_rate_prep <- rep(4,numgroups)
initial_population <- 210000
duration_stages <- 1 / progression
recovery_rate_diag <- 1 / sum(duration_stages)
prob_stages <- duration_stages / sum(duration_stages)
transmissibility_diagnosed <- weighted.mean(infectivity, prob_stages)

# Initialize free parameters
population_fractions_adjusted <- population_fractions
X0 <- 0
flux_on_prep <- c(0,0,0,0)
transmissibility_art <- 0.01
infection_prob <- 0.25
prob_import <- 0
mixing <- 0.5

# Define simulation times
nyears_calibration <- nrow(dat_epi) - 1
years_calibration <- as.integer(dat_epi$Year[-idx_2022])
first_year_fit <- 2017
last_year_fit <- 2021
first_year_fit_idx <- which(years_calibration == first_year_fit)
last_year_fit_idx <- which(years_calibration==last_year_fit)
years_fit <- years_calibration[-c(1:(first_year_fit_idx-1))]
years_fit_idxs <- seq(years_fit[1]-years_calibration[1]+1,years_fit[length(years_fit)]-years_calibration[1]+1)
y_end <- 2024
first_year <- years_calibration[1]
years_sim <- first_year:y_end
dt <- 12
nyears <- y_end - first_year + 1
times <- seq(0,nyears-1,by=1/dt)
nt <- length(times)

#Indexes of fit years: from 2017 to 2021
yf_idxs <- seq(dt*(first_year_fit_idx-1)+1,dt*last_year_fit_idx,by=dt)

# Importing prep data from June 2019 to April 2022
dat_prep <- as.data.frame(read_excel("DATA/EPIDEMIOLOGICAL/HIV_epidemiological_data.xlsx", sheet = "PrEP"))
prep <- c(rep(-1,dat_prep$Month[1]-1),dat_prep$`PrEP users`,rep(-1,12-dat_prep$Month[length(dat_prep$Month)]))
prop3 <- population_fractions[3] / (population_fractions[3] + population_fractions[4])
prop4 <- population_fractions[4] / (population_fractions[3] + population_fractions[4])

# Linear fit of the PrEP users by risk group (3 and 4)
years_prep <- unique(dat_prep$Year)[-length(unique(dat_prep$Year))]
prep_before_y <- c(0,prep[12])
prep_after_y <- c(prep[12],prep[12*length(years_prep)-1])
x_before <- 1:which(years_sim == years_prep[1])
x_after <- 1:(length(years_prep))
prep_after_y3 <- prep_after_y*prop3
prep_before_y3 <- prep_before_y*prop3
prep_after_y4 <- prep_after_y*prop4
prep_before_y4 <- prep_before_y*prop4
m_before3 <- diff(prep_before_y3) / (length(x_before)-1)
m_after3 <- diff(prep_after_y3) / (length(x_after)-1)
m_before4 <- diff(prep_before_y4) / (length(x_before)-1)
m_after4 <- diff(prep_after_y4) / (length(x_after)-1)
prep_fit_idxs <- seq(years_prep[1]-years_calibration[1]+1,years_prep[length(years_prep)]-years_calibration[1]+1)
prep_y <- prep[seq(12,length(years_prep)*12,by=12)]
prep_y <- c(rep(-1,years_prep[1]-years_calibration[1]), prep_y)
prep_y <- c(prep_y, rep(-1,-years_prep[length(years_prep)]+years_calibration[length(years_calibration)]))
t_end_prep <- idx_2022 - 1

# Max capacity of PrEP users
thr_prep <- 10000

# Linear fit of the number of newly imported cases already on ART
totimpdata <- c(newimp,newimp_2022)
totimpdata_fit <- totimpdata - totimpdata[1]
x <- seq_along(totimpdata)
excl_idxs <- c(6,7) # Exclude COVID-19 years
lm_fit_imp <- lm(totimpdata_fit[-excl_idxs] ~ 0 +  x[-excl_idxs])
m_imp <- 0

# Importing diagnostic delay distributions
load("DATA/EPIDEMIOLOGICAL/taus12.RData")
taus12 <- outres
rm(outres)

# Fix the already fitted diagnostic delays
avgtau1 <- mean(taus12$tau1)
avgtau2 <- mean(taus12$tau2)
diag_rate[1] <- avgtau1
diag_rate[2] <- avgtau2

# Build the list of parameters to be passed to the model function
parms  <- c(beta=birth_rate, mu=mortality_rate, gamma = recovery_rate, gamma_diag = recovery_rate_diag, rho=progression, 
            h=infectivity, eps=transmissibility_art, eps_p=transmissibility_prep, eps_d=transmissibility_diagnosed, tau=diag_rate,
            omega=mixing, Omega=effectiveness_prep, NN=initial_population, tau_p=diag_rate_prep, eta=delay_treat,
            lambda=infection_prob, c=contact_rates, q=population_fractions, kon=flux_on_prep,
            kof=flux_off_prep, ur=X0, imp=prob_import, qadj=population_fractions_adjusted
)

# Indexes of free parameters
parms_fit_idxs <- c()
parms_fit_idxs <- c(parms_fit_idxs, which(names(parms) == c("lambda")))
parms_fit_idxs <- c(parms_fit_idxs, which(names(parms) == c("omega")))
parms_fit_idxs <- c(parms_fit_idxs, which(names(parms) == c("eps")))
parms_fit_idxs <- c(parms_fit_idxs, which(names(parms) == c("ur")))
parms_fit_idxs <- c(parms_fit_idxs, which(names(parms) == c("imp")))
parms_fit_idxs <- c(parms_fit_idxs, which(names(parms) == c("qadj1")))
parms_fit_idxs <- c(parms_fit_idxs, which(names(parms) == c("qadj2")))
parms_fit_idxs <- c(parms_fit_idxs, which(names(parms) == c("qadj3")))
parms_fit_idxs <- c(parms_fit_idxs, which(names(parms) == c("qadj4")))
n_pf <- length(parms_fit_idxs)

# Priors of the free parameters
parms_fit_ranges <- matrix(c(0.01, 0.1, #lambda
                             0, 1, #omega
                             0.00001, 0.075, #eps
                             0, currently_infected_est[1] * 30 / 100, #ur
                             0, 0, #imp
                             0.001, 1, #qadj1
                             0.001, 1, #qadj2
                             0.001, 1, #qadj3
                             0.001, 1), #qadj4
                           nrow = n_pf, byrow = T)

# Generate and rescale samples of free parameters based on their priors
num_lhs_samples <- 5000000
lhs_samples <- randomLHS(num_lhs_samples, n_pf)
scaled_lhs_samples <- t(apply(lhs_samples, 1, function(x) {
  parms_fit_ranges[,1] + x * (parms_fit_ranges[,2] - parms_fit_ranges[,1])
}))

# Select only the sets for which Q1 < Q2 < Q3 < Q4
scaled_lhs_samples <- scaled_lhs_samples[which(scaled_lhs_samples[,6] < scaled_lhs_samples[,7] & scaled_lhs_samples[,7] < scaled_lhs_samples[,8] & scaled_lhs_samples[,8] < scaled_lhs_samples[,9]),]
end_abc <- nrow(scaled_lhs_samples)

# Length of accepted distributions
len_post <- 100

# Initialize variables and outputs
S_out <- matrix(0,nt,numgroups)
Sp_out <- matrix(0,nt,numgroups)
Ip_out <- matrix(0,nt,numgroups)
D_out <- matrix(0,nt,numgroups)
A_out <- matrix(0,nt,numgroups)
I_out <- array(0,dim=c(nt,numgroups,numstages))
Ng_out <- matrix(0,nt,numgroups)
prevalence <- rep(0,nt)
impdiag <- rep(0,nt)
impinf <- rep(0,nt)
incidence <- rep(0,nt)
prep_users <- rep(0,nt)
diagnosed <- rep(0,nt)
undiagnosed <- rep(0,nt)
art_users <- rep(0,nt)
newtreat <- rep(0,nt)
newdiag <- rep(0,nt)

accp_list <- matrix(NA, nrow = len_post, ncol = n_pf)

counter_diag <- matrix(NA, nrow = len_post, ncol = nt)
counter_undiag <- matrix(NA, nrow = len_post, ncol = nt)
counter_art <- matrix(NA, nrow = len_post, ncol = nt)
counter_prevalence <- matrix(NA, nrow = len_post, ncol = nt)
counter_impdiag <- matrix(NA, nrow = len_post, ncol = nt)
counter_impundiag <- matrix(NA, nrow = len_post, ncol = nt)
counter_incidence <- matrix(NA, nrow = len_post, ncol = nt)
counter_prep <- matrix(NA, nrow = len_post, ncol = nt)
counter_newtreat <- matrix(NA, nrow = len_post, ncol = nt)
counter_newdiag <- matrix(NA, nrow = len_post, ncol = nt)
counter_pop <- matrix(NA, nrow = len_post, ncol = nt)
counter_fit <- rep(999999, len_post)

rescale_fit <- length(years_fit_idxs)
numaccp <- 1
threshold_fit <- 0.15

#ABC
i <- 1
while(numaccp <= len_post){
  print(c(i,numaccp))
  
  # Propose parameters set
  parms[parms_fit_idxs] <- scaled_lhs_samples[i,]
  
  # Normalize Q
  sum_qadj <- 0
  for(l in 1:numgroups){
    sum_qadj <- sum_qadj + parms[paste0("qadj",l)]
  }
  for(l in 1:numgroups){
    parms[paste0("qadj",l)] <- parms[paste0("qadj",l)] / sum_qadj
  }
  
  # Initialize the dynamics
  AA <- treat_currently[1]
  DD <- init_diag[1]
  II <- runif(1, cilow_undiag[1], cihigh_undiag[1]) + parms["ur"]
  S_0 <- rep(0,numgroups)
  Sp_0 <- rep(0,numgroups)
  I_0 <- matrix(0,numgroups,numstages)
  Ip_0 <- rep(0,numgroups)
  A_0 <- rep(0,numgroups)
  D_0 <- rep(0,numgroups)
  II_groups <- rep(0,numgroups)
  for(l in 1:numgroups){
    A_0[l] <- AA * parms[paste0("qadj",l)]
    D_0[l] <- DD * parms[paste0("qadj",l)]
    II_groups[l] <- II * parms[paste0("qadj",l)]
  }
  I_0 <- matrix(0, numgroups, numstages)
  for (k in 1:numstages) {
    I_0[,k] <- II_groups * prob_stages[k]
  }
  initial_Ng <- initial_population * population_fractions
  for(l in 1:numgroups){
    S_0[l] <- initial_Ng[l] - A_0[l] - D_0[l] - sum(I_0[l,])
  }
  
  # Call the model function and store the output
  xstart <- c(S=S_0, Sp=Sp_0, I=I_0, Ip=Ip_0, A=A_0, D=D_0)
  out <-  as.data.frame(lsoda(xstart, times, model, parms))
  idx <- 1
  for(l in 1:numgroups){
    S_out[,l] <- out[,paste0("S",l)]
    Sp_out[,l] <- out[,paste0("Sp",l)]
    Ip_out[,l] <- out[,paste0("Ip",l)]
    A_out[,l] <- out[,paste0("A",l)]
    D_out[,l] <- out[,paste0("D",l)]
    for(k in 1:numstages){
      I_out[,l,k] <- out[,paste0("I",idx)]
      idx <- idx + 1
    }
  }
  for (t in 1:nt) {
    I_out[t,,] <- t(I_out[t,,])
  }
  for (t in 1:nt) {
    sumIl <- rowSums(I_out[t,,])
    Ng_out[t,] <- S_out[t,] + Sp_out[t,] + Ip_out[t,] + A_out[t,] + D_out[t,] + sumIl
  }
  N_out <- rowSums(Ng_out)
  for (t in 1:nt) {
    M <- matrix(0,numgroups,numgroups)
    denM <- sum(contact_rates*Ng_out[t,])
    for(l in 1:numgroups){
      for (ll in 1:numgroups) {
        M[l,ll] <- parms["omega"] * contact_rates[ll]*Ng_out[t,ll] / denM
        if(l==ll){
          M[l,l] <- M[l,l] + (1 - parms["omega"])
        }
      }
    }
    J <- rep(0,numgroups)
    for(l in 1:numgroups){
      for (ll in 1:numgroups) {
        J[l] <- J[l] + M[l, ll] * (transmissibility_prep * Ip_out[t,ll] / Ng_out[t,ll] + (parms["eps"] * A_out[t,ll] + transmissibility_diagnosed * D_out[t,ll]) / Ng_out[t,ll])
        for (k in 1:numstages) {
          J[l] <- J[l] + M[l, ll] * infectivity[k] * (I_out[t,ll,k]) / Ng_out[t,ll]
        }
      }
      J[l] <- J[l] * (parms["lambda"] * contact_rates[l])
    }
    
    # Calculate the quantities of interest
    ii <- parms["beta"] * initial_population * parms["imp"]
    it <- m_imp * (t-1) / dt + totimpdata[1]
    incidence[t] <- sum(J*S_out[t,]) + effectiveness_prep * sum(J*Sp_out[t,])
    art_users[t] <- sum(A_out[t,])
    diagnosed[t] <- sum(D_out[t,]) + sum(A_out[t,])
    undiagnosed[t] <- sum(Ip_out[t,]) + sum(I_out[t,,])
    prep_users[t] <- sum(Sp_out[t,]) + sum(Ip_out[t,])
    prevalence[t] <- sum(Ip_out[t,]) + sum(I_out[t,,]) + sum(A_out[t,]) + sum(D_out[t,])
    impdiag[t] <- it
    impinf[t] <- ii
    newdiag[t] <- sum(parms["tau1"]*I_out[t,,1]) + sum(parms["tau2"]*I_out[t,,2]) +
      sum(parms["tau3"]*I_out[t,,3]) + sum(parms["tau4"]*I_out[t,,4])
    newtreat[t] <- sum(parms["eta"]*D_out[t,])
  }
  
  # MAPE of newly diagnoses
  fit_newdiag_diff <- sum(abs(newdiag[yf_idxs] - newinc[years_fit_idxs]) / newinc[years_fit_idxs]) / rescale_fit
  
  # MAPE of undiagnosed cases
  fit_undiag_diff <- sum(abs(undiagnosed[yf_idxs] - currently_undiag[years_fit_idxs]) / currently_undiag[years_fit_idxs]) / rescale_fit
  
  # Mean MAPE
  fit_diff <- (fit_newdiag_diff + fit_undiag_diff) / 2
  
  print(fit_diff)
  
  if(fit_diff < threshold_fit){
    accp_list[numaccp, ] <- parms[parms_fit_idxs]
    counter_diag[numaccp, ] <- diagnosed
    counter_undiag[numaccp, ] <- undiagnosed
    counter_art[numaccp, ] <- art_users
    counter_prevalence[numaccp, ] <- prevalence
    counter_impdiag[numaccp, ] <- impdiag
    counter_impundiag[numaccp, ] <- impinf
    counter_incidence[numaccp, ] <- incidence
    counter_prep[numaccp, ] <- prep_users
    counter_newtreat[numaccp, ] <- newtreat
    counter_newdiag[numaccp, ] <- newdiag
    counter_pop[numaccp, ] <- rowSums(Ng_out)
    counter_fit[numaccp] <- fit_diff
    numaccp <- numaccp + 1
  }
  
  i <- i + 1
}
sampled_sets <- i

accp_list <- as.data.frame(accp_list)
if(nrow(accp_list)>0){
  names(accp_list) <- names(parms[parms_fit_idxs])
}

save.image(paste0("RESULTS/fit_threshold",round(100*threshold_fit,0),"_noimp.RData"))

