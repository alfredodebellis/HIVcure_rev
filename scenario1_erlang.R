rm(list=ls())

library(deSolve)
library(lhs)
library(readxl)

load("RESULTS/fit_threshold15.RData")

# Initialize cure parameters and append them to the parameters list of the fit

# Decide whether simulate with baseline hc (infectivity of rebounds) or sensitivity analysis
hf <- 0 # baseline
# hf <- 1 # sensitivity analysis
if(hf == 0){
  # First rebounds stage is chronic
  hcfile <- "MAIN"
  hc <- c(infectivity[2],infectivity[2],infectivity[3],infectivity[4])
}else{
  # First rebounds stage is acute
  hcfile <- "SA"
  hc <- c(infectivity[1],infectivity[2],infectivity[3],infectivity[4])
}
alpha <- 1
parms <- c(parms, alpha)
names(parms)[length(parms)] <- "alpha"

# Decide the value of the cure uptake
ar <- 0.9
# ar <- 0.5
# ar <- 0.1
parms["alpha"] <- -log(1-ar)

c0 <- 0.01
parms <- c(parms, c0)
names(parms)[length(parms)] <- "c0"
r <- 3
parms <- c(parms, r)
names(parms)[length(parms)] <- "r"
eff <- 1
parms <- c(parms, eff)
names(parms)[length(parms)] <- "eff"
# Decide the value of the cure efficacy based on TPP
efficacies <- c(0.2,0.9)
# cure_efficacy_option <- "minimum"
cure_efficacy_option <- "optimal"
if(cure_efficacy_option == "optimal"){
  parms["eff"] <- efficacies[2]
}else{
  parms["eff"] <- efficacies[1]
}
phi <- 1
parms <- c(parms, phi)
names(parms)[length(parms)] <- "phi"
parms <- c(parms, hc)
names(parms)[(length(parms)-4+1):length(parms)] <- sapply(1:numstages, function(k) paste0("hc", k))
tau_c <- rep(1,numstages)
parms <- c(parms, tau_c)
names(parms)[(length(parms)-4+1):length(parms)] <- sapply(1:numstages, function(k) paste0("tau_c", k))

rm(model)

# Define the new times for the simulation
y_end <- 2036
years_sim_extend <- c(years_sim,(years_sim[length(years_sim)]+1):y_end)
nyears <- y_end - first_year + 1
times <- seq(0,nyears-1,by=1/dt)
nt <- length(times)

# Define the cure start year
y_cure <- 2026
t_cure <- which(years_sim_extend==y_cure) - 1

#Erlang (how many additional compartments)
E <- 3

# Model function for scenario 1
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
    omega <- parms["omega"]
    Omega <- parms["Omega"]
    imp <- parms["imp"]
    phi <- parms["phi"]
    alpha <- parms["alpha"]
    NN <- parms["NN"]
    r <- parms["r"]
    c0 <- parms["c0"]
    eff <- parms["eff"]
    rho <- sapply(1:numstages, function(k) parms[paste0("rho", k)])
    h <- sapply(1:numstages, function(k) parms[paste0("h", k)])
    hc <- sapply(1:numstages, function(k) parms[paste0("hc", k)])
    tau_p <- sapply(1:numgroups, function(l) parms[paste0("tau_p", l)])
    tau <- sapply(1:numstages, function(k) parms[paste0("tau", k)])
    q <- sapply(1:numgroups, function(l) parms[paste0("q", l)])
    c <- sapply(1:numgroups, function(l) parms[paste0("c", l)])
    kon <- sapply(1:numgroups, function(l) parms[paste0("kon", l)])
    kof <- sapply(1:numgroups, function(l) parms[paste0("kof", l)])
    tau_c <- sapply(1:numgroups, function(l) parms[paste0("tau_c", l)])
    qadj <- sapply(1:numgroups, function(l) parms[paste0("qadj", l)])
    
    #variables
    S <- sapply(1:numgroups, function(l) x[paste0("S", l)])
    Sp <- sapply(1:numgroups, function(l) x[paste0("Sp", l)])
    Ip <- sapply(1:numgroups, function(l) x[paste0("Ip", l)])
    I <- matrix(sapply(1:(numgroups * numstages), function(i) x[paste0("I", i)]), nrow = numgroups, ncol = numstages)
    D <- sapply(1:numgroups, function(l) x[paste0("D", l)])
    A <- sapply(1:numgroups, function(l) x[paste0("A", l)])
    C <- sapply(1:numgroups, function(l) x[paste0("C", l)])
    R <- sapply(1:numgroups, function(l) x[paste0("R", l)])
    Ic <- matrix(sapply(1:(numgroups * numstages), function(i) x[paste0("Ic", i)]), nrow = numgroups, ncol = numstages)
    Cprimo <- sapply(1:numgroups, function(l) x[paste0("Cprimo", l)])
    Csecondo <- sapply(1:numgroups, function(l) x[paste0("Csecondo", l)])
    
    Ng <- rowSums(cbind(S, Sp, rowSums(I), Ip, D, A, R, C, rowSums(Ic), Cprimo, Csecondo))
    N <- sum(Ng)
    
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
    
    J <- rep(0,numgroups)
    for(l in 1:numgroups){
      for (ll in 1:numgroups) {
        J[l] <- J[l] + M[l, ll] * (eps_p * Ip[ll] / Ng[ll] + eps * (A[ll] + R[ll]) / Ng[ll] + eps_d * D[ll] / Ng[ll])
        for (k in 1:numstages) {
          J[l] <- J[l] + M[l, ll] * (h[k] *I[ll, k] + hc[k] * Ic[ll, k]) / Ng[ll]
        }
      }
      J[l] <- J[l] * (lambda * c[l])
    }
    
    dS <- rep(0,numgroups)
    dSp <- rep(0,numgroups)
    dI <- matrix(0,numgroups,numstages)
    dIp <- rep(0,numgroups)
    dD <- rep(0,numgroups)
    dA <- rep(0,numgroups)
    dC <- rep(0,numgroups)
    dR <- rep(0,numgroups)
    dIc <- matrix(0,numgroups,numstages)
    dCprimo <- rep(0,numgroups)
    dCsecondo <- rep(0,numgroups)
    
    kontime <- rep(0,numgroups)
    
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
    }else if(t>=t_end_prep & t<t_cure){
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
    }else if(t>=t_cure){
      alphatime <- alpha*c0*exp(r*(t-t_cure)) / (1 + c0*(exp(r*(t-t_cure))-1))
      kon[3] <- (m_after3 + kof[3]*Sp[3]) / (S[3])
      kon[4] <- (m_after4 + kof[4]*Sp[4]) / (S[4])
      kontime[3] <- kon[3] / (1 + exp(-1*(-sum(Sp)+thr_prep)) )
      kontime[4] <- kon[4] / (1 + exp(-1*(-sum(Sp)+thr_prep)) )
      dS <- beta * NN * q + kof * Sp - (mu + J + kontime) * S
      dSp <- kontime * S - (mu + Omega * J + kof) * Sp
      dIp <- Omega * J * Sp - (mu + tau_p) * Ip
      dD <- tau_p * Ip - (eta + mu + gamma_diag) * D
      for (k in 1:numstages) {
        dD <- dD + I[,k]*tau[k] + tau_c[k] * Ic[,k]
      }
      dA <- eta * D - (mu + gamma + alphatime) * A + qadj * m_imp * t
      dC <- - (mu + gamma + E*phi) * C + eff*alphatime * A
      dCprimo <- - (mu + gamma) * Cprimo + E*phi * (C - Cprimo)
      dCsecondo <- - (mu + gamma) * Csecondo + E*phi * (Cprimo - Csecondo)
      dR <- (1.-eff)*alphatime * A - (mu + gamma) * R
      dIc[,1] <- phi * E * Csecondo - (mu + rho[1] + tau_c[1]) * Ic[,1]
      dI[,1] <- J * S - (mu + rho[1] + tau[1]) * I[,1] + beta * NN * qadj *(imp * prob_stages[1])
      for (k in 2:numstages) {
        dIc[,k] <- rho[k-1]*Ic[,k-1] - (mu+rho[k]+tau_c[k])*Ic[,k]
        dI[,k] <- rho[k-1]*I[,k-1] - (mu+rho[k]+tau[k])*I[,k] + beta * NN * qadj *(imp * prob_stages[k])
      }
    }
    res <- c(dS, dSp, dI, dIp, dD, dA, dC, dR, dIc, dCprimo, dCsecondo)
    list(res)
  })
}

# Number of accepted free parameters sets to be sampled
len_post <- 100
# len_post <- nrow(accp_list_bests)

# Number of different monitoring strategy for rebounds
n_monitors <- 3

# Define the rebound rates
failures <- c(0, 0.16666666, 0.5)
n_sa <- length(failures)

# Diagnosis uptake for the follow-up monitoring
taulong <- rep(dt*2,numstages)

# Initialize observables
incidence_art_after_new_time_pool <- array(0,dim=c(n_monitors,n_sa,len_post,nt))
incidence_diag_after_new_time_pool <- array(0,dim=c(n_monitors,n_sa,len_post,nt))
incidence_undiag_after_new_time_pool <- array(0,dim=c(n_monitors,n_sa,len_post,nt))
incidence_after_new_time_pool <- array(0,dim=c(n_monitors,n_sa,len_post,nt))
incidence_after_reb_time_pool <- array(0,dim=c(n_monitors,n_sa,len_post,nt))
coverage_time_pool <- array(0,dim=c(n_monitors,n_sa,len_post,nt))
prevalence_after_new_time_pool <- array(0,dim=c(n_monitors,n_sa,len_post,nt))
prevalence_after_reb_time_pool <- array(0,dim=c(n_monitors,n_sa,len_post,nt))
newly_cured_time_pool <- array(0,dim=c(n_monitors,n_sa,len_post,nt))
art_given_time_pool <- array(0,dim=c(n_monitors,n_sa,len_post,nt))
prevalence <- rep(0,nt)
prevalence_cure <- rep(0,nt)
incidence <- rep(0,nt)
incidence_art <- rep(0,nt)
incidence_undiag <- rep(0,nt)
incidence_diag <- rep(0,nt)
coverage <- rep(0,nt)
incidence_cure <- rep(0,nt)
newly_cured <- rep(0,nt)
art_given <- rep(0,nt)
S_out <- matrix(0,nt,numgroups)
Sp_out <- matrix(0,nt,numgroups)
Ip_out <- matrix(0,nt,numgroups)
D_out <- matrix(0,nt,numgroups)
A_out <- matrix(0,nt,numgroups)
I_out <- array(0,dim=c(nt,numgroups,numstages))
C_out <- matrix(0,nt,numgroups)
R_out <- matrix(0,nt,numgroups)
Ic_out <- array(0,dim=c(nt,numgroups,numstages))
Cprimo_out <- matrix(0,nt,numgroups)
Csecondo_out <- matrix(0,nt,numgroups)
Ng_out <- matrix(0,nt,numgroups)

# Simulate scenario
# for (g in 1:n_monitors) {
for (g in 2:2) {
  for(f in failures){
    parms["phi"] <- f
    print(c(g,f))
    for(i in 1:len_post){
      print(i)
      for(j in 1:n_pf){
        parms[parms_fit_idxs][j] <- accp_list[i,j]
        # parms[parms_fit_idxs][j] <- accp_list_bests[i,j]
      }
      
      # Assign monitoring strategy values
      if(g==1){
        for (k in 1:numstages) {
          parms[paste0("tau_c",k)] <- parms[paste0("tau",k)]
        }
      }else if(g==2){
        parms[paste0("tau_c",1)] <- parms[paste0("tau_p",1)]
        parms[paste0("tau_c",2)] <- parms[paste0("tau_p",2)]
        parms[paste0("tau_c",3)] <- parms[paste0("tau",3)]
        parms[paste0("tau_c",4)] <- parms[paste0("tau",4)]
      }else if(g==3){
        for (k in 1:numstages) {
          parms[paste0("tau_c",k)] <- taulong[k]
        }
      }
      
      AA <- treat_currently[1]
      DD <- init_diag[1]
      II <- runif(1, cilow_undiag[1], cihigh_undiag[1]) + parms["ur"]
      S_0 <- rep(0,numgroups)
      Sp_0 <- rep(0,numgroups)
      I_0 <- matrix(0,numgroups,numstages)
      Ip_0 <- rep(0,numgroups)
      Ic_0 <- matrix(0,numgroups,numstages)
      C_0 <- rep(0,numgroups)
      Cprimo_0 <- rep(0,numgroups)
      Csecondo_0 <- rep(0,numgroups)
      R_0 <- rep(0,numgroups)
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
      
      xstart <- c(S=S_0, Sp=Sp_0, I=I_0, Ip=Ip_0, D=D_0, A=A_0, C=C_0, R=R_0, Ic=Ic_0, Cprimo=Cprimo_0, Csecondo=Csecondo_0)
      out <-  as.data.frame(lsoda(xstart, times, model, parms))
      idx <- 1
      for(l in 1:numgroups){
        S_out[,l] <- out[,paste0("S",l)]
        Sp_out[,l] <- out[,paste0("Sp",l)]
        Ip_out[,l] <- out[,paste0("Ip",l)]
        D_out[,l] <- out[,paste0("D",l)]
        A_out[,l] <- out[,paste0("A",l)]
        C_out[,l] <- out[,paste0("C",l)]
        Cprimo_out[,l] <- out[,paste0("Cprimo",l)]
        Csecondo_out[,l] <- out[,paste0("Csecondo",l)]
        R_out[,l] <- out[,paste0("R",l)]
        for(k in 1:numstages){
          I_out[,l,k] <- out[,paste0("I",idx)]
          Ic_out[,l,k] <- out[,paste0("Ic",idx)]
          idx <- idx + 1
        }
      }
      for (t in 1:nt) {
        I_out[t,,] <- t(I_out[t,,])
        Ic_out[t,,] <- t(Ic_out[t,,])
      }
      for (t in 1:nt) {
        sumIl <- rowSums(I_out[t,,])
        sumIcl <- rowSums(Ic_out[t,,])
        Ng_out[t,] <- S_out[t,] + Sp_out[t,] + Ip_out[t,] + D_out[t,] + A_out[t,] + sumIl + C_out[t,] + R_out[t,] + sumIcl + Cprimo_out[t,] + Csecondo_out[t,]
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
        J_undiag <- rep(0,numgroups)
        J_diag <- rep(0,numgroups)
        J_art <- rep(0,numgroups)
        for(l in 1:numgroups){
          for (ll in 1:numgroups) {
            J_art[l] <- J_art[l] +  M[l, ll] * parms["eps"] * (A_out[t,ll] + R_out[t,ll]) / Ng_out[t,ll]
            J_diag[l] <- J_diag[l] +  M[l, ll] * transmissibility_diagnosed * D_out[t,ll] / Ng_out[t,ll]
            J_undiag[l] <- J_undiag[l] +  M[l, ll] * transmissibility_prep * Ip_out[t,ll] / Ng_out[t,ll]
            J[l] <- J[l] + M[l, ll] * (transmissibility_prep * Ip_out[t,ll] / Ng_out[t,ll] + parms["eps"] * (A_out[t,ll] + R_out[t,ll]) / Ng_out[t,ll] + transmissibility_diagnosed * D_out[t,ll] / Ng_out[t,ll])
            for (k in 1:numstages) {
              J_undiag[l] <- J_undiag[l] +  M[l, ll] * (infectivity[k] * I_out[t,ll,k] + hc[k] * Ic_out[t,ll,k]) / Ng_out[t,ll]
              J[l] <- J[l] + M[l, ll] * (infectivity[k] * I_out[t,ll,k] + hc[k] * Ic_out[t,ll,k]) / Ng_out[t,ll]
            }
          }
          J[l] <- J[l] * (parms["lambda"] * contact_rates[l])
          J_art[l] <- J_art[l] * (parms["lambda"] * contact_rates[l])
          J_diag[l] <- J_diag[l] * (parms["lambda"] * contact_rates[l])
          J_undiag[l] <- J_undiag[l] * (parms["lambda"] * contact_rates[l])
        }
        alphatimeext <- parms["alpha"]*parms["c0"]*exp(parms["r"]*(t-dt*t_cure-1)) / (1 + parms["c0"]*(exp(parms["r"]*(t-dt*t_cure-1))-1))
        incidence[t] <- sum(J*S_out[t,]) + effectiveness_prep * sum(J*Sp_out[t,])
        incidence_art[t] <- sum(J_art*S_out[t,]) + effectiveness_prep * sum(J_art*Sp_out[t,])
        incidence_undiag[t] <- sum(J_undiag*S_out[t,]) + effectiveness_prep * sum(J_undiag*Sp_out[t,])
        incidence_diag[t] <- sum(J_diag*S_out[t,]) + effectiveness_prep * sum(J_diag*Sp_out[t,])
        prevalence[t] <- sum(Ip_out[t,]) + sum(I_out[t,,]) + sum(A_out[t,]) + sum(D_out[t,]) + sum(R_out[t,])
        prevalence_cure[t] <- sum(Ic_out[t,,]) + sum(C_out[t,]) + sum(Cprimo_out[t,]) + sum(Csecondo_out[t,])
        incidence_cure[t] <- parms["phi"] * E * (sum(C_out[t,]) + sum(Cprimo_out[t,]) + sum(Csecondo_out[t,]))
        coverage[t] <- (sum(C_out[t,]) + sum(Cprimo_out[t,]) + sum(Csecondo_out[t,])) / ( sum(A_out[t,]) + sum(R_out[t,]) + sum(C_out[t,]) + sum(Cprimo_out[t,]) + sum(Csecondo_out[t,]))
        newly_cured[t] <- sum(parms["eff"] * alphatimeext * A_out[t,])
        art_given[t] <- sum(A_out[t,] + R_out[t,])
      }
      
      coverage_time_pool[g,which(failures==f),i,] <- coverage
      incidence_after_new_time_pool[g,which(failures==f),i,] <- incidence
      incidence_art_after_new_time_pool[g,which(failures==f),i,] <- incidence_art
      incidence_diag_after_new_time_pool[g,which(failures==f),i,] <- incidence_diag
      incidence_undiag_after_new_time_pool[g,which(failures==f),i,] <- incidence_undiag
      incidence_after_reb_time_pool[g,which(failures==f),i,] <- incidence_cure
      prevalence_after_new_time_pool[g,which(failures==f),i,] <- prevalence
      prevalence_after_reb_time_pool[g,which(failures==f),i,] <- prevalence_cure
      newly_cured_time_pool[g,which(failures==f),i,] <- newly_cured
      art_given_time_pool[g,which(failures==f),i,] <- art_given
    }
  }
}

# Calculate means and CrIs of observables
mean_inc_new <- array(apply(incidence_after_new_time_pool, c(1, 2, 4), mean), dim=c(n_monitors, n_sa, nt))
mean_inc_reb <- array(apply(incidence_after_reb_time_pool, c(1, 2, 4), mean), dim=c(n_monitors, n_sa, nt))
q1_inc_new <- array(apply(incidence_after_new_time_pool, c(1, 2, 4), function(x) quantile(x, 0.025)), dim=c(n_monitors, n_sa, nt))
q2_inc_new <- array(apply(incidence_after_new_time_pool, c(1, 2, 4), function(x) quantile(x, 0.975)), dim=c(n_monitors, n_sa, nt))
q1_inc_reb <- array(apply(incidence_after_reb_time_pool, c(1, 2, 4), function(x) quantile(x, 0.025)), dim=c(n_monitors, n_sa, nt))
q2_inc_reb <- array(apply(incidence_after_reb_time_pool, c(1, 2, 4), function(x) quantile(x, 0.975)), dim=c(n_monitors, n_sa, nt))
mean_prev_new <- array(apply(prevalence_after_new_time_pool, c(1, 2, 4), mean), dim=c(n_monitors, n_sa, nt))
mean_prev_reb <- array(apply(prevalence_after_reb_time_pool, c(1, 2, 4), mean), dim=c(n_monitors, n_sa, nt))
q1_prev_new <- array(apply(prevalence_after_new_time_pool, c(1, 2, 4), function(x) quantile(x, 0.025)), dim=c(n_monitors, n_sa, nt))
q2_prev_new <- array(apply(prevalence_after_new_time_pool, c(1, 2, 4), function(x) quantile(x, 0.975)), dim=c(n_monitors, n_sa, nt))
q1_prev_reb <- array(apply(prevalence_after_reb_time_pool, c(1, 2, 4), function(x) quantile(x, 0.025)), dim=c(n_monitors, n_sa, nt))
q2_prev_reb <- array(apply(prevalence_after_reb_time_pool, c(1, 2, 4), function(x) quantile(x, 0.975)), dim=c(n_monitors, n_sa, nt))
mean_coverage <- array(apply(coverage_time_pool, c(1, 2, 4), mean), dim=c(n_monitors, n_sa, nt))
q1_coverage <- array(apply(coverage_time_pool, c(1, 2, 4), function(x) quantile(x, 0.025)), dim=c(n_monitors, n_sa, nt))
q2_coverage <- array(apply(coverage_time_pool, c(1, 2, 4), function(x) quantile(x, 0.975)), dim=c(n_monitors, n_sa, nt))

mean_inc_art_new <- array(apply(incidence_art_after_new_time_pool, c(1, 2, 4), mean), dim=c(n_monitors, n_sa, nt))
mean_inc_diag_new <- array(apply(incidence_diag_after_new_time_pool, c(1, 2, 4), mean), dim=c(n_monitors, n_sa, nt))
mean_inc_undiag_new <- array(apply(incidence_undiag_after_new_time_pool, c(1, 2, 4), mean), dim=c(n_monitors, n_sa, nt))

mean_newly_cured <- array(apply(newly_cured_time_pool, c(1, 2, 4), mean), dim=c(n_monitors, n_sa, nt))
mean_art_given <- array(apply(art_given_time_pool, c(1, 2, 4), mean), dim=c(n_monitors, n_sa, nt))

save.image(paste0("RESULTS/SCENARIOS/SC1_erlangE3_eff", round(parms["eff"]*100,0),"_cu", round(ar*100,0), "_hc", hcfile,"_curestart",y_cure,".RData"))
