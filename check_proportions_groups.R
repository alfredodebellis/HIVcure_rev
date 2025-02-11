rm(list = ls())

library(deSolve)
library(lhs)
library(RColorBrewer)

load("RESULTS/fit_threshold15.RData")
rm(model)

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

# Number of accepted free parameters sets to be sampled
len_post <- 100

prop_inf_groups <- matrix(0,numgroups,nt)
counter_prop_inf1 <- c()
counter_prop_inf2 <- c()
counter_prop_inf3 <- c()
counter_prop_inf4 <- c()

for(i in 1:len_post){
  print(i)
  for(j in 1:n_pf){
    parms[parms_fit_idxs][j] <- accp_list[i,j]
  }
  
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
  initial_Ng <- parms["NN"] * population_fractions
  for(l in 1:numgroups){
    S_0[l] <- initial_Ng[l] - A_0[l] - D_0[l] - sum(I_0[l,])
  }
  
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
    
    for(l in 1:numgroups){
      prop_inf_groups[l,t] <- sum(I_out[t,l,]) / sum((I_out[t,,]))
    }
  }
  counter_prop_inf1 <- rbind(counter_prop_inf1, prop_inf_groups[1,])
  counter_prop_inf2 <- rbind(counter_prop_inf2, prop_inf_groups[2,])
  counter_prop_inf3 <- rbind(counter_prop_inf3, prop_inf_groups[3,])
  counter_prop_inf4 <- rbind(counter_prop_inf4, prop_inf_groups[4,])
}

# Calculate means and credible intervals for the 4 groups
mean_prop_inf1 <- colMeans(counter_prop_inf1)
ci_prop_inf1 <- matrix(0,2,length(mean_prop_inf1))
for (i in 1:length(mean_prop_inf1)) {
  ci_prop_inf1[1,i] <- quantile(counter_prop_inf1[,i], 0.975)
  ci_prop_inf1[2,i] <- quantile(counter_prop_inf1[,i], 0.025)
}
mean_prop_inf2 <- colMeans(counter_prop_inf2)
ci_prop_inf2 <- matrix(0,2,length(mean_prop_inf1))
for (i in 1:length(mean_prop_inf1)) {
  ci_prop_inf2[1,i] <- quantile(counter_prop_inf2[,i], 0.975)
  ci_prop_inf2[2,i] <- quantile(counter_prop_inf2[,i], 0.025)
}
mean_prop_inf3 <- colMeans(counter_prop_inf3)
ci_prop_inf3 <- matrix(0,2,length(mean_prop_inf1))
for (i in 1:length(mean_prop_inf1)) {
  ci_prop_inf3[1,i] <- quantile(counter_prop_inf3[,i], 0.975)
  ci_prop_inf3[2,i] <- quantile(counter_prop_inf3[,i], 0.025)
}
mean_prop_inf4 <- colMeans(counter_prop_inf4)
ci_prop_inf4 <- matrix(0,2,length(mean_prop_inf1))
for (i in 1:length(mean_prop_inf1)) {
  ci_prop_inf4[1,i] <- quantile(counter_prop_inf4[,i], 0.975)
  ci_prop_inf4[2,i] <- quantile(counter_prop_inf4[,i], 0.025)
}

load("PLOTS/colors.Rdata")
save.image(paste0("RESULTS/check_proportions_groups.RData"))
