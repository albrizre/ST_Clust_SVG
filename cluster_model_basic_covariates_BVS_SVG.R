cluster_model_basic_covariates_BVS_SVG <- nimbleCode({
  
  for (j in 1:N_covariates){
    for (k in 1:(N_clusters+1)){
      I[j,k] ~ dbern(0.5)
      # pi[j,k] ~ dbeta(1,9) # Sensitivity analysis on the prior of I_jk (suppl. material)
      # I[j,k] ~ dbern(pi[j,k])
      beta[j,k] ~ dnorm(0, sd = 2)
      beta_I[j,k] <- I[j,k]*beta[j,k]
    }
  }
  
  for (k in 1:(N_clusters+1)){
    sd_alpha[k] ~ dunif(0,5)
  }
  
  for(i in 1:N_genes) {
    y[i, 1:N_comp] ~ dmnorm(mus[1:N_comp, z[i]], wi_tau[1:N_comp, 1:N_comp, i])
    wi_tau[1:N_comp, 1:N_comp, i] <- w[i]*tau[1:N_comp, 1:N_comp]
    w[i] ~ dgamma(2, 2)
    z[i] ~ dcat(omega[i, 1:(N_clusters+1)])
    for (k in 1:(N_clusters+1)){
      alpha[i, k] ~ dnorm(0, sd = sd_alpha[k])
      omega[i, k] <- Phi[i, k]/sum(Phi[i, 1:(N_clusters+1)])
      log(Phi[i, k]) <- alpha[i,k] + inprod(beta_I[1:N_covariates, k],GO_covariates[i, 1:N_covariates])
    }
  }
  for (k in 1:N_clusters){
    mus[1:N_comp, k] ~ dmnorm(mu_0[1:N_comp], tau_0[1:N_comp, 1:N_comp])
  }
  mus[1:N_comp, N_clusters + 1] <- sc_d0[1:N_comp]
  
  tau[1:N_comp, 1:N_comp] ~ dwish(wish_V[1:N_comp, 1:N_comp], N_comp)
}
)
