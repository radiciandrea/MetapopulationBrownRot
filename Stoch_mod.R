if(yy == years[1]){
  id_inoculum = id_inoculum_0 # sample(set_inoculum0,1)
} 

# Initial values - inoculum everywhere
V_0 <- V_t0*rep(1, times = n_cells) # virtual susceptible load, to calculate losses
S_0 <- V_0  #susceptible load
E_0 <- 0*rep(1, times = n_cells) # exposed load, set to 0
I_0 <- 0*rep(1, times = n_cells) # infected load, set to 0  
M_0 <- 0*rep(1, times = n_cells)  # mummifies load, set to 0

# Inoculated cell is completely infected (except a fraction, which is "exposed" to avoid external inputs)
S_0[id_inoculum] = 0
E_0[id_inoculum] = E_t0
I_0[id_inoculum] = V_t0 - E_t0

# load is set to zero in cells where there
V_0[id_NA] <- 0 #
S_0[id_NA] <- 0 # 
E_0[id_NA] <- 0 # 
I_0[id_NA] <- 0 # 
M_0[id_NA] <- 0 # 

# Initialization of the whole matrix (each cell, each day)
V <- matrix(rep(V_0, lt_sim), byrow = T, nrow = lt_sim) 
S <- matrix(rep(S_0, lt_sim), byrow = T, nrow = lt_sim)
E <- matrix(rep(E_0, lt_sim), byrow = T, nrow = lt_sim)
I <- matrix(rep(I_0, lt_sim), byrow = T, nrow = lt_sim)
M <- matrix(rep(M_0, lt_sim), byrow = T, nrow = lt_sim) 

# Initialization of the first day (each cell)
S[1,] <- S_0
E[1,] <- E_0
I[1,] <- I_0
V[1,] <- V_0 
M[1,] <- M_0 

#vector of variables to be integrated
X_0 <- c(S_0, E_0, I_0, V_0, M_0)

# Initialization of the rate of external infection
Rtt <- matrix(0, ncol = n_cells, nrow = (lt_sim-1)) 

# Emply matrix the newly externally infected units
NE <- matrix(NA, ncol = n_cells, nrow = (lt_sim-1))

#Stochastic extraction for external infection
U <- matrix(runif(lt_sim*n_cells), ncol = n_cells)

for(tt in 1:(lt_sim-1)){
  
  #simulation time horizon (1 day)
  t_sim_t <- tt:(tt+Dt)
  
  #integration
  Sim <- ode(X_0, t_sim_t, df, parms)
  
  #splitting of the integrated values inside the matrices
  S[tt+1,] <- as.vector(Sim[2,2:(n_cells+1)])
  E[tt+1,] <- as.vector(Sim[2,(n_cells +2):(2*n_cells +1)])
  I[tt+1,] <- as.vector(Sim[2,(2*n_cells +2):(3*n_cells + 1)])
  V[tt+1,] <- as.vector(Sim[2,(3*n_cells +2):(4*n_cells + 1)])
  M[tt+1,] <- as.vector(Sim[2,(4*n_cells +2):(5*n_cells + 1)])
  
  # Computation of the rate of external infection (Eq. 7)
  Rtt[tt,] = (tt + t_0 -1>t_PH) *  (tt + t_0 -1 < t_H) *(E[tt,] == 0 ) * (S[tt,] > 0)* lambda_b * S[tt,] *list_C[[tt]] %*% (I[tt,] * Ac) 
  
  # Probability of external infection (Eq. 6)
  Pe_t = 1 - exp(-Rtt[tt,]*Dt) 
  
  # Newly infected units
  NE[tt,] = U[tt,]<Pe_t
  
  # Newly infected units are inocultated: A fraction E_t0 becomes exposed
  
  E_t0_eff <- pmin(E_t0, S[tt+1,NE[tt,]]) # check to avoid negative values
  S[tt+1,NE[tt,]] <- S[tt+1,NE[tt,]] - E_t0_eff
  E[tt+1,NE[tt,]] <- E[tt+1,NE[tt,]] + E_t0_eff
  
  X_0 <- c(S[tt+1,], E[tt+1,], I[tt+1,], V[tt+1,], M[tt+1,])
}

# Mummies are computed at the end of the cycle (Eq. 2)
M_tot = M[tt+1, ] + I[tt+1,]

# probability of external infection is updated
ow_possible = pmax(ow_possible, M_tot>0)

# We extract the units what will start the following year as "exposed" bue to overwintering and external infection (Eq. 5)
id_inoculum <- which(runif(n_cells) < ow_possible*Ptilde0 + (1-Ptilde0)*(1-exp(-theta_O*M_tot)))

S[,id_NA] <- NA
V[,id_NA] <- NA