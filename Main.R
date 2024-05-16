# Main code to:
# Run the metapopulation model for brown rot off peaches described in 
# "A metapopulation framework integrating landscape heterogeneity to model an 
# airborne plant pathogen: the case of brown rot of peach in France"
# And compute epidemic isk indices (vulnerbaility, dangerousness)

rm(list = ls())

library(dplyr)
library(sf)
library(Matrix)

setwd("C:/Users/radya/Desktop/Alcuni file permanenti/Dottorato/Codice/07_METAPOP_model/Codice_pulito_rep")

#### LOAD DATA, PARAMETERS AND SETTINGS

domain <- st_read("Data/Domain.shp") #shp of the domain

n_cells = nrow(domain) #number of cells in the domain
reps = 100 #number of stochastic repetitions [100]

### Areas
A <- as.numeric(st_area(domain)) #area of the cells (in m²)
Ac <- domain$Peach_h*10^4 #cultivated surfaces of the cells (in ha)

### times and settings
starting_years = 1991:2010  #set of the possible starting years: min (1981), max (2021 - ly)
l_y = 10 # simulation horizon (number of years)
Dt = 1 # 1 day, integration step of the epidemiological model


#Additional parameters
V_t0 = 15 #initial susceptible fruit load
E_t0 = 0.27 #initial possible exposed fruit load
#C_N8 = 0.053 # average connectivity towards neighbouring cells (will modify connectivity networks C)
Ptilde0 = 0.2 # probability of having an inoculum at t0 due to external sources

# Load calibrated parameters theta
load("Data/thetas.RData")

#### RUN THE METAPOPULATION MODEL IN STOCHASTC MODE
library(doParallel)
library(foreach)
library(parallel)
library(deSolve)

# Cycle over the Monte Carlo repetitions (index kk, from 1 to REP)
doParallel::registerDoParallel(cl)
foreach(kk= 1:reps, .packages=c("dplyr", "Matrix", "pracma", "deSolve")) %dopar%{
  
  # Extract starting years, one of each cell
  starting_years_sampled = sample(starting_years, n_cells, replace = TRUE)
  
  # Extract theta, one of each cell
  theta_E_sampled = sample(thetas$theta_E, n_cells, replace = TRUE) 
  theta_O_sampled = sample(thetas$theta_O, n_cells, replace = TRUE) 
  
  #Load connectivity matrices (a list of 41 lists for the 1981-2021 connectivity matrices)
  # they are already corrected with rain
  load(paste0("Data/List_C_", (1 + (kk-1) %% 10),".RData"))
  
  #Load Parameters (a list of 41 lists for the 1981-2021 parameters)
  load(paste0("Data/List_params_",  (1 + (kk-1) %% 100) ,".RData"))
  # parameters, per each year,
  # eta = matrix of the temperature-dependent spore mortality rate
  # lambda = transmission rate (beta in eq. 1)
  # rho =  infection-related abscission rate (alpha in eq. 1)
  # sigma =  matrix of the rain-dependent infection rate
  # g, d, b = parameters of the natural abscission rate
  # t_B = vector of the blooming dates
  # t_PH = vector of the pit hardening dates
  # t_0 = earliest pit hardening date among cells
  # t_sim = vectors with the days of the year corresponding to the integration horizon
  # cpd = vector of the identifiers of cells successfully accomplishing phenological requirements
  # GDD_c_ec = vector of the Growing Degree Days of the cultivated variety. 678 = early, 1026 = midearly, 1371 = midlate, 1772 = late
  # w_H = weight of the fruit at harvest time
  
  # dataframe in which values of susceptible fruits will be stored
  S_df <- data.frame(matrix(NA, ncol = n_cells+5, nrow = n_cells)) # Only susceptible
  names(S_df) <- c("id_inoculum_0", "kk", "year0", "theta_E", "theta_O", sapply(1:n_cells, function(x){return(paste0("cell_", x))}))
  
  S_df$id_inoculum0 = 1:n_cells # first inoculated cell
  S_df$kk = kk # id of the reperition
  S_df$year0 = starting_years_sampled # first year of the simulation
  S_df$theta_E = theta_E_sampled #tracks the value of the theta_E used for each randomization
  S_df$theta_O = theta_O_sampled #tracks the value of the theta_O used for each randomization
  
  # dataframe in which values of uninfected susceptible fruits will be stored (for losses calculation)
  V_df <- S_df 
  
  for(jj in 1:n_cells){
    
    id_inoculum_0 = jj #inocultated cell
    ow_possible<- rep(0, n_cells) #vector activating the possibility of haveing a random inoculum at the beginning of the ripening season; it is tuned to one once an infection occurs
    
    theta_E = S_df$theta_E[jj] # value of theta_E used in this simulation
    theta_O = S_df$theta_O[jj] # value of theta_O used in this simulation
    
    years = S_df$year0[jj] + 1: l_y - 1 #simulated years

    for(yy in years){
      
      
      # unlist daily connectivity matrices for the simulated year
      list_C <- lapply(X = list_C_81_21[[yy-1980]], FUN = as.matrix)
      
      # unlist parameters for the simulated year
      list2env(list_PAR_81_21[[yy-1980]], environment())
      
      # compute the lambda for external inoculation
      lambda_b = lambda * theta_E
      
      #identify and removes cells where peach cannot be cultivated
      id_NA = which(cpd==FALSE) 
      
      # Set ODEs
      source("Ode_system.R", local = TRUE) 
      
      # lauch simulation within the simulated year
      source("Stoch_mod.R", local = TRUE)

    }
    
    S_df[jj,5+1:n_cells] <-  S[lt_sim,] #susceptible in the last year are stored 
    V_df[jj,5+1:n_cells] <-  V[lt_sim,] #virtual in the last year are stored 

  }

  save(S_df, V_df, file = 
         paste0("Results/Iteration_", kk,".RData"))  

}

#### COMPUTE EPIDEMIC RISK INDICES AND PLOT RESULTS

library(ggplot2)
# library(pracma)
# library(reshape2)
# library(reshape)

file.names = list.files(path = "Results/", pattern = "\\.RData$")

# re-load all saved dataframes in one only
V_tot <- data.frame(matrix(ncol = n_cells+5, nrow = 0)) 
names(V_tot) <-  c("id_inoculum_0", "kk", "year0", "theta_E", "theta_O", sapply(1:n_cells, function(x){return(paste0("cell_", x))}))
S_tot <- V_tot

for(i in 1:length(file.names)){
  load(paste0("Results/", file.names[i]))
  S_tot  <- rbind(S_tot , S_df)
  V_tot  <- rbind(V_tot , V_df)
}

# compute losses in each cell in each simulation (Eq. 8)
L_tot = S_tot #losses
L_tot[, 5+1:n_cells] = 1 - S_tot[, 5+1:n_cells]/V_tot[, 5+1:n_cells] #losses

# summarizes losses provoked in each cell by an infection in any other cell
L_m <- matrix(NA, ncol = n_cells, nrow = n_cells)

for(i in 1:n_cells){
  L_df <- L_tot %>%
    filter(id_inoculum_0 == i)
  
  L_m[i,]  = 100*colSums(as.matrix(L_df[, 5+1:n_cells]), na.rm = T)/colSums(as.matrix(L_df[, 5+1:n_cells])>-1, na.rm = T)
}

# avoid autoinfection
diag(L_m) <- 0 

# losses are aggregated into vulnerability
vulnerability <- colSums(L_m, na.rm = T)/colSums(L_m>-1, na.rm = T)

# Compute area-weighted losses  L_A_m

Ac_v <- domain$Peach_h
Ac_m <- matrix(rep(Ac_v, times = n_cells), nrow = n_cells, byrow = T)
L_A_m <- L_m*Ac_m

# avoid autoinfection
diag(L_A_m) = 0

#dangerousenss
dangerousness <- rowSums(L_A_m, na.rm = T)/(sum(Ac_v)-Ac_v)

#  Plots

ggplot(data = domain)+
  geom_sf(aes(fill = vulnerability))+
  scale_fill_gradient(low = "white",
                      high = "red",
                      name = "% losses",
                      na.value = "grey70") +
  theme(panel.grid.major = element_line(color = gray(.5), 
                                        linetype = "dashed", size = 0.5),
        panel.background = element_rect(fill = "aliceblue"))+
  ggtitle("Vulnerability")

ggplot(data = domain)+
  geom_sf(aes(fill = dangerousness))+
  scale_fill_gradient(low = "white",
                      high = "red",
                      name = "% losses",
                      na.value = "grey70") +
  theme(panel.grid.major = element_line(color = gray(.5), 
                                        linetype = "dashed", size = 0.5),
        panel.background = element_rect(fill = "aliceblue"))+
  ggtitle("Dangerousness")
