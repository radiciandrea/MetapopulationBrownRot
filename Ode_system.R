# parameters to be passed to the ODE system
parms <- list(eta = eta,
              lambda = lambda,
              rho = rho,
              sigma = sigma,
              t_B = t_B,
              t_0 = t_0,
              g = g,
              d = d,
              b = b,
              t_PH = t_PH,
              t_H = t_H) 

df <- function(t, x, parms) {
  
  # initial conditions and paramters
  with(parms, { 
    S <- x[1:(length(x)/5)]
    E <- x[(length(x)/5 +1):(2*length(x)/5)]
    I <- x[(2*length(x)/5+1):(3*length(x)/5)]
    V <- x[(3*length(x)/5+1):(4*length(x)/5)]
    M <- x[(4*length(x)/5+1):length(x)]
  
    # time
    tn <- t_0 + t[1] - 1
    
    # calculate abscission
    N = S + E + I
    mu = pmax((g*(t_0 + t[1] - 1-t_B)^2 + d*(t_0 + t[1] - 1 - t_B))*N^b, 0)
    mu_v = pmax((g*(t_0 + t[1] - 1-t_B)^2 + d*(t_0 + t[1] - 1 - t_B))*V^b, 0)

    # ODE definition 
    dS <- (tn >t_PH) * (tn <t_H) * (- lambda * S * I + eta[,t[1]] * E - mu * S) 
    dE <- (tn >t_PH) * (tn <t_H) * (+ lambda * S * I - (eta[,t[1]] + sigma[,t[1]] + mu) * E)  
    dI <- (tn >t_PH) * (tn <t_H) * (+ sigma[,t[1]] * E - (rho + mu) * I)
    dM <- (tn >t_PH) * (tn <t_H) * (rho + mu) * I
    
    # Virtual class of uninfected fruits (S*) to compute losses
    dV <- (tn >t_PH) * (tn <t_H) * (- mu_v * V)
    dx <- c(dS, dE, dI, dV, dM)
    
    return(list(dx))})
  
}
