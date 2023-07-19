#### Sick-Sicker Markov model in a function ####
run_sick_sicker_markov <- function(v_params) {
  with(as.list(v_params), {
    ## Markov model parameters
    age     <- 25                   # age at baseline
    max_age <- 55                   # maximum age of follow up
    n_t  <- max_age - age           # time horizon, number of cycles
    v_n  <- c("H", "S1", "S2", "D") # the 4 states of the model: Healthy (H), Sick (S1), Sicker (S2), Dead (D)
    n_s <- length(v_n)              # number of health states 
    
    # Transition probabilities and hazard ratios
    p_HD    = 0.005 # probability to die when healthy
    p_HS1   = 0.15  # probability to become sick when healthy
    p_S1H   = 0.5   # probability to become healthy when sick
    # p_S1S2  = 0.105         	# probability to become sicker when sick
    # hr_S1   = 3           	  # hazard ratio of death in sick vs healthy
    # hr_S2   = 10          	  # hazard ratio of death in sicker vs healthy
    # compute internal paramters as a function of external parameters
    r_HD    = - log(1 - p_HD) # rate of death in healthy
    r_S1D   = hr_S1 * r_HD 	  # rate of death in sick
    r_S2D   = hr_S2 * r_HD  	# rate of death in sicker
    p_S1D   = 1 - exp(-r_S1D) # probability to die in sick
    p_S2D   = 1 - exp(-r_S2D) # probability to die in sicker
    
    ####### INITIALIZATION ##########################################
    # create the cohort trace
    m_M <- matrix(NA, nrow = n_t + 1 , 
                  ncol = n_s,
                  dimnames = list(0:n_t, v_n))     # create Markov trace (n_t + 1 because R doesn't understand  Cycle 0)
    
    m_M[1, ] <- c(1, 0, 0, 0)                      # initialize Markov trace
    
    # create transition probability matrix for NO treatment
    m_P <- matrix(0,
                  nrow = n_s, 
                  ncol = n_s,
                  dimnames = list(v_n, v_n))
    # fill in the transition probability array
    ### From Healthy
    m_P["H", "H"]  <- 1 - (p_HS1 + p_HD)
    m_P["H", "S1"] <- p_HS1
    m_P["H", "D"]  <- p_HD
    ### From Sick
    m_P["S1", "H"]  <- p_S1H
    m_P["S1", "S1"] <- 1 - (p_S1H + p_S1S2 + p_S1D)
    m_P["S1", "S2"] <- p_S1S2
    m_P["S1", "D"]  <- p_S1D
    ### From Sicker
    m_P["S2", "S2"] <- 1 - p_S2D
    m_P["S2", "D"]  <- p_S2D
    ### From Dead
    m_P["D", "D"] <- 1
    
    # check rows add up to 1
    if (!isTRUE(all.equal(as.numeric(rowSums(m_P)), as.numeric(rep(1, n_s))))) {
      stop("This is not a valid transition Matrix")
    }
    
    ############# PROCESS ###########################################
    
    for (t in 1:n_t){                              # throughout the number of cycles
      m_M[t + 1, ] <- m_M[t, ] %*% m_P           # estimate the Markov trace for cycle the next cycle (t + 1)
    }
    
    ####### EPIDEMIOLOGICAL OUTPUT  ###########################################
    #### Overall Survival (OS) ####
    v_os <- 1 - m_M[, "D"]                # calculate the overall survival (OS) probability for no treatment
    
    #### Disease prevalence #####
    v_prev <- rowSums(m_M[, c("S1", "S2")])/v_os
    
    #### Proportion of sick in S1 state #####
    v_prop_S1 <- m_M[, "S1"] / v_prev
    
    ####### RETURN OUTPUT  ###########################################
    out <- list(Surv = v_os[-1],
                Prev = v_prev[-1],
                PropSick = v_prop_S1[c(11, 21, 31)])
    
    return(out)
  }
  )
}