
###########################################################################
# This code was created by the DARTH workgroup (www.darthworkgroup.com). 
# When using or modifying this code, please do so with attribution and 
# cite our publications:

# - Alarid-Escudero F, Maclehose RF, Peralta Y, Kuntz KM, Enns EA. 
#   Non-identifiability in model calibration and implications for 
#   medical decision making. Med Decis Making. 2018; 38(7):810-821.

# - Jalal H, Pechlivanoglou P, Krijkamp E, Alarid-Escudero F, Enns E, 
#   Hunink MG. An Overview of R in Health Decision Sciences. 
#   Med Decis Making. 2017; 37(3): 735-746. 

# A walkthrough of the code could be found in the follwing link:
# - https://darth-git.github.io/calibSMDM2018-materials/
###########################################################################


###################  Calibration Specifications  ###################

# Model: Sick-Sicker 4-state Markov Model
# Inputs to be calibrated: p.S1S2, hr.S1, hr.S2
# Targets: Surv, Prev, PropSick

# Search method: Random search using Latin-Hypercube Sampling
# Goodness-of-fit measure: Sum of Log-Likelihood

####################################################################


####################################################################
######  Load packages and function files  ####
####################################################################
# calibration functionality
library(lhs)

# visualization
library(plotrix)
library(psych)
library(scatterplot3d) # now that we have three inputs to estimate, we'll need higher dimension visualization


####################################################################
######  Load target data  ######
####################################################################
load("SickSicker_CalibTargets.RData")
lst_targets <- SickSicker_targets

# Plot the targets

# TARGET 1: Survival ("Surv")
plotrix::plotCI(x = lst_targets$Surv$time, y = lst_targets$Surv$value, 
                ui = lst_targets$Surv$ub,
                li = lst_targets$Surv$lb,
                ylim = c(0, 1), 
                xlab = "Time", ylab = "Pr Survive")

# TARGET 2: Prevalence ("Prev")
plotrix::plotCI(x = lst_targets$Prev$time, y = lst_targets$Prev$value,
                ui = lst_targets$Prev$ub,
                li = lst_targets$Prev$lb,
                ylim = c(0, 1),
                xlab = "Time", ylab = "Prev")

# TARGET 3: Proportion who are Sick ("PropSick"), among all those afflicted (Sick+Sicker)
plotrix::plotCI(x = lst_targets$PropSick$time, y = lst_targets$PropSick$value,
                ui = lst_targets$PropSick$ub,
                li = lst_targets$PropSick$lb,
                ylim = c(0, 1),
                xlab = "Time", ylab = "PropSick")


####################################################################
######  Load model as a function  ######
####################################################################
# - inputs are parameters to be estimated through calibration
# - outputs correspond to the target data

source("SickSicker_MarkovModel_Function.R") # creates the function run_sick_sicker_markov()

# Check that it works
v_params_test <- c(p_S1S2 = 0.105, hr_S1 = 3, hr_S2 = 10)
run_sick_sicker_markov(v_params_test) # It works!


####################################################################
######  Specify calibration parameters  ######
####################################################################
# Specify seed (for reproducible sequence of random numbers)
set.seed(072218)

# number of initial starting points
n_init <- 100

# names and number of input parameters to be calibrated
v_param_names <- c("p_S1S2","hr_S1","hr_S2")
n_param <- length(v_param_names)

# range on input search space
lb <- c(p_S1S2 = 0.01, hr_S1 = 1.0, hr_S2 = 5) # lower bound
ub <- c(p_S1S2 = 0.50, hr_S1 = 4.5, hr_S2 = 15) # upper bound

# number of calibration targets
v_target_names <- c("Surv", "Prev", "PropSick")
n_target <- length(v_target_names)


### Calibration functions

# Write goodness-of-fit function to pass to Nelder-Mead algorithm
f_gof <- function(v_params){
  
  # Run model for parametr set "v_params"
  model_res <- run_sick_sicker_markov(v_params)
  
  # Calculate goodness-of-fit of model outputs to targets
  v_GOF <- numeric(n_target)
  # TARGET 1: Survival ("Surv")
  # log likelihood  
  v_GOF[1] <- sum(dnorm(x = lst_targets$Surv$value,
                        mean = model_res$Surv,
                        sd = lst_targets$Surv$se,
                        log = T))
  
  # TARGET 2: "Prev"
  # log likelihood
  v_GOF[2] <- sum(dnorm(x = lst_targets$Prev$value,
                        mean = model_res$Prev,
                        sd = lst_targets$Prev$se,
                        log = T))
  
  # TARGET 3: "PropSick"
  # log likelihood
  v_GOF[3] <- sum(dnorm(x = lst_targets$PropSick$value,
                        mean = model_res$PropSick,
                        sd = lst_targets$PropSick$se,
                        log = T))
  
  # OVERALL
  # can give different targets different weights
  v_weights <- rep(1,n_target)
  # weighted sum
  GOF_overall <- sum(v_GOF[1:n_target] * v_weights)
  
  # return GOF
  return(GOF_overall)
}


####################################################################
######  Calibrate!  ######
####################################################################
# record start time of calibration
t_init <- Sys.time()

###  Sample multiple random starting values for Nelder-Mead  ###
v_params_init <- matrix(nrow=n_init,ncol=n_param)
for (i in 1:n_param){
  v_params_init[,i] <- runif(n_init,min=lb[i],max=ub[i])
}
colnames(v_params_init) <- v_param_names

###  Run Nelder-Mead for each starting point  ###
m_calib_res <- matrix(nrow = n_init, ncol = n_param+1)
colnames(m_calib_res) <- c(v_param_names, "Overall_fit")
for (j in 1:n_init){
  
  ### use optim() as Nelder-Mead ###
  fit_nm <- optim(v_params_init[j,], f_gof,
                 control = list(fnscale = -1, # switches from minimization to maximization
                                maxit = 1000), hessian = T)
  m_calib_res[j,] <- c(fit_nm$par,fit_nm$value)
  
  ### to use a simulated annealing instead ###
  # fit_sa <- optim(v_params_init[j,], f_gof,
  #                method = c("SANN"),  # switches to using simulated annealing
  #                control = list(temp = 10, tmax = 10, # algorithm tuning parameters
  #                               fnscale = -1, maxit = 1000),
  #                hessian = T)
  # m_calib_res[j,] = c(fit_sa$par,fit_sa$value)
  
  ### to use a genetic algorithm instead ###
  # library(DEoptim)
  # f_fitness <- function(params){
  #   names(params) = v_param_names
  #   return(-f_gof(params))}
  # fit_ga = DEoptim(f_fitness, lower=lb, upper=ub)
  # m_calib_res[j,] = c(fit_ga$optim$bestmem,-1*fit_ga$optim$bestval)
  
}

# Calculate computation time
comp_time <- Sys.time() - t_init


####################################################################
######  Exploring best-fitting input sets  ######
####################################################################

# Arrange parameter sets in order of fit
m_calib_res <- m_calib_res[order(-m_calib_res[,"Overall_fit"]),]

# Examine the top 10 best-fitting sets
m_calib_res[1:10,]

# Plot the top 10 (top 10%)
scatterplot3d(x = m_calib_res[1:10, 1],
              y = m_calib_res[1:10, 2],
              z = m_calib_res[1:10, 3],
              xlim = c(lb[1],ub[1]), ylim = c(lb[2],ub[2]), zlim = c(lb[3],ub[3]),
              xlab = v_param_names[1], ylab = v_param_names[2], zlab = v_param_names[3])

# Pairwise comparison of top 10 sets
pairs.panels(m_calib_res[1:10,v_param_names])

### Plot model-predicted output at best set vs targets ###
v_out_best <- run_sick_sicker_markov(m_calib_res[1,])

# TARGET 1: Survival ("Surv")
plotrix::plotCI(x = lst_targets$Surv$time, y = lst_targets$Surv$value, 
                ui = lst_targets$Surv$ub,
                li = lst_targets$Surv$lb,
                ylim = c(0, 1), 
                xlab = "Time", ylab = "Pr Survive")
points(x = lst_targets$Surv$time, 
       y = v_out_best$Surv, 
       pch = 8, col = "red")
legend("topright", 
       legend = c("Target", "Model-predicted output"),
       col = c("black", "red"), pch = c(1, 8))

# TARGET 2: "Prev"
plotrix::plotCI(x = lst_targets$Prev$time, y = lst_targets$Prev$value,
                ui = lst_targets$Prev$ub,
                li = lst_targets$Prev$lb,
                ylim = c(0, 1),
                xlab = "Time", ylab = "Prev")
points(x = lst_targets$Prev$time,
       y = v_out_best$Prev,
       pch = 8, col = "red")
legend("topright",
       legend = c("Target", "Model-predicted output"),
       col = c("black", "red"), pch = c(1, 8))

# TARGET 3: "PropSick"
plotrix::plotCI(x = lst_targets$PropSick$time, y = lst_targets$PropSick$value,
                ui = lst_targets$PropSick$ub,
                li = lst_targets$PropSick$lb,
                ylim = c(0, 1),
                xlab = "Time", ylab = "PropSick")
points(x = lst_targets$PropSick$time,
       y = v_out_best$PropSick,
       pch = 8, col = "red")
legend("topright",
       legend = c("Target", "Model-predicted output"),
       col = c("black", "red"), pch = c(1, 8))

