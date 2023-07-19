
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

# number of random samples
n_samp <- 1000

# names and number of input parameters to be calibrated
v_param_names <- c("p_S1S2","hr_S1","hr_S2")
n_param <- length(v_param_names)

# range on input search space
lb <- c(p_S1S2 = 0.01, hr_S1 = 1.0, hr_S2 = 5) # lower bound
ub <- c(p_S1S2 = 0.50, hr_S1 = 4.5, hr_S2 = 15) # upper bound

# number of calibration targets
v_target_names <- c("Surv", "Prev", "PropSick")
n_target <- length(v_target_names)


####################################################################
######  Calibrate!  ######
####################################################################
# record start time of calibration
t_init <- Sys.time()

###  Generate a random sample of input values  ###

# Sample unit Latin Hypercube
m_lhs_unit <- randomLHS(n_samp, n_param)

# Rescale to min/max of each parameter
m_param_samp <- matrix(nrow=n_samp,ncol=n_param)
for (i in 1:n_param){
  m_param_samp[,i] <- qunif(m_lhs_unit[,i],
                           min = lb[i],
                           max = ub[i])
}
colnames(m_param_samp) <- v_param_names

# view resulting parameter set samples
pairs.panels(m_param_samp)


###  Run the model for each set of input values ###

# initialize goodness-of-fit vector
m_GOF <- matrix(nrow = n_samp, ncol = n_target)
colnames(m_GOF) <- paste0(v_target_names, "_fit")

# loop through sampled sets of input values
for (j in 1:n_samp){
  
  ###  Run model for a given parameter set  ###
  model_res <- run_sick_sicker_markov(v_params = m_param_samp[j, ])
  
  
  ###  Calculate goodness-of-fit of model outputs to targets  ###

  # TARGET 1: Survival ("Surv")
  # log likelihood  
  m_GOF[j,1] <- sum(dnorm(x = lst_targets$Surv$value,
                       mean = model_res$Surv,
                       sd = lst_targets$Surv$se,
                       log = T))
  
  # weighted sum of squared errors (alternative to log likelihood)
  # w <- 1/(lst_targets$Surv$se^2)
  # m_GOF[j,1] <- -sum(w*(lst_targets$Surv$value - v_res)^2)
  
  
  # TARGET 2: "Prev"
  # log likelihood
  m_GOF[j,2] <- sum(dnorm(x = lst_targets$Prev$value,
                         mean = model_res$Prev,
                         sd = lst_targets$Prev$se,
                         log = T))
  
  # TARGET 3: "PropSick"
  # log likelihood
  m_GOF[j,3] <- sum(dnorm(x = lst_targets$PropSick$value,
                         mean = model_res$PropSick,
                         sd = lst_targets$PropSick$se,
                         log = T))
  
  
} # End loop over sampled parameter sets


###  Combine fits to the different targets into single GOF  ###
# can give different targets different weights
v_weights <- matrix(1, nrow = n_target, ncol = 1)
# matrix multiplication to calculate weight sum of each GOF matrix row
v_GOF_overall <- c(m_GOF%*%v_weights)
# Store in GOF matrix with column name "Overall"
m_GOF <- cbind(m_GOF,Overall_fit=v_GOF_overall)

# Calculate computation time
comp_time <- Sys.time() - t_init


####################################################################
######  Exploring best-fitting input sets  ######
####################################################################

# Arrange parameter sets in order of fit
m_calib_res <- cbind(m_param_samp,m_GOF)
m_calib_res <- m_calib_res[order(-m_calib_res[,"Overall_fit"]),]

# Examine the top 10 best-fitting sets
m_calib_res[1:10,]

# Plot the top 100 (top 10%)
scatterplot3d(x = m_calib_res[1:100, 1],
              y = m_calib_res[1:100, 2],
              z = m_calib_res[1:100, 3],
              xlim = c(lb[1],ub[1]), ylim = c(lb[2],ub[2]), zlim = c(lb[3],ub[3]),
              xlab = v_param_names[1], ylab = v_param_names[2], zlab = v_param_names[3])

# Pairwise comparison of top 100 sets
pairs.panels(m_calib_res[1:100,v_param_names])

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

