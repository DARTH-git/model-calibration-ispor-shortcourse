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

# A walkthrough of the code could be found in the following link:
# - https://darth-git.github.io/calibSMDM2018-materials/
###########################################################################

###################  Calibration Specifications  ###################

# Model: 3-State Cancer Relative Survival (CRS) Markov Model
# Inputs to be calibrated: p_Mets, p_DieMets
# Targets: Surv

# Calibration method: Incremental mixture importance sampling (IMIS)
# Goodness-of-fit measure: Sum of Log-Likelihood

####################################################################

####################################################################
rm(list = ls())

####################################################################
######  Load packages and function files  ####
####################################################################
# calibration functionality
library(lhs)
library(IMIS)
library(matrixStats) # package used for summary statistics

# visualization
library(plotrix)
library(psych)
library(GGally)

# data wrangling
library(dplyr)

# parallel computing
library(doParallel)

####################################################################
######  Load target data  ######
####################################################################
load("ISPOR Calibration Participant Materials/CRS_CalibTargets.RData")
lst_targets <- CRS_targets

# Plot the targets

# TARGET 1: Survival ("Surv")
plotrix::plotCI(x = lst_targets$Surv$time, y = lst_targets$Surv$value, 
                ui = lst_targets$Surv$ub,
                li = lst_targets$Surv$lb,
                ylim = c(0, 1), 
                xlab = "Time", ylab = "Pr Survive")

# TARGET 2: (if you had more...)
# plotrix::plotCI(x = lst_targets$Target2$time, y = lst_targets$Target2$value, 
#                 ui = lst_targets$Target2$ub,
#                 li = lst_targets$Target2$lb,
#                 ylim = c(0, 1), 
#                 xlab = "Time", ylab = "Target 2")


####################################################################
######  Load model as a function  ######
####################################################################
# - inputs are parameters to be estimated through calibration
# - outputs correspond to the target data

source("ISPOR Calibration Participant Materials/CRS_MarkovModel_Function.R") # creates the function run_crs_markov()

# Check that it works
v_params_test <- c(p_Mets = 0.10, p_DieMets = 0.05)
run_crs_markov(v_params_test) # It works!


####################################################################
######  Specify calibration parameters  ######
####################################################################
# Specify seed (for reproducible sequence of random numbers)
set.seed(072218)

# number of random samples
n_resamp <- 1000

# names and number of input parameters to be calibrated
v_param_names <- c("p_Mets","p_DieMets")
n_param <- length(v_param_names)

# range on input search space
lb <- c(p_Mets = 0.04, p_DieMets = 0.04) # lower bound
ub <- c(p_Mets = 0.16, p_DieMets = 0.16) # upper bound

# number of calibration targets
v_target_names <- c("Surv")
n_target     <- length(v_target_names)


### Calibration functions

#  Write function to sample from prior
sample_prior <- function(n_samp){ # n_samp <- 5
  m_lhs_unit   <- randomLHS(n = n_samp, k = n_param)
  m_param_samp <- matrix(nrow = n_samp, ncol = n_param)
  colnames(m_param_samp) <- v_param_names
  for (i in 1:n_param) {
    m_param_samp[, i] <- qunif(m_lhs_unit[,i],
                               min = lb[i],
                               max = ub[i])
    # ALTERNATIVE prior using beta (or other) distributions
    # m_param_samp[, i] <- qbeta(m_lhs_unit[,i],
    #                            shape1 = 1,
    #                            shape2 = 1)
  }
  return(m_param_samp)
}

# view resulting parameter set samples
pairs.panels(sample_prior(1000))


###  PRIOR  ### 
# Write functions to evaluate log-prior and prior

# function that calculates the log-prior
calc_log_prior <- function(v_params){
  if (is.null(dim(v_params))) { # If vector, change to matrix
    v_params <- t(v_params) 
  }
  n_samp <- nrow(v_params)
  colnames(v_params) <- v_param_names
  lprior <- rep(0, n_samp)
  for (i in 1:n_param) {
    lprior <- lprior + dunif(v_params[, i],
                             min = lb[i],
                             max = ub[i], 
                             log = T)
    # ALTERNATIVE prior using beta distributions
    # lprior <- lprior + dbeta(v_params[, i],
    #                          shape1 = 1,
    #                          shape2 = 1, 
    #                          log = T)
  }
  return(lprior)
}
calc_log_prior(v_params = v_params_test)
calc_log_prior(v_params = sample_prior(10))


# function that calculates the (non-log) prior
calc_prior <- function(v_params) { 
  exp(calc_log_prior(v_params)) 
}
calc_prior(v_params = v_params_test)
calc_prior(v_params = sample_prior(10))

###  LIKELIHOOD  ###
# Write functions to evaluate log-likelihood and likelihood

# function to calculate the log-likelihood
calc_log_lik <- function(v_params){
  # par_vector: a vector (or matrix) of model parameters 
  if (is.null(dim(v_params))) { # If vector, change to matrix
    v_params <- t(v_params) 
  }
  n_samp <- nrow(v_params)
  v_llik <- matrix(0, nrow = n_samp, ncol = n_target) 
  llik_overall <- numeric(n_samp)
  for (j in 1:n_samp) { # j=1
    jj <- tryCatch( { 
      ###   Run model for parametr set "v_params" ###
      model_res <- run_crs_markov(v_params[j, ])
      
      ###  Calculate log-likelihood of model outputs to targets  ###
      # TARGET 1: Survival ("Surv")
      # log likelihood  
      v_llik[j, 1] <- sum(dnorm(x = lst_targets$Surv$value,
                                mean = model_res$Surv,
                                sd = lst_targets$Surv$se,
                                log = T))
      
      # TARGET 2: (if you had more...)
      # log likelihood
      # v_llik[j, 2] <- sum(dnorm(x = lst_targets$Target2$value,
      #                        mean = model_res$Target2,
      #                        sd = lst_targets$Target2$se,
      #                        log = T))
      
      # OVERALL 
      llik_overall[j] <- sum(v_llik[j, ])
    }, error = function(e) NA) 
    if (is.na(jj)) { llik_overall <- -Inf }
  } # End loop over sampled parameter sets
  # return LLIK
  return(llik_overall)
}
calc_log_lik(v_params = v_params_test)
calc_log_lik(v_params = sample_prior(10))


# function to calculate the (non-log) likelihood
calc_likelihood <- function(v_params){ 
  exp(calc_log_lik(v_params)) 
}
calc_likelihood(v_params = v_params_test)
calc_likelihood(v_params = sample_prior(10))


###  POSTERIOR  ###
# Write functions to evaluate log-posterior and posterior

# function that calculates the log-posterior
calc_log_post <- function(v_params) { 
  lpost <- calc_log_prior(v_params) + calc_log_lik(v_params)
  return(lpost) 
}
calc_log_post(v_params = v_params_test)
calc_log_post(v_params = sample_prior(10))


# function that calculates the (non-log) posterior
calc_post <- function(v_params) { 
  exp(calc_log_post(v_params)) 
}
calc_post(v_params = v_params_test)
calc_post(v_params = sample_prior(10))

#' Parallel evaluation of log-likelihood function for a sets of parameters
#'
#' \code{log_lik_par} computes a log-likelihood value for one (or multiple) 
#' parameter set(s) using parallel computation.
#'
#' @param v_params Vector (or matrix) of model parameters.
#' @return 
#' A scalar (or vector) with log-likelihood values.
log_lik_par <- function(v_params,
                        ...) { 
  if (is.null(dim(v_params))) { # If vector, change to matrix
    v_params <- t(v_params) 
  }
  
  n_samp <- nrow(v_params)
  
  ### Get OS
  os <- get_os()
  
  no_cores <- parallel::detectCores() - 1
  
  print(paste0("Parallelized Likelihood calculations on ", os, " using ", no_cores, " cores"))
  
  n_time_init_likpar <- Sys.time()
  
  if (os == "macosx") {
    # Initialize cluster object
    cl <- parallel::makeForkCluster(no_cores) 
    doParallel::registerDoParallel(cl)
    v_llk <- foreach::foreach(i = 1:n_samp, .combine = c) %dopar% {
      calc_log_lik(v_params[i, ]) # i = 1
    }
    n_time_end_likpar <- Sys.time()
  }
  if (os == "windows") {
    # Initialize cluster object
    cl <- parallel::makeCluster(no_cores)
    doParallel::registerDoParallel(cl)
    opts <- list(attachExportEnv = TRUE)
    v_llk <- foreach::foreach(i = 1:n_samp, .combine = c,
                              .export = ls(globalenv()),
                              .packages = c(),
                              .options.snow = opts) %dopar% {
                                calc_log_lik(v_params[i, ])
                              }
    n_time_end_likpar <- Sys.time()
  }
  if (os == "linux") {
    # Initialize cluster object
    cl <- parallel::makeCluster(no_cores)
    doMC::registerDoMC(cl)
    v_llk <- foreach::foreach(i = 1:n_samp, .combine = c) %dopar% {
      calc_log_lik(v_params[i, ])
    }
    n_time_end_likpar <- Sys.time()
  }
  
  parallel::stopCluster(cl)
  n_time_total_likpar <- difftime(n_time_end_likpar, n_time_init_likpar, 
                                  units = "hours")
  print(paste0("Runtime: ", round(n_time_total_likpar, 2), " hrs."))
  #-# Try this: # PO
  rm(cl)        # PO
  gc()          # PO
  #-#           # PO
  return(v_llk)
}

#' Likelihood
#'
#' \code{likelihood} computes a likelihood value for one (or multiple) 
#' parameter set(s).
#'
#' @param v_params Vector (or matrix) of model parameters. 
#' @return 
#' A scalar (or vector) with likelihood values.
likelihood <- function(v_params){ 
  v_like <- exp(log_lik_par(v_params)) 
  return(v_like)
}

#' Get operating system
#' 
#' @return 
#' A string with the operating system.
#' @export
get_os <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)) {
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "MacOSX"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
}


####################################################################
######  Calibrate!  ######
####################################################################
# record start time of calibration
t_init <- Sys.time()

###  Bayesian calibration using IMIS  ###
# define three functions needed by IMIS: prior(x), likelihood(x), sample.prior(n)
prior <- calc_prior
# likelihood <- calc_likelihood
sample.prior <- sample_prior

# run IMIS
fit_imis <- IMIS(B = 1000, # the incremental sample size at each iteration of IMIS
                 B.re = n_resamp, # the desired posterior sample size
                 number_k = 10, # the maximum number of iterations in IMIS
                 D = 0) 

# obtain draws from posterior
m_calib_res <- fit_imis$resample

# Calculate log-likelihood (overall fit) and posterior probability of each sample
m_calib_res <- cbind(m_calib_res, 
                     "Overall_fit" = calc_log_lik(m_calib_res[,v_param_names]),
                     "Posterior_prob" = calc_post(m_calib_res[,v_param_names]))

# normalize posterior probability
m_calib_res[,"Posterior_prob"] <- m_calib_res[,"Posterior_prob"]/sum(m_calib_res[,"Posterior_prob"])

# Calculate computation time
comp_time <- Sys.time() - t_init
comp_time

####################################################################
######  Exploring posterior distribution  ######
####################################################################

# Plot the 1000 draws from the posterior
v_post_color <- scales::rescale(m_calib_res[,"Posterior_prob"])
plot(m_calib_res,
     xlim = c(lb[1], ub[1]), ylim = c(lb[2], ub[2]),
     xlab = v_param_names[1], ylab = v_param_names[2],
     col = scales::alpha("black", v_post_color))
# add center of Gaussian components
points(fit_imis$center, col = "red", pch = 8)
legend("topright", c("Draws from posterior", "Center of Gaussian components"),
       col = c("black", "red"), pch = c(1, 8))

# Plot the 1000 draws from the posterior with marginal histograms
pairs.panels(m_calib_res[,v_param_names])

# Compute posterior mean
v_calib_post_mean <- colMeans(m_calib_res[, v_param_names])
v_calib_post_mean

# Compute posterior median and 95% credible interval
m_calib_res_95cr <- colQuantiles(m_calib_res[,v_param_names], probs = c(0.025, 0.5, 0.975))
m_calib_res_95cr

# Compute maximum-a-posteriori (MAP) parameter set
v_calib_map <- m_calib_res[which.max(m_calib_res[,"Posterior_prob"]),]

### Plot model-predicted output at mode vs targets ###
v_out_best <- run_crs_markov(v_calib_map[v_param_names])
v_out_post_mean <- run_crs_markov(v_calib_post_mean)

# set plot margins (helps plotting)
par(mar = c(5, 4, 4, 4)) 
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

points(x = lst_targets$Surv$time, 
       y = v_out_post_mean$Surv, 
       pch = 8, col = "blue")

# TARGET 2: (if you had more...)
# plotrix::plotCI(x = lst_targets$Target2$time, y = lst_targets$Target2$value, 
#                 ui = lst_targets$Target2$ub,
#                 li = lst_targets$Target2$lb,
#                 ylim = c(0, 1), 
#                 xlab = "Time", ylab = "Target 2")
# points(x = lst_targets$Target2$time, 
#        y = v_out_best$Target2, 
#        pch = 8, col = "red")
# legend("topright", 
#        legend = c("Target", "Model-predicted output"),
#        col = c("black", "red"), pch = c(1, 8))

## Fancier pairwise plot ----
gg_post_pairs_corr <- GGally::ggpairs(data.frame(m_calib_res[, v_param_names]),
                                      upper = list(continuous = wrap("cor",
                                                                     color = "black",
                                                                     size = 5)),
                                      diag = list(continuous = wrap("barDiag",
                                                                    alpha = 0.8)),
                                      lower = list(continuous = wrap("points", 
                                                                     alpha = 0.3,
                                                                     size = 0.7)),
                                      columnLabels = v_param_names
) +
  theme_bw(base_size = 18) +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_text(size = 14),
        axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        strip.background = element_rect(fill = "white",
                                        color = "white"),
        strip.text = element_text(hjust = 0))
gg_post_pairs_corr

## Prior vs. posterior ----
m_samp_prior <- sample.prior(n_resamp)
df_samp_prior <- reshape2::melt(cbind(PDF = "Prior", 
                                      as.data.frame(m_samp_prior)), 
                                variable.name = "Parameter")
df_samp_post_imis <- reshape2::melt(cbind(PDF = "Posterior IMIS",
                                          as.data.frame(m_calib_res[, v_param_names])),
                                    variable.name = "Parameter")
df_samp_prior_post <- dplyr::bind_rows(df_samp_prior, 
                                       df_samp_post_imis)

gg_prior_post_imis <- ggplot(df_samp_prior_post, 
                             aes(x = value, y = ..density.., fill = PDF)) +
  facet_wrap(~Parameter, scales = "free", 
             ncol = 4,
             labeller = label_parsed) +
  scale_x_continuous(n.breaks = 6) +
  geom_density(alpha = 0.5) +
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
gg_prior_post_imis

#  Propagate calibrated parameter uncertainty  ----
## Compute IMIS posterior predicted outputs ----
m_out_surv     <- matrix(NA, 
                         nrow = n_resamp, 
                         ncol = length(lst_targets$Surv$value))  

### Run model for each posterior parameter set ----
for (i in 1:n_resamp) { # i = 1
  model_res_temp <- run_crs_markov(m_calib_res[i, ])
  m_out_surv[i, ] <- model_res_temp$Surv
  if (i/100 == round(i/100,0)) { 
    cat('\r', paste(i/n_resamp*100, "% done", sep = ""))
  }
}

## Posterior-predicted mean ----
v_out_surv_post_mean <- colMeans(m_out_surv)

## Posterior-predicted 95% credible interval ----
m_out_surv_95cri <- colQuantiles(m_out_surv, probs = c(0.025, 0.975))

df_out_post <- data.frame(Type = "Model output",
                          dplyr::bind_cols(Outcome = "Survival", 
                                           time = lst_targets[[1]]$time,
                                           value = v_out_surv_post_mean,
                                           lb = m_out_surv_95cri[, 1],
                                           ub = m_out_surv_95cri[, 2]))
df_out_post$Outcome <- ordered(df_out_post$Outcome, 
                               levels = c("Survival"))

## Plot targets vs. model-predicted output ----
df_targets <- data.frame(cbind(Type = "Target", 
                               Outcome = "Survival", 
                               lst_targets[[1]]))
df_targets$Outcome <- ordered(df_targets$Outcome, 
                              levels = c("Survival"))
ggplot(df_targets, aes(x = time, y = value, 
                       ymin = lb, ymax = ub)) +
  # geom_point(shape = 1, size = 2) +
  geom_errorbar() +
  geom_line(data = df_out_post, 
            aes(x = time, y = value), col = "blue") +
  geom_ribbon(data = df_out_post,
              aes(ymin = lb, ymax = ub), alpha = 0.4, fill = "blue") +
  facet_wrap(~ Outcome, scales = "free_x") +
  scale_x_continuous("Time", n.breaks = 8) +
  scale_y_continuous("Proportion", n.breaks = 8) +
  theme_bw(base_size = 16)
