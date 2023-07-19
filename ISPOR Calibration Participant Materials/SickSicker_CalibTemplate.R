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
###########################################################################


###################  Calibration Specifications  ###################

# Model: 4-State Sick-Sicker Markov Model
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




####################################################################
######  Load model as a function  ######
####################################################################
# - inputs are parameters to be estimated through calibration
# - outputs correspond to the target data




####################################################################
######  Specify calibration parameters  ######
####################################################################





####################################################################
######  Calibrate!  ######
####################################################################





####################################################################
######  Exploring best-fitting input sets  ######
####################################################################



# Note: When visualizing the top 100 (top 10%) input sets, use scatterplot3d() instead of plot().
