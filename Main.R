# Author: Roberto Zuniga
# Date: November 12, 2021
# Description: This script partially replicates Blunder, Gowrisankaran and Blundell
# (2020). Specifically, It solver the Bellman equation and computes the nested
# fixed point estimator.


# Set working directories
current_dir <- "./IO_partII"
setwd(current_dir)

outputs_dir <- paste0(current_dir,"Outputs")
bellman_data_file <- paste0(current_dir,"Data/data_for_bellman_computations.csv")
panel_data_file <- paste0(current_dir,"Data/analysis_data.csv")

source("Bellman_functions.R")

# Read Bellman data
data_bellman <- read.csv(bellman_data_file)
data <- subset(data_bellman,naics_recode == 1)
data$index <- 1:nrow(data)

# Read analysis data for lilkelihood computations
panel_data <- read.csv(panel_data_file)
panel_data <- subset(panel_data,naics_recode == 1)

# Define the working parameters
params <- list()
params$data <- data
params$panel_data <- panel_data 
params$N1 <- 161  # number of states per each omega1
params$n_omega1 <- length(unique((data$omega1))) # number of omega1 states
params$DAVgrid <- (0:19)/19*9.5  # DAV grid
params$n_DAV <- length(params$DAVgrid) # number of points in DAV grid
params$DAV_deprate <- 0.1  # DAV depreciation rate
params$beta <- 0.95^(1/4)  # discount factor
params$gamma <- 0.57721566490153*0 # Euler's gamma constant

params$tol <- 1E-6 # stopping condition for the Bellman equation solution
params$nregoutcomes <- 80 # number of regulatory outcomes
params$ntransitions <- 3  # number of states transitions from omega to omega_tilde
################################################################################
# Create structures for fast computation of Bellman equation
params <- AddParams(params)

################################################################################
# Solve the Bellman equation
coeff_BGL <- c(2.872, -0.049, -0.077, -5.980, -0.065) # BGL
resultsBGL <- Bellman(coeff_BGL, params)

# Check resuts
ExpectedResult1 <- paste0(current_dir,"Data/value_part2_problem1_thetaBGL.csv")
CompareValues(resultsBGL,ExpectedResult1,params)

ExpectedResult2 <- paste0(current_dir,"Data/valuetilde_part2_problem1_thetaBGL.csv")
CompareVTilde(resultsBGL,ExpectedResult2,params)

################################################################################
# Compute Loglikelihood
LogLikelihoodBGL <- LogLike(resultsBGL,params)

# Check resuts
ExpectedResult3 <- paste0(current_dir,"Data/likelihood_part2_problem2_thetaBGL.csv")
compTable<-CompareProb(LogLikelihoodBGL$lldata,ExpectedResult3)

################################################################################
# Prob 1.c: Export probabilities for DAV = 2
coeff_1C <- c(2,-0.5,-0.5,-5,-0.1)
results1C <- Bellman(coeff_1C, params)
LogLike1C <- LogLike(results1C,params)
filename1C <- "Output/prob_DAV2_Q1C.csv"
Question1C(results1C,filename1C)
################################################################################
# Prob 2.a: Report quasilikelihood
LogLike1C$ll

################################################################################
# Maximize quasi-likelihood
theta0 <- c(2.872, -0.049, -0.077, -5.980, -0.065)
NFXP <- NestedFixedPoint(theta0,params,3)

# Check that log likelihood is higher at thetaML
theta0 <- c(2.872, -0.049, -0.077, -5.980, -0.065)
loglike_0 <- EvalFunction(theta0,params)
thetaML <- c(2.847, -1.084, -1.751, -3.890, 0.312)
loglike_ML <- EvalFunction(thetaML,params)

################################################################################
# Estimate variance at theta ML
thetaML <- c(2.847, -1.084, -1.751, -3.890, 0.312)
DeltaTheta <- 1e-6
varML <- EstimateMLVariance(thetaML,params,DeltaTheta)
seML <- sqrt(diag(varML))
