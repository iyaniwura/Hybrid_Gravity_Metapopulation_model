library(tidyverse)
library(rstan)
library(shinystan)
library(gridExtra)
library("bayesplot")
library("tidybayes")
library(matrixStats)
library("loo")
library(lubridate)
library(reshape)
library(ggridges)

rstan_options (auto_write = TRUE)   # telling stan not to recompile code that has already been compiled
options (mc.cores = parallel::detectCores() )
set.seed(3)


path = "/home/iyaniwura/Dropbox/BCCDC/MetaPopulationModel_ByNoticeRinga/Manuscript_Files/GithubCodes/DistanceMatrix/"
setwd(path)
source("Fit_Summary.R")   # Loads the function that store the outputs of the samplying in a folder
# the outputs stored by this include: the density plot, trace plot, prediction of cases  plotted 
# together with the actual cases data for each age-group. 


# Model Parameters: PredictedCases_ByRegion_Random_IC.csv

# Defining model parameters
n_region = 13
N_vec = c(161912, 10770, 257926, 105862, 112259, 8931, 161725, 47652, 111502, 82590, 107347, 512436, 254583)  # population of the regions arranged into a vector


incubation_period = 5 # days
Presymptomatic_period = 1 # days
infectious_period = 5 # days
h1 = (incubation_period/7)^(-1) # in weeks
h2 = (Presymptomatic_period/7)^(-1) # in weeks
gamma = ((infectious_period/2)/7)^(-1)
sigma = ((infectious_period/2)/7)^(-1)

# theta = 0.5
# beta_vec = c(1.50, 1.30, 1.20, 1.80, 1.42, 1.51, 1.94, 1.49, 1.79, 1.19, 1.61, 1.55, 1.90)

# reading in time-series mobility data
# TS_Mob_path = "/home/iyaniwura/Dropbox/BCCDC/MetaPopulationModel_ByNoticeRinga/Mobility_Data/TimeSeriesMobility/NonResidents_Only/Normalized_by_LHADevices/"
# TS_Mob_RawData = read.csv(paste(TS_Mob_path, "NonR_TimeSeriesMob_Data.csv", sep="") ) 
TS_Mob_path = "/home/iyaniwura/Dropbox/BCCDC/MetaPopulationModel_ByNoticeRinga/Mobility_Data/Scaled_Unscaled_TS_Mob_data/Scaled/"
TS_Mob_RawData = read.csv(paste(TS_Mob_path, "Scaled_RNonR_TimeSeriesMob_Data.csv", sep="") ) 
TS_Mob_Extract = unname(as.matrix(TS_Mob_RawData[,3:15]))
TS_Mob_Vec = as.vector(t(TS_Mob_Extract))



# computing the number of infected in each region before our study period
Before_data_path = "/home/iyaniwura/Dropbox/BCCDC/MetaPopulationModel_ByNoticeRinga/CasesData/RCodes_DataProcessing/"
Before_Cases_data = read.csv(paste(Before_data_path, "WeeklyCases_Jan2020_June2020.csv", sep="") )
BeforeCasesMat <- unname(as.matrix(Before_Cases_data[,3:15]))
InfectedB4ByLHA = colSums(BeforeCasesMat)
Sus_Pop_LHA = N_vec - InfectedB4ByLHA 

# constructing the distance matrix
LHA_Ditsance <- read.csv("/home/iyaniwura/Dropbox/BCCDC/MetaPopulationModel_ByNoticeRinga/Manuscript_Files/DistanceMatrix/New_Distance_File/Unscaled_DistanceMatrix.csv")  
class(LHA_Ditsance)
D_Mat = unname(as.matrix(LHA_Ditsance))

# constructing the probability of moving from one LHA to another
f_ij = 1 / (rep(1,(n_region*n_region)) + D_Mat); 
pi_ij = f_ij/rowSums(f_ij)
pi_ij_Vec = as.vector(t(pi_ij)) 


# preparing the initial condition vector
Init_infected = BeforeCasesMat[nrow(BeforeCasesMat),]
X0 = rep(0, 7*n_region)  # initializing the initial population vector
X0[1:n_region] = Sus_Pop_LHA  # initial susceptible population
# X0[(3*n_region + 1):(4*n_region)] = Init_infected #rep(2,n_region)
X0[(1*n_region + 1):(5*n_region)] = Init_infected 
sum(X0[(3*n_region + 1):(5*n_region)])

# calculating the initial prevalence fraction
temp = X0[(n_region + 1):(5*n_region)]
state_0 = temp/sum(temp)
sum(state_0)

# Organizing the cases data by region
data_path = "/home/iyaniwura/Dropbox/BCCDC/MetaPopulationModel_ByNoticeRinga/CasesData/RCodes_DataProcessing/"
Raw_Cases = read.csv(paste(data_path, "WeeklyCases_July2020_Jan2021.csv", sep="") )
ReportCases_ByRegion <- Raw_Cases   #[1:20,]



Dates = Raw_Cases$Date
WeeklyCase_Reg_1 <- ReportCases_ByRegion$Abbotsford
WeeklyCase_Reg_2 <-   ReportCases_ByRegion$Agassiz.Harrison
WeeklyCase_Reg_3 <-  ReportCases_ByRegion$Burnaby
WeeklyCase_Reg_4 <- ReportCases_ByRegion$Chilliwack
WeeklyCase_Reg_5 <- ReportCases_ByRegion$Delta
WeeklyCase_Reg_6 <- ReportCases_ByRegion$Hope
WeeklyCase_Reg_7 <- ReportCases_ByRegion$Langley
WeeklyCase_Reg_8 <- ReportCases_ByRegion$Mission
WeeklyCase_Reg_9 <- ReportCases_ByRegion$Maple.Ridge.Pitt.Meadows
WeeklyCase_Reg_10 <-ReportCases_ByRegion$New.Westminster
WeeklyCase_Reg_11 <-ReportCases_ByRegion$South.Surrey.White.Rock
WeeklyCase_Reg_12 <-ReportCases_ByRegion$Surrey
WeeklyCase_Reg_13 <-ReportCases_ByRegion$Tri.Cities



# Preparing the parameters for stan
t0 = 0								# Initial time
n_weeks = nrow(ReportCases_ByRegion)     # number of weeks
ts = 1:n_weeks

data_seir = list(
	n_region = n_region,
	N_vec = N_vec,
	h1 = h1,
	h2 = h2,
	gamma = gamma,
	sigma = sigma,
	pi_ij_Vec = pi_ij_Vec,
	Init_Suscept = Sus_Pop_LHA,  # 
	t0 = t0,
	n_weeks = n_weeks,
	ts = ts,
	inference = 1,
	doprint = 0,
	TS_Mob_Vec = TS_Mob_Vec,
	# Data to fit
	Rep_Reg_1 = WeeklyCase_Reg_1,
	Rep_Reg_2 = WeeklyCase_Reg_2,
	Rep_Reg_3 = WeeklyCase_Reg_3,
	Rep_Reg_4 = WeeklyCase_Reg_4,
	Rep_Reg_5 = WeeklyCase_Reg_5,
	Rep_Reg_6 = WeeklyCase_Reg_6,
	Rep_Reg_7 = WeeklyCase_Reg_7,
	Rep_Reg_8 = WeeklyCase_Reg_8,
	Rep_Reg_9 = WeeklyCase_Reg_9,
	Rep_Reg_10 = WeeklyCase_Reg_10,
	Rep_Reg_11 = WeeklyCase_Reg_11,
	Rep_Reg_12 = WeeklyCase_Reg_12,
	Rep_Reg_13 = WeeklyCase_Reg_13,
	
	# Priors   Note: I have not used them in the computation yet, just created place holders
	# p_i0 = c(log(n_region*2), 1),
	p_i0 = c(log(150), 1),
	p_c0Mean = c(0, 1),
	p_c0Var = c(0.5, 1),
	
	p_dMean = c(0, 1),
	p_dVar = c(0.5, 1),
	
	## priors for theta
	p_theta = c(0.5, 0.1),
	p_phi = c(5),                                # parameters for the prior for phi 
	
	
	# initial condition
	y0_vars = state_0,
	p_underreport_cases=1
)



# initial condition function for sampling
inti_cond <- function(data_seir){
	i0_prior = data_seir$p_i0
	c0Mean_prior = data_seir$p_c0Mean
	c0Var_prior = data_seir$p_c0Var
	
	dMean_prior = data_seir$p_dMean
	dVar_prior = data_seir$p_dVar
	
	## priors for theta
	theta_prior = data_seir$p_theta
	
	
	i0 <- stats::rlnorm(1, i0_prior[1], i0_prior[2]/2)
	c0Mean <- stats::rnorm(1, c0Mean_prior[1], c0Mean_prior[2]/2)
	c0Var <- stats::rnorm(1, c0Var_prior[1], c0Var_prior[2]/2)
	
	p_dMean <- stats::rnorm(1, dMean_prior[1], dMean_prior[2]/2)
	p_dVar <- stats::rnorm(1, dVar_prior[1], dVar_prior[2]/2)
	
	## priors for theta
	theta_1 <- stats::rnorm(1, theta_prior[1], theta_prior[2]/2)
	
	Sample_IC <- list(i0 = i0,  c0Mean = c0Mean, c0Var = c0Var, theta_1 = theta_1,
										p_dMean= p_dMean, p_dVar = p_dVar	)
	Sample_IC
}



AgeModel_Full <- stan_model(paste(path, "SEIR_DistMat.stan", sep="") )
# fit_seir <- sampling(AgeModel_Full,data = data_seir, iter = 5, init = function() inti_cond(data_seir),  chains = 4)
# fit_seir <- sampling(AgeModel_Full, data = data_seir, iter = 5,  chains = 1)

vb_fit_seir  <-  vb(AgeModel_Full, data = data_seir, init = function() inti_cond(data_seir), algorithm = "meanfield", iter = 30000 )
 fit_seir <- vb_fit_seir

 # vb_fit_seir  <-  vb(AgeModel_Full, data = data_seir, algorithm = "meanfield", iter = 30000 )
# fit_seir <- vb_fit_seir


########################## saving stan object #########################
saveRDS(fit_seir, file="Fit_DistMat_Output.RDS")
fit_seir <- readRDS("Fit_DistMat_Output.RDS")



################################## Computing WAIC and loo ##########################################################
log_lik_1 <- extract_log_lik(vb_fit_seir, merge_chains = FALSE)
r_eff <- relative_eff(exp(log_lik_1), cores = 2)
loo_1 <- loo(log_lik_1, r_eff = r_eff, cores = 2)
print(loo_1)
plot(loo_1)

# use the stanfit object directly
loo_exp2 <- loo(vb_fit_seir, pars = "log_lik")
loo_exp2
plot(loo_exp2)

# which is equivalent to
LL <- as.array(vb_fit_seir, pars = "log_lik")
r_eff <- loo::relative_eff(exp(LL))
loo1b <- loo::loo.array(LL, r_eff = r_eff)
print(loo1b)
plot(loo1b)

# WAIC....Note: No applicable method for 'waic' applied to an object of class "stanfit"
WAIC_1 <- waic(log_lik_1)
print(WAIC_1)
plot(WAIC_1)
############################################################################################



# Extracting the posteriors  from the output of the sampling and saving it in a file
Extract <- rstan::extract(fit_seir)
Post  <- data.frame(i0 = Extract$i0, c0_PopMean= Extract$c0_PopMean, c0_PopVar= Extract$c0_PopVar,
										c0_1 = Extract$c0_1, c0_2 = Extract$c0_2, c0_3 = Extract$c0_3, c0_4 = Extract$c0_4, c0_5 = Extract$c0_5,
										c0_6 = Extract$c0_6, c0_7 = Extract$c0_7, c0_8 = Extract$c0_8, c0_9 = Extract$c0_9, c0_10 = Extract$c0_10,
										c0_11 = Extract$c0_11, c0_12 = Extract$c0_12, c0_13 = Extract$c0_13,
										theta_1 = Extract$theta_1, c1 = Extract$c1, d_PopMean = Extract$d_PopMean, d_PopVar = Extract$d_PopVar,
										d1 = Extract$d1, d2 = Extract$d2, d3 = Extract$d3, d4 = Extract$d4, d5 = Extract$d5, d6 = Extract$d6,
										d7 = Extract$d7, phi = Extract$phi)

write.table(Post, file="Sampling_DistMat_Output.csv")


##################################################### displaying fit summary ##########################################################
# checking inference
pars=c( "i0", "c0_PopMean", "c0_PopVar", "c0_1", "c0_2", "c0_3", "c0_4", "c0_5", "c0_6", "c0_7", "c0_8", "c0_9", "c0_10", "c0_11", "c0_12", "c0_13", 'theta_1',
        "c1", "d_PopMean", "d_PopVar", "d1", "d2", "d3", "d4", "d5", "d6", "d7", 'phi') # specifying parameters of interest
print(fit_seir, pars = pars)
stan_dens(fit_seir, pars = pars, separate_chains = TRUE)
traceplot(fit_seir, pars = pars)
pairs(fit_seir, pars = pars)


# plots and saves fit summary in a folder called Fit_Summary
names(ReportCases_ByRegion)[3:15] = c("Reg_1", "Reg_2", "Reg_3", "Reg_4", "Reg_5", "Reg_6",
																			"Reg_7", "Reg_8", "Reg_9", "Reg_10", "Reg_11", "Reg_12", "Reg_13")
Fit_Summary(fit_seir, pars, path, ReportCases_ByRegion)  


########################################################################################################
Post <- read.table(file= paste(path, "Sampling_DistMat_Output.csv", sep="") )
# Post = Post[1:5,]

# Projection codes 
mainDir <- paste(path, "ProjectionFunctions/", sep="")
source(paste(mainDir, "Required_Functions.R", sep=""))
source(paste(mainDir, "SEIR_ODE_function.R", sep=""))


#### stopped here #################
# loading libraries
library(deSolve)
# library.dynam.unload("deSolve", libpath=paste(.libPaths()[1], "//deSolve", sep=""))
# library.dynam("deSolve", package="deSolve", lib.loc=.libPaths()[1])

# forecasting parameters
forecast_wks = 0
tfinal = n_weeks + forecast_wks
tstart = t0            
dt = 1
input_pars = data.frame(n_region=n_region, h1=h1, h2=h2, gamma=gamma, sigma=sigma,
												start_time = tstart, final_time = tfinal, dt = dt)



############################### calling the function that solve ODE system  ############################### 
# calling the function that computed the projection
Proj_Output <- Proj_ODE_Solver(Post, input_pars, Sus_Pop_LHA, state_0, pi_ij_Vec, TS_Mob_Vec)
saveRDS(Proj_Output, file="Proj_DistMat.RDS")


# calling the function that plots the projection
ODE_Proj_Output <- readRDS("Proj_DistMat.RDS")
Proj_ByAge(tfinal, dt, ODE_Proj_Output, path)



# Plotting the beta for each region
Beta_Processor(Post, n_region)













































# 
# #############################################################################
# library(ggridges)
# 
# library(ggplot2)
# 
# d <- data.frame(y = rep(1:5, 3), x = c(rep(0, 5), rep(1, 5), rep(3, 5)),
# 								width = c(0, 1, 3, 4, 0, 1, 2, 3, 5, 4, 0, 5, 4, 4, 1))
# ggplot(d, aes(x, y, width = width, group = x)) + geom_vridgeline(fill="lightblue")
# 
# ggplot(iris, aes(x=Species, y=Sepal.Width, width = ..density.., fill=Species)) +
# 	geom_vridgeline(stat="ydensity", trim=FALSE, alpha = 0.65, scale = 1)
# 
# 
# #############################################################################
# library(tidyverse)
# library(ggridges)
# 
# dfs <-
# 	data.frame(
# 		"estimate" = rnorm(300),
# 		"loading" = factor(rep(1:5, 60)),
# 		"set" = factor(rep(1:3, 100),),
# 		"pop" = rep(c(0.1, 0.15, 0.20, 0.05, 0.7), 60)
# 	)
# 
# ggplot(dfs, aes(x = estimate, y = loading, fill = set)) +
# 	geom_density_ridges(
# 		jittered_points = FALSE,
# 		point_shape = "|", point_size = 2, size = 0.25,
# 		position = position_points_jitter(height = 0), alpha = 0.6, scale = 1.2) 
# #> Picking joint bandwidth of 0.395
# 
# 
# 
# 
# ggplot(dj_longData, aes(x = value, y = variable, fill = Scenario)) +
# 	geom_density_ridges(
# 		jittered_points = FALSE,
# 		point_shape = "|", point_size = 2, size = 0.25,
# 		position = position_points_jitter(height = 0), alpha = 0.6, scale = 1.2) 
# #> Picking joint bandwidth of 0.395
# 
# 
# 
# 
# 
# 
# 
