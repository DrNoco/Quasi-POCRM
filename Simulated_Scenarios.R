################################################################
############# Toxicity Weight matrix Scenarios #################
################################################################
a
ntox <- 3 #Define number of toxicities
W <- array(NA, c(ntox, 4, 3)) #Define Array to contain 3 different toxicity weight matrices 

# Define weight matrices (rows represent toxicities, columns represent grades)

#Weight Matrix 1 with toxic combinations 
W[,,1] <- matrix(c(0.1, 0.35, 0.7, 1.00,
              0.08, 0.23,  0.6, 0.80,
              0.00, 0.15, 0.45, 0.80), nrow=ntox, byrow=T)


#Weight Matrix 2 with combinations not too toxic 
W[,,2] <- matrix(c(0.05, 0.20, 0.7, 1.00,
                   0.05, 0.15,  0.6, 0.80,
                   0.00, 0.15, 0.45, 0.80), nrow=ntox, byrow=T)

#Weight Matrix 3: When non-DLT toxicities don't contribute 
W[,,3] <- matrix(c(0.0, 0, 0.7, 1,
              0.0, 0, 0.6, 0.8,
              0.0, 0, 0, 0.8), nrow=ntox, byrow=T)

################################################################
####  Adjust W such that all possible scores are in [0,1] ######
#### 'tox.weights' function defined in functions script   ######
################################################################
for(i in 1:ntox){
  W[,,i] <- tox.weights(W[,,i]) 
}

#################################################
#######  Define the 3x2 scenarios ###############
#################################################

d <- 6 #define number of dose levels
nps <-6 # number unique scenarios to be defined 
TOX <- array(NA, c(d, 5, ntox, nps)) #Define Array to hold toxicity scenarios

# In the array TOX[a,b,c,d]: rows 'a' represents the dose-level in increasing order, 
# columns 'b represent AE grade (from 0 to 4), position 'c' represents the AE type,
# and position 'd' represents the scenario as a whole. Each scenario (d) has 3 
# AE types (c) considered. The position (a,b) within scenario 'd' 
# represents the probability of a grade 'b' AE of type 'c' at dose level 'a'. 

## Scenario 1: a single best dose-level ##
TOX[,,1,1] <- matrix(c(0.791, 0.172, 0.032, 0.004, 0.001, 
                     0.738, 0.195, 0.043, 0.015, 0.009, 
                     0.555, 0.240, 0.148, 0.044, 0.013,
                     0.602, 0.230, 0.108, 0.046, 0.014,
                     0.305, 0.323, 0.231, 0.091, 0.050,
                     0.200, 0.307, 0.261, 0.154, 0.078), nrow=d, byrow=T)

TOX[,,2,1] <- matrix(c(0.968, 0.029, 0.002, 0.001, 0.000,
                     0.763, 0.192, 0.006, 0.039, 0.000,
                     0.662, 0.233, 0.091, 0.010, 0.004,
                     0.552, 0.255, 0.178, 0.010, 0.005,
                     0.397, 0.258, 0.276, 0.061, 0.008,
                     0.160, 0.377, 0.281, 0.133, 0.049), nrow=d, byrow=T)

TOX[,,3,1] <- matrix(c(0.907, 0.080, 0.008, 0.000, 0.005, 
                     0.552, 0.331, 0.035, 0.040, 0.042,
                     0.306, 0.308, 0.161, 0.091, 0.134,
                     0.015, 0.134, 0.240, 0.335, 0.276, 
                     0.005, 0.052, 0.224, 0.372, 0.347,
                     0.004, 0.022, 0.223, 0.345, 0.406), nrow=d, byrow=T)

## Check to see that the probabilities within a row (dose-level) sum to 1
#apply(TOX[,,1,1], 1, sum)
#apply(TOX[,,2,1], 1, sum)
#apply(TOX[,,3,1], 1, sum)

#Scenario 2 - #Where dose levels 3 and 4 are the same and best
TOX[,,1,2] <- matrix(c(0.791, 0.172, 0.032, 0.004, 0.001, 
                       0.558, 0.245, 0.173, 0.015, 0.009,  
                     0.472, 0.280, 0.188, 0.036, 0.024,
                     0.452, 0.250, 0.178, 0.086, 0.034,
                     0.334, 0.283, 0.232, 0.131, 0.020,
                     0.250, 0.337, 0.251, 0.114, 0.048), nrow=6, byrow=T)

TOX[,,2,2] <- matrix(c(0.868, 0.129, 0.002, 0.001, 0.000,
                       0.603, 0.232, 0.126, 0.039, 0.000,
                     0.542, 0.305, 0.108, 0.040, 0.005,
                     0.552, 0.245, 0.138, 0.050, 0.015,
                     0.396, 0.258, 0.277, 0.061, 0.008,
                     0.260, 0.377, 0.281, 0.073, 0.009), nrow=6, byrow=T)

TOX[,,3,2] <- matrix(c(0.718, 0.220, 0.057, 0.000, 0.005, 
                       0.362, 0.301, 0.135, 0.100, 0.102,
                     0.125, 0.134, 0.210, 0.275, 0.256, 
                     0.045, 0.154, 0.320, 0.295, 0.186, 
                     0.055, 0.102, 0.224, 0.372, 0.247,
                     0.004, 0.023, 0.220, 0.344, 0.409), nrow=6, byrow=T)


#Scenario 3: Dose levels 3 and 5 are the same and best
TOX[,,1,3] <- matrix(c(0.791, 0.172, 0.032, 0.004, 0.001, 
                       0.738, 0.195, 0.043, 0.015, 0.009, 
                     0.552, 0.260, 0.128, 0.046, 0.014,
                     0.655, 0.220, 0.068, 0.044, 0.013,
                     0.562, 0.250, 0.128, 0.046, 0.014,
                     0.390, 0.307, 0.202, 0.073, 0.028), nrow=6, byrow=T)

TOX[,,2,3] <- matrix(c(0.968, 0.029, 0.002, 0.001, 0.000,
                       0.763, 0.192, 0.006, 0.039, 0.000,
                     0.672, 0.205, 0.108, 0.010, 0.005,
                     0.662, 0.233, 0.091, 0.010, 0.004,
                     0.572, 0.305, 0.108, 0.010, 0.005,
                     0.260, 0.377, 0.281, 0.074, 0.008), nrow=6, byrow=T)

TOX[,,3,3] <- matrix(c(0.918, 0.070, 0.007, 0.000, 0.005, 
                       0.602, 0.281, 0.035, 0.040, 0.042,
                     0.015, 0.134, 0.240, 0.335, 0.276,
                     0.296, 0.309, 0.130, 0.131, 0.134,
                     0.015, 0.134, 0.240, 0.335, 0.276,
                     0.004, 0.022, 0.220, 0.395, 0.359), nrow=6, byrow=T)


#Scenario 4: Dose levels 6 is best
TOX[,,1,4] <- matrix(c(0.791, 0.172, 0.032, 0.004, 0.001, 
                       0.738, 0.195, 0.043, 0.015, 0.009, 
                     0.685, 0.190, 0.068, 0.044, 0.013,
                     0.710, 0.200, 0.068, 0.012, 0.010,
                     0.605, 0.230, 0.108, 0.044, 0.013,
                     0.602, 0.200, 0.128, 0.050, 0.020), nrow=6, byrow=T)

TOX[,,2,4] <- matrix(c(0.968, 0.029, 0.002, 0.001, 0.000,
                       0.763, 0.192, 0.006, 0.039, 0.000,
                     0.762, 0.183, 0.041, 0.010, 0.004,
                     0.792, 0.163, 0.035, 0.010, 0.000,
                     0.502, 0.283, 0.161, 0.050, 0.004,
                     0.453, 0.335, 0.147, 0.050, 0.015), nrow=6, byrow=T)


TOX[,,3,4] <- matrix(c(0.917, 0.070, 0.008, 0.000, 0.005, 
                       0.602, 0.281, 0.035, 0.040, 0.042,
                     0.536, 0.209, 0.030, 0.091, 0.134,
                     0.596, 0.219, 0.030, 0.091, 0.064,
                     0.286, 0.308, 0.131, 0.161, 0.114,
                     0.065, 0.134, 0.240, 0.335, 0.226), nrow=6, byrow=T)

#Scenario 5: Dose Level 1 is best
TOX[,,1,5] <- matrix(c(0.662, 0.200, 0.077, 0.046, 0.015, 
                       0.605, 0.223, 0.081, 0.071, 0.020, 
                       0.390, 0.308, 0.201, 0.073, 0.028,
                       0.492, 0.230, 0.169, 0.085, 0.024,
                       0.240, 0.307, 0.301, 0.124, 0.028,
                       0.124, 0.236, 0.469, 0.100, 0.071), nrow=6, byrow=T)

TOX[,,2,5] <- matrix(c(0.582, 0.245, 0.158, 0.010, 0.005,
                       0.397, 0.258, 0.276, 0.061, 0.008,
                       0.260, 0.378, 0.281, 0.073, 0.008,
                       0.352, 0.315, 0.269, 0.050, 0.014,
                       0.260, 0.358, 0.281, 0.073, 0.028,
                       0.134, 0.306, 0.409, 0.100, 0.051), nrow=6, byrow=T)

TOX[,,3,5] <- matrix(c(0.015, 0.134, 0.240, 0.335, 0.276, 
                       0.005, 0.052, 0.224, 0.372, 0.347,
                       0.004, 0.022, 0.221, 0.344, 0.409,
                       0.015, 0.134, 0.240, 0.335, 0.276, 
                       0.084, 0.082, 0.220, 0.284, 0.330,
                       0.004, 0.026, 0.150, 0.400, 0.420), nrow=6, byrow=T)


#Scenario 6 - #Where dose levels 3 and 4 are the same and best (used in sensitivity analysis of Appendix)
TOX[,,1,6] <- matrix(c(0.791, 0.172, 0.032, 0.004, 0.001, 
                       0.608, 0.225, 0.123, 0.035, 0.009,  
                       0.482, 0.280, 0.178, 0.036, 0.024,
                       0.448, 0.295, 0.223, 0.025, 0.009,  
                       0.402, 0.240, 0.238, 0.076, 0.044,
                       0.362, 0.300, 0.278, 0.046, 0.014), nrow=6, byrow=T)

TOX[,,2,6] <- matrix(c(0.868, 0.129, 0.002, 0.001, 0.000,
                       0.603, 0.232, 0.126, 0.039, 0.000,
                       0.542, 0.305, 0.108, 0.040, 0.005,
                       0.503, 0.282, 0.186, 0.029, 0.000,
                       0.422, 0.285, 0.228, 0.050, 0.015,
                       0.302, 0.405, 0.278, 0.010, 0.005), nrow=6, byrow=T)

TOX[,,3,6] <- matrix(c(0.718, 0.220, 0.057, 0.000, 0.005, 
                       0.502, 0.151, 0.085, 0.100, 0.162,
                       0.125, 0.134, 0.210, 0.275, 0.256, 
                       0.202, 0.291, 0.285, 0.110, 0.112,
                       0.045, 0.154, 0.320, 0.295, 0.186,
                       0.004, 0.022, 0.220, 0.395, 0.359), nrow=6, byrow=T)


##############################################################################################
### Assess DLT probabilities and Toxicity Scores for each scenario across dose-levels ########
### doses increase left to right across columns and up rows                           ########

b <- 6 # define which scenario we wish to assess
a <- 1 # define which weight matrix to use for Tox-score calculation
tox.dlt <- c(3,3,4) #define the AE grade that denotes a DLT for each AE type 1-3

# Assess DLT probabilities across dose-levels 
ptox <- round(dlt.prob(TOX[,,,b], ntox, tox.dlt),2)
matrix(c(ptox[4:6], ptox[1:3]),nrow=2,byrow=T)

# Assess the true mean toxicity score based on weights for each dose level
p0 <- round(as.vector(t(get.thresh(ntox, W[,,a], TOX[,,,b]))),2)
matrix(c(p0[4:6], p0[1:3]),nrow=2,byrow=T)

# Assess the true mean TI score at each dose level
ti_scores <- round(as.vector(t(mean.tiscore(ntox, TOX[,,,b], tox.dlt))),3)
matrix(c(ti_scores[4:6], ti_scores[1:3]),nrow=2,byrow=T)


############################################################################
##### Specify the number of possible toxicity orderings/simple-orders ######
s<-5 

### Specify the possible toxicity orderings of the drug combinations
orders<-matrix(nrow=s,ncol=d)
orders[1,]<-c(1,2,3,4,5,6)  
orders[2,]<-c(1,2,4,5,3,6)  
orders[3,]<-c(1,2,4,3,5,6)  
orders[4,]<-c(1,4,2,3,5,6)  
orders[5,]<-c(1,4,2,5,3,6) 

### Specify a set of toxicity skeleton values for Tox Scores and DLT outcomes
init.guess <- 3 #initial guess of position for skeletons
hw <- 0.03 #halfwidth for skeletons of toxicity scores based on weights
hw2 <- 0.04 #halfwidth for skeletons for DLT probabilities
hw3 <- 0.04 #haldwidth for skeletons of TI index scores

### Specify targets for toxicity scores and DLT probability
target <-0.20 # target toxicity score
target2 <- 0.33 # target DLT probability
target3 <- .47 #target TI score

### Get the matrix of model skeletons based on the above specified criteria
p.skel <- get.skeletons(target,  hw, init.guess, d,s, orders)
p.skel2 <- get.skeletons(target2, hw2, init.guess, d,s, orders)
p.skel3 <- get.skeletons(target3, hw3, init.guess, d,s, orders)

#################################################################
######### Run Simulation for a specified scenario ###############
#################################################################

cohortsize=1 ## cohort size for each inclusion
ncohort=36 ## number of cohorts
start.comb=1 ## starting combination
ntrial= 2000 ## number of simulated trials 
seed.no = 11235 ## Set Seed

b <- 1 ## Define which scenario (1 to 6) to test
a <- 1 ## Define which toxicity weight matrix to use (1-3)

## The vector denoting the minimum grade for which each toxicity is a DLT based on specified target scores
tox.dlt <- find.dlt(W[,,a], target)

#The true probability of a DLT for each dose level
ptox <- round(dlt.prob(TOX[,,,b], ntox, tox.dlt),2)

#The true mean toxicity score for each dose level
p0 <- round(as.vector(t(get.thresh(ntox, W[,,a], TOX[,,,b]))),3)

# The TI score 
tiscores <- mean.tiscore(ntox, TOX[,,,b], tox.dlt)

#run Sim WITHOUT overdose control (gives results for QPOCRM-W, QPOCRM-TI, and POCRM)
results <- qpocrm.sim(TOX[,,,b], p.skel, p.skel2, p.skel3, target, target2, target3, tox.dlt, p0, ptox, tiscores, W[,,a], cohortsize,ncohort,ntrial,start.comb, ncores=6)
results

#Run Sim WITH Overdose Control (this runs ONLY the OC model: specify wtype="TI" or wtype="weights")
results <- qpocrm.sim.oc(TOX[,,,b], p.skel3, p.skel2,  target3, target2,  tox.dlt, p0, ptox, tiscores, wtype="TI", W[,,a], cohortsize,ncohort,ntrial,start.comb, ncores=6)
results

#Run Sim WITH Overdose Control w/ reversed order assumption for the overdose constraint (this runs ONLY the OC model: specify wtype="TI" or wtype="weights")
results <- qpocrm.sim.oc2(TOX[,,,b], p.skel3, p.skel2,  target3, target2,  tox.dlt, p0, ptox, tiscores, wtype="TI", W[,,a], cohortsize,ncohort,ntrial,start.comb, ncores=6)
results


##########################################################################################
########################### 4X4 scenarios ################################################
##########################################################################################

#Matrices of probabilities for each toxicity of each grade(columns, for grades 0-4) at each dose level (rows, 6 dose levels)
d <- 16
nps <-2 #number of probability scnenarios to run at once
TOX <- array(NA, c(d, 5, ntox, nps)) #Define Array to hold toxicity weights

#Scenario 1: 16 dose levels for a 4x4 matching nolans pocrm paper example 1
TOX[,,1,1] <- matrix(c(       0.775, 0.172, 0.032, 0.015, 0.006, 
                            0.749, 0.180, 0.043, 0.019, 0.009, 
                            0.662, 0.225, 0.078, 0.025, 0.010,
                            0.612, 0.257, 0.078, 0.042, 0.011,
                            0.555, 0.257, 0.153, 0.025, 0.010, #5
                            0.510, 0.337, 0.111, 0.030, 0.012, #6
                            0.355, 0.298, 0.235, 0.100, 0.012, #7
                            0.305, 0.289, 0.225, 0.130, 0.051, #8
                            0.512, 0.307, 0.128, 0.042, 0.011, #9
                            0.345, 0.308, 0.235, 0.100, 0.012, #10
                            0.305, 0.284, 0.205, 0.150, 0.056, #11
                            0.298, 0.254, 0.185, 0.160, 0.103, #12
                            0.355, 0.234, 0.205, 0.150, 0.056, #13
                            0.285, 0.276, 0.205, 0.160, 0.074, #14
                            0.305, 0.267, 0.165, 0.160, 0.103, #15
                            0.130, 0.204, 0.335, 0.180, 0.151), nrow=d, byrow=T)


TOX[,,2,1] <- matrix(c(       0.775, 0.172, 0.032, 0.015, 0.006, 
                            0.749, 0.180, 0.043, 0.019, 0.009, 
                            0.662, 0.225, 0.078, 0.025, 0.010,
                            0.612, 0.257, 0.078, 0.042, 0.011,
                            0.555, 0.257, 0.153, 0.025, 0.010, #5
                            0.430, 0.337, 0.191, 0.030, 0.012, #6
                            0.305, 0.318, 0.265, 0.100, 0.012, #7
                            0.295, 0.295, 0.229, 0.130, 0.051, #8
                            0.512, 0.307, 0.128, 0.042, 0.011, #9
                            0.345, 0.308, 0.235, 0.100, 0.012, #10
                            0.305, 0.284, 0.205, 0.150, 0.056, #11
                            0.208, 0.254, 0.275, 0.160, 0.103, #12
                            0.275, 0.234, 0.285, 0.150, 0.056, #13
                            0.163, 0.266, 0.331, 0.160, 0.080, #14
                            0.305, 0.267, 0.165, 0.160, 0.103, #15
                            0.130, 0.204, 0.335, 0.180, 0.151), nrow=d, byrow=T)

TOX[,,3,1] <- matrix(c(       0.770, 0.162, 0.032, 0.015, 0.021, 
                            0.730, 0.180, 0.043, 0.019, 0.028, 
                            0.622, 0.220, 0.078, 0.045, 0.035,
                            0.550, 0.227, 0.098, 0.072, 0.053,
                            0.480, 0.277, 0.183, 0.025, 0.035, #5
                            0.450, 0.337, 0.121, 0.050, 0.042, #6
                            0.320, 0.233, 0.205, 0.130, 0.112, #7
                            0.135, 0.229, 0.295, 0.160, 0.181, #8
                            0.550, 0.227, 0.098, 0.072, 0.053, #9
                            0.290, 0.263, 0.235, 0.100, 0.112, #10
                            0.195, 0.204, 0.245, 0.150, 0.206, #11
                            0.058, 0.164, 0.295, 0.220, 0.263, #12
                            0.205, 0.204, 0.195, 0.190, 0.206, #13
                            0.100, 0.181, 0.295, 0.190, 0.234, #14
                            0.030, 0.107, 0.235, 0.365, 0.263, #15
                            0.020, 0.010, 0.239, 0.400, 0.331), nrow=d, byrow=T)


#Scenario 2: 16 dose levels matching POCRM example 2

TOX[,,1,2] <- matrix(c(     0.769, 0.180, 0.023, 0.019, 0.009,
                            0.680, 0.237, 0.041, 0.032, 0.010, 
                            0.592, 0.277, 0.068, 0.052, 0.011,
                            0.382, 0.347, 0.188, 0.052, 0.031,
                            0.655, 0.207, 0.103, 0.025, 0.010, #5
                            0.512, 0.277, 0.158, 0.042, 0.011, #6
                            0.355, 0.328, 0.205, 0.100, 0.012, #7
                            0.305, 0.289, 0.225, 0.130, 0.051, #8
                            0.510, 0.337, 0.111, 0.030, 0.012, #9
                            0.345, 0.308, 0.235, 0.100, 0.012, #10
                            0.305, 0.284, 0.205, 0.150, 0.056, #11
                            0.245, 0.286, 0.235, 0.160, 0.074, #12
                            0.345, 0.308, 0.235, 0.100, 0.012, #13
                            0.305, 0.284, 0.205, 0.150, 0.056, #14
                            0.245, 0.276, 0.245, 0.160, 0.074, #15
                            0.105, 0.327, 0.305, 0.160, 0.103), nrow=d, byrow=T)


TOX[,,2,2] <- matrix(c(     0.749, 0.190, 0.033, 0.019, 0.009, 
                            0.630, 0.267, 0.061, 0.030, 0.012, 
                            0.612, 0.257, 0.078, 0.042, 0.011,
                            0.402, 0.337, 0.208, 0.042, 0.011,
                            0.555, 0.257, 0.153, 0.025, 0.010, #5
                            0.412, 0.337, 0.198, 0.042, 0.011, #6
                            0.365, 0.288, 0.235, 0.100, 0.012, #7
                            0.255, 0.295, 0.269, 0.130, 0.051, #8
                            0.430, 0.337, 0.191, 0.030, 0.012, #9
                            0.345, 0.308, 0.235, 0.100, 0.012, #10
                            0.305, 0.284, 0.205, 0.150, 0.056, #11
                            0.163, 0.266, 0.331, 0.160, 0.080, #12
                            0.345, 0.308, 0.235, 0.100, 0.012, #13
                            0.305, 0.284, 0.205, 0.150, 0.056, #14
                            0.163, 0.316, 0.281, 0.160, 0.080, #15
                            0.105, 0.267, 0.365, 0.160, 0.103), nrow=d, byrow=T)




TOX[,,3,2] <- matrix(c(     0.730, 0.180, 0.043, 0.019, 0.028, 
                            0.600, 0.257, 0.051, 0.050, 0.042, 
                            0.550, 0.227, 0.098, 0.072, 0.053,
                            0.360, 0.337, 0.178, 0.072, 0.053,
                            0.520, 0.297, 0.123, 0.025, 0.035, #5
                            0.400, 0.317, 0.158, 0.072, 0.053, #6
                            0.270, 0.273, 0.215, 0.130, 0.112, #7
                            0.185, 0.219, 0.255, 0.160, 0.181, #8
                            0.450, 0.337, 0.121, 0.050, 0.042, #9
                            0.290, 0.263, 0.235, 0.100, 0.112, #10
                            0.195, 0.204, 0.245, 0.150, 0.206, #11
                            0.100, 0.181, 0.295, 0.190, 0.234, #12
                            0.290, 0.263, 0.235, 0.100, 0.112, #13
                            0.195, 0.204, 0.245, 0.150, 0.206, #14
                            0.100, 0.181, 0.295, 0.190, 0.234, #15
                            0.030, 0.107, 0.235, 0.365, 0.263), nrow=d, byrow=T)

# Check every row sums to 1 
apply(TOX[,,1,1], 1, sum)
apply(TOX[,,2,1], 1, sum)
apply(TOX[,,3,1], 1, sum)
apply(TOX[,,1,2], 1, sum)
apply(TOX[,,2,2], 1, sum)
apply(TOX[,,3,2], 1, sum)

##############################################################################################
### Assess DLT probabilities and Toxicity Scores for each scenario across dose-levels ########
### doses increase left to right across columns and up rows                           ########

b <- 2 # define which scenario we wish to assess
a <- 1 # define which weight matrix to use for Tox-score calculation
tox.dlt <- c(3,3,4) #define the AE grade that denotes a DLT for each AE type 1-3

# Assess DLT probabilities across dose-levels 
ptox <- round(dlt.prob(TOX[,,,b], ntox, tox.dlt),2)
matrix(c(ptox[13:16], ptox[9:12], ptox[5:8], ptox[1:4]),nrow=4,byrow=T)

# Assess the true mean toxicity score based on weights for each dose level
p0 <- round(as.vector(t(get.thresh(ntox, W[,,a], TOX[,,,b]))),2)
matrix(c(p0[13:16], p0[9:12], p0[5:8], p0[1:4]),nrow=4,byrow=T)

# Assess the true mean TI score at each dose level
ti_scores <- round(as.vector(t(mean.tiscore(ntox, TOX[,,,b], tox.dlt))),3)
matrix(c(ti_scores[13:16], ti_scores[9:12], ti_scores[5:8], ti_scores[1:4]),nrow=4,byrow=T)


##############################################################
#####Specify the number of possible toxicity orderings
#Orders for 4x4 trial
d=16
s=3
orders<-matrix(nrow=s,ncol=d)
orders[1,]<-c(1,2,5,3,6,9,4,7,10,13,8,11,14,12,15,16)
#orders[1,]<-c(1,2,5,9,6,3,7,4,10,13,14,11,8,12,15,16)
orders[2,]<-c(1,5,2,3,6,9,10,13,7,4,8,11,14,15,12,16)
orders[3,]<-c(1,5,2,9,6,3,13,10,7,4,14,11,8,15,12,16)

### Specify a set of toxicity skeleton values for Tox Scores and DLT outcomes
init.guess <- 8 #ininital guess of position for skeletons
hw <- .03 # halfwidth for toxicity scores (weights) skeleton
hw2 <- .04 # halfwidth for DLT probability skeleton
hw3 <- .04 #halfwidth for TI score skeleton

### Specify targets for toxicity scores and DLT probability
target <-0.20 # target toxicity score
target2 <- 0.30 # target DLT probability
target3 <- .47 #target TI score

### Get the matrix of model skeletons based on the above specified criteria
p.skel <- get.skeletons(target,  hw, init.guess, d,s, orders)
p.skel2 <- get.skeletons(target2, hw2, init.guess, d,s, orders)
p.skel3 <- get.skeletons(target3, hw3, init.guess, d,s, orders)

#################################################################
######### Run Simulation for a specified scenario ###############
#################################################################

cohortsize=1 ## cohort size for each inclusion
ncohort=36 ## number of cohorts
start.comb=1 ## starting combination
ntrial= 2000 ## number of simulated trials 
seed.no = 11235 ## Set Seed

b <- 1 ## Define which scenario (1 to 2) to test
a <- 1 ## Define which toxicity weight matrix to use (1-3)

## The vector denoting the minimum grade for which each toxicity is a DLT based on specified target scores
tox.dlt <- find.dlt(W[,,a], target)

#The true probability of a DLT for each dose level
ptox <- round(dlt.prob(TOX[,,,b], ntox, tox.dlt),2)

#The true mean toxicity score for each dose level
p0 <- round(as.vector(t(get.thresh(ntox, W[,,a], TOX[,,,b]))),3)

# The TI score 
tiscores <- mean.tiscore(ntox, TOX[,,,b], tox.dlt)

#run Sim WITHOUT overdose control (gives results for QPOCRM-W, QPOCRM-TI, and POCRM)
results <- qpocrm.sim(TOX[,,,b], p.skel, p.skel2, p.skel3, target, target2, target3, tox.dlt, p0, ptox, tiscores, W[,,a], cohortsize,ncohort,ntrial,start.comb, ncores=6)
results

#Run Sim WITH Overdose Control (this runs ONLY the OC model: specify wtype="TI" or wtype="weights")
results <- qpocrm.sim.oc(TOX[,,,b], p.skel3, p.skel2,  target3, target2,  tox.dlt, p0, ptox, tiscores, wtype="TI", W[,,a], cohortsize,ncohort,ntrial,start.comb, ncores=6)
results

#Run Sim WITH Overdose Control w/ reversed order assumption for the overdose constraint (this runs ONLY the OC model: specify wtype="TI" or wtype="weights")
results <- qpocrm.sim.oc2(TOX[,,,b], p.skel3, p.skel2,  target3, target2,  tox.dlt, p0, ptox, tiscores, wtype="TI", W[,,a], cohortsize,ncohort,ntrial,start.comb, ncores=6)
results


##########################################################
########################### 3X5 scenarios ###############
############################################################
#Matrices of probabilities for each toxicity of each grade(columns, for grades 0-4) at each dose level (rows, 6 dose levels)
d <- 15
nps <-2 #number of probability scnenarios to run at once
TOX <- array(NA, c(d, 5, ntox, nps)) #Define Array to hold toxicity weights

#Scenario 3: 15 dose levels for a 3x5 matching  pocrm paper example 3
TOX[,,1,1] <- matrix(c(     0.775, 0.172, 0.032, 0.015, 0.006, 
                            0.719, 0.200, 0.053, 0.019, 0.009, 
                            0.662, 0.225, 0.078, 0.025, 0.010,
                            0.355, 0.298, 0.235, 0.100, 0.012,
                            0.305, 0.284, 0.205, 0.150, 0.056, #5
                            
                            0.510, 0.337, 0.111, 0.030, 0.012, #6
                            0.512, 0.247, 0.178, 0.052, 0.011, #7
                            0.355, 0.298, 0.235, 0.100, 0.012, #8
                            0.305, 0.284, 0.205, 0.150, 0.056, #9
                            0.285, 0.276, 0.205, 0.160, 0.074,  #10
                            
                            0.512, 0.257, 0.178, 0.042, 0.011, #11
                            0.355, 0.298, 0.235, 0.100, 0.012, #12
                            0.305, 0.284, 0.205, 0.150, 0.056, #13
                            0.285, 0.276, 0.205, 0.160, 0.074, #14
                            0.255, 0.267, 0.215, 0.160, 0.103), nrow=d, byrow=T)




TOX[,,2,1] <- matrix(c(     0.775, 0.172, 0.032, 0.015, 0.006, 
                            0.749, 0.180, 0.043, 0.019, 0.009, 
                            0.662, 0.225, 0.078, 0.025, 0.010,
                            0.305, 0.318, 0.265, 0.100, 0.012,
                            0.305, 0.284, 0.205, 0.150, 0.056, #5
                            
                            0.430, 0.337, 0.191, 0.030, 0.012, #6
                            0.412, 0.357, 0.178, 0.042, 0.011, #7
                            0.335, 0.298, 0.255, 0.100, 0.012, #8
                            0.305, 0.284, 0.205, 0.150, 0.056, #9
                            0.163, 0.266, 0.331, 0.160, 0.080, #10
                            
                            0.462, 0.307, 0.178, 0.042, 0.011, #11
                            0.305, 0.318, 0.265, 0.100, 0.012, #12
                            0.305, 0.284, 0.205, 0.150, 0.056, #13
                            0.203, 0.266, 0.291, 0.160, 0.080, #14
                            0.155, 0.267, 0.315, 0.160, 0.103), nrow=d, byrow=T)


TOX[,,3,1] <- matrix(c(     0.770, 0.162, 0.032, 0.015, 0.021, 
                            0.730, 0.180, 0.043, 0.019, 0.028, 
                            0.622, 0.220, 0.078, 0.045, 0.035,
                            0.330, 0.233, 0.195, 0.130, 0.112,
                            0.195, 0.204, 0.245, 0.150, 0.206, #5
                            
                            0.450, 0.337, 0.121, 0.050, 0.042, #6
                            0.550, 0.227, 0.098, 0.072, 0.053, #7
                            0.290, 0.233, 0.235, 0.130, 0.112, #8
                            0.195, 0.204, 0.245, 0.150, 0.206, #9
                            0.100, 0.181, 0.295, 0.190, 0.234, #10
                            
                            0.350, 0.327, 0.198, 0.072, 0.053, #11
                            0.250, 0.303, 0.205, 0.130, 0.112, #12
                            0.195, 0.204, 0.245, 0.150, 0.206, #13
                            0.100, 0.181, 0.295, 0.190, 0.234, #14
                            0.030, 0.107, 0.235, 0.365, 0.263), nrow=d, byrow=T)


#Scenario 4: 15 dose levels for a 3x5 matching pocrm paper example 4
TOX[,,1,2] <- matrix(c(       0.775, 0.172, 0.032, 0.015, 0.006, 
                            0.662, 0.225, 0.078, 0.025, 0.010, 
                            0.612, 0.257, 0.078, 0.042, 0.011,
                            0.355, 0.298, 0.235, 0.100, 0.012,
                            0.305, 0.284, 0.205, 0.150, 0.056, #5
                            
                            0.662, 0.225, 0.078, 0.025, 0.010, #6
                            0.355, 0.298, 0.235, 0.100, 0.012, #7
                            0.305, 0.284, 0.205, 0.150, 0.056, #8
                            0.130, 0.204, 0.335, 0.180, 0.151, #9
                            0.070, 0.174, 0.325, 0.250, 0.181,  #10
                            
                            0.305, 0.284, 0.205, 0.150, 0.056,  #11
                            0.305, 0.267, 0.165, 0.160, 0.103, #12
                            0.130, 0.204, 0.335, 0.180, 0.151, #13
                            0.070, 0.174, 0.325, 0.250, 0.181, #14
                            0.020, 0.154, 0.285, 0.320, 0.221), nrow=d, byrow=T)


TOX[,,2,2] <- matrix(c(      0.775, 0.172, 0.032, 0.015, 0.006, 
                            0.662, 0.225, 0.078, 0.025, 0.010, 
                            0.612, 0.257, 0.078, 0.042, 0.011,
                            0.305, 0.318, 0.265, 0.100, 0.012,
                            0.305, 0.284, 0.205, 0.150, 0.056, #5
                            
                            0.662, 0.225, 0.078, 0.025, 0.010, #6
                            0.305, 0.318, 0.265, 0.100, 0.012, #7
                            0.305, 0.284, 0.205, 0.150, 0.056, #8
                            0.130, 0.204, 0.335, 0.180, 0.151, #9
                            0.090, 0.204, 0.305, 0.220, 0.181, #10
                            
                            0.305, 0.284, 0.205, 0.150, 0.056, #11
                            0.305, 0.267, 0.165, 0.160, 0.103, #12
                            0.130, 0.204, 0.335, 0.180, 0.151, #13
                            0.090, 0.204, 0.305, 0.220, 0.181, #14
                            0.020, 0.124, 0.285, 0.310, 0.261), nrow=d, byrow=T)


TOX[,,3,2] <- matrix(c(      0.770, 0.162, 0.032, 0.015, 0.021, 
                            0.622, 0.220, 0.078, 0.045, 0.035, 
                            0.550, 0.227, 0.098, 0.072, 0.053,
                            0.340, 0.233, 0.185, 0.130, 0.112,
                            0.195, 0.204, 0.245, 0.150, 0.206, #5
                            
                            0.622, 0.220, 0.078, 0.045, 0.035, #6
                            0.320, 0.253, 0.215, 0.100, 0.112, #7
                            0.195, 0.204, 0.245, 0.150, 0.206, #8
                            0.020, 0.010, 0.239, 0.400, 0.331, #9
                            0.020, 0.080, 0.189, 0.310, 0.401, #10
                            
                            0.195, 0.204, 0.245, 0.150, 0.206, #11
                            0.030, 0.107, 0.235, 0.365, 0.263, #12
                            0.020, 0.010, 0.239, 0.400, 0.331, #13
                            0.020, 0.080, 0.189, 0.310, 0.401, #14
                            0.010, 0.040, 0.149, 0.330, 0.471), nrow=d, byrow=T)

# Check every row sums to 1 
apply(TOX[,,1,1], 1, sum)
apply(TOX[,,2,1], 1, sum)
apply(TOX[,,3,1], 1, sum)
apply(TOX[,,1,2], 1, sum)
apply(TOX[,,2,2], 1, sum)
apply(TOX[,,3,2], 1, sum)

##############################################################################################
### Assess DLT probabilities and Toxicity Scores for each scenario across dose-levels ########
### doses increase left to right across columns and up rows                           ########

b <- 2 # define which scenario we wish to assess
a <- 1 # define which weight matrix to use for Tox-score calculation
tox.dlt <- c(3,3,4) #define the AE grade that denotes a DLT for each AE type 1-3

# Assess DLT probabilities across dose-levels 
ptox <- round(dlt.prob(TOX[,,,b], ntox, tox.dlt),2)
matrix(c(ptox[11:15], ptox[6:10], ptox[1:5]),nrow=3,byrow=T)

# Assess the true mean toxicity score based on weights for each dose level
p0 <- round(as.vector(t(get.thresh(ntox, W[,,a], TOX[,,,b]))),2)
matrix(c(p0[11:15], p0[6:10], p0[1:5]),nrow=3,byrow=T)

# Assess the true mean TI score at each dose level
ti_scores <- round(as.vector(t(mean.tiscore(ntox, TOX[,,,b], tox.dlt))),3)
matrix(c(ti_scores[11:15], ti_scores[6:10], ti_scores[1:5]),nrow=3,byrow=T)


##############################################################
#####Specify the number of possible toxicity orderings
#Orders for 3x5 trial
d=15
s=3
orders<-matrix(nrow=s,ncol=d)
orders[1,]<-c(1,2,6,3,7,11,4,8,12,5,9,13,10,14,15)
orders[2,]<-c(1,6,2,11,7,3,12,8,4,13,9,5,14,10,15)
orders[3,]<-c(1,6,2,7,11,3,8,12,4,9,13,5,14,10,15)

### Specify a set of toxicity skeleton values for Tox Scores and DLT outcomes
init.guess <- 8 #ininital guess of position for skeletons
hw <- .03 # halfwidth for toxicity scores (weights) skeleton
hw2 <- .04 # halfwidth for DLT probability skeleton
hw3 <- .04 #halfwidth for TI score skeleton

### Specify targets for toxicity scores and DLT probability
target <-0.20 # target toxicity score
target2 <- 0.30 # target DLT probability
target3 <- .47 #target TI score

### Get the matrix of model skeletons based on the above specified criteria
p.skel <- get.skeletons(target,  hw, init.guess, d,s, orders)
p.skel2 <- get.skeletons(target2, hw2, init.guess, d,s, orders)
p.skel3 <- get.skeletons(target3, hw3, init.guess, d,s, orders)



#################################################################
######### Run Simulation for a specified scenario ###############
#################################################################

cohortsize=1 ## cohort size for each inclusion
ncohort=36 ## number of cohorts
start.comb=1 ## starting combination
ntrial= 2000 ## number of simulated trials 
seed.no = 11235 ## Set Seed

b <- 1 ## Define which scenario (1 to 2) to test
a <- 1 ## Define which toxicity weight matrix to use (1-3)

## The vector denoting the minimum grade for which each toxicity is a DLT based on specified target scores
tox.dlt <- find.dlt(W[,,a], target)

#The true probability of a DLT for each dose level
ptox <- round(dlt.prob(TOX[,,,b], ntox, tox.dlt),2)

#The true mean toxicity score for each dose level
p0 <- round(as.vector(t(get.thresh(ntox, W[,,a], TOX[,,,b]))),3)

# The TI score 
tiscores <- mean.tiscore(ntox, TOX[,,,b], tox.dlt)

#run Sim WITHOUT overdose control (gives results for QPOCRM-W, QPOCRM-TI, and POCRM)
results <- qpocrm.sim(TOX[,,,b], p.skel, p.skel2, p.skel3, target, target2, target3, tox.dlt, p0, ptox, tiscores, W[,,a], cohortsize,ncohort,ntrial,start.comb, ncores=6)
results

#Run Sim WITH Overdose Control (this runs ONLY the OC model: specify wtype="TI" or wtype="weights")
results <- qpocrm.sim.oc(TOX[,,,b], p.skel3, p.skel2,  target3, target2,  tox.dlt, p0, ptox, tiscores, wtype="TI", W[,,a], cohortsize,ncohort,ntrial,start.comb, ncores=6)
results

#Run Sim WITH Overdose Control w/ reversed order assumption for the overdose constraint (this runs ONLY the OC model: specify wtype="TI" or wtype="weights")
results <- qpocrm.sim.oc2(TOX[,,,b], p.skel3, p.skel2,  target3, target2,  tox.dlt, p0, ptox, tiscores, wtype="TI", W[,,a], cohortsize,ncohort,ntrial,start.comb, ncores=6)
results


########################################################################################
########################### 4X5 scenarios ##############################################
########################################################################################
#Matrices of probabilities for each toxicity of each grade(columns, for grades 0-4) at each dose level (rows, 6 dose levels)
d <- 20
nps <-2 #number of probability scnenarios to run at once
TOX <- array(NA, c(d, 5, ntox, nps)) #Define Array to hold toxicity weights

#Scenario 5: 20 dose levels for a 4x5 matching  pocrm paper example 5
TOX[,,1,1] <- matrix(c(     0.495, 0.268, 0.165, 0.052, 0.020, 
                            0.380, 0.303, 0.167, 0.110, 0.040, 
                            0.305, 0.274, 0.215, 0.160, 0.046,
                            0.255, 0.286, 0.205, 0.170, 0.084,
                            0.195, 0.267, 0.235, 0.180, 0.123, #5
                            
                            0.370, 0.299, 0.175, 0.112, 0.039, #6
                            0.335, 0.284, 0.195, 0.130, 0.056, #7
                            0.315, 0.257, 0.185, 0.153, 0.090, #8
                            0.205, 0.207, 0.235, 0.230, 0.123, #9
                            0.120, 0.184, 0.335, 0.200, 0.161,  #10
                            
                            0.305, 0.314, 0.195, 0.140, 0.046, #11
                            0.225, 0.307, 0.205, 0.180, 0.083, #12
                            0.155, 0.217, 0.265, 0.260, 0.103, #13
                            0.110, 0.164, 0.335, 0.240, 0.151, #14
                            0.050, 0.134, 0.365, 0.260, 0.191, #15
                            
                            0.305, 0.246, 0.195, 0.170, 0.084, #16
                            0.205, 0.217, 0.225, 0.230, 0.123, #17
                            0.160, 0.164, 0.285, 0.220, 0.171, #18
                            0.080, 0.134, 0.335, 0.260, 0.191, #19
                            0.050, 0.124, 0.365, 0.280, 0.181), nrow=d, byrow=T)


TOX[,,2,1] <- matrix(c(       0.405, 0.298, 0.185, 0.090, 0.022, 
                              0.334, 0.295, 0.180, 0.146, 0.045, 
                              0.265, 0.294, 0.215, 0.160, 0.066,
                              0.203, 0.266, 0.281, 0.170, 0.080,
                              0.195, 0.267, 0.245, 0.180, 0.113, #5
                              
                              0.325, 0.265, 0.189, 0.180, 0.041, #6
                              0.265, 0.240, 0.265, 0.160, 0.070, #7
                              0.235, 0.237, 0.265, 0.160, 0.103, #8
                              0.205, 0.267, 0.225, 0.170, 0.133, #9
                              0.120, 0.184, 0.315, 0.200, 0.181, #10
                              
                              0.305, 0.274, 0.185, 0.174, 0.062, #11
                              0.285, 0.277, 0.185, 0.160, 0.093, #12
                              0.195, 0.267, 0.225, 0.200, 0.113, #13
                              0.120, 0.184, 0.315, 0.220, 0.161, #14
                              0.090, 0.194, 0.305, 0.230, 0.181, #15
                              
                              0.163, 0.266, 0.341, 0.170, 0.060, #16
                              0.205, 0.267, 0.225, 0.190, 0.113, #17
                              0.120, 0.184, 0.315, 0.220, 0.161, #18
                              0.090, 0.174, 0.305, 0.250, 0.181, #19
                              0.030, 0.104, 0.355, 0.310, 0.201), nrow=d, byrow=T)


TOX[,,3,1] <- matrix(c(       0.400, 0.265, 0.155, 0.100, 0.080, 
                              0.315, 0.244, 0.175, 0.145, 0.121, 
                              0.275, 0.204, 0.195, 0.160, 0.166,
                              0.150, 0.201, 0.245, 0.190, 0.214,
                              0.070, 0.107, 0.195, 0.365, 0.263, #5
                              
                              0.285, 0.309, 0.185, 0.130, 0.091, #6
                              0.205, 0.264, 0.225, 0.180, 0.126, #7
                              0.120, 0.147, 0.215, 0.255, 0.263, #8
                              0.080, 0.107, 0.235, 0.315, 0.263, #9
                              0.030, 0.010, 0.249, 0.380, 0.331, #10
                              
                              0.204, 0.205, 0.235, 0.200, 0.156, #11
                              0.170, 0.127, 0.175, 0.265, 0.263, #12
                              0.060, 0.107, 0.235, 0.335, 0.263, #13
                              0.040, 0.010, 0.239, 0.380, 0.331, #14
                              0.020, 0.080, 0.189, 0.310, 0.401, #15
                              
                              0.100, 0.201, 0.295, 0.210, 0.194, #16
                              0.070, 0.107, 0.235, 0.325, 0.263, #17
                              0.050, 0.010, 0.239, 0.370, 0.331, #18
                              0.030, 0.080, 0.179, 0.310, 0.401, #19
                              0.020, 0.080, 0.149, 0.300, 0.451), nrow=d, byrow=T)


#Scenario 6: 20 dose levels for a 4x5 matching pocrm paper example 4
TOX[,,1,2] <- matrix(c(     0.582, 0.217, 0.118, 0.062, 0.021, 
                            0.555, 0.229, 0.105, 0.085, 0.026, 
                            0.405, 0.239, 0.205, 0.121, 0.030, 
                            0.325, 0.254, 0.235, 0.140, 0.046,
                            0.295, 0.226, 0.225, 0.170, 0.084, #5
                            
                            0.505, 0.258, 0.155, 0.070, 0.012, #6
                            0.435, 0.229, 0.185, 0.122, 0.029, #7
                            0.345, 0.244, 0.225, 0.140, 0.046, #8
                            0.285, 0.307, 0.165, 0.160, 0.103, #9
                            0.155, 0.227, 0.265, 0.230, 0.123,  #10
                            
                            0.444, 0.230, 0.175, 0.130, 0.031, #11
                            0.355, 0.244, 0.205, 0.150, 0.046, #12
                            0.305, 0.257, 0.165, 0.180, 0.093, #13
                            0.185, 0.217, 0.235, 0.240, 0.123, #14
                            0.110, 0.164, 0.335, 0.220, 0.171, #15
                            
                            0.325, 0.284, 0.205, 0.140, 0.046, #16
                            0.305, 0.267, 0.165, 0.160, 0.103, #17
                            0.180, 0.204, 0.245, 0.220, 0.151, #18
                            0.155, 0.197, 0.285, 0.240, 0.123, #19
                            0.050, 0.134, 0.365, 0.260, 0.191), nrow=d, byrow=T)


TOX[,,2,2] <- matrix(c(       0.512, 0.357, 0.078, 0.042, 0.011, 
                              0.455, 0.265, 0.149, 0.111, 0.020, 
                              0.375, 0.275, 0.159, 0.161, 0.030, 
                              0.305, 0.264, 0.205, 0.170, 0.056,
                              0.153, 0.266, 0.331, 0.170, 0.080, #5
                              
                              0.425, 0.268, 0.195, 0.110, 0.002, #6
                              0.375, 0.245, 0.159, 0.180, 0.041, #7
                              0.295, 0.254, 0.215, 0.170, 0.066, #8
                              0.255, 0.277, 0.215, 0.160, 0.093, #9
                              0.205, 0.267, 0.225, 0.190, 0.113, #10
                              
                              0.424, 0.205, 0.170, 0.170, 0.031, #11
                              0.345, 0.224, 0.185, 0.190, 0.056, #12
                              0.305, 0.277, 0.165, 0.160, 0.093, #13
                              0.195, 0.267, 0.225, 0.200, 0.113, #14
                              0.120, 0.184, 0.315, 0.220, 0.161, #15
                              
                              0.325, 0.234, 0.205, 0.180, 0.056, #16
                              0.285, 0.277, 0.165, 0.180, 0.093, #17
                              0.195, 0.307, 0.205, 0.180, 0.113, #18
                              0.120, 0.184, 0.305, 0.230, 0.161, #19
                              0.090, 0.194, 0.305, 0.230, 0.181), nrow=d, byrow=T)

TOX[,,3,2] <- matrix(c(       0.550, 0.227, 0.098, 0.072, 0.053, 
                              0.425, 0.229, 0.165, 0.100, 0.081,  
                              0.355, 0.199, 0.175, 0.150, 0.121, 
                              0.215, 0.204, 0.205, 0.220, 0.156,
                              0.100, 0.201, 0.295, 0.190, 0.214, #5
                              
                              0.410, 0.233, 0.175, 0.100, 0.082, #6
                              0.335, 0.269, 0.185, 0.120, 0.091, #7
                              0.275, 0.244, 0.215, 0.150, 0.116, #8
                              0.160, 0.107, 0.235, 0.235, 0.263, #9
                              0.030, 0.107, 0.235, 0.365, 0.263, #10
                              
                              0.320, 0.264, 0.165, 0.140, 0.111, #11
                              0.245, 0.254, 0.185, 0.160, 0.156, #12
                              0.150, 0.107, 0.185, 0.305, 0.253, #13
                              0.100, 0.107, 0.185, 0.345, 0.263, #14
                              0.020, 0.010, 0.239, 0.400, 0.331, #15
                              
                              0.265, 0.194, 0.185, 0.180, 0.176, #16
                              0.120, 0.107, 0.215, 0.335, 0.223, #17
                              0.060, 0.077, 0.245, 0.345, 0.273, #18
                              0.030, 0.010, 0.219, 0.390, 0.351, #19
                              0.020, 0.080, 0.189, 0.310, 0.401), nrow=d, byrow=T)

# Check every row sums to 1 
apply(TOX[,,1,1], 1, sum)
apply(TOX[,,2,1], 1, sum)
apply(TOX[,,3,1], 1, sum)
apply(TOX[,,1,2], 1, sum)
apply(TOX[,,2,2], 1, sum)
apply(TOX[,,3,2], 1, sum)

##############################################################################################
### Assess DLT probabilities and Toxicity Scores for each scenario across dose-levels ########
### doses increase left to right across columns and up rows                           ########

b <- 2 # define which scenario we wish to assess
a <- 1 # define which weight matrix to use for Tox-score calculation
tox.dlt <- c(3,3,4) #define the AE grade that denotes a DLT for each AE type 1-3

# Assess DLT probabilities across dose-levels 
ptox <- round(dlt.prob(TOX[,,,b], ntox, tox.dlt),2)
matrix(c(ptox[16:20], ptox[11:15], ptox[6:10], ptox[1:5]),nrow=4,byrow=T)

# Assess the true mean toxicity score based on weights for each dose level
p0 <- round(as.vector(t(get.thresh(ntox, W[,,a], TOX[,,,b]))),2)
matrix(c(p0[16:20], p0[11:15], p0[6:10], p0[1:5]),nrow=4,byrow=T)

# Assess the true mean TI score at each dose level
tiscores <- round(as.vector(t(mean.tiscore(ntox, TOX[,,,b], tox.dlt))),3)
matrix(c(tiscores[16:20], tiscores[11:15], tiscores[6:10], tiscores[1:5]),nrow=4,byrow=T)



##############################################################
#####Specify the number of possible toxicity orderings
#Orders for 4x5 trial
d=20
s=3
orders<-matrix(nrow=s,ncol=d)
orders[1,]<-c(1,2,6,3,7,11,4,8,12,16,5,9,13,17,10,14,18,15,19,20)
#orders[1,]<-c(1,2,6,7,3,11,4,8,12,16,5,9,13,17,10,14,18,15,19,20)
orders[2,]<-c(1,6,2,11,7,3,16,12,8,4,17,13,9,5,18,14,10,19,15,20)
orders[3,]<-c(1,6,2,7,11,3,12,8,16,4,13,9,17,5,14,18,10,19,15,20)

### Specify a set of toxicity skeleton values for Tox Scores and DLT outcomes
init.guess <- 8 #ininital guess of position for skeletons
hw <- .03 # halfwidth for toxicity scores (weights) skeleton
hw2 <- .04 # halfwidth for DLT probability skeleton
hw3 <- .04 #halfwidth for TI score skeleton

### Specify targets for toxicity scores and DLT probability
target <-0.21 # target toxicity score
target2 <- 0.40 # target DLT probability
target3 <- .49 #target TI score

### Get the matrix of model skeletons based on the above specified criteria
p.skel <- get.skeletons(target,  hw, init.guess, d,s, orders)
p.skel2 <- get.skeletons(target2, hw2, init.guess, d,s, orders)
p.skel3 <- get.skeletons(target3, hw3, init.guess, d,s, orders)

#################################################################
######### Run Simulation for a specified scenario ###############
#################################################################

cohortsize=1 ## cohort size for each inclusion
ncohort=36 ## number of cohorts
start.comb=1 ## starting combination
ntrial= 2000 ## number of simulated trials 
seed.no = 11235 ## Set Seed

b <- 1 ## Define which scenario (1 to 2) to test
a <- 1 ## Define which toxicity weight matrix to use (1-3)

## The vector denoting the minimum grade for which each toxicity is a DLT based on specified target scores
tox.dlt <- find.dlt(W[,,a], target)

#The true probability of a DLT for each dose level
ptox <- round(dlt.prob(TOX[,,,b], ntox, tox.dlt),2)

#The true mean toxicity score for each dose level
p0 <- round(as.vector(t(get.thresh(ntox, W[,,a], TOX[,,,b]))),3)

# The TI score 
tiscores <- mean.tiscore(ntox, TOX[,,,b], tox.dlt)

#run Sim WITHOUT overdose control (gives results for QPOCRM-W, QPOCRM-TI, and POCRM)
results <- qpocrm.sim(TOX[,,,b], p.skel, p.skel2, p.skel3, target, target2, target3, tox.dlt, p0, ptox, tiscores, W[,,a], cohortsize,ncohort,ntrial,start.comb, ncores=6)
results

#Run Sim WITH Overdose Control (this runs ONLY the OC model: specify wtype="TI" or wtype="weights")
results <- qpocrm.sim.oc(TOX[,,,b], p.skel3, p.skel2,  target3, target2,  tox.dlt, p0, ptox, tiscores, wtype="TI", W[,,a], cohortsize,ncohort,ntrial,start.comb, ncores=6)
results

#Run Sim WITH Overdose Control w/ reversed order assumption for the overdose constraint (this runs ONLY the OC model: specify wtype="TI" or wtype="weights")
results <- qpocrm.sim.oc2(TOX[,,,b], p.skel3, p.skel2,  target3, target2,  tox.dlt, p0, ptox, tiscores, wtype="TI", W[,,a], cohortsize,ncohort,ntrial,start.comb, ncores=6)
results