##### Load required R packages
library(nnet)
library(dfcrm)
library(R.utils)
library(foreach)
library(doSNOW)
library(doRNG)
library(parallel)

####### Function to obtain average toxicity scores at each dose level, 
####### which is the weighted sum of toxicity weights and the probability  
####### of toxicity at each dose level.    

get.thresh <- function(ntox, weights, TOX){
  tox.type.names <- paste("tox.type", seq(1:ntox), sep="")
  tox_type <- seq(from=0, to=4,by=1)
  tox.type <- matrix(rep(tox_type, each=ntox), ncol=ntox, byrow=T)
  tox.type <- as.data.frame(tox.type)
  possible_outcomes <- expand.grid(tox.type[,1:ntox]) #Permutes all possible toxicity profiles into a data set
  #map toxicity profiles to weights
  mapped.weight <- NA
  v <- NULL
  for(i in 1:ntox){
    for(k in 1:nrow(possible_outcomes)){
      mapped.weight[k] <- ifelse(possible_outcomes[k,i] > 0, weights[i,possible_outcomes[k,i]], 0)
    }
    v <- cbind(v, mapped.weight)
  }
  mapped_data <- matrix(v, ncol=3, byrow=F)
  scores <- apply(mapped_data, 1, sum) #Calulate toxicity scores for each profile in data set
  prob.data <- array(NA, c(nrow(possible_outcomes), ntox, d))
  for(j in 1:d){
    for(i in 1:nrow(possible_outcomes)){
      for(k in 1:ntox){
        tox.prob <- ifelse(possible_outcomes[i,k] > 0, TOX[j, possible_outcomes[i,k]+1, k], TOX[j,1,k])
        prob.data[i,k,j] <- tox.prob
      }
    }
  }
  prob.tox <- apply(prob.data[,,],c(1,3),prod)
  prob.scores <- scores*prob.tox
  thresh <- apply(prob.scores, 2, sum)
  return(thresh)
}

##### Function to compute Toxicity Index scores standardized to [0,1] given a toxicity profile
get_ti_score <- function(df, dlts){
  df <- as.numeric(df)
  
  for(k in 1:length(df)){
    df[k] <- ifelse(dlts[k]==4 & df[k]>1, df[k]-1, 
                    ifelse(dlts[k]==4 & df[k]==1, df[k]-.5, df[k]))
  }
  
  #order toxicities
  df <- df[order(-df)]
  
  #calculate weights
  ti_weights <- 1
  for(k in 2:(length(df))){
    ti_weights[k] <- (1+df[k-1])*prod(ti_weights[k-1])
  }
  
  #calculate ti_score
  ti_score <- sum(df/ti_weights)
  ti_score <- ti_score/5
  ti_score
}

##### Computes the mean Toxicity Index score at each dose level given a specified scenario and vector denoting 
#####  the AE grade that define a DLT for each AE type considered
mean.tiscore <- function(ntox, TOX, dlts){
  tox.type.names <- paste("tox.type", seq(1:ntox), sep="")
  tox_type <- seq(from=0, to=4,by=1)
  tox.type <- matrix(rep(tox_type, each=ntox), ncol=ntox, byrow=T)
  tox.type <- as.data.frame(tox.type)
  possible_outcomes <- expand.grid(tox.type[,1:ntox]) #Permutes all possible toxicity profiles into a data set
  dlts <- dlts
  #calculate TI Score
  scores <- apply(possible_outcomes, 1, get_ti_score, dlts=dlts) #Calulate toxicity scores for each profile in data set
  prob.data <- array(NA, c(nrow(possible_outcomes), ntox, d)) #Create an empty array indexed by (number of toxicity profiles, number of toxicities, number of dose levels)
  #fill in array with probability of no toxicity for each possible toxicity at each dose
  for(j in 1:d){
    for(i in 1:nrow(possible_outcomes)){
      for(k in 1:ntox){
        tox.prob <- ifelse(possible_outcomes[i,k] > 0, TOX[j, possible_outcomes[i,k]+1, k], TOX[j,1,k])
        prob.data[i,k,j] <- tox.prob
      }
    }
  }
  prob.tox <- apply(prob.data[,,],c(1,3),prod) #yields probability of each toxicity combination (rows) occurring at each dose level (columns)
  prob.scores <- scores*prob.tox #yields expected score of toxicity combination at each dose level 
  thresh <- apply(prob.scores, 2, sum) #Sum across columns to get expected toxicity score at each dose level 
  return(thresh)
}


#### Function to establish which grade of each AE type distinguishes a DLT based on a specified toxicity threshold
#### and specified toxicity weights. This ensures the DLT definition will be consistent with toxicity weights 
find.dlt <- function(W2, tul){
  dlts <- NA
  for(i in 1:ntox){
    dlts[i] <- ifelse(W2[i,1] >=tul, 1, 
                      ifelse(W2[i,2] >=tul, 2,
                             ifelse(W2[i,3] >=tul,3,4)))
  }
  return(dlts)
}

#### Function to calculate the probability of a DLT at each dose level using 
#### the individual toxicity probabilities across grades and doses
dlt.prob <- function(TOX, ntox, tox.dlt){
  dlt.matrix <- matrix(0, nrow=ntox, ncol=nrow(TOX))
  for(i in 1:ntox){
    xxx <- TOX[,-(1:tox.dlt[i]),i]
    xxx <- cbind(xxx, 0)
    dlt.matrix[i,] <- apply(xxx, 1, sum)
  }
  ptox <- NULL
  for(i in 1:d){
    ptox[i] <- 1-prod(1-dlt.matrix[,i]) #Ptox is probability of DLT at each dose level
  }
  return(ptox)
}

#### Function to obtain adjusted toxicity weights to range of 0-1
tox.weights <- function(W){
  thetamax <- sum(W[,4])
  W2 <- W/thetamax
  return(W2)
}


#### function to get skeletons for simulations, which uses the Cheung `dfcrm' getprior function
get.skeletons <- function(target1, hw, init.guess, d, s, orders){
  skeleton1 <- getprior(hw, target1, init.guess, d, model="empiric") #Get Skeleton for toxicity scores
  #For Tox Score Model
  p.skel<-matrix(0,nrow=s,ncol=d)
  for(j in 1:s){
    p.skel[j,]<-skeleton1[order(orders[j,])]
  }
  p.skel
}

#############################################################################################################################
#### qpocrm: Run a single trial that yields results for based on toxicity scores, TI scores, and binary DLT probability 
####         Note, this runs each trial design independently. This functions simply simulates all three trial designs at once
####         allowing for a direct comparison of results when averaged over multiple trial simulations. 
####         TOX: specified the model scenario being assessed. This should be an array with 3 dimensions. 
####         p.skel: model skeleton for weights-based toxicity scores
####         p.skel2: model skeleton for binary DLT model
####         p.skel3: model skeleton for TI model
####         tul: target toxicity score based on weights
####         tul2: target DLT probability for conventional model
####         tul3: target toxicity score baesd on Toxicity Index (TI)
####         tox.dlt: vector denoting the AE grades that distinguish a DLT for each AE type
####         ptox: mean toxicity probability at each dose-level (used for conventional model)
####         weights: user-defined weight matrix of dimension Lx5, where L represents the # of AEs considered
####         cohortsize: size of each enrolled cohort
####         ncohort: number of cohorts to run trial
####         start.comb: starting dose-level. 
#####################################################################################################################
qpocrm<-function(TOX, p.skel, p.skel2, p.skel3, tul, tul2, tul3, tox.dlt, ptox, weights, cohortsize,ncohort,start.comb){
  if(is.vector(p.skel)) p.skel=t(as.matrix(p.skel)); 
  nord.tox = nrow(p.skel); 
  mprior.tox = rep(1/nord.tox, nord.tox);  # prior probability for each ordering is 1/M (M is total number of orders)
  
  #Establish Normal Prior for scaling parameter based on 'Empiric' model dose response function 
  post.tox<-function(a,p,y,n){
    s2=1.34                                  # Prior variance is 1.34
    lik=1
    for(j in 1:length(p)){
      pj=p[j]**exp(a)                      # "Empiric" model dose reponse function for parameter 'a'          
      lik=lik*pj^y[j]*(1-pj)^(n[j]-y[j]);  # Bernoulli Likelihood (Same functional form for quasi-Bernoulli)
    }
    return(lik*exp(-0.5*a*a/s2));            # product of Bernoulli Likelihood and exponent from Normal prior 
    # non-exponent: (1/sqrt(2*pi*s2)) cancels out in posterior #
    # prior mean=0 and prior variance=1.34 (Oquigley 1996, Cheung 2011)
  }               
  
  #Function for calulation of the the posterior mean of parameter alpha which will be integrated over using Gauss-Hermite Quadrature
  posttoxf <- function(a, p, y, n, j) { p[j]^(exp(a))*post.tox(a, p, y, n); }
  
  ### Run Trial ###
  
  #Define initial parameters 
  ncomb = ncol(p.skel);   #number of drug combinations (i.e. dose levels)
  y=y2=y3= rep(0,ncomb);  #cumulative toxicity score and DLTs at each dose level for both models
  n=n2=n3= rep(0,ncomb);  #number of treated patients at each dose level for both models
  comb.curr = comb.curr2 = comb.curr3= start.comb;  # current dose level at each iteration of a trial	 
  ptox.hat = ptox.hat2 = ptox.hat3 = numeric(ncomb); # estimate of toxicity probability at each iteration of a trial 
  comb.select=comb.select2= comb.select3 = rep(0,ncomb); # a vector indicating selected dose for each cohort/patient in a trial
  Tox = Tox2 =  NA #Initializes matrix to hold toxicity values. Columns 1 = Tox A, column 2 = Tox B
  doses <- doses2 <- doses3 <- NULL
  i=1	
  while(i <= ncohort)
  {
    doses <- c(doses, comb.curr)
    doses2 <- c(doses2, comb.curr2)
    doses3 <- c(doses3, comb.curr3)
    
    for(z in 1:cohortsize){
      #Random occurrence of each toxicity and grade based on dose-level and corresponding tox probabilities for that dose level
      random.tox <- runif(ntox)
      for(k in 1:ntox){
        Tox[k] <- ifelse(random.tox[k] < TOX[comb.curr, 1, k] || TOX[comb.curr,1,k]==0, 0, 
                         ifelse(random.tox[k] < TOX[comb.curr, 2,k]+TOX[comb.curr,1,k], 1, 
                                ifelse(random.tox[k] < TOX[comb.curr,2,k]+TOX[comb.curr,1,k]+TOX[comb.curr, 3,k], 2,
                                       ifelse(random.tox[k] < TOX[comb.curr,2,k]+TOX[comb.curr,1,k]+TOX[comb.curr, 3,k]+TOX[comb.curr,4,k], 3, 4))))
      }
      
      #Calculate Toxicity Scores
      toxscores <- NA
      for(k in 1:ntox){
        toxscores[k] = ifelse(max(Tox[k]) == 0, 0, weights[k,max(Tox[k])]) 
      }
      toxscore = sum(toxscores) #The toxicity score for patient i
      
      y[comb.curr] = y[comb.curr] + toxscore; #gives cumulative toxicity score for each group
      n[comb.curr] = n[comb.curr] + 1;	#number of patients treated at each dose for toxicity score model
  
      ######## calculate toxicity index  #############
      random.tox <- runif(ntox)
      for(k in 1:ntox){
        Tox2[k] = ifelse(random.tox[k] < TOX[comb.curr3, 1, k] || TOX[comb.curr3,1,k]==0, 0, 
                         ifelse(random.tox[k] < TOX[comb.curr3, 2,k]+TOX[comb.curr3,1,k], 1, 
                                ifelse(random.tox[k] < TOX[comb.curr3,2,k]+TOX[comb.curr3,1,k]+TOX[comb.curr3, 3,k], 2,
                                       ifelse(random.tox[k] < TOX[comb.curr3,2,k]+TOX[comb.curr3,1,k]+TOX[comb.curr3, 3,k]+TOX[comb.curr3,4,k], 3, 4))))
      }
      ti_score = get_ti_score(Tox2, tox.dlt)
      y3[comb.curr3] = y3[comb.curr3] + ti_score; #gives cumulative toxicity score for each group
      n3[comb.curr3] = n3[comb.curr3] + 1;	#number of patients treated at each dose for toxicity score model
      
      #Calculate DLTs
      y2[comb.curr2] = y2[comb.curr2] + rbinom(1,1,ptox[comb.curr2]); #gives cumulative number of DLTs for each group
      n2[comb.curr2] = n2[comb.curr2] + 1;	#number of patients treated at each dose for binary model
    }
    
    ########################################################
    # Posterior probability estimates for each order using #
    # 'integrate' function to integrate over the parameter #
    #  space for 'a' [-inf, inf]                           # 
    ########################################################
    
    #For Toxicity Score Model
    marginal.tox = rep(0, nord.tox);
    for(k in 1:nord.tox){
      marginal.tox[k] = integrate(post.tox,lower=-Inf,upper=Inf, p=p.skel[k,], y=y,n=n)$value;
    }
    postprob.tox = (marginal.tox*mprior.tox)/sum(marginal.tox*mprior.tox); #probability of each order for tox score model
    
    #For Binary-DLT model
    marginal.tox2 = rep(0, nord.tox);
    for(k in 1:nord.tox){
      marginal.tox2[k] = integrate(post.tox,lower=-Inf,upper=Inf, p=p.skel2[k,], y=y2,n=n2)$value;
    }
    
    postprob.tox2 = (marginal.tox2*mprior.tox)/sum(marginal.tox2*mprior.tox); #probability of each order for binary DLT model
    
    #For TI Score model
    marginal.tox3 = rep(0, nord.tox);
    for(k in 1:nord.tox){
      marginal.tox3[k] = integrate(post.tox,lower=-Inf,upper=Inf, p=p.skel3[k,], y=y3,n=n3)$value;
    }
    
    postprob.tox3 = (marginal.tox3*mprior.tox)/sum(marginal.tox3*mprior.tox); #probability of each order for binary DLT model
    
    ########################################################
    # Using poseterior probabilities for each model order, #
    # select the assumed order. For the first 1/3 patients #
    # order is selected based off the posterior probability#
    # of each order. After the first 1/3 patients, the     #
    # order with the highest posterior probability is      #
    # selected for the next patient.                       #
    ########################################################
    
    #Select order for toxicity score model
    if(nord.tox>1 && i <= ncohort/3){ 
      mtox.sel = sample(1:nord.tox, size=1, p=postprob.tox)
    } else if(nord.tox>1 && i > ncohort/3){
      mtox.sel = which.is.max(postprob.tox); 
    } else{
      mtox.sel = 1;
    }
    
    #For order for bianry DLT model 
    if(nord.tox>1 && i <= ncohort/3){ 
      mtox.sel2 =sample(1:nord.tox, size=1, p=postprob.tox2)
    } else if(nord.tox>1 && i > ncohort/3){
      mtox.sel2 = which.is.max(postprob.tox2); 
    } else{
      mtox.sel2 = 1;
    }
    
    #For order for TI score model 
    if(nord.tox>1 && i <= ncohort/3){ 
      mtox.sel3 =sample(1:nord.tox, size=1, p=postprob.tox3)
    } else if(nord.tox>1 && i > ncohort/3){
      mtox.sel3 = which.is.max(postprob.tox3); 
    } else{
      mtox.sel3 = 1;
    }
    
    ##################################################################
    # Calculate posterior mean of toxicity score and toxicity        #
    # probability for each model using via integration for each      #
    # drug combination assuming the dose order selected in the step  #
    # above.                                                         #
    ##################################################################
    
    # For Toxicity Score Model
    for(j in 1:ncomb){
      ptox.hat[j] = integrate(posttoxf,lower=-Inf,upper=Inf, p.skel[mtox.sel,], y, n,j)$value/marginal.tox[mtox.sel]; 
    }
    
    # For probabiltiy of DLT Model
    for(j in 1:ncomb){
      ptox.hat2[j] = integrate(posttoxf,lower=-Inf,upper=Inf, p.skel2[mtox.sel2,], y2, n2,j)$value/marginal.tox2[mtox.sel2]; 
    }
    
    # For TI score model
    for(j in 1:ncomb){
      ptox.hat3[j] = integrate(posttoxf,lower=-Inf,upper=Inf, p.skel3[mtox.sel3,], y3, n3,j)$value/marginal.tox3[mtox.sel3]; 
    }
    
    ##################################################################
    ############ Select the next dose for each model #################
    ##################################################################
    distance=abs(ptox.hat-tul) 
    distance2 = abs(ptox.hat2-tul2)
    distance3 = abs(ptox.hat3 - tul3)
    comb.best= which.is.max(-distance) #Selects dose closest to target score for toxicity score model
    comb.best2 = which.is.max(-distance2) #Selects dose closest to target DLT rate for binary CRM model
    comb.best3 = which.is.max(-distance3)

    comb.curr=min(comb.best,comb.curr+1) #Selects next dose, but doesn't allow dose-skipping of untried dose-levels
    comb.curr2=min(comb.best2, comb.curr2+1) #selects next dose, but doesn't allow dose skipping of untried dose-levels
    comb.curr3=min(comb.best3, comb.curr3+1)
    i=i+1
  }
  comb.select[comb.curr]=comb.select[comb.curr]+1;
  comb.select2[comb.curr2]=comb.select2[comb.curr2]+1;
  comb.select3[comb.curr3]=comb.select3[comb.curr3]+1;
  return(list(comb.select=comb.select, comb.select2=comb.select2, comb.select3=comb.select3, 
              tox.data=y, tox.data2 = y2, tox.data3 = y3, 
              pt.allocation=n, pt.allocation2 = n2, pt.allocation3 = n3, msel=mtox.sel, msel2= mtox.sel2, msel3= mtox.sel3,
              doses = doses, doses2=doses2, doses3=doses3))
}

#######################################################################################################################################
####### qpocrm.oc: Run a single trial using the Overdose constraint, specifying wtype='TI" for toxicity index or wtype="weights" for 
#######            scores based on weights. If wtype="TI", the 'weights' argument is ignored. 
#######            p.skel: skeleton for toxicity score model being used
#######            p.skel2: skeleton for binary DLT model
#######            tul: target toxicity score
#######            tul2: target DLT probability
#######            *other parameters are same as in function 'qpocrm' defined above
########################################################################################################################################
qpocrm.oc<-function(TOX, p.skel, p.skel2, tul, tul2, tox.dlt, wtype="TI", weights=NULL, cohortsize,ncohort,start.comb){
  if(is.vector(p.skel)) p.skel=t(as.matrix(p.skel)); 
  
  nord.tox = nrow(p.skel); 
  mprior.tox = rep(1/nord.tox, nord.tox);  # prior for each toxicity ordering 
  
  post.tox<-function(a,p,y,n){
    s2=1.34
    lik=1
    for(j in 1:length(p)){
      pj=p[j]**exp(a)
      lik=lik*pj^y[j]*(1-pj)^(n[j]-y[j]);
    }
    return(lik*exp(-0.5*a*a/s2));
  }
  
    #the posterior mean of ptox
  posttoxf <- function(a, p, y, n, j) { p[j]^(exp(a))*post.tox(a, p, y, n); }
  ################## POWER MODEL ###########################
  
  ###run a trial 	
  ncomb = ncol(p.skel);   #number of combos
  y=y2 = rep(0,ncomb);  #cumulative toxicity score at each dose level
  n=n2 = rep(0,ncomb);  #number of treated patients at each dose level
  comb.curr = comb.curr2 = start.comb;  # current dose level	 
  ptox.hat = ptox.hat2  = numeric(ncomb); # estimate of toxicity prob
  comb.select=comb.select2= rep(0,ncomb); # a vector of indicators for dose selection
  Tox = NA #Defines matrix to hold toxicity values. Columns 1 = Tox A, column 2 = Tox B
  doses <- NULL
  i=1	
  while(i <= ncohort)
  { 
    doses <- c(doses, comb.curr)
    #generate data for a new cohort of patients
    #Toxicities given dose under tox scores
    #Model 1 with tox scores and no safety
    for(z in 1:cohortsize){
      random.tox <- runif(ntox)
      for(k in 1:ntox){
        Tox[k] = ifelse(random.tox[k] < TOX[comb.curr, 1, k] || TOX[comb.curr,1,k]==0, 0, 
                         ifelse(random.tox[k] < TOX[comb.curr, 2,k]+TOX[comb.curr,1,k], 1, 
                                ifelse(random.tox[k] < TOX[comb.curr,2,k]+TOX[comb.curr,1,k]+TOX[comb.curr, 3,k], 2,
                                       ifelse(random.tox[k] < TOX[comb.curr,2,k]+TOX[comb.curr,1,k]+TOX[comb.curr, 3,k]+TOX[comb.curr,4,k], 3, 4))))
      }
      
      
      #Get DLTs for binary model
      DLTs = NULL
      for(k in 1:ntox){
        DLTs[k] = ifelse(Tox[k]>=tox.dlt[k], 1,0)
      }
      DLT = ifelse(sum(DLTs)>0, 1, 0)
      
      y2[comb.curr] = y2[comb.curr] + DLT;
      
      #Get toicity scores based on TI or weights
      toxscores = NA
      if(wtype=="weights"){
        for(k in 1:ntox){
          toxscores[k] <- ifelse(max(Tox[k]) == 0, 0, weights[k,max(Tox[k])])
        }
        toxscore <- sum(toxscores)
  
        y[comb.curr] = y[comb.curr] + toxscore; #gives the total toxicity score for each group
        n[comb.curr] = n[comb.curr] + 1;
      } else if(wtype=="TI"){
      
        ######## calculate toxicity index  #############
        ti_score = get_ti_score(Tox, tox.dlt)
  
        y[comb.curr] = y[comb.curr] + ti_score; #gives the total toxicity score for each group
        n[comb.curr] = n[comb.curr] + 1;
      } 
    }
    
    #For tox scores
    marginal.tox = NA;
    for(k in 1:nord.tox){
      marginal.tox[k] = integrate(post.tox,lower=-Inf,upper=Inf, p=p.skel[k,], y=y,n=n)$value;
    }
    postprob.tox = (marginal.tox*mprior.tox)/sum(marginal.tox*mprior.tox); #probability of each order
    
    marginal.tox2 = NA;
    for(k in 1:nord.tox)
    {
      marginal.tox2[k] = integrate(post.tox,lower=-Inf,upper=Inf, p=p.skel2[k,], y=y2,n=n)$value;
    }
    postprob.tox2 = (marginal.tox2*mprior.tox)/sum(marginal.tox2*mprior.tox); #probability of each order
    
    ##Bayesian model selection, identify the model with the highest posterior prob:
    #For Tox scores
    if(nord.tox>1 && i <= ncohort/3){ 
      mtox.sel = sample(1:nord.tox, size=1, p=postprob.tox)
    } else if(nord.tox>1 && i > ncohort/3){
      mtox.sel = which.is.max(postprob.tox); 
    } else{
      mtox.sel = 1;
    }
    
    ##calculate posterior mean of toxicity probability at each combo
    for(j in 1:ncomb){
      ptox.hat[j] = integrate(posttoxf,lower=-Inf,upper=Inf, p.skel[mtox.sel,], y, n,j)$value/marginal.tox[mtox.sel]; 
    }
    
    for(j in 1:ncomb){
      ptox.hat2[j] = integrate(posttoxf,lower=-Inf,upper=Inf, p.skel2[mtox.sel,], y2, n,j)$value/marginal.tox2[mtox.sel]; 
    }
    
    distance=abs(ptox.hat-tul)
    distance2 = abs(ptox.hat2-tul2)
    comb.best= which.is.max(-distance)
    comb.best2 = which.is.max(-distance2)
    comb.best = ifelse(ptox.hat2[comb.best2] < ptox.hat2[comb.best], comb.best2, comb.best)
    #comb.best = ifelse(ptox.hat2[comb.best] < tul2+.05, comb.best, comb.best2)
    #nexthighest=ifelse(ptox.hat[comb.curr]==max(ptox.hat), which.max(ptox.hat),which(ptox.hat==min(ptox.hat[(ptox.hat>ptox.hat[comb.curr])])))
    comb.curr=min(comb.best,comb.curr+1)
    i=i+1
  }
  comb.select[comb.curr]=comb.select[comb.curr]+1;
  return(list(comb.select=comb.select, tox.data=y, tox.data2= y2, pt.allocation=n, msel=mtox.sel))
}

###########################################################################################################
#### qpocrm.oc2: This is nearly the same as 'qpocrm.oc' except that it flips the overdose constraint 
####             such that the comparison of selected doses with respect to DLT probability are done 
####             so under the order established by the DLT model (as discussed in Appendix F)
#########################################################################################################
qpocrm.oc2 <-function(TOX, p.skel, p.skel2, tul, tul2, tox.dlt, wtype="TI", weights=NULL, cohortsize,ncohort,start.comb){
  if(is.vector(p.skel)) p.skel=t(as.matrix(p.skel)); 
  
  nord.tox = nrow(p.skel); 
  mprior.tox = rep(1/nord.tox, nord.tox);  # prior for each toxicity ordering 
  
  post.tox<-function(a,p,y,n){
    s2=1.34
    lik=1
    for(j in 1:length(p)){
      pj=p[j]**exp(a)
      lik=lik*pj^y[j]*(1-pj)^(n[j]-y[j]);
    }
    return(lik*exp(-0.5*a*a/s2));
  }
  
  #the posterior mean of ptox
  posttoxf <- function(a, p, y, n, j) { p[j]^(exp(a))*post.tox(a, p, y, n); }
  ################## POWER MODEL ###########################
  
  ###run a trial 	
  ncomb = ncol(p.skel);   #number of combos
  y=y2 = rep(0,ncomb);  #cumulative toxicity score at each dose level
  n=n2 = rep(0,ncomb);  #number of treated patients at each dose level
  comb.curr = comb.curr2 = start.comb;  # current dose level	 
  ptox.hat = ptox.hat2  = ptox.hat3 = numeric(ncomb); # estimate of toxicity prob
  comb.select=comb.select2= rep(0,ncomb); # a vector of indicators for dose selection
  Tox = NA #Defines matrix to hold toxicity values. Columns 1 = Tox A, column 2 = Tox B
  i=1	
  while(i <= ncohort)
  { 
    #generate data for a new cohort of patients
    #Toxicities given dose under tox scores
    #Model 1 with tox scores and no safety
    for(z in 1:cohortsize){
      random.tox <- runif(ntox)
      for(k in 1:ntox){
        Tox[k] = ifelse(random.tox[k] < TOX[comb.curr, 1, k] || TOX[comb.curr,1,k]==0, 0, 
                        ifelse(random.tox[k] < TOX[comb.curr, 2,k]+TOX[comb.curr,1,k], 1, 
                               ifelse(random.tox[k] < TOX[comb.curr,2,k]+TOX[comb.curr,1,k]+TOX[comb.curr, 3,k], 2,
                                      ifelse(random.tox[k] < TOX[comb.curr,2,k]+TOX[comb.curr,1,k]+TOX[comb.curr, 3,k]+TOX[comb.curr,4,k], 3, 4))))
      }
      
      
      #Get DLTs for binary model
      DLTs = NULL
      for(k in 1:ntox){
        DLTs[k] = ifelse(Tox[k]>=tox.dlt[k], 1,0)
      }
      DLT = ifelse(sum(DLTs)>0, 1, 0)
      
      y2[comb.curr] = y2[comb.curr] + DLT;
      
      #Get toicity scores based on TI or weights
      toxscores = NA
      if(wtype=="weights"){
        for(k in 1:ntox){
          toxscores[k] <- ifelse(max(Tox[k]) == 0, 0, weights[k,max(Tox[k])])
        }
        toxscore <- sum(toxscores)
        
        y[comb.curr] = y[comb.curr] + toxscore; #gives the total toxicity score for each group
        n[comb.curr] = n[comb.curr] + 1;
      } else if(wtype=="TI"){
        
        ######## calculate toxicity index  #############
        ti_score = get_ti_score(Tox, tox.dlt)
        
        y[comb.curr] = y[comb.curr] + ti_score; #gives the total toxicity score for each group
        n[comb.curr] = n[comb.curr] + 1;
      } 
    }
    
    #For tox scores
    marginal.tox = NA;
    for(k in 1:nord.tox){
      marginal.tox[k] = integrate(post.tox,lower=-Inf,upper=Inf, p=p.skel[k,], y=y,n=n)$value;
    }
    postprob.tox = (marginal.tox*mprior.tox)/sum(marginal.tox*mprior.tox); #probability of each order
    
    marginal.tox2 = NA;
    for(k in 1:nord.tox)
    {
      marginal.tox2[k] = integrate(post.tox,lower=-Inf,upper=Inf, p=p.skel2[k,], y=y2,n=n)$value;
    }
    postprob.tox2 = (marginal.tox2*mprior.tox)/sum(marginal.tox2*mprior.tox); #probability of each order
    
    ##Bayesian model selection, identify the model with the highest posterior prob:
    #For Tox scores
    if(nord.tox>1 && i <= ncohort/3){ 
      mtox.sel = sample(1:nord.tox, size=1, p=postprob.tox)
    } else if(nord.tox>1 && i > ncohort/3){
      mtox.sel = which.is.max(postprob.tox); 
    } else{
      mtox.sel = 1;
    }
    
    #for DLT prob
    if(nord.tox>1 && i <= ncohort/3){ 
      mtox.sel2 = sample(1:nord.tox, size=1, p=postprob.tox2)
    } else if(nord.tox>1 && i > ncohort/3){
      mtox.sel2 = which.is.max(postprob.tox2); 
    } else{
      mtox.sel2 = 1;
    }
    
    ##calculate posterior mean of toxicity probability at each combo
    for(j in 1:ncomb){
      ptox.hat[j] = integrate(posttoxf,lower=-Inf,upper=Inf, p.skel[mtox.sel,], y, n,j)$value/marginal.tox[mtox.sel]; 
    }
    
    for(j in 1:ncomb){
      ptox.hat2[j] = integrate(posttoxf,lower=-Inf,upper=Inf, p.skel2[mtox.sel,], y2, n,j)$value/marginal.tox2[mtox.sel]; 
    }
    
    for(j in 1:ncomb){
      ptox.hat3[j] = integrate(posttoxf,lower=-Inf,upper=Inf, p.skel2[mtox.sel2,], y2, n,j)$value/marginal.tox2[mtox.sel2]; 
    }
    
    distance=abs(ptox.hat-tul)
    distance2 = abs(ptox.hat2-tul2)
    comb.best= which.is.max(-distance)
    comb.best2 = which.is.max(-distance2)
    # if(ptox.hat2[comb.best2] < ptox.hat2[comb.best]){
    #   distance2 = abs(ptox.hat3-tul2)
    #   comb.best <- which.is.max(-distance2)
    # } else{
    #   comb.best <- comb.best
    # }
    
    comb.best = ifelse(ptox.hat2[comb.best2] < ptox.hat2[comb.best], comb.best2, comb.best)
    comb.curr=min(comb.best,comb.curr+1)
    i=i+1
  }
  comb.select[comb.curr]=comb.select[comb.curr]+1;
  return(list(comb.select=comb.select, tox.data=y, tox.data2= y2, pt.allocation=n, msel=mtox.sel))
}

# THe following set up functions 'qp', 'qp.oc', and 'qp.oc2" call the model functions above
# and relabel the results output so that they can be more easily summarized and ran in parallel
qp <- function(TOX, p.skel, p.skel2, p.skel3, tul, tul2, tul3, tox.dlt, p0, ptox,tiscores, weights, cohortsize, ncohort, start.comb){    
  result<-qpocrm(TOX, p.skel, p.skel2, p.skel3, tul, tul2, tul3, tox.dlt, ptox, weights, cohortsize, ncohort, start.comb)
  comb.select=result$comb.select
  comb.select2=result$comb.select2
  comb.select3 = result$comb.select3
  y=result$tox.data
  y2=result$tox.data2
  y3 = result$tox.data3
  n=result$pt.allocation
  n2=result$pt.allocation2
  n3 = result$pt.allocation3
  dose = c(1:d)
  msel=result$msel
  msel2 = result$msel2
  msel3 = result$msel3
  doses = result$doses
  doses2 = result$doses2
  doses3 = result$doses3
  dat <- cbind(dose, p0, ptox, tiscores, comb.select, comb.select2, comb.select3, y, y2, y3, n, n2, n3, msel, msel2, msel3)
} 

qp.oc <- function(TOX, p.skel, p.skel2, tul, tul2, tox.dlt, p0, ptox, tiscores, wtype, weights, cohortsize, ncohort, start.comb){    
  result<-qpocrm.oc(TOX, p.skel, p.skel2, tul, tul2, tox.dlt, wtype, weights, cohortsize, ncohort, start.comb)
  comb.select=result$comb.select
  y=result$tox.data
  y2=result$tox.data2
  n=result$pt.allocation
  dose = c(1:d)
  msel=result$msel
  dat <- cbind(dose, p0, ptox, tiscores, comb.select, y, y2, n, msel)
} 

qp.oc2 <- function(TOX, p.skel, p.skel2, tul, tul2, tox.dlt, p0, ptox, tiscores, wtype, weights, cohortsize, ncohort, start.comb){    
  result<-qpocrm.oc2(TOX, p.skel, p.skel2, tul, tul2, tox.dlt, wtype, weights, cohortsize, ncohort, start.comb)
  comb.select=result$comb.select
  y=result$tox.data
  y2=result$tox.data2
  n=result$pt.allocation
  msel=result$msel
  dose = c(1:d)
  dat <- cbind(dose, p0, ptox, tiscores, comb.select, y, y2, n, msel)
} 


######################################################################################################
####  The following functions run the trials in parallel tos summarize across multiple trials. The additional 
####  parameters are added beyond those defined for the 'qpocrm' function above. 
####  p0: the mean toxicity scores at each dose level (these are used for presentation of results)
####  tiscores: the mean TI scores at each dose level (these are used for presentation of results)
####  ntrial: number of trials to run and average over
####  ncores: the number of CPU cores to run trials over in parallel
####  seed: the random seed that can be specified for reproducible results
####################################################################################################
qpocrm.sim<-function(TOX, p.skel, p.skel2, p.skel3, tul, tul2, tul3, tox.dlt, p0, ptox, tiscores, weights, cohortsize,ncohort, ntrial,start.comb, ncores=1, seed=NA){
  ncomb=nrow(TOX)
  c=ncores
  cl = makeSOCKcluster(c) #make Cluster
  registerDoSNOW(cl) #register cluseter
  pb2 <- txtProgressBar(min=1, max=ntrial, style=3)
  progress = function(n) setTxtProgressBar(pb2, n)
  opts =list(progress=progress)
  if(is.na(seed)==T){
    res <-  foreach(i = 1:ntrial, .combine='rbind', .packages="nnet", .options.snow=opts, 
                    .export=c("qpocrm", "qp", "get_ti_score")) %dopar% {
                      ntox <- length(tox.dlt)
                      d <- ncomb
                      qp(TOX, p.skel, p.skel2, p.skel3, tul, tul2, tul3, tox.dlt, p0, ptox, tiscores, weights, cohortsize,ncohort,start.comb)
                    }
  } else {
    set.seed(seed)
    res <-  foreach(i = 1:ntrial, .combine='rbind', .packages="nnet", .options.snow=opts, 
                    .export=c("qpocrm", "qp", "get_ti_score")) %dorng% {
                      ntox <- length(tox.dlt)
                      d <- ncomb
                      qp(TOX, p.skel, p.skel2, p.skel3, tul, tul2, tul3, tox.dlt, p0, ptox, tiscores, weights, cohortsize,ncohort,start.comb)
                    }
  }
  close(pb2)
  stopCluster(cl)
  res <- data.frame(res)
  comb.select2 = comb.select = comb.select3 = y3 =  y2 = y = n3 = n2 = n = msel=msel2=msel3= NULL
  for(k in 1:ncomb){
    comb.select=cbind(comb.select, res[which(res$dose==k),]$comb.select)
    comb.select2=cbind(comb.select2, res[which(res$dose==k),]$comb.select2)
    comb.select3=cbind(comb.select3, res[which(res$dose==k),]$comb.select3)
    y=cbind(y, res[which(res$dose==k),]$y)
    y2=cbind(y2, res[which(res$dose==k),]$y2)
    y3=cbind(y3, res[which(res$dose==k),]$y3)
    n=cbind(n, res[which(res$dose==k),]$n)
    n2=cbind(n2, res[which(res$dose==k),]$n2)
    n3=cbind(n3, res[which(res$dose==k),]$n3)
  }
  msel <- table(res$msel)/(ncomb*ntrial)
  msel2 <- table(res$msel2)/(ncomb*ntrial)
  msel3 <- table(res$msel3)/(ncomb*ntrial)
  return(list(`True Mean Scores (weights)`=round(p0,2), `True DLT Probabilities`=ptox, `True Mean TI scores` = round(tiscores,2), 
              `MTD Selection % based on toxscores (weights)`=colMeans(comb.select), `MTD Selection % based on DLT Prob`=colMeans(comb.select2), `MTD Selection % based on TI`=colMeans(comb.select3),
              #tox.data=colMeans(y), tox.data2 = colMeans(y2), tox.data3 =colMeans(y3), 
              `PT Allocation by dose base on scores (weights)`=round(colMeans(n),1), `PT Allocation by dose base on DLT prob` = round(colMeans(n2),1), 
              `PT Allocation by dose base on TI` = round(colMeans(n3),1),
              `Order selection % (scores)`=msel, `Order selection % (DLT)` =msel2, `Order selection % (TI)`=msel3))
}



qpocrm.sim.oc<-function(TOX, p.skel, p.skel2, tul, tul2, tox.dlt, p0, ptox, tiscores, wtype, weights, cohortsize,ncohort,ntrial,start.comb, ncores=1,  seed=NA){
  ncomb=d
  c=ncores
  cl <- makeSOCKcluster(c) #make Cluster
  registerDoSNOW(cl) #register cluseter
  pb2 <- txtProgressBar(min=1, max=ntrial, style=3)
  progress <- function(n) setTxtProgressBar(pb2, n)
  opts <- list(progress=progress)
  if(is.na(seed)==T){
    res <-  foreach(i = 1:ntrial, .combine='rbind', .packages="nnet", .options.snow=opts, 
                    .export=c("qpocrm.oc", "qp.oc", "get_ti_score")) %dopar% {
                      ntox <- length(tox.dlt)
                      d <- ncomb
                      qp.oc(TOX, p.skel, p.skel2, tul, tul2, tox.dlt, p0, ptox, tiscores, wtype, weights, cohortsize,ncohort,start.comb)
                    }
  } else {
    set.seed(seed)
    res <-  foreach(i = 1:ntrial, .combine='rbind', .packages="nnet", .options.snow=opts, 
                    .export=c("qpocrm.oc", "qp.oc", "get_ti_score")) %dorng% {
                      ntox <- length(tox.dlt)
                      d <- ncomb
                      qp.oc(TOX, p.skel, p.skel2, tul, tul2, tox.dlt, p0, ptox, tiscores, wtype, weights, cohortsize,ncohort,start.comb)
                    }
  }
  close(pb2)
  stopCluster(cl)
  res <- data.frame(res)
  comb.select = y2 = y = n = NULL
  for(k in 1:ncomb){
    comb.select=cbind(comb.select, res[which(res$dose==k),]$comb.select)
    y=cbind(y, res[which(res$dose==k),]$y)
    y2=cbind(y2, res[which(res$dose==k),]$y2)
    n=cbind(n, res[which(res$dose==k),]$n)
    
  }
  msel <- table(res$msel)/(ncomb*ntrial)
  return(list(`True Mean Scores (weights)`=round(p0,2), `True DLT Probabilities`=ptox, `True Mean TI scores` = round(tiscores,2),
              `MTD Selection %`=colMeans(comb.select), 
              #tox.data=colMeans(y), tox.data2 = colMeans(y2), 
              `PT Allocation by dose`=round(colMeans(n),1), 
              `Order Selection %`=msel))
}


qpocrm.sim.oc2<-function(TOX, p.skel, p.skel2, tul, tul2, tox.dlt, p0, ptox, tiscores, wtype, weights, cohortsize,ncohort,ntrial,start.comb, ncores=1,  seed=NA){
  ncomb=d
  c=ncores
  cl <- makeSOCKcluster(c) #make Cluster
  registerDoSNOW(cl) #register cluseter
  pb2 <- txtProgressBar(min=1, max=ntrial, style=3)
  progress <- function(n) setTxtProgressBar(pb2, n)
  opts <- list(progress=progress)
  if(is.na(seed)==T){
    res <-  foreach(i = 1:ntrial, .combine='rbind', .packages="nnet", .options.snow=opts, 
                    .export=c("qpocrm.oc2", "qp.oc2", "get_ti_score")) %dopar% {
                      ntox <- length(tox.dlt)
                      d <- ncomb
                      qp.oc2(TOX, p.skel, p.skel2, tul, tul2, tox.dlt, p0, ptox, tiscores, wtype, weights, cohortsize,ncohort,start.comb)
                    }
  } else {
    set.seed(seed)
    res <-  foreach(i = 1:ntrial, .combine='rbind', .packages="nnet", .options.snow=opts, 
                    .export=c("qpocrm.oc2", "qp.oc2", "get_ti_score")) %dorng% {
                      ntox <- length(tox.dlt)
                      d <- ncomb
                      qp.oc2(TOX, p.skel, p.skel2, tul, tul2, tox.dlt, p0, ptox, tiscores, wtype, weights, cohortsize,ncohort,start.comb)
                    }
  }
  close(pb2)
  stopCluster(cl)
  res <- data.frame(res)
  comb.select = y2 = y = n = NULL
  for(k in 1:ncomb){
    comb.select=cbind(comb.select, res[which(res$dose==k),]$comb.select)
    y=cbind(y, res[which(res$dose==k),]$y)
    y2=cbind(y2, res[which(res$dose==k),]$y2)
    n=cbind(n, res[which(res$dose==k),]$n)
  }
  msel <- table(res$msel)/(ncomb*ntrial)
  return(list(`True Mean Scores (weights)`=round(p0,2), `True DLT Probabilities`=ptox, `True Mean TI scores` = round(tiscores,2),
              `MTD Selection %`=colMeans(comb.select), 
              #tox.data=colMeans(y), tox.data2 = colMeans(y2), 
              `PT Allocation by dose`=round(colMeans(n),1),
              `Order Selection %`=msel))
}

