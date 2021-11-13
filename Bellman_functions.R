Bellman <- function(coeff, params,V0){
  # Calls BellmanParams to create the required matrices for fast computation
  # of the Bellman equation. Then calls BellmanLoop for the actual fixed point
  # solution
  
  # Construct matrices to compute the Bellman equation
  if (!exists("BellmanParams",where = params)){
    params$BellmanParams <- PackBellmanParams(params)
  }

  # Call BellmanLoop
  if(!missing(V0)) {
    return(BellmanLoop(coeff, params,V0))
  }
  else {
    return(BellmanLoop(coeff, params))
  }
  
}

PackBellmanParams <- function(params){
  # Create matrices that allow fast computation of the Bellman equation
  
  # Other parameters
  n_omega1 <- params$n_omega1          # number of omega1 states
  N1 <- params$N1                      # number of states for each omega1
  nstates <- n_omega1*N1               # Number of states
  n_regoutcomes <- params$nregoutcomes # Number of regulatory outcomes (80)
  n_transitions <- params$ntransitions # violator transtions (3)
  dep <- params$DAV_deprate            # depreciation rate
  beta <- params$beta                  # discount factor
  gamma <- params$gamma                # euler's constant
  
  # Initialize matrices 
  Vtilde <- matrix(0,nstates,1)
  NewV   <- matrix(0,nstates,1)
  OldV   <- matrix(0,nstates,1)
  Investprob <- matrix(0,nstates,1)
  
  # Get ordered values if the state for each index 
  # states are (omega1,lagInv1,lagInv2,orderViol,DAV) for NewV or OldV
  # or         (omega1,lagInv1,currentV,orderViol,DAV) for vtilde
  index <- c(1:nstates)
  states <- t(sapply(index,IndexToValues)) 
  # get states when the plant invests (from omegaTilde to omega)
  statesI <- states
  statesI[,2] <- 1 # lag1 of investment
  statesI[,3] <- states[,2] # lag2 of investment
  # get states when the plant does not invest
  statesNI <- states
  statesNI[,2] <- 0 # lag1 of investment
  statesNI[,3] <- states[,2] # lag2 of investment
  # get index when the plant is in compliance
  statesComp <- states
  statesComp[,2:4]<-0
  statesComp[,5]<-1
  indexComp <- apply(statesComp,1,ValuesToIndex)
  
  # Get DAV values
  DAVgrid <- params$DAVgrid
  DAV <- DAVgrid[states[,5]]
  npointsDAV  <- length(DAVgrid)
  maxpointDAV <- DAVgrid[npointsDAV]
  widthDAV <- DAVgrid[2]-DAVgrid[1]
  
  # compute updated DAV
  newDAV <- (1-dep)*DAV + states[,3]
  
  # Interpolate newDAV using the grid
  #indexes
  i_below <- pmin(1 + floor(newDAV*(npointsDAV-1)/maxpointDAV),npointsDAV)
  i_above <- pmin(2 + floor(newDAV*(npointsDAV-1)/maxpointDAV),npointsDAV)
  # weights
  w_below <- pmax((DAVgrid[i_above]-newDAV)/widthDAV,0)
  w_above <- 1-w_below
  
  # Solve OldV if investment
  statesI_below <- statesI
  statesI_below[,5] <- i_below
  indexI_below <- apply(statesI_below,1,ValuesToIndex)
  statesI_above <- statesI
  statesI_above[,5] <- i_above
  indexI_above <- apply(statesI_above,1,ValuesToIndex)
  
  # Solve OldV if no investment
  statesNI_below <- statesNI
  statesNI_below[,5] <- i_below
  indexNI_below <- apply(statesNI_below,1,ValuesToIndex)
  statesNI_above <- statesNI
  statesNI_above[,5] <- i_above
  indexNI_above <- apply(statesNI_above,1,ValuesToIndex)
  
  
  # Compute U, CaseProb and indexes for each case
  # create lists of names
  ins_list <- c()
  vio_list <- c()
  fine_list <- c()
  prob_list <- c()
  for (j_ro in 1:n_regoutcomes){
    
    ins_name <- paste0("inspection",j_ro)
    ins_list <- c(ins_list,ins_name)
    
    vio_name <- paste0("violation",j_ro)
    vio_list <- c(vio_list,vio_name)
    
    fine_name <- paste0("fine",j_ro)
    fine_list <- c(fine_list,fine_name)
    
    prob_name <- paste0("prob",j_ro)
    prob_list <- c(prob_list,prob_name)
  }
  
  # create j_ro and j_tr indexes
  jRo <- rep(NA, n_regoutcomes*n_transitions)
  jTr <- rep(NA, n_regoutcomes*n_transitions)
  for (j_ro in 1:n_regoutcomes){
    for (j_tr in 1:n_transitions){
      jRo[(j_tr-1)*n_regoutcomes+j_ro] <- j_ro
      jTr[(j_tr-1)*n_regoutcomes+j_ro] <- j_tr
    }
  }
  
  # Compute U and CaseProb (probability of each case)
  
  insCase  <- params$data[index,ins_list[jRo]]
  vioCase  <- params$data[index,vio_list[jRo]]
  fineCase <- params$data[index,fine_list[jRo]]
  HPVCase  <- matrix(rep(jTr,nstates),nstates,n_regoutcomes*n_transitions,byrow = T)
  HPVCase  <- 1*(HPVCase == 3) 
  
  # Compute probability of each case
  prob_ro <- params$data[index,prob_list[jRo]]
  transition_name <- paste0("transition",jTr-1,"_",jRo)
  prob_tr  <- params$data[index,transition_name]
  CaseProb <- prob_ro*prob_tr
  
  # Reshape as vectors
  CaseProb <- c(as.matrix(CaseProb))
  
  # compute index for violation cases
  indexVioCase <- matrix(rep(index,n_regoutcomes*n_transitions),nstates,n_regoutcomes*n_transitions)
  # loop through columns (cases) to compute the transition
  for (i in 1:(n_regoutcomes*n_transitions)){
    # get state values
    stateCase <- t(sapply(indexVioCase[,i],IndexToValues))
    
    # change current violation
    stateCase[,3] <- vioCase[,i]
    
    # change HPV status
    stateCase[,4] <- jTr[i]-1
    
    # recompute the index
    indexVioCase[,i] <- apply(stateCase,1,ValuesToIndex)
  }
  
  # Pack Bellman parameters
  paramsBellman <- list()
  paramsBellman$NewV              <- NewV          
  paramsBellman$OldV              <- OldV          
  paramsBellman$w_below           <- w_below       
  paramsBellman$w_above           <- w_above       
  paramsBellman$indexI_below      <- indexI_below  
  paramsBellman$indexNI_below     <- indexNI_below 
  paramsBellman$indexI_above      <- indexI_above  
  paramsBellman$indexNI_above     <- indexNI_above 
  paramsBellman$states            <- states        
  paramsBellman$indexVioCase      <- indexVioCase  
  paramsBellman$CaseProb          <- CaseProb      
  paramsBellman$insCase           <- insCase       
  paramsBellman$vioCase           <- vioCase       
  paramsBellman$fineCase          <- fineCase      
  paramsBellman$HPVCase           <- HPVCase
  paramsBellman$indexComp         <- indexComp
  
  return(paramsBellman)
}

BellmanLoop <- function(coeff, params,V0){
  # Estimates the fixed point of the Bellman equation
  
  # Retrieve Bellman parameters
  paramsBellman <- params$BellmanParams
  
  # Unpack Bellman parameters
  NewV <- paramsBellman$NewV
  OldV <- paramsBellman$OldV
  w_below <- paramsBellman$w_below
  w_above <- paramsBellman$w_above
  indexI_below <- paramsBellman$indexI_below
  indexNI_below <- paramsBellman$indexNI_below
  indexI_above <- paramsBellman$indexI_above
  indexNI_above <- paramsBellman$indexNI_above
  states <- paramsBellman$states
  indexVioCase <- paramsBellman$indexVioCase
  CaseProb <- paramsBellman$CaseProb
  insCase  <- paramsBellman$insCase
  vioCase  <- paramsBellman$vioCase
  fineCase <- paramsBellman$fineCase
  HPVCase  <- paramsBellman$HPVCase
  indexComp <- paramsBellman$indexComp
  beta <- params$beta                  # discount factor
  gamma <- params$gamma                # euler's constant
  n_omega1 <- params$n_omega1          # number of omega1 states
  N1 <- params$N1                      # number of states for each omega1
  nstates <- n_omega1*N1               # Number of states
  n_regoutcomes <- params$nregoutcomes # Number of regulatory outcomes (80)
  n_transitions <- params$ntransitions # violator transtions (3)
  
  if(!missing(V0)) {
    OldV <- V0
  } 
  
  # Unpack coefficient (theta)
  theta_X <- coeff[1]
  theta_I <- coeff[2]
  theta_V <- coeff[3]
  theta_F <- coeff[4]
  theta_H <- coeff[5]
  
  # Per period utility
  U <- theta_I*insCase + theta_V*vioCase + theta_F*fineCase +
    theta_H*HPVCase
  U <- c(as.matrix(U))
  
  # Loop until NewV converges to the fixed point
  
  iconv = 0 # counts the number of iterations
  norm_Vdiff <- 1
  while (norm_Vdiff>=params$tol){
    iconv = iconv + 1
    
    # 1) Update Vtilde and investment probabilities ############################

    # Solve OldV if investment
    OldV_I <- w_below*OldV[indexI_below] + w_above*OldV[indexI_above]
    
    # Solve OldV if no investment
    OldV_NI <- w_below*OldV[indexNI_below] + w_above*OldV[indexNI_above]
    
    # Use logit value to find Vtilde
    V_C  <- gamma + beta*OldV[indexComp] # value if in compliance
    V_NC <- gamma + log(exp(beta*OldV_NI)+exp(-theta_X+beta*OldV_I)) # not in compliance
    
    # Update Vtilde
    Vtilde <- (states[,4] == 0)*V_C + (states[,4] > 0)*V_NC # states[,4] -> ordered violator status
  
    # 2) Update NewV ###########################################################
    
    VCase <- (U+Vtilde[indexVioCase])*CaseProb
    VCase <- matrix(VCase,nstates,n_regoutcomes*n_transitions)
    NewV <- rowSums(VCase)
    
    # 3) Calculate norm of diff between OldV and NewV ##########################
    norm_Vdiff <- max(abs(OldV-NewV))
    OldV <- NewV
    #print(paste0("Iteration # ",iconv,", norm = ",norm_Vdiff))
  }
  
  # 4) Pack and return main results ############################################
  results <- list()
  
  # Save prob. of investment
  Investprob <- ((states[,4] > 0))*
    exp(-theta_X+beta*OldV_I)/(exp(beta*OldV_NI)+exp(-theta_X+beta*OldV_I))
  
  # If V contains NA (exploding case) just return 1
  if (any(is.na(NewV))){
    NewV <- paramsBellman$OldV+1
    Vtilde <- paramsBellman$OldV+1
    Investprob <- paramsBellman$OldV*0
  }
  
  results$NewV <- NewV
  results$Vtilde <- Vtilde
  results$Investprob <- Investprob
  
  return(results)
}

CompareValues <- function(outBellman,file,params){
  # Function to check value function with provided file
  
  library(ggplot2)
  library(tidyr)
  
  N1 <- params$N1
  n_omega1 <- params$n_omega1
  nstates <- n_omega1*N1
  myResults <- outBellman$NewV
  
  expectedResults <- read.csv(file)
  expectedResults <- expectedResults[1:nstates,"value"]
  
  diff <- max(abs(myResults-expectedResults))
  
  print(paste0("the difference is ",diff))
  
  
  index <- c(1:nstates)
  plotData <- data.frame(myResults,expectedResults,index)
  
  plotData %>%
    gather(key,value, myResults, expectedResults) %>%
    ggplot(aes(x=index, y=value, colour=key)) +
    geom_line()
  
}

CompareVTilde <- function(outBellman,file,params){
  # Function to check vtilde function with provided file
  
  library(ggplot2)
  library(tidyr)
  
  N1 <- params$N1
  n_omega1 <- params$n_omega1
  nstates <- n_omega1*N1
  myResults <- outBellman$Vtilde
  
  expectedResults <- read.csv(file)
  expectedResults$compliance <- 1*(expectedResults$orderedvio>0)
  expectedResults <- expectedResults[order(expectedResults$NAICS,
                                           expectedResults$region,
                                           expectedResults$gravity,
                                           expectedResults$compliance,
                                           expectedResults$laginv,
                                           expectedResults$violation,
                                           expectedResults$orderedvio,
                                           expectedResults$DAVgrid),]
  expectedResults <- expectedResults[1:nstates,"value"]
  
  diff <- max(abs(myResults-expectedResults))
  
  print(paste0("the difference is ",diff))
  
  
  index <- c(1:nstates)
  plotData <- data.frame(myResults,expectedResults,index)
  
  plotData %>%
    gather(key,value, myResults, expectedResults) %>%
    ggplot(aes(x=index, y=value, colour=key)) +
    geom_line()
  
}

CompareProb <- function(LLdata,file){
  # Function to check investment probabilites and likelihood
  library(ggplot2)
  library(tidyr)
  
  LLdata$obs2 <- c(1:nrow(LLdata))
  LLdata <- LLdata[LLdata$viostatus>0,]
  
  
  expectedResults <- read.csv(file)
 
  merged <- cbind(LLdata,expectedResults)
  tol <- 1e-3
  #merged <- merge(x = expectedResults, y = LLdata, by = "obs", all.x = TRUE)
  merged$probdiff <- abs(merged$invprob - merged$prob)>tol
  merged <- merged[,c("obs","obs2","o1","li1","li2","viostatus","violation","DAV","x",
                      "invprob","prob","probdiff","loglike","l_i")]
  
  merged$DAV0 <- 1*(merged$DAV == 0)
  nobs <- nrow(merged)
  ndiff <- sum(merged$probdiff)             
  
  
  diff <- max(abs(merged$invprob-merged$prob))
  
  print(paste0(ndiff," observations out of ",nobs," have different prob. values"))
  print(paste0("the difference in probailities is ",diff))
  diff <- max(abs(merged$loglike-merged$l_i))
  print(paste0("the difference in log-likelihood is ",diff))
  
  p<-ggplot(merged,aes(x=invprob, y=prob)) +
    geom_point()+xlab("Investment probability (expected result)")+
    ylab("Investment probability (my result)") + xlim(0,1)+ylim(0,1)
  show(p)
  return(merged)
}

Question1C <- function(outBellman,file){
  # Recover the liklihood for states with DAV = 2
  index <- c(1:length(outBellman$Investprob))
  states <- t(sapply(index,IndexToValues)) 
  dataout <- data.frame(states)
  colnames(dataout) <- c("omega1","lag1_inv","violation","viostatus","intDAV")
  dataout$DAV <- (dataout$intDAV-1)*0.5 
  dataout$invprob <- outBellman$Investprob
  dataout <- dataout[dataout$DAV==2,]
  write.csv(dataout,file, row.names = FALSE)
}

AddParams <- function(params){
  # Add BellmanParams and LogLikeParams for faster computations
  params$BellmanParams <- PackBellmanParams(params)
  params$LogLikeParams <- LogLikeParams(params)
  return(params)
}

EstimateMLVariance <- function(thetaML, params,DeltaTheta,V0){
  # Estimates the variance of the ML estimator usin the outer product approximation
  
  # Construct matrices to compute the Bellman equation
  if (!exists("BellmanParams",where = params)){
    params$BellmanParams <- PackBellmanParams(params)
  }
  
  # Estimate likehood at thetaML
  if (missing(V0)){
    outML <- BellmanLoop(thetaML, params)  
  }
  else {
    outML <- BellmanLoop(thetaML, params,V0)
  }
  
  LLOutput <- LogLike(outML,params)
  LLData <- LLOutput$lldata
  LLData$LL_ML <- LLData$l_i  
  
  # Estimate likelihood at deviations
  for (i in 1:length(thetaML)){
    theta_i <- thetaML
    theta_i[i] <- theta_i[i] - DeltaTheta # Delta-deviation
    out <- BellmanLoop(theta_i, params,outML$NewV)
    LLOutput <- LogLike(out,params)
    varname <- paste0("dlog_",i)
    LLData[,varname] <- LLData$LL_ML-LLOutput$lldata$l_i
  }
  
  # Compute the outer product of the gradient of the log likelihood function
  GradProd <- matrix(0,length(thetaML),length(thetaML))
  Ncols <- ncol(LLData)
  LLData<- as.matrix(LLData[,c((Ncols + 1 - length(thetaML)):Ncols)])
  
  for (i in 1:nrow(LLData)){
    GradProd = GradProd + (t(t(LLData[i,]))%*%LLData[i,])/(DeltaTheta^2)
  }
  return(solve(GradProd))
}
 
LogLike <- function(BellmanOut,params){
  # Computes the Log likelihood at all observations. Also returns investment
  # probabilities
  
  if (!exists("LogLikeParams",where = params)){
    params$LogLikeParams <- LogLikeParams(params)
  }
  
  i1 <- params$LogLikeParams$i1
  i2 <- params$LogLikeParams$i2
  w_below <- params$LogLikeParams$w_below
  w_above <- params$LogLikeParams$w_above
  x <- params$LogLikeParams$x
  o1 <- params$LogLikeParams$o1
  li1 <- params$LogLikeParams$li1
  li2 <- params$LogLikeParams$li2
  viostatus <- params$LogLikeParams$viostatus
  violation <- params$LogLikeParams$violation
  DAV <- params$LogLikeParams$DAV

  p1 <- BellmanOut$Investprob[i1]
  p2 <- BellmanOut$Investprob[i2]
  
  prob <- w_below*p1 + w_above*p2
  
  l_i <- log(x*prob + (1-x)*(1-prob))
  ll <- sum(l_i)

  out <- list()
  out$ll <- ll
  
  lldata <- data.frame(o1,li1,li2,viostatus,violation,DAV,x,prob,l_i)
  
  out$lldata <- lldata
  return(out)
}


LogLikeParams <- function(params){
  # Defines matrices that allow faster repeated computation of the likelihood
  
  DAVgrid <-params$DAVgrid
  data <-params$panel_data
  n_obs <-nrow(data)
  
  DAVgrid <- params$DAVgrid
  npointsDAV  <- length(DAVgrid)
  maxpointDAV <- DAVgrid[npointsDAV]
  widthDAV <- DAVgrid[2]-DAVgrid[1]
  
  DAV <- data$DAV
  x   <- data$investment
  o1  <- data$omega1
  li1 <- data$lag_investment
  li2 <- data$lag2_investment
  viostatus <- data$ordered_violator
  violation <- data$violation
  
  # Interpolate newDAV using the grid
  # indexes in DAV grid
  i_below <- pmin(1 + floor(DAV*(npointsDAV-1)/maxpointDAV),npointsDAV)
  i_above <- pmin(2 + floor(DAV*(npointsDAV-1)/maxpointDAV),npointsDAV)
  # weights
  w_below <- pmax((DAVgrid[i_above]-DAV)/widthDAV,0)
  w_above <- 1-w_below
  # states
  state1 <- cbind(o1,li1,violation,viostatus,i_below)
  state2 <- cbind(o1,li1,violation,viostatus,i_above)
  # indexes in probability matrix
  i1 <- apply(state1,1,ValuesToIndex)
  i2 <- apply(state2,1,ValuesToIndex)
  
  llparams <- list()
  llparams$i1 <- i1
  llparams$i2 <- i2
  llparams$w_below <- w_below
  llparams$w_above <- w_above
  llparams$x <- x
  llparams$o1 <- o1
  llparams$li1 <- li1
  llparams$li2 <- li2
  llparams$viostatus <- viostatus
  llparams$violation <- violation
  llparams$DAV <- DAV
  
  
  return(llparams)
}

LogLike2 <- function(BellmanOut,params){
  # Faster computation of the total likelihood

  i1 <- params$LogLikeParams$i1
  i2 <- params$LogLikeParams$i2
  w_below <- params$LogLikeParams$w_below
  w_above <- params$LogLikeParams$w_above
  x <- params$LogLikeParams$x
  
  p1 <- BellmanOut$Investprob[i1]
  p2 <- BellmanOut$Investprob[i2]
  
  prob <- w_below*p1 + w_above*p2
  
  l_i <- log(x*prob + (1-x)*(1-prob))
  ll <- sum(l_i)
  return(ll)
}


EvalFunction <- function(theta,params){
  # Function to be maximized in ML estimation. First computes the solution to the
  # Bellman equation, then computes the log-likelihood
  
  print(theta)
  
  if (exists("V0")){
    outBellman <- BellmanLoop(theta, params,V0)
  }
  else {
    outBellman <- BellmanLoop(theta, params)
  }
  assign("V0",outBellman$NewV,envir = .GlobalEnv)
  LL <- LogLike2(outBellman,params)
  return(LL)
}

NestedFixedPoint<-function(theta0,params,solver){
  # Nested fixed point estimation
  
  # Construct matrices to compute the Bellman equation
  if (!exists("BellmanParams",where = params)){
    params$BellmanParams <- PackBellmanParams(params)
  }
  
  if (!exists("LogLikeParams",where = params)){
    params$LogLikeParams <- LogLikeParams(params)
  }
  
  # Solver selection
  if (solver== 1){
    library(pracma)
    out<-fminsearch(EvalFunction, theta0, params = params,
                    minimize= FALSE, method = "Hooke-Jeeves", maxiter = 500, tol = 1e-6)
  }
  else if (solver == 2){
    library(pracma)
    out<-fminsearch(EvalFunction, theta0, params = params,
                    minimize= FALSE, method = "Nelder-Mead", maxiter = 500, tol = 1e-6)
  }
  else if (solver == 3){
    library(optimx)
    out <- optimx(theta0, EvalFunction, gr=NULL, hess=NULL, lower=c(0,-Inf,-Inf,-Inf,-Inf),
                  upper=c(Inf,0,0,0,0), 
           method="CG", itnmax=500, hessian=FALSE,
           control=list(maximize=TRUE),
           params = params)
  }
  else if (solver == 4){
    library(optimx)
    out <- optimx(theta0, EvalFunction, gr=NULL, hess=NULL, lower=c(0,-Inf,-Inf,-Inf,-Inf),
                  upper=c(Inf,0,0,0,0), 
                  method="L-BFGS-B", itnmax=500, hessian=FALSE,
                  control=list(maximize=TRUE),
                  params = params)
  }
  else{
    library(maxLik)
    
    A <- matrix(c(0,0,0,0,-1), 1, 5)
    B <- 0
    out<- maxNR(EvalFunction, grad = NULL, hess = NULL, start = theta0, print.level = 1,
                tol = 1e-06, reltol=sqrt(.Machine$double.eps), gradtol = 1e-06,
                steptol = 1e-10, lambdatol = 1e-06, qrtol = 1e-10,
                iterlim = 500,
                constraints=list(eqA=A, eqB=B),
                params = params)
  }
  
  return(out)
}

IndexToValues <- function(index){
  # Function that takes an index of the V matrix and return the corresponding
  # omega (or omega_tilde) states
  
  N1 <- 161
  N2 <- 80
  N3 <- 40
  N4 <- 20
  
  omega1 <- 1 + (index-1)%/%N1
  index <- (index-1)%%N1-1 # index goes from -1 to N1-2
  violator <- index>-1
  lagi1 <- (index%/%N2)*violator
  index <- index%%N2  # index goes from 0 to N2-1
  lagi2 <- (index%/%N3)*violator
  index <- index%%N3
  vio   <- (index%/%N4+1)*violator
  index <- index%%N4+1
  dav   <- index*violator + (violator==0)
  return(c(omega1,lagi1,lagi2,vio,dav))
}

ValuesToIndex <- function(x){
  # Function that takes an omega state and returns its corresponding index in the
  # V matrix
  
  N1 <- 161
  N2 <- 80
  N3 <- 40
  N4 <- 20
  
  i_o1  <- x[1]
  i_l1  <- x[2]
  i_l2  <- x[3]
  i_vio <- x[4]
  i_DAV <- x[5]
  
  i <- (i_o1-1)*N1 + 1 + (i_vio>=1)*((i_l1)*N2 + (i_l2)*N3 + (i_vio-1)*N4 + i_DAV)
  return(i)
}

