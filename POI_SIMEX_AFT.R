# Functions used to simulate and fit an Accelerated Failure time (AFT) model using
#  the POI_SIMEX method.  The package 'simexaft' functions were modified for the
#  conditional Poisson distributed covariates.

# REFERENCE:   Xiong J, He W, Yi GY (2019). simexaft: R package version 1.0.7.1,
#    <https://CRAN.R-project.org/package=simexaft>.

# REFERENCE:  "POI-SIMEX for conditionally 
#    Poisson distributed biomarkers from tissue microarrays" (2024)
#    Aijun Yang (aijunyan@uvic.ca), Phineas T. Hamilton, Brad H. Nelson, Julian J. Lum, 
#    Mary Lesperance, Farouk S. Nathoo

sim.lognormT<-function(a,b,n,unif.a0=0.5,unif.b0=9,seed=12345,rcensor=0){
 # Simulate data from log-normal AFT model with:
  #   intercept b0; one covariate measured without error, Z~unif(unif.a0, unif.b0), with coef b2;
  #   one covariate measured with error ~ Poisson(true mean = lambda*A) with coef b1, 
  #    with lambda~gamma(shape=a, scale=b), and A = area.
  #   The random normal component has standard deviation sig.e.
  # Arguments
  # a, b are gamma shape and scale parameters of distribution for true covariate, lambda
  # n, the number of observations
  # b0, b1, b2 are the model coefficients
  # sig.e is standard deviation of random normal
  # unif.a0, unif.b0 are two parameters, Z~unif(unif.a0, unif.b0)
  # seed: random seed
  # rcensor: % of censoring, values in [0, 1)
  # output dataframe: simdata, 
  
  set.seed(seed)
  simdata<-c()
  #area in POI process--set to 1 for simplicity
  A<-rep(1,n)
  ##covariate without ME
  Z<-runif(n,unif.a0,unif.b0)
  
  ##a=gamma shape b=gamma scale
  
  lambda<-rgamma(n,shape=a,scale=b)  #rate=1/scale
    
    W<- rpois(n, lambda*A) #generate Poissons, mis-measured covariate
    lambda.naive<-W/A
    
    log.Y<-b0+b1*lambda+b2*Z+sig.e*rnorm(n)  ##set standard normal dist rnorm(0,1)
    Y<-exp(log.Y)
    delta<- rep(1, n) #event indicator
    delta[sample(n, rcensor*n)] <- 0 # random censorship
    
    simdata<-data.frame(Y=Y,Z=Z,lambda=lambda,lambda.naive=lambda.naive,A=A,delta=delta)
  return (simdata)
}

# REFERENCE:   Xiong J, He W, Yi GY (2019). simexaft: R package version 1.0.7.1,
#    <https://CRAN.R-project.org/package=simexaft>.

linearextrapolation <-
  function(A1,A2,A3,lambda){
   # A1 estimates obtained from each level of labmda.
   # A2 variances estimates obtained from each level of lambda.
   # A3 scale estimates obtained from each level of lambda.
   # lambda vector of lambdas, the grids for the extrapolation step.
    
    #extrapolation back to lambda=-1 yield the SIMEX estimates
    reg1<-numeric()
    
    #extrapolation back to lambda=-1 yield the SIMEX estimates of variances
    reg2<-numeric() 
    
    #extrapolation back to lambda=-1 yield the SIMEX estimates of scale
    scalereg<-numeric()
    
    D<-ncol(A1)                         
    
    #regression on estimates,tau#
    for(i in 1:D)
    {    
      e1<-coef(lm(A1[,i]~lambda))
      a1<- e1[1] - e1[2]
      reg1 <- c(reg1,a1)
      
      e2 <-coef(lm(A2[,i]~lambda))        
      a2 <- e2[1] - e2[2]        
      reg2 <- c(reg2, a2)
      
    }
    
    #regression on scale#
    e3 <- coef(lm(A3[,1]~lambda))
    a3 <-  e3[1] - e3[2]
    scalereg<- c(scalereg,a3)
    
    return(list("reg1"=reg1,"reg2"=reg2,"scalereg"=scalereg))
  }


quadraticextrapolation <-
  function(A1,A2,A3,lambda){
   # A1 estimates obtained from each level of labmda.
   # A2 variances estimates obtained from each level of lambda.
   # A3 scale estimates obtained from each level of lambda.
   # lambda vector of lambdas, the grids for the extrapolation step.
    
    #extrapolation back to lambda=-1 yield the SIMEX estimates
    reg1<-numeric()
    
    #extrapolation back to lambda=-1 yield the SIMEX estimates of variances
    reg2<-numeric() 
    
    #extrapolation back to lambda=-1 yield the SIMEX estimates of scale
    scalereg<-numeric()
    
    D<-ncol(A1)
    
    #regression on estimates,tau#
    for(i in 1:D)
    {                      
      lambda2<-lambda^2                                
      e1<-coef(lm(A1[,i]~lambda + lambda2))
      a1 <- e1[1] - e1[2] + e1[3]
      reg1 <- c(reg1,a1)
      
      e2 <- coef(lm(A2[,i]~lambda + lambda2))
      a2 <- e2[1] - e2[2] + e2[3]
      reg2 <- c(reg2, a2)                                     
    }
    
    
    #regression on scale#
    e3<-coef(lm(A3[,1]~lambda+lambda2))
    a3<- e3[1]-e3[2]+e3[3]
    scalereg<- c(scalereg,a3)
    
    return(list("reg1"=reg1,"reg2"=reg2,"scalereg"=scalereg))
    
  }

POI.simexaft <-
  function(formula=formula(data),data=parent.frame(),SIMEXvariable,areaVariable,B=50,lambda=seq(0,2,0.1),extrapolation="quadratic",dist="lognormal")
  {
   # formula --the model used to fit the observed covariate with measurement error
   # data --input data, must be in data.frame format
   # SIMEXvariable, covariate with measurement error and need SIMEX for bias correction
   # areaVariable,  area or length corresponding to the conditional Poisson distributed covariate
   # B, the number of simulated data containing additional measurement error for each lambda value
   # lambda, grid of values in POI-SIMEX simulation step
   # extrapolation, indicating extrapolation method, linear or quadratic
   # dist, distribution used in AFT model, current implmentation is only for lognormal
    
    ############################ check the input of SIMEXvariable #######################################
    colname<-colnames(data)
    SIMEXvariable<-unique(SIMEXvariable)
    nSIMEXvariable<-length(SIMEXvariable)
    
    xbar<-mean(data[,SIMEXvariable])
    ##divided by area
    sig<-xbar/data[,areaVariable]
    
    
    if(!is.character(SIMEXvariable) | nSIMEXvariable>length(colname)){
      
      stop("Invalid SIMEXvariable object")
      
    }
    if(!all(SIMEXvariable %in% colname)){
      
      stop("SIMEXvariable must selected from the data")
      
    }
    
    
    ############################  lambda ################################
    if(!is.vector(lambda) |!is.numeric(lambda)){
      
      stop(":Invalide lambda object")
      
    }
    
    if (any(lambda < 0)) {
      
      warning("Lambda should be positive values. Negative values will be ignored",call. = FALSE)
      
      lambda <- lambda[lambda >= 0]
      
    }
    
    
    ##### specifies the regression form used in the extrpolation step #####
    extrapolation = substr(extrapolation, 1, 4)
    
    if(!is.character(extrapolation) | length(extrapolation)!=1){
      
      warning("Invalid extrapolation object. Using: quadatic\n\n",call.=FALSE)
      
      extrapolation<-"quad"
    }
    
    
    ##############################  prepare reading data from data.frame############
    
    
    ndata<-nrow(data)

    # nformula=length(attr(terms(formula),"term.labels"))+1
    # update the above code for the covariate having more than 2 levels
    
    ###fit the naive model
    ncoef<-survreg(formula=formula,data=data,dist=dist,robust=T)$coefficients
    nformula<-length(ncoef)
    
    nlambda=length(lambda)
    
    #the matrixes to save the estimates of the coefficients, variance #
    A1<-matrix(data=NA,nlambda,nformula)
    
    A2<-matrix(data=NA,nlambda,nformula)
    
    A3<-matrix(data=NA,nlambda,nformula)
    
    #the matrix of estimate#
    theta<-matrix(data=NA,B,nformula)
    
    colnames(theta)<-names(ncoef)
    
    p.names<-colnames(theta)
    
    #all estimates corresponding to each lambda#
    theta.all<-vector(mode="list",nlambda)
    
    
    #####################use SIMEX method to analysis the data#########
    k<-1
    while(k <= length(lambda))
    {
      ## the coefficients estimates of the kth sample##
      w<-numeric()
      
      ##the variance of the coefficient estimate of the kth sample##
      v<-numeric()
      
      ##the variance estimate of the kth sample##
      omega<-numeric()
      
      temp<-data
      
      #the matrix of variance estimate#
      estivarB<-matrix(data=NA,B,nformula)
      
      #the vector of the estimate scale#
      estiscaleB<-matrix(data=NA,B,ncol=1)
      
      ### simulation step of the SIMEX algorithm ###
      for(r in 1:B)
      { 
        
        temp[SIMEXvariable]<-data[SIMEXvariable]+sqrt(lambda[k]*sig)*rnorm(ndata,0,1)
   
        ### use survreg to fit the artifical simulated data to the the AFT model ###
        re <- survreg(formula=formula,data=temp,dist=dist,robust=T)
        
        ##obtain the coefficients estimates w, associated variance estimate omega and the scale estimate scale##
        scale<-re$scale
        w<-re$coefficients
        omega<-diag(re$var)[1:nformula]
        theta[r,]<-w
        estivarB[r,]<-omega
        estiscaleB[r,]<-scale
      }
      
      ###for each lambda value calculate the mean of the each estimates of B samples ###
      ##the coefficients estimates and the variance of the estiamtes##
      w<-apply(theta,2,FUN=mean)
      v<-apply(theta,2,FUN=var)
      
      
      ##the variance estimates##
      omega <-apply(estivarB,2,FUN=mean)
      
      
      ##the difference of the variance estimates and the variance of the coefficient estimates, use this to estimate the variance##
      tau<-omega-v
      
      ### remove the sample with negative tau###
      if(all(tau>0)){
        A1[k,] <- w
        A2[k,] <- tau
        A3[k,] <- apply(estiscaleB,2,FUN=mean)
        theta.all[[k]]<-theta
        k <- k + 1
      }

    }
    
    ### save all the estimates in the theta.all###
    theta<-matrix(unlist(theta.all),nrow=B)
    theta.all<-list()
    
    for (i in 1:nformula){
      
      theta.all[[p.names[i]]] <- data.frame(theta[,seq(i, nformula * nlambda, by = nformula)])
    }
    
    
    #####################extrapolation step of the SIMEX algorithm, fit the regression model#####################
    if(extrapolation=="line"){
      
      result1<-linearextrapolation(A1,A2,A3,lambda)
      
      
    }
    
    
    else if(extrapolation=="quad"){
      
      result1<-quadraticextrapolation(A1,A2,A3,lambda)
      
    }
    
    
    else stop("extrapolation method must be linear or quadratic")
    
    
    ####### outputs###########################
    estimate<-result1$reg1
    names(estimate)<-p.names
    se<-sqrt(result1$reg2)
    names(se)<-p.names
    scalereg<-result1$scalereg
    pvalue<-2*(1-pnorm(abs(estimate/se)))

    ##drop repind=repind in the output list as just for one replicate
    erg<-list(coefficients=estimate,se=se,scalereg=scalereg,pvalue=pvalue,lambda=lambda, B=B, formula=formula, extrapolation=extrapolation,SIMEXvariable=SIMEXvariable,theta=theta.all)

    return(erg)
    
  }
