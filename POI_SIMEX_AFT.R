
linearextrapolation <-
  function(A1,A2,A3,lambda){
    
    #the fitted value of linear extrapolation at lambda=-1#
    reg1<-numeric()
    
    #the fitted value of linear extrapolation on diff at lambda=-1#
    reg2<-numeric() 
    
    #the fitted value of linear extrapolation on scale at lambda=-1#
    scalereg<-numeric()
    
    D=ncol(A1)                         
    
    #regression on estimates,tau#
    for(i in 1:D)
    {    
      e1=coef(lm(A1[,i]~lambda))
      a1= e1[1] - e1[2]
      reg1 = c(reg1,a1)
      
      e2 = coef(lm(A2[,i]~lambda))        
      a2 = e2[1] - e2[2]        
      reg2 = c(reg2, a2)
      
    }
    
    #regression on scale#
    e3 = coef(lm(A3[,1]~lambda))
    a3 =  e3[1] - e3[2]
    scalereg= c(scalereg,a3)
    
    return(list("reg1"=reg1,"reg2"=reg2,"scalereg"=scalereg))
  }


quadraticextrapolation <-
  function(A1,A2,A3,lambda){
    #the fitted value of quadratic extrapolation at lambda=-1#
    reg1<-numeric()
    
    #the fitted value of quadratic extrapolation on diff at lambda=-1#
    reg2<-numeric() 
    
    #the fitted value of quadratic extrapolation on scale at lambda=-1#
    scalereg<-numeric()
    
    D=ncol(A1)
    
    #regression on estimates,tau#
    for(i in 1:D)
    {                      
      lambda2=lambda^2                                
      e1=coef(lm(A1[,i]~lambda + lambda2))
      a1 = e1[1] - e1[2] + e1[3]
      reg1 = c(reg1,a1)
      
      e2 = coef(lm(A2[,i]~lambda + lambda2))
      a2 = e2[1] - e2[2] + e2[3]
      reg2 = c(reg2, a2)                                     
    }
    
    
    #regression on scale#
    e3=coef(lm(A3[,1]~lambda+lambda2))
    a3= e3[1]-e3[2]+e3[3]
    scalereg= c(scalereg,a3)
    
    return(list("reg1"=reg1,"reg2"=reg2,"scalereg"=scalereg))
    
  }

POI.simexaft <-
  function(formula=formula(data),data=parent.frame(),SIMEXvariable,areaVariable,repind=list(),B=50,lambda=seq(0,2,0.1),extrapolation="quadratic",dist="weibull")
  {
    #for debug
    #data=bio_surv_fm;SIMEXvariable="bio";areaVariable="area";repind=list();B=50;lambda=seq(0,2,0.1);extrapolation="quad";dist="lognormal";
    
 
    ############################ check the input of SIMEXvariable #######################################
    colname=colnames(data)
    SIMEXvariable=unique(SIMEXvariable)
    nSIMEXvariable=length(SIMEXvariable)
    
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
      
      extrapolation="quad"
    }
    
    
    ##############################  prepare reading data from data.frame############
    
    
    ndata=nrow(data)
    
    ###fit the naive model
    ncoef<-survreg(formula=formula,data=data,dist=dist,robust=T)$coefficients
    nformula=length(ncoef)
    
    #nformula=length(attr(terms(formula),"term.labels"))+1
    
    nlambda=length(lambda)
    
    #the matrixes to save the estimates of the coefficents, variance #
    A1=matrix(data=NA,nlambda,nformula)
    
    A2=matrix(data=NA,nlambda,nformula)
    
    A3=matrix(data=NA,nlambda,nformula)
    
    #the matrix of estimate#
    theta=matrix(data=NA,B,nformula)
    
    #colnames(theta)=c("Intercept",attr(terms(formula),"term.labels"))
    
    colnames(theta)=names(ncoef)
    
    p.names=colnames(theta)
    
    #all estimates corresponding to each lambda#
    theta.all=vector(mode="list",nlambda)
    
    
    #####################use SIMEX method to analysis the data#########
    k=1
    while(k <= length(lambda))
    {
      ## the coefficients estimates of the kth sample##
      w=numeric()
      
      ##the variance of the coefficent estimate of the kth sample##
      v=numeric()
      
      ##the variance estimate of the kth sample##
      omega=numeric()
      
      
      temp=data
      
      #the matrix of variance estimate#
      estivarB=matrix(data=NA,B,nformula)
      
      
      #the vector of the estimate scale#
      estiscaleB=matrix(data=NA,B,ncol=1)
      
      ### simulation step of the SIMEX algorithm ###
      for(r in 1:B)
      { 
        
        temp[SIMEXvariable]=data[SIMEXvariable]+sqrt(lambda[k]*sig)*rnorm(ndata,0,1)
   
        ### use survreg to fit the artifical simulated data to the the AFT model ###
        re = survreg(formula=formula,data=temp,dist=dist,robust=T)
        
        ##obtain the coefficients estimates w, associated variance estimate omega and the scale estimate scale##
        scale=re$scale
        w=re$coefficients
        omega=diag(re$var)[1:nformula]
        theta[r,]=w
        estivarB[r,]=omega
        estiscaleB[r,]=scale
      }
      
      ###for each lambda value calculate the mean of the each estimates of B samples ###
      ##the coefficients estimates and the variance of the estiamtes##
      w=apply(theta,2,FUN=mean)
      v=apply(theta,2,FUN=var)
      
      
      ##the variance estimates##
      omega =apply(estivarB,2,FUN=mean)
      
      
      ##the difference of the variance estimates and the variance of the coefficient estimates, use this to estimate the variance##
      tau=omega-v
      
      ### remove the sample with negative tau###
      if(all(tau>0)){
        A1[k,] = w
        A2[k,] = tau
        A3[k,] = apply(estiscaleB,2,FUN=mean)
        theta.all[[k]]=theta
        k = k + 1
      }

    }
    
    ### save all the estimates in the theta.all###
    theta=matrix(unlist(theta.all),nrow=B)
    theta.all=list()
    
    for (i in 1:nformula){
      
      theta.all[[p.names[i]]] <- data.frame(theta[,seq(i, nformula * nlambda, by = nformula)])
    }
    
    
    #####################extrapolation step of the SIMEX algorithm, fit the regression model#####################
    if(extrapolation=="line"){
      
      result1=linearextrapolation(A1,A2,A3,lambda)
      
      
    }
    
    
    else if(extrapolation=="quad"){
      
      result1=quadraticextrapolation(A1,A2,A3,lambda)
      
    }
    
    
    else stop("extrapolation method must be linear or quadratic")
    
    
    ####### outputs###########################
    estimate=result1$reg1
    names(estimate)=p.names
    se=sqrt(result1$reg2)
    names(se)=p.names
    scalereg=result1$scalereg
    pvalue=2*(1-pnorm(abs(estimate/se)))
    
    erg=list(coefficients=estimate,se=se,scalereg=scalereg,pvalue=pvalue,lambda=lambda, B=B, formula=formula, extrapolation=extrapolation,SIMEXvariable=SIMEXvariable,repind=repind,theta=theta.all)
      
    #class(erg)<-("simexaft")
    
    return(erg)
    
  }