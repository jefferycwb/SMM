pxt <- function(qxt, valyear, ageStart, omegaAge){
  t_start <- which(colnames(qxt) == toString(valyear))
  x_start <- which(rownames(qxt) == toString(ageStart))
  qxt <- diag(qxt[x_start:(x_start+(omegaAge-ageStart)), t_start:(t_start+(omegaAge-ageStart))])
  pxt <- c(1, head(cumprod(1- qxt), -1))
  pxt
}

createActuarialTable<-function(yearOfBirth,rate){
  qxcoh <- rate[1:nrow(rate), (yearOfBirth-min(years)+1):ncol(rate)]
  cohort.qx <- diag(qxcoh)
  #get projected Px
  fittedQx = cohort.qx #add px to table
  qx4Completion = seq(from = cohort.qx[length(fittedQx)], to = 0, length = 20)
  totalQx = c(fittedQx,qx4Completion[2:length(qx4Completion)])
  #create life table
  irate=1.04/1.02-1
  
  cohortLt=probs2lifetable(probs=totalQx, radix=100000,type="qx",
                           name=paste("Cohort",yearOfBirth))
  cohortAct=new("actuarialtable",x=cohortLt@x, lx=cohortLt@lx,
                interest=irate, name=cohortLt@name)
  return(cohortAct)
}

# Kannisto extrapolation of mortality rates for old ages
kannistoExtrapolation = function(qx, ages, years, max_age=120, nObs=15){
  ## calculate ages to be extrapolated
  agesExtrap <- (tail(ages,1)+1):max_age
  
  ## extract the ages to be used for the regression
  obsAges <- tail(ages, nObs)
  
  #extrapolate for each year
  extrapolateQxt <- sapply(years, function(x){
    ## extract the mortality rates and ages which are going to be used for the regression
    obsQx <- pmin(tail(qx[,toString(x)], nObs),1)
    obsQxt  <- obsQx 
    ## transform the mortality rates
    obsQx <- -log(1 - obsQx)
    obsQx <- log(obsQx/(1 - obsQx))
    
    if (length(na.omit(obsQx))>2){
      ## apply the regression and obtain phi1 and phi2
      model <- lm(formula = obsQx~obsAges)
      phi1 <- exp(model$coefficient[1])
      phi2 <- model$coefficient[2]
      
      #extrapolate mortality rates and add to the overall matrix
      1-exp(-phi1*exp(phi2*agesExtrap)/(1+phi1*exp(phi2*agesExtrap)))
    }else{
      seq(obsQxt[nObs], to = 1, length = 30)
    }
  })
  
  qxt <- rbind(qx,  extrapolateQxt)
  ages <- c(ages, agesExtrap)
  
  list(qxt = qxt, ages.fit = ages)
}


DFcashflow = function(qxt, ageStart, omegaAge, pensionAge, valyear, ir, type = 1){
  if (type == 1){
    CF <- c(rep(0,max(0,(pensionAge-ageStart+1))), rep(1,(omegaAge-max(pensionAge, ageStart))))
  }else{
    CF <- c(rep(1,max(0,(pensionAge-ageStart+1))), rep(0,(omegaAge-max(pensionAge, ageStart))))
  }
  
  t_start <- which(colnames(qxt) == toString(valyear))
  x_start <- which(rownames(qxt) == toString(ageStart))
  qx      <- diag(qxt[x_start:(x_start+(omegaAge-ageStart)), t_start:(t_start+(omegaAge-ageStart))])
  pxt <- cumprod(1- qx)
  k <- 0:(omegaAge-ageStart)
  interestRates <- rep(ir, (omegaAge-ageStart)+1)
  v <- (1 + interestRates)^-k
  dfCF <- sum(pxt*v*CF)
  dfCF
}