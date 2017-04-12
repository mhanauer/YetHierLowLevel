# YetHierLowLevel
# Jags-Ymet-XmetMultiHier-Mrobust.R 
# Accompanies the book:
#  Kruschke, J. K. (2015). Doing Bayesian Data Analysis, Second Edition: 
#  A Tutorial with R, JAGS, and Stan. Academic Press / Elsevier.

source("DBDA2E-utilities.R")

#===============================================================================

genMCMC = function( data , xName="x" , yName="y" , sName = "s", 
                    numSavedSteps=10000 , thinSteps=10 , saveName=NULL  ,
                    runjagsMethod=runjagsMethodDefault , 
                    nChains=nChainsDefault ) { 
  require(runjags)
  #-----------------------------------------------------------------------------
  # THE DATA.
  y = data[,yName]
  x = as.matrix(data[,xName],ncol=length(xName))
  s = as.numeric(factor(data[,sName]))
  # Do some checking that data make sense:
  if ( any( !is.finite(y) ) ) { stop("All y values must be finite.") }
  if ( any( !is.finite(x) ) ) { stop("All x values must be finite.") }
  cat("\nCORRELATION MATRIX OF PREDICTORS:\n ")
  show( round(cor(x),3) )
  cat("\n")
  flush.console()
  # Specify the data in a list, for later shipment to JAGS:
  dataList = list(
    x = x ,
    y = y ,
    s = s,
    Nsubj = max(s),  # should equal length(unique(s))
    Nx = dim(x)[2] ,
    Ntotal = dim(x)[1]
  )
  #-----------------------------------------------------------------------------
  # THE MODEL.
  modelString = "
  # Standardize the data:
  data {
    ym <- mean(y)
    ysd <- sd(y)
    for ( i in 1:Ntotal ) {
      zy[i] <- ( y[i] - ym ) / ysd
    }
    for ( j in 1:Nx ) {
      xm[j]  <- mean(x[,j])
      xsd[j] <-   sd(x[,j])
      for ( i in 1:Ntotal ) {
        zx[i,j] <- ( x[i,j] - xm[j] ) / xsd[j]
      }
    }
  }
  # Specify the model for standardized data:
  model {
    for ( i in 1:Ntotal ) {
      zy[i] ~ dt( zbeta0[s[i]] + sum( zbeta[s[1:Nx]] * zx[s[i,1:Nx]] ) , 1/zsigma^2 , nu )
    }
  for (i in 1:Nsubj){
    for (j in 1:Nx ){
      zbeta0[j] ~ dnorm( 0 , 1/2^2 ) # One intercept for each subject  
      zbeta[i,j] ~ dnorm( 0 , 1/2^2 ) # One for each covariate and subject
    }  
  }  
  # Priors vague on standardized scale. May need to fix this later :
    zbeta0Mu ~ dnorm( 0 , 1/2^2 ) 
    for ( j in 1:Nx ) {
      zbetaMu[j] ~ dnorm( 0 , 1/2^2 )
    }
    zsigma ~ dunif( 1.0E-5 , 1.0E+1 )
    nu ~ dexp(1/30.0)
    # Transform to original scale:
    nu ~ dexp(1/30.0)
  # Transform to original scale:
  for ( i in 1:Nsubj ) {
    for (j in 1:Nx) {
  beta[i,j] <- ( zbeta[i,j] / xsd[1:Nx] )*ysd
  beta0[j] <-  zbeta0[j]*ysd  + ym - sum( zbeta[i,j] * xm[1:Nx] / xsd[1:Nx] )*ysd
    }
  }


    betaMu[1:Nx] <- ( zbetaMu[1:Nx] / xsd[1:Nx] )*ysd
    beta0Mu <- zbeta0*ysd  + ym - sum( zbeta[1:Nx] * xm[1:Nx] / xsd[1:Nx] )*ysd
    sigma <- zsigma*ysd
  }
  " # close quote for modelString
  # Write out modelString to a text file
  writeLines( modelString , con="TEMPmodel.txt" )
  #-----------------------------------------------------------------------------
  # RUN THE CHAINS
  parameters = c(  "betoMu","betaMu","beta0" ,  "beta" ,  "sigma","zbeta0Mu","zbetaMu", 
                  "zbeta0" , "zbeta" , "zsigma", "nu" )
  adaptSteps = 500  # Number of steps to "tune" the samplers
  burnInSteps = 1000
  runJagsOut <- run.jags( method=runjagsMethod ,
                          model="TEMPmodel.txt" , 
                          monitor=parameters , 
                          data=dataList ,  
                          #inits=initsList , 
                          n.chains=nChains ,
                          adapt=adaptSteps ,
                          burnin=burnInSteps , 
                          sample=ceiling(numSavedSteps/nChains) ,
                          thin=thinSteps ,
                          summarise=FALSE ,
                          plots=FALSE )
  codaSamples = as.mcmc.list( runJagsOut )
  # resulting codaSamples object has these indices: 
  #   codaSamples[[ chainIdx ]][ stepIdx , paramIdx ]
  if ( !is.null(saveName) ) {
    save( codaSamples , file=paste(saveName,"Mcmc.Rdata",sep="") )
  }
  return( codaSamples )
} # end function

#===============================================================================

smryMCMC = function(  codaSamples , 
                      saveName=NULL ) {
  summaryInfo = NULL
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  paramName = colnames(mcmcMat)
  for ( pName in paramName ) {
    summaryInfo = rbind( summaryInfo , summarizePost( mcmcMat[,pName] ) )
  }
  rownames(summaryInfo) = paramName
  summaryInfo = rbind( summaryInfo , 
                       "log10(nu)" = summarizePost( log10(mcmcMat[,"nu"]) ) )
  if ( !is.null(saveName) ) {
    write.csv( summaryInfo , file=paste(saveName,"SummaryInfo.csv",sep="") )
  }
  return( summaryInfo )
}

#===============================================================================

