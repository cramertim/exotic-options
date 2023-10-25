#### functions for calculation---------------------------------------------------

## return matrix with only the first positive value of each row and rest 0
firstValueRow <-
  function(x) {
    cumSumMat <- matrix(NA, nrow = dim(x)[1], ncol = dim(x)[2])
    for(i in 1:(dim(x)[1])) {
      cumSumMat[i,] <- cumsum(x[i,])
    }
    
    cumSumMat2 <- cbind(matrix(0, nrow = dim(x)[1], ncol = 1), cumSumMat[, -(dim(x)[2])])
    ResultMat<-matrix(NA, nrow = dim(x)[1], ncol = dim(x)[2])
    for(i in 1:dim(x)[2]) {
      ResultMat[, i] <- ifelse(cumSumMat2[, i] > 0, 0, x[,i])
    }
    return(ResultMat)	
}

## heston simulation which creates stock paths
simulate.Heston <- function(S0, sqrt_V0, # sqrt is the initial volatility not variance!
                            kappa, theta, rho, epsilon, r, dr, N, M, dt){
  
  
  St <- S0
  V <- sqrt_V0^2 
  
  St_matrix <- matrix(NA, nrow = N, ncol = M)
  mu <- c(0,0)
  cov <- matrix(c(1, rho, rho, 1), nrow = 2, ncol = 2, byrow = FALSE)
  
  # calculate value of the volatility and price for N timesteps
  for(i in 1:N){
    
    # if M is equal to one the matrix has to be defined and transposed seperatly
    # otherwise the dimension and the data type are not correct
    if(M == 1){
      W <- t(matrix(mvrnorm(M, mu, cov)))
    } else {
      W <- mvrnorm(M, mu, cov)
    }
    
    # stochastic process for the variance
    V <- V + kappa*(theta-V)*dt + as.vector(epsilon*sqrt(V*dt) * W[,2])
    V[V < 0] <- 0
    
    # stochastic process for the stock price
    St = St * exp(((r-dr)-0.5 * V) * dt + as.vector(sqrt(V*dt) * W[,1]))
    
    St_matrix[i,] <- St
    
    # increment for the progressbar
    incProgress(1/N)
  }
  
  
  return(St_matrix)
  
}

## generate random stock paths with either Black-Scholes or Heston model
create.StockPaths <- function(simulationModel,
                              initialPrice,  maturity, 
                              volatility,    interestRate,
                              dividendRate = NULL,
                              
                              sampleSize,    numberTimeSteps = NULL,
                              skalarTimeSteps = NULL, #how many timesteps per day
                              startDate = NULL, # startDate for forwardStart
                              forwardStart = FALSE, # flag to indicate forward start
                              
                              # heston-model specific input paramters
                              kappa = NULL, # volatility mean-reversion speed
                              epsilon = NULL, # volatility of volatility
                              rho = NULL, # Correlation between stoch. vol and spot prices
                              theta = NULL  #long-term variance <- variance of the asset
){
  
  
  if(is.convertible.to.date(maturity)){
    
    currentDate <- Sys.Date()
    endDate <- as.Date(maturity)
    
    # if ... then forwardstart
    if(forwardStart) {
      
      startDate <- as.Date(startDate)
      
      # T2 - T1
      T <- as.numeric(difftime(as.Date(endDate), startDate, unit="weeks"))/52.25
      
    } else {
      # T2 - T0 
      T <- as.numeric(difftime(as.Date(endDate), as.Date(currentDate), unit="weeks"))/52.25
      
      # if numbertimesteps has to be an skalar of the days to maturity...
      if(!is.null(skalarTimeSteps) && is.null(numberTimeSteps)){
        
        numberTimeSteps <- as.numeric(difftime(as.Date(endDate), as.Date(currentDate), unit="days")) * skalarTimeSteps
        
      } else if(!is.null(skalarTimeSteps) && !is.null(numberTimeSteps)){
        
        stop("skalarTimeSteps and numberTimeSteps in function!")
        
      } else {
        
      }
      
    }
    
  } 
  else {
    # if maturity is numeric (time in years)
    T <- maturity
    
    # if numbertimesteps has to be an skalar of the days to maturity...
    if(!is.null(skalarTimeSteps) && is.null(numberTimeSteps)){
      
      currentDate <- Sys.Date()
      maturityDate <- as.Date(Sys.Date() + floor(365 * T))
      
      numberTimeSteps <- as.numeric(difftime(as.Date(maturityDate), as.Date(currentDate), unit="days")) * skalarTimeSteps
      
    } else if(!is.null(skalarTimeSteps) && !is.null(numberTimeSteps)){
      
      stop("skalarTimeSteps and numberTimeSteps in function!")
      
    } else {
      
    }
    
  }
  
  # define general parameters
  S0 <- initialPrice        # initial stockprice
  sigma <- volatility       # in Percent
  r <- interestRate         # risk free interest rate
  M <- sampleSize 
  N <- numberTimeSteps      
  dt <- T/N                 # length of subintervall
  
  if(is.null(dividendRate)){
    dr <- 0
  } else{
    dr <- dividendRate
  }
  
  discF <- exp(-r*T)        # discount factor
  
  
  
  
  switch (simulationModel,
          "Black-Scholes"={
            
            # Monte Carlo via Back-Scholes-model
            epsilon_black <- matrix(rnorm(N * M), N, M)
            
            term <- matrix(NA_real_, nrow=N, ncol=M)
            
            term <-  exp(dt * ((r-dr) - ((sigma)^2)/2) + sqrt(dt) * epsilon_black * sigma)
            term2 <- exp(dt * ((r-dr) - ((sigma)^2)/2) - sqrt(dt) * epsilon_black * sigma)
            
            St <- S0 * colCumprods(term)
            
            St2 <- S0 * colCumprods(term2)
            
          },
          "Heston"={
            
            # Monte Carlo via the Heston-model
            St <- simulate.Heston(S0 = S0, sqrt_V0 = sigma,
                                  r = r, dr = dr,
                                  N = N, dt = dt,
                                  M = M,
                                  kappa = kappa, theta = theta, 
                                  rho = rho, epsilon = epsilon)
            St2 <- St
            
          }
  )
  
  
  
  returnlist <- list(St = St, St2 = St2, maturity = T, 
                     interestRate = r, sampleSize = M, 
                     numberTimeSteps = N, skalarTimeSteps = skalarTimeSteps,
                     startDate = startDate,
                     stockPrice = S0, dividendRate = dividendRate,
                     volatility = sigma,
                     
                     kappa = kappa,
                     epsilon = epsilon,
                     rho = rho,
                     theta = theta,
                     
                     simulationModel = simulationModel
  )
}


## checks if handed over string is convertible to an date object
is.convertible.to.date <- function(x) !is.na(as.Date(as.character(x), format = '%Y-%m-%d'))

## function to exclude bank holidays from calendar
create.Holidays <- function(){
  year <- getRmetricsOption("currentYear")
  years <- numeric()
  for(i in 0:15){
    years[i+1] <- year + i
  }
  hol <- holidayNYSE(years, type = "standard")
  holidays <- as.Date(hol@Data)
  return(holidays)
}

# function for option pricing
price.Option <- function(type,              # "american", "asian", "european" etc.
                         optionType = NULL, # "Call", "Put"
                         subType = NULL,    # not for all types, e.g. "downAndOut"
                         
                         stockPaths,        # stockpath object containing all needed variables 
                         strikePrice = NULL,       # strikeprice of the option
                         
                         # asian
                         asianRange = NULL, #"Total" - every timestep, "interval" -> intervalObs, "period" -> periodObs
                         intervalObs_asian = NULL, # x days between each observation
                         periodObs_asian = NULL, # last x days of timeseries observed
                         
                         # barrier option
                         barrier = NULL, # barrier as number
                         
                         # binary option cashOrNothing
                         binaryPayoff = NULL,
                         
                         # chooser
                         chooser = FALSE, # flag to indicate chooser
                         choosingDate = NULL,
                         chooserOption = NULL, # option where put/call will be choosed
                         chooserSubType = NULL, # subtype of chooser option
                         
                         # compound
                         secondSampleSize = NULL, # second... -> refers to parameters of 2. option
                         secondNumberTimeSteps = NULL,
                         secondSkalarTimeSteps = NULL,
                         secondExerciseDate = NULL,
                         type2 = NULL,
                         secondStrikePrice  = NULL,
                         secondSubType = NULL,
                         
                         # european
                         # -
                         
                         # forward start
                         startDate = NULL,
                         
                         # gap
                         strikePrice1 = NULL, # strike price of the option
                         strikePrice2 = NULL # trigger price of the option
                         
                         # lookback
                         # -
){
  
  
  # define all general variables
  T <- stockPaths$maturity

  r <- stockPaths$interestRate
  M <- stockPaths$sampleSize
  N <- NULL
  sigma <- stockPaths$volatility
  
  
  
  
  ## chooser specific
  # calculate choosing date as number in dependency of the skalatTimesteps
  if(type == "chooser" && is.convertible.to.date(choosingDate)){
    
    currentDate <- Sys.Date()
    choosingDate <- as.Date(choosingDate)
    
    cD <- as.numeric(difftime(as.Date(choosingDate), as.Date(currentDate), unit="days")) * stockPaths$skalarTimeSteps
  } else if(type != "chooser"){
    
    # this is for the call/put of the options in the chooser option, cD has been already calculated from date to numeric
    
  } else {
    stop("Choosing date not convertible to date.")
  }
  
  
  # if chooser TRUE an option with maturity choosingDate in the chooser option is calculated
  if(chooser == TRUE){
    N <- choosingDate
  }else if(chooser == FALSE){
    N <- stockPaths$numberTimeSteps
  }
  
  
  dr <- stockPaths$dividendRate 
  
  K <- strikePrice
  
  
  # discount factor
  discF <- exp(-r*T)
  
  St <- stockPaths$St
  St2 <- stockPaths$St2
  
  St <- head(St, N)
  St2 <- head(St2, N)
  
  endPrice <- St[N ,]
  endPrice2 <- St2[N ,]
  
  
  ## define option specific variables
  # american
  # -
  
  # asian
  a_range <- asianRange
  a_period <- NULL #corresponding to periodObs_asian and used in switch loop
  a_interval <- NULL #corresponding to intervalObs_asian and used in switch loop
  periodObs_asian <- periodObs_asian # for chooser
  intervalObs_asian <- intervalObs_asian
  
  
  # barrier option
  B <- barrier 
  
  # binary option cashOrNothing
  BP <- binaryPayoff
  
  # european
  # -
  
  # chooser
  cO <- chooserOption
  cST <- chooserSubType
  
  
  # compound
  M2 <- secondSampleSize
  N2 <- secondNumberTimeSteps
  N2_skalar <- secondSkalarTimeSteps
  t2 <- type2
  K2.com <- secondStrikePrice 
  T2 <- secondExerciseDate 
  subType2 = secondSubType
  
  
  # european
  # -
  
  # forward start
  startDate <- as.Date(startDate)
  
  # gap
  K1 <- strikePrice1
  K2 <- strikePrice2
  
  
  # lookback
  # -

  
  
  
  switch(type,
         "american"={
           
           if(optionType == "Call"){
             
             #transpose St matrix in order to fit with rest of code with override of variables
             m <- N
             n <- M
             St <- t(St)
             
             #t est: if each entry bigger than Strike price 
             X<-ifelse(St>strikePrice,St,NA)
             
             CFL<-matrix(pmax(0,St - strikePrice), nrow=n, ncol=m)
             # param for regressions
             Xsh<-X[,-m]
             X2sh<-Xsh*Xsh
             # discount
             Y1<-CFL*exp(-1*r*(T/m))
             Y2<-cbind((matrix(NA, nrow=n, ncol=m-1)), Y1[,m])
             continuationVal<-matrix(NA, nrow=n, ncol=m-1)
             
             # regression and  iterative process
             try(for(i in (m-1):1) {
               reg1<-lm(Y2[,i+1]~Xsh[,i]+X2sh[,i])
               continuationVal[,i]<-(matrix(reg1$coefficients)[1,1])+((matrix(reg1$coefficients)[2,1])*Xsh[,i])+((matrix(reg1$coefficients)[3,1])*X2sh[,i])
               continuationVal[,i]<-(ifelse(is.na(continuationVal[,i]),0,continuationVal[,i]))
               Y2[,i]<-ifelse(CFL[,i]>continuationVal[,i], Y1[,i], Y2[,i+1]*exp(-1*r*(T/m)))
               incProgress(1 /(1.1*m))
             }
             , silent = TRUE)
             
             continuationVal<-ifelse(is.na(continuationVal),0,continuationVal)
             CVp<-cbind(continuationVal, (matrix(0, nrow=n, ncol=1)))
             POF<-ifelse(CVp>CFL,0,CFL)
             FPOF<-firstValueRow(POF)
             dFPOF<-matrix(NA, nrow=n, ncol=m)
             for(i in 1:m) {
               dFPOF[,i]<-FPOF[,i]*exp(-1*T/m*r*i)
               incProgress(m/10)
             }
             payoff <- rowSums(dFPOF)
             
             
           } else if(optionType == "Put"){
             
             #transpose St matrix in order to fit with rest of code with override of variables
             m <- N
             n <- M
             St <- t(St)
             
             #test: if each entry smaller than Strike price 
             X<-ifelse(St<strikePrice,St,NA)
             
             CFL<-matrix(pmax(0,strikePrice-St), nrow=n, ncol=m)
             # param for regression
             Xsh<-X[,-m]
             X2sh<-Xsh*Xsh
             
             # discount
             Y1<-CFL*exp(-1*r*(T/m))
             Y2<-cbind((matrix(NA, nrow=n, ncol=m-1)), Y1[,m])
             continuationVal<-matrix(NA, nrow=n, ncol=m-1)
             
             # regression and iterative process
             try(for(i in (m-1):1) {
               reg1<-lm(Y2[,i+1]~Xsh[,i]+X2sh[,i])
               continuationVal[,i]<-(matrix(reg1$coefficients)[1,1])+((matrix(reg1$coefficients)[2,1])*Xsh[,i])+((matrix(reg1$coefficients)[3,1])*X2sh[,i])
               continuationVal[,i]<-(ifelse(is.na(continuationVal[,i]),0,continuationVal[,i]))
               Y2[,i]<-ifelse(CFL[,i]>continuationVal[,i], Y1[,i], Y2[,i+1]*exp(-1*r*(T/m)))
               incProgress(1 /(1.1*m))
             }
             , silent = TRUE)
             
             continuationVal<-ifelse(is.na(continuationVal),0,continuationVal)
             CVp<-cbind(continuationVal, (matrix(0, nrow=n, ncol=1)))
             POF<-ifelse(CVp>CFL,0,CFL)
             FPOF<-firstValueRow(POF)
             dFPOF<-matrix(NA, nrow=n, ncol=m)
             for(i in 1:m) {
               dFPOF[,i]<-FPOF[,i]*exp(-1*T/m*r*i)
               incProgress(m/10)
             }
             payoff <- rowSums(dFPOF)
           }
         },
         "asian"={
           if(a_range == "Interval"){
             #get indices for the closing price of each day defined by the interval
             obs <- seq(1 ,N ,by = intervalObs_asian*stockPaths$skalarTimeSteps) + stockPaths$skalarTimeSteps - 1
             
             St <- matrix(St[obs,], ncol = M)
             St2 <- matrix(St2[obs,], ncol = M)
             
             meanSt <- colMeans2(St)
             meanSt2 <- colMeans2(St2)
             
             geometricMean <- apply(St, 2, function(x) exp(mean(log(x[x>0]))))  #geometric Mean
             geometricMean2 <- apply(St2, 2, function(x) exp(mean(log(x[x>0]))))  #geometric Mean
             
           }
           else if(a_range == "Period"){   
             # only the last x days are considered
             a_period <- periodObs_asian * stockPaths$skalarTimeSteps
             
             St <- tail(St, a_period)
             St2 <-tail(St2, a_period)
             
             meanSt <- colMeans2(St)
             meanSt2 <- colMeans2(St2)
             
             geometricMean <- apply(St, 2, function(x) exp(mean(log(x[x>0]))))  
             geometricMean2 <- apply(St2, 2, function(x) exp(mean(log(x[x>0])))) 
           }
           else{
             # default - all datapoints are considered
             meanSt <- colMeans2(St)
             meanSt2 <- colMeans2(St2)
             
             geometricMean <- apply(St, 2, function(x) exp(mean(log(x[x>0]))))  
             geometricMean2 <- apply(St2, 2, function(x) exp(mean(log(x[x>0]))))
           }
           
           switch(subType,
                  "arithmetic"={
                    if (optionType == "Call") {
                      payoff <- ifelse(meanSt > K, meanSt - K, 0)
                      payoff2 <- ifelse(meanSt2 > K, meanSt2 - K, 0)
                    }
                    else if(optionType == "Put") {
                      payoff <- ifelse(K > meanSt, K - meanSt, 0) 
                      payoff2 <- ifelse(K > meanSt2, K - meanSt2, 0) 
                    }
                  },
                  "geometric"={
                    if (optionType == "Call") {
                      payoff <- ifelse(geometricMean > K, geometricMean - K, 0)
                      payoff2 <-  ifelse(geometricMean2 > K, geometricMean2 - K, 0)
                    }
                    else if(optionType == "Put") {
                      payoff <- ifelse(K > geometricMean, K - geometricMean, 0)
                      payoff2 <- ifelse(K > geometricMean2, K - geometricMean2, 0) 
                    }
                  })
         },
         "barrier"={
           maxPrice <- colMaxs(St)
           maxPrice2 <- colMaxs(St2)
           minPrice <- colMins(St)
           minPrice2 <- colMins(St2)
           
           
           switch (subType,
                   "Down-And-Out"={
                     if (optionType == "Call") {
                       # price stays over the barrier -> 1
                       mask <- matrix(ifelse(minPrice > B, 1, 0), ncol = M)
                       mask2 <- matrix(ifelse(minPrice2 > B, 1, 0), ncol = M)
                       
                       payoff <- (matrix(ifelse(St[N, ] > K, St[N, ] - K, 0), ncol = M) * mask)
                       payoff2 <- (matrix(ifelse(St2[N, ] > K, St2[N, ] - K, 0), ncol = M) * mask2)
                       
                     }
                     else if(optionType == "Put") {
                       mask <- matrix(ifelse(minPrice > B, 1, 0), ncol = M)
                       mask2 <- matrix(ifelse(minPrice2 > B, 1, 0), ncol = M)
                       
                       payoff <- (matrix(ifelse(K > St[N, ], K - St[N, ], 0), ncol = M) * mask)
                       payoff2 <- (matrix(ifelse(K > St2[N, ], K - St2[N, ], 0), ncol = M) * mask2)
                     }
                   },
                   "Down-And-In"={
                     if (optionType == "Call") {
                       mask <- matrix(ifelse(B > minPrice, 1, 0), ncol = M)
                       mask2 <- matrix(ifelse(B > minPrice2, 1, 0), ncol = M)
                       
                       payoff <- (matrix(ifelse(St[N, ] > K, St[N, ] - K, 0), ncol = M) * mask)
                       payoff2 <- (matrix(ifelse(St2[N, ] > K, St2[N, ] - K, 0), ncol = M) * mask2) 
                       
                     }
                     else if(optionType == "Put") {
                       mask <- matrix(ifelse(B > minPrice, 1, 0), ncol = M)
                       mask2 <- matrix(ifelse(B > minPrice2, 1, 0), ncol = M)
                       
                       payoff <- (matrix(ifelse(K > St[N, ], K - St[N, ], 0), ncol = M) * mask)
                       payoff2 <-  (matrix(ifelse(K > St2[N, ], K - St2[N, ], 0), ncol = M) * mask2)
                     }
                   },
                   "Up-And-Out"={
                     if (optionType == "Call") {
                       mask <- matrix(ifelse(B > maxPrice, 1, 0), ncol = M)
                       mask2 <- matrix(ifelse(B > maxPrice2, 1, 0), ncol = M)
                       
                       payoff <- (matrix(ifelse(St[N, ] > K, St[N, ] - K, 0), ncol = M) * mask)
                       payoff2 <-  (matrix(ifelse(St2[N, ] > K, St2[N, ] - K, 0), ncol = M) * mask2)
                       
                     }
                     else if(optionType == "Put") {
                       mask <- matrix(ifelse(B > maxPrice, 1, 0), ncol = M)
                       mask2 <- matrix(ifelse(B > maxPrice2, 1, 0), ncol = M)
                       
                       payoff <- (matrix(ifelse(K > St[N, ], K - St[N, ], 0), ncol = M) * mask)
                       payoff2 <- (matrix(ifelse(K > St2[N, ], K - St2[N, ], 0), ncol = M) * mask2)
                     }
                   },
                   "Up-And-In"={
                     if (optionType == "Call") {
                       mask <- matrix(ifelse(maxPrice > B, 1, 0), ncol = M)
                       mask2 <- matrix(ifelse(maxPrice2 > B, 1, 0), ncol = M)
                       
                       payoff <- (matrix(ifelse(St[N, ] > K, St[N, ] - K, 0), ncol = M) * mask)
                       payoff2 <- (matrix(ifelse(St2[N, ] > K, St2[N, ] - K, 0), ncol = M) * mask2)
                       
                     }
                     else if(optionType == "Put") {
                       mask <- matrix(ifelse(maxPrice > B, 1, 0), ncol = M)
                       mask2 <- matrix(ifelse(maxPrice2 > B, 1, 0), ncol = M)
                       
                       payoff <- (matrix(ifelse(K > St[N, ], K - St[N, ], 0), ncol = M) * mask)
                       payoff2 <- (matrix(ifelse(K > St2[N, ], K - St2[N, ], 0), ncol = M) * mask2)
                     }
                   }
           )
           
         },
         "binary"={
           switch(subType,
                  "Cash-or-nothing"={
                    if (optionType == "Call") {
                      payoff <- ifelse(St[N, ] > K, BP, 0)
                      payoff2 <- ifelse(St2[N, ] > K, BP, 0)
                    }
                    else if(optionType == "Put") {
                      payoff <- ifelse(K > St[N, ], BP, 0) 
                      payoff2 <- ifelse(K > St2[N, ], BP, 0)
                    }
                  },
                  "Asset-or-nothing"={
                    if (optionType == "Call") {
                      payoff <- ifelse(St[N, ] > K, St[N, ], 0)
                      payoff2 <- + ifelse(St2[N, ] > K, St2[N, ], 0)
                    }
                    else if(optionType == "Put") {
                      payoff <- ifelse(K > St[N, ], St[N, ], 0)
                      payoff2 <- ifelse(K > St2[N, ], St2[N, ], 0) 
                    }
                  }
           )
         },
         "chooser"={
           # See if Call or Put on the choosingDate is worth more and decide for it
           # Then calculate the price of the option with this decision at T
           
           ## Payoff at ChoosingDate - call
           call <- price.Option(type = cO,
                                optionType= "Call",
                                subType = cST,
                                
                                stockPaths = c(stockPaths),        # stockpath object containing all needed variables
                                strikePrice = K,       # strikeprice of the option
                                
                                # asian
                                asianRange = a_range,
                                intervalObs_asian = intervalObs_asian, 
                                periodObs_asian = (periodObs_asian*stockPaths$skalarTimeSteps - (N-cD))/stockPaths$skalarTimeSteps, # hier änderung
                                
                                # barrier option
                                barrier = B,
                                
                                # binary option cashOrNothing
                                binaryPayoff = BP,
                                
                                # chooser -> Flag to check determined day
                                chooser = TRUE,
                                choosingDate = cD, 
                                
                                # forward start
                                startDate = startDate,
                                
                                # gap
                                strikePrice1 = K1,
                                strikePrice2 = K2
                                
           )
           
           ## put
           put <- price.Option(type = cO,
                               optionType= "Put",
                               subType = cST,
                               
                               stockPaths = c(stockPaths),        # stockpath object containing all needed variables
                               strikePrice = K,       # strikeprice of the option
                               
                               # asian
                               asianRange = a_range,
                               intervalObs_asian = intervalObs_asian,
                               periodObs_asian = (periodObs_asian*stockPaths$skalarTimeSteps - (N-cD))/stockPaths$skalarTimeSteps,
                               
                               # barrier option
                               barrier = B,
                               
                               # binary option cashOrNothing
                               binaryPayoff = BP,
                               
                               # chooser -> Flag to check determined day
                               chooser = TRUE,
                               choosingDate = cD, 
                               
                               # forward start
                               startDate = startDate,
                               
                               # gap
                               strikePrice1 = K1,
                               strikePrice2 = K2
           )
           
           #Choose the higher Payoff Option at ChoosingDate
           maskCall <- matrix(ifelse(call$payoff > put$payoff, 1, 0), ncol = M)
           maskPut  <- matrix(ifelse(put$payoff > call$payoff, 1, 0), ncol = M)

           #Sometimes both Call and Put reward the same Payoff (Zero), in this Case choose random optionType
           #By adding maskCall to maskPut the Zeros will identify this cases
           Mask <- matrix(maskCall + maskPut, ncol = M)
           random <- matrix(runif(M, min = -1, max = 1), ncol = M)
           
           #Check for Zero and then add the decision for Call/ Put
           for(i in 1:M) {
             if(Mask[1, i] < 1){
               if(random[i] > 0) {
                 maskCall[1, i] <- 1
               }else {
                 maskPut[1, i] <- 1
               } 
             }
           }
           
           ## Payoff at the end of the timeframe - call
           call <- price.Option(type = cO,
                                optionType= "Call",
                                subType = cST,
                                
                                stockPaths = c(stockPaths),        # stockpath object containing all needed variables
                                strikePrice = K,       # strikeprice of the option
                                
                                # asian
                                asianRange = a_range,
                                intervalObs_asian = intervalObs_asian,
                                periodObs_asian = periodObs_asian,
                                
                                # barrier option
                                barrier = B,
                                
                                # binary option cashOrNothing
                                binaryPayoff = BP,
                                
                                # forward start
                                startDate = startDate,
                                
                                # gap
                                strikePrice1 = K1,
                                strikePrice2 = K2
           )
           
           ## put
           put <- price.Option(type = cO,
                               optionType= "Put",
                               subType = cST,
                               
                               stockPaths = c(stockPaths),        # stockpath object containing all needed variables
                               strikePrice = K,       # strikeprice of the option
                               
                               # asian
                               asianRange = a_range,
                               intervalObs_asian = intervalObs_asian,
                               periodObs_asian = periodObs_asian,
                               
                               # barrier option
                               barrier = B,
                               
                               # binary option cashOrNothing
                               binaryPayoff = BP,
                               
                               # forward start
                               startDate = startDate,
                               
                               # gap
                               strikePrice1 = K1,
                               strikePrice2 = K2
                               
           )
           
           # Multiply with mask to eliminate the Paths that are not chosen
           call.payoff <- call$payoff_plot * maskCall  #*2 because otherwise /4 due to antithetic payoff
           call.payoff2 <- call$payoff_anti * maskCall #*2 because otherwise /4 due to antithetic payoff
           
           put.payoff <- put$payoff_plot * maskPut  
           put.payoff2 <- put$payoff_anti * maskPut 
           
           payoff <- call.payoff + put.payoff
           payoff2 <- call.payoff2 + put.payoff2
           
           
         },
         "compound"={
           switch(subType,
                  #1. Option is always an european option
                  #2. Option can be anything except compound
                  
                  "CoC"={
                    # Call on Call
                    
                    #calculate the payoff of the first option
                    first <- price.Option(type = "european",
                                          optionType = "Call",
                                          strikePrice = K,
                                          stockPaths = c(stockPaths))
                    payoff <- first$payoff
                    payoff2 <- first$payoff_anti
                    
                    if(is.convertible.to.date(secondExerciseDate)){
                      
                      T2 <- as.numeric(difftime(as.Date(secondExerciseDate), as.Date(Sys.Date() + floor(first[["stockPaths"]][["maturity"]] * 365)), unit="weeks"))/52.25
                      
                    }
                    
                    # if the value of the option is greater than the first strike price, execute the option
                    for(i in 1:M) {
                      incProgress(1/M)
                      
                      if(type2 == "asian" || type2 == "lookback"){
                        
                        #for every executed path calculate the value of the second option
                        stockPathsSecond <- create.StockPaths(simulationModel = "Black-Scholes",
                                                              initialPrice = endPrice[i],       maturity = T2, #T2 - T1
                                                              volatility = sigma,               interestRate = r, 
                                                              dividendRate = dr,
                                                              sampleSize = M2,                  skalarTimeSteps = N2_skalar)
                        
                      } else {
                        
                        #for every executed path calculate the value of the second option
                        stockPathsSecond <- create.StockPaths(simulationModel = "Black-Scholes",
                                                              initialPrice = endPrice[i],       maturity = T2, #T2 - T1
                                                              volatility = sigma,               interestRate = r, 
                                                              dividendRate = dr,
                                                              sampleSize = M2,                  numberTimeSteps = N2)
                        
                      }
                      
                      ## payoff second
                      second <- price.Option(type = t2,
                                             optionType= "Call",
                                             subType = subType2,
                                             strikePrice = K2.com,
                                             stockPaths = c(stockPathsSecond),
                                             
                                             # asian
                                             asianRange = a_range,
                                             intervalObs_asian = intervalObs_asian,
                                             periodObs_asian = periodObs_asian,
                                             
                                             # barrier option
                                             barrier = B,
                                             
                                             # binary option cashOrNothing
                                             binaryPayoff = BP,
                                             
                                             # chooser 
                                             chooser = FALSE,
                                             choosingDate = NULL,
                                             chooserOption = NULL,
                                             chooserSubType = NULL,
                                             
                                             # forward start
                                             startDate = startDate,
                                             
                                             # gap
                                             strikePrice1 = K1,
                                             strikePrice2 = K2
                      )
                      
                      #the payoff at T1 is the value of the second option minus the price to execute it
                      payoff[i] <- max(second$optionValue - K, 0)
                      payoff2[i] <- max(second$optionValue - K, 0)
                      
                    }
                    
                  },
                  "CoP"={
                    #calculate the payoff of the first option
                    first <- price.Option(type = "european",
                                          optionType = "Call",
                                          strikePrice = K,
                                          stockPaths = c(stockPaths))
                    payoff <- first$payoff
                    payoff2 <- first$payoff_anti
                    
                    if(is.convertible.to.date(secondExerciseDate)){
                      
                      T2 <- as.numeric(difftime(as.Date(secondExerciseDate), as.Date(Sys.Date() + floor(first[["stockPaths"]][["maturity"]] * 365)), unit="weeks"))/52.25
                      
                    }
                    
                    # if the value of the option is greater than the first strike price, execute the option
                    for(i in 1:M) {
                     
                      if(type2 == "asian" || type2 == "lookback"){
                        
                        #for every executed path calculate the value of the second option
                        stockPathsSecond <- create.StockPaths(simulationModel = "Black-Scholes",
                                                              initialPrice = endPrice[i],       maturity = T2, #T2 - T1
                                                              volatility = sigma,               interestRate = r, 
                                                              dividendRate = dr,
                                                              sampleSize = M2,                  skalarTimeSteps = N2_skalar)
                        
                      } else {
                        
                        #for every executed path calculate the value of the second option
                        stockPathsSecond <- create.StockPaths(simulationModel = "Black-Scholes",
                                                              initialPrice = endPrice[i],       maturity = T2, #T2 - T1
                                                              volatility = sigma,               interestRate = r, 
                                                              dividendRate = dr,
                                                              sampleSize = M2,                  numberTimeSteps = N2)
                        
                      }
                      
                      ## payoff second
                      second <- price.Option(type = t2,
                                             optionType= "Put",
                                             subType = subType2,
                                             strikePrice = K2.com,
                                             stockPaths = c(stockPathsSecond),
                                             
                                             # asian
                                             asianRange = a_range,
                                             intervalObs_asian = intervalObs_asian,
                                             periodObs_asian = periodObs_asian,
                                             
                                             # barrier option
                                             barrier = B,
                                             
                                             # binary option cashOrNothing
                                             binaryPayoff = BP,
                                             
                                             # chooser 
                                             chooser = FALSE,
                                             choosingDate = NULL,
                                             chooserOption = NULL,
                                             chooserSubType = NULL,
                                             
                                             # forward start
                                             startDate = startDate,
                                             
                                             # gap
                                             strikePrice1 = K1,
                                             strikePrice2 = K2
                                            
                      )
                      
                      #the payoff at T1 is the value of the second option minus the price to execute it
                      payoff[i] <- max(second$optionValue - K, 0)
                      payoff2[i] <- max(second$optionValue - K, 0)
                    }
                    
                  },
                  "PoC"={
                    #calculate the payoff of the first option
                    first <- price.Option(type = "european",
                                          optionType = "Put",
                                          strikePrice = K,
                                          stockPaths = c(stockPaths))
                    payoff <- first$payoff
                    payoff2 <- first$payoff_anti
                    
                    if(is.convertible.to.date(secondExerciseDate)){
                      
                      T2 <- as.numeric(difftime(as.Date(secondExerciseDate), as.Date(Sys.Date() + floor(first[["stockPaths"]][["maturity"]] * 365)), unit="weeks"))/52.25
                      
                    }
                    # if the value of the option is greater than the first strike price, execute the option
                    for(i in 1:M) {
                      if(type2 == "asian" || type2 == "lookback"){
                        
                        #for everey executed paath calculate the value of the second option
                        stockPathsSecond <- create.StockPaths(simulationModel = "Black-Scholes",
                                                              initialPrice = endPrice[i],       maturity = T2, #T2 - T1
                                                              volatility = sigma,               interestRate = r, 
                                                              dividendRate = dr,
                                                              sampleSize = M2,                  skalarTimeSteps = N2_skalar)
                        
                      } else {
                        
                        #for everey executed paath calculate the value of the second option
                        stockPathsSecond <- create.StockPaths(simulationModel = "Black-Scholes",
                                                              initialPrice = endPrice[i],       maturity = T2, #T2 - T1
                                                              volatility = sigma,               interestRate = r, 
                                                              dividendRate = dr,
                                                              sampleSize = M2,                  numberTimeSteps = N2)
                        
                      }
                      
                      ## payoff second
                      second <- price.Option(type = t2,
                                             optionType= "Call",
                                             subType = subType2,
                                             strikePrice = K2.com,
                                             stockPaths = c(stockPathsSecond),
                                             
                                             # asian
                                             asianRange = a_range,
                                             intervalObs_asian = intervalObs_asian,
                                             periodObs_asian = periodObs_asian,
                                             
                                             # barrier option
                                             barrier = B,
                                             
                                             # binary option cashOrNothing
                                             binaryPayoff = BP,
                                             
                                             # chooser 
                                             chooser = FALSE,
                                             choosingDate = NULL,
                                             chooserOption = NULL,
                                             chooserSubType = NULL,
                                          
                                             # forward start
                                             startDate = startDate,
                                             
                                             # gap
                                             strikePrice1 = K1,
                                             strikePrice2 = K2
                      )
                      
                      #the payoff at T1 is the value of the second option minus the price to execute it
                      payoff[i] <- max(K - second$optionValue, 0) 
                      payoff2[i] <- max(K - second$optionValue, 0) 
                      
                    }
                    
                  },
                  "PoP"={
                    #calculate the payoff of the first option
                    first <- price.Option(type = "european",
                                          optionType = "Put",
                                          strikePrice = K,
                                          stockPaths = c(stockPaths))
                    payoff <- first$payoff
                    payoff2 <- first$payoff_anti
                    
                    if(is.convertible.to.date(secondExerciseDate)){
                      
                      T2 <- as.numeric(difftime(as.Date(secondExerciseDate), as.Date(Sys.Date() + floor(first[["stockPaths"]][["maturity"]] * 365)), unit="weeks"))/52.25
                      
                    }
                    # if the value of the option is greater than the first strike price, execute the option
                    for(i in 1:M) {
                      
                      #for everey executed paath calculate the value of the second option
                      if(type2 == "asian" || type2 == "lookback"){
                        
                        #for every executed path calculate the value of the second option
                        stockPathsSecond <- create.StockPaths(simulationModel = "Black-Scholes",
                                                              initialPrice = endPrice[i],       maturity = T2, #T2 - T1
                                                              volatility = sigma,               interestRate = r, 
                                                              dividendRate = dr,
                                                              sampleSize = M2,                  skalarTimeSteps = N2_skalar)
                        
                      } else {
                        
                        #for every executed path calculate the value of the second option
                        stockPathsSecond <- create.StockPaths(simulationModel = "Black-Scholes",
                                                              initialPrice = endPrice[i],       maturity = T2, #T2 - T1
                                                              volatility = sigma,               interestRate = r, 
                                                              dividendRate = dr,
                                                              sampleSize = M2,                  numberTimeSteps = N2)
                        
                      }
                      
                      ## payoff second
                      second <- price.Option(type = t2,
                                             optionType= "Put",
                                             subType = subType2,
                                             strikePrice = K2.com,
                                             stockPaths = c(stockPathsSecond),
                                             
                                             # asian
                                             asianRange = a_range,
                                             intervalObs_asian = intervalObs_asian,
                                             periodObs_asian = periodObs_asian,
                                             
                                             # barrier option
                                             barrier = B,
                                             
                                             # binary option cashOrNothing
                                             binaryPayoff = BP,
                                             
                                             # chooser
                                             chooser = FALSE,
                                             choosingDate = NULL,
                                             chooserOption = NULL,
                                             chooserSubType = NULL,
                                            
                                             # forward start
                                             startDate = startDate,
                                             
                                             # gap
                                             strikePrice1 = K1,
                                             strikePrice2 = K2
                      )
                      
                      #the payoff at T1 is the value of the second option minus the price to execute it
                      payoff[i] <- max(K - second$optionValue, 0)
                      payoff2[i] <- max(K - second$optionValue, 0)
                    }
                    
                  }
           )
         },
         "european"={
           if (optionType == "Call") {
             payoff <- ifelse(St[N, ] > K, St[N, ] - K, 0)
             payoff2 <- ifelse(St2[N, ] > K, St2[N, ] - K, 0)
           }
           else if(optionType == "Put") {
             payoff <- ifelse(K > St[N, ], K - St[N, ], 0)
             payoff2 <- ifelse(K > St2[N, ], K - St2[N, ], 0)
           }
         },
         "forwardStart"={
           if(!is.null(startDate)){
             
             # define T1 (in years)
             T1 <- as.numeric(difftime(startDate, Sys.Date(), unit="weeks"))/52.25
             discF <- exp(-r*T1)
             
           } else {
             stop("no startdate defined")
           }
           
           
           # calculating payoff
           if (optionType == "Call") {
             
             c <- ifelse(St[N, ] > K, St[N, ] - K, 0)  #old european payoff
             c2 <- ifelse(St2[N, ] > K, St2[N, ] - K, 0)
             # exp((r-q)T1) with dividend q
             payoff <- c * exp(r*T1)
             payoff2 <- c2 * exp(r*T1)
             
             
           }
           else if(optionType == "Put") {
             
             c <- ifelse(K > St[N, ], K - St[N, ], 0)  #old european payoff
             c2 <- ifelse(K > St2[N, ], K - St2[N, ], 0)
             payoff <- c * exp(r*T1)
             payoff2 <- c2 * exp(r*T1)
             
           }
         },
         "gap"={
           
           if (optionType == "Call") {
             payoff <- ifelse(St[N, ] > K2, St[N, ] - K1, 0)
             payoff2 <- ifelse(St2[N, ] > K2, St2[N, ] - K1, 0)
           }
           else if(optionType == "Put") {
             payoff <- ifelse(K2 > St[N, ], K1 - St[N, ], 0)
             payoff2 <- ifelse(K2 > St2[N, ], K1 - St2[N, ], 0)
           }
           
         },
         "lookback"={
           maxPrice <- colMaxs(St)
           maxPrice2 <- colMaxs(St2)
           minPrice <- colMins(St)
           minPrice2 <- colMins(St2)
           
           switch(subType,
                  "Fixed"={
                    if (optionType == "Call") {
                      payoff <- ifelse(maxPrice > K, maxPrice - K, 0)
                      payoff2 <- ifelse(maxPrice2 > K, maxPrice2 - K, 0)
                    }
                    else if(optionType == "Put") {
                      payoff <- ifelse(K > minPrice, K - minPrice, 0)
                      payoff2 <- ifelse(K > minPrice2, K - minPrice2, 0)
                    }
                  },
                  "Floating"={
                    if (optionType == "Call") {
                      payoff <- ifelse(endPrice > minPrice, endPrice - minPrice, 0)
                      payoff2 <- ifelse(endPrice2 > minPrice2, endPrice2 - minPrice2, 0)
                    }
                    else if(optionType == "Put") {
                      payoff <- ifelse(endPrice < maxPrice, maxPrice - endPrice, 0)
                      payoff2 <- ifelse(endPrice2 < maxPrice2, maxPrice2 - endPrice2, 0)
                    }
                  }       
           )
           
         }
  )
  
  #discount the payoff
  if(type != "american"){
    payoff_plot <- payoff
    payoff_anti <- payoff2
    
    #payoff <- payoff/2                # antithetic
    payoff <- (payoff + payoff2)/2     # alternative plot
    
    C0 <- discF * payoff
    C0_plot <- discF * payoff_plot
    
  }else if(type == "american"){
    
    #already discounted in switch function, therefore overwritten with payoff
    C0 <- payoff
    C0_plot <- payoff
    
    payoff_plot <- NULL
    payoff_anti <- NULL
    
  }
  
  # the two variables are saves as data for the x option -> for plotting
  # numeric is needed
  if(type == "barrier" || type == "chooser"){
    C0 <- as.numeric(C0)
    C0_plot <- as.numeric(C0_plot)
    
    payoff <- as.numeric(payoff)
  }
  
  
  optionValue <- mean(C0)
  optionSD <- ifelse(M <= 1, 0, sd(C0))
  optionSE <- ifelse(M <= 1, 0, sd(C0)/sqrt(M))
  
  # return list
  returnlist <- list(C,
                     optionValue = optionValue, 
                     optionSD = optionSD,
                     optionSE = optionSE,
                     type = list(type,subType,optionType), 
                     C0 = C0,
                     C0_plot = C0_plot,
                     payoff = payoff,
                     payoff_plot = payoff_plot,
                     payoff_anti = payoff_anti,
                     stockPaths = c(stockPaths),
                     strikePrice = strikePrice,
                     # store all option specific input parameters for comparison
                     inputparams = list(asianRange = asianRange,
                                        intervalObs_asian = intervalObs_asian, 
                                        periodObs_asian = periodObs_asian,
                                        
                                        # barrier option
                                        barrier = barrier, 
                                        
                                        # binary option cashOrNothing
                                        binaryPayoff = binaryPayoff,
                                        
                                        # chooser
                                        chooser = chooser,
                                        choosingDate = choosingDate,
                                        chooserOption = chooserOption,
                                        chooserSubType = chooserSubType,
                                        
                                        # compound
                                        secondSampleSize = secondSampleSize,
                                        secondNumberTimeSteps = secondNumberTimeSteps,
                                        secondSkalarTimeSteps = secondSkalarTimeSteps,
                                        secondExerciseDate = T2,
                                        type2 = type2,
                                        secondStrikePrice  = secondStrikePrice,
                                        secondSubType = secondSubType,
                                        
                                        # forward start
                                        startDate = startDate,
                                        
                                        # gap
                                        strikePrice1 = strikePrice1,
                                        strikePrice2 = strikePrice2)
  )
  
  return(returnlist)  
}




#### functions not directly for calculation of the optionprices------------------

# set variable to check which tab is active to use in "check.Input"
set.active.option <- function(reactives, type){
  
  # only the list element handed over in type is set to true
  for(i in 1:length(reactives$check_active_option)){
    
    reactives$check_active_option[i] <- FALSE
    
  }
  
  reactives$check_active_option[type] <- TRUE
  
}

# check if the paramters used for the heston path plot are NA
check.Heston_na <- function(reactives, input){
  
  check <- TRUE
  warning <- NULL
  
  # Check general parameters that ar used for the heston plot
  if(is.na(input$stockPrice)){
    
    check <- FALSE
    warning <- "Stockprice mustn't be empty."
    
  } else if(is.na(input$volatility)){
    
    check <- FALSE
    warning <- "Volatility mustn't be empty."
    
  } else if(is.na(input$riskFreeRate)){
    
    check <- FALSE
    warning <- "Risk free rate mustn't be empty."
    
  } else if(input$dividend_YesOrNo == TRUE){
    
    if(is.na(input$dividendRate)){
      
      check <- FALSE
      warning <- "Dividend yield mustn't be empty."
      
    }
    
  }
  
  # Check Heston specific parameters
  if(is.na(input$kappa)){
    
    check <- FALSE
    warning <- "Kappa mustn't be empty."
    
  } else if(is.na(input$epsilon)){
    
    check <- FALSE
    warning <- "Epsilon mustn't be empty."
    
  } else if(is.na(input$rho)){
    
    check <- FALSE
    warning <- "Correlation mustn't be empty."
    
  } else if(is.na(input$theta)){
    
    check <- FALSE
    warning <- "Long time variance mustn't be empty."
    
  }
  
  # return list with logical check and the warning message
  returnlist <- list(check = check,
                     warning = warning)
  
  return(returnlist)
  
}

# check if the paramters used for the heston path plot are correct
check.Heston <- function(reactives, input){
  
  check <- TRUE
  warning <- NULL
  
  # try catches errors due to missing values in the if loop
  # this happens if one parameter is empty (deleted by hand in the app input)
  try(
    # Check general parameters that are used for the heston plot
    if(input$stockPrice <= 0){
      
      check <- FALSE
      warning <- "Stockprice has to be greater than zero."
      
    } else if(input$volatility < 0 || input$volatility > 100){
      
      check <- FALSE
      warning <- "Volatility has to be between zero and onehundred percent."
      
    } else if(input$riskFreeRate < -100 || input$riskFreeRate > 100){
      
      check <- FALSE
      warning <- "Risk free rate has to be between minus and plus onehundred percent."
      
    } else if(input$dividend_YesOrNo == TRUE){
      
      if(input$dividendRate < 0){
        
        check <- FALSE
        warning <- "The dividend yield must be greater or equal to zero."
        
      } else if(input$dividendRate > 100){
        
        check <- FALSE
        warning <- "The dividend yield must be below onehundred."
        
      }
      
    }
    , silent = TRUE
  )
  
  
  try(
    # Check Heston specific parameters
    if(input$kappa < 0){
      
      check <- FALSE
      warning <- "Kappa has to be greater than zero."
      
    } else if(input$epsilon < 0 || input$epsilon > 1){
      
      check <- FALSE
      warning <- "Epsilon has to be between zero and one."
      
    } else if(input$rho < -1 || input$rho > 1){
      
      check <- FALSE
      warning <- "Correlation hast to be between minus one and one."
      
    } else if(input$theta < 0 || input$theta > 100){
      
      check <- FALSE
      warning <- "Long time variance has to be between zero and onehundred percent."
      
    }
    , silent = TRUE
  )
  
  try(
    # check that the date of maturity is in the future
    if((input$hestonPathDate - Sys.Date()) <= 0){
      
      check <- FALSE
      warning <- "Date of maturity has to be in the future."
      
    } else{
    } 
    , silent = TRUE
  )
  
  
  # return list with logical check and the warning message
  returnlist <- list(check = check,
                     warning = warning)
  
  return(returnlist)
  
  
}


# function to detect if input is NA
check.Input_na <- function(reactives, input){
  
  check <- TRUE
  warning <- NULL
  
  #### American
  if(reactives$check_active_option[["american"]] == TRUE){
    
    if(is.na(input$timeSteps)){
      
      check <- FALSE
      warning <- "Number of timesteps mustn't be empty."
      
    }
    
    # if loop for date is similar for each option
   if(input$american_dateOrYear == "years"){
      
      if(is.na(input$timeToMaturity_american)){
        
        check <- FALSE
        warning <- "The time to maturity mustn't be empty."
        
      } else{
      }
      
    }
    
    # if loop with the other option specific input parameters
    if(is.na(input$strikePrice_american)){
      
      check <- FALSE
      warning <- "Strikeprice mustn't be empty."
      
    } else{
    }
    
    #### Asian
  } else if(reactives$check_active_option[["asian"]] == TRUE){
    
    if(is.na(input$skalar_timeSteps)){
      
      check <- FALSE
      warning <- "Number of timesteps mustn't be empty."
      
    }
    
    # if loop for date is similar for each option
    if(input$asian_dateOrYear == "years"){
      
      if(is.na(input$timeToMaturity_asian)){
        
        check <- FALSE
        warning <- "The time to maturity mustn't be empty."
        
      } else{
      }
      
    }
    # error handling of the Asian total
    if(input$observation_asian == "Total"){
      if(is.na(input$strikePrice_asian_total)){
        check <- FALSE
        warning <- "Strike price mustn't be empty."
      }
    } 
    # error handling Interval Asian
    else if(input$observation_asian == "Interval"){
      if(is.na(input$strikePrice_asian_interval)){
        check <- FALSE
        warning <- "Strike price mustn't be empty."
      }
      if(is.na(input$daysBetweenObs_asian)){
        check <- FALSE
        warning <- "Positive number mustn't be empty."
      }
    } 
    # error handling Period Asian
    else if(input$observation_asian == "Period"){
      if(is.na(input$strikePrice_asian_period)){
        check <- FALSE
        warning <- "Strike price mustn't be empty."
      }
      if(is.na(input$periodOfObs_asian)){
        check <- FALSE
        warning <- "Period of observation mustn't be empty."
      }
    }
    
    
    
    
    #### Barrier
  } else if(reactives$check_active_option[["barrier"]] == TRUE){
    
    if(is.na(input$timeSteps)){
      
      check <- FALSE
      warning <- "Number of timesteps mustn't be empty."
      
    }
    
    # if loop for date is similar for each option
    if(input$barrier_dateOrYear == "years"){
      
      if(is.na(input$timeToMaturity_barrier)){
        
        check <- FALSE
        warning <- "The time to maturity mustn't be empty."
        
      } else{
      }
      
    }
    
    # if loop with the other option specific input parameters
    if(is.na(input$strikePrice_barrier)){
      
      check <- FALSE
      warning <- "Strikeprice mustn't be empty."
      
    } else if(is.na(input$barrierValue)){
      
      check <- FALSE
      warning <- "Barrier mustn't be empty."
      
    } else{
      
    }
    
    #### Binary
  } else if(reactives$check_active_option[["binary"]] == TRUE){
    
    if(is.na(input$timeSteps)){
      
      check <- FALSE
      warning <- "Number of timesteps mustn't be empty."
      
    }
    
    # if loop for date is similar for each option
    if(input$binary_dateOrYear == "years"){
      
      if(is.na(input$timeToMaturity_binary)){
        
        check <- FALSE
        warning <- "The time to maturity mustn't be empty.."
        
      } else{
      }
      
    }
    
    # if loop with the other option specific input parameters
    if(input$subType_binary == "Asset-or-nothing"){
      
      if(is.na(input$strikePrice_binary_asset)){
        
        check <- FALSE
        warning <- "Strikeprice mustn't be empty."
        
      } else{
        
      }
      
      
    } else if(input$subType_binary == "Cash-or-nothing"){
      
      if(is.na(input$binaryPayoff)){
        
        check <- FALSE
        warning <- "Binary payoff mustn't be empty."
        
      } else if(is.na(input$strikePrice_binary_cash)){
        
        check <- FALSE
        warning <- "Strikeprice mustn't be empty."
        
      } else{
        
      }
      
    } else{
      
    }
    
    
    #### Chooser
  } else if(reactives$check_active_option[["chooser"]] == TRUE){
    
    # check if skalar timesteps is greater than zero
    if(is.na(input$skalar_timeSteps)){
      
      check <- FALSE
      warning <- "Number of timesteps mustn't be empty."
      
    }
    
    # if loop for maturity date and decision date
    if(input$chooser_dateOrYear == "date"){
      
      if(input$chooser_dateOrYear == "years"){
      
      reference_years <- as.numeric(difftime(input$chooserDecisionDate, Sys.Date(), unit="weeks"))/52.25
      
      if(is.na(input$timeToMaturity_chooser)){
        
        check <- FALSE
        warning <- "The time to maturity mustn't be empty."
        
      }
    }
    
    # if loop with the other option specific input parameters
    if(is.na(input$strikePrice_chooser)){
      
      check <- FALSE
      warning <- "Strikeprice mustn't be empty."
      
    } else if(input$chooser_option == "asian"){
      
      # error handling Interval Asian
      if(input$observation_asian_chooser == "Interval"){
        if(is.na(input$daysBetweenObs_asian_chooser)){
          check <- FALSE
          warning <- "Observations mustn't be empty."
        }
        
      } 
      # error handling Period Asian
      else if(input$observation_asian_chooser == "Period"){
        if(is.na(input$periodOfObs_asian_chooser)){
          check <- FALSE
          warning <- "Period of observation mustn't be empty."
        }
      }
      
    } else if(input$chooser_option == "barrier"){
      
      # if loop with the other option specific input parameters
      if(is.na(input$barrierValue_chooser)){
        
        check <- FALSE
        warning <- "Barrier mustn't be empty."
        
      } else{
        
      }
      
    } else if(input$chooser_option == "binary"){
      
      # if loop with the other option specific input parameters
      if(input$chooser_subtype_binary == "Cash-or-nothing"){
        
        if(is.na(input$binaryPayoff_chooser)){
          
          check <- FALSE
          warning <- "Binary payoff mustn't be empty."
          
        } else{
          
        }
        
      } else{
        
      }
      
    } else if(input$chooser_option == "european"){
      
      # no parameters to check
      
    } else if(input$chooser_option == "forwardStart"){
      
      # if loop for maturity date and start date
      if(input$chooser_dateOrYear == "date"){
        
        if((input$maturityDate_chooser - Sys.Date()) <= 0){
          
          check <- FALSE
          warning <- "Date of maturity has to be in the future."
          
        } else if((input$maturityDate_chooser - Sys.Date()) > 0){
          
          if((input$forwardStartDate_chooser-input$maturityDate_chooser) >= 0){
            
            check <- FALSE
            warning <- "Start date must be before the maturity date."
            
          } else if((input$forwardStartDate_chooser-Sys.Date()) <= 0){
            
            check <- FALSE
            warning <- "Start date must be after today."
            
          } else {
            
          }
          
        }
        
      } else if(input$chooser_dateOrYear == "years"){
        
        reference_years <- as.numeric(difftime(input$forwardStartDate_chooser, Sys.Date(), unit="weeks"))/52.25
        
        if(is.na(input$timeToMaturity_chooser)){
          
          check <- FALSE
          warning <- "The time to maturity mustn't be empty."
          
        } 
      }
      
      
      
    } else if(input$chooser_option == "gap"){
      
      if(is.na(input$triggerPrice_gap_chooser)){
        
        check <- FALSE
        warning <- "Trigger price mustn't be empty."
        
      } else{
      }
      
    } else if(input$chooser_option == "lookback"){
      
      # no parameters to check
      
    } else {
      
    }
      
    }
    
    #### Compound
  } else if(reactives$check_active_option[["compound"]] == TRUE){
    
    # if loop for date of the overlying option
    if(input$compound_first_dateOrYear == "years"){
      
      if(is.na(input$timeToMaturity_compound_first)){
        
        check <- FALSE
        warning <- "The time to maturity mustn't be empty."
        
      } else{
      }
      
    }
    
    # if loop for date of the underlying option
   if(input$compound_second_dateOrYear == "years"){
      
      if(is.na(input$timeToMaturity_compound_second)){
        
        check <- FALSE
        warning <- "The time to maturity mustn't be empty."
        
      } else{
      }
      
    }
    
    
    # if loop with the other option specific input parameters
    if(is.na(input$strikePrice_compound_first)){
      
      check <- FALSE
      warning <- "First strikeprice mustn't be empty."
      
    } else if(is.na(input$strikePrice_compound_second)){
      
      check <- FALSE
      warning <- "Second strikeprice mustn't be empty."
      
    } else if(is.na(input$sampleSize_compound)){
      
      check <- FALSE
      warning <- "The number of Monte Carlo samples for the overlying option mustn't be empty."
      
    } else if(is.na(input$sampleSize_compound_second)){
      
      check <- FALSE
      warning <- "The number of Monte Carlo samples for the underlying option mustn't be empty."
      
    } else if(input$compound_option_second == "asian"){
      
      # check if skalar timesteps is greater than zero
      if(is.na(input$skalar_timeSteps)){
        
        check <- FALSE
        warning <- "Number of timesteps mustn't be empty."
        
      }
      
      
      # error handling interval asian
      if(input$observation_asian_compound == "Interval"){
        if(is.na(input$daysBetweenObs_asian_compound)){
          check <- FALSE
          warning <- "Observations mustn't be empty."
        }
      } 
      # error handling period asian
      else if(input$observation_asian_compound == "Period"){
        if(is.na(input$periodOfObs_asian_compound)){
          check <- FALSE
          warning <- "Period of observation mustn't be empty."
        }
        #' @add period > maturity check einfügen
        
      }
      
    } else if(input$compound_option_second == "barrier"){
      
      if(is.na(input$timeSteps)){
        
        check <- FALSE
        warning <- "Number of timesteps mustn't be empty."
        
      }
      
      if(is.na(input$barrierValue_compound)){
        
        check <- FALSE
        warning <- "Barrier mustn't be empty."
        
      } else{
        
      }
      
    } else if(input$compound_option_second == "binary"){
      
      if(is.na(input$timeSteps)){
        
        check <- FALSE
        warning <- "Number of timesteps mustn't be empty."
        
      }
      
      # if loop with the other option specific input parameters
      if(input$compound_subtype_binary == "Cash-or-nothing"){
        
        if(is.na(input$binaryPayoff_compound)){
          
          check <- FALSE
          warning <- "Binary payoff mustn't be empty."
          
        } else{
          
        }
        
      } else{
        
      }
      
    } else if(input$compound_option_second == "european"){
      
      # no parameters to check
      
    } else if(input$compound_option_second == "forwardStart"){
      
      if(is.na(input$timeSteps)){
        
        check <- FALSE
        warning <- "Number of timesteps mustn't be empty."
        
      }
      
      
      
    } else if(input$compound_option_second == "gap"){
      
      if(is.na(input$timeSteps)){
        
        check <- FALSE
        warning <- "Number of timesteps mustn't be empty."
        
      }
      
      # if loop with the other option specific input parameters
      if(is.na(input$triggerPrice_gap_compound)){
        
        check <- FALSE
        warning <- "Trigger price mustn't be empty."
        
      } else{
      }
      
    } else if(input$compound_option_second == "lookback"){
      
      if(is.na(input$skalar_timeSteps)){
        
        check <- FALSE
        warning <- "Number of timesteps mustn't be empty."
        
      }
      
    } else {
      
    }
    
    
    #### European
  } else if(reactives$check_active_option[["european"]] == TRUE){
    
    # if loop for date is similar for each option
    if(input$european_dateOrYear == "years"){
      
      if(is.na(input$timeToMaturity_european)){
        
        check <- FALSE
        warning <- "The time to maturity mustn't be empty."
        
      } else{
      }
      
    }
    
    # if loop with the other option specific input parameters
    if(is.na(input$strikePrice_european)){
      
      check <- FALSE
      warning <- "Strikeprice mustn't be empty."
      
    } else{
    }
    
    
    #### Forward Start
  } else if(reactives$check_active_option[["forwardStart"]] == TRUE){
    
    if(is.na(input$timeSteps)){
      
      check <- FALSE
      warning <- "Number of timesteps mustn't be empty."
      
    }
    
    # if loop for maturity date and start date
    if(input$forwardStart_dateOrYear == "years"){
      
      reference_years <- as.numeric(difftime(input$forwardStartDate, Sys.Date(), unit="weeks"))/52.25
      
      if(is.na(input$timeToMaturity_forwardStart)){
        
        check <- FALSE
        warning <- "The time to maturity mustn't be empty."
        
      }
    }
    
    # if loop with the other option specific input parameters
    if(is.na(input$strikePrice_forwardStart)){
      
      check <- FALSE
      warning <- "Strikeprice mustn't be empty."
      
    } else{
    }
    
    #### Gap
  } else if(reactives$check_active_option[["gap"]] == TRUE){
    
    if(is.na(input$timeSteps)){
      
      check <- FALSE
      warning <- "Number of timesteps mustn't be empty."
      
    }
    
    # if loop for date is similar for each option
    if(input$gap_dateOrYear == "years"){
      
      if(is.na(input$timeToMaturity_gap)){
        
        check <- FALSE
        warning <- "The time to maturity mustn't be empty."
        
      } else{
      }
      
    }
    
    # if loop with the other option specific input parameters
    if(is.na(input$strikePrice_gap1)){
      
      check <- FALSE
      warning <- "First strikeprice mustn't be empty."
      
    } else if(is.na(input$strikePrice_gap2)){
      
      check <- FALSE
      warning <- "Trigger price mustn't be empty."
      
    } else{
    }
    
    #### Lookback
  } else if(reactives$check_active_option[["lookback"]] == TRUE){
    
    if(is.na(input$timeSteps)){
      
      check <- FALSE
      warning <- "Number of timesteps mustn't be empty."
      
    }
    
    # if loop for date is similar for each option
    if(input$lookback_dateOrYear == "years"){
      
      if(is.na(input$timeToMaturity_lookback)){
        
        check <- FALSE
        warning <- "The time to maturity mustn't be empty."
        
      } else{
      }
      
    }
    
    # if loop with the other option specific input parameters
    if(is.na(input$strikePrice_lookback)){
      
      check <- FALSE
      warning <- "Strikeprice mustn't be empty."
      
    } else{
      
    }
    
  }
  
  
  
  
  # Check for parameters that ar similar for Black-Scholes and Heston
  if(is.na(input$stockPrice)){
    
    check <- FALSE
    warning <- "Stockprice mustn't be empty."
    
  } else if(is.na(input$volatility)){
    
    check <- FALSE
    warning <- "Volatility mustn't be empty."
    
  } else if(is.na(input$riskFreeRate)){
    
    check <- FALSE
    warning <- "Risk free rate mustn't be empty."
    
  } else if(input$dividend_YesOrNo == TRUE){
    
    if(is.na(input$dividendRate)){
      
      check <- FALSE
      warning <- "Dividend yield mustn't be empty."
      
    }
    
  }
  
  # Check Heston specific parameters
  if(input$simulationModel == "Heston"){
    if(is.na(input$kappa)){
      
      check <- FALSE
      warning <- "Kappa mustn't be empty."
      
    } else if(is.na(input$epsilon)){
      
      check <- FALSE
      warning <- "Epsilon mustn't be empty."
      
    } else if(is.na(input$rho)){
      
      check <- FALSE
      warning <- "Correlation mustn't be empty."
      
    } else if(is.na(input$theta)){
      
      check <- FALSE
      warning <- "Long time variance mustn't be empty."
      
    }
  }
  
  # return list with logical check and the warning message
  returnlist <- list(check = check,
                     warning = warning)
  
  return(returnlist)
  
}


# check if all input parameters are correct
check.Input <- function(reactives,
                        input){
  
  check <- TRUE
  warning <- NULL
  
  
  #### American
  if(reactives$check_active_option[["american"]] == TRUE){
    
    # try catches errors due to missing values in the if loop
    # this happens if one parameter is empty (deleted by hand in the app input)
    try ( 
    # check if timesteps is greater than zero
    if(input$timeSteps <= 0){
      
      check <- FALSE
      warning <- "Number of timesteps must be greater than zero."
      
    }
    , silent = TRUE
    )
    
    try(
    
    # if loop for date is similar for each option
    if(input$american_dateOrYear == "date"){
      
      if((input$maturityDate_american - Sys.Date()) <= 0){
        
        check <- FALSE
        warning <- "Date of maturity has to be in the future."
        
      } else{
      }
      
    } else if(input$american_dateOrYear == "years"){
      
      if(input$timeToMaturity_american <= 0){
        
        check <- FALSE
        warning <- "The time to maturity has to be greater than zero."
        
      } else{
      }
      
    }
    , silent = TRUE
    )
    
    try(
    
    # if loop with the other option specific input parameters
    if(input$strikePrice_american <= 0){
      
      check <- FALSE
      warning <- "Strikeprice has to be greater than zero."
      
    } else{
    }
    , silent = TRUE
    )
    
    #### Asian
  } else if(reactives$check_active_option[["asian"]] == TRUE){
    try(
    # check skalar timesteps
    if(input$skalar_timeSteps <= 0){
      
      check <- FALSE
      warning <- "Number of timesteps must be greater than zero."
      
    }
    , silent = TRUE
    )
    try(
    if(input$skalar_timeSteps %%1 != 0){
      check <- FALSE
      warning <- "Number of time steps per day has to be an integer."
    }
    , silent = TRUE
    )
    try(
    # if loop for date is similar for each option
    if(input$asian_dateOrYear == "date"){
      
      if((input$maturityDate_asian - Sys.Date()) <= 0){
        
        check <- FALSE
        warning <- "Date of maturity has to be in the future."
        
      } else{
      }
    } else if(input$asian_dateOrYear == "years"){
      
      if(input$timeToMaturity_asian <= 0){
        
        check <- FALSE
        warning <- "The time to maturity has to be greater than zero."
        
      } else{
      }
      
    }
    , silent = TRUE
    )
    try(
    # error handling asian total
    if(input$observation_asian == "Total"){
      if(input$strikePrice_asian_total <= 0){
        check <- FALSE
        warning <- "Strike price has to be greater than zero."
      }
    } 
    # error handling asian interval
    else if(input$observation_asian == "Interval"){
      if(input$strikePrice_asian_interval <= 0){
        check <- FALSE
        warning <- "Strike price has to be greater than zero."
      }
      if(input$daysBetweenObs_asian <= 0){
        check <- FALSE
        warning <- "Positive number of days between observations required."
      }
      if(input$daysBetweenObs_asian %% 1 != 0){
        check <- FALSE
        warning <- "Interval between observation needs to be an integer."
      }
    } 
    # error handling asian period
    else if(input$observation_asian == "Period"){
      if(input$strikePrice_asian_period <= 0){
        check <- FALSE
        warning <- "Strike price has to be greater than zero."
      }
      if(input$periodOfObs_asian <= 1){
        check <- FALSE
        warning <- "Period of observation must be greater than one."
      }
      if(input$periodOfObs_asian %% 1 != 0){
        check <- FALSE
        warning <- "Period of observation needs to be an integer."
      }
      if(input$asian_dateOrYear == "date"){
        if((input$maturityDate_asian - Sys.Date()) < input$periodOfObs_asian){
          check <- FALSE
          warning <- "Period of observation is longer than time to expiration."
        }
      } else if (input$asian_dateOrYear == "years"){
        if(floor(365 * input$timeToMaturity_asian) < input$periodOfObs_asian){
          check <- FALSE
          warning <- "Period of observation is longer than time to expiration."
        }
      }
    }
    , silent = TRUE
    )
    
    
    
    #### Barrier
  } else if(reactives$check_active_option[["barrier"]] == TRUE){
    try(
    # check if timesteps is greater than zero
    if(input$timeSteps <= 0){
      
      check <- FALSE
      warning <- "Number of timesteps must be greater than zero."
      
    }
    , silent = TRUE
    )
    try(
    
    # if loop for date is similar for each option
    if(input$barrier_dateOrYear == "date"){
      
      if((input$maturityDate_barrier - Sys.Date()) <= 0){
        
        check <- FALSE
        warning <- "Date of maturity has to be in the future."
        
      } else{
      }
      
    } else if(input$barrier_dateOrYear == "years"){
      
      if(input$timeToMaturity_barrier <= 0){
        
        check <- FALSE
        warning <- "The time to maturity has to be greater than zero."
        
      } else{
      }
      
    }
    , silent = TRUE
    )
    try(
    
    # if loop with the other option specific input parameters
    if(input$strikePrice_barrier <= 0){
      
      check <- FALSE
      warning <- "Strikeprice has to be greater than zero."
      
    } else if(input$barrierValue <= 0){
      
      check <- FALSE
      warning <- "Barrier has to be greater than zero."
      
    } else{
      
    }
    , silent = TRUE
    )
    
    #### Binary
  } else if(reactives$check_active_option[["binary"]] == TRUE){
    
    try(
    # check if timesteps is greater than zero
    if(input$timeSteps <= 0){
      
      check <- FALSE
      warning <- "Number of timesteps must be greater than zero."
      
    }
    , silent = TRUE
    )
    try(
    # if loop for date is similar for each option
    if(input$binary_dateOrYear == "date"){
      
      if((input$maturityDate_binary - Sys.Date()) <= 0){
        
        check <- FALSE
        warning <- "Date of maturity has to be in the future."
        
      } else{
      }
      
    } else if(input$binary_dateOrYear == "years"){
      
      if(input$timeToMaturity_binary <= 0){
        
        check <- FALSE
        warning <- "The time to maturity has to be greater than zero."
        
      } else{
      }
      
    }
    , silent = TRUE
    )
    
    try(
    # if loop with the other option specific input parameters
    if(input$subType_binary == "Asset-or-nothing"){
      
      if(input$strikePrice_binary_asset <= 0){
        
        check <- FALSE
        warning <- "Strikeprice has to be greater than zero."
        
      } else{
        
      }
      
      
    } else if(input$subType_binary == "Cash-or-nothing"){
      
      if(input$binaryPayoff <= 0){
        
        check <- FALSE
        warning <- "Binary payoff must be greater than zero."
        
      } else if(input$strikePrice_binary_cash <= 0){
        
        check <- FALSE
        warning <- "Strikeprice has to be greater than zero."
        
      } else{
        
      }
      
    } else{
      
    }
    , silent = TRUE
    )
    
    
    #### Chooser
  } else if(reactives$check_active_option[["chooser"]] == TRUE){
    
    try(
    # check if skalar timesteps is greater than zero
    if(input$skalar_timeSteps <= 0){
      
      check <- FALSE
      warning <- "Number of timesteps must be greater than zero."
      
    }
    , silent = TRUE
    )
    try(
    if(input$skalar_timeSteps %%1 != 0){
      check <- FALSE
      warning <- "Number of time steps per day has to be an integer."
    }
    , silent = TRUE
    )
    try(
    # if loop for maturity date and decision date
    if(input$chooser_dateOrYear == "date"){
      
      if((input$maturityDate_chooser - Sys.Date()) <= 0){
        
        check <- FALSE
        warning <- "Date of maturity has to be in the future."
        
      } else if((input$maturityDate_chooser- Sys.Date()) > 0){
        
        if((input$chooserDecisionDate-input$maturityDate_chooser) >= 0){
          
          check <- FALSE
          warning <- "Decision date must be before the maturity date."
          
        } else if((input$chooserDecisionDate-Sys.Date()) <= 0){
          
          check <- FALSE
          warning <- "Decision date must be after today."
          
        } else{
          
        }
        
      }
      
    } else if(input$chooser_dateOrYear == "years"){
      
      reference_years <- as.numeric(difftime(input$chooserDecisionDate, Sys.Date(), unit="weeks"))/52.25
      
      if(input$timeToMaturity_chooser <= 0){
        
        check <- FALSE
        warning <- "The time to maturity has to be greater than zero."
        
      } else if(input$timeToMaturity_chooser > 0){
        
        if((reference_years - input$timeToMaturity_chooser) >= 0){
          
          check <- FALSE
          warning <- "Decision date must be before the maturity date."
          
        } else if((input$chooserDecisionDate-Sys.Date()) <= 0){
          
          check <- FALSE
          warning <- "Decision date must be after today."
          
        } else{
          
        }
      }
    }
    , silent = TRUE
    )
    
    try(
    # if loop with the other option specific input parameters
    if(input$strikePrice_chooser <= 0){
      
      check <- FALSE
      warning <- "Strikeprice has to be greater than zero."
      
    } else if(input$chooser_option == "asian"){
      
      # Behandlung Interval Asian
      if(input$observation_asian_chooser == "Interval"){
        if(input$daysBetweenObs_asian_chooser <= 0){
          check <- FALSE
          warning <- "Positive number of days between observations required."
        }
        if(input$daysBetweenObs_asian_chooser %% 1 != 0){
          check <- FALSE
          warning <- "Interval between observation needs to be an integer."
        }
      } 
      # Behandlung Period Asian
      else if(input$observation_asian_chooser == "Period"){
        if(input$periodOfObs_asian_chooser <= 1){
          check <- FALSE
          warning <- "Period of observation must be greater than one."
        }
        if(input$periodOfObs_asian_chooser %% 1 != 0){
          check <- FALSE
          warning <- "Period of observation needs to be an integer."
        }
        if(input$chooser_dateOrYear == "date"){
          if((input$maturityDate_chooser - Sys.Date()) < input$periodOfObs_asian_chooser){
            check <- FALSE
            warning <- "Period of observation is longer than time to expiration."
          }
          if(input$chooserDecisionDate - as.Date((as.numeric(input$maturityDate_chooser) - input$periodOfObs_asian_chooser)) < 0){
            check <- FALSE
            warning <- "Decision date must be within the period."
          }
        } else if (input$chooser_dateOrYear == "years"){
          if(floor(365 * input$timeToMaturity_chooser) < input$periodOfObs_asian_chooser){
            check <- FALSE
            warning <- "Period of observation is longer than time to expiration."
          }
          if(input$chooserDecisionDate - as.Date(as.numeric(Sys.Date()) + (floor(365 * input$timeToMaturity_chooser) - input$periodOfObs_asian_chooser)) < 0){
            check <- FALSE
            warning <- "Decision date must be within the period."
          }
        }
      }
      
    } else if(input$chooser_option == "barrier"){
      
      # if loop with the other option specific input parameters
      if(input$barrierValue_chooser <= 0){
        
        check <- FALSE
        warning <- "Barrier has to be greater than zero."
        
      } else{
        
      }
      
    } else if(input$chooser_option == "binary"){
      
      # if loop with the other option specific input parameters
      if(input$chooser_subtype_binary == "Cash-or-nothing"){
        
        if(input$binaryPayoff_chooser <= 0){
          
          check <- FALSE
          warning <- "Binary payoff must be greater than zero."
          
        } else{
          
        }
        
      } else{
        
      }
      
    } else if(input$chooser_option == "european"){
      
      # no parameters to check
      
    } else if(input$chooser_option == "forwardStart"){
      
      # if loop for maturity date and start date
      if(input$chooser_dateOrYear == "date"){

        if((input$maturityDate_chooser - Sys.Date()) <= 0){

          check <- FALSE
          warning <- "Date of maturity has to be in the future."

        } else if((input$maturityDate_chooser - Sys.Date()) > 0){

          if((input$forwardStartDate_chooser-input$maturityDate_chooser) >= 0){

            check <- FALSE
            warning <- "Start date must be before the maturity date."

          } else if((input$forwardStartDate_chooser-Sys.Date()) <= 0){

            check <- FALSE
            warning <- "Start date must be after today."

          } else {

          }

        }

      } else if(input$chooser_dateOrYear == "years"){

        reference_years <- as.numeric(difftime(input$forwardStartDate_chooser, Sys.Date(), unit="weeks"))/52.25

        if(input$timeToMaturity_chooser <= 0){

          check <- FALSE
          warning <- "The time to maturity has to be greater than zero."

        } else if(input$timeToMaturity_chooser > 0){

          if((reference_years - input$timeToMaturity_chooser) >= 0){

            check <- FALSE
            warning <- "Start date must be before the maturity date."

          } else if((input$forwardStartDate_chooser-Sys.Date()) <= 0){

            check <- FALSE
            warning <- "Start date must be after today."

          } else {

          }
        }
      }
      
      if(input$chooserDecisionDate < input$forwardStartDate_chooser){
        
        check <- FALSE
        warning <- "Decision date must be after start date."
        
      }
      
      
    } else if(input$chooser_option == "gap"){
      
      # if loop with the other option specific input parameters
      if(input$triggerPrice_gap_chooser <= 0){
        
        check <- FALSE
        warning <- "Trigger price has to be greater than zero."
        
      } else{
      }
      
    } else if(input$chooser_option == "lookback"){
      
      # no parameters to check
      
    } else {
      
    }
    , silent = TRUE
    )
    
    
    #### Compound
  } else if(reactives$check_active_option[["compound"]] == TRUE){
    try(
    # if loop date overlying option
    if(input$compound_first_dateOrYear == "date"){
      
      if((input$maturityDate_compound_first - Sys.Date()) <= 0){
        
        check <- FALSE
        warning <- "Date of maturity has to be in the future."
        
      } else{
      }
      
    } else if(input$compound_first_dateOrYear == "years"){
      
      if(input$timeToMaturity_compound_first <= 0){
        
        check <- FALSE
        warning <- "The time to maturity has to be greater than zero."
        
      } else{
      }
      
    }
    , silent = TRUE
    )
    try(
    # if loop date underlying option
    if(input$compound_second_dateOrYear == "date"){

      if(input$compound_first_dateOrYear == "date"){

        if((input$maturityDate_compound_second - input$maturityDate_compound_first) <= 0){

          check <- FALSE
          warning <- "Date of maturity has to be past the maturity date of the overlying option."

        } else{
        }

      } else if(input$compound_first_dateOrYear == "years"){

        maturityDate_reference_compound <- as.Date(Sys.Date() + floor(365 * input$timeToMaturity_compound_first))

        if((input$maturityDate_compound_second - maturityDate_reference_compound) <= 0){

          check <- FALSE
          warning <- "Date of maturity has to be past the maturity date of the overlying option."

        }

      }

    } else if(input$compound_second_dateOrYear == "years"){

      if(input$timeToMaturity_compound_second <= 0){

        check <- FALSE
        warning <- "The time to maturity has to be greater than zero."

      } else{
      }

    }
    , silent = TRUE
    )
    try(
    
    # if loop with the other option specific input parameters
    if(input$strikePrice_compound_first <= 0){
      
      check <- FALSE
      warning <- "First strikeprice has to be greater than zero."
      
    } else if(input$strikePrice_compound_second <= 0){
      
      check <- FALSE
      warning <- "Second strikeprice has to be greater than zero."
      
    } else if(input$sampleSize_compound <= 0){
      
      check <- FALSE
      warning <- "The number of Monte Carlo samples for the overlying option is zero or less."
      
    } else if(input$sampleSize_compound > 10000){
      
      check <- FALSE
      warning <- "The number of Monte Carlo samples for the overlying option is bigger than 10000."
      
    } else if(input$sampleSize_compound_second <= 0){
      
      check <- FALSE
      warning <- "The number of Monte Carlo samples for the underlying option is zero or less."
      
    }  else if(input$sampleSize_compound_second > 10000){
      
      check <- FALSE
      warning <- "The number of Monte Carlo samples for the underlying option is bigger than 10000."
      
    } else if(input$compound_option_second == "asian"){
      
      # check if skalar timesteps is greater than zero
      if(input$skalar_timeSteps <= 0){
        
        check <- FALSE
        warning <- "Number of timesteps must be greater than zero."
        
      } else if(input$skalar_timeSteps %%1 != 0){
        check <- FALSE
        warning <- "Number of time steps per day has to be an integer."
      }
      
      
      # error handling asian interval
      if(input$observation_asian_compound == "Interval"){
        if(input$daysBetweenObs_asian_compound <= 0){
          check <- FALSE
          warning <- "Positive number of days between observations required."
        }
        if(input$daysBetweenObs_asian_compound %% 1 != 0){
          check <- FALSE
          warning <- "Interval between observation needs to be an integer."
        }
      } 
      # error handling asian period
      else if(input$observation_asian_compound == "Period"){
        if(input$periodOfObs_asian_compound <= 1){
          check <- FALSE
          warning <- "Period of observation must be greater than one."
        }
        if(input$periodOfObs_asian_compound %% 1 != 0){
          check <- FALSE
          warning <- "Period of observation needs to be an integer."
        }
        
        if(input$compound_second_dateOrYear == "date"){
          if((input$maturityDate_compound_second - Sys.Date()) < input$periodOfObs_asian_compound){
            check <- FALSE
            warning <- "Period of observation is longer than time to expiration."
          }
        } else if (input$compound_second_dateOrYear == "years"){
          if(floor(365 * input$timeToMaturity_compound_second) < input$periodOfObs_asian_compound){
            check <- FALSE
            warning <- "Period of observation is longer than time to expiration."
          }
        }
      }
      
    } else if(input$compound_option_second == "barrier"){
      
      # check if timesteps is greater than zero
      if(input$timeSteps <= 0){
        
        check <- FALSE
        warning <- "Number of timesteps must be greater than zero."
        
      }
      
      # if loop with the other option specific input parameters
      if(input$barrierValue_compound <= 0){
        
        check <- FALSE
        warning <- "Barrier has to be greater than zero."
        
      } else{
        
      }
      
    } else if(input$compound_option_second == "binary"){
      
      # check if timesteps is greater than zero
      if(input$timeSteps <= 0){
        
        check <- FALSE
        warning <- "Number of timesteps must be greater than zero."
        
      }
      
      # if loop with the other option specific input parameters
      if(input$compound_subtype_binary == "Cash-or-nothing"){
        
        if(input$binaryPayoff_compound <= 0){
          
          check <- FALSE
          warning <- "Binary payoff must be greater than zero."
          
        } else{
          
        }
        
      } else{
        
      }
      
    } else if(input$compound_option_second == "european"){
      
      # no parameters to check
      
    } else if(input$compound_option_second == "forwardStart"){
      
      # check if timesteps is greater than zero
      if(input$timeSteps <= 0){
        
        check <- FALSE
        warning <- "Number of timesteps must be greater than zero."
        
      }
      
      # if loop for forwardstartdate and maturity of the first and second option
      if(input$compound_first_dateOrYear == "date"){
        
        if((input$forwardStartDate_compound-input$maturityDate_compound_first) < 0){
          
          check <- FALSE
          warning <- "Start date must be after the maturity date of the overlying option."
          
        }
        
        if(input$compound_second_dateOrYear == "date"){
          
          if((input$forwardStartDate_compound-input$maturityDate_compound_second) >= 0){
            
            check <- FALSE
            warning <- "Start date must be before the maturity date."
            
          }
          
        } else if(input$compound_second_dateOrYear == "years"){
          
          reference_years_compound_forwardStart <- as.numeric(difftime(input$forwardStartDate_compound, input$maturityDate_compound_first, unit="weeks"))/52.25
          
          if((reference_years_compound_forwardStart - input$timeToMaturity_compound_second) >= 0){
            
            check <- FALSE
            warning <- "Start date must be before the maturity date of the underlying option."
            
          }
          
        }
        
      } else if(input$compound_first_dateOrYear == "years"){
        
        maturityDate_reference_compound_forwardStart <- as.Date(Sys.Date() + floor(365 * input$timeToMaturity_compound_first))
        
        if((input$forwardStartDate_compound-maturityDate_reference_compound_forwardStart) < 0){
          
          check <- FALSE
          warning <- "Start date must be after the maturity date of the overlying option."
          
        }
        
        if(input$compound_second_dateOrYear == "date"){
          
          if((input$forwardStartDate_compound-input$maturityDate_compound_second) >= 0){
            
            check <- FALSE
            warning <- "Start date must be before the maturity date."
            
          }
          
        } else if(input$compound_second_dateOrYear == "years"){
          
          reference_years_compound_forwardStart <- as.numeric(difftime(input$forwardStartDate_compound, maturityDate_reference_compound_forwardStart, unit="weeks"))/52.25
          
          if((reference_years_compound_forwardStart - input$timeToMaturity_compound_second) >= 0){
            
            check <- FALSE
            warning <- "Start date must be before the maturity date of the underlying option."
            
          }
          
        }
        
      }
      
      
    } else if(input$compound_option_second == "gap"){
      
      # check if timesteps is greater than zero
      if(input$timeSteps <= 0){
        
        check <- FALSE
        warning <- "Number of timesteps must be greater than zero."
        
      }
      
      # if loop with the other option specific input parameters
      if(input$triggerPrice_gap_compound <= 0){
        
        check <- FALSE
        warning <- "Trigger price has to be greater than zero."
        
      } else{
      }
      
    } else if(input$compound_option_second == "lookback"){
      
      # check if skalar timesteps is greater than zero
      if(input$skalar_timeSteps <= 0){
        
        check <- FALSE
        warning <- "Number of timesteps must be greater than zero."
        
      } else if(input$skalar_timeSteps %%1 != 0){
        check <- FALSE
        warning <- "Number of time steps per day has to be an integer."
      }
      
    } else {
      
    }
    , silent = TRUE
    )
    
    #### European
  } else if(reactives$check_active_option[["european"]] == TRUE){
    try(
    # if loop for date is similar for each option
    if(input$european_dateOrYear == "date"){
      
      if((input$maturityDate_european - Sys.Date()) <= 0){
        
        check <- FALSE
        warning <- "Date of maturity has to be in the future."
        
      } else{
      }
      
    } else if(input$european_dateOrYear == "years"){
      
      if(input$timeToMaturity_european <= 0){
        
        check <- FALSE
        warning <- "The time to maturity has to be greater than zero."
        
      } else{
      }
      
    }
    , silent = TRUE
    )
    try(
    # if loop with the other option specific input parameters
    if(input$strikePrice_european <= 0){
      
      check <- FALSE
      warning <- "Strikeprice has to be greater than zero."
      
    } else{
    }
    , silent = TRUE
    )
    
    #### Forward Start
  } else if(reactives$check_active_option[["forwardStart"]] == TRUE){
    try(
    # check if timesteps is greater than zero
    if(input$timeSteps <= 0){
      
      check <- FALSE
      warning <- "Number of timesteps must be greater than zero."
      
    }
    , silent = TRUE
    )
    try(
    
    # if loop for maturity date and start date
    if(input$forwardStart_dateOrYear == "date"){
      
      if((input$maturityDate_forwardStart - Sys.Date()) <= 0){
        
        check <- FALSE
        warning <- "Date of maturity has to be in the future."
        
      } else if((input$maturityDate_forwardStart - Sys.Date()) > 0){
        
        if((input$forwardStartDate-input$maturityDate_forwardStart) >= 0){
          
          check <- FALSE
          warning <- "Start date must be before the maturity date."
          
        } else if((input$forwardStartDate-Sys.Date()) <= 0){
          
          check <- FALSE
          warning <- "Start date must be after today."
          
        } else{
          
        }
        
      }
      
    } else if(input$forwardStart_dateOrYear == "years"){
      
      reference_years <- as.numeric(difftime(input$forwardStartDate, Sys.Date(), unit="weeks"))/52.25
      
      if(input$timeToMaturity_forwardStart <= 0){
        
        check <- FALSE
        warning <- "The time to maturity has to be greater than zero."
        
      } else if(input$timeToMaturity_forwardStart > 0){
        
        if((reference_years - input$timeToMaturity_forwardStart) >= 0){
          
          check <- FALSE
          warning <- "Start date must be before the maturity date."
          
        } else if((input$forwardStartDate-Sys.Date()) <= 0){
          
          check <- FALSE
          warning <- "Start date must be after today."
          
        } else{
          
        }
      }
    }
    , silent = TRUE
    )
    try(
    # if loop with the other option specific input parameters
    if(input$strikePrice_forwardStart <= 0){
      
      check <- FALSE
      warning <- "Strikeprice has to be greater than zero."
      
    } else{
    }
    , silent = TRUE
    )
    #### Gap
  } else if(reactives$check_active_option[["gap"]] == TRUE){
    try(
    # check if timesteps is greater than zero
    if(input$timeSteps <= 0){
      
      check <- FALSE
      warning <- "Number of timesteps must be greater than zero."
      
    }
    , silent = TRUE
    )
    try(
    # if loop for date is similar for each option
    if(input$gap_dateOrYear == "date"){
      
      if((input$maturityDate_gap - Sys.Date()) <= 0){
        
        check <- FALSE
        warning <- "Date of maturity has to be in the future."
        
      } else{
      }
      
    } else if(input$gap_dateOrYear == "years"){
      
      if(input$timeToMaturity_gap <= 0){
        
        check <- FALSE
        warning <- "The time to maturity has to be greater than zero."
        
      } else{
      }
      
    }
    , silent = TRUE
    )
    try(
    # if loop with the other option specific input parameters
    if(input$strikePrice_gap1 <= 0){
      
      check <- FALSE
      warning <- "First strikeprice has to be greater than zero."
      
    } else if(input$strikePrice_gap2 <= 0){
      
      check <- FALSE
      warning <- "Trigger price has to be greater than zero."
      
    } else{
    }
    , silent = TRUE
    )
    
    #### Lookback
  } else if(reactives$check_active_option[["lookback"]] == TRUE){
    try(
    # check if timesteps is greater than zero
    if(input$timeSteps <= 0){
      
      check <- FALSE
      warning <- "Number of timesteps must be greater than zero."
      
    }
    , silent = TRUE
    )
    try(
    
    # if loop for date is similar for each option
    if(input$lookback_dateOrYear == "date"){
      
      if((input$maturityDate_lookback - Sys.Date()) <= 0){
        
        check <- FALSE
        warning <- "Date of maturity has to be in the future."
        
      } else{
      }
      
    } else if(input$lookback_dateOrYear == "years"){
      
      if(input$timeToMaturity_lookback <= 0){
        
        check <- FALSE
        warning <- "The time to maturity has to be greater than zero."
        
      } else{
      }
      
    }
    , silent = TRUE
    )
    try(
    # if loop with the other option specific input parameters
    if(input$strikePrice_lookback <= 0){
      
      check <- FALSE
      warning <- "Strikeprice has to be greater than zero."
      
    } else{
      
    }
    , silent = TRUE
    )
    
  }
  
  
  
  try(
  # Check for parameters that ar similar for Black-Scholes and Heston
  if(input$stockPrice <= 0){
    
    check <- FALSE
    warning <- "Stockprice has to be greater than zero."
    
  } else if(input$volatility < 0 || input$volatility > 100){
    
    check <- FALSE
    warning <- "Volatility has to be between zero and onehundred percent."
    
  } else if(input$riskFreeRate < -100 || input$riskFreeRate > 100){
    
    check <- FALSE
    warning <- "Risk free rate has to be between minus and plus onehundred percent."
    
  } else if(input$dividend_YesOrNo == TRUE){
    
    if(input$dividendRate < 0){
      
      check <- FALSE
      warning <- "The dividend yield must be greater or equal to zero."
      
    } else if(input$dividendRate > 100){
      
      check <- FALSE
      warning <- "The dividend yield must be below onehundred."
      
    }
    
  }
  , silent = TRUE
  )
  
  
  try(
    
    # in case of the compound option the "input$sampleSize is not used"
    if(input$tab != "compound"){
      
      if(input$sampleSize == 0){
        
        check <- FALSE
        warning <- "The number of Monte Carlo samples is zero."
        
      }
      
    } 
    , silent = TRUE
  )
  
  
  
  try(
  # Check Heston specific parameters
  if(input$simulationModel == "Heston"){
    if(input$kappa < 0){
      
      check <- FALSE
      warning <- "Kappa has to be greater than zero."
      
    } else if(input$epsilon < 0 || input$epsilon > 1){
      
      check <- FALSE
      warning <- "Epsilon has to be between zero and one."
      
    } else if(input$rho < -1 || input$rho > 1){
      
      check <- FALSE
      warning <- "Correlation hast to be between minus one and one."
      
    } else if(input$theta < 0 || input$theta > 100){
      
      check <- FALSE
      warning <- "Long time variance has to be between zero and onehundred percent."
      
    }
  }
  , silent = TRUE
  )
  
  # return list with logical check and the warning message
  returnlist <- list(check = check,
                     warning = warning)
  
  return(returnlist)
  
}



# function to create warning if the input parameters have changed and the displayed info/plot is not up to update
create.warning <- function(reactives, input, type){
  
  warning <- NULL
  
  warning_dummy <- HTML(paste(tags$span(
    "Attention: The option calculation is not up to date. Please recalculate.",
      style = "color:#ff0000; font-size:20px; font-weight: bold;"
    )
  ))
  
  # variable to help setting the reactive warning variable to "TRUE"
  # because of the reactive behavior the reactive variable herself cannot be used
  reactive_warning_dummy <- FALSE
  
  if(!is.null(reactives$result[[type]])){
    
    switch(type,
           "american"={
             
             # try catches errors due to missing values in the if loop
             # this happens if one parameter is empty (deleted by hand in the app input)
             try (
             # save current maturity
             if(input$american_dateOrYear == "date"){
               maturity_reference <- as.numeric(difftime(input$maturityDate_american, Sys.Date(), unit="weeks"))/52.25 
             } else if(input$american_dateOrYear == "years"){
               if(input$timeToMaturity_american) {
                 maturity_reference <- input$timeToMaturity_american
               }
             }
             , silent = TRUE
             )
             
             try(
             # check for option specific input parameters
             if(input$strikePrice_american != reactives[["result"]][[type]][["strikePrice"]] ||
                input$timeSteps != reactives[["result"]][[type]][["stockPaths"]][["numberTimeSteps"]] ||
                maturity_reference != reactives[["result"]][[type]][["stockPaths"]][["maturity"]] ||
                input$optionType_american != reactives[["result"]][[type]][["type"]][[3]]){
               
               reactive_warning_dummy <- TRUE
               warning <- warning_dummy
               
             }
             , silent = TRUE
             )
             
           },
           "asian"={
             try (
             # save current maturity
             if(input$asian_dateOrYear == "date"){
               maturity_reference <- as.numeric(difftime(input$maturityDate_asian, Sys.Date(), unit="weeks"))/52.25 
             } else if(input$asian_dateOrYear == "years"){
               maturity_reference <- input$timeToMaturity_asian
             }
             , silent = TRUE
             )
             
             try (
             # check for option specific input parameters
             if(maturity_reference != reactives[["result"]][[type]][["stockPaths"]][["maturity"]] ||
                input$optionType_asian != reactives[["result"]][[type]][["type"]][[3]] || 
                input$subType_asian != reactives[["result"]][[type]][["type"]][[2]] ||
                input$skalar_timeSteps != reactives[["result"]][[type]][["stockPaths"]][["skalarTimeSteps"]]){
               
               reactive_warning_dummy <- TRUE
               warning <- warning_dummy
             }
             , silent = TRUE
             )
             
             try (
             # check for different types of observation
             if(input$observation_asian != reactives[["result"]][[type]][["inputparams"]][["asianRange"]]){
               reactive_warning_dummy <- TRUE
               warning <- warning_dummy
             }else if(input$observation_asian == reactives[["result"]][[type]][["inputparams"]][["asianRange"]]){
               if(input$observation_asian == "Total"){
                 if(input$strikePrice_asian_total != reactives[["result"]][[type]][["strikePrice"]]){
                   reactive_warning_dummy <- TRUE
                   warning <- warning_dummy
                 }
               }else if(input$observation_asian == "Period"){
                 if(input$strikePrice_asian_period != reactives[["result"]][[type]][["strikePrice"]] || 
                    input$periodOfObs_asian != reactives[["result"]][[type]][["inputparams"]][["periodObs_asian"]]){
                   reactive_warning_dummy <- TRUE
                   warning <- warning_dummy
                 }
                               
               }else if(input$observation_asian == "Interval"){
                 if(input$strikePrice_asian_interval != reactives[["result"]][[type]][["strikePrice"]] || 
                    input$daysBetweenObs_asian != reactives[["result"]][[type]][["inputparams"]][["intervalObs_asian"]]){
                   reactive_warning_dummy <- TRUE
                   warning <- warning_dummy
                 }
               }
             }
             , silent = TRUE
             )
             
             
           },
           "barrier"={
             try (
             # save current maturity
             if(input$barrier_dateOrYear == "date"){
               maturity_reference <- as.numeric(difftime(input$maturityDate_barrier, Sys.Date(), unit="weeks"))/52.25 
             } else if(input$barrier_dateOrYear == "years"){
               maturity_reference <- input$timeToMaturity_barrier
             }
             , silent = TRUE
             )
             
             
             try(
             # check for option specific input parameters
             if(input$strikePrice_barrier != reactives[["result"]][[type]][["strikePrice"]] ||
                input$subType_barrier != reactives[["result"]][[type]][["type"]][[2]] ||
                input$barrierValue != reactives[["result"]][[type]][["inputparams"]][["barrier"]] ||
                input$timeSteps != reactives[["result"]][[type]][["stockPaths"]][["numberTimeSteps"]] ||
                maturity_reference != reactives[["result"]][[type]][["stockPaths"]][["maturity"]] ||
                input$optionType_barrier != reactives[["result"]][[type]][["type"]][[3]]){
               
               reactive_warning_dummy <- TRUE
               warning <- warning_dummy
             }
             , silent = TRUE
             )
             
           },
           "binary"={
             try ( 
             # save current maturity
             if(input$binary_dateOrYear == "date"){
               maturity_reference <- as.numeric(difftime(input$maturityDate_binary, Sys.Date(), unit="weeks"))/52.25 
             } else if(input$binary_dateOrYear == "years"){
               maturity_reference <- input$timeToMaturity_binary
             }
             , silent = TRUE
             )
             
             try(
             # check for option specific input parameters
             if(input$subType_binary != reactives[["result"]][[type]][["type"]][[2]] ||
                input$timeSteps != reactives[["result"]][[type]][["stockPaths"]][["numberTimeSteps"]] ||
                maturity_reference != reactives[["result"]][[type]][["stockPaths"]][["maturity"]] ||
                input$optionType_binary != reactives[["result"]][[type]][["type"]][[3]]){
               
               reactive_warning_dummy <- TRUE
               warning <- warning_dummy
             }
             , silent = TRUE
             )
             
             
             try(
             if(input$subType_binary == "Asset-or-nothing"){
               
               if(input$strikePrice_binary_asset != reactives[["result"]][[type]][["strikePrice"]]){
                 
                 reactive_warning_dummy <- TRUE
                 warning <- warning_dummy
                 
               } else{
                 
               }
               
               
             } else if(input$subType_binary == "Cash-or-nothing"){
               
               if(input$strikePrice_binary_cash != reactives[["result"]][[type]][["strikePrice"]]){
                 
                 reactive_warning_dummy <- TRUE
                 warning <- warning_dummy
                 
               } else{
                 
               }
               
             }
             , silent = TRUE
             )
             
             
             try(
             if(input$subType_binary == "Cash-or-nothing" && !is.null(reactives[["result"]][[type]])){
               
               if(input$binaryPayoff != reactives[["result"]][[type]][["inputparams"]][["binaryPayoff"]]){
                 
                 reactive_warning_dummy <- TRUE
                 warning <- warning_dummy
                 
               }
               
             }
             , silent = TRUE
             )
             
           },
           "chooser"={
             try(
             # save current maturity
             if(input$chooser_dateOrYear == "date"){
               maturity_reference <- as.numeric(difftime(input$maturityDate_chooser, Sys.Date(), unit="weeks"))/52.25 
             } else if(input$chooser_dateOrYear == "years"){
               maturity_reference <- input$timeToMaturity_chooser
             }
             , silent = TRUE
             )
             
             
             try(
             # check for option specific input parameters
             if(input$strikePrice_chooser != reactives[["result"]][[type]][["strikePrice"]] ||
                input$skalar_timeSteps != reactives[["result"]][[type]][["stockPaths"]][["skalarTimeSteps"]] ||
                maturity_reference != reactives[["result"]][[type]][["stockPaths"]][["maturity"]] ||
                input$chooserDecisionDate != reactives[["result"]][[type]][["inputparams"]][["choosingDate"]] ||
                input$chooser_option != reactives[["result"]][[type]][["inputparams"]][["chooserOption"]]){

               reactive_warning_dummy <- TRUE
               warning <- warning_dummy

             }
             , silent = TRUE
             )
             
             try(
             if(input$chooser_option == "asian"){
               
               if(input$chooser_subtype_asian != reactives[["result"]][[type]][["inputparams"]][["chooserSubType"]]){
                 
                 reactive_warning_dummy <- TRUE
                 warning <- warning_dummy
                 
               }
               
               if(!is.null(reactives[["result"]][[type]][["inputparams"]][["asianRange"]])){
                 
                 # check for different types of observation
                 if(input$observation_asian_chooser != reactives[["result"]][[type]][["inputparams"]][["asianRange"]]){
                   reactive_warning_dummy <- TRUE
                   warning <- warning_dummy
                 }
                 
               }
               
               if(!is.null(reactives[["result"]][[type]][["inputparams"]][["asianRange"]]) && !is.null(reactives[["result"]][[type]][["inputparams"]][["periodObs_asian"]])){
                 
                 if(input$observation_asian_chooser == "Period"){
                   if(input$periodOfObs_asian_chooser != reactives[["result"]][[type]][["inputparams"]][["periodObs_asian"]]){
                     reactive_warning_dummy <- TRUE
                     warning <- warning_dummy
                   }
                   
                 }
                 
               }
               
               if(!is.null(reactives[["result"]][[type]][["inputparams"]][["asianRange"]]) && !is.null(reactives[["result"]][[type]][["inputparams"]][["intervalObs_asian"]])){
                 
                 if(input$observation_asian_chooser == "Interval"){
                   if(input$daysBetweenObs_asian_chooser != reactives[["result"]][[type]][["inputparams"]][["intervalObs_asian"]]){
                     reactive_warning_dummy <- TRUE
                     warning <- warning_dummy
                   }
                   
                 }
                 
               }
               
               
             } else if(input$chooser_option == "barrier"){
               
               if(input$chooser_subtype_barrier != reactives[["result"]][[type]][["inputparams"]][["chooserSubType"]]){
                 
                 reactive_warning_dummy <- TRUE
                 warning <- warning_dummy
                 
               }
               
               if(!is.null(reactives[["result"]][[type]][["inputparams"]][["barrier"]])){
                 
                 # check for option specific input parameters
                 if(input$barrierValue_chooser != reactives[["result"]][[type]][["inputparams"]][["barrier"]]){
                   
                   reactive_warning_dummy <- TRUE
                   warning <- warning_dummy
                   
                 }
                 
               }
               
               
               
               
             } else if(input$chooser_option == "binary"){
               
               if(input$chooser_subtype_binary != reactives[["result"]][[type]][["inputparams"]][["chooserSubType"]]){
                 
                 reactive_warning_dummy <- TRUE
                 warning <- warning_dummy
                 
               }
               
               if(!is.null(reactives[["result"]][[type]][["inputparams"]][["binaryPayoff"]])){
                 
                 if(input$chooser_subtype_binary == "Cash-or-nothing" && !is.null(reactives[["result"]][[type]])){
                   
                   if(input$binaryPayoff_chooser != reactives[["result"]][[type]][["inputparams"]][["binaryPayoff"]]){
                     
                     reactive_warning_dummy <- TRUE
                     warning <- warning_dummy
                     
                   }
                   
                 }
                 
               }
               
               
               
               
             } else if(input$chooser_option == "european"){
               
               # no parameters to check
               
             } else if(input$chooser_option == "forwardStart"){
               
               if(input$chooser_subtype_forwardStart != reactives[["result"]][[type]][["inputparams"]][["chooserSubType"]]){
                 
                 reactive_warning_dummy <- TRUE
                 warning <- warning_dummy
                 
               }
               
               if(!is.null(reactives[["result"]][[type]][["inputparams"]][["startDate"]])){
                 
                 # check for option specific input parameters
                 if(input$forwardStartDate_chooser != reactives[["result"]][[type]][["inputparams"]][["startDate"]]){
                   
                   reactive_warning_dummy <- TRUE
                   warning <- warning_dummy
                   
                 } else {
                   
                 }
                 
               }
               
               
             } else if(input$chooser_option == "gap"){
               
               if(input$chooser_subtype_gap != reactives[["result"]][[type]][["inputparams"]][["chooserSubType"]]){
                 
                 reactive_warning_dummy <- TRUE
                 warning <- warning_dummy
                 
               }
               
               if(!is.null(reactives[["result"]][[type]][["inputparams"]][["strikePrice2"]])){
                 
                 # check for option specific input parameters
                 if(input$triggerPrice_gap_chooser != reactives[["result"]][[type]][["inputparams"]][["strikePrice2"]]){
                   
                   reactive_warning_dummy <- TRUE
                   warning <- warning_dummy
                   
                 } else {
                   
                 }
                 
               }
               
               
             } else if(input$chooser_option == "lookback"){
               
               # no parameters to check
               
             }
             , silent = TRUE
             )
             
           },
           "compound"={
             
             try(
               # save current maturity for the overlying option
               if(input$compound_first_dateOrYear == "date"){
                 maturity_reference_first <- as.numeric(difftime(input$maturityDate_compound_first, Sys.Date(), unit="weeks"))/52.25 
               } else if(input$compound_first_dateOrYear == "years"){
                 maturity_reference_first <- input$timeToMaturity_compound_first
               }
               , silent = TRUE
             )
             
             try(
               # save current maturity for the underlying option
               if(input$compound_second_dateOrYear == "date"){
                 
                 if(input$compound_first_dateOrYear == "date"){
                   
                   maturity_reference_second <- as.numeric(difftime(input$maturityDate_compound_second, as.Date(Sys.Date() + floor(reactives[["result"]][[type]][["stockPaths"]][["maturity"]] * 365)), unit="weeks"))/52.25
                   
                   
                 } else if(input$compound_first_dateOrYear == "years"){
                   
                   maturity_reference_second <- as.numeric(difftime(input$maturityDate_compound_second, as.Date(Sys.Date() + floor(365 * input$timeToMaturity_compound_first)), unit="weeks"))/52.25
                   
                 }
                 
                 
               } else if(input$compound_second_dateOrYear == "years"){
                 maturity_reference_second <- input$timeToMaturity_compound_second
               }
               , silent = TRUE
             )
             
             try(
               # check for option specific input parameters
               if(input$strikePrice_compound_first != reactives[["result"]][[type]][["strikePrice"]] ||
                  input$strikePrice_compound_second != reactives[["result"]][[type]][["inputparams"]][["secondStrikePrice"]] ||

                  input$skalar_timeSteps != reactives[["result"]][[type]][["inputparams"]][["secondSkalarTimeSteps"]] ||

                  maturity_reference_first != reactives[["result"]][[type]][["stockPaths"]][["maturity"]] ||
                  maturity_reference_second != reactives[["result"]][[type]][["inputparams"]][["secondExerciseDate"]] ||

                  input$sampleSize_compound_second != reactives[["result"]][[type]][["inputparams"]][["secondSampleSize"]] ||

                  input$optionType_compound != reactives[["result"]][[type]][["type"]][[2]] ||
                  input$compound_option_second != reactives[["result"]][[type]][["inputparams"]][["type2"]]){

                 reactive_warning_dummy <- TRUE
                 warning <- warning_dummy

               }
               , silent = TRUE
             )
             
             
             try(

               if(input$compound_option_second != "european"){

                 if(input$timeSteps != reactives[["result"]][[type]][["inputparams"]][["secondNumberTimeSteps"]]){

                   reactive_warning_dummy <- TRUE
                   warning <- warning_dummy

                 }

               }
               , silent = TRUE
             )
             
             
             try(
               if(input$compound_option_second == "asian"){
                 
                 if(input$compound_subtype_asian != reactives[["result"]][[type]][["inputparams"]][["secondSubType"]]){
                   
                   reactive_warning_dummy <- TRUE
                   warning <- warning_dummy
                   
                 }
                 
                 if(!is.null(reactives[["result"]][[type]][["inputparams"]][["asianRange"]])){
                   
                   # check for different types of observation
                   if(input$observation_asian_compound != reactives[["result"]][[type]][["inputparams"]][["asianRange"]]){
                     reactive_warning_dummy <- TRUE
                     warning <- warning_dummy
                   }
                   
                 }
                 
                 if(!is.null(reactives[["result"]][[type]][["inputparams"]][["asianRange"]]) && !is.null(reactives[["result"]][[type]][["inputparams"]][["periodObs_asian"]])){
                   
                   if(input$observation_asian_compound == "Period"){
                     if(input$periodOfObs_asian_compound != reactives[["result"]][[type]][["inputparams"]][["periodObs_asian"]]){
                       reactive_warning_dummy <- TRUE
                       warning <- warning_dummy
                     }
                     
                   }
                   
                 }
                 
                 if(!is.null(reactives[["result"]][[type]][["inputparams"]][["asianRange"]]) && !is.null(reactives[["result"]][[type]][["inputparams"]][["intervalObs_asian"]])){
                   
                   if(input$observation_asian_compound == "Interval"){
                     if(input$daysBetweenObs_asian_compound != reactives[["result"]][[type]][["inputparams"]][["intervalObs_asian"]]){
                       reactive_warning_dummy <- TRUE
                       warning <- warning_dummy
                     }
                     
                   }
                   
                 }
                 
                 
               } else if(input$compound_option_second == "barrier"){
                 
                 if(input$compound_subtype_barrier != reactives[["result"]][[type]][["inputparams"]][["secondSubType"]]){
                   
                   reactive_warning_dummy <- TRUE
                   warning <- warning_dummy
                   
                 }
                 
                 if(!is.null(reactives[["result"]][[type]][["inputparams"]][["barrier"]])){
                   
                   # check for option specific input parameters
                   if(input$barrierValue_compound != reactives[["result"]][[type]][["inputparams"]][["barrier"]]){
                     
                     reactive_warning_dummy <- TRUE
                     warning <- warning_dummy
                     
                   }
                   
                 }
                 
               } else if(input$compound_option_second == "binary"){
                 
                 if(input$compound_subtype_binary != reactives[["result"]][[type]][["inputparams"]][["secondSubType"]]){
                   
                   reactive_warning_dummy <- TRUE
                   warning <- warning_dummy
                   
                 }
                 
                 if(!is.null(reactives[["result"]][[type]][["inputparams"]][["binaryPayoff"]])){
                   
                   if(input$compound_subtype_binary == "Cash-or-nothing" && !is.null(reactives[["result"]][[type]])){
                     
                     if(input$binaryPayoff_compound != reactives[["result"]][[type]][["inputparams"]][["binaryPayoff"]]){
                       
                       reactive_warning_dummy <- TRUE
                       warning <- warning_dummy
                       
                     }
                     
                   }
                   
                 }
                 
                 
               } else if(input$compound_option_second == "european"){
                 
                 # no parameters to check
                 
               } else if(input$compound_option_second == "forwardStart"){
                 
                 if(!is.null(reactives[["result"]][[type]][["inputparams"]][["startDate"]])){
                   
                   # check for option specific input parameters
                   if(input$forwardStartDate_compound != reactives[["result"]][[type]][["inputparams"]][["startDate"]]){
                     
                     reactive_warning_dummy <- TRUE
                     warning <- warning_dummy
                     
                   } else {
                     
                   }
                   
                 }
                 
                 
               } else if(input$compound_option_second == "gap"){
                 
                 if(!is.null(reactives[["result"]][[type]][["inputparams"]][["strikePrice2"]])){
                   
                   # check for option specific input parameters
                   if(input$triggerPrice_gap_compound != reactives[["result"]][[type]][["inputparams"]][["strikePrice2"]]){
                     
                     reactive_warning_dummy <- TRUE
                     warning <- warning_dummy
                     
                   } else {
                     
                   }
                   
                 }
                 
                 
               } else if(input$compound_option_second == "lookback"){
                 
                 # no parameters to check
                 
               }
               , silent = TRUE
             )
             
             
           },
           "european"={
             
             try(
               
               # save current maturity
               if(input$european_dateOrYear == "date"){
                 maturity_reference <- as.numeric(difftime(input$maturityDate_european, Sys.Date(), unit="weeks"))/52.25
               } else if(input$european_dateOrYear == "years"){
                 maturity_reference <- input$timeToMaturity_european
               }
               
               , silent = TRUE
             )
             
             try(
               
               # check for option specific input parameters
               if(input$strikePrice_european != reactives[["result"]][[type]][["strikePrice"]] ||
                  maturity_reference != reactives[["result"]][[type]][["stockPaths"]][["maturity"]] ||
                  input$optionType_european != reactives[["result"]][[type]][["type"]][[3]]){
                 
                 reactive_warning_dummy <- TRUE
                 warning <- warning_dummy
                 
               } else {
                 
               }
               , silent = TRUE
             )
             
           },
           "forwardStart"={
             try(
             # save current maturity
             if(input$forwardStart_dateOrYear == "date"){
               maturity_reference <- as.numeric(difftime(input$maturityDate_forwardStart, Sys.Date(), unit="weeks"))/52.25 
             } else if(input$forwardStart_dateOrYear == "years"){
               maturity_reference <- input$timeToMaturity_forwardStart
             }
             , silent = TRUE
             )
             
             try(
             # check for option specific input parameters
             if(input$forwardStartDate != reactives[["result"]][[type]][["inputparams"]][["startDate"]] ||
                input$strikePrice_forwardStart != reactives[["result"]][[type]][["strikePrice"]] ||
                maturity_reference != reactives[["result"]][[type]][["stockPaths"]][["maturity"]] ||
                input$optionType_forwardStart != reactives[["result"]][[type]][["type"]][[3]] ||
                input$timeSteps != reactives[["result"]][[type]][["stockPaths"]][["numberTimeSteps"]]
             ){
               
               reactive_warning_dummy <- TRUE
               warning <- warning_dummy
               
             } else {
               
             }
             , silent = TRUE
             )
             
           },
           "gap"={
             try(
             # save current maturity
             if(input$gap_dateOrYear == "date"){
               maturity_reference <- as.numeric(difftime(input$maturityDate_gap, Sys.Date(), unit="weeks"))/52.25 
             } else if(input$gap_dateOrYear == "years"){
               maturity_reference <- input$timeToMaturity_gap
             }
             , silent = TRUE
             )
             
             try(
             # check for option specific input parameters
             if(input$strikePrice_gap1 != reactives[["result"]][[type]][["inputparams"]][["strikePrice1"]] ||
                input$strikePrice_gap2 != reactives[["result"]][[type]][["inputparams"]][["strikePrice2"]] ||
                maturity_reference != reactives[["result"]][[type]][["stockPaths"]][["maturity"]] ||
                input$optionType_gap != reactives[["result"]][[type]][["type"]][[3]] ||
                input$timeSteps != reactives[["result"]][[type]][["stockPaths"]][["numberTimeSteps"]]
                ){
               
               reactive_warning_dummy <- TRUE
               warning <- warning_dummy
               
             } else {
               
             }
             , silent = TRUE
             )
             
           },
           "lookback"={
             try(
             # save current maturity
             if(input$lookback_dateOrYear == "date"){
               maturity_reference <- as.numeric(difftime(input$maturityDate_lookback, Sys.Date(), unit="weeks"))/52.25 
             } else if(input$lookback_dateOrYear == "years"){
               maturity_reference <- input$timeToMaturity_lookback
             }
             , silent = TRUE
             )
             
             try(
             # check for option specific input parameters
             if(input$strikePrice_lookback != reactives[["result"]][[type]][["strikePrice"]] ||
                input$skalar_timeSteps != reactives[["result"]][[type]][["stockPaths"]][["skalarTimeSteps"]] ||
                input$subType_lookback != reactives[["result"]][[type]][["type"]][[2]] ||
                maturity_reference != reactives[["result"]][[type]][["stockPaths"]][["maturity"]] ||
                input$optionType_lookback != reactives[["result"]][[type]][["type"]][[3]]){
               
               reactive_warning_dummy <- TRUE
               warning <- warning_dummy
               
             }
             , silent = TRUE
             )
             
           })
    
    try(
    # check for dividend
    if(!is.null(reactives[["result"]][[type]][["stockPaths"]][["dividendRate"]])){

      if(input$dividend_YesOrNo == TRUE){

        if(input$dividendRate/100 != reactives[["result"]][[type]][["stockPaths"]][["dividendRate"]]){

          reactive_warning_dummy <- TRUE
          warning <- warning_dummy

        } else {

        }

      } else if(input$dividend_YesOrNo == FALSE){

        reactive_warning_dummy <- TRUE
        warning <- warning_dummy

      }

    } else if(is.null(reactives[["result"]][[type]][["stockPaths"]][["dividendRate"]])){

      if(input$dividend_YesOrNo == TRUE){

        reactive_warning_dummy <- TRUE
        warning <- warning_dummy

      } else if(input$dividend_YesOrNo == FALSE){

      }

    }
    , silent = TRUE
    )

    
    try(
    # check for simulationmodel type and the corresponding input parameters
    if(reactives$result[[type]][["stockPaths"]][["simulationModel"]] == "Black-Scholes"){

      if(input$stockPrice != reactives[["result"]][[type]][["stockPaths"]][["stockPrice"]]||
         input$volatility/100 != reactives[["result"]][[type]][["stockPaths"]][["volatility"]]||
         input$riskFreeRate/100 != reactives[["result"]][[type]][["stockPaths"]][["interestRate"]]){
        
        reactive_warning_dummy <- TRUE
        warning <- warning_dummy
        
        
      } else{
        
      }

    } else if(reactives$result[[type]][["stockPaths"]][["simulationModel"]] == "Heston"){

      if(input$stockPrice != reactives[["result"]][[type]][["stockPaths"]][["stockPrice"]]||
         input$volatility/100 != reactives[["result"]][[type]][["stockPaths"]][["volatility"]]||
         input$riskFreeRate/100 != reactives[["result"]][[type]][["stockPaths"]][["interestRate"]] ||
         input$kappa != reactives[["result"]][[type]][["stockPaths"]][["kappa"]] ||
         input$epsilon != reactives[["result"]][[type]][["stockPaths"]][["epsilon"]] ||
         input$rho != reactives[["result"]][[type]][["stockPaths"]][["rho"]] ||
         input$theta != reactives[["result"]][[type]][["stockPaths"]][["theta"]]){
        
        reactive_warning_dummy <- TRUE
        warning <- warning_dummy
        
      }
      
    } else{
      
    }
    , silent = TRUE
    )
    
    
    try(
      
      if(input$simulationModel != reactives$result[[type]][["stockPaths"]][["simulationModel"]]){
        
        reactive_warning_dummy <- TRUE
        warning <- warning_dummy
        
      } else {
        
      }
      , silent = TRUE
    )
    
    
    try(

      if(type == "compound"){

        if(input$sampleSize_compound != reactives[["result"]][[type]][["stockPaths"]][["sampleSize"]]){

          reactive_warning_dummy <- TRUE
          warning <- warning_dummy

        }

      } else {

        if(input$sampleSize*1000 != reactives[["result"]][[type]][["stockPaths"]][["sampleSize"]]){

          reactive_warning_dummy <- TRUE
          warning <- warning_dummy

        }
      }
      , silent = TRUE
    )
    
    reactives[["warning"]][[type]] <- reactive_warning_dummy
    
  }
  
  return(warning)
   
}


# function to detect if maturity as date or in years
choose.maturity.Input <- function(dateOrYear, maturityDate, timeToMaturity){
  
  if(dateOrYear == "date"){
    
    maturity <- maturityDate
    
  } else if(dateOrYear == "years"){
    
    maturity <- timeToMaturity
    
  }
  
  return(maturity)
  
}

## function to set default for risk free rate using 3 month US T bill
get.riskFreeRate <- function(){
  getSymbols('DTB3', src = 'FRED')
  mat <- data.matrix(as.data.frame(last(DTB3, '7 days')))
  r <- c()
  
  for(i in 7:1){
    if(is.na(mat[i,1])){
      
    }else{
      r <- mat[i,1]
      break
    }
  }
  return(as.numeric(r))
}




## create the progressboxes
# price
create.Progressbox_price <- function(reactives, type, boxtype, color, icon){
  
  box <- valueBox(
    
    if(!is.null(reactives$result[[type]])){
      
      paste0(round(reactives$result[[type]][[boxtype]][1], digits = 3))
      
    } else {
      
      paste0(0)
      
    }, "Price", icon = icon(icon, lib = "glyphicon"),
    color = color
  )
  
  
  return(box)
  
}

# sd
create.Progressbox_SD <- function(reactives, type, boxtype, color, icon){
  
  box <- valueBox(
    
    if(!is.null(reactives$result[[type]])){
      
      paste0(round(reactives$result[[type]][[boxtype]][1], digits = 3))
      
    } else {
      
      paste0(0)
      
    }, "Standard deviation", icon = icon(icon, lib = "glyphicon"),
    color = color
  )
  
  
  return(box)
  
}


# se
create.Progressbox_SE <- function(reactives, type, boxtype, color, icon){
  
  box <- valueBox(
    
    if(!is.null(reactives$result[[type]])){
      
      paste0(round(reactives$result[[type]][[boxtype]][1], digits = 3))
      
    } else {
      
      paste0(0)
      
    }, "Standard error", icon = icon(icon, lib = "glyphicon"),
    color = color
  )
  
  
  return(box)
  
}

# two conditional panels for maturity
create.Conditional.Maturity <- function(type){
  dateString <- paste("input.",type,"_dateOrYear == 'date'", sep = "")
  yearString <- paste("input.",type,"_dateOrYear == 'years'", sep = "")
  
  dateId <- paste("maturityDate_",type, sep = "")
  yearId <- paste("timeToMaturity_",type, sep = "")
  
  panelList <- list(
    conditionalPanel(
      condition = dateString,
      
      dateInput(inputId = dateId,
                label = "",
                value = (Sys.Date() + 7),
                daysofweekdisabled = c(0,6),
                datesdisabled = create.Holidays())
    ),
    conditionalPanel(
      condition = yearString,
      
      numericInput(inputId = yearId,
                   label = "",
                   value = 0.5)
    )
  )
  
  return(panelList)
  
}

# generate plot (histogram)
create.Histogram <- function(reactives, type){
  
  # default empty histogram
  histogram <- ggplot(data = data.frame(c(0:1))) + theme_bw() + annotate(geom = "text", x = 0.5, y = 0.1, label = "Option price not calculated yet.") + 
    labs(y = "frequency", 
         title = "",
         subtitle = "")
  
  # set the header for the plot
  plot_string <- ""
  
  if(!is.null(reactives$result[[type]])){
    
    if(!is.null(reactives[["result"]][[type]][["type"]][[3]])){
      
      if(reactives[["result"]][[type]][["type"]][[3]] == "Call"){
        
        call_or_put <- "call"
        
      } else if(reactives[["result"]][[type]][["type"]][[3]] == "Put"){
        
        call_or_put <- "put"
        
      }
      
    }
    
    switch(type,
           "american"={
             plot_string <- paste("American", call_or_put, "option")
           },
           "asian"={
             plot_string <- paste("Asian", call_or_put, "option")
           },
           "barrier"={
             plot_string <- paste("Barrier", call_or_put, "option")
           },
           "binary"={
             plot_string <- paste("Binary", call_or_put, "option")
           },
           "chooser"={
             plot_string <- paste("Chooser option -", reactives[["result"]][[type]][["inputparams"]][["chooserOption"]])
           },
           "compound"={
             plot_string <- paste("Compound option -", reactives[["result"]][[type]][["inputparams"]][["type2"]])
           },
           "european"={
             plot_string <- paste("European", call_or_put, "option")
           },
           "forwardStart"={
             plot_string <- paste("Forward Start", call_or_put, "option")
           },
           "gap"={
             plot_string <- paste("Gap", call_or_put, "option")
           },
           "lookback"={
             plot_string <- paste("Lookback", call_or_put, "option")
           }
    )
    
    # if the x axis range is to big scale it by log10
    if((max(unlist(reactives$payoff_data[[type]]["payoff"])) - quantile(unlist(reactives$payoff_data[[type]]["payoff"]), probs = 0.9)[[1]]) > 300){
      
      histogram <- ggplot(data = reactives$payoff_data[[type]]["payoff"], aes(x=payoff))+ geom_histogram(binwidth = 0.05, 
                                                                                                         colour = "darkolivegreen", 
                                                                                                         fill= "lightgreen", 
                                                                                                         linetype = 1) + 
        theme_light() + labs(title = plot_string,
                             subtitle = paste("Price =  ", round(reactives$result[[type]][["optionValue"]][1], digits = 3))) + 
        scale_x_log10(name="payoff", labels = scales::label_comma(scale = 1/1)) +
        scale_y_log10(name="frequency", labels = scales::label_comma(scale = 1/1)) 
      
    } else {
      
      histogram <- ggplot(data = reactives$payoff_data[[type]]["payoff"], aes(x=payoff))+ geom_histogram(binwidth = .3, 
                                                                                                         colour = "darkolivegreen", 
                                                                                                         fill= "lightgreen", 
                                                                                                         linetype = 1) + 
        theme_light() + labs(title = plot_string,
                             subtitle = paste("Price =  ", round(reactives$result[[type]][["optionValue"]][1], digits = 3))) + 
        scale_x_continuous(name="payoff", labels = scales::label_comma(scale = 1/1)) +
        scale_y_log10(name="frequency", labels = scales::label_comma(scale = 1/1)) 
      
    } 
    
    
  } else {
    
  }
  
  
  return(histogram)
}