## load packages
library(ggplot2)
library(gridExtra)
library(tidyr) 
library(timeDate)
library(matrixStats)
library(MASS)


# Source functions ----
source("functions.R")

# ignore warnings from the ggplot an other warnings
options(warn=-2)

##### Server ----------------------------------------------------------------------

# Server definieren
server <- function(input, output, session){
  
  ## store all reactive elements in the reactives list
  # those variables are visible for all wrapped reactive environments like "observeEvent"
  reactives <- reactiveValues(
    
    # "price.Option"  results
    result = list(),
    
    # data for the histogram for each option
    payoff_data = list(),
    
    # list that stores all options and only one of them can be true (last calculated option)
    # functionality implemented by the function "set.active.option"
    check_active_option = list(american = FALSE, asian = FALSE,
                               barrier = FALSE, binary = FALSE,
                               chooser = FALSE,
                               compound = FALSE, european = FALSE,
                               forwardStart = FALSE, gap = FALSE,
                               lookback = FALSE, shout = FALSE),
    
    # list for the warning, if option calculation is not up to date
    warning = list(american = FALSE, asian = FALSE,
                   barrier = FALSE, binary = FALSE,
                   chooser = FALSE,
                   compound = FALSE, european = FALSE,
                   forwardStart = FALSE, gap = FALSE,
                   lookback = FALSE, shout = FALSE)
    
  )
  
  ## events and ouput independent of specific options ----------------------
  
  ## conditional panel for the simulation model: if heston is chosen an additional tab is shown
  output$heston_active <- renderUI({
    
    if(input$simulationModel == "Black-Scholes"){
      
      sidebar <- sidebarMenu()
      
    } else if(input$simulationModel == "Heston"){
      
      sidebar <- sidebarMenu(
        
        shinyFeedback::useShinyFeedback(),
        
        waiter::use_waitress(),
        
        menuItem("Adjust Wiener processes", tabName="hestonTab"),
        
        id = "tab2"
      )
      
    }
    
    return(sidebar)
    
  })
  
  
  ## conditional panel for the dividend
  output$yield_chosen <- renderUI({
    
    if(input$dividend_YesOrNo == TRUE){
      
      sidebar <- sidebarMenu(
        
        shinyFeedback::useShinyFeedback(),
        waiter::use_waitress(),
        
        # dividend yield
        numericInput(inputId = "dividendRate",
                     label ="Dividend yield in percent",
                     value = 1.1,
                     min = 0,
                     max = 100
        ),
        id = "tab3"
      )
      
    } else {
      sidebar <- sidebarMenu()
    }
    
    return(sidebar)
    
  })
  
  
  ## conditional panel for the number of time steps: for options where the path has not to 
  # be an multiple of the number of days an integer greater than zero can be handed over
  output$timesteps <- renderUI({
    
    if(input$tab == "about" ||
       input$tab == "american" ||
       input$tab == "barrier" ||
       input$tab == "shout" ||
       input$tab == "forward" ||
       input$tab == "binary" ||
       input$tab == "gap" ||
       (input$tab == "compound" && input$compound_option_second != "asian" 
       && input$simulationModel == "Heston" && input$compound_option_second != "lookback") ||
       (input$tab == "compound" && input$compound_option_second != "asian" && input$compound_option_second != "lookback" 
        && input$simulationModel == "Black-Scholes" && input$compound_option_second != "european")){
      
      sidebar <- sidebarMenu(
        
        numericInput(inputId = "timeSteps", 
                     label = "Time steps",
                     value = 100, min = 1, max = 1500),
        id = "tab4"
      )
     }else {
        sidebar <- sidebarMenu()
     }
    return(sidebar)
    
  })
  
  
  ## conditionalpanel for the number of timesteps if Asian, Chooser, Lookback:
  # the timesteps are an multiple of the days of the maturity time days
  output$skalar_timesteps <- renderUI({
    if (input$tab == "asian" || input$tab == "chooser" || input$tab == "lookback" ||
        (input$tab == "compound" && input$compound_option_second == "asian") ||
        (input$tab == "compound" && input$compound_option_second == "lookback")){
      sidebar <- sidebarMenu(
        
        numericInput(inputId = "skalar_timeSteps", 
                     label = "Daily time steps",
                     value = 2, min = 1, max = 14),
        id = "tab5"
      )
    }else {
      sidebar <- sidebarMenu()
    }
  })
  
  
  
  ## heston tab: shown if heston is chosen
  # inital path plot after app has been started
  output$path_plot <- renderPlot({
    
    ggplot(data = data.frame(c(0:1))) + theme_bw() + annotate(geom = "text", x = 0.5, y = 0.1, label = "Stock path not calculated yet.") + 
      labs(y = "h(x)", 
           title = paste(""))
    
  },
  height = "auto",width = "auto")
  
  
  # plot of the Heston Path after button is pressed
  observeEvent(input$hestonPlot, {
    
    # check if there are missing or false values required for the path plot
    checkHeston_na <- check.Heston_na(reactives, input)
    checkHeston <- check.Heston(reactives, input)
    
    if(checkHeston$check && checkHeston_na$check){
      
      start_date <- Sys.Date()
      end_date <- input$hestonPathDate
      range <- seq(start_date, end_date,"days")
      allDates <- length(range)
      N <- length(range)-1
      
      # create an black-scholes as well as an heston path for the plot 
      if(input$dividend_YesOrNo == TRUE){
        
        path_Black <- create.StockPaths(simulationModel = "Black-Scholes",
                                        initialPrice = input$stockPrice, maturity = end_date,
                                        volatility = input$volatility/100, interestRate = input$riskFreeRate/100,
                                        sampleSize = 1, numberTimeSteps = N,
                                        dividendRate = input$dividendRate/100)
        
        path_Heston <- create.StockPaths(simulationModel = "Heston",
                                         initialPrice = input$stockPrice, maturity = end_date,
                                         volatility = input$volatility/100, interestRate = input$riskFreeRate/100,
                                         sampleSize = 1, numberTimeSteps = N,
                                         dividendRate = input$dividendRate/100,
                                         
                                         # heston-model specific input paramters
                                         kappa = input$kappa, epsilon = input$epsilon,
                                         rho = input$rho, theta = input$theta)
        
      } else if(input$dividend_YesOrNo == FALSE){
        
        path_Black <- create.StockPaths(simulationModel = "Black-Scholes",
                                        initialPrice = input$stockPrice, maturity = end_date,
                                        volatility = input$volatility/100, interestRate = input$riskFreeRate/100,
                                        dividendRate = 0, sampleSize = 1, numberTimeSteps = N)
        
        
        path_Heston <- create.StockPaths(simulationModel = "Heston",
                                         initialPrice = input$stockPrice, maturity = end_date,
                                         volatility = input$volatility/100, interestRate = input$riskFreeRate/100,
                                         dividendRate = 0, sampleSize = 1, numberTimeSteps = N,
                                         
                                         # heston-model specific input paramters
                                         kappa = input$kappa, epsilon = input$epsilon,
                                         rho = input$rho, theta = input$theta)
        
      }
      
      # the initial value of the asset shold be displayed as well
      hestonPath <- c(input$stockPrice, path_Heston[["St"]])
      blackPath <- c(input$stockPrice, path_Black[["St"]])
      
      # create the data frame
      plot_data <- data.frame(Heston=hestonPath, BlackScholes = blackPath, dates = range)
      
      stock_prices.tidy <- gather(plot_data, Symbol, Prices, -dates)
      colnames(stock_prices.tidy) <- c("Date", "Model", "Prices")
      
      # update the path plot
      output$path_plot <- renderPlot({
        
        ggplot(stock_prices.tidy, aes(x = Date, y = Prices, color = Model)) +
          geom_line(size = .9) + theme_light() + labs(title = "")
        
      },
      height = "auto",width = "auto")
      
    } else if(checkHeston$check == FALSE) { # if false parameters then show warning
      
      showNotification(checkHeston$warning, type = "error", duration = 3)
      
    } else if(checkHeston_na$check == FALSE){ # if missing parameters then show warning
      
      showNotification(checkHeston_na$warning, type = "error", duration = 3)
      
    }
    
    
  })
  
  # text for the description of the heston path plot
  output$hestonplotText <- renderText({
    paste('Choose the date of maturity and press the button in the box "Set parameters of the heston pricepath plot" to generate
           the plot of the pricepaths. It shows one Heston and one Black-Scholes path
           depending on the current input parameters in the sidebar- and heston input parameter panel. Only the sample size is not used.')
  })
  
  
  ### NOTICE: in the following, the option implementation is pretty mutch the same
  # for all options. Because of that the american option is commented detailed and all following options are
  # only commented if something special for the option occurs
  
  
  ## american ------------------------------------------------------------------
  # event that is triggered if the mc start button is pressed
  observeEvent(input$americanSimulation, {
    
    # if not fully loaded no calculation is made
    req(input$maturityDate_american)
    
    # set the active option to american: this is used in the "check.active.option
    set.active.option(reactives, "american")
    
    # check if the option specific parameters are correct
    checkInput <- check.Input(reactives, input)
    # check if the option specific parameters are not NA
    checkInput_na <- check.Input_na(reactives, input) 
    
    # initalize the warning variable: means that the option calculation is up to date
    reactives[["warning"]][["american"]] <- FALSE
    
    # if the input params are correct and not NA start the calculation
    if(checkInput$check && checkInput_na$check){
      
      # funcitn detects if the selected maturity is "date" or "years"
      maturity <- choose.maturity.Input(input$american_dateOrYear, 
                                        input$maturityDate_american, 
                                        input$timeToMaturity_american)
      
      
        # create the stockpaths in dependecy if the model contains an dividend yield or not
        if(input$dividend_YesOrNo == TRUE){
          
          if(input$simulationModel == "Black-Scholes"){
            withProgress(message = 'Generating stock paths',{
            # calling the path and option function
            stockPaths <- create.StockPaths(simulationModel= input$simulationModel,
                                            initialPrice = input$stockPrice,     maturity = maturity,
                                            volatility = input$volatility/100,   interestRate = input$riskFreeRate/100,
                                            dividendRate = input$dividendRate/100,
                                            sampleSize = input$sampleSize*1000,  numberTimeSteps = input$timeSteps)
            })
          } else if(input$simulationModel == "Heston"){
            withProgress(message = 'Generating stock paths',{
            # calling the path and option function
            stockPaths <- create.StockPaths(simulationModel= input$simulationModel,
                                            initialPrice = input$stockPrice,     maturity = maturity,
                                            volatility = input$volatility/100,   interestRate = input$riskFreeRate/100,
                                            dividendRate = input$dividendRate/100,
                                            sampleSize = input$sampleSize*1000,  numberTimeSteps = input$timeSteps,
                                            kappa = input$kappa, epsilon = input$epsilon,
                                            rho = input$rho, theta = input$theta)
            })
          }
          
        } else if(input$dividend_YesOrNo == FALSE){
          
          if(input$simulationModel == "Black-Scholes"){
            withProgress(message = 'Generating stock paths',{
            # calling the path and option function
            stockPaths <- create.StockPaths(simulationModel= input$simulationModel,
                                            initialPrice = input$stockPrice,     maturity = maturity,
                                            volatility = input$volatility/100,   interestRate = input$riskFreeRate/100,
                                            sampleSize = input$sampleSize*1000,  numberTimeSteps = input$timeSteps)
            })
          } else if(input$simulationModel == "Heston"){
            withProgress(message = 'Generating stock paths',{
            # calling the path and option function
            stockPaths <- create.StockPaths(simulationModel= input$simulationModel,
                                            initialPrice = input$stockPrice,     maturity = maturity,
                                            volatility = input$volatility/100,   interestRate = input$riskFreeRate/100,
                                            sampleSize = input$sampleSize*1000,  numberTimeSteps = input$timeSteps,
                                            kappa = input$kappa, epsilon = input$epsilon,
                                            rho = input$rho, theta = input$theta)
            })
          }
          
        }
      
      # create an progressbar that indicates the calculation progress
      withProgress(message = 'Calculating option value',{
        
        # option price is calculated
        result <- price.Option(type = "american",
                               optionType = input$optionType_american,
                               strikePrice = input$strikePrice_american,
                               stockPaths = c(stockPaths))
      
      # store the result in the reactives list
      reactives$result[["american"]] <- result
      # create data frame for the histogram plot
      reactives$payoff_data[["american"]] <- data.frame(payoff = result[["C0_plot"]])
      })
      
      # create a warning if the option calculation is not up to date
      output$warning_american <- renderUI({
        
        # because of the reactive behavior the function is checked each time a input or reactive variable changes
        create.warning(reactives, input, "american")
        
      })
      
      
      # if one or more input params are not correct show a notification and dont calculate
    } else if(checkInput$check == FALSE){
      
      showNotification(checkInput$warning, type = "error", duration = 3)
      
      reactives$result[["american"]] <- NULL
      reactives$payoff_data[["american"]] <- NULL
      
      # if one or more input params are NA show a notification and dont calculate
    } else if(checkInput_na$check == FALSE){
      
      showNotification(checkInput_na$warning, type = "error", duration = 3)
      
      reactives$result[["american"]] <- NULL
      reactives$payoff_data[["american"]] <- NULL
      
    }
    
  })
  
  
  ## progressboxes
  # american price progress box
  output$progressBox_american_price <- renderValueBox({
    
    # creates a progressbox for the result
    create.Progressbox_price(reactives= reactives, "american", "optionValue",
                             "blue", "euro")
    
  })
  
  # american option sd progress box
  output$progressBox_american_sd <- renderValueBox({
    
    # creates a progressbox for the standard deviation
    create.Progressbox_SD(reactives= reactives, "american", "optionSD",
                          "light-blue", "resize-horizontal")
    
  })
  
  # american option se progress box
  output$progressBox_american_se <- renderValueBox({
    
    # creates a progressbox for the standard error
    create.Progressbox_SE(reactives= reactives, "american", "optionSE",
                          "aqua", "exclamation-sign")
    
  })
  
  ## histogram of the discounted payoff
  observeEvent(reactives[["warning"]][["american"]], {
    
    # as long as the option calculation is up to date the histogram is shown
    # if not the plot is set to null so that the empty plot is shown
    if(reactives[["warning"]][["american"]] == FALSE){
      
      output$american_plot <- renderPlot({
        
        # creates an histogram of the discounted payoff with the package ggplot2
        create.Histogram(reactives = reactives, "american")
        
      },
      height = "auto",width = "auto")
      
    } else if(reactives[["warning"]][["american"]] == TRUE){
      
      output$american_plot <- NULL
      
    }
    
  })
  
  
  # conditional panel for maturity in date or years
  output$maturity_american <- renderUI({
    
    # function creates the panel with the maturity: either date or year
    create.Conditional.Maturity("american")
    
  })
  
  
  
  ## asian ---------------------------------------------------------------------
  observeEvent(input$asianSimulation, {
    
    req(input$maturityDate_asian)
    
    set.active.option(reactives, "asian")
    
    checkInput <- check.Input(reactives, input)
    checkInput_na <- check.Input_na(reactives, input) 
    
    reactives[["warning"]][["asian"]] <- FALSE
    
    if(checkInput$check && checkInput_na$check){
      
      maturity <- choose.maturity.Input(input$asian_dateOrYear, 
                                        input$maturityDate_asian, 
                                        input$timeToMaturity_asian)
      
      
      if(input$dividend_YesOrNo == TRUE){
        
        if(input$simulationModel == "Black-Scholes"){
          withProgress(message = 'Generating stock paths',{
          # calling the path and option function
          stockPaths <- create.StockPaths(simulationModel= input$simulationModel,
                                          initialPrice = input$stockPrice,     maturity = maturity,
                                          volatility = input$volatility/100,   interestRate = input$riskFreeRate/100,
                                          dividendRate = input$dividendRate/100,
                                          sampleSize = input$sampleSize*1000,  skalarTimeSteps = input$skalar_timeSteps)
          })
        } else if(input$simulationModel == "Heston"){
          withProgress(message = 'Generating stock paths',{
          # calling the path and option function
          stockPaths <- create.StockPaths(simulationModel= input$simulationModel,
                                          initialPrice = input$stockPrice,     maturity = maturity,
                                          volatility = input$volatility/100,   interestRate = input$riskFreeRate/100,
                                          dividendRate = input$dividendRate/100,
                                          sampleSize = input$sampleSize*1000,  skalarTimeSteps = input$skalar_timeSteps,
                                          kappa = input$kappa, epsilon = input$epsilon,
                                          rho = input$rho, theta = input$theta)
          })
        }
        
      } else if(input$dividend_YesOrNo == FALSE){
        
        if(input$simulationModel == "Black-Scholes"){
          withProgress(message = 'Generating stock paths',{
          # calling the path and option function
          stockPaths <- create.StockPaths(simulationModel= input$simulationModel,
                                          initialPrice = input$stockPrice,     maturity = maturity,
                                          volatility = input$volatility/100,   interestRate = input$riskFreeRate/100,
                                          sampleSize = input$sampleSize*1000,  skalarTimeSteps = input$skalar_timeSteps)
          })
        } else if(input$simulationModel == "Heston"){
          withProgress(message = 'Generating stock paths',{
          # calling the path and option function
          stockPaths <- create.StockPaths(simulationModel= input$simulationModel,
                                          initialPrice = input$stockPrice,     maturity = maturity,
                                          volatility = input$volatility/100,   interestRate = input$riskFreeRate/100,
                                          sampleSize = input$sampleSize*1000,  skalarTimeSteps = input$skalar_timeSteps,
                                          kappa = input$kappa, epsilon = input$epsilon,
                                          rho = input$rho, theta = input$theta)
          })
        }
        
      }
      
      withProgress(message = 'Calculating option value', {
      # seperation via if regarded because of check if caclulation is not up to date
      # and because of three different strike prices (otherwise reactive does not work)
      if(input$observation_asian == "Total"){
        
        result <- price.Option(type = "asian",
                               optionType = input$optionType_asian,
                               subType = input$subType_asian,
                               strikePrice = input$strikePrice_asian_total,
                               asianRange = "Total",
                               stockPaths = c(stockPaths))
        
      } else if (input$observation_asian == "Interval"){
        
        result <- price.Option(type = "asian",
                               optionType = input$optionType_asian,
                               subType = input$subType_asian,
                               strikePrice = input$strikePrice_asian_interval,
                               asianRange = "Interval",
                               intervalObs_asian = input$daysBetweenObs_asian,
                               stockPaths = c(stockPaths))
        
      } else if(input$observation_asian == "Period"){
        
        result <- price.Option(type = "asian",
                               optionType = input$optionType_asian,
                               subType = input$subType_asian,
                               strikePrice = input$strikePrice_asian_period,
                               asianRange = "Period",
                               periodObs_asian = input$periodOfObs_asian,
                               stockPaths = c(stockPaths))
        
      }
      
      reactives$result[["asian"]] <- result
      reactives$payoff_data[["asian"]] <- data.frame(payoff = result[["C0_plot"]])
      
      })
      
      output$warning_asian <- renderUI({
        
        create.warning(reactives, input, "asian")
        
      })
      
      
    } else if(checkInput$check == FALSE){
      
      showNotification(checkInput$warning, type = "error", duration = 3)
      
      reactives$result[["asian"]] <- NULL
      reactives$payoff_data[["asian"]] <- NULL
      
    }else if(checkInput_na$check == FALSE){
      
      showNotification(checkInput_na$warning, type = "error", duration = 3)
      
      reactives$result[["asian"]] <- NULL
      reactives$payoff_data[["asian"]] <- NULL
      
    }
    
  })
  
  
  # progressboxes
  # asian price progress box
  output$progressBox_asian_price <- renderValueBox({
    
    create.Progressbox_price(reactives= reactives, "asian", "optionValue",
                             "blue", "euro")
    
  })
  
  # asian option sd progress box
  output$progressBox_asian_sd <- renderValueBox({
    
    create.Progressbox_SD(reactives= reactives, "asian", "optionSD",
                          "light-blue", "resize-horizontal")
    
  })
  
  # asian option se progress box
  output$progressBox_asian_se <- renderValueBox({
    
    create.Progressbox_SE(reactives= reactives, "asian", "optionSE",
                          "aqua", "exclamation-sign")
    
  })
  
  
  ## histogram of the discounted payoff
  observeEvent(reactives[["warning"]][["asian"]], {
    
    if(reactives[["warning"]][["asian"]] == FALSE){
      
      output$asian_plot <- renderPlot({
        
        create.Histogram(reactives = reactives, "asian")
        
      },
      height = "auto",width = "auto")
      
    } else if(reactives[["warning"]][["asian"]] == TRUE){
      
      output$asian_plot <- NULL
      
    }
    
  })
  
  
  ## conditional panel for maturity in date or years
  output$maturity_asian <- renderUI({
    
    create.Conditional.Maturity("asian")
    
  })
  
  ## conditional panel for total range or partial range
  output$date_range_asian <- renderUI({
    
    list(
      conditionalPanel(
        condition = "input.observation_asian == 'Total'",
        
        numericInput(inputId = "strikePrice_asian_total",
                     label = "Strike price",
                     value = 50,
                     width = "99%")
      ),
      conditionalPanel(
        condition = "input.observation_asian == 'Interval'",
        
        splitLayout(
          cellWidths = c("30%","70%"),
          
          numericInput(inputId = "strikePrice_asian_interval",
                       label = "Strike price",
                       value = 50,
                       width = "99%"),
          
          numericInput(inputId = "daysBetweenObs_asian",
                       label = tags$span(
                         "Days between observations ", 
                         tags$i(
                           class = "glyphicon glyphicon-info-sign", 
                           style = "color:#0072B2;",
                           title = 'Interval between observations in days (incl. weekends).' 
                         )
                       ),
                       value = 1,
                       min = 1,
                       width = "95%")
          
        )
        
      ),
      conditionalPanel(
        condition = "input.observation_asian == 'Period'",
        
        splitLayout(
          cellWidths = c("30%","70%"),
          
          numericInput(inputId = "strikePrice_asian_period",
                       label = "Strike price",
                       value = 50,
                       width = "99%"),
          
          numericInput(inputId = "periodOfObs_asian",
                       label = tags$span(
                         "Days before expiry ", 
                         tags$i(
                           class = "glyphicon glyphicon-info-sign", 
                           style = "color:#0072B2;",
                           title = 'Number of days from the period of observation.' 
                         )
                       ),
                       value = 90,
                       min = 2,
                       width = "95%")
          
        )
        
        
      )
    )
  })
  
  
  ## barrier -------------------------------------------------------------------
  observeEvent(input$barrierSimulation, {
    
    # if not fully loaded no calculation is made
    req(input$maturityDate_barrier)
    
    set.active.option(reactives, "barrier")
    
    checkInput <- check.Input(reactives, input)
    checkInput_na <- check.Input_na(reactives, input) 
    
    reactives[["warning"]][["barrier"]] <- FALSE
    
    if(checkInput$check && checkInput_na$check){
      
      maturity <- choose.maturity.Input(input$barrier_dateOrYear, 
                                        input$maturityDate_barrier, 
                                        input$timeToMaturity_barrier)
      
      if(input$dividend_YesOrNo == TRUE){
        
        if(input$simulationModel == "Black-Scholes"){
          withProgress(message = 'Generating stock paths',{
          # calling the path and option function
          stockPaths <- create.StockPaths(simulationModel= input$simulationModel,
                                          initialPrice = input$stockPrice,     maturity = maturity,
                                          volatility = input$volatility/100,   interestRate = input$riskFreeRate/100,
                                          dividendRate = input$dividendRate/100,
                                          sampleSize = input$sampleSize*1000,  numberTimeSteps = input$timeSteps)
          })
        } else if(input$simulationModel == "Heston"){
          withProgress(message = 'Generating stock paths',{
          # calling the path and option function
          stockPaths <- create.StockPaths(simulationModel= input$simulationModel,
                                          initialPrice = input$stockPrice,     maturity = maturity,
                                          volatility = input$volatility/100,   interestRate = input$riskFreeRate/100,
                                          dividendRate = input$dividendRate/100,
                                          sampleSize = input$sampleSize*1000,  numberTimeSteps = input$timeSteps,
                                          kappa = input$kappa, epsilon = input$epsilon,
                                          rho = input$rho, theta = input$theta)
          })
        }
        
      } else if(input$dividend_YesOrNo == FALSE){
        
        if(input$simulationModel == "Black-Scholes"){
          withProgress(message = 'Generating stock paths',{
          # calling the path and option function
          stockPaths <- create.StockPaths(simulationModel= input$simulationModel,
                                          initialPrice = input$stockPrice,     maturity = maturity,
                                          volatility = input$volatility/100,   interestRate = input$riskFreeRate/100,
                                          sampleSize = input$sampleSize*1000,  numberTimeSteps = input$timeSteps)
          })
        } else if(input$simulationModel == "Heston"){
          withProgress(message = 'Generating stock paths',{
          # calling the path and option function
          stockPaths <- create.StockPaths(simulationModel= input$simulationModel,
                                          initialPrice = input$stockPrice,     maturity = maturity,
                                          volatility = input$volatility/100,   interestRate = input$riskFreeRate/100,
                                          sampleSize = input$sampleSize*1000,  numberTimeSteps = input$timeSteps,
                                          kappa = input$kappa, epsilon = input$epsilon,
                                          rho = input$rho, theta = input$theta)
          })
        }
        
      }
      
      withProgress(message = 'Calculating option value', {
      result <- price.Option(type = "barrier",
                             optionType = input$optionType_barrier,
                             subType = input$subType_barrier,
                             barrier = input$barrierValue,
                             strikePrice = input$strikePrice_barrier,
                             stockPaths = c(stockPaths))
      
      reactives$result[["barrier"]] <- result
      reactives$payoff_data[["barrier"]] <- data.frame(payoff = result[["C0_plot"]])
      })
      
      output$warning_barrier <- renderUI({
        
        create.warning(reactives, input, "barrier")
        
      })
      
      
    } else if(checkInput$check == FALSE){
      
      showNotification(checkInput$warning, type = "error", duration = 3)
      
      reactives$result[["barrier"]] <- NULL
      reactives$payoff_data[["barrier"]] <- NULL
      
    }else if(checkInput_na$check == FALSE){
      
      showNotification(checkInput_na$warning, type = "error", duration = 3)
      
      reactives$result[["barrier"]] <- NULL
      reactives$payoff_data[["barrier"]] <- NULL
      
    }
    
  })
  
  
  ## progressboxes
  # barrier price progress box
  output$progressBox_barrier_price <- renderValueBox({
    
    create.Progressbox_price(reactives= reactives, "barrier", "optionValue",
                             "blue", "euro")
    
  })
  
  # barrier option sd progress box
  output$progressBox_barrier_sd <- renderValueBox({
    
    create.Progressbox_SD(reactives= reactives, "barrier", "optionSD",
                          "light-blue", "resize-horizontal")
    
  })
  
  # barrier option se progress box
  output$progressBox_barrier_se <- renderValueBox({
    
    create.Progressbox_SE(reactives= reactives, "barrier", "optionSE",
                          "aqua", "exclamation-sign")
    
  })
  
  
  ## histogram of the discounted payoff
  observeEvent(reactives[["warning"]][["barrier"]], {
    
    if(reactives[["warning"]][["barrier"]] == FALSE){
      
      output$barrier_plot <- renderPlot({
        
        create.Histogram(reactives = reactives, "barrier")
        
      },
      height = "auto",width = "auto")
      
    } else if(reactives[["warning"]][["barrier"]] == TRUE){
      
      output$barrier_plot <- NULL
      
    }
    
  })
  
  
  ## conditional panel for maturity in date or years
  output$maturity_barrier <- renderUI({
    
    create.Conditional.Maturity("barrier")
    
  })
  
  
  
  ## binary --------------------------------------------------------------------
  observeEvent(input$binarySimulation, {
    
    # if not fully loaded no calculation is made
    req(input$maturityDate_binary)
    
    set.active.option(reactives, "binary")
    
    checkInput <- check.Input(reactives, input)
    checkInput_na <- check.Input_na(reactives, input) 
    
    reactives[["warning"]][["binary"]] <- FALSE
    
    if(checkInput$check && checkInput_na$check){
      
      maturity <- choose.maturity.Input(input$binary_dateOrYear, 
                                        input$maturityDate_binary, 
                                        input$timeToMaturity_binary)
      
      if(input$dividend_YesOrNo == TRUE){
        
        if(input$simulationModel == "Black-Scholes"){
          withProgress(message = 'Generating stock paths',{
          # calling the path and option function
          stockPaths <- create.StockPaths(simulationModel= input$simulationModel,
                                          initialPrice = input$stockPrice,     maturity = maturity,
                                          volatility = input$volatility/100,   interestRate = input$riskFreeRate/100,
                                          dividendRate = input$dividendRate/100,
                                          sampleSize = input$sampleSize*1000,  numberTimeSteps = input$timeSteps)
          })
        } else if(input$simulationModel == "Heston"){
          withProgress(message = 'Generating stock paths',{
          # calling the path and option function
          stockPaths <- create.StockPaths(simulationModel= input$simulationModel,
                                          initialPrice = input$stockPrice,     maturity = maturity,
                                          volatility = input$volatility/100,   interestRate = input$riskFreeRate/100,
                                          dividendRate = input$dividendRate/100,
                                          sampleSize = input$sampleSize*1000,  numberTimeSteps = input$timeSteps,
                                          kappa = input$kappa, epsilon = input$epsilon,
                                          rho = input$rho, theta = input$theta)
          })
        }
        
      } else if(input$dividend_YesOrNo == FALSE){
        
        if(input$simulationModel == "Black-Scholes"){
          withProgress(message = 'Generating stock paths',{  
          # calling the path and option function
          stockPaths <- create.StockPaths(simulationModel= input$simulationModel,
                                          initialPrice = input$stockPrice,     maturity = maturity,
                                          volatility = input$volatility/100,   interestRate = input$riskFreeRate/100,
                                          sampleSize = input$sampleSize*1000,  numberTimeSteps = input$timeSteps)
          })
        } else if(input$simulationModel == "Heston"){
          withProgress(message = 'Generating stock paths',{
          # calling the path and option function
          stockPaths <- create.StockPaths(simulationModel= input$simulationModel,
                                          initialPrice = input$stockPrice,     maturity = maturity,
                                          volatility = input$volatility/100,   interestRate = input$riskFreeRate/100,
                                          sampleSize = input$sampleSize*1000,  numberTimeSteps = input$timeSteps,
                                          kappa = input$kappa, epsilon = input$epsilon,
                                          rho = input$rho, theta = input$theta)
          })
        }
        
      }
      
      withProgress(message = 'Calculating option value', {
      if(input$subType_binary == "Asset-or-nothing"){
        
        result <- price.Option(type = "binary",
                               optionType = input$optionType_binary,
                               subType = input$subType_binary,
                               binaryPayoff = 0,
                               strikePrice = input$strikePrice_binary_asset,
                               stockPaths = c(stockPaths))
        
      } else if(input$subType_binary == "Cash-or-nothing"){
        
        result <- price.Option(type = "binary",
                               optionType = input$optionType_binary,
                               subType = input$subType_binary,
                               binaryPayoff = input$binaryPayoff,
                               strikePrice = input$strikePrice_binary_cash,
                               stockPaths = c(stockPaths))
        
      }
      
      
      
      reactives$result[["binary"]] <- result
      reactives$payoff_data[["binary"]] <- data.frame(payoff = result[["C0_plot"]])
      })
      
      output$warning_binary <- renderUI({
        
        create.warning(reactives, input, "binary")
        
      })
      
      
    } else if(checkInput$check == FALSE){
      
      showNotification(checkInput$warning, type = "error", duration = 3)
      
      reactives$result[["binary"]] <- NULL
      reactives$payoff_data[["binary"]] <- NULL
      
    }else if(checkInput_na$check == FALSE){
      
      showNotification(checkInput_na$warning, type = "error", duration = 3)
      
      reactives$result[["binary"]] <- NULL
      reactives$payoff_data[["binary"]] <- NULL
      
    }
    
  })
  
  
  ## progressboxes
  # binary price progress box
  output$progressBox_binary_price <- renderValueBox({
    
    create.Progressbox_price(reactives= reactives, "binary", "optionValue",
                             "blue", "euro")
    
  })
  
  # binary option sd progress box
  output$progressBox_binary_sd <- renderValueBox({
    
    create.Progressbox_SD(reactives= reactives, "binary", "optionSD",
                          "light-blue", "resize-horizontal")
    
  })
  
  # binary option se progress box
  output$progressBox_binary_se <- renderValueBox({
    
    create.Progressbox_SE(reactives= reactives, "binary", "optionSE",
                          "aqua", "exclamation-sign")
    
  })
  
  
  ## histogram of the discounted payoff
  observeEvent(reactives[["warning"]][["binary"]], {
    
    if(reactives[["warning"]][["binary"]] == FALSE){
      
      output$binary_plot <- renderPlot({
        
        create.Histogram(reactives = reactives, "binary")
        
      },
      height = "auto",width = "auto")
      
    } else if(reactives[["warning"]][["binary"]] == TRUE){
      
      output$binary_plot <- NULL
      
    }
    
  })
  
  ## conditional panel for maturity in date or years
  output$maturity_binary <- renderUI({
    
    create.Conditional.Maturity("binary")
    
  })
  
  
  ## conditional panel for cash if cashornothing
  output$cash_or_nothing_input <- renderUI({
    
    list(
      conditionalPanel(
        condition = "input.subType_binary == 'Asset-or-nothing'",
        
        numericInput(inputId = "strikePrice_binary_asset",
                     label = "Strike price",
                     value = 40)
        
        
      ),
      conditionalPanel(
        condition = "input.subType_binary == 'Cash-or-nothing'",
        
        splitLayout(
          
          numericInput(inputId = "strikePrice_binary_cash",
                       label = "Strike price",
                       value = 50),
          
          numericInput(inputId = "binaryPayoff",
                       label = "Binary payoff",
                       value = 20,
                       width = "97%")
        )
        
        
        
      )
    )
  })
  
  
  
  ## chooser -------------------------------------------------------------------
  observeEvent(input$chooserSimulation, {
    
    # if not fully loaded no calculation is made
    req(input$chooser_subtype_asian)
    
    set.active.option(reactives, "chooser")
    
    checkInput <- check.Input(reactives, input)
    checkInput_na <- check.Input_na(reactives, input)
    
    reactives[["warning"]][["chooser"]] <- FALSE
    
    if(checkInput$check && checkInput_na$check){

      maturity <- choose.maturity.Input(input$chooser_dateOrYear,
                                        input$maturityDate_chooser,
                                        input$timeToMaturity_chooser)

      
      
      if(input$dividend_YesOrNo == TRUE){

        if(input$simulationModel == "Black-Scholes"){
          withProgress(message = 'Generating stock paths',{
          # calling the path and option function
          stockPaths <- create.StockPaths(simulationModel= input$simulationModel,
                                          initialPrice = input$stockPrice,     maturity = maturity,
                                          volatility = input$volatility/100,   interestRate = input$riskFreeRate/100,
                                          dividendRate = input$dividendRate/100,
                                          sampleSize = input$sampleSize*1000,  skalarTimeSteps = input$skalar_timeSteps)
          })
        } else if(input$simulationModel == "Heston"){
          withProgress(message = 'Generating stock paths',{
          # calling the path and option function
          stockPaths <- create.StockPaths(simulationModel= input$simulationModel,
                                          initialPrice = input$stockPrice,     maturity = maturity,
                                          volatility = input$volatility/100,   interestRate = input$riskFreeRate/100,
                                          dividendRate = input$dividendRate/100,
                                          sampleSize = input$sampleSize*1000,  skalarTimeSteps = input$skalar_timeSteps,
                                          kappa = input$kappa, epsilon = input$epsilon,
                                          rho = input$rho, theta = input$theta)
          })
        }

      } else if(input$dividend_YesOrNo == FALSE){
        
        if(input$simulationModel == "Black-Scholes"){
          withProgress(message = 'Generating stock paths',{
          # calling the path and option function
          stockPaths <- create.StockPaths(simulationModel= input$simulationModel,
                                          initialPrice = input$stockPrice,     maturity = maturity,
                                          volatility = input$volatility/100,   interestRate = input$riskFreeRate/100,
                                          sampleSize = input$sampleSize*1000,  skalarTimeSteps = input$skalar_timeSteps)
          })
        } else if(input$simulationModel == "Heston"){
          withProgress(message = 'Generating stock paths',{
          # calling the path and option function
          stockPaths <- create.StockPaths(simulationModel= input$simulationModel,
                                          initialPrice = input$stockPrice,     maturity = maturity,
                                          volatility = input$volatility/100,   interestRate = input$riskFreeRate/100,
                                          sampleSize = input$sampleSize*1000,  skalarTimeSteps = input$skalar_timeSteps,
                                          kappa = input$kappa, epsilon = input$epsilon,
                                          rho = input$rho, theta = input$theta)
          })
        }

      }
      
      # here the loop for the calculation is pretty big because of the number of cases
      withProgress(message = 'Calculating option value', {
      if(input$chooser_option == "asian"){
        
        result <- price.Option(type="chooser",
                               strikePrice = input$strikePrice_chooser,
                               stockPaths = c(stockPaths),
                               choosingDate = input$chooserDecisionDate,
                               
                               chooserOption ="asian",
                               asianRange = input$observation_asian_chooser,
                               intervalObs_asian = input$daysBetweenObs_asian_chooser,
                               periodObs_asian = input$periodOfObs_asian_chooser,
                               chooserSubType = input$chooser_subtype_asian
        )
        
        
      } else if(input$chooser_option == "barrier"){
        
        result <- price.Option(type="chooser",
                               strikePrice = input$strikePrice_chooser,
                               stockPaths = c(stockPaths),
                               choosingDate = input$chooserDecisionDate,
                               
                               chooserOption ="barrier",
                               chooserSubType = input$chooser_subtype_barrier,
                               barrier = input$barrierValue_chooser
        )
        
        
      } else if(input$chooser_option == "binary"){
        
        result <- price.Option(type="chooser",
                               strikePrice = input$strikePrice_chooser,
                               stockPaths = c(stockPaths),
                               choosingDate = input$chooserDecisionDate,
                               
                               chooserOption ="binary",
                               chooserSubType = input$chooser_subtype_binary,
                               binaryPayoff = input$binaryPayoff_chooser
        )
        
      } else if(input$chooser_option == "european"){
        
        result <- price.Option(type="chooser",
                               strikePrice = input$strikePrice_chooser,
                               stockPaths = c(stockPaths),
                               choosingDate = input$chooserDecisionDate,
                               
                               chooserOption ="european"
        )
        
        
      } else if(input$chooser_option == "forwardStart"){
        
        result <- price.Option(type="chooser",
                               strikePrice = input$strikePrice_chooser,
                               stockPaths = c(stockPaths),
                               choosingDate = input$chooserDecisionDate,
                               
                               chooserOption ="forwardStart",
                               startDate = input$forwardStartDate_chooser
        )
        
      } else if(input$chooser_option == "gap"){
        
        result <- price.Option(type="chooser",
                               stockPaths = c(stockPaths),
                               choosingDate = input$chooserDecisionDate,

                               chooserOption ="gap",
                               strikePrice1 = input$strikePrice_chooser,
                               strikePrice2 = input$triggerPrice_gap_chooser
        )
        
        
      } else if(input$chooser_option == "lookback"){
        
        result <- price.Option(type="chooser",
                               strikePrice = input$strikePrice_chooser,
                               stockPaths = c(stockPaths),
                               choosingDate = input$chooserDecisionDate,
                               
                               chooserOption ="lookback",
                               chooserSubType = input$chooser_subtype_lookback
        )
        
      }
      
      reactives$result[["chooser"]] <- result
      reactives$payoff_data[["chooser"]] <- data.frame(payoff = result[["C0_plot"]])
      })
      
      output$warning_chooser <- renderUI({
        
        create.warning(reactives, input, "chooser")
        
      })
      
    } else if(checkInput$check == FALSE){
      
      showNotification(checkInput$warning, type = "error", duration = 3)
      
      reactives$result[["chooser"]] <- NULL
      reactives$payoff_data[["chooser"]] <- NULL
      
    }else if(checkInput_na$check == FALSE){
      
      showNotification(checkInput_na$warning, type = "error", duration = 3)
      
      reactives$result[["chooser"]] <- NULL
      reactives$payoff_data[["chooser"]] <- NULL
      
    }
    
  })
  
  
  
  ## progressboxes
  # chooser price progress box
  output$progressBox_chooser_price <- renderValueBox({
    
    create.Progressbox_price(reactives= reactives, "chooser", "optionValue",
                             "blue", "euro")
    
  })
  
  # chooser option sd progress box
  output$progressBox_chooser_sd <- renderValueBox({
    
    create.Progressbox_SD(reactives= reactives, "chooser", "optionSD",
                          "light-blue", "resize-horizontal")
    
  })
  
  # chooser option se progress box
  output$progressBox_chooser_se <- renderValueBox({
    
    create.Progressbox_SE(reactives= reactives, "chooser", "optionSE",
                          "aqua", "exclamation-sign")
    
  })
  
  
  ## histogram of the discounted payoff
  observeEvent(reactives[["warning"]][["chooser"]], {
    
    if(reactives[["warning"]][["chooser"]] == FALSE){
      
      output$chooser_plot <- renderPlot({
        
        create.Histogram(reactives = reactives, "chooser")
        
      },
      height = "auto",width = "auto")
      
    } else if(reactives[["warning"]][["chooser"]] == TRUE){
      
      output$chooser_plot <- NULL
      
    }
    
  })
  
  
  ## conditional panel for maturity in date or years
  output$maturity_chooser <- renderUI({
    
    create.Conditional.Maturity("chooser")
    
  })
  
  
  
  ## conditional panel for the different option subtypes
  output$optionPanel_chooser_subtype <- renderUI({

    list(
      conditionalPanel(
        condition = "input.chooser_option == 'asian'",

        selectInput(inputId = "chooser_subtype_asian",
                    label = "Subtype",
                    choices = c("Arithmetic" = "arithmetic", "Geometric" = "geometric"),
                    width = "97%"
                    )

      ),
      conditionalPanel(
        condition = "input.chooser_option == 'barrier'",
        
        selectInput(inputId = "chooser_subtype_barrier",
                    label = "Subtype",
                    choices = c("Up-And-In", "Up-And-Out", "Down-And-In",
                                "Down-And-Out"),
                    width = "97%"
        )

      ),
      conditionalPanel(
        condition = "input.chooser_option == 'binary'",
        
        selectInput(inputId = "chooser_subtype_binary",
                    label = "Subtype",
                    choices = c("Asset-or-nothing", "Cash-or-nothing"),
                    width = "97%"
        )

      ),
      conditionalPanel(
        condition = "input.chooser_option == 'european'",
        
        selectInput(inputId = "chooser_subtype_european",
                    label = "Subtype",
                    choices = c("No subtype"),
                    width = "97%"
        )

      ),
      conditionalPanel(
        condition = "input.chooser_option == 'forwardStart'",
        
        selectInput(inputId = "chooser_subtype_forwardStart",
                    label = "Subtype",
                    choices = c("No subtype"),
                    width = "97%"
        )

      ),
      conditionalPanel(
        condition = "input.chooser_option == 'gap'",
        
        selectInput(inputId = "chooser_subtype_gap",
                    label = "Subtype",
                    choices = c("No subtype"),
                    width = "97%"
        )

      ),
      conditionalPanel(
        condition = "input.chooser_option == 'lookback'",
        
        selectInput(inputId = "chooser_subtype_lookback",
                    label = "Subtype",
                    choices = c("Fixed", "Floating"),
                    width = "97%"
        )

      )
    )
  })
  
  
  ## conditional panel for the different option subtype parameters
  output$optionPanel_chooser_optionvalues <- renderUI({

    list(
      conditionalPanel(
        condition = "input.chooser_option == 'asian'",

        splitLayout(
          
          radioButtons(
            inputId = "observation_asian_chooser",
            label = "Range",
            choices = c("Total","Interval",
                        "Period"),
            inline = FALSE
          ),
          
          uiOutput("date_range_asian_chooser")
          
        )
        

      ),
      conditionalPanel(
        condition = "input.chooser_option == 'barrier'",

        numericInput(inputId = "barrierValue_chooser",
                     label = "Barrier",
                     value = 55)

      ),
      conditionalPanel(
        condition = "input.chooser_option == 'binary'",
        
        numericInput(inputId = "binaryPayoff_chooser",
                     label = tags$span(
                       "Binary payoff ", 
                       tags$i(
                         class = "glyphicon glyphicon-info-sign", 
                         style = "color:#0072B2;",
                         title = 'Is only used if "Cash-or-nothing" is selected.' 
                       )
                     ),
                     value = 20,
                     width = "100%")

      ),
      conditionalPanel(
        condition = "input.chooser_option == 'forwardStart'",

        dateInput(inputId = "forwardStartDate_chooser",
                  label = "Start date",
                  value = (Sys.Date() + 6),
                  daysofweekdisabled = c(0,6),
                  datesdisabled = create.Holidays()
        )

      ),
      conditionalPanel(
        condition = "input.chooser_option == 'gap'",

        numericInput(inputId = "triggerPrice_gap_chooser",
                     label = "Trigger price",
                     value = 50)

      )
    )
  })
  
  
  ## conditional panel for total range or partial range of the asian chooser
  output$date_range_asian_chooser <- renderUI({

    list(
      conditionalPanel(
        condition = "input.observation_asian_chooser == 'Interval'",

        numericInput(inputId = "daysBetweenObs_asian_chooser",
                     label = tags$span(
                       "Days between observations ", 
                       tags$i(
                         class = "glyphicon glyphicon-info-sign", 
                         style = "color:#0072B2;",
                         title = 'Interval between observations in days (incl. weekends).' 
                       )
                     ),
                     value = 1,
                     min = 1,
                     width = "95%")

      ),
      conditionalPanel(
        condition = "input.observation_asian_chooser == 'Period'",

        numericInput(inputId = "periodOfObs_asian_chooser",
                     label = tags$span(
                       "Days before expiry ", 
                       tags$i(
                         class = "glyphicon glyphicon-info-sign", 
                         style = "color:#0072B2;",
                         title = 'Number of days from the period of observation.' 
                       )
                     ),
                     value = 90,
                     min = 2,
                     width = "95%")

      )
    )
  })
  
  ## compound ------------------------------------------------------------------
  observeEvent(input$compoundSimulation, {
    
    # if not fully loaded no calculation is made
    req(input$compound_subtype_asian)
    
    set.active.option(reactives, "compound")

    checkInput <- check.Input(reactives, input)
    checkInput_na <- check.Input_na(reactives, input)
    
    reactives[["warning"]][["compound"]] <- FALSE
    
    if(checkInput$check && checkInput_na$check){

      maturity_first <- choose.maturity.Input(input$compound_first_dateOrYear,
                                              input$maturityDate_compound_first,
                                              input$timeToMaturity_compound_first)

      maturity_second <- choose.maturity.Input(input$compound_second_dateOrYear,
                                               input$maturityDate_compound_second,
                                               input$timeToMaturity_compound_second)
      
      withProgress(message = 'Generating  stock paths', {
      
        if(input$dividend_YesOrNo == TRUE){
          
          if(input$simulationModel == "Black-Scholes"){
            
            # calling the path and option function
            stockPaths <- create.StockPaths(simulationModel= input$simulationModel,
                                            initialPrice = input$stockPrice,     maturity = maturity_first,
                                            volatility = input$volatility/100,   interestRate = input$riskFreeRate/100,
                                            dividendRate = input$dividendRate/100,
                                            sampleSize = input$sampleSize_compound,  numberTimeSteps = 1)
            
          } else if(input$simulationModel == "Heston"){
            
            # calling the path and option function
            stockPaths <- create.StockPaths(simulationModel= input$simulationModel,
                                            initialPrice = input$stockPrice,     maturity = maturity_first,
                                            volatility = input$volatility/100,   interestRate = input$riskFreeRate/100,
                                            sampleSize = input$sampleSize_compound, 
                                            dividendRate = input$dividendRate/100,
                                            numberTimeSteps = 1,
                                            kappa = input$kappa, epsilon = input$epsilon,
                                            rho = input$rho, theta = input$theta)
            
          }
          
        } else if(input$dividend_YesOrNo == FALSE){
          
          if(input$simulationModel == "Black-Scholes"){
            
            # calling the path and option function
            stockPaths <- create.StockPaths(simulationModel= input$simulationModel,
                                            initialPrice = input$stockPrice,     maturity = maturity_first,
                                            volatility = input$volatility/100,   interestRate = input$riskFreeRate/100,
                                            sampleSize = input$sampleSize_compound,  numberTimeSteps = 1)
            
            
            
          } else if(input$simulationModel == "Heston"){
            
            # calling the path and option function
            stockPaths <- create.StockPaths(simulationModel= input$simulationModel,
                                            initialPrice = input$stockPrice,     maturity = maturity_first,
                                            volatility = input$volatility/100,   interestRate = input$riskFreeRate/100,
                                            sampleSize = input$sampleSize_compound,  numberTimeSteps = 1,
                                            kappa = input$kappa, epsilon = input$epsilon,
                                            rho = input$rho, theta = input$theta)
          }
          
        }
      })
      
      withProgress(message = 'Calculating option value',{
        
        # as in the chooser option the loop for the option price is pretty big due to the different options
        if(input$compound_option_second == "asian"){
          
          result <- price.Option(type="compound",
                                 type2="asian",
                                 subType = input$optionType_compound,
                                 strikePrice = input$strikePrice_compound_first,
                                 stockPaths = c(stockPaths),
                                 
                                 #Second option parameters
                                 secondSampleSize = input$sampleSize_compound_second,
                                 secondSkalarTimeSteps = input$skalar_timeSteps,
                                 secondExerciseDate = maturity_second,
                                 secondSubType = input$compound_subtype_asian,
                                 
                                 secondStrikePrice = input$strikePrice_compound_second,
                                 
                                 asianRange = input$observation_asian_compound,
                                 intervalObs_asian = input$daysBetweenObs_asian_compound,
                                 periodObs_asian = input$periodOfObs_asian_compound)
          
        } 
        else if(input$compound_option_second == "barrier"){
          
          result <- price.Option(type="compound",
                                 type2="barrier",
                                 subType = input$optionType_compound,
                                 strikePrice = input$strikePrice_compound_first,
                                 stockPaths = c(stockPaths),
                                 
                                 #Second option parameters
                                 secondSampleSize = input$sampleSize_compound_second,
                                 secondNumberTimeSteps = input$timeSteps,
                                 secondExerciseDate = maturity_second,
                                 secondSubType = input$compound_subtype_barrier,
                                 
                                 secondStrikePrice = input$strikePrice_compound_second,
                                 
                                 barrier = input$barrierValue_compound)
          
        } else if(input$compound_option_second == "binary"){
          
          result <- price.Option(type="compound",
                                 type2="binary",
                                 subType = input$optionType_compound,
                                 strikePrice = input$strikePrice_compound_first,
                                 stockPaths = c(stockPaths),
                                 
                                 #Second option parameters
                                 secondSampleSize = input$sampleSize_compound_second,
                                 secondNumberTimeSteps = input$timeSteps,
                                 secondExerciseDate = maturity_second,
                                 secondSubType = input$compound_subtype_binary,
                                 
                                 secondStrikePrice = input$strikePrice_compound_second,
                                 
                                 binaryPayoff = input$binaryPayoff_compound)
          
        } else if(input$compound_option_second == "european"){
          
          result <- price.Option(type="compound",
                                 type2="european",
                                 subType = input$optionType_compound,
                                 strikePrice = input$strikePrice_compound_first,
                                 stockPaths = c(stockPaths),
                                 
                                 #Second option parameters
                                 secondSampleSize = input$sampleSize_compound_second,
                                 secondNumberTimeSteps = 1,
                                 secondExerciseDate = maturity_second,
                                 secondSubType = input$compound_subtype_european,
                                 
                                 secondStrikePrice = input$strikePrice_compound_second
          )
          
        } else if(input$compound_option_second == "forwardStart"){
          
          result <- price.Option(type="compound",
                                 type2="forwardStart",
                                 subType = input$optionType_compound,
                                 strikePrice = input$strikePrice_compound_first,
                                 stockPaths = c(stockPaths),
                                 
                                 #Second option parameters
                                 secondSampleSize = input$sampleSize_compound_second,
                                 secondNumberTimeSteps = input$timeSteps,
                                 secondExerciseDate = maturity_second,
                                 secondSubType = input$compound_subtype_forwardStart,
                                 
                                 secondStrikePrice = input$strikePrice_compound_second,
                                 
                                 startDate = input$forwardStartDate_compound)
          
        } else if(input$compound_option_second == "gap"){
          
          result <- price.Option(type="compound",
                                 type2="gap",
                                 subType = input$optionType_compound,
                                 strikePrice = input$strikePrice_compound_first,
                                 stockPaths = c(stockPaths),
                                 
                                 #Second option parameters
                                 secondSampleSize = input$sampleSize_compound_second,
                                 secondNumberTimeSteps = input$timeSteps,
                                 secondExerciseDate = maturity_second,
                                 secondSubType = input$compound_subtype_gap,
                                 
                                 # secondStrikePrice = input$strikePrice_compound_second,
                                 
                                 strikePrice1 = input$strikePrice_compound_second,
                                 strikePrice2 = input$triggerPrice_gap_compound)
          
          
          
        } else if(input$compound_option_second == "lookback"){
          
          result <- price.Option(type="compound",
                                 type2="lookback",
                                 subType = input$optionType_compound,
                                 strikePrice = input$strikePrice_compound_first,
                                 stockPaths = c(stockPaths),
                                 
                                 #Second option parameters
                                 secondSampleSize = input$sampleSize_compound_second,
                                 secondSkalarTimeSteps = input$skalar_timeSteps,
                                 secondExerciseDate = maturity_second,
                                 secondSubType = input$compound_subtype_lookback,
                                 
                                 secondStrikePrice = input$strikePrice_compound_second)
          
          
          
        }
      })
      reactives$result[["compound"]] <- result
      reactives$payoff_data[["compound"]] <- data.frame(payoff = result[["C0_plot"]])


      output$warning_compound <- renderUI({

        create.warning(reactives, input, "compound")

      })


    } else if(checkInput$check == FALSE){

      showNotification(checkInput$warning, type = "error", duration = 3)

      reactives$result[["compound"]] <- NULL
      reactives$payoff_data[["compound"]] <- NULL

    }else if(checkInput_na$check == FALSE){
      
      showNotification(checkInput_na$warning, type = "error", duration = 3)
      
      reactives$result[["compound"]] <- NULL
      reactives$payoff_data[["compound"]] <- NULL
      
    }
    
  })
  
  
  ## progressboxes
  # compound price progress box
  output$progressBox_compound_price <- renderValueBox({
    
    create.Progressbox_price(reactives= reactives, "compound", "optionValue",
                             "blue", "euro")
    
  })
  
  # compound option sd progress box
  output$progressBox_compound_sd <- renderValueBox({
    
    create.Progressbox_SD(reactives= reactives, "compound", "optionSD",
                          "light-blue", "resize-horizontal")
    
  })
  
  # compound option se progress box
  output$progressBox_compound_se <- renderValueBox({
    
    create.Progressbox_SE(reactives= reactives, "compound", "optionSE",
                          "aqua", "exclamation-sign")
    
  })
  
  
  ## histogram of the discounted payoff
  observeEvent(reactives[["warning"]][["compound"]], {
    
    if(reactives[["warning"]][["compound"]] == FALSE){
      
      output$compound_plot <- renderPlot({
        
        create.Histogram(reactives = reactives, "compound")
        
      },
      height = "auto",width = "auto")
      
    } else if(reactives[["warning"]][["compound"]] == TRUE){
      
      output$compound_plot <- NULL
      
    }
    
  })
  
  
  # conditional panel for maturity in date or years
  output$maturity_compound_first <- renderUI({
    
    create.Conditional.Maturity("compound_first")
    
  })
  
  # conditional panel for maturity in date or years
  output$maturity_compound_second <- renderUI({
    
    create.Conditional.Maturity("compound_second")
    
  })
  
  
  
  ## conditional panel for the different option subtypes
  output$optionPanel_compound_subtype <- renderUI({
    
    list(
      conditionalPanel(
        condition = "input.compound_option_second == 'asian'",
        
        selectInput(inputId = "compound_subtype_asian",
                    label = "Subtype",
                    choices = c("Arithmetic" = "arithmetic", "Geometric" = "geometric"),
                    width = "97%"
        )
        
      ),
      conditionalPanel(
        condition = "input.compound_option_second == 'barrier'",
        
        selectInput(inputId = "compound_subtype_barrier",
                    label = "Subtype",
                    choices = c("Up-And-In", "Up-And-Out", "Down-And-In",
                                "Down-And-Out"),
                    width = "97%"
        )
        
      ),
      conditionalPanel(
        condition = "input.compound_option_second == 'binary'",
        
        selectInput(inputId = "compound_subtype_binary",
                    label = "Subtype",
                    choices = c("Asset-or-nothing", "Cash-or-nothing"),
                    width = "97%"
        )
        
      ),
      conditionalPanel(
        condition = "input.compound_option_second == 'european'",
        
        selectInput(inputId = "compound_subtype_european",
                    label = "Subtype",
                    choices = c("No subtype"),
                    width = "97%"
        )
        
      ),
      conditionalPanel(
        condition = "input.compound_option_second == 'forwardStart'",
        
        selectInput(inputId = "compound_subtype_forwardStart",
                    label = "Subtype",
                    choices = c("No subtype"),
                    width = "97%"
        )
        
      ),
      conditionalPanel(
        condition = "input.compound_option_second == 'gap'",
        
        selectInput(inputId = "compound_subtype_gap",
                    label = "Subtype",
                    choices = c("No subtype"),
                    width = "97%"
        )
        
      ),
      conditionalPanel(
        condition = "input.compound_option_second == 'lookback'",
        
        selectInput(inputId = "compound_subtype_lookback",
                    label = "Subtype",
                    choices = c("Fixed", "Floating"),
                    width = "97%"
        )
        
      )
    )
  })
  
  
  ## conditional panel for the different option subtype parameters
  output$optionPanel_compound_optionvalues <- renderUI({
    
    list(
      conditionalPanel(
        condition = "input.compound_option_second == 'asian'",
        
        splitLayout(
          
          radioButtons(
            inputId = "observation_asian_compound",
            label = "Range",
            choices = c("Total","Interval",
                        "Period"),
            inline = FALSE
          ),
          
          uiOutput("date_range_asian_compound")
          
        )
        
        
      ),
      conditionalPanel(
        condition = "input.compound_option_second == 'barrier'",
        
        numericInput(inputId = "barrierValue_compound",
                     label = "Barrier",
                     value = 55)
        
      ),
      conditionalPanel(
        condition = "input.compound_option_second == 'binary'",
        
        numericInput(inputId = "binaryPayoff_compound",
                     label = tags$span(
                       "Binary payoff ", 
                       tags$i(
                         class = "glyphicon glyphicon-info-sign", 
                         style = "color:#0072B2;",
                         title = 'Is only used if "Cash-or-nothing" is selected.' 
                       )
                     ),
                     value = 20,
                     width = "100%")
        
      ),
      conditionalPanel(
        condition = "input.compound_option_second == 'forwardStart'",
        
        dateInput(inputId = "forwardStartDate_compound",
                  label = "Start date",
                  value = (Sys.Date() + 6),
                  daysofweekdisabled = c(0,6),
                  datesdisabled = create.Holidays()
        )
        
      ),
      conditionalPanel(
        condition = "input.compound_option_second == 'gap'",
        
        numericInput(inputId = "triggerPrice_gap_compound",
                     label = "Trigger price",
                     value = 50)
        
      )
    )
  })
  
  
  ## conditional panel for total range or partial range of the asian compound
  output$date_range_asian_compound <- renderUI({
    
    list(
      conditionalPanel(
        condition = "input.observation_asian_compound == 'Interval'",
        
        numericInput(inputId = "daysBetweenObs_asian_compound",
                     label = tags$span(
                       "Days between observations ", 
                       tags$i(
                         class = "glyphicon glyphicon-info-sign", 
                         style = "color:#0072B2;",
                         title = 'Interval between observations in days (incl. weekends).' 
                       )
                     ),
                     value = 1,
                     min = 1,
                     width = "95%")
        
      ),
      conditionalPanel(
        condition = "input.observation_asian_compound == 'Period'",
        
        numericInput(inputId = "periodOfObs_asian_compound",
                     label = tags$span(
                       "Days before expiry ", 
                       tags$i(
                         class = "glyphicon glyphicon-info-sign", 
                         style = "color:#0072B2;",
                         title = 'Number of days from the period of observation.' 
                       )
                     ),
                     value = 90,
                     min = 2,
                     width = "95%")
        
      )
    )
  })
  
  
  ## european  -----------------------------------------------------------
  observeEvent(input$europeanSimulation, {
    
    # if not fully loaded no calculation is made
    req(input$maturityDate_european)
    
    set.active.option(reactives, "european")
    
    checkInput <- check.Input(reactives, input)
    checkInput_na <- check.Input_na(reactives, input) 
    
    reactives[["warning"]][["european"]] <- FALSE
    
    if(checkInput$check && checkInput_na$check){
      
      maturity <- choose.maturity.Input(input$european_dateOrYear, 
                                        input$maturityDate_european, 
                                        input$timeToMaturity_european)
      withProgress(message = 'Generating stock paths',{
        if(input$dividend_YesOrNo == TRUE){
          
          if(input$simulationModel == "Black-Scholes"){
            
            # calling the path and option function
            stockPaths <- create.StockPaths(simulationModel= input$simulationModel,
                                            initialPrice = input$stockPrice,     maturity = maturity,
                                            volatility = input$volatility/100,   interestRate = input$riskFreeRate/100,
                                            dividendRate = input$dividendRate/100,
                                            sampleSize = input$sampleSize*1000,  numberTimeSteps = 1)
            
          } else if(input$simulationModel == "Heston"){
            
            # calling the path and option function
            stockPaths <- create.StockPaths(simulationModel= input$simulationModel,
                                            initialPrice = input$stockPrice,     maturity = maturity,
                                            volatility = input$volatility/100,   interestRate = input$riskFreeRate/100,
                                            dividendRate = input$dividendRate/100,
                                            sampleSize = input$sampleSize*1000,  numberTimeSteps = 1,
                                            kappa = input$kappa, epsilon = input$epsilon,
                                            rho = input$rho, theta = input$theta)
            
          }
          
        } else if(input$dividend_YesOrNo == FALSE){
          
          if(input$simulationModel == "Black-Scholes"){
            
            # calling the path and option function
            stockPaths <- create.StockPaths(simulationModel= input$simulationModel,
                                            initialPrice = input$stockPrice,     maturity = maturity,
                                            volatility = input$volatility/100,   interestRate = input$riskFreeRate/100,
                                            sampleSize = input$sampleSize*1000,  numberTimeSteps = 1)
            
          } else if(input$simulationModel == "Heston"){
            
            # calling the path and option function
            stockPaths <- create.StockPaths(simulationModel= input$simulationModel,
                                            initialPrice = input$stockPrice,     maturity = maturity,
                                            volatility = input$volatility/100,   interestRate = input$riskFreeRate/100,
                                            sampleSize = input$sampleSize*1000,  numberTimeSteps = 1,
                                            kappa = input$kappa, epsilon = input$epsilon,
                                            rho = input$rho, theta = input$theta)
            
          }
          
        }
      })
      
      withProgress(message = 'Calculating option value', {
      result <- price.Option(type = "european",
                             optionType = input$optionType_european,
                             strikePrice = input$strikePrice_european,
                             stockPaths = c(stockPaths))
      
      reactives$result[["european"]] <- result
      reactives$payoff_data[["european"]] <- data.frame(payoff = result[["C0_plot"]])
      })
      
      output$warning_european <- renderUI({
        
        create.warning(reactives, input, "european")
        
      })
        
      
    } else if(checkInput$check == FALSE){
      
      showNotification(checkInput$warning, type = "error", duration = 3)
      
      reactives$result[["european"]] <- NULL
      reactives$payoff_data[["european"]] <- NULL
      
    }else if(checkInput_na$check == FALSE){
      
      showNotification(checkInput_na$warning, type = "error", duration = 3)
      
      reactives$result[["european"]] <- NULL
      reactives$payoff_data[["european"]] <- NULL
      
    }
    
  })
  
  ## progressboxes
  # european price progress box
  output$progressBox_european_price <- renderValueBox({
    
    create.Progressbox_price(reactives= reactives, "european", "optionValue",
                             "blue", "euro")
    
  })
  
  # european option sd progress box
  output$progressBox_european_sd <- renderValueBox({
    
    create.Progressbox_SD(reactives= reactives, "european", "optionSD",
                          "light-blue", "resize-horizontal")
    
  })
  
  # european option se progress box
  output$progressBox_european_se <- renderValueBox({
    
    create.Progressbox_SE(reactives= reactives, "european", "optionSE",
                          "aqua", "exclamation-sign")
    
  })
  
  
  ## histogram of the discounted payoff
  observeEvent(reactives[["warning"]][["european"]], {
    
    if(reactives[["warning"]][["european"]] == FALSE){
      
      output$european_plot <- renderPlot({
        
        create.Histogram(reactives = reactives, "european")
        
      },
      height = "auto",width = "auto")
      
    } else if(reactives[["warning"]][["european"]] == TRUE){
      
      output$european_plot <- NULL
      
    }
    
  })
  
  
  ## conditional panel for maturity in date or years
  output$maturity_european <- renderUI({
    
    create.Conditional.Maturity("european")
    
  })
  
  
  ## forwardStart --------------------------------------------------------------
  observeEvent(input$forwardStartSimulation, {
    
    # if not fully loaded no calculation is made
    req(input$maturityDate_forwardStart)
    
    set.active.option(reactives, "forwardStart")
    
    checkInput <- check.Input(reactives, input)
    checkInput_na <- check.Input_na(reactives, input) 
    
    reactives[["warning"]][["forwardStart"]] <- FALSE
    
    if(checkInput$check && checkInput_na$check){
      
      maturity <- choose.maturity.Input(input$forwardStart_dateOrYear, 
                                        input$maturityDate_forwardStart, 
                                        input$timeToMaturity_forwardStart)
      
      withProgress(message = 'Generating stock paths',{
      if(input$dividend_YesOrNo == TRUE){
        
        if(input$simulationModel == "Black-Scholes"){
          
          # calling the path and option function
          stockPaths <- create.StockPaths(simulationModel= input$simulationModel,
                                          initialPrice = input$stockPrice,     maturity = maturity,
                                          volatility = input$volatility/100,   interestRate = input$riskFreeRate/100,
                                          dividendRate = input$dividendRate/100,
                                          sampleSize = input$sampleSize*1000,  numberTimeSteps = input$timeSteps)
          
        } else if(input$simulationModel == "Heston"){
          
          # calling the path and option function
          stockPaths <- create.StockPaths(simulationModel= input$simulationModel,
                                          initialPrice = input$stockPrice,     maturity = maturity,
                                          volatility = input$volatility/100,   interestRate = input$riskFreeRate/100,
                                          dividendRate = input$dividendRate/100,
                                          sampleSize = input$sampleSize*1000,  numberTimeSteps = input$timeSteps,
                                          kappa = input$kappa, epsilon = input$epsilon,
                                          rho = input$rho, theta = input$theta)
          
        }
        
      } else if(input$dividend_YesOrNo == FALSE){
        
        if(input$simulationModel == "Black-Scholes"){
          
          # calling the path and option function
          stockPaths <- create.StockPaths(simulationModel= input$simulationModel,
                                          initialPrice = input$stockPrice,     maturity = maturity,
                                          volatility = input$volatility/100,   interestRate = input$riskFreeRate/100,
                                          sampleSize = input$sampleSize*1000,  numberTimeSteps = input$timeSteps)
          
        } else if(input$simulationModel == "Heston"){
          
          # calling the path and option function
          stockPaths <- create.StockPaths(simulationModel= input$simulationModel,
                                          initialPrice = input$stockPrice,     maturity = maturity,
                                          volatility = input$volatility/100,   interestRate = input$riskFreeRate/100,
                                          sampleSize = input$sampleSize*1000,  numberTimeSteps = input$timeSteps,
                                          kappa = input$kappa, epsilon = input$epsilon,
                                          rho = input$rho, theta = input$theta)
          
        }
        
      }
      })
      
      withProgress(message = 'Calculating option value', {
      result <- price.Option(type = "forwardStart",
                             optionType = input$optionType_forwardStart,
                             strikePrice = input$strikePrice_forwardStart,
                             startDate = input$forwardStartDate,
                             stockPaths = c(stockPaths))
      
      reactives$result[["forwardStart"]] <- result
      reactives$payoff_data[["forwardStart"]] <- data.frame(payoff = result[["C0_plot"]])
      })
      
      output$warning_forwardStart <- renderUI({
        
        create.warning(reactives, input, "forwardStart")
        
      })
      
      
    } else if(checkInput$check == FALSE){
      
      showNotification(checkInput$warning, type = "error", duration = 3)
      
      reactives$result[["forwardStart"]] <- NULL
      reactives$payoff_data[["forwardStart"]] <- NULL
      
    }else if(checkInput_na$check == FALSE){
      
      showNotification(checkInput_na$warning, type = "error", duration = 3)
      
      reactives$result[["forwardStart"]] <- NULL
      reactives$payoff_data[["forwardStart"]] <- NULL
      
    }
    
  })
  
  
  ## progressboxes
  # forwardStart price progress box
  output$progressBox_forwardStart_price <- renderValueBox({
    
    create.Progressbox_price(reactives= reactives, "forwardStart", "optionValue",
                             "blue", "euro")
    
  })
  
  # forwardStart option sd progress box
  output$progressBox_forwardStart_sd <- renderValueBox({
    
    create.Progressbox_SD(reactives= reactives, "forwardStart", "optionSD",
                          "light-blue", "resize-horizontal")
    
  })
  
  # forwardStart option se progress box
  output$progressBox_forwardStart_se <- renderValueBox({
    
    create.Progressbox_SE(reactives= reactives, "forwardStart", "optionSE",
                          "aqua", "exclamation-sign")
    
  })
  
  ## histogram of the discounted payoff
  observeEvent(reactives[["warning"]][["forwardStart"]], {
    
    if(reactives[["warning"]][["forwardStart"]] == FALSE){
      
      output$forwardStart_plot <- renderPlot({
        
        create.Histogram(reactives = reactives, "forwardStart")
        
      },
      height = "auto",width = "auto")
      
    } else if(reactives[["warning"]][["forwardStart"]] == TRUE){
      
      output$forwardStart_plot <- NULL
      
    }
    
  })
  
  ## conditional panel for maturity in date or years
  output$maturity_forwardStart <- renderUI({
    
    create.Conditional.Maturity("forwardStart")
    
  })
  
  
  
  ## gap -----------------------------------------------------------------------
  observeEvent(input$gapSimulation, {
    
    # if not fully loaded no calculation is made
    req(input$maturityDate_gap)
    
    set.active.option(reactives, "gap")
    
    checkInput <- check.Input(reactives, input)
    checkInput_na <- check.Input_na(reactives, input) 
    
    reactives[["warning"]][["gap"]] <- FALSE
    
    if(checkInput$check && checkInput_na$check){
      
      maturity <- choose.maturity.Input(input$gap_dateOrYear, 
                                        input$maturityDate_gap, 
                                        input$timeToMaturity_gap)
      
      withProgress(message = 'Generating stock paths',{
        
      if(input$dividend_YesOrNo == TRUE){
        
        if(input$simulationModel == "Black-Scholes"){
          
          # calling the path and option function
          stockPaths <- create.StockPaths(simulationModel= input$simulationModel,
                                          initialPrice = input$stockPrice,     maturity = maturity,
                                          volatility = input$volatility/100,   interestRate = input$riskFreeRate/100,
                                          dividendRate = input$dividendRate/100,
                                          sampleSize = input$sampleSize*1000,  numberTimeSteps = input$timeSteps)
          
        } else if(input$simulationModel == "Heston"){
          
          # calling the path and option function
          stockPaths <- create.StockPaths(simulationModel= input$simulationModel,
                                          initialPrice = input$stockPrice,     maturity = maturity,
                                          volatility = input$volatility/100,   interestRate = input$riskFreeRate/100,
                                          dividendRate = input$dividendRate/100,
                                          sampleSize = input$sampleSize*1000,  numberTimeSteps = input$timeSteps,
                                          kappa = input$kappa, epsilon = input$epsilon,
                                          rho = input$rho, theta = input$theta)
          
        }
        
      } else if(input$dividend_YesOrNo == FALSE){
        
        if(input$simulationModel == "Black-Scholes"){
          
          # calling the path and option function
          stockPaths <- create.StockPaths(simulationModel= input$simulationModel,
                                          initialPrice = input$stockPrice,     maturity = maturity,
                                          volatility = input$volatility/100,   interestRate = input$riskFreeRate/100,
                                          sampleSize = input$sampleSize*1000,  numberTimeSteps = input$timeSteps)
          
        } else if(input$simulationModel == "Heston"){
          
          # calling the path and option function
          stockPaths <- create.StockPaths(simulationModel= input$simulationModel,
                                          initialPrice = input$stockPrice,     maturity = maturity,
                                          volatility = input$volatility/100,   interestRate = input$riskFreeRate/100,
                                          sampleSize = input$sampleSize*1000,  numberTimeSteps = input$timeSteps,
                                          kappa = input$kappa, epsilon = input$epsilon,
                                          rho = input$rho, theta = input$theta)
          
        }
        
      }
      })
      
      withProgress(message = 'Calculating option value', {
      result <- price.Option(type = "gap",
                             optionType = input$optionType_gap,
                             strikePrice1 = input$strikePrice_gap1,
                             strikePrice2 = input$strikePrice_gap2,
                             stockPaths = c(stockPaths))
      
      reactives$result[["gap"]] <- result
      reactives$payoff_data[["gap"]] <- data.frame(payoff = result[["C0_plot"]])
      })
      
      output$warning_gap <- renderUI({
        
        create.warning(reactives, input, "gap")
        
      })
      
      
    } else if(checkInput$check == FALSE){
      
      showNotification(checkInput$warning, type = "error", duration = 3)
      
      reactives$result[["gap"]] <- NULL
      reactives$payoff_data[["gap"]] <- NULL
      
    }else if(checkInput_na$check == FALSE){
      
      showNotification(checkInput_na$warning, type = "error", duration = 3)
      
      reactives$result[["gap"]] <- NULL
      reactives$payoff_data[["gap"]] <- NULL
      
    }
    
  })
  
  
  
  ## progressboxes
  # gap price progress box
  output$progressBox_gap_price <- renderValueBox({
    
    create.Progressbox_price(reactives= reactives, "gap", "optionValue",
                             "blue", "euro")
    
  })
  
  # gap option sd progress box
  output$progressBox_gap_sd <- renderValueBox({
    
    create.Progressbox_SD(reactives= reactives, "gap", "optionSD",
                          "light-blue", "resize-horizontal")
    
  })
  
  # gap option se progress box
  output$progressBox_gap_se <- renderValueBox({
    
    create.Progressbox_SE(reactives= reactives, "gap", "optionSE",
                          "aqua", "exclamation-sign")
    
  })
  
  
  ## histogram of the discounted payoff
  observeEvent(reactives[["warning"]][["gap"]], {
    
    if(reactives[["warning"]][["gap"]] == FALSE){
      
      output$gap_plot <- renderPlot({
        
        create.Histogram(reactives = reactives, "gap")
        
      },
      height = "auto",width = "auto")
      
    } else if(reactives[["warning"]][["gap"]] == TRUE){
      
      output$gap_plot <- NULL
      
    }
    
  })
  
  
  ## conditional panel for maturity in date or years
  output$maturity_gap <- renderUI({
    
    create.Conditional.Maturity("gap")
    
  })
  
  
  
  ## lookback ------------------------------------------------------------------
  observeEvent(input$lookbackSimulation, {
    
    # if not fully loaded no calculation is made
    req(input$maturityDate_lookback)
    
    set.active.option(reactives, "lookback")
    
    checkInput <- check.Input(reactives, input)
    checkInput_na <- check.Input_na(reactives, input) 
    
    reactives[["warning"]][["lookback"]] <- FALSE
    
    if(checkInput$check && checkInput_na$check){
      
      maturity <- choose.maturity.Input(input$lookback_dateOrYear, 
                                        input$maturityDate_lookback, 
                                        input$timeToMaturity_lookback)
      
      withProgress(message = 'Generating stock paths',{
        
      if(input$dividend_YesOrNo == TRUE){
        
        if(input$simulationModel == "Black-Scholes"){
          
          # calling the path and option function
          stockPaths <- create.StockPaths(simulationModel= input$simulationModel,
                                          initialPrice = input$stockPrice,     maturity = maturity,
                                          volatility = input$volatility/100,   interestRate = input$riskFreeRate/100,
                                          dividendRate = input$dividendRate/100,
                                          sampleSize = input$sampleSize*1000,  skalarTimeSteps = input$skalar_timeSteps)
          
        } else if(input$simulationModel == "Heston"){
          
          # calling the path and option function
          stockPaths <- create.StockPaths(simulationModel= input$simulationModel,
                                          initialPrice = input$stockPrice,     maturity = maturity,
                                          volatility = input$volatility/100,   interestRate = input$riskFreeRate/100,
                                          dividendRate = input$dividendRate/100,
                                          sampleSize = input$sampleSize*1000,  skalarTimeSteps = input$skalar_timeSteps,
                                          kappa = input$kappa, epsilon = input$epsilon,
                                          rho = input$rho, theta = input$theta)
          
        }
        
      } else if(input$dividend_YesOrNo == FALSE){
        
        if(input$simulationModel == "Black-Scholes"){
          
          # calling the path and option function
          stockPaths <- create.StockPaths(simulationModel= input$simulationModel,
                                          initialPrice = input$stockPrice,     maturity = maturity,
                                          volatility = input$volatility/100,   interestRate = input$riskFreeRate/100,
                                          sampleSize = input$sampleSize*1000,  skalarTimeSteps = input$skalar_timeSteps)
          
        } else if(input$simulationModel == "Heston"){
          
          # calling the path and option function
          stockPaths <- create.StockPaths(simulationModel= input$simulationModel,
                                          initialPrice = input$stockPrice,     maturity = maturity,
                                          volatility = input$volatility/100,   interestRate = input$riskFreeRate/100,
                                          sampleSize = input$sampleSize*1000,  skalarTimeSteps = input$skalar_timeSteps,
                                          kappa = input$kappa, epsilon = input$epsilon,
                                          rho = input$rho, theta = input$theta)
          
        }
        
      }
      })
      
      withProgress(message = 'Calculating option value', {
      result <- price.Option(type = "lookback",
                             optionType = input$optionType_lookback,
                             subType = input$subType_lookback,
                             strikePrice = input$strikePrice_lookback,
                             stockPaths = c(stockPaths))
      
      reactives$result[["lookback"]] <- result
      reactives$payoff_data[["lookback"]] <- data.frame(payoff = result[["C0_plot"]])
      })
      
      output$warning_lookback <- renderUI({
        
        create.warning(reactives, input, "lookback")
        
      })
      
      
    } else if(checkInput$check == FALSE){
      
      showNotification(checkInput$warning, type = "error", duration = 3)
      
      reactives$result[["lookback"]] <- NULL
      reactives$payoff_data[["lookback"]] <- NULL
      
    }else if(checkInput_na$check == FALSE){
      
      showNotification(checkInput_na$warning, type = "error", duration = 3)
      
      reactives$result[["lookback"]] <- NULL
      reactives$payoff_data[["lookback"]] <- NULL
      
    }
    
  })
  
  
  ## progressboxes
  # lookback price progress box
  output$progressBox_lookback_price <- renderValueBox({
    
    create.Progressbox_price(reactives= reactives, "lookback", "optionValue",
                             "blue", "euro")
    
  })
  
  # lookback option sd progress box
  output$progressBox_lookback_sd <- renderValueBox({
    
    create.Progressbox_SD(reactives= reactives, "lookback", "optionSD",
                          "light-blue", "resize-horizontal")
    
  })
  
  # lookback option se progress box
  output$progressBox_lookback_se <- renderValueBox({
    
    create.Progressbox_SE(reactives= reactives, "lookback", "optionSE",
                          "aqua", "exclamation-sign")
    
  })
  
  
  ## histogram of the discounted payoff
  observeEvent(reactives[["warning"]][["lookback"]], {
    
    if(reactives[["warning"]][["lookback"]] == FALSE){
      
      output$lookback_plot <- renderPlot({
        
        create.Histogram(reactives = reactives, "lookback")
        
      },
      height = "auto",width = "auto")
      
    } else if(reactives[["warning"]][["lookback"]] == TRUE){
      
      output$lookback_plot <- NULL
      
    }
    
  })
  
  ## conditional panel for maturity in date or years
  output$maturity_lookback <- renderUI({
    
    create.Conditional.Maturity("lookback")
    
  })
  
  
  #image output----
  output$tutorial <- renderImage({
    filename <- 'images for UI/tutorial.png'
    list(
      src = filename, 
      height = "auto",
      width = "100%"
    )
  }	, deleteFile = FALSE)

  
}  