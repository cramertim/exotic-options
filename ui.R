## load packages
library(shiny)
library(shinydashboard)
library(waiter)
library(quantmod)
library(gridExtra)
library(shinyjs)
library(shinycssloaders)
library(timeDate)
library(quantmod)

# Source functions ----
source("functions.R")

# UI content ----
# header
header <-  dashboardHeader(title = "MC option pricing")


# Sidebar content
sidebar <-  dashboardSidebar(
  
  sidebarMenu(
    
    id = "tab",
    
    shinyFeedback::useShinyFeedback(),
    
    waiter::use_waitress(),
    
    menuItem("Vanilla options",
             #Tabs für Vanilla Option
             menuSubItem("European", tabName = "european", selected = TRUE),
             menuSubItem("American", tabName = "american"),
             startExpanded = TRUE
    ),
    
    #Tabs für Exotic Option
    menuItem("Exotic options",
             #Path-dependent option
             menuItem("Path-dependent options",
                      menuSubItem("Asian", tabName = "asian"),
                      menuSubItem("Barrier", tabName = "barrier"),
                      menuSubItem("Lookback", tabName = "lookback")),
             
             #Time-dependent Option
             menuItem("Time-dependent options",
                      menuSubItem("Forward Start", tabName = "forward"),
                      menuSubItem("Chooser", tabName = "chooser")),
             
             #Other Exotic Option
             menuItem("Other exotic options",
                      menuSubItem("Binary", tabName = "binary"),
                      menuSubItem("Gap", tabName = "gap"),
                      menuSubItem("Compound", tabName = "compound")
             )
    ),
    
    menuItem("About",
             tabName = "about"),
    
    #Choose if dividend yield is used
    checkboxInput(
      inputId = "dividend_YesOrNo",
      label = "Asset with dividend",
      value = FALSE
    ),
    
    #Choose model for Simulation
    selectInput(inputId = "simulationModel",
                label = "Simulation Model",
                choices = c("Black-Scholes","Heston"),
                selected = "Black-Scholes"),
    
    
    ## input of stockprice
    #@default = 0
    numericInput(inputId = "stockPrice",  
                 label = "Stock price",
                 value = 50
    ),
    ## volatility of the underlying asset
    # @default = 0
    numericInput(inputId = "volatility",
                 label = tags$span(
                   "Volatility ", 
                   tags$i(
                     class = "glyphicon glyphicon-info-sign", 
                     style = "color:#0072B2;",
                     title = "Volatility of the underlying asset at time zero (in percent)"
                   )
                 ),
                 value = 20)
    
    
  ),
  
  uiOutput("heston_active"),
  
  sidebarMenu(
    
    #input of risk free interest rate
    #@default = 3 month US Treasury
    numericInput(inputId = "riskFreeRate", 
                 label =  tags$span(
                   "Risk free rate ", 
                   tags$i(
                     class = "glyphicon glyphicon-info-sign", 
                     style = "color:#0072B2;",
                     title = "latest 10 Year Treasury Rate (in percent)"
                   )
                 ), 
                 #current risk free rate 
                 value = get.riskFreeRate()
    ),
    
    
    uiOutput("yield_chosen"),
    
    uiOutput("timesteps"),
    
    uiOutput("skalar_timesteps"),
    
    ## Anzahl der Ziehungen
    sliderInput(
      inputId = "sampleSize",
      label = "Sample size in 1.000",
      min = 0, max = 100,
      value = 1, step = 2
    )
    #uiOutput("sampleSize")
    
  )
  
)


# Body content
body <-  dashboardBody(
  tabItems(
    
    ## Start Tab--------------------------------------------------------------------     
    tabItem(tabName = "about",
            
            tabsetPanel(
              type = "tabs",
              tabPanel("Theory",
                       
                       fluidRow(
                         
                         infoBox(title = strong("About page"),
                                 subtitle = "This page provides an informational overview about the app, 
                                   the options and the used calculation formulas. This tab contains the theory behind the calculation 
                                   of the option prices.", 
                                 width = 12)
                         
                       ),
                       
                       fluidRow(
                         
                         box(
                           title = NULL, 
                           width = 12, 
                           solidHeader = FALSE, 
                           status = "info", 
                           height = "3000",
                           withMathJax(),
                           tags$script("MathJax.Hub.Config({tex2jax: {inlineMath: [['$','$'],
                                ['\\(','\\)']],processEscapes: true}});"),
                           
                           h1("Theoretical background of the calculation"),
                           h4("This section provides an overview and explanation of the procedure and methods used to calculate the price of the exotic option."),
                           p("First, stock paths are generated based on provided input parameters using either the ",strong("Black-Scholes"), " or the ", strong("Heston"),"model. 
                               Based on the simulated paths, the payoff is calculated using the specific ",strong("formula for each option"),", which can be found ",strong("next to")," the 
                               ",strong("plot of each option"),". After calculating the payoff for each paths, the average is used to value the option as a whole. This process is 
                               part of the ",strong("Monte-Carlo Simulation"), "and explained in more detail below. To reduce variance of the payoff parameters, the", strong("antithetic 
                               variable technique")," is implemented."),
                           
                           h3("Black-Scholes stock model"),
                           
                           p("The Black-Scholes model is a mathematical model for pricing options. It is based on the
assumption of constant volatility for the underlying asset and the existence of a risk-free
rate of return. The model uses five key inputs: the current price of the underlying asset,
the exercise price of the option, the time to expiration, the risk-free interest rate, and the
volatility of the underlying asset’s price. We use the following notation"),
                           withMathJax(helpText(
                             '$\\mu$ average growth rate of the underlying stock' )),
                           withMathJax(helpText(
                             '$W$ Wiener process $\\text{(}$one-dimensional Brownian motion$\\text{)}$')),
                           withMathJax(helpText(
                             '$\\sigma$ volatility' )),
                           withMathJax(helpText(
                             '$Z$ standard normal random variable' )),
                           p("The Black Scholes model assumes that the 
                               stock price evolves according to the following stochastic differential equation:"),
                           withMathJax("$$dS_t = \\mu S_tdt + \\sigma S_tdW_t$$"),
                           p("This implies that we can simulate the stock price at time t+1 using the formula:"),
                           withMathJax("$$S_t = S_0 exp ((\\mu-\\frac{1}{2}\\sigma^2)t  + \\sigma\\sqrt{t}Z)$$"),
                           
                           h3("Heston stock model"),
                           p("The Heston model is an extension of the Black-Scholes model and allows the volatility 
                               of the underlying to vary over time. The model uses nine key inputs: the current price of the underlying asset,
the exercise price of the option, the time to expiration, the risk-free interest rate, the
volatility of the underlying asset’s price, the long run average variance, the rate of reversion in the variance, 
the volatility of volatility, and the correlation. We use the following notation:"),
                           withMathJax(helpText(
                             '$\\mu$ average growth rate' )),
                           withMathJax(helpText(
                             '$\\theta$ long run average variance' )),
                           withMathJax(helpText(
                             '$\\kappa$ rate of mean reversion in the variance' )),
                           withMathJax(helpText(
                             '$\\xi$ volatility of volatility' )),
                           withMathJax(helpText(
                             '$r$ risk free rate' )),
                           withMathJax(helpText(
                             '$dW_t^S$ and $dW_t^{\\upsilon}$ are Wiener processes with instantaneous correlation $\\rho$' )),
                           p("The Heston model assumes that the stock price and its variance evolve according to the following
                               two stochastic differential equations:"),
                           withMathJax("$$dS_t = \\mu S_tdt + \\sqrt{\\upsilon_t}S_tdW_t^S$$"),
                           withMathJax("$$d\\upsilon_t = \\kappa(\\theta - \\upsilon_t)dt + \\xi\\sqrt{\\upsilon_t}dW_t^\\upsilon$$"),
                           p("This implies that we can simulate the stock price at time t using the formula:"),
                           withMathJax("$$S_{t+i} = S_te^{(r-\\frac{\\upsilon_t}{2})dt}+\\sqrt{\\upsilon_i}dW_t^S$$"),
                           withMathJax("$$\\upsilon_{t+1} = \\upsilon_t + \\kappa(\\theta - \\upsilon_t)dt + \\xi\\sqrt{\\upsilon_t}dW_t^\\upsilon$$"),
                           
                           h3("Option pricing with Monte-Carlo Simulation"),
                           p("Monte Carlo simulation is a method of valuing options that uses random sampling to simulate
the possible future prices of the underlying asset. The simulation is run many times, with
each iteration producing a different potential outcome. The average of these outcomes is
used to estimate the option’s value."),
                           p("1. Using a risk-neutral probability distribution to simulate S"),
                           p("2. Calculate the payoff of the derivative"),
                           p("3. Repeat steps 1 and 2 to get multiple sample values of the derivative's payoff."),
                           p("4. Calculate the mean of the sample payoffs to get an estimate of the expected payoff"),
                           p("5. Discount this expected payoff at the risk-free rate to get an estimate of the value of the
derivative."),
                           
                           
                           
                           
                           h3("Variance reduction"),
                           p("The approach we use to decrease the speed of convergence of the Monte Carlo method is the antithetic variable technique. 
The key idea of this approach is that every simulation trial involves calculating two values of the payoff of the derivative. 
The first value (f_1) is calculated in the usual way; the second value (f_2) is calculated by switching the sign of all the random 
samples from standard normal distributions. The sample value of the derivative calculated from a simulation trial
is the average of (f_1) and (f_2)
"),        
                           p("Denote (f_{ave}) as the average of (f_1) and (f_2)"),
                           withMathJax("$$f_{ave} = \\frac{f_1 + f_2}{2}$$")
                           
                         )
                         
                       )
                       
                       
              ),
              
              tabPanel("Tutorial",
                       
                       fluidRow(
                         
                         infoBox(title = strong("About page"),
                                 subtitle = "This page provides an informational overview about the app, 
                                   the options and the used calculation formulas. This tab contains an tutorial for the usage of the app.", 
                                 width = 12),
                         box(
                           title = NULL, 
                           status = "info",
                           width = 12, 
                           height = 2500,
                           solidHeader = FALSE, 
                           imageOutput('tutorial')
                         )
                       )
                       
              )
              
              
            )
            
    ),
    
    ## Heston Tab--------------------------------
    tabItem(tabName = "hestonTab",
            fluidRow(
              column(width = 4,
                     
                     box(
                       title = "Adjust Heston model parameters", 
                       width = NULL, 
                       solidHeader = TRUE, 
                       status = "warning", 
                       height = 610,
                       
                       h4("Parameters of the volatility"),
                       p("The volatility is a stochastic, mean-reverting process."),
                       wellPanel(style = "background: white",
                                 
                                 ## input of kappa
                                 #@default = 0
                                 numericInput(inputId = "kappa",
                                              label = tags$span(
                                                "Rate of mean reversion in the variance", 
                                                tags$i(
                                                  class = "glyphicon glyphicon-info-sign", 
                                                  style = "color:#0072B2;",
                                                  title = "Kappa"
                                                )
                                              ),
                                              value = 3.8),
                                 
                                 ## input of epsilon
                                 #@default = 0
                                 numericInput(inputId = "epsilon",
                                              label = tags$span(
                                                "Volatility of volatility", 
                                                tags$i(
                                                  class = "glyphicon glyphicon-info-sign", 
                                                  style = "color:#0072B2;",
                                                  title = "Epsilon"
                                                )
                                              ),
                                              value = 0.9288),
                                 
                                 ## input of theta
                                 #@default = 0
                                 numericInput(inputId = "theta",
                                              label = tags$span(
                                                "Long run average variance", 
                                                tags$i(
                                                  class = "glyphicon glyphicon-info-sign", 
                                                  style = "color:#0072B2;",
                                                  title = "Theta"
                                                )
                                              ),
                                              value = 0.04)
                                 
                                 
                       ),
                       
                       h4("Correlation"),
                       wellPanel(style = "background: white",
                                 
                                 ## input of rho
                                 #@default = 0
                                 numericInput(inputId = "rho",
                                              label = tags$span(
                                                "Correlation between volatility and prices", 
                                                tags$i(
                                                  class = "glyphicon glyphicon-info-sign", 
                                                  style = "color:#0072B2;",
                                                  title = "Rho: instantaneous correlation of the two Wiener processes"
                                                )
                                              ),
                                              value = 0)
                                 
                       )
                       
                       
                     ),
                     
                     box(title = "Set parameters of the heston pricepath plot", width = NULL, solidHeader = TRUE, status = "primary", height = 200,
                         
                         ## input of expiration date
                         #@default = today
                         dateInput(inputId = "hestonPathDate",
                                   label = "Maturity",
                                   value = (Sys.Date() + 100),
                                   daysofweekdisabled = c(0,6),
                                   datesdisabled = create.Holidays()),
                         
                         ## ActionButton to start simulation
                         actionButton(inputId = "hestonPlot",
                                      label = strong("Calculate pricepath"),
                                      width = "100%",
                                      style = 'background: MediumAquamarine; height:40px; font-size:110%')
                     )
                     
                     
                     
              ),
              
              
              box(title = "Heston pricepath plot", width = 8, solidHeader = TRUE, status = "primary", height = 830,
                  
                  textOutput("hestonplotText"),
                  
                  withSpinner(plotOutput("path_plot", height = 670))
                  
              )   
            ),
            fluidRow(
              mainPanel(
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br()
              )
            )
    ),
    
    
    ## American Tab-----------------------------------------------------------------
    tabItem(tabName = "american",
            fluidRow(
              infoBox(title = strong("American option"),
                      subtitle = "An American option is a financial derivative that gives the holder the right to buy or sell an underlying asset at a predetermined price at any time before the expiration date.
                        This type of option is more flexible than a European option, and is therefore often traded at a small premium in comparison to an equivalent European option.",
                      width = 12)
            ),
            
            fluidRow(
              
              box(
                title = "Option specific input parameters",
                width = 4,
                solidHeader = FALSE,
                status = "primary",
                height = 470,
                
                wellPanel(style = "background: white",
                          splitLayout(
                            cellWidths = c("40%","60%"),
                            
                            #choose wether call or put
                            radioButtons(
                              inputId = "optionType_american",
                              label = "Type",
                              choices = c("Call","Put"),
                              inline = FALSE
                            ),
                            
                            ## input of strike Price
                            #@default = 0
                            numericInput(inputId = "strikePrice_american",
                                         label = "Strike price",
                                         value = 50))), br(),
                
                wellPanel(style = "background: white",
                          splitLayout(
                            cellWidths = c("40%","60%"),
                            radioButtons("american_dateOrYear",
                                         label = "Maturity as",
                                         c("Date" = "date",
                                           "Years" = "years")),
                            uiOutput("maturity_american"))), br(),
                
                ## ActionButton to start simulation
                actionButton(inputId = "americanSimulation",
                             label = strong("Start MC simulation"),
                             width = "100%",
                             style = 'background: MediumAquamarine; height:40px; font-size:110%')
              ),
              tabBox(
                title = "",
                side = "left",
                width = 8,
                height = 470,
                tabPanel(
                  title = tags$span(
                    "Payoff histogram", 
                    tags$i(
                      class = "glyphicon glyphicon-info-sign", 
                      style = "color:#0072B2;",
                      title = "Histogram considering only the positive paths of variance reduction"
                    )
                  ), width = 8, solidHeader = FALSE, status = "primary", height = 550,
                  htmlOutput("warning_american"),
                  
                  withSpinner(plotOutput("american_plot", height = 400))
                )
                ,
                tabPanel(title = "Formula",
                         h3("American option calculation"),
                         p("Implementation of Least-Square method proposed by Longstaff & Schwart (2001)"),
                         tags$a(href="https://www.researchgate.net/publication/5216848_Valuing_American_Options_by_Simulation_A_Simple_Least-Squares_Approach",
                                "Click here to view")
                         
                )
              )
              
              
            ),
            
            fluidRow(
              (valueBoxOutput("progressBox_american_price")),
              (valueBoxOutput("progressBox_american_sd")),
              (valueBoxOutput("progressBox_american_se"))
            ),
            fluidRow(
              mainPanel(
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br()
              )
            )
    ),
    
    
    
    ## Asian Tab--------------------------------------------------------------------
    tabItem(tabName = "asian",
            fluidRow(
              infoBox(title = strong("Asian option"),
                      subtitle = "An Asian option is a type of exotic option, where the payoff is determined 
                      by the average price of the underlying asset during a certain period of time. They are 
                      commonly used in the commodity markets and are used as a way to hedge against price volatility.", 
                      width = 12)
            ),
            
            fluidRow(
              tags$head(
                tags$style(HTML('.shiny-split-layout>div {overflow: hidden;}')),
              ),
              
              box(
                title = "Option specific input parameters", 
                width = 4, 
                solidHeader = FALSE, 
                status = "info", 
                height = 590,
                
                
                selectInput(inputId = "subType_asian",
                            label = "Method of calculation",
                            choices = c("Arithmetic"="arithmetic", "Geometric"="geometric")
                ),
                
                wellPanel(style = "background: white",
                          
                          verticalLayout(
                            
                            splitLayout(
                              cellWidths = c("40%","60%"),
                              
                              #choose wether call or put
                              radioButtons(
                                inputId = "optionType_asian",
                                label = "Type",
                                choices = c("Call","Put"),
                                inline = FALSE
                              ),
                              
                              #choose whether call or put
                              radioButtons(
                                inputId = "observation_asian",
                                label = "Range",
                                choices = c("Total","Interval",
                                            "Period"),
                                inline = FALSE
                              )
                              
                            ),
                            
                            uiOutput("date_range_asian")
                            
                          )
                          
                ),
                
                wellPanel(style = "background: white",
                          splitLayout(
                            cellWidths = c("40%","60%"),
                            radioButtons("asian_dateOrYear",
                                         label = "Maturity as",
                                         c("Date" = "date",
                                           "Years" = "years"),),
                            uiOutput("maturity_asian"))), 
                
                ## ActionButton to start simulation
                actionButton(inputId = "asianSimulation",
                             label = strong("Start MC simulation"),
                             width = "100%",
                             style = 'background: MediumAquamarine; height:40px; font-size:110%')
              ),
              
              tabBox(
                title = "",
                side = "left",
                width = 8, 
                height = 590,
                tabPanel(title = tags$span(
                  "Payoff histogram", 
                  tags$i(
                    class = "glyphicon glyphicon-info-sign", 
                    style = "color:#0072B2;",
                    title = "Histogram considering only the positive paths of variance reduction"
                  )
                ), width = 8, solidHeader = FALSE, status = "primary", height = 580,
                htmlOutput("warning_asian"),
                
                withSpinner(plotOutput("asian_plot", height = 520))),
                tabPanel("Formula",
                         h3("Asian option calculation"),
                         withMathJax(),
                         tags$script("MathJax.Hub.Config({tex2jax: {inlineMath: [['$','$'],
                                ['\\(','\\)']],processEscapes: true}});"),
                         p("(S_{ave}) average price across observations"),
                         p("(K) strike price"),
                         withMathJax("$$\\text{Arithmetic average } S_{ave} = \\frac{\\sum_{i=0}^NSi}{N}$$"),
                         withMathJax("$$\\text{Geometric average } S_{ave} = \\sqrt{\\prod_{i} S_i}$$"),
                         withMathJax("$$\\text{payoff of call option} = max(S_{ave} − K, 0)$$"),
                         withMathJax("$$\\text{payoff of put option} = max(K − S_{ave}, 0)$$"))
              )
              
            ),
            
            fluidRow(
              (valueBoxOutput("progressBox_asian_price")),
              (valueBoxOutput("progressBox_asian_sd")),
              (valueBoxOutput("progressBox_asian_se"))
            ),
            fluidRow(
              mainPanel(
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br()
              )
            )
    ),
    
    ## Barrier Tab------------------------------------------------------------------
    tabItem(tabName = "barrier",
            fluidRow(
              infoBox(title = strong("Barrier option"),
                      subtitle = "Barrier options are a type of financial derivative that gives the holder the right to buy or sell an underlying asset at a predetermined price 
                    if a certain price level is reached or breached. They can be either knock-in, which means the option becomes active if the barrier is 
                    breached, or knock-out, which means the option becomes invalid if the barrier is breached.", 
                      width = 12)
            ),
            
            fluidRow(
              
              box(
                title = "Option specific input parameters", 
                width = 4, 
                solidHeader = FALSE, 
                status = "info", 
                height = 490,
                
                ## Choose barrier type
                # @default = Up-And-In
                selectInput(inputId = "subType_barrier",
                            label = "Type of barrier option",
                            choices = c("Up-And-In", "Up-And-Out", "Down-And-In",
                                        "Down-And-Out")
                ),
                
                
                wellPanel(style = "background: white",
                          splitLayout(
                            cellWidths = c("33%","33%","33%"),
                            
                            #choose wether call or put
                            radioButtons(
                              inputId = "optionType_barrier",
                              label = "Type",
                              choices = c("Call","Put"),
                              inline = FALSE
                            ),
                            
                            ## input of strike price
                            #@default = 0
                            numericInput(inputId = "strikePrice_barrier",
                                         label = "Strike price",
                                         value = 50),
                            
                            #input of barrier
                            #@default = 0
                            numericInput(inputId = "barrierValue",
                                         label = "Barrier",
                                         value = 55)
                          )
                ), 
                
                wellPanel(style = "background: white",
                          splitLayout(
                            cellWidths = c("40%","60%"),
                            radioButtons("barrier_dateOrYear",
                                         label = "Maturity as",
                                         c("Date" = "date",
                                           "Years" = "years")),
                            uiOutput("maturity_barrier"))), 
                
                ## ActionButton to start simulation
                actionButton(inputId = "barrierSimulation",
                             label = strong("Start MC simulation"),
                             width = "100%",
                             style = 'background: MediumAquamarine; height:40px; font-size:110%')
              ),
              
              tabBox(
                title = "",
                side = "left",
                width = 8, 
                height = 490,
                tabPanel(title = tags$span(
                  "Payoff histogram", 
                  tags$i(
                    class = "glyphicon glyphicon-info-sign", 
                    style = "color:#0072B2;",
                    title = "Histogram considering only the positive paths of variance reduction"
                  )
                ), width = 8, solidHeader = FALSE, status = "primary", height = 520,
                htmlOutput("warning_barrier"),
                
                withSpinner(plotOutput("barrier_plot", height = 420))),
                tabPanel("Formula",
                         h3("Barrier option calculation"),
                         withMathJax(),
                         tags$script("MathJax.Hub.Config({tex2jax: {inlineMath: [['$','$'],
                                ['\\(','\\)']],processEscapes: true}});"),
                         p("(B) barrier"),
                         p("(S) price of underlying asset at maturity"),
                         p("(K) strike price"),
                         p("(S_{min/max}) minimal/maximal (S_t) in 0 ⩽ t ⩽ T"),
                         splitLayout(
                           cellWidths = c("50%","50%"),
                           box(title =  "Payoff of call option", width = 12,
                               withMathJax("$$\\text{Down-And-Out: if } (B<S_{min}) \\text{ payoff } = max(S-K)\\text{, else 0}$$"),
                               withMathJax("$$\\text{Down-And-In: if } (B>S_{min}) \\text{ payoff } = max(S-K)\\text{, else 0}$$"),
                               withMathJax("$$\\text{Up-And-Out: if } (B>S_{max}) \\text{ payoff } = max(S-K)\\text{, else 0}$$"),
                               withMathJax("$$\\text{Up-And-In: if } (B<S_{max}) \\text{ payoff } = max(S-K)\\text{, else 0}$$"),),
                           box(title =  "Payoff of put option", width = 12,
                               withMathJax("$$\\text{Down-And-Out: if } (B<S_{min}) \\text{ payoff } = max(K-S)\\text{, else 0}$$"),
                               withMathJax("$$\\text{Down-And-In: if } (B>S_{min}) \\text{ payoff } = max(K-S)\\text{, else 0}$$"),
                               withMathJax("$$\\text{Up-And-Out: if } (B>S_{max}) \\text{ payoff } = max(K-S)\\text{, else 0}$$"),
                               withMathJax("$$\\text{Up-And-In: if } (B<S_{max}) \\text{ payoff } = max(K-S)\\text{, else 0}$$"),)
                         )
                         
                )
              )
            ),
            
            fluidRow(
              (valueBoxOutput("progressBox_barrier_price")),
              (valueBoxOutput("progressBox_barrier_sd")),
              (valueBoxOutput("progressBox_barrier_se"))
            ),
            fluidRow(
              mainPanel(
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br()
              )
            )
    ),
    ## Binary Tab------------------------------------------------------------------
    tabItem(tabName = "binary",
            fluidRow(
              infoBox(title = strong("Binary option"),
                      subtitle = "Binary options are a type of financial derivative that allow traders to speculate on whether the
                    price of an underlying asset will be above or below a certain level at a specific time in the future.
                    They may have a fixed payout or are dependent on the price at expiration date.. A cash-or-nothing binary option pay the difference between strike price and price at expiration,
                      whereas asset-or-nothing options pay the underlying asset if they are in in the money.", 

                      width = 12)
            ),
            
            fluidRow(
              
              box(
                title = "Option specific input parameters", 
                width = 4, 
                solidHeader = FALSE, 
                status = "info", 
                height = 490,
                
                ## choose binary type
                # @default = Up-And-In
                selectInput(inputId = "subType_binary",
                            label = "Type of binary option",
                            choices = c("Asset-or-nothing", "Cash-or-nothing")
                ),
                
                wellPanel(style = "background: white",
                          splitLayout(
                            cellWidths = c("40%","60%"),
                            
                            #choose wether call or put
                            radioButtons(
                              inputId = "optionType_binary",
                              label = "Type",
                              choices = c("Call","Put"),
                              inline = FALSE
                            ),
                            
                            uiOutput("cash_or_nothing_input")
                            
                          )
                ),
                
                wellPanel(style = "background: white",
                          splitLayout(
                            cellWidths = c("40%","60%"),
                            radioButtons("binary_dateOrYear",
                                         label = "Maturity as",
                                         c("Date" = "date",
                                           "Years" = "years")),
                            uiOutput("maturity_binary"))), 
                
                ## ActionButton to start simulation
                actionButton(inputId = "binarySimulation",
                             label = strong("Start MC simulation"),
                             width = "100%",
                             style = 'background: MediumAquamarine; height:40px; font-size:110%')
              ),
              
              tabBox(
                title = "",
                side = "left",
                width = 8, 
                height = 490,
                tabPanel(
                  title = tags$span(
                    "Payoff histogram", 
                    tags$i(
                      class = "glyphicon glyphicon-info-sign", 
                      style = "color:#0072B2;",
                      title = "Histogram considering only the positive paths of variance reduction"
                    )
                  ), width = 8, solidHeader = FALSE, status = "primary", height = 540,
                  htmlOutput("warning_binary"),
                  
                  withSpinner(plotOutput("binary_plot", height = 420))),
                tabPanel("Formula",
                         h3("Binary option calculation"),
                         withMathJax(),
                         tags$script("MathJax.Hub.Config({tex2jax: {inlineMath: [['$','$'],
                                ['\\(','\\)']],processEscapes: true}});"),
                         p("(S) price of underlying asset at maturity"),
                         p("(K) strike price"),
                         p("(Q) fixed amount payoff"),
                         br(),
                         strong("Cash-or-nothing "),
                         withMathJax("$$\\text{payoff of call option: if } (S>K) \\text{ payoff } = Q\\text{, else 0}$$"),
                         withMathJax("$$\\text{payoff of put option: if } (S<K) \\text{ payoff } = Q\\text{, else 0}$$"),
                         br(),
                         strong("Asset-or-nothing "),
                         withMathJax("$$\\text{payoff of call option: if } (S>K) \\text{ payoff } = S\\text{, else 0}$$"),
                         withMathJax("$$\\text{payoff of put option: if } (S<K) \\text{ payoff } = S\\text{, else 0}$$")
                )
              )
            ),
            
            fluidRow(
              (valueBoxOutput("progressBox_binary_price")),
              (valueBoxOutput("progressBox_binary_sd")),
              (valueBoxOutput("progressBox_binary_se"))
            ),
            fluidRow(
              mainPanel(
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br()
              )
            )
    ),
    ## Chooser Tab-----------------------------------------------------------------
    tabItem(tabName = "chooser",
            fluidRow(
              infoBox(title = strong("Chooser option"),
                      subtitle = "A chooser option is a financial derivative that gives the holder the right to choose whether the option 
                    is a call or a put at a specific time in the future. This decides wether the holder has the right to buy the underlying asset at a 
                    predetermined price, or the right to sell the underlying asset at a predetermined price.", 
                      width = 12)
            ),
            
            fluidRow(
              
              box(
                title = "Option specific input parameters", 
                width = 4, 
                solidHeader = FALSE, 
                status = "info", 
                height = 580,
                
                splitLayout(
                  cellWidths = c("0%","50%", "50%"),
                  tags$head(tags$style(HTML("
                              .shiny-split-layout > div {
                                overflow: visible;
                              }
                              "))),
                  
                  ## choose chooser type
                  selectInput(inputId = "chooser_option",
                              label = "Type of chooser option",
                              choices = c("Asian" = "asian", 
                                          "Barrier" = "barrier",
                                          "Binary" = "binary",
                                          "European" = "european",
                                          "Forward Start" = "forwardStart",
                                          "Gap" = "gap",
                                          "Lookback" = "lookback"
                              )
                  ),
                  ## choose subType
                  uiOutput("optionPanel_chooser_subtype")
                ),
                
                
                
                
                wellPanel(style = "background: white",
                          splitLayout(
                            cellWidths = c("50%","50%"),
                            
                            ## input of decisionDate
                            #@default = today + 7
                            dateInput(inputId = "chooserDecisionDate",
                                      label = "Decision time",
                                      value = (Sys.Date() + 6),
                                      daysofweekdisabled = c(0,6),
                                      datesdisabled = create.Holidays()
                            ),
                            
                            ##input of strikePrice
                            #@default = 0
                            numericInput(inputId = "strikePrice_chooser",
                                         label = "Strike price",
                                         value = 50)),
                          
                          
                          
                          uiOutput("optionPanel_chooser_optionvalues")
                          
                ), 
                
                
                wellPanel(style = "background: white",
                          splitLayout(
                            cellWidths = c("40%","60%"),
                            radioButtons("chooser_dateOrYear",
                                         label = "Maturity as",
                                         c("Date" = "date",
                                           "Years" = "years")),
                            uiOutput("maturity_chooser"))), 
                
                ## ActionButton to start simulation
                actionButton(inputId = "chooserSimulation",
                             label = strong("Start MC simulation"),
                             width = "100%",
                             style = 'background: MediumAquamarine; height:40px; font-size:110%')
                
              ),
              
              tabBox(
                title = "",
                side = "left",
                width = 8, 
                height = 580,
                tabPanel(
                  title = tags$span(
                    "Payoff histogram", 
                    tags$i(
                      class = "glyphicon glyphicon-info-sign", 
                      style = "color:#0072B2;",
                      title = "Histogram considering only the positive paths of variance reduction"
                    )
                  ), width = 8, solidHeader = FALSE, status = "primary", height = 550,
                  
                  htmlOutput("warning_chooser"),
                  
                  withSpinner(plotOutput("chooser_plot", height = 510))),
                tabPanel("Formula",
                         h3("Chooser option calculation"),
                         withMathJax(),
                         tags$script("MathJax.Hub.Config({tex2jax: {inlineMath: [['$','$'],
                                ['\\(','\\)']],processEscapes: true}});"),
                         p("(C) call option value"),
                         p("(P) put option value"),
                         p("(T_1) decision time $\\text{(}$in the future$\\text{)}$"),
                         withMathJax("$$\\text{payoff at }T_1 = max(C,P)$$")
                )
              )
            ),
            
            fluidRow(
              (valueBoxOutput("progressBox_chooser_price")),
              (valueBoxOutput("progressBox_chooser_sd")),
              (valueBoxOutput("progressBox_chooser_se"))
            ),
            fluidRow(
              mainPanel(
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br()
              )
            )
    ),
    ## Compound Tab-----------------------------------------------------------------
    tabItem(tabName = "compound",
            fluidRow(
              infoBox(title = strong("Compound option"),
                      subtitle = "A compound option is an option for which its underlying security is another option. Therefore, there are two strike prices and two exercise dates. They are available for any combination of calls and puts. For example, a put where the underlying is a call option or a call where the underlying is a put option.", 
                      width = 12)
            ),
            
            fluidRow(
              tabBox(
                title = "",
                side = "left",
                width = 4, 
                height = 580,
                
                tabPanel("Overlying option",
                         
                         strong("Attention: The overlying option is always a standard European option."),
                         br(),
                         br(),
                         
                         wellPanel(style = "background: white",
                                   
                                   
                                   splitLayout(
                                     cellWidths = c("40%","60%"),
                                     
                                     #choose wether call or put
                                     radioButtons(
                                       inputId = "optionType_compound",
                                       label = "Type",
                                       choices = c("Call on call" = "CoC",
                                                   "Call on put" = "CoP",
                                                   "Put on call" = "PoC",
                                                   "Put on put" = "PoP"),
                                       inline = FALSE
                                     ),
                                     
                                     verticalLayout(
                                       
                                       #@default = 0
                                       numericInput(inputId = "strikePrice_compound_first",
                                                    label = "Strike price",
                                                    value = 50),
                                       
                                       #@default = 0
                                       numericInput(inputId = "sampleSize_compound",
                                                    label =  tags$span(
                                                      "Sample size ", 
                                                      tags$i(
                                                        class = "glyphicon glyphicon-info-sign", 
                                                        style = "color:#0072B2;",
                                                        title = "Sample size of the overlying option. 
                                                        Sample size in the sidebar is not used for this option."
                                                      )
                                                    ),
                                                    value = 1000, min = 1, max = 10000)
                                       
                                       
                                     )
                                     )), 
                         br(),
                         
                         wellPanel(style = "background: white",
                                   splitLayout(
                                     cellWidths = c("40%","60%"),
                                     radioButtons("compound_first_dateOrYear",
                                                  label = "Maturity as",
                                                  c("Date" = "date",
                                                    "Years" = "years")),
                                     uiOutput("maturity_compound_first")))),
                
                tabPanel("Underlying option", 
                         
                         
                         splitLayout(
                           cellWidths = c("0%","50%", "50%"),
                           tags$head(tags$style(HTML("
                              .shiny-split-layout > div {
                                overflow: visible;
                              }
                              "))),
                           
                           ## choose compound type
                           selectInput(inputId = "compound_option_second",
                                       label = "Type of underlying option",
                                       choices = c("Asian" = "asian", 
                                                   "Barrier" = "barrier",
                                                   "Binary" = "binary",
                                                   "European" = "european",
                                                   "Forward Start" = "forwardStart",
                                                   "Gap" = "gap",
                                                   "Lookback" = "lookback"
                                       )
                           ),
                           ## choose subType
                           uiOutput("optionPanel_compound_subtype")
                         ),
                         
                         
                         
                         wellPanel(style = "background: white",
                                   
                                   splitLayout(
                                     cellWidths = c("50%", "50%"),
                                     
                                     ## input of strike Price second
                                     #@default = 0
                                     numericInput(inputId = "strikePrice_compound_second",
                                                  label = "Strike price",
                                                  value = 50),
                                     
                                     ## input of second sample Size
                                     #@default = 0
                                     numericInput(inputId = "sampleSize_compound_second",
                                                  label = "Sample size",
                                                  value = 100, min = 1, max = 10000)
                                   ),
                                   
                                   uiOutput("optionPanel_compound_optionvalues")
                                   
                         ), 
                         
                         
                         wellPanel(style = "background: white",
                                   splitLayout(
                                     cellWidths = c("40%","60%"),
                                     radioButtons("compound_second_dateOrYear",
                                                  label = "Maturity as",
                                                  c("Date" = "date",
                                                    "Years" = "years")),
                                     uiOutput("maturity_compound_second"))),
                         
                         ## ActionButton um Simulation zu starten
                         actionButton(inputId = "compoundSimulation",
                                      label = strong("Start MC simulation"),
                                      width = "100%",
                                      style = 'background: MediumAquamarine; height:40px; font-size:110%')
                         
                )
                
                
              ),
              
              tabBox(
                title = "",
                side = "left",
                width = 8, 
                height = 580,
                tabPanel(title = tags$span(
                  "Payoff histogram", 
                  tags$i(
                    class = "glyphicon glyphicon-info-sign", 
                    style = "color:#0072B2;",
                    title = "Histogram considering only the positive paths of variance reduction"
                  )
                ), width = 8, solidHeader = FALSE, status = "primary", height = 550,
                
                htmlOutput("warning_compound"),
                
                withSpinner(plotOutput("compound_plot", height = 510))
                ),
                tabPanel("Formula",
                         h3("Compound option calculation"),
                         withMathJax(),
                         tags$script("MathJax.Hub.Config({tex2jax: {inlineMath: [['$','$'],
                                ['\\(','\\)']],processEscapes: true}});"),
                         p("(S) price of underlying asset at maturity"),
                         p("(K_1) first strike price"),
                         p("(K_2) second strike strike price"),
                         p("(T_1) first expiration date"),
                         p("(T_2) second expiration date"),
                         p("(c)$\\text{()}$ call option value"),
                         p("(p)$\\text{()}$ put option value"),
                         withMathJax("$$\\text{price for first option} = c(S, T_1 , K_1) \\text{ or } p(S, T_1 , K_1)$$"),
                         withMathJax("$$\\text{price at $T_1$} = max[ c(S, T_2 - T_1; K_2) , p(S, T_2 - T_1; K_2)]$$")
                )
              )
            ),
            
            fluidRow(
              (valueBoxOutput("progressBox_compound_price")),
              (valueBoxOutput("progressBox_compound_sd")),
              (valueBoxOutput("progressBox_compound_se"))
            ),
            fluidRow(
              mainPanel(
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br()
              )
            )
    ),
    
    ## European Tab-----------------------------------------------------------------
    tabItem(tabName = "european",
            fluidRow(
              infoBox(title = strong("European option"),
                      subtitle = "A European option is a version of an options contract that gives the holder the right to buy or sell an underlying asset at a predetermined price on a certain date in the future. 
                             An investor would not be able to exercise the option early and take delivery of or sell the shares. 
                             Instead, the call or put action will only take place on the date of option maturity.", 
                      width = 12)
            ),
            
            fluidRow(
              
              box(
                title = "Option specific input parameters", 
                width = 4, 
                solidHeader = FALSE, 
                status = "info", 
                height = 470,
                
                wellPanel(style = "background: white",
                          splitLayout(
                            cellWidths = c("40%","60%"),
                            
                            #choose whether call or put
                            radioButtons(
                              inputId = "optionType_european",
                              label = "Type",
                              choices = c("Call","Put"),
                              inline = FALSE
                            ),
                            
                            ## input of strikePrice
                            #@default = 0
                            numericInput(inputId = "strikePrice_european",
                                         label = "Strike price",
                                         value = 50))), br(),
                
                wellPanel(style = "background: white",
                          splitLayout(
                            cellWidths = c("40%","60%"),
                            radioButtons("european_dateOrYear",
                                         label = "Maturity as",
                                         c("Date" = "date",
                                           "Years" = "years")),
                            uiOutput("maturity_european"))), br(),
                
                actionButton(inputId = "europeanSimulation",
                             label = strong("Start MC simulation"),
                             width = "100%",
                             style = 'background: MediumAquamarine; height:40px; font-size:110%')
              ),
              tabBox(
                title = "",
                side = "left",
                width = 8, 
                height = 470,
                tabPanel(title = tags$span(
                  "Payoff histogram", 
                  tags$i(
                    class = "glyphicon glyphicon-info-sign", 
                    style = "color:#0072B2;",
                    title = "Histogram considering only the positive paths of variance reduction"
                  )
                ), width = 8, solidHeader = FALSE, status = "primary", height = 450,
                
                htmlOutput("warning_european"),
                
                withSpinner(plotOutput("european_plot", height = 400))),
                tabPanel("Formula",
                         h3("European option calculation"),
                         withMathJax(),
                         tags$script("MathJax.Hub.Config({tex2jax: {inlineMath: [['$','$'],
                                ['\\(','\\)']],processEscapes: true}});"),
                         p("(S) price of underlying asset at maturity"),
                         p("(K) strike price"),
                         withMathJax("$$\\text{payoff of call option} = max(S − K, 0)$$"),
                         withMathJax("$$\\text{payoff of put option} = max(K − S, 0)$$")
                ),
              )
            ),
            
            fluidRow(
              (valueBoxOutput("progressBox_european_price")),
              (valueBoxOutput("progressBox_european_sd")),
              (valueBoxOutput("progressBox_european_se"))
            ),
            fluidRow(
              mainPanel(
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br()
              )
            )
    ),
    
    ## Forward Start Tab-----------------------------------------------------------------
    tabItem(tabName = "forward",
            fluidRow(
              infoBox(title = strong("Forward start option"),
                      subtitle = "A forward start option is a financial derivative that gives the holder the right to buy or sell an underlying 
                    asset at a predetermined price at a specific time in the future, known as the 'start date.' The start date for 
                    a forward start option is later than the date the option is purchased, and the option has an expiration date 
                    that is set at the time of purchase.", 
                      width = 12)
            ),
            
            fluidRow(
              
              box(
                title = "Option specific input parameters", 
                width = 4, 
                solidHeader = FALSE, 
                status = "info", 
                height = 480,
                
                ## input of start date
                #@default = today + 7
                dateInput(inputId = "forwardStartDate",
                          label = "Start date",
                          value = (Sys.Date() + 6),
                          daysofweekdisabled = c(0,6),
                          datesdisabled = create.Holidays()
                ),
                
                wellPanel(style = "background: white",
                          splitLayout(
                            cellWidths = c("40%","60%"),
                            
                            #choose wether call or put
                            radioButtons(
                              inputId = "optionType_forwardStart",
                              label = "Type",
                              choices = c("Call","Put"),
                              inline = FALSE
                            ),
                            
                            #@default = 0
                            numericInput(inputId = "strikePrice_forwardStart",
                                         label = "Strike price",
                                         value = 50))), 
                
                wellPanel(style = "background: white",
                          splitLayout(
                            cellWidths = c("40%","60%"),
                            radioButtons("forwardStart_dateOrYear",
                                         label = "Maturity as",
                                         c("Date" = "date",
                                           "Years" = "years")),
                            uiOutput("maturity_forwardStart"))), 
                
                actionButton(inputId = "forwardStartSimulation",
                             label = strong("Start MC simulation"),
                             width = "100%",
                             style = 'background: MediumAquamarine; height:40px; font-size:110%')
              ),
              
              tabBox(
                title = "",
                side = "left",
                width = 8, 
                height = 480,
                tabPanel(title = tags$span(
                  "Payoff histogram", 
                  tags$i(
                    class = "glyphicon glyphicon-info-sign", 
                    style = "color:#0072B2;",
                    title = "Histogram considering only the positive paths of variance reduction"
                  )
                ), width = 8, solidHeader = FALSE, status = "primary", height = 470,
                
                htmlOutput("warning_forwardStart"),
                
                withSpinner(plotOutput("forwardStart_plot", height = 410))
                ),
                tabPanel("Formula",
                         h3("Forward start option calculation"),
                         withMathJax(),
                         tags$script("MathJax.Hub.Config({tex2jax: {inlineMath: [['$','$'],
                                ['\\(','\\)']],processEscapes: true}});"),
                         p("$\\widehat{E}$ expected value under the risk-neutral measure"),
                         p("(T_1) start date $\\text{(}$in the future$\\text{)}$"),
                         p("(c) value at time zero of an at-the-money option"),
                         p("(S_0) price of underlying at time zero"),
                         p("(S_1) price of underlying at time (T_1)"),
                         withMathJax("$$\\text{price} = e^{-rT_1}\\widehat{E} [c\\frac{S_0}{S_1}]$$")
                )
              )
            ),
            
            fluidRow(
              (valueBoxOutput("progressBox_forwardStart_price")),
              (valueBoxOutput("progressBox_forwardStart_sd")),
              (valueBoxOutput("progressBox_forwardStart_se"))
            ),
            fluidRow(
              mainPanel(
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br()
              )
            )
            
    ),
    ## Gap Tab----------------------------------------------------------------------
    tabItem(tabName = "gap",
            fluidRow(
              infoBox(title = strong("Gap Option"),
                      subtitle = "Gap options are a type of financial derivative that gives the holder the right to 
                    buy or sell an underlying asset at a predetermined price if the price of the underlying 
                    asset 'gaps' beyond a certain level. Gap options are typically used to hedge 
                    against the risk of large price movements in the underlying asset.", 
                      width = 12)
            ),
            
            fluidRow(
              
              box(
                title = "Option specific input parameters", 
                width = 4, 
                solidHeader = FALSE, 
                status = "info", 
                height = 450,
                
                wellPanel(style = "background: white",
                          splitLayout(
                            cellWidths = c("33%","33%","33%"),
                            
                            #choose wether call or put
                            radioButtons(
                              inputId = "optionType_gap",
                              label = "Type",
                              choices = c("Call","Put"),
                              inline = FALSE
                            ),
                            
                            #@default = 0
                            numericInput(inputId = "strikePrice_gap1",
                                         label = "Strike price",
                                         value = 40),
                            #@default = 0
                            numericInput(inputId = "strikePrice_gap2",
                                         label = "Trigger price",
                                         value = 50)
                          )
                ), br(),
                
                wellPanel(style = "background: white",
                          splitLayout(
                            cellWidths = c("40%","60%"),
                            radioButtons("gap_dateOrYear",
                                         label = "Maturity as",
                                         c("Date" = "date",
                                           "Years" = "years")),
                            uiOutput("maturity_gap")
                          )
                ), br(),
                
                actionButton(inputId = "gapSimulation",
                             label = strong("Start MC simulation"),
                             width = "100%",
                             style = 'background: MediumAquamarine; height:40px; font-size:110%')
              ),
              tabBox(
                title = "",
                side = "left",
                width = 8, 
                height = 450,
                tabPanel(
                  title = tags$span(
                    "Payoff histogram", 
                    tags$i(
                      class = "glyphicon glyphicon-info-sign", 
                      style = "color:#0072B2;",
                      title = "Histogram considering only the positive paths of variance reduction"
                    )
                  ), width = 8, solidHeader = FALSE, status = "primary", height = 450,
                  
                  htmlOutput("warning_gap"),
                  
                  withSpinner(plotOutput("gap_plot", height = 380))
                ),
                tabPanel("Formula",
                         h3("Gap option calculation"),
                         withMathJax(),
                         tags$script("MathJax.Hub.Config({tex2jax: {inlineMath: [['$','$'],
                                ['\\(','\\)']],processEscapes: true}});"),
                         p("(S) price of underlying asset at maturity"),
                         p("(K_1) strike price"),
                         p("(K_2) trigger price"),
                         splitLayout(
                           cellWidths = c("50%","50%"),
                           box(
                             title = tags$span(
                             "Payoff of call option", 
                             tags$i(
                               class = "glyphicon glyphicon-info-sign", 
                               style = "color:#0072B2;",
                               title = "if K_1 = K_2, the gap option payoff will be the same as that of a european option"
                             )
                           ), width = 12,
                           p(strong("if "),"((K_2) > (K_1))",strong("then")),
                           p(HTML('&nbsp;'),HTML('&nbsp;'),"(S) > (K_2): payoff = (S) - (K_1)"),
                           p(HTML('&nbsp;'),HTML('&nbsp;'),"(S) < (K_2): payoff = 0"),
                           p(strong("else if "),"( (K_2) < (K_1))",strong("then")),
                           p(HTML('&nbsp;'),HTML('&nbsp;'),"payoff = (K_2) - (K_1)")),
                           box(title = tags$span(
                             "Payoff of put option", 
                             tags$i(
                               class = "glyphicon glyphicon-info-sign", 
                               style = "color:#0072B2;",
                               title = "if K_1 = K_2, the gap option payoff will be the same as that of a european option."
                             )
                           ), width = 12,
                           p(strong("if "),"((K_2) < (K_1))",strong("then")),
                           p(HTML('&nbsp;'),HTML('&nbsp;'),"(S) < (K_2): payoff = (K_1) - (S)"),
                           p(HTML('&nbsp;'),HTML('&nbsp;'),"(S) > (K_2): payoff = 0"),
                           p(strong("else if "),"( (K_2) > (K_1))",strong("then")),
                           p(HTML('&nbsp;'),HTML('&nbsp;'),"payoff = (K_1) - (K_2)"))
                         )
                )
              )
            ),
            
            fluidRow(
              (valueBoxOutput("progressBox_gap_price")),
              (valueBoxOutput("progressBox_gap_sd")),
              (valueBoxOutput("progressBox_gap_se"))
            ),
            fluidRow(
              mainPanel(
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br()
              )
            )
            
    ),
    ## Lookback Tab-----------------------------------------------------------------
    tabItem(tabName = "lookback",
            fluidRow(
              infoBox(title = strong("Lookback option"),
                      subtitle = "Lookback options come in two types: 'floating strike' options, in which the strike price is set at the expiration 
                    date based on the maximum or minimum price achieved by the underlying asset, and 'fixed strike' options, 
                    in which the strike price is set at the time the option is purchased.", 
                      width = 12)
            ),
            
            fluidRow(
              
              box(
                title = "Option specific input parameters", 
                width = 4, 
                solidHeader = FALSE, 
                status = "info", 
                height = 490,
                
                # @default = Up-And-In
                selectInput(inputId = "subType_lookback",
                            label = "Type of lookback option",
                            choices = c("Fixed", "Floating")
                ),
                
                wellPanel(style = "background: white",
                          splitLayout(
                            cellWidths = c("40%","60%"),
                            
                            #choose wether call or put
                            radioButtons(
                              inputId = "optionType_lookback",
                              label = "Type",
                              choices = c("Call","Put"),
                              inline = FALSE
                            ),
                            
                            ## numerische Eingabe des Ausübungspreises
                            #@default = 0
                            numericInput(inputId = "strikePrice_lookback",
                                         label = "Strike price",
                                         value = 50))), 
                
                wellPanel(style = "background: white",
                          splitLayout(
                            cellWidths = c("40%","60%"),
                            radioButtons("lookback_dateOrYear",
                                         label = "Maturity as",
                                         c("Date" = "date",
                                           "Years" = "years")),
                            uiOutput("maturity_lookback"))), 
                
                ## ActionButton um Simulation zu starten
                actionButton(inputId = "lookbackSimulation",
                             label = strong("Start MC simulation"),
                             width = "100%",
                             style = 'background: MediumAquamarine; height:40px; font-size:110%')
              ),
              
              tabBox(
                title = "",
                side = "left",
                width = 8, 
                height = 490,
                tabPanel(
                  title = tags$span(
                    "Payoff histogram", 
                    tags$i(
                      class = "glyphicon glyphicon-info-sign", 
                      style = "color:#0072B2;",
                      title = "Histogram considering only the positive paths of variance reduction"
                    )
                  ), width = 8, solidHeader = FALSE, status = "primary", height = 470,
                  htmlOutput("warning_lookback"),
                  
                  withSpinner(plotOutput("lookback_plot", height = 420))
                ),
                tabPanel("Formula",
                         h3("Lookback option calculation"),
                         withMathJax(),
                         tags$script("MathJax.Hub.Config({tex2jax: {inlineMath: [['$','$'],
                                ['\\(','\\)']],processEscapes: true}});"),
                         p("(S) price of underlying asset at maturity"),
                         p("(K) strike price"),
                         p("(S_{min/max}) minimum/maximum (S_t) in 0 ⩽ t ⩽ T"),
                         br(),
                         strong("Fixed lookback"),
                         withMathJax("$$\\text{payoff of call option} = max (S_{max} − K, 0)$$"),
                         withMathJax("$$\\text{payoff of put option} = max (K − S_{min}, 0)$$"),
                         br(),
                         strong("Floating lookback"),
                         withMathJax("$$\\text{payoff of call option} = max (S − S_{min}, 0)$$"),
                         withMathJax("$$\\text{payoff of put option} = max (S_{max} − S, 0)$$")
                )
              )
            ),
            
            fluidRow(
              (valueBoxOutput("progressBox_lookback_price")),
              (valueBoxOutput("progressBox_lookback_sd")),
              (valueBoxOutput("progressBox_lookback_se"))
            ),
            fluidRow(
              mainPanel(
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br()
              )
            )
    )
    
    
    
  )
)

#kombinieren
ui <- dashboardPage(header, sidebar, body, skin = "black")

