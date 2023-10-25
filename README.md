# exoticOptions

App for the valuation of exotic options using the Monte Carlo simulation in different model environments.

Online version: 
https://exoticoptions.shinyapps.io/exoticoptions/?_ga=2.136252210.2064526663.1675284358-1995696272.1673430467<br>
Made with "shinyapps.io". Tutorial here: https://shiny.rstudio.com/articles/shinyapps.html


## List of the used packages
Packages **beyond the default packages** are listed here. Those packages are **mandatory** for execution
of the program.

1. shiny
2. shinydashboard
3. waiter
4. quantmod
5. gridExtra
6. shinyjs
7. shinycssloaders
8. timeDate
9. ggplot2
10. matrixStats
11. MASS

## Code
This is a small overview of the used funcitons inside the files.

### functions.R
**Functions for calculation of the option prices** <br>
* firstValueRow <br>
*Returns matrix with only the first value of each row replacing all other values by 0. Used after iterative process for American options.*
* simulate.Heston <br>
*Heston simulation which creates stock paths.*
* create.StockPaths <br>
*Generate random stock paths with either Black-Scholes or Heston model.*
* is.convertible.to.date <br>
*Checks if handed over string is convertible to an date object.*
* create.Holidays <br>
*Excludes bank holidays from calendar.*
* price.Option <br>
*Price the passed option.*

**Functions not directly for calculation of the option prices** <br>
* set.active.option <br>
*Function that detects, which option is calculated and stores this value, so that the functions "check.Input_na" and "check.Input" know, which option is active.*
* check.Heston_na <br>
*Check if the paramters used for the heston path plot are NA.*
* check.Heston <br>
*Check if the paramters used for the heston path plot are correct.*
* check.Input_na <br>
*Check if the option specific parameters are not NA.*
* check.Input <br>
*Check if the option specific parameters are correct.*
* create.warning <br>
*Checks if the option specific calculation is up to date. Works because of the reactive behavior. The function is checked each time a input or reactive variable changes.*
* choose.maturity.Input <br>
*Detects if the selected maturity is an "date" or "years" input.*
* get.riskFreeRate <br>
*Calls the current risk free rate from the web.*
* create.Progressbox_price <br>
*Creates a progressbox for the price of the option.*
* create.Progressbox_SD <br>
*Creates a progressbox for the standard deviation of the option.*
* create.Progressbox_SE <br>
*Creates a progressbox for the standard error of the option.*
* create.Conditional.Maturity <br>
*Creates the panel with the maturity as either date or year.*
* create.Histogram <br>
*Creates an histogram of the discounted payoff with the package ggplot2.*


### Server.R
All functions and techniques are used pretty straightforward and consistent to the shiny tutorials. Additionally, the comments help to understand what happens in the according lines.

### ui.R
The UI is built with the additional package "shinydashboard". A tutorial can be found in the web.

