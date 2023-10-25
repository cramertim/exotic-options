## note: the working directory is automatically the folder which contains the app.R file

# source server ----
source("server.R")

# Source ui ----
source("ui.R")

# start app ----
shinyApp(ui = ui, server = server)