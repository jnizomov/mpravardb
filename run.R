library(shiny)

runApp(
  appDir = getwd(),
  launch.browser = FALSE,
  workerId = "",
  quiet = FALSE,
  display.mode = c("auto", "normal", "showcase"),
  test.mode = getOption("shiny.testmode", FALSE),
  host = "0.0.0.0",
  port = 19401
)