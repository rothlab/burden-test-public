### library
library(logger)

# Handle INFO message
printInfo = function(message) {
  log_info(message, "\n")
}

# Handle ERROR message
throwError = function(message) {
  stop(message)
}

# Setup logger
log_appender(appender_console)
log_errors()
