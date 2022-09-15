library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(tools)
library(reshape2)
library(odbc)
library(RMariaDB)
library(yaml)

# Load database (DB) parameters
dbParams = read_yaml("../variant-effect-predictor-assessment_450k/db_connect.yaml")

# Load phenotypes

# Connect to DB
con <- dbConnect(
  drv = MariaDB(), 
  username = dbParams$username,
  password = dbParams$password, 
  host = dbParams$host,
  port = dbParams$port,
  dbname = dbParams$dbname
)
