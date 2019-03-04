#Load the libraries and install them if needed
if (!require(devtools)) install.packages("devtools")
if (!require(foreign))  install.packages("foreign")
if (!require(doRNG))    install.packages("doRNG")
if (!require(xtable))   install.packages("xtable")
if (!require(Hmisc))    install.packages("Hmisc")
install_github("bmelly/discreteQ", force = TRUE)
library(discreteQ)
library(foreign)
library(doParallel)
library(doRNG)
library(xtable)
library(Hmisc)

#Define the directory containing the data and codes
#(leave as it as "." if the currentworking directory contains the files)
code.dir <- "."
#This directory must contain a folder 'Data' containing both data files
#and a folder 'Codes' containing the .R files.
if (!dir.exists(paste0(code.dir, "/Data")))
  stop("Data folder not found")
if (!dir.exists(paste0(code.dir, "/Codes")))
  stop("Codes folder not found")
if (!dir.exists(paste0(code.dir, "/Results")))
  dir.create(paste0(code.dir, "/Results"))

#Create a parallel socket cluster
#(or leave it as NULL if you don't want it)
cl <- NULL

#Figure 1: histograms
#Indicative computing time: less than 1 minute
source(paste0(code.dir, "/Codes/figure1.R"))

#Figure 2 and 3: construction of the bands
#Indicative computing time: less than 1 minute
source(paste0(code.dir, "/Codes/figure2and3.R"))

#Figure 4: Insurance coverage and health care utilization
#Indicative computing time: 4 hours
source(paste0(code.dir, "/Codes/figure4.R"))

#Figure 5: Racial differences in mental ability of young children at 8 months
#Indicative computing time: 3 days
source(paste0(code.dir, "/Codes/figure5.R"))

#Figure 6: Racial differences in mental ability of young children at 7 years
#Indicative computing time: 10 days
source(paste0(code.dir, "/Codes/figure6.R"))

#Simulations (Tables 1 to 4)
#Indicative computing time: 1 day
source(paste0(code.dir, "/Codes/simulations.R"))
