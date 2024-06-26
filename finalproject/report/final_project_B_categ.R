library(dplyr)

source("pair_plot.R")

dube <- read.csv("proteomics_dube_2023.csv")
n <- nrow(dube)
p <- ncol(dube)
