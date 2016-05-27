#!/usr/bin/env Rscript

library(devtools)
library(testthat)
load_all(".")
test_package("fgsea")
