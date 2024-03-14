#script to install R packages
#setwd("~/Desktop/demo/basicR")

library(readxl)

#1.using CRAN
install.packages('readxl')

#install packages from github(public repositories)
#https://github.com/rstudio/shiny
install.packages('remotes')
library(remotes)
remotes::install_github('rstudio/shiny')
library(shiny)

#install from bioconductor
#we need a package called BiocManager which is available on CRAN
install.packages('BiocManager')
library(BiocManager)
BiocManager::install('Biostrings')
library(Biostrings)

dnastring=DNAString('AATCGAACTGG')
reverseComplement(dnastring)

#where are my packages saved?
.libPaths()
