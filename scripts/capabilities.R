#### Dependencies ####
library(plFA)
library(tidyverse)
library(RcppClock)
library(latex2exp)
ifelse(!dir.exists('output/'), dir.create('output/'), "Folder already exists")

#### import data ####
dataRL<-read.table("data/capabilities.TXT", header=T)
names(dataRL) <- c( "cc1","cc2","cc3","cc4","cc5",

                    "ic1","ic2","ic3",

                    "rt1","rt2","rt3","rt4","rt5",

                    "eu1","eu2","eu3","eu4","eu5")
dataRL
summary(dataRL)
newData <- as.matrix(dataRL-1)
summary(newData)
p <- ncol(newData)
K = p*(p-1)/2; K
n <- nrow(newData)
A <- matrix(0, p, 4)
A[1:5, 1] <- A[6:8, 2] <- A[9:13, 3] <- A[14:18, 4] <-  1
A

RLnumFit <- fit_plFA(
  DATA = newData,
  CONSTR_LIST = list('CONSTRMAT' = A, 'CORRFLAG' = 1),
  METHOD = 'ucminf'
)

getPar(RLnumFit, 'list')
numSe <- computeVar(RLnumFit, DATA = newData, NUMDERIV = T, OPTION = 'transformed')$asymptotic_variance |> sqrt()
numSe[109:length(numSe)]

ppi <- 8
maxt <- 5000
burn <- 500
my_seq <-  seq(0, 1, by = .1)*maxt
RLstoFit <- fit_plFA(
  DATA = newData,
  CONSTR_LIST = list('CONSTRMAT' = A, 'CORRFLAG' = 1),
  METHOD = 'hyper',
  CONTROL = list(PAIRS_PER_ITERATION = ppi,
                 MAXT = maxt, BURN = burn,
                 ETA = .005, EACHCLOCK = 250,
                 PAR2 = .001, TOL = 1e-3, TOLCOUNT = 1, EACHCHECK = 500, CHECKCONV = 1),
  ITERATIONS_SUBSET = my_seq,
  NCORES = 1
)
RLstoFit
RLstoFit@stoFit@pathValNll
getPar(RLstoFit, 'list')

stoVa <- computeVar(RLstoFit, DATA = newData, NUMDERIV = T, OPTION = 'transformed')
stoSe <- sqrt(stoVa$asymptotic_variance+stoVa$optimisation_noise)

tibble(num = getPar(RLnumFit, 'transformed')[109:length(numSe)],
       num_se = round(numSe[109:length(numSe)],3),
       sto = getPar(RLstoFit, 'transformed')[109:length(numSe)],
       sto_se = round(stoSe[109:length(stoSe)],3)) |> print(n=100)



