#### Dependencies ####
library(plFA)
library(dplyr)
library(tidyr)
library(purrr)
library(RcppClock)
ifelse(!dir.exists('output/'), dir.create('output/'), "Folder already exists")

#### Settings setup ####
settings <- list(
  I = list(
    p = 40,
    q = 4,
    n =  c(1000, 2000, 5000)
  ),
  II = list(
    p = 40,
    q = 8,
    n =  c(1000, 2000, 5000)
  )
)
set.seed(1)
for (set in 1:length(settings)) {
  settings[[set]]$cnstr <- build_constrMat(P = settings[[set]]$p, Q = settings[[set]]$q, STRUCT = 'crossed')
  settings[[set]]$loads <- gen_loadings(CONSTRMAT = settings[[set]]$cnstr, LB = .2, UB = .8)
  settings[[set]]$lat_corr <- get_S(THETA = runif(settings[[set]]$q*(settings[[set]]$q-1)/2, -.8, .8), Q = settings[[set]]$q)
  for (j in 1:settings[[set]]$p) {
    if(t(settings[[set]]$loads[j,])%*%settings[[set]]$lat_corr%*%(settings[[set]]$loads[j,])>=1){
      settings[[set]]$loads[j,] <- settings[[set]]$loads[j,]/  as.numeric((t(settings[[set]]$loads[j,])%*%settings[[set]]$lat_corr%*%(settings[[set]]$loads[j,])) + 1e-2)
    }
  }
  settings[[set]]$thr <- c(-1.2, 0, 1.2)
  settings[[set]]$theta <- get_theta(rep(settings[[set]]$thr, settings[[set]]$p),
                                     settings[[set]]$loads, settings[[set]]$lat_corr,
                                     rep(length(settings[[set]]$thr)+1, settings[[set]]$p),
                                     settings[[set]]$cnstr)
  settings[[set]]$true_theta <- getPar(settings[[set]]$theta,
                                       OPTION = 'transformed',
                                       C = settings[[set]]$thr* settings[[set]]$p,
                                       P = settings[[set]]$p,
                                       Q = settings[[set]]$q,
                                       CONSTRMAT = settings[[set]]$cnstr)

}

save(settings, file = 'output/settings.rda')

#### SIMULATION ####
nsim <- 1000

#### setting I ####
#### sto ####
settingLab <- 'I'
setting_grid <- expand_grid(id = 1:nsim,
                   n = settings[[settingLab]]$n,
                   p = settings[[settingLab]]$p,
                   q = settings[[settingLab]]$q)
path <- paste0('output/', settingLab)
ifelse(!dir.exists(path), dir.create(path), "Folder already exists")

lab <- 'st8'


cat(lab,':\n')

est_obj_list <- pbapply::pblapply(purrr::transpose(.l = setting_grid), function(x){


  D <- sim_data(
    SAMPLE_SIZE = x$n ,
    LOADINGS = settings[[settingLab]]$loads,
    THRESHOLDS = settings[[settingLab]]$thr,
    LATENT_COV = settings[[settingLab]]$lat_corr,
    SEED = x$id)

  ppi <- 8
  maxt <- 3*x$n
  burn <- 1*x$n
  my_seq <-  seq(0, 4, by = .5)*x$n
  suppressMessages(
    fit <- fit_plFA(
      DATA = D,
      CONSTR_LIST = list('CONSTRMAT' = settings[[settingLab]]$cnstr, 'CORRFLAG'=1),
      METHOD = 'hyper',
      CONTROL = list(PAIRS_PER_ITERATION = ppi, MAXT = maxt, BURN = burn,  ETA = .01),
      ITERATIONS_SUBSET = my_seq,
      NCORES = 1
    )
  )


  return(fit)
}, cl = 8)


res <- setting_grid
res$lab <- lab
res$fit <- est_obj_list
saveRDS(res, file = paste0(path,'/', lab,'.rds'))
#### num ####
settingLab <- 'I'
setting_grid <- expand_grid(id = 1:nsim,
                            n = settings[[settingLab]]$n,
                            p = settings[[settingLab]]$p,
                            q = settings[[settingLab]]$q)
path <- paste0('output/', settingLab)
ifelse(!dir.exists(path), dir.create(path), "Folder already exists")

lab <- 'num'


cat(lab,':\n')

est_obj_list <- pbapply::pblapply(purrr::transpose(.l = setting_grid), function(x){


  D <- sim_data(
    SAMPLE_SIZE = x$n ,
    LOADINGS = settings[[settingLab]]$loads,
    THRESHOLDS = settings[[settingLab]]$thr,
    LATENT_COV = settings[[settingLab]]$lat_corr,
    SEED = x$id)

  suppressMessages(
    fit <- fit_plFA(
      DATA = D,
      CONSTR_LIST = list('CONSTRMAT' = settings[[settingLab]]$cnstr, 'CORRFLAG'=1),
      METHOD = 'ucminf',
    )
  )


  return(fit)
}, cl = 8)


res <- setting_grid
res$lab <- lab
res$fit <- est_obj_list
saveRDS(res, file = paste0(path,'/', lab,'.rds'))


#### Setting II ####
#### sto ####
settingLab <- 'II'
setting_grid <- expand_grid(id = 1:nsim,
                            n = settings[[settingLab]]$n,
                            p = settings[[settingLab]]$p,
                            q = settings[[settingLab]]$q)
path <- paste0('output/', settingLab)
ifelse(!dir.exists(path), dir.create(path), "Folder already exists")

lab <- 'st8'

cat(lab,':\n')

est_obj_list <- pbapply::pblapply(purrr::transpose(.l = setting_grid), function(x){


  D <- sim_data(
    SAMPLE_SIZE = x$n ,
    LOADINGS = settings[[settingLab]]$loads,
    THRESHOLDS = settings[[settingLab]]$thr,
    LATENT_COV = settings[[settingLab]]$lat_corr,
    SEED = x$id)

  ppi <- 8
  maxt <- 3*x$n
  burn <- 1*x$n
  my_seq <-  seq(0, 4, by = .5)*x$n
  suppressMessages(
    fit <- fit_plFA(
      DATA = D,
      CONSTR_LIST = list('CONSTRMAT' = settings[[settingLab]]$cnstr, 'CORRFLAG'=1),
      METHOD = 'hyper',
      CONTROL = list(PAIRS_PER_ITERATION = ppi, MAXT = maxt, BURN = burn,  ETA = .01),#cpp_ctrl,
      ITERATIONS_SUBSET = my_seq,
      NCORES = 1
    )
  )


  return(fit)
}, cl = 8)

res <- setting_grid
res$lab <- lab
res$fit <- est_obj_list
saveRDS(res, file = paste0(path,'/', lab,'.rds'))

#### num ####
settingLab <- 'II'
setting_grid <- expand_grid(id = 1:nsim,
                            n = settings[[settingLab]]$n,
                            p = settings[[settingLab]]$p,
                            q = settings[[settingLab]]$q)
path <- paste0('output/', settingLab)
ifelse(!dir.exists(path), dir.create(path), "Folder already exists")

lab <- 'num'


cat(lab,':\n')

est_obj_list <- pbapply::pblapply(purrr::transpose(.l = setting_grid), function(x){


  D <- sim_data(
    SAMPLE_SIZE = x$n ,
    LOADINGS = settings[[settingLab]]$loads,
    THRESHOLDS = settings[[settingLab]]$thr,
    LATENT_COV = settings[[settingLab]]$lat_corr,
    SEED = x$id)

  suppressMessages(
    fit <- fit_plFA(
      DATA = D,
      CONSTR_LIST = list('CONSTRMAT' = settings[[settingLab]]$cnstr, 'CORRFLAG'=1),
      METHOD = 'ucminf',
    )
  )


  return(fit)
}, cl = 8)


res <- setting_grid
res$lab <- lab
res$fit <- est_obj_list
saveRDS(res, file = paste0(path,'/', lab,'.rds'))



#### Pointwise performance ####
# rm(list = ls())
#
library(ggplot2)
load("output/settings.rda")
resI <- readRDS("output/I/st8.rds")
resII <- readRDS("output/II/st8.rds")
numI <- readRDS("output/I/num.rds")
numII <- readRDS("output/II/num.rds")


settingLab <- 'I'
nthr <- length(settings[[settingLab]]$thr)*settings[[settingLab]]$p
nload <- sum(settings[[settingLab]]$cnstr)
ncorr <- settings[[settingLab]]$q*(settings[[settingLab]]$q-1)/2
theta <- settings[[settingLab]]$true_theta
mseI <- resI %>%
  mutate(
    path = map(fit, ~getThetaPath(.x))
  ) %>%
  unnest(path) %>%
  mutate(
    thresholds = map_dbl(par, ~mean((.x[1:nthr]-theta[1:nthr])^2)),
    loadings = map_dbl(par, ~mean((.x[(nthr+1):(nthr+nload)]-theta[(nthr+1):(nthr+nload)])^2)),
    correlations = map_dbl(par, ~mean((.x[(nthr+nload+1):(nthr+nload+ncorr)]-theta[(nthr+nload+1):(nthr+nload+ncorr)])^2))
  ) %>%
  gather(key = par_type, value = mse, thresholds, loadings, correlations)
nummseI <- numI |>
  mutate(
    par = map(fit, ~getPar(.x))
  ) %>%
  mutate(
    thresholds = map_dbl(par, ~mean((.x[1:nthr]-theta[1:nthr])^2)),
    loadings = map_dbl(par, ~mean((.x[(nthr+1):(nthr+nload)]-theta[(nthr+1):(nthr+nload)])^2)),
    correlations = map_dbl(par, ~mean((.x[(nthr+nload+1):(nthr+nload+ncorr)]-theta[(nthr+nload+1):(nthr+nload+ncorr)])^2))
  ) %>%
  gather(key = par_type, value = mse, thresholds, loadings, correlations)




settingLab <- 'II'
nthr <- length(settings[[settingLab]]$thr)*settings[[settingLab]]$p
nload <- sum(settings[[settingLab]]$cnstr)
ncorr <- settings[[settingLab]]$q*(settings[[settingLab]]$q-1)/2
theta <- settings[[settingLab]]$true_theta
mseII <- resII %>%
  mutate(
    path = map(fit, ~getThetaPath(.x))
  ) %>%
  filter(n %in% c(1000, 2000, 5000)) |>
  unnest(path) %>%
  mutate(
    thresholds = map_dbl(par, ~mean((.x[1:nthr]-theta[1:nthr])^2)),
    loadings = map_dbl(par, ~mean((.x[(nthr+1):(nthr+nload)]-theta[(nthr+1):(nthr+nload)])^2)),
    correlations = map_dbl(par, ~mean((.x[(nthr+nload+1):(nthr+nload+ncorr)]-theta[(nthr+nload+1):(nthr+nload+ncorr)])^2))
  ) %>%
  gather(key = par_type, value = mse, thresholds, loadings, correlations)
nummseII <- numII |>
  mutate(
    par = map(fit, ~getPar(.x))
  ) %>%
  mutate(
    thresholds = map_dbl(par, ~mean((.x[1:nthr]-theta[1:nthr])^2)),
    loadings = map_dbl(par, ~mean((.x[(nthr+1):(nthr+nload)]-theta[(nthr+1):(nthr+nload)])^2)),
    correlations = map_dbl(par, ~mean((.x[(nthr+nload+1):(nthr+nload+ncorr)]-theta[(nthr+nload+1):(nthr+nload+ncorr)])^2))
  ) %>%
  gather(key = par_type, value = mse, thresholds, loadings, correlations)

mse <- mseI |>
  bind_rows(mseII) |>
  select(-fit,-par) |>
  group_by(n, p, q, lab, iter, par_type) |>
  summarise(mse = mean(mse)) |>
  mutate(
    Iter = factor(iter/n, levels = seq(0,4,by=.25), labels = paste0(seq(0,4,by=.25),'n'), ordered = T),
    q = factor(q, levels = c(4, 8), labels = c('q = 4', 'q = 8'), ordered = T),
    n2 = factor(n, levels = settings[[settingLab]]$n, labels = paste0('n = ', settings[[settingLab]]$n), ordered = T)) |>
  filter(par_type!='thresholds')
nummse <- nummseI |>
  bind_rows(nummseII) |>
  select(-fit,-par) |>
  group_by(n, p, q, lab, par_type) |>
  summarise(mse = mean(mse)) |>
  mutate(
    q = factor(q, levels = c(4, 8), labels = c('q = 4', 'q = 8'), ordered = T),
    n2 = factor(n, levels = settings[[settingLab]]$n, labels = paste0('n = ', settings[[settingLab]]$n), ordered = T)) |>
  filter(par_type!='thresholds')

gg <- mse |>
  ggplot(aes(x = Iter, y = mse, group = par_type)) +
  ggh4x::facet_nested(par_type~q+n2, scales = 'free') +
  scale_y_log10()+
  scale_x_discrete(breaks = paste0(seq(0,4,by=1),'n'), labels = c('0', 'n', '2n', '3n', '4n'))+
  geom_rect( aes(xmin=-Inf,xmax=3,ymin=0,ymax=Inf),
            fill="lightgrey", color =NA)+
  geom_hline(data = nummse, aes(yintercept=mse), linetype = 'dashed', col = 'darkgrey', alpha = .8)+
  theme_bw() +
  geom_line() +
  labs(y = 'Mean square error', x = 'Iterations') +
  theme(legend.position = 'none')
gg
ggsave(gg, file = 'output/mse.pdf', height = 4, width = 10)


#### Time ####
timing <- resI |>
  bind_rows(resII) |>
  bind_rows(numI) |>
  bind_rows(numII) |>
  mutate(time = map_dbl(fit, ~.x@RTime-.x@freqTime)) |>
  select(-fit)

timing |>
  mutate(
    q = factor(q, levels = c(4, 8), labels = c('q = 4', 'q = 8'), ordered = T),
    n2 = factor(n, levels = settings[[1]]$n, labels = paste0('n = ', settings[[1]]$n), ordered = T)) |>
  ggplot(aes(x = lab, y = time)) +
  ggh4x::facet_nested(~q+n2) +
  theme_bw() +
  geom_boxplot()
#### Variability ####
# rm(list = ls())

library(ggplot2)
load("output/settings.rda")
resI <- readRDS("output/I/st8.rds")
resII <- readRDS("output/II/st8.rds")


#### Setting II ####
settingLab <- 'II'
nthr <- length(settings[[settingLab]]$thr)*settings[[settingLab]]$p
nload <- sum(settings[[settingLab]]$cnstr)
ncorr <- settings[[settingLab]]$q*(settings[[settingLab]]$q-1)/2
theta <- settings[[settingLab]]$theta

varII <- resII |>
  mutate(
    path = map(fit, ~getThetaPath(.x, OPTION = 'raw'))
  ) %>%
  unnest(path)
DI <- sim_data(
  SAMPLE_SIZE = 5000,
  LOADINGS = settings[[settingLab]]$loads,
  THRESHOLDS = settings[[settingLab]]$thr,
  LATENT_COV = settings[[settingLab]]$lat_corr,
  SEED = 123)

fr <- pairs_freq(DI, C_VEC = rep(length(settings[[settingLab]]$thr)+1, settings[[settingLab]]$p))
J <- estimate_J(
  Y = DI,
  C_VEC = rep(length(settings[[settingLab]]$thr)+1, settings[[settingLab]]$p),
  A = settings[[settingLab]]$cnstr,
  THETA = theta,
  CORRFLAG = 1
)$est_J
H <- estimate_H(
  C_VEC = rep(length(settings[[settingLab]]$thr)+1, settings[[settingLab]]$p),
  A = settings[[settingLab]]$cnstr,
  THETA = theta,
  FREQ = fr,
  N = nrow(DI),
  CORRFLAG = 1
)$est_H
Hinv <- solve(H)
sum(is.na(J))


est_varII <- pbapply::pblapply(purrr::transpose(.l = varII), function(x){

  n <- x$fit@dims@n
  Tn <- x$iter - x$fit@stoFit@control$BURN
  p <- x$fit@dims@p
  pairs <- p*(p-1)/2
  ppi <- x$fit@stoFit@control$PAIRS_PER_ITERATION

  if(Tn > 0){
    sandwich <- Hinv%*%J%*%Hinv
    asy_var <- diag(sandwich/n)
    a1 <- pairs*(pairs - ppi)/(ppi*(pairs-1))
    a2 <- (pairs-ppi)/(ppi*(pairs-1))
    opt_noise <- diag(Hinv%*%(a1*H - a2*J)%*%Hinv/(n*Tn))




    out <- list( dims = list(p = p, ppi = ppi, pairs = pairs, n = n, Tn = Tn), asy_var = asy_var, opt_noise = opt_noise)
  }else{
    out <- list( dims = list(p = p, ppi = ppi, pairs = pairs, n = n, Tn = Tn), asy_var = NA, opt_noise = NA)
  }


  return(out)
}, cl = 8)
varII$est_var <- est_varII

covII <- varII |>
  mutate(
    est = map2(par, est_var, ~{

      tibble(
        parId = 1:length(.x),
        par = getPar(OBJ = .x,
                     OPTION = 'raw',
                     C = settings[[settingLab]]$thr* settings[[settingLab]]$p,
                     P = settings[[settingLab]]$p,
                     Q = settings[[settingLab]]$q,
                     CONSTRMAT = settings[[settingLab]]$cnstr),
        asy_var = .y$asy_var,
        opt_noise = .y$opt_noise,
        asy_lb = getPar(OBJ = .x - 1.96*sqrt(.y$asy_var),
                        OPTION = 'raw',
                        C = settings[[settingLab]]$thr* settings[[settingLab]]$p,
                        P = settings[[settingLab]]$p,
                        Q = settings[[settingLab]]$q,
                        CONSTRMAT = settings[[settingLab]]$cnstr),
        asy_ub = getPar(OBJ = .x + 1.96*sqrt(.y$asy_var),
                        OPTION = 'raw',
                        C = settings[[settingLab]]$thr* settings[[settingLab]]$p,
                        P = settings[[settingLab]]$p,
                        Q = settings[[settingLab]]$q,
                        CONSTRMAT = settings[[settingLab]]$cnstr),
        fin_lb = getPar(OBJ = .x - 1.96*sqrt(.y$asy_var+.y$opt_noise),
                        OPTION = 'raw',
                        C = settings[[settingLab]]$thr* settings[[settingLab]]$p,
                        P = settings[[settingLab]]$p,
                        Q = settings[[settingLab]]$q,
                        CONSTRMAT = settings[[settingLab]]$cnstr),
        fin_ub = getPar(OBJ = .x + 1.96*sqrt(.y$asy_var+.y$opt_noise),
                        OPTION = 'raw',
                        C = settings[[settingLab]]$thr* settings[[settingLab]]$p,
                        P = settings[[settingLab]]$p,
                        Q = settings[[settingLab]]$q,
                        CONSTRMAT = settings[[settingLab]]$cnstr)
      )
    }

    )) %>%
  select(-fit, -est_var, -par) |>
  unnest(est) %>%
  mutate(
    asy_cov = if_else(settings[[settingLab]]$theta[parId] <= asy_ub & settings[[settingLab]]$theta[parId] >= asy_lb, 1, 0),
    fin_cov = if_else(settings[[settingLab]]$theta[parId] <= fin_ub & settings[[settingLab]]$theta[parId] >= fin_lb, 1, 0)
  )


#### Setting I ####
settingLab <- 'I'
nthr <- length(settings[[settingLab]]$thr)*settings[[settingLab]]$p
nload <- sum(settings[[settingLab]]$cnstr)
ncorr <- settings[[settingLab]]$q*(settings[[settingLab]]$q-1)/2
theta <- settings[[settingLab]]$theta

varI <- resI |>
  mutate(
    path = map(fit, ~getThetaPath(.x, OPTION = 'raw'))
  ) %>%
  unnest(path)

DI <- sim_data(
  SAMPLE_SIZE = 5000,
  LOADINGS = settings[[settingLab]]$loads,
  THRESHOLDS = settings[[settingLab]]$thr,
  LATENT_COV = settings[[settingLab]]$lat_corr,
  SEED = 123)


fr <- pairs_freq(DI, C_VEC = rep(length(settings[[settingLab]]$thr)+1, settings[[settingLab]]$p))
J <- estimate_J(
  Y = DI,
  C_VEC = rep(length(settings[[settingLab]]$thr)+1, settings[[settingLab]]$p),
  A = settings[[settingLab]]$cnstr,
  THETA = theta,
  CORRFLAG = 1
)$est_J
tmp <- estimate_H(
  C_VEC = rep(length(settings[[settingLab]]$thr)+1, settings[[settingLab]]$p),
  A = settings[[settingLab]]$cnstr,
  THETA = theta,
  FREQ = fr,
  N = nrow(DI),
  CORRFLAG = 1
)

H <- tmp$est_H


est_varI <- pbapply::pblapply(purrr::transpose(.l = varI), function(x){

  n <- x$fit@dims@n
  Tn <- x$iter - x$fit@stoFit@control$BURN
  p <- x$fit@dims@p
  pairs <- p*(p-1)/2
  ppi <- x$fit@stoFit@control$PAIRS_PER_ITERATION

  if(Tn > 0){
    sandwich <- Hinv%*%J%*%Hinv
    asy_var <- diag(sandwich/n)
    a1 <- pairs*(pairs - ppi)/(ppi*(pairs-1))
    a2 <- (pairs-ppi)/(ppi*(pairs-1))
    opt_noise <- diag(Hinv%*%(a1*H - a2*J)%*%Hinv/(n*Tn))




    out <- list( dims = list(p = p, ppi = ppi, pairs = pairs, n = n, Tn = Tn), asy_var = asy_var, opt_noise = opt_noise)
  }else{
    out <- list( dims = list(p = p, ppi = ppi, pairs = pairs, n = n, Tn = Tn), asy_var = NA, opt_noise = NA)
  }

  return(out)
}, cl = 8)
varI$est_var <- est_varI

covI <- varI |>
  mutate(
    est = map2(par, est_var, ~{

      tibble(
        parId = 1:length(.x),
        par = getPar(OBJ = .x,
                     OPTION = 'raw',
                     C = settings[[settingLab]]$thr* settings[[settingLab]]$p,
                     P = settings[[settingLab]]$p,
                     Q = settings[[settingLab]]$q,
                     CONSTRMAT = settings[[settingLab]]$cnstr),
        asy_var = .y$asy_var,
        opt_noise = .y$opt_noise,
        asy_lb = getPar(OBJ = .x - 1.96*sqrt(.y$asy_var),
                        OPTION = 'raw',
                        C = settings[[settingLab]]$thr* settings[[settingLab]]$p,
                        P = settings[[settingLab]]$p,
                        Q = settings[[settingLab]]$q,
                        CONSTRMAT = settings[[settingLab]]$cnstr),
        asy_ub = getPar(OBJ = .x + 1.96*sqrt(.y$asy_var),
                        OPTION = 'raw',
                        C = settings[[settingLab]]$thr* settings[[settingLab]]$p,
                        P = settings[[settingLab]]$p,
                        Q = settings[[settingLab]]$q,
                        CONSTRMAT = settings[[settingLab]]$cnstr),
        fin_lb = getPar(OBJ = .x - 1.96*sqrt(.y$asy_var+.y$opt_noise),
                        OPTION = 'raw',
                        C = settings[[settingLab]]$thr* settings[[settingLab]]$p,
                        P = settings[[settingLab]]$p,
                        Q = settings[[settingLab]]$q,
                        CONSTRMAT = settings[[settingLab]]$cnstr),
        fin_ub = getPar(OBJ = .x + 1.96*sqrt(.y$asy_var+.y$opt_noise),
                        OPTION = 'raw',
                        C = settings[[settingLab]]$thr* settings[[settingLab]]$p,
                        P = settings[[settingLab]]$p,
                        Q = settings[[settingLab]]$q,
                        CONSTRMAT = settings[[settingLab]]$cnstr)
      )
    }

    )) %>%
  select(-fit, -est_var, -par) |>
  unnest(est) %>%
  mutate(
    asy_cov = if_else(settings[[settingLab]]$theta[parId] <= asy_ub & settings[[settingLab]]$theta[parId] >= asy_lb, 1, 0),
    fin_cov = if_else(settings[[settingLab]]$theta[parId] <= fin_ub & settings[[settingLab]]$theta[parId] >= fin_lb, 1, 0)
  )


#### TOT ####
emp_cov <- covI |>
  bind_rows(covII) |>
  select(-id) |>
  group_by(n, p, q, lab, parId, iter) |>
  summarise_all(mean)

summary(emp_cov)


gg <- emp_cov |>
  filter(!is.na(asy_cov)) |>
  filter(n %in% c(1000, 2000, 5000)) |>
  mutate(par_type = 'thresholds',
         par_type = if_else(parId>nthr, 'loadings', par_type),
         par_type = if_else(parId>nthr+nload, 'reparametrised correlations', par_type)) |>
  filter(par_type!='thresholds') |>
  select(n, p, q, lab, parId, par_type, iter, asymptotic = asy_cov, corrected = fin_cov) |>
  gather(key = 'var', val = 'cov', asymptotic, corrected) |>
  mutate(Iter = factor(iter/n, levels = seq(0,5,by=.25), labels = paste0(seq(0,5,by=.25),'n'), ordered = T),
         n = factor(n, levels = settings[[settingLab]]$n, labels = paste0('n = ', settings[[settingLab]]$n), ordered = T),
         q = factor(q, levels = c(4, 8), labels = c('q = 4', 'q = 8'), ordered = T)
  ) |>
  ggplot(aes(y = cov, x = Iter, group = interaction(iter, var))) +
  ggh4x::facet_nested(par_type~q+n, scales = 'free')+
  scale_y_continuous(breaks = c(.5, .75, .9, .94, .95,.96, 1))+
  labs(y = 'Empirical coverage', x = 'Iterations', col = '') +
  theme_bw() +
  geom_hline(yintercept = .95, linetype = 'dashed') +
  geom_boxplot(aes(col = var)) +
  scale_color_grey(start = .7, end = .3) +
  theme(legend.position = 'top')
gg
ggsave(gg, file = 'output/gg_coverage.pdf', height = 6, width = 10)

#### settings plots####
loadsI <- reshape::melt(settings$I$loads)
loadsI$lab <- 'q = 4'
loadsII <- reshape::melt(settings$II$loads)
loadsII$lab <- 'q = 8'
gg_load <- loadsI |>
  bind_rows(loadsII) |>
  ggplot(aes(x = X2, y = X1, fill = value)) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = 'right') +
  facet_grid(~lab, scales = 'free') +
  labs(x = 'Latent factors', y = 'Items', fill = ' ') +
  geom_tile() +
  lims(fill = c(0,1)) +
  scale_fill_gradient(low = 'white', high = 'black', limits = c(0,1.05))
ggsave(gg_load, file = 'output/true_loadings.pdf', width = 10, height = 5)


corrI <- reshape::melt(settings$I$lat_corr)
corrI$lab <- 'q = 4'
corrII <- reshape::melt(settings$II$lat_corr)
corrII$lab <- 'q = 8'
gg_corr <- corrI |>
  bind_rows(corrII) |>
  mutate(X1 = factor(X1, levels = 8:1, ordered = T)) |>
  mutate(sign = if_else(value >= 0, 'positive', 'negative')) |>
  ggplot(aes(x = X2, y = X1, col = sign)) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = 'right') +
  facet_wrap(~lab, scales = 'free') +
  labs(x = '', y = '', size = '', col = '') +
  geom_point(aes(size = abs(value), col = sign)) +
  scale_color_grey(start = .7, end = .3)

ggsave(gg_corr, file = 'output/true_corr.pdf', width = 6, height = 3)
