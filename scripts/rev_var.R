# Proper asymptotic variances computed at the stochastic estimator
# Asymptotic covariance matrices in sims.R are evaluated at the true value
# of the parameter vector


library(plFA)
library(dplyr)
library(tidyr)
library(purrr)
ifelse(!dir.exists('output/'), dir.create('output/'), "Folder already exists")
library(ggplot2)
load("output/settings.rda")
resII <- readRDS("output/II/st8.rds")


#### Setting II ####
settingLab <- 'II'
nthr <- length(settings[[settingLab]]$thr)*settings[[settingLab]]$p
nload <- sum(settings[[settingLab]]$cnstr)
ncorr <- settings[[settingLab]]$q*(settings[[settingLab]]$q-1)/2
theta <- settings[[settingLab]]$theta

varII_1000 <- resII %>%
  mutate(
    path = map(fit, ~getThetaPath(.x, OPTION = 'raw'))
  ) %>%
  unnest(path) %>%
  filter(id %in% 1:2, n == 1000)

est_varII <- pbapply::pblapply(purrr::transpose(.l = varII_1000), function(x){

  n <- x$fit@dims@n
  Tn <- x$iter - x$fit@stoFit@control$BURN
  p <- x$fit@dims@p
  pairs <- p*(p-1)/2
  ppi <- x$fit@stoFit@control$PAIRS_PER_ITERATION
  freq <- x$fit@freq
  bartheta <- as.numeric(x$fit@stoFit@pathAvTheta[which(x$fit@stoFit@trajSubset==x$iter)[1], ])



  if(Tn > 0){
    D <- sim_data(
      SAMPLE_SIZE = x$n ,
      LOADINGS = settings[[settingLab]]$loads,
      THRESHOLDS = settings[[settingLab]]$thr,
      LATENT_COV = settings[[settingLab]]$lat_corr,
      SEED = x$id)
    H <- estimate_H(
      C_VEC = rep(length(settings[[settingLab]]$thr)+1, settings[[settingLab]]$p),
      A = settings[[settingLab]]$cnstr,
      THETA = bartheta,
      FREQ = freq,
      N = nrow(D),
      CORRFLAG = 1
    )$est_H
    Hinv <- solve(H)
    J <- estimate_J(
      Y = D,
      C_VEC = rep(length(settings[[settingLab]]$thr)+1, settings[[settingLab]]$p),
      A = settings[[settingLab]]$cnstr,
      THETA = bartheta,
      CORRFLAG = 1
    )$est_J

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
varII_1000$est_var <- est_varII
save(varII_1000, file = 'output/varII_1000.rds')

varII_2000 <- resII %>%
  mutate(
    path = map(fit, ~getThetaPath(.x, OPTION = 'raw'))
  ) %>%
  unnest(path) %>%
  filter(id %in% 1:2, n == 2000)
est_varII <- pbapply::pblapply(purrr::transpose(.l = varII_2000), function(x){

  n <- x$fit@dims@n
  Tn <- x$iter - x$fit@stoFit@control$BURN
  p <- x$fit@dims@p
  pairs <- p*(p-1)/2
  ppi <- x$fit@stoFit@control$PAIRS_PER_ITERATION
  freq <- x$fit@freq
  bartheta <- as.numeric(x$fit@stoFit@pathAvTheta[which(x$fit@stoFit@trajSubset==x$iter)[1], ])



  if(Tn > 0){
    D <- sim_data(
      SAMPLE_SIZE = x$n ,
      LOADINGS = settings[[settingLab]]$loads,
      THRESHOLDS = settings[[settingLab]]$thr,
      LATENT_COV = settings[[settingLab]]$lat_corr,
      SEED = x$id)
    H <- estimate_H(
      C_VEC = rep(length(settings[[settingLab]]$thr)+1, settings[[settingLab]]$p),
      A = settings[[settingLab]]$cnstr,
      THETA = bartheta,
      FREQ = freq,
      N = nrow(D),
      CORRFLAG = 1
    )$est_H
    Hinv <- solve(H)
    J <- estimate_J(
      Y = D,
      C_VEC = rep(length(settings[[settingLab]]$thr)+1, settings[[settingLab]]$p),
      A = settings[[settingLab]]$cnstr,
      THETA = bartheta,
      CORRFLAG = 1
    )$est_J

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
varII_2000$est_var <- est_varII
save(varII_2000, file = 'output/varII_2000.rds')

varII_5000 <- resII %>%
  mutate(
    path = map(fit, ~getThetaPath(.x, OPTION = 'raw'))
  ) %>%
  unnest(path) %>%
  filter(id %in% 1:2, n ==5000)
est_varII <- pbapply::pblapply(purrr::transpose(.l = varII_5000), function(x){

  n <- x$fit@dims@n
  Tn <- x$iter - x$fit@stoFit@control$BURN
  p <- x$fit@dims@p
  pairs <- p*(p-1)/2
  ppi <- x$fit@stoFit@control$PAIRS_PER_ITERATION
  freq <- x$fit@freq
  bartheta <- as.numeric(x$fit@stoFit@pathAvTheta[which(x$fit@stoFit@trajSubset==x$iter)[1], ])



  if(Tn > 0){
    D <- sim_data(
      SAMPLE_SIZE = x$n ,
      LOADINGS = settings[[settingLab]]$loads,
      THRESHOLDS = settings[[settingLab]]$thr,
      LATENT_COV = settings[[settingLab]]$lat_corr,
      SEED = x$id)
    H <- estimate_H(
      C_VEC = rep(length(settings[[settingLab]]$thr)+1, settings[[settingLab]]$p),
      A = settings[[settingLab]]$cnstr,
      THETA = bartheta,
      FREQ = freq,
      N = nrow(D),
      CORRFLAG = 1
    )$est_H
    Hinv <- solve(H)
    J <- estimate_J(
      Y = D,
      C_VEC = rep(length(settings[[settingLab]]$thr)+1, settings[[settingLab]]$p),
      A = settings[[settingLab]]$cnstr,
      THETA = bartheta,
      CORRFLAG = 1
    )$est_J

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
varII_5000$est_var <- est_varII
save(varII_5000, file = 'output/varII_5000.rds')



#### read and compare ####
load("output/settings.rda")

# I
settingLab <- 'I'
nthr <- length(settings[[settingLab]]$thr)*settings[[settingLab]]$p
nload <- sum(settings[[settingLab]]$cnstr)
ncorr <- settings[[settingLab]]$q*(settings[[settingLab]]$q-1)/2

emp_cov_bar_I <-  covI_1000 |>
  bind_rows(covI_2000) |>
  bind_rows(covI_5000) |>
  select(-id) |>
  group_by(n, p, q, lab, parId, iter) |>
  summarise_all(mean) |>
  mutate(par_type = 'thresholds',
         par_type = if_else(parId>nthr, 'loadings', par_type),
         par_type = if_else(parId>nthr+nload, 'correlations', par_type))

# II
settingLab <- 'II'
nthr <- length(settings[[settingLab]]$thr)*settings[[settingLab]]$p
nload <- sum(settings[[settingLab]]$cnstr)
ncorr <- settings[[settingLab]]$q*(settings[[settingLab]]$q-1)/2

emp_cov_bar_II <-  covII_1000 |>
  bind_rows(covII_2000) |>
  bind_rows(covII_5000) |>
  select(-id) |>
  group_by(n, p, q, lab, parId, iter) |>
  summarise_all(mean) |>
  mutate(par_type = 'thresholds',
         par_type = if_else(parId>nthr, 'loadings', par_type),
         par_type = if_else(parId>nthr+nload, 'correlations', par_type))




emp_cov_bar <- emp_cov_bar_I |>
  bind_rows(emp_cov_bar_II)


gg <- emp_cov_bar |>
  filter(!is.na(asy_cov)) |>
  filter(par_type!='thresholds') |>
  select(n, p, q, lab, parId, par_type,  iter, asymptotic = asy_cov, corrected = fin_cov) |>
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
  scale_fill_grey(start = .7, end = 1) +
  theme(legend.position = 'top')
gg

ggsave(gg, file = 'output/new_coverage.pdf', width = 10, height = 6)





emp_cov_tot <- emp_cov |>
  mutate(ev = 'theta_star') |>
  bind_rows(emp_cov_bar |>
              mutate(ev = 'theta_bar'))
gg <- emp_cov_tot |>
  filter(!is.na(asy_cov)) |>
  filter(n %in% c(1000, 2000, 5000)) |>
  mutate(par_type = 'thresholds',
         par_type = if_else(parId>nthr, 'loadings', par_type),
         par_type = if_else(parId>nthr+nload, 'correlations', par_type)) |>
  filter(par_type!='thresholds') |>
  select(n, p, q, lab, parId, par_type, ev, iter, asymptotic = asy_cov, corrected = fin_cov) |>
  gather(key = 'var', val = 'cov', asymptotic, corrected) |>
  mutate(Iter = factor(iter/n, levels = seq(0,5,by=.25), labels = paste0(seq(0,5,by=.25),'n'), ordered = T),
         n = factor(n, levels = settings[[settingLab]]$n, labels = paste0('n = ', settings[[settingLab]]$n), ordered = T),
         q = factor(q, levels = c(4, 8), labels = c('q = 4', 'q = 8'), ordered = T)
  ) |>
  ggplot(aes(y = cov, x = Iter, group = interaction(iter, var, ev))) +
  ggh4x::facet_nested(par_type+ev~q+n, scales = 'free')+
  scale_y_continuous(breaks = c(.5, .75, .9, .94, .95,.96, 1))+
  labs(y = 'Empirical coverage', x = 'Iterations', col = '') +
  theme_bw() +
  geom_hline(yintercept = .95, linetype = 'dashed') +
  geom_boxplot(aes(col = var)) +
  scale_color_grey(start = .7, end = .3) +
  scale_fill_grey(start = .7, end = 1) +
  theme(legend.position = 'top')
gg

saveRDS(emp_cov_tot, file = 'output/emp_cov_tot.rds')
ggsave(gg, file = 'output/coverage_comparison.pdf', width = 10, height = 8)











