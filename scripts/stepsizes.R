#### Setup ####
library(dplyr)
library(tidyr)
library(plFA)
load("output/settings.rda")

#### read setting ####
nsim <- 1000
settingLab <- 'II'

path <- paste0('output/', settingLab)
ifelse(!dir.exists(path), dir.create(path), "Folder already exists")
nthr <- length(settings[[settingLab]]$thr)*settings[[settingLab]]$p
nload <- sum(settings[[settingLab]]$cnstr)
ncorr <- settings[[settingLab]]$q*(settings[[settingLab]]$q-1)/2
theta <- settings[[settingLab]]$theta
true_theta <- settings[[settingLab]]$true_theta


#### stepsizes ####
setting_grid <- expand_grid(id = 1:nsim,
                            n = settings[[settingLab]]$n,
                            p = settings[[settingLab]]$p,
                            q = settings[[settingLab]]$q,
                            eta = c(.02, .01, .005, .0025, .00125))


lab <- 'st8_stepsizes'

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
      CONTROL = list(PAIRS_PER_ITERATION = ppi, MAXT = maxt, BURN = burn,  ETA = x$eta),#cpp_ctrl,
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

#### mse ####
library(ggplot2)
library(purrr)
library(plFA)
res <- readRDS(paste0(path,'/st8_stepsizes.rds'))

mse <- res |>
  mutate(
    path = map(fit, ~try(getThetaPath(.x, OPTION = 'transformed'))),
    err = map_lgl(path, ~if_else(typeof(.x)=='list', F, T))
  ) |>
  filter(!err) |>
  unnest(path) |>
  filter(iter == 1*n) |>
  mutate(
    mse = map_dbl(par, ~mean((.x-true_theta)^2))
    ) |>
  select(-fit,-par) |>
  group_by(n, p, q, lab, eta) |>
  summarise(mse = mean(mse)) |>
  mutate(
    q = factor(q, levels = c(4, 8), labels = c('q = 4', 'q = 8'), ordered = T),
    n2 = factor(n, levels = settings[[settingLab]]$n, labels = paste0('n = ', settings[[settingLab]]$n), ordered = T))

library(latex2exp)
gg <- mse |>
  ggplot(aes(x = eta, y = mse)) +
  ggh4x::facet_nested(~q+n2, scales = 'free') +
  scale_y_log10()+
  theme_bw() +
  geom_point() +
  geom_line() +
  labs(y = 'Mean square error ( after 0.5n iterations)', x = TeX(r'($\eta_0$)')) +
  scale_color_grey(start = .5, end = .2)
gg
ggsave(gg, file = 'output/stepsizes.pdf', width = 10, height = 3)

plotly::ggplotly(gg, dynamicTicks = T)


msetype <- res |>
  mutate(
    path = map(fit, ~try(getThetaPath(.x, OPTION = 'transformed'))),
    err = map_lgl(path, ~if_else(typeof(.x)=='list', F, T))
  ) |>
  filter(!err) |>
  filter(n %in% c(1000, 2000, 5000)) |>
  unnest(path) |>
  mutate(
    thresholds = map_dbl(par, ~mean((.x[1:nthr]-true_theta[1:nthr])^2)),
    loadings = map_dbl(par, ~mean((.x[(nthr+1):(nthr+nload)]-true_theta[(nthr+1):(nthr+nload)])^2)),
    correlations = map_dbl(par, ~mean((.x[(nthr+nload+1):(nthr+nload+ncorr)]-true_theta[(nthr+nload+1):(nthr+nload+ncorr)])^2))
  ) |>
  gather(key = par_type, value = mse, thresholds, loadings, correlations) |>
  select(-fit,-par) |>
  group_by(n, p, q, lab, iter, par_type, eta) |>
  summarise(mse = mean(mse)) |>
  mutate(
    Iter = factor(iter/n, levels = seq(0,4,by=.25), labels = paste0(seq(0,4,by=.25),'n'), ordered = T),
    q = factor(q, levels = c(4, 8), labels = c('q = 4', 'q = 8'), ordered = T),
    n2 = factor(n, levels = settings[[settingLab]]$n, labels = paste0('n = ', settings[[settingLab]]$n), ordered = T))

msetype |>
  ggplot(aes(x = Iter, y = mse, group = interaction(par_type, eta))) +
  ggh4x::facet_nested(par_type~q+n2, scales = 'free') +
  scale_y_log10()+
  scale_x_discrete(breaks = paste0(seq(0,4,by=1),'n'), labels = c('0', 'n', '2n', '3n', '4n'))+
  geom_rect( aes(xmin=-Inf,xmax=3,ymin=0,ymax=Inf),
             fill="lightgrey", color =NA)+
  theme_bw() +
  geom_line(alpha = .8, aes(col = as.factor(eta))) +
  labs(y = 'Mean square error', x = 'Iterations', col = '') +
  # scale_color_grey(start = .5, end = .2) +
  theme(legend.position = 'top')
