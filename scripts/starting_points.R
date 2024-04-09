#### Setup ####
library(dplyr)
library(tidyr)
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

#### starting points ####
lambda0_init <- c()
s <- 0

start1 <- c(rep(c(-1, 0, 1), settings[[settingLab]]$p), rep(0.5, nload),rep(0, ncorr))
start2 <- c(rep(c(-2, 0, 2), settings[[settingLab]]$p), rep(0.5, nload),rep(0, ncorr))
start3 <- c(rep(c(-.5, 0, .5), settings[[settingLab]]$p), rep(0.5, nload),rep(0, ncorr))
start4 <- c(rep(c(-2, 0, .5), settings[[settingLab]]$p), rep(0.5, nload),rep(0, ncorr))

start_tib <- tribble(
  ~start_lab, ~start_val,
  'start_1', start1,
  'start_2', start2,
  'start_3', start3,
  'start_4', start4,
)

####
setting_grid <- expand_grid(id = 1:nsim,
                            n = settings[[settingLab]]$n,
                            p = settings[[settingLab]]$p,
                            q = settings[[settingLab]]$q,
                            start_tib)


lab <- 'st8_starting_points'

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
      INIT = x$start_val,
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


#### mse ####
library(ggplot2)
library(purrr)
library(plFA)
res <- readRDS(paste0(path,'/st8_starting_points.rds'))

numII <- readRDS(paste0(path,'/num.rds'))
mse <- res %>%
  mutate(
    path = map(fit, ~getThetaPath(.x, OPTION = 'transformed'))
  ) %>%
  filter(n %in% c(1000, 2000, 5000)) |>
  unnest(path) %>%
  mutate(
    thresholds = map_dbl(par, ~mean((.x[1:nthr]-true_theta[1:nthr])^2)),
    loadings = map_dbl(par, ~mean((.x[(nthr+1):(nthr+nload)]-true_theta[(nthr+1):(nthr+nload)])^2)),
    correlations = map_dbl(par, ~mean((.x[(nthr+nload+1):(nthr+nload+ncorr)]-true_theta[(nthr+nload+1):(nthr+nload+ncorr)])^2))
  ) %>%
  gather(key = par_type, value = mse, thresholds, loadings, correlations) |>
  select(-fit,-par) |>
  group_by(n, p, q, lab, iter, start_lab, par_type) |>
  summarise(mse = mean(mse)) |>
  mutate(
    Iter = factor(iter/n, levels = seq(0,4,by=.25), labels = paste0(seq(0,4,by=.25),'n'), ordered = T),
    q = factor(q, levels = c(4, 8), labels = c('q = 4', 'q = 8'), ordered = T),
    n2 = factor(n, levels = settings[[settingLab]]$n, labels = paste0('n = ', settings[[settingLab]]$n), ordered = T))

nummse <- numII |>
  mutate(
    par = map(fit, ~getPar(.x, OPTION = 'transformed'))
  ) %>%
  mutate(
    thresholds = map_dbl(par, ~mean((.x[1:nthr]-true_theta[1:nthr])^2)),
    loadings = map_dbl(par, ~mean((.x[(nthr+1):(nthr+nload)]-true_theta[(nthr+1):(nthr+nload)])^2)),
    correlations = map_dbl(par, ~mean((.x[(nthr+nload+1):(nthr+nload+ncorr)]-true_theta[(nthr+nload+1):(nthr+nload+ncorr)])^2))
  ) %>%
  gather(key = par_type, value = mse, thresholds, loadings, correlations) |>
  select(-fit,-par) |>
  group_by(n, p, q, lab, par_type) |>
  summarise(mse = mean(mse)) |>
  mutate(
    q = factor(q, levels = c(4, 8), labels = c('q = 4', 'q = 8'), ordered = T),
    n2 = factor(n, levels = settings[[settingLab]]$n, labels = paste0('n = ', settings[[settingLab]]$n), ordered = T))

gg <- mse |>
  ggplot(aes(x = Iter, y = mse, group = interaction(par_type, start_lab))) +
  ggh4x::facet_nested(par_type~q+n2, scales = 'free') +
  scale_y_log10()+
  scale_x_discrete(breaks = paste0(seq(0,4,by=1),'n'), labels = c('0', 'n', '2n', '3n', '4n'))+
  geom_rect( aes(xmin=-Inf,xmax=3,ymin=0,ymax=Inf),
             fill="lightgrey", color =NA)+
  geom_hline(data = nummse, aes(yintercept=mse), linetype = 'dashed', col = 'darkgrey', alpha = .8)+
  theme_bw() +
  geom_line(alpha = .8, aes(col = start_lab)) +
  labs(y = 'Mean square error', x = 'Iterations', col = '') +
  scale_color_grey(start = .5, end = .2) +
  theme(legend.position = 'top')
gg
ggsave(gg, file = 'output/starting_points.pdf', width = 10, height = 6)

plotly::ggplotly(gg, dynamicTicks = T)
