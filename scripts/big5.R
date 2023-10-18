#### Dependencies ####
library(haven)
library(plFA)
library(tidyverse)
library(RcppClock)
ifelse(!dir.exists('output/'), dir.create('output/'), "Folder already exists")
ifelse(!dir.exists('data/'), dir.create('data/'), "Folder already exists")

url <- "https://osf.io/download/t4rms/"
file_name <- "IPIP120.por.zip"
file_path <- "data/"
download.file(url, paste(file_path, file_name, sep = ""), mode = "wb")
unzip("data/IPIP120.por.zip", exdir = 'data/')

#######
B5 <-  read_spss("data/IPIP120.por")
B5 <- B5[, c(11:130)] - 1
B5m <- B5[rowSums(B5 < 0)==0,]
rowSums(B5m < 0)

#### Dictionary ####
# Neuroticism
factorj <- seq(from = 1, to = 91, by = 30)
itemOrder <- factorj
for ( lat in 1:5){
  factorj <- factorj+5
  itemOrder <- c(itemOrder, factorj)
}
facets <- c('Anxiety', 'Anger' , 'Depression', 'Self-Consciousness','Immoderation','Vulnerability')
facets <- rep(facets, each = 4)
dictN <- tibble(
  item = itemOrder,
  group = rep('Neuroticism', 24),
  facet = facets,
  lab = paste0('n', 1:24)
)

# Extraversion
factorj <- seq(from = 2, to = 92, by = 30)
itemOrder <- factorj
for ( lat in 1:5){
  factorj <- factorj+5
  itemOrder <- c(itemOrder, factorj)
}
facets <- c('Friendliness', 'Gregariousness' , 'Assertiveness', 'Activity level','Excitement seeking','Cheerfulness')
facets <- rep(facets, each = 4)
dictE <- tibble(
  item = itemOrder,
  group = rep('Extraversion', 24),
  facet = facets,
  lab = paste0('e', 1:24)
)

# Openness to experience
factorj <- seq(from = 3, to = 93, by = 30)
itemOrder <- factorj
for ( lat in 1:5){
  factorj <- factorj+5
  itemOrder <- c(itemOrder, factorj)
}
facets <- c('Imagination', 'Artistic interests' , 'Emotionality', 'Adventurousness','Intellect','Liberalism')
facets <- rep(facets, each = 4)
dictO <- tibble(
  item = itemOrder,
  group = rep('Openness to experience', 24),
  facet = facets,
  lab = paste0('o', 1:24)
)

# Agreeableness
factorj <- seq(from = 4, to = 94, by = 30)
itemOrder <- factorj
for ( lat in 1:5){
  factorj <- factorj+5
  itemOrder <- c(itemOrder, factorj)
}
facets <- c('Trust', 'Morality' , 'Altruism', 'Cooperation','Modesty','Sympathy')
facets <- rep(facets, each = 4)
dictA <- tibble(
  item = itemOrder,
  group = rep('Agreeableness', 24),
  facet = facets,
  lab = paste0('a', 1:24)
)

# Conscientiousness
factorj <- seq(from = 5, to = 95, by = 30)
itemOrder <- factorj
for ( lat in 1:5){
  factorj <- factorj+5
  itemOrder <- c(itemOrder, factorj)
}
facets <- c('Self-efficacy', 'Orderliness' , 'Dutilfulness', 'Achievement-striving','Self-discipline','Cautiousness')
facets <- rep(facets, each = 4)
dictC <- tibble(
  item = itemOrder,
  group = rep('Conscientiousness', 24),
  facet = facets,
  lab = paste0('c', 1:24)
)

# Join
dict <- bind_rows(dictN, dictE, dictO, dictA, dictC) %>%
  mutate_all(as.factor)

save(B5m, dict, file = 'data/Big5.rda')

#### Estimate model #####
rm(list=ls())
load(file = 'data/Big5.rda')
p <- ncol(B5m); q = 30
A <- build_constrMat(P = p, Q = q, STRUCT = 'simple')
itemOrder <- dict$item
train.id <- sample(1:nrow(B5m), .6*nrow(B5m))
B5train <- as.matrix(B5m[train.id, itemOrder])
B5val <- as.matrix(B5m[!(1:nrow(B5m)%in%train.id), itemOrder])
ppi <- 1
maxt <- 50000
burn <- 5000
my_seq <-  seq(0, 1, by = .1)*maxt
fit <- fit_plFA(
  DATA = B5train,
  CONSTR_LIST = list('CONSTRMAT' = A, 'CORRFLAG' = 1),
  VALDATA = B5val,
  METHOD = 'hyper',
  CONTROL = list(PAIRS_PER_ITERATION = ppi,
                 MAXT = maxt, BURN = burn,
                 ETA = .0001, EACHCLOCK = 1000,
                 PAR2 = .0001, EACHCHECK = 5000,
                 CHECKCONV = 1,
                 TOL = 1e-3, TOLCOUNT = 1, SEED = 1),#cpp_ctrl,
  ITERATIONS_SUBSET = my_seq,
  NCORES = 1
)
fit
fit@stoFit@pathValNll


#### Trajectories ####
nthr <- 4*p
nload <- sum(A)
ncorr <- q*(q-1)/2
path <- getThetaPath(fit)

path <- getThetaPath(fit) |>
  mutate(par = map(par, ~tibble(parId = 1:length(.x), par = .x))) |>
  unnest(par) |>
  mutate(par_type = map_chr(parId, ~{
    out <- 'correlations'
    if(.x %in% 1:nthr) out <- 'thresholds'
    if(.x %in% (nthr+1):(nthr+nload)) out <- 'loadings'
    return(out)
  }))

gg_traj <- path |>
  ggplot(aes(x = iter, y = par)) +
  geom_rect( aes(xmin=-Inf,xmax=burn,ymin=-Inf,ymax=Inf),
             fill="lightgrey", color =NA)+
  theme_bw() +
  facet_wrap(~par_type, scales = 'free') +
  geom_line(aes(group = parId), alpha = .5) +
  labs(x = 'Iterations', y = '')
gg_traj
gg_vnll <- as_tibble(fit@stoFit@pathValNll) |>
  ggplot(aes(x = V1, y = V2)) +
  theme_bw() +
  labs(x = '', y = 'Negative complete pairwise likelihood on validation data') +
  geom_line()
gg_checks <- ggpubr::ggarrange(gg_traj, gg_vnll, nrow = 1, widths = c(3, 1))
ggsave(gg_checks, file = 'output/gg_big5_checks.pdf', height = 5, width = 10)

#### Estimated parameters ####
facet_names <- c(paste0('N', 1:6), paste0('E', 1:6), paste0('O', 1:6), paste0('A', 1:6), paste0('C', 1:6))
items_names <- c(paste0('n', 1:24), paste0('e', 1:24), paste0('o', 1:24), paste0('a', 1:24), paste0('c', 1:24))

loads <- getPar(fit, OPTION = 'list')$loadings
loads.mat <- reshape::melt(loads) |>
  mutate(X1 = factor(X1, levels = 120:1, labels = rev(items_names), ordered = T),
         X2 = factor(X2, levels = 1:30, labels = facet_names, ordered = T))

latcorr <- getPar(fit, OPTION = 'list')$latent_correlations
latcorr.mat <- reshape::melt(latcorr) |>
  mutate(X1 = factor(X1, levels = 1:30, labels = facet_names, ordered = T),
         X2 = factor(X2, levels = 30:1, labels = rev(facet_names), ordered = T))


gg_estload <- loads.mat |>
  filter(value!=0)%>%
  ggplot(aes(x = X2, y = X1, fill = value))+
  geom_tile() +
  labs(x = '', y = '', fill = '')+
  theme_minimal() +
  theme( legend.position = 'top', plot.title=element_text(hjust=0.5), strip.text.x = element_blank(), panel.grid = element_blank()) +
  facet_wrap(vars(X2), scales = 'free') +
  scale_fill_gradient(low = 'white', high = 'black', limits = c(0,1.05))
# ggsave(gg_estload, file = 'output/gg_estload.pdf', height = 8, width = 10)

linecol <-'darkgrey'#  '#424242'
linety <- 'dashed'
linesi <- .5


gg_estcorr <- latcorr.mat |>
  mutate(sign = if_else(value >= 0, 'positive', 'negative')) |>
  ggplot(aes(x = X1, y = X2)) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = 'top', strip.text.x = element_blank(), strip.text.y = element_blank()) +
  labs(x = '', y = '', size = '', col = '') +
  geom_point(aes(size = abs(value), col = sign)) +
  scale_color_grey(start = .7, end = .3)  +
  geom_vline(xintercept = c(6.5, 12.5, 18.5, 24.5), col = linecol, linetype = linety, linewidth = linesi) +
  geom_hline(yintercept = c(6.5, 12.5, 18.5, 24.5), col = linecol, linetype = linety, linewidth = linesi)

#ggsave(gg_estcorr, file = 'output/gg_estcorr.pdf', height = 8, width = 10)

gg <- ggpubr::ggarrange(gg_estload, gg_estcorr, nrow = 1, widths = c(1,2))
gg
ggsave(gg, file = 'output/gg_big5_estpar.pdf', height = 7, width = 10)

