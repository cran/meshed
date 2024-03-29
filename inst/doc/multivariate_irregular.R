## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align="center",
  message=FALSE
)

## -----------------------------------------------------------------------------
library(magrittr)
library(dplyr)
library(ggplot2)
library(meshed)

set.seed(2021)

SS <- 25 # coord values for jth dimension 
dd <- 2 # spatial dimension
n <- SS^2 # number of locations
q <- 3 # number of outcomes
k <- 2 # true number of spatial factors used to make the outcomes
p <- 3 # number of covariates

coords <- cbind(runif(n), runif(n)) %>% 
  as.data.frame() 
colnames(coords) <- c("Var1", "Var2")

clist <- 1:q %>% lapply(function(i) coords %>% 
                          mutate(mv_id=i) %>% 
                          as.matrix()) 

philist <- c(5, 10)

# cholesky decomp of covariance matrix
LClist <- 1:k %>% lapply(function(i) t(chol(
  exp(-philist[i] * as.matrix(dist(clist[[i]])) ))))

# generating the factors
wlist <- 1:k %>% lapply(function(i) LClist[[i]] %*% rnorm(n))

# factor matrix
WW <- do.call(cbind, wlist)

# factor loadings
Lambda <- matrix(0, q, ncol(WW))
diag(Lambda) <- runif(k, 1, 2)
Lambda[lower.tri(Lambda)] <- runif(sum(lower.tri(Lambda)), -1, 1)

# nuggets
tau.sq <- rep(.05, q)
TTsq <- matrix(1, nrow=n) %x% matrix(tau.sq, ncol=length(tau.sq))
# measurement errors
EE <- ( rnorm(n*length(tau.sq)) %>% matrix(ncol=length(tau.sq)) ) * TTsq^.5

XX <- matrix(rnorm(n*p), ncol=p)
Beta <- matrix(rnorm(p*q), ncol=q)

# outcome matrix, fully observed
YY_full <- XX %*% Beta + WW %*% t(Lambda) + EE

# .. introduce some NA values in the outcomes
# all at different locations
YY <- YY_full
for(i in 1:q){
  YY[sample(1:n, n/5, replace=FALSE), i] <- NA
}

## ---- fig.show="hold", fig.width=2.5, fig.height=2.5--------------------------
simdata <- coords %>%
  cbind(data.frame(Outcome_full=YY_full, 
                   Outcome_obs = YY)) 

simdata %>%
  ggplot(aes(Var1, Var2, color=Outcome_obs.2)) + 
  geom_point() + scale_color_viridis_c() +
  theme_minimal() + theme(legend.position="none")

## -----------------------------------------------------------------------------
mcmc_keep <- 200 # too small! this is just a vignette.
mcmc_burn <- 400
mcmc_thin <- 2

mesh_total_time <- system.time({
  meshout <- spmeshed(y=YY, x=XX, coords=coords, k = 2,
                    grid_size = c(20, 20), 
                    block_size = 16, 
                    n_samples = mcmc_keep, 
                    n_burn = mcmc_burn, 
                    n_thin = mcmc_thin, 
                    n_threads = 2,
                    prior = list(phi=c(2, 20)),
                    verbose=0
  )})

## ---- echo=FALSE, include=FALSE-----------------------------------------------
plot_cube <- function(cube_mcmc, q, k, Ptrue, Pname, full=F){
  oldpar <- par(mar=c(2.5,2,1,1), mfrow=c(q,k))
  for(i in 1:q){
    for(j in 1:k){
      if(full|(j<=i)){
        plot(density(cube_mcmc[i, j,], bw=.05), main=paste0(Pname, i,j))
        abline(v=Ptrue[i,j], col="red")
      } else {
        plot(c(0))
      }
    }
  }
  par(oldpar)
}

## -----------------------------------------------------------------------------
plot_cube(meshout$lambda_mcmc, ncol(YY), 2, Lambda, "Lambda")

## -----------------------------------------------------------------------------
plot_cube(meshout$beta_mcmc, ncol(YY), p, Beta, "Beta", T)

## -----------------------------------------------------------------------------
mcmc_summary <- function(x) c(quantile(x, .025), mean(x), quantile(x, .975))
y_post_sample <- meshout$yhat_mcmc %>% 
  abind::abind(along=3) %>%
  apply(1:2, mcmc_summary)

# posterior mean for 3rd outcome:
meshout$coordsdata %>% cbind(y_pm_3 = y_post_sample[2,,3]) %>%
  filter(forced_grid==0) %>%
  ggplot(aes(Var1, Var2, color=y_pm_3)) +
  geom_point() +
  scale_color_viridis_c() +
  theme_minimal() + theme(legend.position="none")

## -----------------------------------------------------------------------------
perf <- meshout$coordsdata %>% 
  cbind(y_pm_3 = y_post_sample[2,,3]) %>%
  left_join(simdata, by = c("Var1", "Var2"))

perf %>% filter(forced_grid==0, !complete.cases(Outcome_obs.3)) %>% 
  with(cor(y_pm_3, Outcome_full.3))

