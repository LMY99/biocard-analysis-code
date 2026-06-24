set.seed(2)
source("functions_S.R")
usePackage("splines2")
usePackage("TruncatedNormal")
usePackage("mvtnorm")
usePackage("matrixStats")
usePackage("foreach")
library(progress)
library(doParallel)
library(doRNG)
library(rstan)
options(warn = 0)
cv_exclude <- as.integer(commandArgs(trailingOnly = T)[1])

# Loading Data ---------------------------------------------------
load("biocard.rda")
load('cv_split.rda')

# Create consecutive pseudo-IDs for each individual for easy coding
unique.IDs <- sort(unique(df$Study.ID))
df$ID <- match(df$Study.ID, unique.IDs)
df_test <- df[cv_index[df$ID]==cv_exclude,]
df <- df[cv_index[df$ID]!=cv_exclude,]
K <- 11 # Number of biomarkers
# All non-age covariates + Intercept
X <- as.matrix(cbind(intercept = 1, df[, c("apoe", "SEX", "education")]))
Y <- df[, 1:K] # Biomarkers array
t <- df$ageori # Age in original scale
unique.IDs <- sort(unique(df$Study.ID))
df$ID <- match(df$Study.ID, unique.IDs)

# All non-age covariates
Y <- as.matrix(Y, ncol = K) # Biomarkers array
y_obs <- matrix(nrow = nrow(Y), ncol = ncol(Y))
for (i in 1:nrow(Y)) {
  for (j in 1:ncol(Y)) {
    y_obs[i, j] <- as.numeric(!is.na(Y[i, j]))
  }
}
Y[is.na(Y)] <- 0
dat <- list(
  Nobs = nrow(X), Npred = ncol(X),
  Nout = ncol(Y), Nind = max(df$ID),
  x = X, t = df$ageori, y = Y, y_obs = y_obs,
  jj = df$ID
)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = FALSE)
R <- 1e4

stan.fit <- stan(file = "logistic.stan", data = dat, iter = R, chains = 1, seed = 419)
stan.array <- rstan::extract(stan.fit)
rm(stan.fit)
gc()

X_test <- as.matrix(cbind(intercept = 1, df_test[, c("apoe", "SEX", "education")]))
FE_fit <- X_test %*% t(apply(stan.array$beta,c(2,3),mean))
scale_fit <- colMeans(stan.array$lscale)
pos_fit <- colMeans(stan.array$lpos)
amp_fit <- colMeans(stan.array$lamp)
logit_fit <- matrix(ncol=ncol(FE_fit),nrow=nrow(FE_fit))
for(i in 1:ncol(FE_fit))
  logit_fit[,i] <- plogis((df_test$ageori - pos_fit[i])/scale_fit[i])*amp_fit[i]
Y_test_fit <- FE_fit + logit_fit
Y_test <- df_test[,1:11]
cv_error <- mean(colMeans((Y_test - Y_test_fit)^2, na.rm = TRUE))
var_noise <- mean(stan.array$sigmaerror)
var_between <- mean(stan.array$sigmarandom)
save(cv_error, var_noise, var_between, 
     Y_test, Y_test_fit,
     file = paste('para_cv',cv_exclude,'.rda',sep=''))