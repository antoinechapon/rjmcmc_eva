library(tidyverse)
library(reshape2)
library(lubridate)
library(rlist)
library(xtable)
library(nloptr)
library(zoo)
library(tictoc)

library(extRemes)
library(quantreg)
library(VineCopula)

library(hdrcde)

load("minip_surge.Rdata")
load("t_vec.Rdata")

mean_imp <- apply(minip_surge, 1, mean)
mean_imp <- na.approx(mean_imp)

data <- bind_cols(t = t_vec, surge = mean_imp)


load("hs_50_1.Rdata")


t_common <- as.POSIXct(intersect(t_vec, hs_50_1$t), tz = "UTC", origin = "1970-01-01")


data <- data.frame(
  t = t_common,
  surge = mean_imp[t_vec %in% t_common],
  hs = hs_50_1$y[hs_50_1$t %in% t_common]
)

data$cos <- cos(yday(data$t) / (365 + leap_year(data$t)) * 2 * pi)
data$sin <- sin(yday(data$t) / (365 + leap_year(data$t)) * 2 * pi)

# cop <- BiCopSelect(pobs(data$surge),
#                    pobs(data$hs))
# contour(cop)



# CDE avec cos sin de year en covar
# - permet genere var1 en conction de year
# - get dyncop en fonction de var1 et year
# - get var2 en fonction de dyncop et year
#
# Manière semi-paramétrique de générer des obs en 2D sans casser le cycle annuel
# commu aux deux varialbe.
# A partir de ces valeurs generes B fois, GEV avec block maxima annuel.
# 
# Bien justfier qu'on voulait faire du resampling pour que ça corresponde au cours
# mais makheureursement on a du faire du parametric bootstrap, pour que ce soit
# adapté aux extremes.
# 
# Et ca permet de faire CI en 2D, vu que CI c'est dans le cours
# (en utilisant le truc de Serinaldi ? ou alors suffit de sampler en Volpi Fiori
# dans le resultat du bootstrap pour CI à la Serinaldi ?)





# str(data)
# 
# 
# 
# library(gamlss)
# 
# gam_surge <- gamlss(formula = surge ~ yr_cos + yr_sin,
#                     sigma.formula = surge ~ yr_cos + yr_sin,
#                     data = data)
# 
# gam_hs <- gamlss(formula = hs ~ yr_cos + yr_sin,
#                  sigma.formula = hs ~ yr_cos + yr_sin,
#                  data = data[, c(3, 4, 5)])
# 
# # In both cases the global deviance is lower with varying sigma.
# 
# predict(gam_hs,
#         newdata = data.frame(yr_cos = c(0.1, 0.1),
#                              yr_sin = c(0.1, 0.5)),
#         type = "response")
# 


######################################################

# Nonstationary bivariate bayesian bootstrap for EVA

library(extRemes)

# yr_max <- data.frame(
#   t = unique(year(data$t))
# )
# temp_surge <- aggregate(data$surge, list(year(data$t)), FUN = max)
# temp_hs <- aggregate(data$hs, list(year(data$t)), FUN = max)
# yr_max$surge <- temp_surge$x
# yr_max$hs <- temp_hs$x
# 
# 
# cop <- BiCopSelect(pobs(yr_max$surge),
#                    pobs(yr_max$hs))
# 
# contour(cop)
# BiCopPar2TailDep(cop)



#########################

# BDMCMC pour margins POT et dyncop
# avec Jumps pour les covars : season, linear, quadratic
# et Jumps de la famille de la pair-cop.

thresh_surge <- quantreg::rq(surge ~ cos * sin, data = data, tau = .99)
thresh_hs <- quantreg::rq(hs ~ cos * sin, data = data, tau = .99)

dec_surge <- decluster(data$surge, thresh_surge$fitted.values, run = 2)
dec_hs <- decluster(data$hs, thresh_hs$fitted.values, run = 2)

pot_surge <- data.frame(
  t = data$t,
  y = as.numeric(dec_surge),
  cos = data$cos,
  sin = data$sin
)
pot_hs <- data.frame(
  t = data$t,
  y = as.numeric(dec_hs),
  cos = data$cos,
  sin = data$sin
)


init_surge <- fevd(x = y,
                   data = pot_surge,
                   use.phi = TRUE,
                   threshold = thresh_surge$fitted.values)
init_hs <- fevd(x = y,
                data = pot_hs,
                use.phi = TRUE,
                threshold = thresh_hs$fitted.values)


bdmcmc_marg <- function(input,
                        n_mc = 1e4) {
  # BDMCMC with mu and sigma dependent on linear or quadratic time and season.
  # using phi instead of sigma in the algo.
  
  # Matrices of the models space. One for each covariate.
  # Dim 1 is for mu, dim 2 for phi (xi could be added with 3d arrays).
  mat_time <- matrix(FALSE, nrow = 3, ncol = 3)
  mat_seas <- matrix(FALSE, nrow = 2, ncol = 2)
  
  mats <- list(time = mat_time, seas = mat_seas)
  
  t_vec <- 1:nrow(input)
  
  for (i in 1:n_mc) {
    # Pick one covariate.
    
  }
}

bdmcmc_marg(pot_surge)