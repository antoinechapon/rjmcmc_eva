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
                   threshold = thresh_surge$fitted.values,
                   type = "PP")
init_surge
# init_hs <- fevd(x = y,
#                 data = pot_hs,
#                 use.phi = TRUE,
#                 threshold = thresh_hs$fitted.values)



test <- fevd(x = y,
             data = pot_surge,
             use.phi = TRUE,
             type = "PP",
             method = "GMLE",
             threshold = thresh_surge$fitted.values,
             location.fun = ~ cbind(pot_surge$cos, pot_surge$sin),
             scale.fun = ~ cbind(pot_surge$cos, pot_surge$sin))
test

levd(x = pot_surge$y,
     threshold = thresh_surge$fitted.values,
     location = init_surge$results$par[1],
     scale = exp(init_surge$results$par[2]),
     shape = init_surge$results$par[3],
     type = "PP")

init_surge$results$par

bdmcmc_nhpp <- function(input,
                        thresh,
                        fevd_init,
                        n_mc = 1e4) {
  # BDMCMC with mu and sigma dependent on linear or quadratic time and season.
  # Using phi instead of sigma.
  
  t_vec <- scale(1:nrow(input))

  # Lists of vectors indication which covariate is currently active.
  locat_l <- list(
    trend = rep(FALSE, 2),
    season = F
  )
  scale_l <- list(
    trend = rep(FALSE, 2),
    season = FALSE
  )
  par_l <- list(locat_l, scale_l)
  # Matrices of covariates.
  trend <- cbind(t_vec, t_vec^2)
  season <- cbind(input$cos, input$sin)

  # Initializing the parameter vector and the loglikelihood.
  par_old <- list(
    mu = fevd_init$results$par["location"],
    phi = fevd_init$results$par["log.scale"],
    xi = fevd_init$results$par["shape"]
  )
  lik_old <- exp(-fevd_init$results$value)
  
  n_par <- unlist(fevd_init$results$num.pars)
  
  out <- list()
  
  for (i in 1:n_mc) {
    # Stay step: update parameters with Metropolis-Hastings (actually only
    # Metropolis for mu and phi).
    
    for (j in 1:length(par_old)) {
      # Iterating EVD parameters.
      for (k in length(par_old[[j]])) {
        # Iterating the jth parameter hyperparameters.
        par_new <- par_old
        
        par_new[[j]][k] <- rnorm(1, par_old[[j]][k], 1e-1)
        
        # Assemble the covariate matrix.
        cov_mu <- cbind(
          rep(1, length(nrow(input))),
          trend[, par_l[[1]]$trend],
          season[, par_l[[1]]$season]
        )
        cov_phi <- cbind(
          rep(1, length(nrow(input))),
          trend[, par_l[[2]]$trend],
          season[, par_l[[2]]$season]
        )
        # Compute de vector of the parameters.
        mu_vec <- cov_mu %*% par_new[[1]]
        phi_vec <- cov_phi %*% par_new[[2]]
        xi_vec <- rep(par_new[[3]], nrow(input)) # Just in case varying xi is implemented later.
        
        # Likelihood of new model.
        lik_new <- levd(x = input$y,
                        threshold = thresh,
                        location = mu_vec,
                        scale = exp(phi_vec),
                        shape = xi_vec,
                        type = "PP",
                        log = FALSE,
                        negative = FALSE)
        
        alpha <- lik_new / lik_old
        
        if (runif(1) < alpha) {
          par_old[[j]][k] <- par_new[[j]][k]
          lik_old <- lik_new
        }
        out <- list.append(out, unlist(par_old)) # ITS BAD TO GROW LIST ON MCMC ?
        model_id <- paste(as.numeric(unlist(par_l)), collapse = "")
        names(out)[length(out)] <- model_id
      }
    }
  }
  return(out)
}

zut <- bdmcmc_nhpp(pot_surge,
                   thresh = thresh_surge$fitted.values,
                   fevd_init = init_surge,
                   n_mc = 1e2)

zut2 <- bind_rows(zut)

hist(zut2$mu0, 30)
hist(zut2$mu1, 30)
hist(zut2$log.scale, 30)
