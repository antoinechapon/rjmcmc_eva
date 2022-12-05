library(tidyverse)
library(reshape2)
library(lubridate)
library(rlist)
library(xtable)
library(nloptr)
library(zoo)
library(tictoc)
library(latex2exp)

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


toplot <- data.frame(x = yday(data$t),
                     y = data$surge)
ggplot(toplot) +
  geom_point(aes(x = x, y = y), alpha = .5) +
  geom_line(aes(x = x, y = thresh_surge$fitted.values), col = "red") +
  theme_grey(base_size = 12) +
  labs(x = "jour de l'année", y = "surcote (m)") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


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
             location.fun = ~ scale(1:nrow(pot_surge)))
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
  library(mvtnorm)
  tic()
  t_vec <- scale(1:nrow(input))
  # Lists of vectors indication which covariate is currently active.
  locat_l <- list(
    trend = rep(FALSE, 2),
    season = FALSE
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
    mu = c(fevd_init$results$par["location"], rep(NA, 4)),
    phi = c(fevd_init$results$par["log.scale"], rep(NA, 4)),
    xi = fevd_init$results$par["shape"]
  )
  lik_old <- exp(-fevd_init$results$value)
  # n_par <- unlist(fevd_init$results$num.pars) # useless line ?
  # Function to optimize to get mean(s) of the jump proposals.
  optim_fevd <- function(x, par_l, n_hyp) {
    cov_mu <- cbind(
      rep(1, length(nrow(input))),
      trend[, par_l[[1]]$trend],
      season[, rep(par_l[[1]]$season, 2)]
    )
    cov_phi <- cbind(
      rep(1, length(nrow(input))),
      trend[, par_l[[2]]$trend],
      season[, rep(par_l[[2]]$season, 2)]
    )
    mu_vec <- cov_mu %*% x[which(n_hyper == 1)]
    phi_vec <- cov_phi %*% x[which(n_hyper == 2)]
    xi_vec <- x[which(n_hyper == 3)]
    levd(x = input$y,
         threshold = thresh,
         location = mu_vec,
         scale = exp(phi_vec),
         shape = xi_vec,
         type = "PP")
  }
  out <- list()
  # Main RJMCMC loop.
  for (i in 1:n_mc) {
    if (i %in% seq(100, n_mc, 100)) {
      print(paste(i / n_mc * 100, "% completed."))
    }
    # STAY step:
    # Update every parameter with Metropolis. (NEED BOUNDED MH ON [-.5,.5] FOR XI)
    for (j in 1:length(par_old)) {
      # Iterating EVD parameters.
      for (k in 1:length(par_old[[j]])) {
        # Iterating the jth parameter hyperparameters.
        if (!is.na(par_old[[j]][k])) {
          par_new <- par_old
          par_new[[j]][k] <- rnorm(1, par_old[[j]][k], 1e-1)
          # Assemble the covariate matrix.
          cov_mu <- cbind(
            rep(1, length(nrow(input))),
            trend[, par_l[[1]]$trend],
            season[, rep(par_l[[1]]$season, 2)] # Need to take both columns.
          )
          cov_phi <- cbind(
            rep(1, length(nrow(input))),
            trend[, par_l[[2]]$trend],
            season[, rep(par_l[[2]]$season, 2)]
          )
          # Compute de vector of the parameters.
          mu_vec <- cov_mu %*% na.omit(par_new[[1]])
          phi_vec <- cov_phi %*% na.omit(par_new[[2]])
          xi_vec <- rep(par_new[[3]], nrow(input)) # Just in case varying xi is implemented later.
          # Likelihood of the new model.
          lik_new <- levd(x = input$y,
                          threshold = thresh,
                          location = mu_vec,
                          scale = exp(phi_vec),
                          shape = xi_vec,
                          type = "PP",
                          log = FALSE,
                          negative = FALSE)
          # Acceptance probability for stay step.
          alpha <- lik_new / lik_old
          if (runif(1) < alpha) {
            par_old[[j]][k] <- par_new[[j]][k]
            lik_old <- lik_new
          }
        }
      }
    }
    # JUMP step:
    # Pick a parameter for the jump.
    idx_par <- floor(runif(1, 1, length(par_l) + 1))
    # Pick a covariate of this parameter.
    idx_cov <- floor(runif(1, 1, length(par_l[[idx_par]]) + 1))
    # "idx_cov" also corresponds to the number of hyperparameter to change for
    # each covariate when 1=trend and 2= season, so used as such later (but bad
    # for generalisation of function to other covariates).
    # Select move (switching last T to F or first F to T).
    len_cov <- length(par_l[[idx_par]][[idx_cov]])
    sum_cov <- sum(par_l[[idx_par]][[idx_cov]])
    if (sum_cov == 0) {
      sum_new <- 1
    } else if (sum_cov == len_cov) {
      sum_new <- len_cov - 1
    } else {
      sum_new <- sum_cov + sign(rnorm(1))
    }
    # Compute jump and reverse jump relative probabilities.
    if (sum_cov == 0 | sum_cov == len_cov) {
      p_j <- 1
    } else {
      p_j <- .5
    }
    if (sum_new == 0 | sum_new == len_cov) {
      p_rj <- 1
    } else {
      p_rj <- .5
    }
    vec_new <- c(rep(TRUE, sum_new), rep(FALSE, len_cov - sum_new))
    par_l_new <- par_l
    par_l_new[[idx_par]][[idx_cov]] <- vec_new
    diff_dir <- sum(par_l_new[[idx_par]][[idx_cov]]) - sum_cov
    par_new <- par_old
    par_bool <- list(
      c(par_l_new[[1]][[1]], rep(par_l_new[[1]][[2]], 2)),
      c(par_l_new[[2]][[1]], rep(par_l_new[[2]][[2]], 2)))
    if (diff_dir > 0) {
      # Adding hyperparameter(s). Need to optimize their values to get mean
      # of the proposal (as in Min & Czado 2011 for the vine copula).
      par_new[[idx_par]][c(FALSE, par_bool[[idx_par]])] <- Inf
      vec <- na.omit(unlist(par_new))
      x0 <- vec
      lb <- vec
      ub <- vec
      n_hyper <- c(
        rep(1, length(na.omit(par_new[[1]]))),
        rep(2, length(na.omit(par_new[[2]]))),
        3
      )
      x0[which(x0 == Inf)] <- 0
      lb[which(lb == Inf)] <- -10
      ub[which(ub == Inf)] <- 10
      opt <- nloptr(
        x0 = x0,
        eval_f = optim_fevd,
        lb = lb,
        ub = ub,
        opts = list(
          algorithm = "NLOPT_LN_SBPLX", # bug with "NLOPT_LN_COBYLA"
          ftol_rel = 1e-2,
          max_eval = -1,
          print_level = 0
        ),
        par_l = par_l_new,
        n_hyp = n_hyper
      )
      # Removing the first hyperparameter that is never concerned by jumps.
      prop_mean <- opt$solution[n_hyper == idx_par][-1]
      draw <- rmvnorm(1,
                      mean = prop_mean,
                      sigma = diag(1e-2,
                                   length(prop_mean),
                                   length(prop_mean)))
      modif_par <- which(par_new[[idx_par]] == Inf)
      par_new[[idx_par]][modif_par] <- draw
    } else {
      # Removing hyperparameter.
      par_new[[idx_par]][
        c(FALSE, !c(par_l_new[[idx_par]][[1]], rep(par_l_new[[idx_par]][[2]], 2)))
          ] <- NA
    }
    # Assemble the covariate matrix.
    cov_mu <- cbind(
      rep(1, length(nrow(input))),
      trend[, par_l_new[[1]]$trend],
      season[, rep(par_l_new[[1]]$season, 2)]
    )
    cov_phi <- cbind(
      rep(1, length(nrow(input))),
      trend[, par_l_new[[2]]$trend],
      season[, rep(par_l_new[[2]]$season, 2)]
    )
    
    # Compute de vector of the parameters.
    mu_vec <- cov_mu %*% na.omit(par_new[[1]])
    phi_vec <- cov_phi %*% na.omit(par_new[[2]])
    xi_vec <- rep(par_new[[3]], nrow(input)) # Just in case varying xi is implemented later.
    # Likelihood of the new model.
    
    lik_new <- levd(x = input$y,
                    threshold = thresh,
                    location = mu_vec,
                    scale = exp(phi_vec),
                    shape = xi_vec,
                    type = "PP",
                    log = FALSE,
                    negative = FALSE)
    # Density of the proposal for the old and new hyperparameters.
    hyp_old <- par_old[[idx_par]][modif_par]
    if (all(is.na(hyp_old))) {
      dens_old <- 1
    } else {
      dval <- na.omit(par_old[[idx_par]][modif_par])
      dmean <- prop_mean[!is.na(par_old[[idx_par]][modif_par])]
      dens_old <- dnorm(dval[1], mean = dmean[1], sd = 1e-2)
      if (length(dval) == 2) {
        dens_old <- dens_old * dnorm(dval[2], mean = dmean[2], sd = 1e-2)
      }
    }
    hyp_new <- par_new[[idx_par]][modif_par]
    if (all(is.na(hyp_new))) {
      dens_new <- 1
    } else {
      dval <- na.omit(par_new[[idx_par]][modif_par])
      dmean <- prop_mean[!is.na(par_new[[idx_par]][modif_par])]
      dens_new <- dnorm(dval[1], mean = dmean[1], sd = 1e-2)
      if (length(dval) == 2) {
        dens_new <- dens_new * dnorm(dval[2], mean = dmean[2], sd = 1e-2)
      }
    }
    # Acceptance probability for jump step.
    # Priors cancel so are not included. (BUT ADD THEM IF SHAPE IMPLEMENTED)
    alpha <- lik_new / lik_old * p_rj / p_j * dens_old / dens_new
    if (!is.nan(alpha)) {
      if (runif(1) < alpha) {
        par_old <- par_new
        lik_old <- lik_new
        par_l <- par_l_new
      }
    }
    # Output result for this RJMCMC iteration.
    out <- list.append(out, unlist(par_old))
    model_id <- paste(as.numeric(!is.na(unlist(par_old))), collapse = "")
    names(out)[length(out)] <- model_id
  }
  # Remove burn-in period (taking first 10% by default).
  out <- out[-(1:(n_mc * .1))]
  # Vectors of the models visited. PUTTING INDEX IN mod_vec$id BY HAND FOR SIMPLICITY
  mod_vec <- list(id = names(out))
  mod_vec <- list.append(
    mod_vec, mu1 = as.factor(substr(mod_vec$id, 1, 3)))
  mod_vec <- list.append(
    mod_vec, mu2 = as.factor(substr(mod_vec$id, 4, 5)))
  mod_vec <- list.append(
    mod_vec, sigma1 = as.factor(substr(mod_vec$id, 6, 8)))
  mod_vec <- list.append(
    mod_vec, sigma2 = as.factor(substr(mod_vec$id, 9, 10)))
  # Split they output by model.
  out <- split(bind_rows(out), names(out))
  toc()
  return(list(models = out, ids = mod_vec))
}

run <- bdmcmc_nhpp(pot_surge,
                   thresh = thresh_surge$fitted.values,
                   fevd_init = init_surge,
                   n_mc = 2e4)
# save(run, file = "run.RData")
load("run.RData")


# REMOVE BURN-IN ####

toplot <- bind_cols(x = 1:length(run$ids[[1]]), run$ids[-1])
names(toplot)[-1] <- c("mu trend", "mu season",
                       "sigma trend", "sigma season")
toplot <- melt(toplot, measure.vars = -1)
toplot$value <- factor(toplot$value)

# Model jump trace plot.
ggplot(toplot) +
  geom_point(aes(x = x, y = value, group = variable)) +
  facet_grid(rows = vars(variable), scales = "free_y") +
  theme_grey(base_size = 12) +
  labs(x = "itération", y = "modèle") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


visit <- data.frame(mod = factor(names(run$models)))
visit$n <- NA
for (i in 1:nrow(visit)) {
  visit$n[i] <- nrow(run$models[[i]])
}
ggplot(visit) +
  geom_col(aes(x = n, y = mod)) +
  labs(x = "itérations", y = "modèle") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# MAKE SURE 10011111111 IS THE MOST VISITED ####

tohist <- bind_cols(x = 1:nrow(run$models$`10011111111`), run$models$`10011111111`)
tohist <- tohist[, -(3:4)]
tohist <- melt(tohist, measure.vars = -1)

# Not real trace plot since other models are visited in between.
ggplot(tohist) +
  geom_line(aes(x = x, y = value, group = variable)) +
  facet_wrap(facets = vars(variable), scales = "free_y") +
  labs(x = "itérations", y = "hyperparamètre") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggplot(tohist) +
  geom_histogram(aes(x = value, group = variable), bins = 30) +
  facet_wrap(facets = vars(variable), scales = "free_x") +
  labs(x = "itérations", y = "hyperparamètre") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

