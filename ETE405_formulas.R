# SIMPLE ####

# _Moyenne de la pop ####

sigma_ybar <- function(sigma, n, N) {
  sigma / sqrt(n) * sqrt(1 - n / N)
}

s2 <- function(y, n) {
  sum( (y - mean(y))^2 ) / (n - 1)
}

s_ybar <- function(s, n, N) {
  s / sqrt(n) * sqrt(1 - n / N)
}


# _Somme des elements de la pop ####

sigma_Yhat <- function(sigma, n, N) {
  N * sigma / sqrt(n) * sqrt(1 - n / N)
}

s_Yhat <- function(s, n, N) {
  N * s / sqrt(n) * sqrt(1 - n / N)
}


# _Estimation d'un rapport ####

sigma_Rhat <- function(Rhat, Xbar, y, x, n, N) {
  sqrt(1 - n / N) / (sqrt(n) * Xbar) * sqrt( sum((y - Rhat * x)^2) / N - 1 )
}

s_Rhat <- function(Rhat, Xbar, y, x, n, N) {
  sqrt(1 - n / N) / (sqrt(n) * Xbar) * sqrt( sum((y - Rhat * x)^2) / N - 1 )
}


# _Determination de la taille de l'echantillon ####

# Moyenne

n_mean <- function(u_alpha2, sigma, L, N) {
  n_0 <- (u_alpha2 * sigma / L)^2
  n_0 / (1 + n_0 / N)
}

# Proportion

n_prop <- function(u_alpha2, sigma, L, N, p) {
  n_0 <- p * (1 - p) / L^2 * u_alpha2^2
  n_0 / (1 + n_0 / N)
}

# Correlation spatiale et temporelle

# rho_l : temporelle, rho_c : spatiale

var_ybar <- function(n, ns, sigma, rho_l, rho_c) {
  sigma^2 / (n * ns) * ( 1 + 2/n * sum(n - l) * rho_l ) * ( 1 + rho_c * (ns - 1) )
}

n_cor <- function(u_alpha2, sigma, L, rho_l, rho_c, ns) {
  l <- 1:length(rho_l)
  D <- (u_alpha2 * sigma / L)^2
  D / (2 * ns) * (1 + 2*sum(rho_l) + rho_c*(ns - 1) * (1 + 2*sum(rho_l))) +
    D / 2 * ( 1 / ns^2 * (1 + 2*sum(rho_l) + rho_c*(ns - 1) * (1 + 2*sum(rho_l)))^2 -
                8 / (ns*D) * (sum(l * rho_l) - rho_c * (ns - 1) * sum(l * rho_l)) )^(1 / 2)
}



# _Echantillonnage aleatoire simple avec variable supplementaire

# estimation par un rapport

var_Yhat_r <- function(x, y, n, N, Rhat) {
  N^2/n * (1 - n/N) * sum((y - Rhat*x)^2) / (N - 1)
}

var_Yhatbar_r <- function(x, y, n, N, Rhat) {
  1/n * (1 - n/N) * sum((y - Rhat*x)^2) / (N - 1)
}

# estimation par regression

ybar_r <- function(ybar, b, Xbar, xbar){
  ybar + b * (Xbar - xbar)
}
# coef b fixe a priori
var_y_r_fixe <- function(b, x, y, n, N) {
  s_x2 <- 1/(n-1) * sum((x - mean(x))^2)
  s_y2 <- 1/(n-1) * sum((y - mean(y))^2)
  s_xy2 <- 1/(n-1) * sum((x - mean(x))* (y - mean(y)))
  var <- (1 - n/N) / n * (s_y2 - 2*b*s_xy2 + b^2 * s_x2)
  return(var)
}
# coef b determine a partir de l'echantillon
coef_b <- function(x, y) {
  sum((x - mean(x))* (y - mean(y))) / sum((x - mean(x))^2)
}
var_y_r_dete <- function(x, y, n, N) {
  s_xy2 <- 1/(n-1) * sum((x - mean(x))* (y - mean(y)))
  rho <- s_xy2 / (sqrt(s2(x)) * sqrt(s2(y)))
  var <- (1 - n/N) / n * s2(y) * (1 - rho^2)
}



# STRATIFIE ####

# _Stratifie proportionnel ####

nh_prop <- function(Nh, n, N){
  Nh * n / N
}

ybar_strat <- function(N, Nh, yh) {
  1/N * sum(Nh * yh)
}

var_ybar_strat <- function(N, Nh, nh, yh_list) { # generale
  yh_vec <- rep(NA, length(yh_list))
  for (i in 1:length(yh_list)) {
    yh_vec[[i]] <- s2(yh_list[[i]], nh[i])
  }
  sum( (Nh/N)^2 * (1 - nh/Nh) * yh_vec / nh )
}

var_ybar_strat2 <- function(N, Nh, nh, sh) { # avec sh^2
  sum( (Nh/N)^2 * (1 - nh/Nh) * sh^2 / nh )
}


var_Yhat_strat <- function(Nh, nh, yh_list){ # generale
  sh2 <- rep(NA, length(nh))
  for (i in 1:lenght(nh)) {
    sh2[i] <- s2(yh_list[[i]], nh[i])
  }
  sum( Nh^2 * (1 - nh/Nh) * sh2 / nh )
}


# _Stratifie optimum ####

n_optim <- function(Nh, C, C0, ch, sigmah) {
  (C - C0) * sum(Nh * sigmah / sqrt(ch)) / sum(Nh * sigmah * sqrt(ch))
}

nh_optim <- function(n, Nh, ch, sigmah) {
  n * Nh * sigmah / (sqrt(ch) * sum(Nh * sigmah / sqrt(ch)))
}

var_yst_min <- function(N, n, sigmah) {
  1/N^2 * (1/n * sum(Nh * sigmah)^2 - sum(Nh * sigmah^2))
}


# _Stratitife pour estimation de proportions ####

p_st <- function(ph, Nh, N) {
  sum(Nh * ph / N)
}

var_P_st <- function(Nh, N, nh, Ph) {
  sum( (Nh/N)^2 * Ph*(1-Ph)/nh * (1 - nh/Nh) )
}

var_p_st <- function(Nh, N, nh, Ph) {
  sum( (Nh/N)^2 * Ph*(1-Ph)/nh )
}

# estim proportion avec strate optimales (p. 46)
nh_optim_estprop <- function(n, Nh, Ph, ch) {
  n * (Nh * sqrt(Ph*(1-Ph)/ch)) / sum(Nh * sqrt(Ph*(1-Ph)/ch))
}



# SYSTEMATIQUE ####







#######################

# exemple 13, page 42

Nh <- c(394, 461, 391, 334, 169, 113, 148)
sh <- c(8.3, 13.3, 15.1, 19.8, 24.5, 26, 35.2)

var_ybar_strat2(N = sum(Nh), Nh = Nh, nh = nh, sh = sh)

round(nh_optim(n = 100, Nh = Nh, ch = rep(1, length(Nh)), sigmah = sh))
  
var_yst_min(N = sum(Nh), n = 100, sigmah = sh)

