# Echantillonnage aleatoire simple ####

# _Moyenne de la pop ####

sigma_ybar <- function(sigma, n, N) {
  sigma / sqrt(n) * sqrt(1 - n / N)
}

s2 <- function(y, ybar, n) {
  sum( (y - ybar)^2 ) / (n - 1)
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