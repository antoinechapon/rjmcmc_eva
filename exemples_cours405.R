# exemple 1

fish <- c(rep(42, 23),
          rep(41, 4),
          36,
          32,
          29,
          27, 27,
          23,
          19,
          16, 16,
          15, 15,
          14,
          11,
          10,
          9,
          7,
          rep(6, 3),
          5, 5,
          4,
          3)

fish_mean <- mean(fish)

s2_fish <- s2(fish, fish_mean, 50)

fish_mean + qt(0.1, 49) * s_ybar(sqrt(s2_fish), 50, 676)
fish_mean + qt(0.9, 49) * s_ybar(sqrt(s2_fish), 50, 676)


# exemple 2

sum(fish)

Yhat <- 676 / 50 * sum(fish)

Yhat + qt(0.1, 49) * s_Yhat(sqrt(s2_fish), 50, 676)
Yhat + qt(0.9, 49) * s_Yhat(sqrt(s2_fish), 50, 676)


