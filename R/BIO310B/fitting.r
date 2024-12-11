install.packages("extraDistr")
install.packages("viridis")
install.packages("ggthemes")
install.packages("npreg")
library(ggplot2)


path <- "/home/vice/Documents/git/BIO310B/out/precision-loo"
data <- read.csv(file.path(path, "e-loocv_ecfp4_precision.csv"))

path <- "/home/vice/Documents/git/BIO310B/out/precision-loo-log-x"
data <- read.csv(file.path(path, "Chembl28CCandD-loocv_ecfp4_precision_log_x.csv"))
data <- read.csv(file.path(path, "nr-loocv_ecfp4_precision_log_x.csv"))

# --- Models
model <- lm(Precision ~ Score, data = data) # lineal
model <- lm(Precision ~ Score + I(Score^2), data = data) # logistico
model <- lm(Precision ~ poly(Score, 2), data = data) # polinomial grado 2
model <- lm(Precision ~ poly(Score, 3), data = data) # polinomial grado 3
model <- lm(Precision ~ poly(Score, 4), data = data) # polinomial grado 4
model <- lm(Precision ~ poly(Score, 15), data = data) # polinomial grado 5
model <- lm(Precision ~ poly(Score, 8), data = data) # polinomial grado 6

# --- Plots
summary(model)
ggplot(data, aes(x=Score, y=Precision)) + 
  geom_point(color='#2980B9', size = 1) +
  geom_smooth(method=lm, formula = y ~ poly(x, 15), color='#2C3E50')

# --- Test
library(extraDistr)
library(viridis)
library(ggthemes)
mean_precision <- mean(data$Precision)
size <- sd(data$Score) * 10

beta_mu_phi <- ggplot() +
  geom_function(fun = dprop, args = list(mean = mean_precision, size = size),
                aes(color = paste0("dprop(mean = ", round(mean_precision, 3), 
                                   ", size = ", round(size, 1), ")")),
                size = 1) +
  scale_color_viridis_d(option = "plasma", end = 0.8, name = "",
                        guide = guide_legend(ncol = 1)) +
  labs(title = "Mean- and precision-based beta distributions",
       x = "Precision", y = "Density") +
  theme_clean() +
  theme(legend.position = "bottom")

print(beta_mu_phi)
# ---
n <- 101
x <- seq(0, 1, length.out = n)
fx <- sin(2 * pi * x)
set.seed(1)
y <- fx + rnorm(n, sd = 0.5)
plot(x, y)             # data
lines(x, fx, lwd = 2)  # f(x)

library(npreg)
x <- data$Score
y <- data$Precision
mod.ss <- ss(x, y, nknots = 10)
print(mod.ss)
plot(x, y, main = "Splines", col = "blue", pch = 16,
     xlab = "Score", ylab = "Precision")
lines(x, predict(mod.ss, x), col = "red", lwd = 2)
