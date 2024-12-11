install.packages("car")
library(betareg)

path <- "/home/vice/Documents/git/BIO310B/out/precision-beta-reg"
data <- read.csv(file.path(path, "Chembl28CCandD-loocv_ecfp4_for_beta_reg.csv"))
data <- read.csv(file.path(path, "Chembl28CCandD-loocv_fcfp4_for_beta_reg.csv"))
data <- read.csv(file.path(path, "Chembl28CCandD-loocv_maccs_for_beta_reg.csv"))
data <- read.csv(file.path(path, "e-loocv_ecfp4_for_beta_reg.csv"))
data <- read.csv(file.path(path, "e-loocv_fcfp4_for_beta_reg.csv"))
data <- read.csv(file.path(path, "e-loocv_maccs_for_beta_reg.csv"))
data <- read.csv(file.path(path, "global-loocv_ecfp4_for_beta_reg.csv"))
data <- read.csv(file.path(path, "global-loocv_fcfp4_for_beta_reg.csv"))
data <- read.csv(file.path(path, "global-loocv_maccs_for_beta_reg.csv"))
data <- read.csv(file.path(path, "gpcr-loocv_ecfp4_for_beta_reg.csv"))
data <- read.csv(file.path(path, "gpcr-loocv_fcfp4_for_beta_reg.csv"))
data <- read.csv(file.path(path, "gpcr-loocv_maccs_for_beta_reg.csv"))
data <- read.csv(file.path(path, "ic-loocv_ecfp4_for_beta_reg.csv"))
data <- read.csv(file.path(path, "ic-loocv_fcfp4_for_beta_reg.csv"))
data <- read.csv(file.path(path, "ic-loocv_maccs_for_beta_reg.csv"))
data <- read.csv(file.path(path, "nr-loocv_ecfp4_for_beta_reg.csv"))
data <- read.csv(file.path(path, "nr-loocv_fcfp4_for_beta_reg.csv"))
data <- read.csv(file.path(path, "nr-loocv_maccs_for_beta_reg.csv"))


sm <- lm(Precision ~ Score, data = data)
logit <- lm(Precision ~ Score + I(Score^2), data = data)

glmodel <- glm(Precision ~ Score, family = binomial(link = "logit"), data = data)



fit.train = lm(Precision ~ poly(Score,2), data=data)
summary(fit.train)
plot(fit.train)

library(fitdistrplus)
fit <- fitdist(data$Precision, "beta")
summary(fit)
# Kolmogorov-Smirnov
ks_test <- ks.test(data$Precision, "pbeta", fit$estimate[1], fit$estimate[2], 0, 1)
ks_test
# Anderson-Darling
library(goftest)
ad_test <- ad.test(data$Precision, "beta", shape1 = fit$estimate[1], shape2 = fit$estimate[2])
ad_test

modelo_beta <- betareg(Precision ~ Score, data = data, 
                       control = betareg.control(trace = TRUE)
                       )
summary(modelo_beta)
plot(modelo_beta)

plot(modelo_beta, which = 1,
     main = "Residuals vs Fitted")
plot(modelo_beta, which = 2,
     main = "Normal Q-Q Plot")
plot(modelo_beta, which = 3,
     main = "Scale-Location Plot")
plot(modelo_beta, which = 4,
     main = "Cook's Distance")


# ---
data("GasolineYield", package = "betareg")
gy_logit <- betareg(yield ~ temp, data = GasolineYield)
summary(gy_logit)
plot(gy_logit)
library(fitdistrplus)
fit <- fitdist(GasolineYield$yield, "beta")
summary(fit)
# Kolmogorov-Smirnov
ks_test <- ks.test(GasolineYield$yield, "pbeta", fit$estimate[1], fit$estimate[2], 0, 1)
ks_test
# Anderson-Darling
library(goftest)
ad_test <- ad.test(data$Precision, "beta", shape1 = fit$estimate[1], shape2 = fit$estimate[2])
ad_test
# ---
