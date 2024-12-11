#install.packages("minpack.lm")
library(ggplot2)
library(minpack.lm)

path <- "/home/vice/Documents/git/BIO310B/out/precision-loo"

data <- read.csv(file.path(path, "Chembl28CCandD-loocv_ecfp4_precision.csv"))
data <- read.csv(file.path(path, "Chembl28CCandD-loocv_fcfp4_precision.csv"))
data <- read.csv(file.path(path, "Chembl28CCandD-loocv_maccs_precision.csv"))
data <- read.csv(file.path(path, "e-loocv_ecfp4_precision.csv"))
data <- read.csv(file.path(path, "e-loocv_fcfp4_precision.csv"))
data <- read.csv(file.path(path, "e-loocv_maccs_precision.csv"))
data <- read.csv(file.path(path, "global-loocv_ecfp4_precision.csv"))
data <- read.csv(file.path(path, "global-loocv_fcfp4_precision.csv"))
data <- read.csv(file.path(path, "global-loocv_maccs_precision.csv"))
data <- read.csv(file.path(path, "gpcr-loocv_ecfp4_precision.csv"))
data <- read.csv(file.path(path, "gpcr-loocv_fcfp4_precision.csv"))
data <- read.csv(file.path(path, "gpcr-loocv_maccs_precision.csv"))
data <- read.csv(file.path(path, "ic-loocv_ecfp4_precision.csv"))
data <- read.csv(file.path(path, "ic-loocv_fcfp4_precision.csv"))
data <- read.csv(file.path(path, "ic-loocv_maccs_precision.csv"))
data <- read.csv(file.path(path, "nr-loocv_ecfp4_precision.csv"))
data <- read.csv(file.path(path, "nr-loocv_fcfp4_precision.csv"))
data <- read.csv(file.path(path, "nr-loocv_maccs_precision.csv"))

logistic_model <- nlsLM(
  Precision ~ L / (1 + exp(-k * (Score - x0))),
  data = data,
  start = list(L = max(data$Precision), k = 1, x0 = median(data$Score))
)

data$Predicted <- predict(logistic_model, newdata = data)

ggplot(data, aes(x = Score, y = Precision)) +
  geom_point() + geom_line(aes(y = Predicted)) +  
  theme_minimal() + labs(title = "Fit Logístico")

# ---
logm <- function(x, a, b, c) {
  c / (1 + exp(-a * (x - b)))
}

fit <- nlsLM(Precision ~ logm(Score, a, b, c), 
                data = data, 
                start = list(a = 1, b = mean(data$Score), c = 1))

cf <- coef(fit)
cat("Ecuación ajustada: Precision = ", cf["c"], "/ (1 + exp(-", cf["a"], " * (Score - ", cf["b"], ")))\n")

x_pred <- seq(min(data$Score), max(data$Score), length.out = 500)
y_pred <- logm(x_pred, cf["a"], cf["b"], cf["c"])

plot(data$Score, data$Precision, 
     xlab = "Score", ylab = "Precision", 
     pch = 1, cex = 0.5,
     main = "fit")
lines(x_pred, y_pred)

# ---
file_names <- c(
  "nr-loocv",
  "Chembl28CCandD-loocv",
  "e-loocv",
  "global-loocv",
  "gpcr-loocv",
  "ic-loocv"
)

generar_grafico <- function(path, nombre_archivo, output_path) {
  archivo_completo <- file.path(path, nombre_archivo)
  data <- read.csv(archivo_completo)
  
  fit <- glm(Precision ~ Score, 
             data = data, 
             family = binomial(link = "logit"))
  
  cf <- coef(fit)
  cat("logit(y) = ", cf["(Intercept)"], " + ", cf["Score"], " * x\n")
  
  x_pred <- seq(min(data$Score), max(data$Score), length.out = 500)
  y_pred <- 1 / (1 + exp(-(cf["(Intercept)"] + cf["Score"] * x_pred)))
  
  ecuacion <- paste0("logit(y) = ", 
                     round(cf["(Intercept)"], 3), 
                     " + ", 
                     round(cf["Score"], 3), 
                     " * x")
  
  archivo_salida <- file.path(output_path, paste0(tools::file_path_sans_ext(nombre_archivo), "_logx.png"))
  
  png(filename = archivo_salida, width = 1000, height = 800)
  par(mar = c(5, 4, 4, 2))
  
  plot(data$Score, data$Precision, 
       xlab = "Score", ylab = "Precision", 
       pch = 1, cex = 0.5,
       main = paste("fit log_X", nombre_archivo))
  lines(x_pred, y_pred)
  text(x = min(data$Score) + (max(data$Score) - min(data$Score)) * 0.05,
       y = max(data$Precision) * 0.8,
       labels = ecuacion, 
       cex = 0.8, 
       pos = 4)
  
  dev.off()
}

path <- "/home/vice/Documents/git/BIO310B/out/precision-loo-log-x"
output_path <- "/home/vice/Documents/git/BIO310B/img/R/logx"
if (!dir.exists(output_path)) {
  dir.create(output_path)
}

for (archivo in file_names) {
  generar_grafico(path, archivo, output_path)
}


# ---
# Función para generar el gráfico combinado por base de datos
generar_grafico_combinado <- function(path, base_datos, output_path) {
  # Identificar los archivos correspondientes a la base de datos
  archivos <- list.files(path, pattern = paste0("^", base_datos, "_.*\\.csv$"), full.names = TRUE)
  
  # Colores para los modelos
  colores <- c("red", "blue", "green")
  nombres_modelos <- c()
  curvas <- list()
  
  # Procesar cada archivo
  for (i in seq_along(archivos)) {
    data <- read.csv(archivos[i])
    
    # Ajustar el modelo logístico
    fit <- glm(Precision ~ Score, 
               data = data, 
               family = binomial(link = "logit"))
    
    cf <- coef(fit)
    cat("Archivo:", archivos[i], "\nlogit(y) = ", cf["(Intercept)"], " + ", cf["Score"], " * x\n")
    
    # Calcular las predicciones
    x_pred <- seq(min(data$Score), max(data$Score), length.out = 500)
    y_pred <- 1 / (1 + exp(-(cf["(Intercept)"] + cf["Score"] * x_pred)))
    
    # Guardar las curvas y el nombre del modelo
    curvas[[i]] <- list(x = x_pred, y = y_pred)
    nombres_modelos[i] <- tools::file_path_sans_ext(basename(archivos[i]))
  }
  
  # Generar el gráfico
  archivo_salida <- file.path(output_path, paste0(base_datos, "_combinado.png"))
  png(filename = archivo_salida, width = 600, height = 00)
  
  par(mar = c(5, 4, 4, 2))
  plot(NA, NA, 
       xlim = range(unlist(lapply(curvas, function(c) c$x))), 
       ylim = range(unlist(lapply(curvas, function(c) c$y))),
       xlab = "Score", ylab = "Precision",
       main = paste("Ajustes combinados -", base_datos))
  
  # Dibujar las curvas en el mismo gráfico
  for (i in seq_along(curvas)) {
    lines(curvas[[i]]$x, curvas[[i]]$y, col = colores[i], lty = 1)
  }
  
  # Añadir leyenda
  legend("topright", legend = nombres_modelos, col = colores, lty = 1, cex = 0.8)
  
  dev.off()
}

# Directorios
path <- "/home/vice/Documents/git/BIO310B/out/precision-loo"
output_path <- "/home/vice/Documents/git/BIO310B"
if (!dir.exists(output_path)) {
  dir.create(output_path)
}

# Nombres base de datos
file_names <- c("nr-loocv", "Chembl28CCandD-loocv", "e-loocv", 
                "global-loocv", "gpcr-loocv", "ic-loocv")

# Generar gráficos combinados para cada base de datos
for (base_datos in file_names) {
  generar_grafico_combinado(path, base_datos, output_path)
}
