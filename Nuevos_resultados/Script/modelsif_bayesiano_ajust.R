pacman::p_load(rstan, rio, tidyverse)

# Cargar los datos

vih <- readRDS("vih.rds")
# Definir variables
N <- nrow(vih)
t1 <- vih$reactivo_sif1
x1 <- vih$edad_gr
x2 <- vih$sex
x3 <- vih$vih_1
x4 <- vih$reactivo_hepb
x5 <- vih$tip_pob_HSH
x6 <- vih$tip_pob_LGTBIQ
x7 <- vih$`tip_pob_OTRAS POBLACIONES`
x8 <- vih$`tip_pob_TRABAJADOR SEXUAL`

SN <- 0.964
SP <- 0.974

# Estimación inicial de betas con regresión logística
# zz <- glm(t1 ~ x1 + x2 + x3 + x4, data = data1, family = "binomial")
# betas_init <- as.numeric(zz$coefficients)

# Generar valores iniciales aleatorios para el estado latente d1
set.seed(123)
#d1_init <- rbinom(N, 1, 0.5)  # Supongamos que la mitad de la población está enferma

# Estructurar los datos para Stan
stan_data <- list(
  N = N,
  x1 = x1,
  x2 = x2,
  x3 = x3,
  x4 = x4,
  x5 = x5,
  x6 = x6,
  x7 = x7,
  x8 = x8,
  t1 = t1,
  SN = SN,
  SP = SP
)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

stanc(file = "Nuevos_resultados/Modelo_logistico_Bayesiano/Model_bayesiano1.stan")
# Especificar la ruta al archivo .stan
stan_file <- "Nuevos_resultados/Modelo_logistico_Bayesiano/Model_bayesiano1.stan"
#COMPILAR MODELO 
stan_model <- stan_model(file = stan_file, auto_write = TRUE)

sif_bayesiano_ajust <- sampling(stan_model, data = stan_data, chains = 4,
                                 iter = 4000)

saveRDS(sif_bayesiano_ajust,"Nuevos_resultados/Modelos/sif_bayesiano_ajust.rds")

# Ver resultados
print(sif_bayesiano_ajust, pars = c("betas"))

traceplot_sif_ajust <- traceplot(sif_bayesiano_ajust)
ggsave(plot = traceplot_sif_ajust, filename = "traceplot_sif_ajust.png",device = "png",
       path = "Nuevos_resultados/Modelos",dpi = 600, units = "cm",width = 17 ,height = 12)


o <- summary(sif_bayesiano_ajust)

# coeficientes
coef <- o$summary[,1][-10]

# errores estandar
se <- o$summary[,3][-10]

coef_models_sif <- coef_models_sif[,-8]
# adicion a tabla de coeficientes
coef_models_sif <- add_column(coef_models_sif) %>% 
  data.frame(coef) %>% 
  rename(sif_bayesiano_ajust = coef)

coef_models_sif <- add_column(coef_models_sif) %>% 
  data.frame(se) %>% 
  rename(se_sif_bayesiano_ajust = se)


saveRDS(coef_models_sif, "Nuevos_resultados/Resultados/coef_models_sif.rds")

nom_coeficientes <- names(mod_sif_std$coefficients)

names(coef) <- nom_coeficientes

coef

##### modelo dummy para predicciones
modelo_dummy <- glm(reactivo_sif1 ~ edad_gr+sex+vih_1+
                      reactivo_hepb+tip_pob_HSH+tip_pob_LGTBIQ+`tip_pob_OTRAS POBLACIONES`+
                      `tip_pob_TRABAJADOR SEXUAL`, family = "binomial", data = vih)

# Reemplazamos los coeficientes con los estimados con el Modelo de Liu

modelo_dummy$coefficients <- coef

prev_global <- predictions(modelo_dummy, newdata = new_datasif) %>% 
  select(2,5:6)

prevalencias[11,1] <- "Sifilis -Bayesiano_ajust"

prevalencias[11,"Se"] <- SN

prevalencias[11,"Sp"] <- SP

prevalencias[11,"P_ajustada"] <- prev_global[1]

prevalencias[11,"lower_adj"] <- prev_global[2]

prevalencias[11,"upper_adj"] <- prev_global[3]

saveRDS(prevalencias,"Nuevos_resultados/Resultados/prevalencias.rds")
