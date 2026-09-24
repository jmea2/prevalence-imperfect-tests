
pacman::p_load(tidyverse, rio, binom, marginaleffects, rstan)

vih <- readRDS("vih.rds")


x <- vih[,c("edad_gr","sex","reactivo_sif1","reactivo_hepb",
            "tip_pob_HSH", "tip_pob_LGTBIQ", "tip_pob_OTRAS POBLACIONES","tip_pob_TRABAJADOR SEXUAL")]
N <- 11452  # Número de observaciones
P <- 8    # Número de covariables
X <- as.matrix.data.frame(x)
# Asignar sensibilidad y especificidad
Se <- 0.9999  # Sensibilidad
Sp <- 0.983066  # Especificidad
y_obs <- vih$vih_2 # variable dependiente con error de medida
# Preparar datos para Stan
data_list <- list(N = N, P = P, X = X, y_obs = y_obs, Se = Se, Sp = Sp)

#
stanc(file = "Nuevos_resultados/Modelo_logistico_Bayesiano/VIH2.stan")
# Especificar la ruta al archivo .stan
stan_file <- "Nuevos_resultados/Modelo_logistico_Bayesiano/VIH2.stan"
#COMPILAR MODELO 
stan_model <- stan_model(file = stan_file, auto_write = TRUE)
# Configuración para mejorar el rendimiento en RStan
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
# Ajustar el modelo en Stan
mod_vih2_bayesiano <- sampling(stan_model, data = data_list, chains = 4,
                               iter = 2000)


saveRDS(mod_vih2_bayesiano,"Nuevos_resultados/Modelos/mod_vih2_bayesiano.rds")
# Resumen de los resultados
# Visualizar los resultados

mod_vhi2_bayesiano <- stan_glm(vih_2 ~ edad_gr+sex+reactivo_sif1+
                reactivo_hepb+tip_pob_HSH+tip_pob_LGTBIQ+`tip_pob_OTRAS POBLACIONES`+
                `tip_pob_TRABAJADOR SEXUAL`, data = vih,
              family = binomial(link = "logit"),cores = 2, seed = 12345,prior = normal(0,12))

mod_vih2_bayesiano$coefficients
saveRDS(mod_vih2_bayesiano,"Nuevos_resultados/Modelos/mod_vih2_bayesiano.rds")

# calculo de prevalencia global con modelo regresion logistica 

prev_global <- predictions(mod_vih2_bayesiano, newdata = new_datavih) %>% 
  select(2:4)

#### intervalo de prevalencia aparente global y ajustada por Se y Sp REGRESION LOGISTICA
#### ESTANDAR

prevalencias <- prevalencias %>% 
  add_case(tibble(Disease = "VIH 2-Bayesiano",
                  prev_global,
                  Se = Se, Sp = Sp) %>% 
             mutate(P_ajustada = estimate-(1-Sp)/(1-(1-Se)-(1-Sp)),
                    lower_adj = conf.low-(1-Sp)/(1-(1-Se)-(1-Sp)), # corregida por Se y Sp del test
                    upper_adj = conf.high-(1-Sp)/(1-(1-Se)-(1-Sp)))) # orignal combinado 







