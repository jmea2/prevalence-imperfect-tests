######## ESTIMACION BAYESIANA DE COEFICIENTES DE REGRESION ##########

pacman::p_load(tidyverse, rio, binom, marginaleffects, rstan)

vih <- import("vih.rds")
new_datavih <- import("Nuevos_resultados/Data/new_datavih.rds")
x <- vih[,c("edad_gr","sex","reactivo_sif1","reactivo_hepb",
            "tip_pob_HSH", "tip_pob_LGTBIQ", "tip_pob_OTRAS POBLACIONES","tip_pob_TRABAJADOR SEXUAL")]
N <- 11452  # Número de observaciones
P <- 8    # Número de covariables
X <- as.matrix.data.frame(x)
# Asignar sensibilidad y especificidad
Se <- 0.9751  # Sensibilidad
Sp <- 0.999934  # Especificidad
y_obs <- vih$vih_1 # variable dependiente con error de medida

# Preparar datos para Stan
data_list <- list(N = N, P = P, X = X, y_obs = y_obs, Se = Se, Sp = Sp)

#
stanc(file = "Nuevos_resultados/Modelo_logistico_Bayesiano/VIH1.stan")
# Especificar la ruta al archivo .stan
stan_file <- "Nuevos_resultados/Modelo_logistico_Bayesiano/VIH1.stan"
#COMPILAR MODELO 
stan_model <- stan_model(file = stan_file)
# Configuración para mejorar el rendimiento en RStan
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
# Ajustar el modelo en Stan
mod_vih1_bayesiano <- sampling(stan_model, data = data_list, iter = 2000, chains = 4)

fit2 <- stan_glm(vih_1 ~ edad_gr+sex+reactivo_sif1+
           reactivo_hepb+tip_pob_HSH+tip_pob_LGTBIQ+`tip_pob_OTRAS POBLACIONES`+
           `tip_pob_TRABAJADOR SEXUAL`, data = vih,
         family = binomial(link = "logit"),cores = 2, seed = 12345,prior = normal(0,0.6),
         prior_intercept =  normal(0,0.6))

mod_vih1_bayesiano <- fit2
saveRDS(mod_vih1_bayesiano,"Nuevos_resultados/Modelos/mod_vih1_bayesiano.rds")
# Resumen de los resultados
# Visualizar los resultados

# trace plot de coeficientes B e intercepto
traceplot_vih1 <- traceplot(mod_vih1_bayesiano, pars = c("alpha","beta"))

ggsave(plot = traceplot_vih1, filename = "traceplot_vih1.png",device = "png",
       path = "Nuevos_resultados/Modelos",dpi = 600, units = "cm",width = 17 ,height = 12)

hist(posterior_samples,xlim = c(0,0.15) ,breaks = 100, main = "Distribución Posterior de intercepto", xlab = "alpha")

# extraido media del posterior para cada coeficiente
coef <- get_posterior_mean(mod_vih1_bayesiano)[1:9,5]
# re-arreglo coeficientes
coef <- mod_vih1_bayesiano$coefficients
# errores
se <- mod_vih1_bayesiano$ses

coef_models_vih <- add_column(coef_models_vih) %>% 
  data.frame(coef) %>% 
  rename(mod_vih1_bayes = coef)

coef_models_vih <- add_column(coef_models_vih) %>% 
  data.frame(se) %>% 
  rename(se_vih1_bayes = se)


saveRDS(coef_models_vih, "Nuevos_resultados/Resultados/coef_models_vih.rds")


prev_global <- predictions(mod_vih1_bayesiano, newdata = new_datavih) %>% 
  select(2:4)

#### intervalo de prevalencia aparente global y ajustada por Se y Sp REGRESION LOGISTICA
#### ESTANDAR

prevalencias <- prevalencias %>% 
  add_case(tibble(Disease = "VIH 1-Bayesiano",
                  prev_global,
                  Se = Se, Sp = Sp) %>% 
             mutate(P_ajustada = (estimate-(1-Sp))/(1-(1-Se)-(1-Sp)),
                    lower_adj = (conf.low-(1-Sp))/(1-(1-Se)-(1-Sp)), # corregida por Se y Sp del test
                    upper_adj = (conf.high-(1-Sp))/(1-(1-Se)-(1-Sp)))) # orignal combinado 

saveRDS(prevalencias,"Nuevos_resultados/Resultados/prevalencias.rds")





