install.packages("pacman")

pacman::p_load(tidyverse, rio, binom,
               gtsummary, janitor, fastDummies,
               marginaleffects, rstanarm, bayesplot)


# carga de datos

# carga de datos de combinaciones para marginalizacion estandarizada

new_datasif <- readRDS("Nuevos_resultados/Data/new_datasif.rds")

# saveRDS(new_datasif,"Nuevos_resultados/Data/new_datasif.rds")

######################## Se Sp del test utilizado #################

# PRUEBA SIFILIS   Marca Bioline Syphilis 3.0 seleccion de limite inferior de IC95%
# sensibilidad: 0.964
# especificidad: 0.974


# sensibilidades y especificidad combinadas para test Syphilis 3.0
se <- 0.964
sp <- 0.974


# MODELO DE REGRESION LOGISRICA CON STAN
mod_sif_bayesiano <- stan_glm(reactivo_sif1 ~ edad_gr+sex+vih_1+
                   reactivo_hepb+tip_pob_HSH+tip_pob_LGTBIQ+`tip_pob_OTRAS POBLACIONES`+
                   `tip_pob_TRABAJADOR SEXUAL`, data = vih,
                 family = binomial(link = "logit"),cores = 2, seed = 12345,prior = normal(0,0.6),
                 prior_intercept =  normal(0,0.6))

saveRDS(mod_sif_bayesiano,"Nuevos_resultados/Modelos/mod_sif_bayesiano.rds")

posterior <- as.matrix(mod_sif_bayesiano)

plot_title <- ggtitle("Posterior distributions",
                      "with medians and 95% intervals")

# distribuciones posteriores
mcmc_areas(posterior,
           pars = c("reactivo_hepb","tip_pob_HSH","tip_pob_LGTBIQ"),
           prob = 0.95) + plot_title

mcmc_areas(posterior,
           pars = c("edad_gr","sex","vih_1"),
           prob = 0.95) + plot_title

# trace plot
plot_title <- ggtitle("Trace plot coefficients for syphilis Model")
traceplot_sifilis <- mcmc_trace(posterior, n_warmup = 300)+plot_title

ggsave(filename = "traceplot_sifilis.png", plot = traceplot_sifilis,device = "png",
       units = "cm",path = "Nuevos_resultados/Modelos",width = 25,height = 10,dpi = 300)

# inetrvals
mcmc_intervals(
  posterior,
  pars = character(),
  regex_pars = character(),
  transformations = list(),
  prob = 0.5,
  prob_outer = 0.9,
  point_est = c("median"),
  outer_size = 0.5,
  inner_size = 2,
  point_size = 4,
  rhat = numeric()
)
# re-arreglo coeficientes
coef <- mod_sif_bayesiano$coefficients
# errores standar
se <- mod_sif_bayesiano$ses

coef_models_sif <- add_column(coef_models_sif) %>% 
  data.frame(coef) %>% 
  rename(mod_sif_bayes = coef)

coef_models_sif <- add_column(coef_models_sif) %>% 
  data.frame(se) %>% 
  rename(se_sif_bayes = se)

saveRDS(coef_models_sif, "Nuevos_resultados/Resultados/coef_models_sif.rds")

# Resumen de los resultados
prev_global <- predictions(mod_sif_bayesiano, newdata = new_datasif) %>% 
  select(2:4)

#### intervalo de prevalencia aparente global y ajustada por Se y Sp REGRESION LOGISTICA
#### ESTANDAR

prevalencias <- prevalencias %>% 
  add_case(tibble(Disease = "Sifilis - Bayesiano",
                  prev_global,
                  Se = se, Sp = sp) %>% 
             mutate(P_ajustada = (estimate-(1-Sp))/(1-(1-Se)-(1-Sp)),
                    lower_adj = (conf.low-(1-Sp))/(1-(1-Se)-(1-Sp)), # corregida por Se y Sp del test
                    upper_adj = (conf.high-(1-Sp))/(1-(1-Se)-(1-Sp)))) # orignal combinado 

saveRDS(prevalencias,"Nuevos_resultados/Resultados/prevalencias.rds")

