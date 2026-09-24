
install.packages("pacman")

pacman::p_load(logistic4p, tidyverse, rio, binom,
               gtsummary, janitor, fastDummies,
               marginaleffects)


# carga de datos

vih <- readRDS("vih.rds")

# carga de datos de combinaciones para marginalizacion estandarizada

# new_datasif <- vih %>% 
#   select(edad_gr,sex,vih_1,reactivo_hepb,tip_pob_HSH,
#          tip_pob_LGTBIQ, `tip_pob_OTRAS POBLACIONES`,
#          `tip_pob_TRABAJADOR SEXUAL`)
# 
# new_datasif <- datagrid(newdata = new_datasif)

new_datasif <- readRDs("Nuevos_resultados/Data/new_datasif.rds")

# saveRDS(new_datasif,"Nuevos_resultados/Data/new_datasif.rds")

######################## Se Sp del test utilizado #################

# PRUEBA SIFILIS   Marca Bioline Syphilis 3.0 seleccion de limite inferior de IC95%
# sensibilidad: 0.964
# especificidad: 0.974


# sensibilidades y especificidad combinadas para test Syphilis 3.0
se <- 0.964
sp <- 0.974
# calculo de r0 (Falsos positivos, FP)y r1 (falsos negativos, FN)

r0 <- 1-sp  # FP
r1 <- 1-se  # FN

####### ESTIMACION DE PREVALENCIA APARENTE CON  con ro y r1 desconocido ###########
##### metodo de Liu et cols. ESTIAMCION DE PREVALENCIA AJUSTADA(Covariables) Y AJUSTADA CORREGIDA por Se y Sp

x <- vih[,c("edad_gr","sex","vih_1","reactivo_hepb",
            "tip_pob_HSH", "tip_pob_LGTBIQ", "tip_pob_OTRAS POBLACIONES","tip_pob_TRABAJADOR SEXUAL")] # predictors
y <- vih[,"reactivo_sif1"]      # outcome

mod_sif_fpfn <- logistic4p(x, y, model='fp.fn', max.iter = 1000, detail = TRUE)
mod_sif_eq <- logistic4p(x, y, model='equal', max.iter = 1000, detail = TRUE)

print(mod_sif_fpfn)
print(mod_sif_eq)

# estiamdos de FP y FN por Liu en sifilis

fn <- 0  # fp = fn = -0.06394875 p-value = 0.00567
fp <- 0

saveRDS(mod_sif_eq, "Nuevos_resultados/Modelos/mod_sif_eq.rds")
saveRDS(mod_sif_fpfn, "Nuevos_resultados/Modelos/mod_sif_fpfn.rds")

### obtener coeficientes del modelo de Liu eq

coef <- mod_sif_eq$estimates[-1,1]


# creamos modelo dummy con coeficionete del modelo de Liu
modelo_dummy <- glm(reactivo_sif1 ~ edad_gr+sex+vih_1+
                      reactivo_hepb+tip_pob_HSH+tip_pob_LGTBIQ+`tip_pob_OTRAS POBLACIONES`+
                      `tip_pob_TRABAJADOR SEXUAL`, family = "binomial", data = vih)
# Reemplazamos los coeficientes
modelo_dummy$coefficients <- coef

# guarda coeficientes estiamdos

coef_models_sif <- data.frame(coef) %>% 
  rename(coef_sif_eq = coef)

saveRDS(coef_models_sif, "Nuevos_resultados/Resultados/coef_models_sif.rds")

# calculo de prevalencia global y IC95%

prev_global <- predictions(modelo_dummy, newdata = new_datasif) %>% 
  select(2,4:5)

#### intervalo de prevalencia aparente global y ajustada por Se y Sp

prevalencias <- prevalencias %>% 
  add_case(
  tibble(Disease = "Sifilis - Liu.eq",
                       prev_global,
                       Se = 1-fn, Sp = 1-fp) %>% 
  mutate(P_ajustada = (estimate-fp)/(1-fn-fp),
         lower_adj = (conf.low-fp)/(1-fn-fp),
         upper_adj = (conf.high-fp)/(1-fn-fp))
  )

saveRDS(prevalencias,"Nuevos_resultados/Resultados/prevalencias.rds")


############# RESULTADOS POR  REGRESION LOGISITCA ESTANDAR ##############

mod_sif_std <- glm(reactivo_sif1 ~ edad_gr+sex+vih_1+
      reactivo_hepb+tip_pob_HSH+tip_pob_LGTBIQ+`tip_pob_OTRAS POBLACIONES`+
        `tip_pob_TRABAJADOR SEXUAL`, family = "binomial", data = vih)


summary(mod_sif_std)
saveRDS(mod_sif_std,"Nuevos_resultados/Modelos/mod_sif_std.rds")
### obtener coeficientes del modelo de Liu

coef <- mod_sif_std$coefficients

# guarda coeficientes estiamdos

coef_models_sif <- add_column(coef_models_sif) %>% 
data.frame(coef) %>% 
  rename(coef_sif_std = coef)

saveRDS(coef_models_sif, "Nuevos_resultados/Resultados/coef_models_sif.rds")

# calculo de prevalencia global con modelo regresion logistica 

prev_global <- predictions(mod_sif_std, newdata = new_datasif) %>% 
  select(2,4:5)

#### intervalo de prevalencia aparente global y ajustada por Se y Sp REGRESION LOGISTICA
#### ESTANDAR

prevalencias <- prevalencias %>% 
  add_case(tibble(Disease = "Sifilis-STD",
                  prev_global,
                  Se = se, Sp = sp) %>% 
             mutate(P_ajustada = estimate-(1-Sp)/(1-(1-Se)-(1-Sp)),
                    lower_adj = conf.low-(1-Sp)/(1-(1-Se)-(1-Sp)), # corregida por Se y Sp del test
                    upper_adj = conf.high-(1-Sp)/(1-(1-Se)-(1-Sp)))) # orignal combinado 

saveRDS(prevalencias,"Nuevos_resultados/Resultados/prevalencias.rds")
