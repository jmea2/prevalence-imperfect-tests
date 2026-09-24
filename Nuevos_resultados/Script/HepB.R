
install.packages("pacman")

pacman::p_load(logistic4p, tidyverse, rio, binom, gtsummary,
               janitor, fastDummies, marginaleffects)


# carga de datos

vih <- readRDS("vih.rds")

# carga de datos de combinaciones para marginalizacion estandarizada

new_datahepb <- vih %>%
  select(edad_gr,sex,vih_1,reactivo_sif1,tip_pob_HSH,
         tip_pob_LGTBIQ, `tip_pob_OTRAS POBLACIONES`,
         `tip_pob_TRABAJADOR SEXUAL`)

new_datahepb <- datagrid(newdata = new_datahepb)
saveRDS(new_datahepb,"Nuevos_resultados/Data/new_datahepb.rds")
new_dataHepB <- readRDS("Nuevos_resultados/Data/new_datahepb.rds")

######################## Se Sp del test utilizado #################

# PRUEBA HEPATITIS B   DETERMINE HBsAg 2 
# sensibilidad: 0.967
# especificidad:  0.991


# sensibilidades y especificidad combinadas para test Syphilis 3.0
se <- 0.967
sp <- 0.991
# calculo de r0 (Falsos positivos, FP)y r1 (falsos negativos, FN)

r0 <- 1-sp  # FP
r1 <- 1-se  # FN

####### ESTIMACION DE PREVALENCIA APARENTE CON  con ro y r1 desconocido ###########
##### metodo de Liu et cols. ESTIAMCION DE PREVALENCIA AJUSTADA(Covariables) Y AJUSTADA CORREGIDA por Se y Sp
x <- vih[,c("edad_gr","sex","vih_1","reactivo_sif1",
            "tip_pob_HSH", "tip_pob_LGTBIQ", "tip_pob_OTRAS POBLACIONES",
            "tip_pob_TRABAJADOR SEXUAL")] # predictors
y <- vih[,"reactivo_hepb"]      # outcome

mod_hepb_eq <- logistic4p(x, y, initial = ini, model='fp.fn', max.iter = 1000, detail = TRUE)
# # modelo fn.fp Warning: The algorithm does not converge after 1000 iterations.
# The current parameter estimates are 
# -10.07084 0.03077804 0.7266397 1.824813 2.09813 0.4629579 2.081515
# 1.736641 -3.004855 

# datos para modelo reducido, menos covariables y con valores iniciales programados
x <- vih[,c("edad_gr","sex","vih_1","reactivo_sif1")] # predictors
y <- vih[,"reactivo_hepb"]      # outcome

mod_hepb_eq <- logistic4p(x, y, initial = ini, model='equal', max.iter = 1000, detail = TRUE)
# Warning: The algorithm does not converge after 1000 iterations.
# The current parameter estimates are  0.001405356 0.1735142 -0.3729698 2.16583 0.157321 1.084756 
# Please increase the maximum of iterations, change the starting value or fit a different model.
print(mod_hepb_eq)

fp <- 
fn <- 

saveRDS(mod_hepb_eq, "Nuevos_resultados/Modelos/mod_hepb_eq.rds")

### obtener coeficientes del modelo de Liu

coef <- mod_hepb_eq$estimates[,1][-1]

# guarda coeficientes estiamdos

coef_models_hepb <- data.frame(coef) %>% 
  rename(coef_hepb_eq = coef)

saveRDS(coef_models_hepb, "Nuevos_resultados/Resultados/coef_models_hepb.rds")

# creamos modelo dummy con coeficiente del modelo de Liu
modelo_dummy <- glm(reactivo_hepb ~ edad_gr+sex+reactivo_sif1+
                      vih_1+tip_pob_HSH+tip_pob_LGTBIQ+`tip_pob_OTRAS POBLACIONES`+
                      `tip_pob_TRABAJADOR SEXUAL`, family = "binomial", data = vih)
# Reemplazamos los coeficientes
modelo_dummy$coefficients <- coef


# calculo de prevalencia global a la media e IC95
prev_global <- predictions(modelo_dummy, newdata = new_datahepb) %>% 
  select(2,4:5)
#### intervalo de prevalencia aparente global y ajustada por Se y Sp

prevalencias <- prevalencias %>% 
  add_case(
    tibble(Disease = "Hepatitis B - Liu.eq",
           prev_global,
           Se = 1-fn, Sp = 1-fp) %>% 
      mutate(P_ajustada = (estimate-fp)/(1-fn-fp),
             lower_adj = (conf.low-fp)/(1-fn-fp),
             upper_adj = (conf.high-fp)/(1-fn-fp))
  )
saveRDS(prevalencias,"Nuevos_resultados/Resultados/prevalencias.rds")

############# RESULTADOS POR  REGRESION LOGISITCA ESTANDAR ##############

mod_hepb_std <- glm(reactivo_hepb ~ edad_gr+sex+reactivo_sif1+
      vih_1+tip_pob_HSH+tip_pob_LGTBIQ+`tip_pob_OTRAS POBLACIONES`+
        `tip_pob_TRABAJADOR SEXUAL`, family = "binomial", data = vih)

summary(mod_hepb_std)

### obtener coeficientes del modelo de Liu

coef <- mod_hepb_std$coefficients

# guarda coeficientes estiamdos

coef_models_hepb <- data.frame(coef) %>% 
  rename(coef_hepb_std = coef)

saveRDS(coef_models_hepb, "Nuevos_resultados/Resultados/coef_models_hepb.rds")

# calculo de prevalencia global y varianza del estimador
prev_global <- predictions(mod_hepb_std, newdata = new_dataHepB) %>% 
  select(2,4:5)

#### intervalo de prevalencia aparente global y ajustada por Se y Sp REGRESION LOGISTICA
#### ESTANDAR

prevalencias <- prevalencias %>% add_case(
  tibble(Disease = "Hepatitis B-STD",
                  prev_global,
                  Se = se, Sp = sp) %>% 
             mutate(P_ajustada = estimate-(1-Sp)/(1-(1-Se)-(1-Sp)),
                    lower_adj = conf.low-(1-Sp)/(1-(1-Se)-(1-Sp)), # corregida por Se y Sp del test
                    upper_adj = conf.high-(1-Sp)/(1-(1-Se)-(1-Sp)))) # orignal combinado 


saveRDS(prevalencias,"Nuevos_resultados/Resultados/prevalencias.rds")
