
install.packages("pacman")

pacman::p_load(logistic4p, tidyverse, rio, binom,
               gtsummary, janitor, fastDummies,
               marginaleffects)


# carga de datos

vih <- readRDS("vih.rds")

# carga de datos de combinaciones para marginalizacion estandarizada

new_data <- readRDS("Nuevos_resultados/Data/new_datavih.rds")
######################## Se Sp de los test utilizados #################

# PRUEBA VIH 1   (ANTICUERPOS Y ANTIGENO) Marca determine Early detect 4a generacion
# sensibilidad: 0.995
# especificidad:  .994

# PRUEBA VIH 2   (ANTICUERPOS Y ANTIGENO) Bioline HIV 1/2 3.0  3a generacion

# sensibilidad: 0.980
# especificidad:  .989

# ####### Expresiones de Se y Sp de combinacion para el caso actual:###########
# sensitivity: (A)sen + (B)sen - [(A)sen x (B)sen]                            #
# specificity: (A)spec x (B)spec                                              #
###############################################################################

se_1 = 0.995
sp_1 = 0.994
se_2 = 0.980
sp_2 = 0.989

# Calculo de la Se y Sp de los test combinados, teniendo que si Test A es Positivo Ò
# el test B es positivo, otras cominaciones son negativas
# sensibilidades y especificidad combinadas para test vih
se <- se_1 + se_2 - (se_1*se_2) 
sp <- sp_1*sp_2 
# calculo de r0 (Falsos positivos, FP)y r1 (falsos negativos, FN)

r0 <- 1-sp  # FP
r1 <- 1-se  # FN

####### ESTIMACION DE PREVALENCIA APARENTE CON  con ro y r1 desconocido ###########
##### metodo de Liu et cols. ESTIAMCION DE PREVALENCIA AJUSTADA(Covariables) Y AJUSTADA CORREGIDA por Se y Sp

x <- vih[,c("edad_gr","sex","reactivo_sif1","reactivo_hepb",
            "tip_pob_HSH", "tip_pob_LGTBIQ", "tip_pob_OTRAS POBLACIONES","tip_pob_TRABAJADOR SEXUAL")] # predictors
y <- vih[,"vih_2"]      # outcome

mod_vih2_fpfn <- logistic4p(x, y, model='fp.fn', max.iter = 1000, detail = TRUE)
mod_vih2_eq <- logistic4p(x, y, model='equal', max.iter = 1000, detail = TRUE)

print(mod_vih2_eq)


### prueba de hipoteiss para modelo fp.fn Vs fp=fn

d <- -2*(-725.2589-(-724.923))

pvalue_prueba_chi2 <- stats::pchisq(d, 1,lower.tail = FALSE)

# no diferencia entre modelos es mejor el reducido Eq.

fp <- 0  # por modelo liu fp = fn = -0.007979 p-value = 0.258
fn <- 0  

saveRDS(mod_vih2_eq, "Nuevos_resultados/Modelos/mod_vih2_eq.rds")
saveRDS(mod_vih2_fpfn, "Nuevos_resultados/Modelos/mod_vih2_fpfn.rds")



# crear funcion de prediccion para el modelo de Liu
source("function_predict.R")

### obtener coeficientes del modelo de Liu

coef <- mod_vih2_eq$estimates[-1,1]

# modelo dummy para predicciones
modelo_dummy$coefficients <- coef

# guarda coeficientes estiamdos

coef_models_vih <- add_column(coef_models_vih) %>% 
  data.frame(coef) %>% 
  rename(mod_vih2eq = coef)

saveRDS(coef_models_vih, "Nuevos_resultados/Resultados/coef_models_vih.rds")

# calculo de prevalencia global y IC95%


prev_global <- predictions(modelo_dummy, newdata = new_data) %>% 
  select(2,4:5)

#### intervalo de prevalencia aparente global y ajustada por Se y Sp por Liu
prevalencias <- prevalencias %>% 
  add_case(
  tibble(Disease = "VIH 2 - Liu.eq",
                       prev_global,
                       Se = 1-fn, Sp = 1-fp) %>% 
  mutate(P_ajustada = (estimate-fp)/(1-fn-fp),
         lower_adj = (conf.low-fp)/(1-fn-fp),
         upper_adj = (conf.high-fp)/(1-fn-fp)))


saveRDS(prevalencias,"Nuevos_resultados/Resultados/prevalencias.rds")


#########################################################################
############# RESULTADOS POR  REGRESION LOGISITCA ESTANDAR ##############
##########################################################################

mod_vih2_std <- glm(vih_2 ~ edad_gr+sex+reactivo_sif1+
      reactivo_hepb+tip_pob_HSH+tip_pob_LGTBIQ+`tip_pob_OTRAS POBLACIONES`+
        `tip_pob_TRABAJADOR SEXUAL`, family = "binomial", data = vih)

summary(mod_vih2_std)

saveRDS(mod_vih2_std, "Nuevos_resultados/Modelos/mod_vih2_std.rds")

### obtener coeficientes del modelo de Liu

coef <- mod_vih2_std$coefficients

# guarda coeficientes estiamdos

coef_models_vih <- add_column(coef_models_vih) %>% 
data.frame(coef) %>% 
  rename(mod_vih2_std = coef)

saveRDS(coef_models_vih, "Nuevos_resultados/Resultados/coef_models_vih.rds")

# calculo de prevalencia global y varianza del estimador


prev_global <- predictions(mod_vih2_std, newdata = new_data) %>% 
  select(2,4:5)
  
#### intervalo de prevalencia aparente global y ajustada por Se y Sp REGRESION LOGISTICA
#### ESTANDAR

prevalencias <- prevalencias %>% 
  add_case(tibble(Disease = "VIH 2-STD",
                  prev_global,
                  Se = se, Sp = sp) %>% 
             mutate(P_ajustada = estimate-(1-Sp)/(1-(1-Se)-(1-Sp)),
                    lower_adj = conf.low-(1-Sp)/(1-(1-Se)-(1-Sp)), # corregida por Se y Sp del test
                    upper_adj = conf.high-(1-Sp)/(1-(1-Se)-(1-Sp)))) # orignal combinado 





saveRDS(prevalencias,"Nuevos_resultados/Resultados/prevalencias.rds")

