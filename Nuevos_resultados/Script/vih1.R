
install.packages("pacman")

pacman::p_load(logistic4p, tidyverse, rio, binom, gtsummary, janitor, fastDummies, marginaleffects)


vih <- readRDS("vih.rds")
# combinaciones para marginalizacion estandarizada

new_data <- vih %>% 
  select(edad_gr,sex,reactivo_sif1,reactivo_hepb,tip_pob_HSH,
        tip_pob_LGTBIQ, `tip_pob_OTRAS POBLACIONES`, `tip_pob_TRABAJADOR SEXUAL`) %>% 
 datagrid(newdata = new_data)

######################## Se Sp de los test utilizados #################

# PRUEBA VIH 1   (ANTICUERPOS Y ANTIGENO) Marca determine Early detect 4a generacion
# limite inferior del IC95%
# sensibilidad: 0.995
# especificidad:  0.994

# PRUEBA VIH 2   (ANTICUERPOS Y ANTIGENO) Bioline HIV 1/2 3.0  3a generacion

# sensibilidad: 0.98
# especificidad:  .989


####### estimacion de prevalencia de VIH aparente con ro y r1 desconocido #######

# Calculo de la Se y Sp de los test combinados, teniendo que si Test A es negativo
# el test B no requiere ser aplicado.
# el primero debe dar positivo para aplicar el segundo (ambos positivos definen un caso)

# ####### Expresiones: #############################
# sensitivity: (A)sen x (B)sen                     #
# specificity: (A)spec + [1 - (A)spec] x (B)spec   #
####################################################

 se_1 = 0.995
 sp_1 = 0.994
 se_2 = 0.980
 sp_2 = 0.989

# sensibilidades y especificidad combinadas para test vih
se <- se_1*se_2   
sp <- sp_1+(1 - sp_1)*sp_2 
# calculo de r0 y r1

r0 <- 1-sp # fp
r1 <- 1-se # fn

##### metodo de Liu et cols. ESTIAMCION DE PREVALENCIA AJUSTADA Y AJUSTADA CORREGIDA

x <- vih[,c("edad_gr","sex","reactivo_sif1","reactivo_hepb",
            "tip_pob_HSH", "tip_pob_LGTBIQ", "tip_pob_OTRAS POBLACIONES","tip_pob_TRABAJADOR SEXUAL")] # predictors
y <- vih[,"vih_1"]      # outcome

# modelo Liu con fp = fp  SI converge
mod_vih1_eq <- logistic4p(x, y,model='equal', max.iter = 1000, detail = TRUE)

# fp y fn con valores iniciales no converge
mod_vih1_fpfn <- logistic4p(x, y, initial = c(r0,r1,model_vih1_std$coefficients),model='fp.fn', max.iter = 1000, detail = TRUE)


print(mod_vih1_eq) # trabajo con este para las estimaciones de prevalencia

fp <- 0  # por modelo liu fp = fn 
fn <- 0  

saveRDS(mod_vih1_eq, "Nuevos_resultados/Modelos/mod_vih1_eq.rds")

# crear funcion de prediccion para el modelo de Liu
source("function_predict.R")

### obtener coeficientes del modelo de Liu

coef <- mod_vih1_eq$estimates[,1][-1]

# guarda coeficientes estiamdos

coef_models_vih <- data.frame(coef) %>% 
  rename(modvih1_liu_eq = coef)

saveRDS(coef_models_vih, "Nuevos_resultados/Resultados/coef_models_vih.rds")

######### calculo de prevalencia global y varianza del estimador ########

# para el caso de modelo de Liu al no ser un objeto tipo modelo libreria
# marginal effects no lo acepta entonces hacemos un modelo dummy

# Ajustamos un modelo logístico dummy

modelo_dummy <- glm(vih_1 ~ edad_gr+sex+reactivo_sif1+
                      reactivo_hepb+tip_pob_HSH+tip_pob_LGTBIQ+`tip_pob_OTRAS POBLACIONES`+
                      `tip_pob_TRABAJADOR SEXUAL`, family = "binomial", data = vih)

# Reemplazamos los coeficientes con los estimados con el Modelo de Liu

modelo_dummy$coefficients <- coef
prev_global <- predictions(modelo_dummy, newdata = new_data) %>% 
  select(2,4:5)

#### intervalo de prevalencia aajustado por ro y r1

prevalencias <- tibble(Disease = "VIH 1 - Liu.eq",
                    prev_global,
       Se = 1-fn, Sp = 1-fp) %>% 
  mutate(P_ajustada = (estimate-fp)/(1-fn-fp),
         lower_adj = (conf.low-fp)/(1-fn-fp),
         upper_adj = (conf.high-fp)/(1-fn-fp))

saveRDS(prevalencias,"Nuevos_resultados/Resultados/prevalencias.rds")

############# RESULTADOS POR  REGRESION LOGISITCA ESTANDAR ##############

mod_vih1_std <- glm(vih_1 ~ edad_gr+sex+reactivo_sif1+
      reactivo_hepb+tip_pob_HSH+tip_pob_LGTBIQ+`tip_pob_OTRAS POBLACIONES`+
        `tip_pob_TRABAJADOR SEXUAL`, family = "binomial", data = vih)


summary(mod_vih1_std)

saveRDS(mod_vih1_std,"Nuevos_resultados/Modelos/mod_vih1_std.rds")

### obtener coeficientes del Modelo estandar

coef <- mod_vih1_std$coefficients

# guarda coeficientes estiamdos

coef_models_vih <- add_column(coef_models_vih) %>% 
data.frame(coef) %>% 
  rename(modvih1_std = coef)

saveRDS(coef_models_vih, "Nuevos_resultados/Resultados/coef_models_vih.rds")

# calculo de prevalencia global con modelo regresion logistica 

prev_global <- predictions(mod_vih1_std, newdata = new_data) %>% 
  select(2,4:5)

#### intervalo de prevalencia aparente global y ajustada por Se y Sp REGRESION LOGISTICA
#### ESTANDAR

prevalencias <- prevalencias %>% 
  add_case(tibble(Disease = "VIH 1-STD",
                    prev_global,
                    Se = se, Sp = sp) %>% 
  mutate(P_ajustada = estimate-(1-Sp)/(1-(1-Se)-(1-Sp)),
         lower_adj = conf.low-(1-Sp)/(1-(1-Se)-(1-Sp)), # corregida por Se y Sp del test
         upper_adj = conf.high-(1-Sp)/(1-(1-Se)-(1-Sp)))) # orignal combinado 

saveRDS(prevalencias,"Nuevos_resultados/Resultados/prevalencias.rds")

################# ESTIMACION PREVALENCIA POR MODELO LOGISTICO #############
################################### BAYESIANO #############################


