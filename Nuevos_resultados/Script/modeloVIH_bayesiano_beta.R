pacman::p_load(rstan, rio, marginaleffects, tidyverse)

# datos
vih <- readRDS("vih.rds")

# Definir variables
N <- nrow(vih)
t1 <- vih$vih_1
x1 <- vih$edad_gr
x2 <- vih$sex
x3 <- vih$reactivo_sif1
x4 <- vih$reactivo_hepb
x5 <- vih$tip_pob_HSH
x6 <- vih$tip_pob_LGTBIQ
x7 <- vih$`tip_pob_OTRAS POBLACIONES`
x8 <- vih$`tip_pob_TRABAJADOR SEXUAL`

# sensibilidad y especificidad del test
Se <- 0.9751
Sp <- 0.999934

# Parámetros del gold standard y test empleado
nplus <- 682   # detectados by gold standard method
tplus <- 665  # numero de individuos positivos entre los confirmados por gold
nminus <- 2469 # negativos by gold standard method
tminus <- 2468 # negativos entre descartados por gold

set.seed(123)

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
  nplus = nplus, 
  tplus = tplus,
  nminus = nminus,
  tminus = tminus)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

stanc(file = "Nuevos_resultados/Modelo_logistico_Bayesiano/Model_bayesiano2.stan")
# Especificar la ruta al archivo .stan
stan_file <- "Nuevos_resultados/Modelo_logistico_Bayesiano/Model_bayesiano2.stan"
#COMPILAR MODELO 
vih_bayesiano_ajustBeta <- stan_model(file = stan_file, auto_write = TRUE)

vih_bayesiano_ajustBeta <- sampling(stan_model, data = stan_data, chains = 7,
                                 iter = 6000)


# Ver resultados
print(vih_bayesiano_ajustBeta, pars = c("betas", "SN", "SP"))
