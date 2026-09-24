############################################################################
# REANÁLISIS PARA LA REVISIÓN DEL MANUSCRITO (SMMR)
#
# Se corre desde la raíz del proyecto (donde está vih.rds), igual que los
# demás scripts de Nuevos_resultados/Script:
#   Rscript Nuevos_resultados/Script/reanalisis_revision.R vih
#   Rscript Nuevos_resultados/Script/reanalisis_revision.R sif
#
# Cambios frente a los scripts originales:
#  1. Rogan–Gladen con el paréntesis corregido: (p - (1 - Sp)) / (Se + Sp - 1).
#  2. Liu ajustado por máxima verosimilitud propia (logistic4p ya no está en
#     CRAN): versión "equal" sin restringir (reproduce lo publicado), "equal"
#     restringida a r en [0, 0.5) y "fp.fn" restringida.
#  3. Intervalos de Liu y de los modelos bayesianos obtenidos de su propia
#     incertidumbre (Hessiana / posterior), no del glm "dummy".
#  4. BIC con cuatro especificaciones de priors (análisis de sensibilidad).
#  5. Diagnósticos de convergencia para todos los parámetros.
#  6. Prevalencia para el perfil medio de covariables (como en el artículo) y
#     prevalencia marginal estandarizada (promedio sobre la muestra).
############################################################################

pacman::p_load(tidyverse, rstan, rstanarm, posterior, numDeriv, MASS)

args <- commandArgs(trailingOnly = TRUE)
outcome <- if (length(args)) args[1] else "vih"

################################ CONFIGURACIÓN ##############################

covs_base <- c("edad_gr", "sex", "reactivo_hepb", "tip_pob_HSH", "tip_pob_LGTBIQ",
               "tip_pob_OTRAS POBLACIONES", "tip_pob_TRABAJADOR SEXUAL")

cfg <- switch(outcome,
  vih = list(
    y = "vih_1",
    covs = append(covs_base, "reactivo_sif1", after = 2),
    # algoritmo serial Determine -> Bioline (límites inferiores IC95 %)
    Se = 0.995 * 0.980,
    Sp = 0.994 + (1 - 0.994) * 0.989,
    # validación frente a estándar de oro (Model_bayesiano2.stan)
    a_se = c(665 + 1, 682 - 665 + 1),
    a_sp = c(2468 + 1, 2469 - 2468 + 1)
  ),
  sif = list(
    y = "reactivo_sif1",
    covs = append(covs_base, "vih_1", after = 2),
    Se = 0.964,
    Sp = 0.974,
    # sin conteos de validación disponibles: Beta con la media del fabricante
    # y un tamaño efectivo de 500 (supuesto, declararlo en el texto)
    a_se = 500 * c(0.964, 1 - 0.964),
    a_sp = 500 * c(0.974, 1 - 0.974)
  ),
  stop("El argumento debe ser 'vih' o 'sif'")
)

sd_newman <- sqrt(pi^2 / (3 * (length(cfg$covs) + 1)))   # varianza 0.366
iter <- 2000; chains <- 4; seed <- 2026
reusar <- TRUE   # TRUE: carga los ajustes Stan ya guardados en Modelos/revision

dir_res <- "Nuevos_resultados/Resultados/revision"
dir_mod <- "Nuevos_resultados/Modelos/revision"
dir.create(dir_res, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_mod, showWarnings = FALSE, recursive = TRUE)

options(mc.cores = min(chains, parallel::detectCores()))
rstan_options(auto_write = TRUE)
set.seed(seed)

################################ DATOS ######################################

vih <- readRDS("vih.rds")
y <- vih[[cfg$y]]
X <- as.matrix(vih[, cfg$covs])
N <- nrow(X); P <- ncol(X)
x_mean <- colMeans(X)                      # perfil medio (datagrid del original)
coef_names <- c("(Intercept)", cfg$covs)

rg <- function(p, Se, Sp) (p - (1 - Sp)) / (Se + Sp - 1)
# estimación puntual: media posterior (bayesianos) o valor en el MLE (point)
qs <- function(x, point = mean(x)) c(estimate = point,
                                     lower = unname(quantile(x, .025)),
                                     upper = unname(quantile(x, .975)))
# prevalencia (perfil medio y marginal) a partir de una matriz de draws [alpha, beta]
prev_draws <- function(B) {
  Xa <- cbind(1, X)
  list(perfil = plogis(cbind(1, t(x_mean)) %*% t(B))[1, ],
       marginal = colMeans(plogis(Xa %*% t(B))))
}

prev_tab <- list(); coef_tab <- list(); diag_tab <- list()
add_prev <- function(model, tipo, v) {
  prev_tab[[length(prev_tab) + 1]] <<- tibble(model, tipo, !!!as.list(v))
}

############################ 1. PREVALENCIA CRUDA ###########################
# Rogan–Gladen con intervalo Agresti–Coull (IC_agresti_se_sp.R)

z <- qnorm(.975); x <- sum(y)
AP <- (x + z^2 / 2) / (N + z^2)
p_rg <- max(rg(AP, cfg$Se, cfg$Sp), 0)
ee <- z / (cfg$Se + cfg$Sp - 1) * sqrt(AP * (1 - AP) / (N + z^2))
add_prev("Observada", "cruda", c(estimate = x / N,
                                 lower = binom::binom.agresti.coull(x, N)$lower,
                                 upper = binom::binom.agresti.coull(x, N)$upper))
add_prev("Crude (RG)", "cruda", c(estimate = p_rg, lower = max(p_rg - ee, 0),
                                  upper = p_rg + ee))

############################ 2. MODELO ESTÁNDAR (STD) ########################

d <- data.frame(y = y, X, check.names = FALSE)
mod_std <- glm(y ~ ., family = binomial, data = d)
b_std <- coef(mod_std)
B_std <- MASS::mvrnorm(4000, b_std, vcov(mod_std))
pd <- prev_draws(B_std)
p0 <- prev_draws(t(b_std))
for (t in names(pd)) {
  add_prev("STD", t, qs(rg(pd[[t]], cfg$Se, cfg$Sp), rg(p0[[t]], cfg$Se, cfg$Sp)))
  add_prev("STD sin corregir", t, qs(pd[[t]], p0[[t]]))
}
coef_tab$STD <- tibble(term = coef_names, estimate = b_std,
                       se = sqrt(diag(vcov(mod_std))))

############################ 3. MODELO DE LIU ###############################
# P(Y=1|X) = r0 + (1 - r0 - r1) * expit(eta)
# tipo "libre": r sin restringir (lo que hace logistic4p, admite r < 0)
# tipo "restr": r = 0.5 * expit(theta), es decir r en (0, 0.5)

liu_fit <- function(model = c("equal", "fp.fn"), restr = TRUE) {
  model <- match.arg(model)
  k <- if (model == "equal") 1 else 2
  Xa <- cbind(1, X)
  tr  <- if (restr) function(t) 0.5 * plogis(t) else identity
  dtr <- if (restr) function(t) 0.5 * plogis(t) * (1 - plogis(t)) else function(t) 1
  nll <- function(par) {
    r <- tr(par[1:k]); r0 <- r[1]; r1 <- r[k]
    p <- r0 + (1 - r0 - r1) * plogis(drop(Xa %*% par[-(1:k)]))
    if (any(p <= 0 | p >= 1)) return(Inf)
    -sum(dbinom(y, 1, p, log = TRUE))
  }
  grad <- function(par) {
    th <- par[1:k]; r <- tr(th); r0 <- r[1]; r1 <- r[k]
    pi <- plogis(drop(Xa %*% par[-(1:k)])); p <- r0 + (1 - r0 - r1) * pi
    w <- y / p - (1 - y) / (1 - p)
    g_r <- if (k == 1) sum(w * (1 - 2 * pi)) * dtr(th) else
      c(sum(w * (1 - pi)) * dtr(th[1]), sum(-w * pi) * dtr(th[2]))
    -c(g_r, colSums(w * (1 - r0 - r1) * pi * (1 - pi) * Xa))
  }
  # varios valores iniciales de r (la verosimilitud es multimodal); se
  # conserva la solución con mayor log-verosimilitud finita
  r_ini <- if (restr) qlogis(2 * c(0.001, 0.01, 0.05)) else c(0.01, 0, -0.005)
  fits <- lapply(r_ini, function(r)
    nlminb(c(rep(r, k), b_std), nll, grad,
           control = list(iter.max = 2000, eval.max = 4000)))
  op <- fits[[which.min(sapply(fits, `[[`, "objective"))]]
  H <- numDeriv::jacobian(grad, op$par)
  H <- (H + t(H)) / 2
  V <- tryCatch({
    if (any(eigen(H, only.values = TRUE)$values <= 0)) stop("no def. pos.")
    solve(H)
  }, error = function(e) NULL)
  pd_hessian <- !is.null(V)
  if (!pd_hessian) {
    # r en el borde / no identificada: varianza de beta condicionada a r estimada
    V <- matrix(0, length(op$par), length(op$par))
    Hb <- H[-(1:k), -(1:k)]
    V[-(1:k), -(1:k)] <- tryCatch(solve(Hb), error = function(e) MASS::ginv(Hb))
  }
  list(par = op$par, V = V, loglik = -op$objective, conv = op$message,
       r = tr(op$par[1:k]), k = k, tr = tr, pd_hessian = pd_hessian)
}

liu <- list(`Liu equal (libre)` = liu_fit("equal", FALSE),
            `Liu equal (restr.)` = liu_fit("equal", TRUE),
            `Liu fp.fn (restr.)` = liu_fit("fp.fn", TRUE))

liu_resumen <- imap_dfr(liu, function(f, nm) {
  tibble(model = nm, r0 = f$r[1], r1 = f$r[f$k], loglik = f$loglik,
         convergencia = f$conv, hessiana_def_pos = f$pd_hessian)
})

for (nm in names(liu)) {
  f <- liu[[nm]]
  B <- MASS::mvrnorm(4000, f$par, f$V)[, -(1:f$k)]
  pd <- prev_draws(B)
  p0 <- prev_draws(t(f$par[-(1:f$k)]))
  for (t in names(pd)) add_prev(nm, t, qs(pd[[t]], p0[[t]]))
  coef_tab[[nm]] <- tibble(term = coef_names, estimate = f$par[-(1:f$k)],
                           se = sqrt(pmax(diag(f$V), 0))[-(1:f$k)])
}

#################### 4. BAYESIANO CON CORRECCIÓN EXTERNA (BEC) ##############

f_bec <- file.path(dir_mod, paste0("bec_", outcome, ".rds"))
if (reusar && file.exists(f_bec)) {
  mod_bec <- readRDS(f_bec)
} else {
  mod_bec <- stan_glm(y ~ ., data = d, family = binomial(link = "logit"),
                      prior = normal(0, 0.6, autoscale = FALSE),
                      prior_intercept = normal(0, 0.6, autoscale = FALSE),
                      chains = chains, iter = iter, seed = seed, refresh = 0)
  saveRDS(mod_bec, f_bec)
}
B <- as.matrix(mod_bec)[, 1:(P + 1)]
pd <- prev_draws(B)
for (t in names(pd)) add_prev("BEC", t, qs(rg(pd[[t]], cfg$Se, cfg$Sp)))
coef_tab$BEC <- tibble(term = coef_names, estimate = colMeans(B),
                       se = apply(B, 2, sd))
diag_tab$BEC <- summarise_draws(as_draws_array(as.array(mod_bec)),
                                "mean", "sd", "rhat", "ess_bulk", "ess_tail")

#################### 5. BAYESIANO CON CORRECCIÓN INTERNA (BIC) ##############

stan_bic <- stan_model("Nuevos_resultados/Modelo_logistico_Bayesiano/bic_logit.stan")

escenarios <- list(
  # prior del artículo: N(0, 0.366) en intercepto y pendientes, X sin centrar
  `BIC original` = list(center = FALSE, sd_int = sd_newman, sd_b = sd_newman, sesp = 0),
  # Newman en pendientes; intercepto (X centrada) con prior débil
  `BIC intercepto débil` = list(center = TRUE, sd_int = 2.5, sd_b = sd_newman, sesp = 0),
  # prior vaga en todos los coeficientes
  `BIC vaga` = list(center = FALSE, sd_int = sqrt(10), sd_b = sqrt(10), sesp = 0),
  # como "intercepto débil", con Se y Sp inciertas (prior Beta)
  `BIC Se/Sp inciertas` = list(center = TRUE, sd_int = 2.5, sd_b = sd_newman, sesp = 1)
)

bic_sampler <- list()
for (nm in names(escenarios)) {
  e <- escenarios[[nm]]
  Xs <- if (e$center) sweep(X, 2, x_mean) else X
  f_bic <- file.path(dir_mod, paste0("bic_", outcome, "_", make.names(nm), ".rds"))
  if (reusar && file.exists(f_bic)) {
    fit <- readRDS(f_bic)
  } else {
    fit <- sampling(stan_bic, seed = seed, chains = chains, iter = iter, refresh = 0,
                    data = list(N = N, P = P, X = Xs, y = y,
                                sd_intercept = e$sd_int, sd_beta = e$sd_b,
                                estimate_sesp = e$sesp,
                                Se_fixed = cfg$Se, Sp_fixed = cfg$Sp,
                                a_se = cfg$a_se, a_sp = cfg$a_sp))
    saveRDS(fit, f_bic)
  }
  dr <- rstan::extract(fit, pars = c("alpha", "beta"))
  alpha <- dr$alpha
  if (e$center) alpha <- alpha - drop(dr$beta %*% x_mean)  # a escala original
  B <- cbind(alpha, dr$beta)
  pd <- prev_draws(B)
  for (t in names(pd)) add_prev(nm, t, qs(pd[[t]]))
  coef_tab[[nm]] <- tibble(term = coef_names, estimate = colMeans(B),
                           se = apply(B, 2, sd))
  pars <- c("alpha", "beta", if (e$sesp) c("Se", "Sp"), "lp__")
  diag_tab[[nm]] <- summarise_draws(as_draws_array(as.array(fit, pars = pars)),
                                    "mean", "sd", "rhat", "ess_bulk", "ess_tail")
  bic_sampler[[nm]] <- tibble(model = nm,
    divergencias = get_num_divergent(fit),
    max_treedepth = get_num_max_treedepth(fit),
    min_ebfmi = min(get_bfmi(fit)))
}

################################ SALIDAS ####################################

prevalencias <- bind_rows(prev_tab) %>%
  mutate(ancho_IC = upper - lower,
         cambio_vs_cruda = 100 * (estimate / estimate[model == "Crude (RG)"] - 1))
prevalencias <- prevalencias %>%
  group_by(tipo) %>%
  mutate(cambio_vs_STD = if (any(model == "STD"))
    100 * (estimate / estimate[model == "STD"] - 1) else NA_real_) %>%
  ungroup()

coeficientes <- bind_rows(coef_tab, .id = "model")
diagnosticos <- bind_rows(diag_tab, .id = "model")
sampler <- bind_rows(bic_sampler)

write_csv(prevalencias, file.path(dir_res, paste0("prevalencias_", outcome, ".csv")))
write_csv(coeficientes, file.path(dir_res, paste0("coeficientes_", outcome, ".csv")))
write_csv(diagnosticos, file.path(dir_res, paste0("diagnosticos_", outcome, ".csv")))
write_csv(sampler, file.path(dir_res, paste0("sampler_", outcome, ".csv")))
write_csv(liu_resumen, file.path(dir_res, paste0("liu_", outcome, ".csv")))

print(liu_resumen)
print(prevalencias, n = Inf, width = Inf)
print(coeficientes %>% dplyr::select(-se) %>%
        pivot_wider(names_from = model, values_from = estimate), width = Inf)
print(diagnosticos %>% group_by(model) %>%
        summarise(max_rhat = max(rhat), min_ess_bulk = min(ess_bulk),
                  min_ess_tail = min(ess_tail)))
print(sampler)
