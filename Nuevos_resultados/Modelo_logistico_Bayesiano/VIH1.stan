data {
  int<lower=0> N;  // Numero de observaciones
  int<lower=0> P;  // Numero de covariables
  matrix[N, P] X;  // Matriz de covariables
  int<lower=0, upper=1> y_obs[N];  // Variable observada (con error)
  real<lower=0, upper=1> Se;  // Sensibilidad
  real<lower=0, upper=1> Sp;  // Especificidad
}

parameters {
  vector[P] beta;  // Coeficientes de regresion
  real alpha;  // Intercepto
}

transformed parameters {
  vector[N] p_true; // Probabilidad de la variable latente
  for (i in 1:N) {
    p_true[i] = inv_logit(alpha + dot_product(X[i], beta));
  }
}

model {
  // Priores
  alpha ~ normal(0,0.365541);
  beta ~ normal(0,0.365541);

  // Verosimilitud ajustada por error de medicion
  for (i in 1:N) {
    real p_obs_1 = Se * p_true[i] + (1 - Sp) * (1 - p_true[i]);
    y_obs[i] ~ bernoulli(p_obs_1);
  }
}
