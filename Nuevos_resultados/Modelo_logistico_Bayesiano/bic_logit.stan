// Regresión logística bayesiana con corrección interna del error de clasificación (BIC).
// Admite Se/Sp fijas o con prior Beta, y una prior normal distinta para el intercepto
// y para las pendientes. Las covariables pueden entrar centradas o sin centrar.
data {
  int<lower=1> N;
  int<lower=1> P;
  matrix[N, P] X;
  array[N] int<lower=0, upper=1> y;
  real<lower=0> sd_intercept;
  real<lower=0> sd_beta;
  int<lower=0, upper=1> estimate_sesp;   // 0 = Se/Sp fijas; 1 = prior Beta
  real<lower=0, upper=1> Se_fixed;
  real<lower=0, upper=1> Sp_fixed;
  vector<lower=0>[2] a_se;                // Beta(a_se[1], a_se[2])
  vector<lower=0>[2] a_sp;                // Beta(a_sp[1], a_sp[2])
}
parameters {
  real alpha;
  vector[P] beta;
  array[estimate_sesp] real<lower=0.5, upper=1> Se_par;
  array[estimate_sesp] real<lower=0.5, upper=1> Sp_par;
}
transformed parameters {
  real Se = estimate_sesp ? Se_par[1] : Se_fixed;
  real Sp = estimate_sesp ? Sp_par[1] : Sp_fixed;
}
model {
  vector[N] p_true = inv_logit(alpha + X * beta);
  alpha ~ normal(0, sd_intercept);
  beta ~ normal(0, sd_beta);
  if (estimate_sesp) {
    Se_par[1] ~ beta(a_se[1], a_se[2]);
    Sp_par[1] ~ beta(a_sp[1], a_sp[2]);
  }
  y ~ bernoulli(Se * p_true + (1 - Sp) * (1 - p_true));
}
