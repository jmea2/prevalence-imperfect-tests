data {
  int<lower=1> N;       // Número de observaciones
  vector[N] x1;         // Covariable 1
  vector[N] x2;
  vector[N] x3;
  vector[N] x4;
  vector[N] x5;
  vector[N] x6;
  vector[N] x7;
  vector[N] x8;
  int<lower=0,upper=1> t1[N]; // Resultado del test diagnóstico
  real<lower=0,upper=1> SN;   // Sensibilidad fija
  real<lower=0,upper=1> SP;   // Especificidad fija
}

parameters {
  vector[9] betas; // Coeficientes de regresión logística
}

model {
  vector[N] logit_p1;
  vector[N] p1;
  vector[N] p2;

  // Priores para los betas
  betas ~ normal(0, sqrt(0.366));

  // Modelo de probabilidad de enfermedad
  for (i in 1:N) {
    logit_p1[i] = betas[1] + betas[2] * x1[i] + betas[3] * x2[i] + betas[4] * x3[i] + betas[5] * x4[i] + betas[6] * x5[i] + betas[7] * x6[i] + betas[8] * x7[i] + betas[9] * x8[i];
    p1[i] = inv_logit(logit_p1[i]); // Probabilidad de estar enfermo

    // Modelamos el estado latente `d1[i]` implícitamente como una variable oculta
    p2[i] = p1[i] * SN + (1 - p1[i]) * (1 - SP);

    // Modelo de observación: test diagnóstico
    t1[i] ~ bernoulli(p2[i]);
  }
}
