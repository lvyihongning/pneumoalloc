// =====================================================================
// pneumoalloc: presence-based carriage
// inst/stan/presence.stan
//
// Identical (apart from this header) to Sensi_presence_based_20261002.stan
// of the PCV10 Nha Trang analysis. Outcomes are ordered as overall
// carriage, VT, NVT, multiple carriage, then each serotype category.
// Presence-based sensitivity analysis.
//
// Every outcome is a child-level binary indicator with its own logistic
// regression on the same design matrix:
//   y[c, n] = 1 if outcome c is present in child n
// The R script orders the outcomes as overall carriage, VT, NVT,
// multiple serotype carriage, then each retained serotype.
//
// A child carrying 1 VT + 2 NVT serotypes is positive for overall, VT,
// NVT and multiple carriage, and for each of the three serotypes. VT and
// NVT prevalences are therefore NOT constrained to sum to overall
// carriage, unlike the allocation approaches.
//
// Prevalences are marginally standardised to the observed covariate
// distribution: X_pred[t] sets every child to survey round t.

data {
  int<lower=1> T;                          // survey rounds
  int<lower=1> N;                          // children
  int<lower=1> K;                          // design-matrix columns
  int<lower=1> C;                          // binary outcomes

  array[C, N] int<lower=0, upper=1> y;
  matrix[N, K] X;                          // ~ 0 + survey (+ covariates)
  array[T] matrix[N, K] X_pred;            // counterfactual designs

  int<lower=1, upper=T> baseline;
  real<lower=0> prior_sd;                  // [2.5], as prior_sd_car in M1to4
}

parameters {
  matrix[K, C] beta;
}

model {
  to_vector(beta) ~ normal(0, prior_sd);

  for (c in 1:C)
    y[c] ~ bernoulli_logit_glm(X, 0.0, beta[, c]);
}

generated quantities {
  array[T, C] real<lower=0, upper=1> P;    // standardised prevalence
  array[T, C] real<lower=0> PR;            // prevalence ratio vs baseline

  for (t in 1:T) {
    matrix[N, C] eta = X_pred[t] * beta;
    for (c in 1:C)
      P[t, c] = mean(inv_logit(eta[, c]));
  }

  for (t in 1:T)
    for (c in 1:C)
      PR[t, c] = P[t, c] / P[baseline, c];
}
