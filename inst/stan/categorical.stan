// =====================================================================
// pneumoalloc: multinomial-logit model with posterior standardisation
// inst/stan/categorical.stan
//
// Used for
//   * carriage states (no carriage, 1, 2, 3, ... serotypes),
//     reference = category 1 (no carriage);
//   * types of multiple carriage (PCV10-only, non-PCV10-only, mixed
//     PCV10-dominant, mixed non-PCV10-dominant), reference = last.
//
// Generalises 1_Co-colonisation_20261002.stan (G = 4, ref = 1) and
// Proportion in co-colonisation_20261002.stan (G = 4, ref = 4) to G
// categories. The parameters (beta_raw, K x (G - 1), the non-reference
// columns in order), priors and likelihood are the same as in those
// files.
// =====================================================================

data {
  int<lower=1> N;                       // children
  int<lower=1> K;                       // design matrix columns
  int<lower=1> T;                       // survey rounds
  int<lower=2> G;                       // outcome categories
  int<lower=1, upper=G> ref;            // reference category
  real<lower=0> prior_sd;               // normal(0, prior_sd) on all coefficients

  array[N] int<lower=1, upper=G> y;     // outcome category
  matrix[N, K] X;                       // first T columns = survey indicators
  array[T] matrix[N, K] X_pred;         // counterfactual designs
}

transformed data {
  if (K < T) reject("X must start with the T survey indicator columns");
}

parameters {
  // Coefficients of the non-reference categories. The design has one
  // column per survey round (~ 0 + survey + ...), so there is no separate
  // intercept.
  matrix[K, G - 1] beta_raw;
}

transformed parameters {
  vector[G] alpha = rep_vector(0.0, G);
  matrix[K, G] beta;
  {
    int j = 1;
    for (g in 1:G) {
      if (g == ref) {
        beta[:, g] = rep_vector(0.0, K);
      } else {
        beta[:, g] = beta_raw[:, j];
        j += 1;
      }
    }
  }
}

model {
  to_vector(beta_raw) ~ normal(0, prior_sd);
  y ~ categorical_logit_glm(X, alpha, beta);
}

generated quantities {
  // Population-standardised probability of each category by survey:
  // every child is set to survey t and the predicted probabilities are
  // averaged over children (g-computation).
  array[T] vector[G] P_class;

  for (t in 1:T) {
    matrix[N, G] logits = rep_matrix(alpha', N) + X_pred[t] * beta;
    vector[G] acc = rep_vector(0.0, G);
    for (n in 1:N) acc += softmax(logits[n]');
    P_class[t] = acc / N;
  }
}
