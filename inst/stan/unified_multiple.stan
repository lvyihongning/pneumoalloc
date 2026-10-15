// =====================================================================
// pneumoalloc: unified Bayesian framework among children with multiple
// serotype carriage
// inst/stan/unified_multiple.stan
//
// Identical (apart from this header) to Multi_M1to4_unified_20261002.stan,
// except that m_sh (allocation categories) may exceed m_obs: m_obs is the
// number of detected serotypes capped at 3 ("3 or more"), so the check
// m_sh <= m_obs was removed, and the commune random effects of the
// published file were removed (pneumoalloc fits no cluster effects).
// Neither change alters the model fitted in the analysis (use_re = 0).
// Population: children with >= 2 detected serotypes (m = 2 or 3).
// There is no carriage component: every child here is a multiple
// carrier, so the estimand is the composition of multiple carriage,
//   share_s(t), with sum_s share_s(t) = 1.
//
// Three components, fitted jointly:
//   (b) serotype composition : fractional multinomial logit,
//                              reference = last category
//   (c) density              : lognormal, shift for m = 3 vs m = 2
//   (d) multiplicity         : Bernoulli logit for m = 3 vs m = 2, so
//                              that E[D | X, t] can be marginalised
//
// As in the primary model, the four approaches differ ONLY in the data
// passed to (b):
//
//   approach              within-child allocation q_ks    use_dens_weight
//   i   dominant          one-hot on the dominant serotype        0
//   ii  proportional      1 / m_i for each detected serotype      0
//   iii relative abund.   relative abundance theta_is             0
//   iv  density           relative abundance theta_is             1
//
// Child weight in (b): h_k = 1 (i-iii) or D_k / Dbar (iv), and the same
// weight in the standardisation (1 or E[D_i | X_i, t]).
//
// (c) and (d) share no parameters with (b), so meanD, P_m3 and
// dens_fold_3vs2 have the same posterior under all four approaches (up
// to Monte Carlo error); only the composition differs.
// =====================================================================

data {
  int<lower=1> T;                       // survey rounds
  int<lower=2> C;                       // serotype categories
  int<lower=1> N;                       // children with multiple carriage
  int<lower=1> K;                       // design matrix columns
  int<lower=1, upper=T> baseline;

  array[C] int<lower=0, upper=1> vt_mask_sero;

  matrix[N, K] X;                       // first T columns = survey indicators
  array[T] matrix[N, K] X_pred;         // counterfactual designs

  array[N] int<lower=2, upper=3> m_obs; // detected serotypes

  // Within-child allocation, sparse; active sh_w entries sum to 1.
  int<lower=1> m_sh_max;
  array[N] int<lower=1, upper=m_sh_max> m_sh;
  array[N, m_sh_max] int<lower=1, upper=C> sh_idx;
  array[N, m_sh_max] real<lower=0> sh_w;

  vector[N] logD;                       // log(dens), natural log
  vector<lower=0>[N] dens;

  // 0 = approaches i-iii (h_k = 1); 1 = approach iv (h_k = D_k / Dbar)
  int<lower=0, upper=1> use_dens_weight;

  // Prior scales (defaults as in the primary model)
  real<lower=0> prior_sd_share;         // [2]
  real<lower=0> prior_sd_mult;          // [2]
  real<lower=0> prior_sd_dens_survey;   // [3]
  real<lower=0> prior_sd_dens_cov;      // [1]
  real<lower=0> prior_sd_gamma;         // [1]
  real<lower=0> prior_rate_sigma_dens;  // [0.5] exponential rate
  real dens_prior_center;               // [log(1e6)]

}

transformed data {
  vector[N] w_dens;
  vector[N] m3_vec;
  array[N] int m3;
  vector[C] vt_mask_vec;

  // Dbar is computed here, so sum(w_dens) = N under every approach.
  real dens_scale = mean(dens);
  real kish_ess_dens_weights;

  if (K < T) reject("X must start with the T survey indicator columns");

  for (i in 1:N) {
    w_dens[i] = use_dens_weight ? dens[i] / dens_scale : 1.0;
    m3[i]     = m_obs[i] == 3;
    m3_vec[i] = m3[i];
  }
  for (c in 1:C) vt_mask_vec[c] = vt_mask_sero[c];

  // Kish effective sample size of the child weights; N when unweighted.
  kish_ess_dens_weights = square(sum(w_dens)) / dot_self(w_dens);

  for (i in 1:N) {
    real sw = 0;
    for (j in 1:m_sh[i]) sw += sh_w[i, j];
    if (abs(sw - 1.0) > 1e-8)
      reject("Active sh_w values do not sum to 1 at child i = ", i);
  }

  print("use_dens_weight = ", use_dens_weight, " ; Dbar = ", dens_scale,
        " ; Kish effective sample size = ", kish_ess_dens_weights, " of ", N);
}

parameters {
  // ---- (b) Serotype composition; category C is the reference ----
  matrix[K, C-1] beta_share_raw;

  // ---- (c) Density ----
  vector[K] beta_dens;
  real<lower=1e-6, upper=10> sigma_dens;
  real gamma_ge3;                       // log fold-change, m = 3 vs m = 2

  // ---- (d) Multiplicity: m = 3 vs m = 2 ----
  vector[K] beta_m3;
}

transformed parameters {
  matrix[K, C] beta_share = append_col(beta_share_raw, rep_vector(0.0, K));
}

model {
  // ---------------- Priors ----------------
  to_vector(beta_share_raw) ~ normal(0, prior_sd_share);
  beta_m3 ~ normal(0, prior_sd_mult);

  // Survey coefficients are log-density levels; the rest are slopes.
  beta_dens[1:T] ~ normal(dens_prior_center, prior_sd_dens_survey);
  if (K > T) beta_dens[(T + 1):K] ~ normal(0, prior_sd_dens_cov);
  sigma_dens ~ exponential(prior_rate_sigma_dens);
  gamma_ge3 ~ normal(0, prior_sd_gamma);


  // ---------------- Likelihood ----------------
  // (c) Density (natural-log lognormal)
  logD ~ normal(X * beta_dens + gamma_ge3 * m3_vec,
                sigma_dens);

  // (d) Multiplicity
  m3 ~ bernoulli_logit(X * beta_m3);

  // (b) Composition: lp += h_k * sum_s q_ks * log rho_ks.
  //     w_dens is the OBSERVED D_k / Dbar for approach iv.
  {
    matrix[N, C] eta = X * beta_share;
    for (i in 1:N) {
      real lse = log_sum_exp(eta[i]);
      real acc = 0;
      for (j in 1:m_sh[i]) acc += sh_w[i, j] * (eta[i, sh_idx[i, j]] - lse);
      target += w_dens[i] * acc;
    }
  }
}

generated quantities {
  real dens_weight_kish_ess = kish_ess_dens_weights;
  real dens_fold_3vs2 = exp(gamma_ge3);

  // ---- Composition of multiple carriage (standardised) ----
  array[T, C] real<lower=0, upper=1> share;
  array[T] real<lower=0, upper=1> share_VT;
  array[T] real<lower=0, upper=1> share_NVT;

  array[T, C] real<lower=0> PR_sero;    // ratio vs baseline
  array[T] real<lower=0> PR_VT;
  array[T] real<lower=0> PR_NVT;
  array[T, C] real PD_sero;             // difference vs baseline

  // ---- Density and multiplicity among multiple carriers ----
  array[T] real<lower=0> meanD;         // marginal over multiplicity
  array[T] real<lower=0> PR_meanD;
  array[T] real<lower=0, upper=1> P_m3; // proportion with 3 serotypes

  {
    real eps = 1e-12;
    real half_v = 0.5 * square(sigma_dens);

    for (t in 1:T) {
      matrix[N, C] eta_t = X_pred[t] * beta_share;
      vector[N] p3 = inv_logit(X_pred[t] * beta_m3);
      // E[D | X_i, t], marginalised over m = 2 / 3
      vector[N] Dhat = exp(X_pred[t] * beta_dens + half_v)
                       .* (1 - p3 + p3 * exp(gamma_ge3));

      vector[C] numer = rep_vector(0.0, C);
      real denom = 0;
      for (i in 1:N) {
        // child weight: 1 (i-iii) or E[D_i] (iv)
        real wgt = use_dens_weight ? Dhat[i] : 1.0;
        numer += wgt * softmax(eta_t[i]');
        denom += wgt;
      }

      vector[C] sv = numer / fmax(denom, eps);
      for (c in 1:C) share[t, c] = sv[c];
      share_VT[t]  = fmin(1.0, fmax(0.0, dot_product(sv, vt_mask_vec)));
      share_NVT[t] = 1.0 - share_VT[t];

      meanD[t] = mean(Dhat);
      P_m3[t]  = mean(p3);
    }

    for (t in 1:T) {
      PR_VT[t]    = share_VT[t]  / fmax(share_VT[baseline],  eps);
      PR_NVT[t]   = share_NVT[t] / fmax(share_NVT[baseline], eps);
      PR_meanD[t] = meanD[t] / fmax(meanD[baseline], eps);
      for (c in 1:C) {
        PR_sero[t, c] = share[t, c] / fmax(share[baseline, c], eps);
        PD_sero[t, c] = share[t, c] - share[baseline, c];
      }
    }
  }
}
