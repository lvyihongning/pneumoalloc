// =====================================================================
// pneumoalloc: unified Bayesian framework, approaches i-iv, all children
// inst/stan/unified.stan
//
// The model of M1to4_unified_v4_20261002.stan (PCV10 Nha Trang analysis):
// same data, parameters, priors, likelihood and generated quantities, with
// a faster likelihood computation:
//   * carriage, density and multiplicity use Stan's GLM functions
//     (bernoulli_logit_glm, normal_id_glm, categorical_logit_glm), which
//     have analytic gradients;
//   * composition: with W = diag(h) A (allocation weights A, carrier
//     weights h) and rows of A summing to one,
//       sum_k h_k sum_s A_ks log softmax(eta_k)_s
//         = sum(beta_share .* X_pos' W) - sum_k h_k log_sum_exp(eta_k),
//     where eta_k = X_pos[k] * beta_share. The first term is linear in
//     beta_share and X_pos' W is computed once (transformed data). For
//     h_k = 1 (approaches i-iii), sum_k log_sum_exp(eta_k) =
//       sum_k eta_k[1] - categorical_logit_glm_lpmf(1 | X_pos, 0, beta_share).
// Its log density equals that of the published file up to floating-point
// rounding (tests/testthat/fixtures/unified_published.stan and
// test-unified_equivalence.R). The commune random effects of the published
// file were removed (pneumoalloc fits no cluster effects).
// Approach iv weights are normalised WITHIN each survey round:
//   h_k = Theta_k / Theta_bar_t, with Theta_bar_t the mean total density
// among the carriers of survey round t (the published file used one mean
// over all rounds). The weights of each round then add up to its number of
// carriers, so no round gains or loses information through the weighting.
// m_obs is the number of detected serotypes, capped at M
// (prepare_carriage_data(max_serotypes = )); m_sh, the number of
// allocation categories, may exceed it.
// =====================================================================
// One model, three components, fitted jointly:
//   (a) overall carriage        : Bernoulli logit, all N children
//   (b) serotype composition    : fractional multinomial logit among
//                                 carriers, reference = last category
//   (c) density among carriers  : lognormal, with a multiplicity
//                                 (1/2/3 serotypes) sub-model so that
//                                 E[D | X, C = 1, t] can be marginalised
//
// The four approaches differ ONLY in the data passed to (b):
//
//   approach              within-carrier allocation q_ks    use_dens_weight
//   i   dominant          one-hot on the dominant serotype        0
//   ii  proportional      1 / m_i for each detected serotype      0
//   iii relative abund.   relative abundance theta_is             0
//   iv  density           relative abundance theta_is             1
//
// Carrier weight in (b):  h_k = 1                         (i-iii)
//                         h_k = Theta_k / Theta_bar_t     (iv)
// Theta_k: total density of carrier k; theta_ks: relative abundance of
// serotype s; Theta_bar_t: mean total density among carriers of round t.
// and the same weight is used in the standardisation in generated
// quantities (psi_i for i-iii, psi_i * E[D] for iv).
//
// Because (a) and (c) share no parameters with (b), P_overall, P_mult
// and meanD_pos have the same posterior under all four approaches (up
// to Monte Carlo error); only the composition differs.
//
// Changes relative to M4_psiweighted_v3.stan:
//   1. use_dens_weight flag (data): h_k = 1 when 0, Theta_k / Theta_bar_t when 1,
//      applied in the share likelihood AND in both standardisations.
//   2. m_sh and m_obs are separated:
//        m_obs[k] = number of detected serotypes (multiplicity, density)
//        m_sh[k]  = number of non-zero allocation weights (share loop)
//      They differ under approach i (one-hot) and whenever two
//      detected serotypes collapse into the same category
//      (e.g. "Other NVT"). v3 used m_sh for both.
//   3. meanD_pos has its own accumulator, so it is the carriage-weighted
//      mean density under every approach (v3 reused the share
//      denominator, which is only correct when use_dens_weight = 1).
//
// Estimand (main text):
//   P_s(t) = P_overall(t) * share_s(t),  so sum_s P_s(t) = P_overall(t).
// =====================================================================

functions {
  // Multiplicity shift on the natural-log density scale.
  // gamma_ge2 : log fold-change for m >= 2 versus m = 1
  // gamma_ge3 : additional log fold-change for m = 3 versus m = 2
  real multiplicity_shift(int g, real gamma_ge2, real gamma_ge3) {
    real out = 0;
    if (g >= 2) out += gamma_ge2;
    if (g >= 3) out += gamma_ge3;
    return out;
  }
}

data {
  int<lower=1> T;                       // survey rounds
  int<lower=2> C;                       // serotype categories
  int<lower=1> N;                       // all sampled children
  int<lower=1> N_pos;                   // carriers
  int<lower=1> K;                       // design matrix columns
  int<lower=2, upper=5> M;              // max detected serotypes (3 here)
  int<lower=1, upper=T> baseline;

  array[C] int<lower=0, upper=1> vt_mask_sero;

  matrix[N, K] X;
  array[T] matrix[N, K] X_pred;         // counterfactual designs

  // Carriage outcome: 1 if pneumococcus detected, else 0.
  // carrier[i] == 1 iff i is in pos_idx (checked below).
  array[N] int<lower=0, upper=1> carrier;

  // Detected serotypes per carrier (multiplicity + density models).
  array[N_pos] int<lower=1, upper=M> m_obs;

  // Within-carrier allocation, sparse. Categories are 1..C and the
  // active entries of each sh_w row sum to 1.
  int<lower=1> m_sh_max;
  array[N_pos] int<lower=1, upper=m_sh_max> m_sh;  // share loop only
  array[N_pos, m_sh_max] int<lower=1, upper=C> sh_idx;
  array[N_pos, m_sh_max] real<lower=0> sh_w;

  array[N_pos] int<lower=1, upper=N> pos_idx;

  vector[N_pos] logD_pos;                 // log(dens_pos), natural log
  vector<lower=0>[N_pos] dens_pos;

  // 0 = approaches i-iii (h_k = 1); 1 = approach iv (h_k = Theta_k / Theta_bar_t)
  int<lower=0, upper=1> use_dens_weight;

  // Prior scales, supplied as data so sensitivity analyses need no
  // recompilation. Suggested defaults in brackets.
  real<lower=0> prior_sd_car;             // [2.5]
  real<lower=0> prior_sd_share;           // [2] primary; [5] sensitivity
  real<lower=0> prior_sd_mult;            // [2]
  real<lower=0> prior_sd_dens_survey;     // [3]
  real<lower=0> prior_sd_dens_cov;        // [1]
  real<lower=0> prior_sd_gamma;           // [1]
  real<lower=0> prior_rate_sigma_dens;    // [0.5] exponential rate

  real dens_prior_center;                 // [log(1e6)]

  int<lower=0, upper=1> compute_log_lik;  // 0 to save output size
}

transformed data {
  matrix[N_pos, K] X_pos;

  vector[N_pos] w_dens;
  vector[C] vt_mask_vec;

  // Theta_bar_t is computed here rather than passed in, so that the
  // weights of each round add up to its number of carriers by construction.
  real dens_scale = mean(dens_pos);   // overall mean, reported only

  real w_dens_sum;
  matrix[K, C] XtWA;
  vector[K] x_pos_colsum;
  array[N_pos] int ones_pos;
  vector[N_pos] ge2;
  vector[N_pos] ge3;
  // Kish effective sample size of the carrier weights,
  // (sum w)^2 / sum(w^2). Equals N_pos when use_dens_weight = 0.
  real kish_ess_dens_weights;
  // by survey round: sum of the carrier weights (= carriers in the round)
  // and their Kish effective sample size
  vector[T] w_sum_round = rep_vector(0.0, T);
  vector[T] kish_ess_round;

  for (k in 1:N_pos) X_pos[k] = X[pos_idx[k]];

  // Approach iv weights: Theta_k / Theta_bar_t, normalised within each
  // survey round (round = the survey indicator among the first T columns).
  {
    array[N_pos] int t_pos;
    vector[T] w2_round = rep_vector(0.0, T);
    vector[T] sum_t = rep_vector(0.0, T);
    vector[T] n_t = rep_vector(0.0, T);
    for (k in 1:N_pos) {
      t_pos[k] = 0;
      for (t in 1:T) if (X_pos[k, t] == 1) t_pos[k] = t;
      if (t_pos[k] == 0) reject("no survey indicator for carrier k = ", k);
      sum_t[t_pos[k]] += dens_pos[k];
      n_t[t_pos[k]] += 1;
    }
    for (k in 1:N_pos) {
      w_dens[k] = use_dens_weight
                  ? dens_pos[k] / (sum_t[t_pos[k]] / n_t[t_pos[k]]) : 1.0;
      w_sum_round[t_pos[k]] += w_dens[k];
      w2_round[t_pos[k]] += square(w_dens[k]);
    }
    for (t in 1:T)
      kish_ess_round[t] = w2_round[t] > 0 ? square(w_sum_round[t]) / w2_round[t] : 0;
  }
  for (c in 1:C) vt_mask_vec[c] = vt_mask_sero[c];

  w_dens_sum = sum(w_dens);
  kish_ess_dens_weights = square(w_dens_sum) / dot_self(w_dens);

  // Consistency check: carrier flag must agree with pos_idx.
  {
    array[N] int chk = rep_array(0, N);
    for (k in 1:N_pos) chk[pos_idx[k]] = 1;
    for (i in 1:N) {
      if (chk[i] != carrier[i])
        reject("carrier[] disagrees with pos_idx at i = ", i);
    }
  }

  // Allocation weights must sum to one for every carrier.
  for (k in 1:N_pos) {
    real sw = 0;
    for (j in 1:m_sh[k]) sw += sh_w[k, j];
    if (abs(sw - 1.0) > 1e-8)
      reject("Active sh_w values do not sum to 1 at carrier k = ", k,
             " (sum = ", sw, ")");
  }

  // Linear part of the composition likelihood: X_pos' W, K x C
  {
    matrix[N_pos, C] WA = rep_matrix(0.0, N_pos, C);
    for (k in 1:N_pos)
      for (j in 1:m_sh[k])
        WA[k, sh_idx[k, j]] += w_dens[k] * sh_w[k, j];
    XtWA = X_pos' * WA;
  }
  x_pos_colsum = (rep_row_vector(1.0, N_pos) * X_pos)';
  ones_pos = rep_array(1, N_pos);
  for (k in 1:N_pos) {
    ge2[k] = m_obs[k] >= 2;
    ge3[k] = m_obs[k] >= 3;
  }

  print("use_dens_weight = ", use_dens_weight,
        " ; mean density (all rounds) = ", dens_scale,
        " ; sum(w_dens) = ", w_dens_sum, " (should equal N_pos = ", N_pos, ")",
        " ; Kish effective sample size of carrier weights = ",
        kish_ess_dens_weights, " of ", N_pos, " carriers");
}

parameters {
  // ---- Carriage (Bernoulli) ----
  vector[K] beta_car;

  // ---- Density among carriers ----
  vector[K] beta_dens;
  real<lower=1e-6, upper=10> sigma_dens;
  real gamma_ge2;
  real gamma_ge3;

  // ---- Multiplicity among carriers ----
  // Category 1 (single serotype) is the reference.
  matrix[K, M-1] beta_mult_raw;

  // ---- Serotype composition among carriers ----
  // Category C is the reference.
  matrix[K, C-1] beta_share_raw;
}

transformed parameters {
  matrix[K, C] beta_share =
    append_col(beta_share_raw, rep_vector(0.0, K));

  matrix[K, M] beta_mult =
    append_col(rep_vector(0.0, K), beta_mult_raw);

  // Per-child linear predictors are deliberately NOT declared here
  // (they would be written to the CSV on every draw).
}

model {
  // ---------------- Priors ----------------
  beta_car ~ normal(0, prior_sd_car);
  to_vector(beta_share_raw) ~ normal(0, prior_sd_share);
  to_vector(beta_mult_raw)  ~ normal(0, prior_sd_mult);

  // Density design is ~ 0 + survey + covariates, so the first T
  // coefficients are survey-specific log-density levels.
  for (kk in 1:K) {
    if (kk <= T) beta_dens[kk] ~ normal(dens_prior_center, prior_sd_dens_survey);
    else         beta_dens[kk] ~ normal(0, prior_sd_dens_cov);
  }

  sigma_dens ~ exponential(prior_rate_sigma_dens);
  gamma_ge2 ~ normal(0, prior_sd_gamma);
  gamma_ge3 ~ normal(0, prior_sd_gamma);

  // ---------------- Likelihood ----------------

  // (a) Carriage among all children
  carrier ~ bernoulli_logit_glm(X, 0.0, beta_car);

  // (c) Density among carriers (natural-log lognormal)
  logD_pos ~ normal_id_glm(X_pos, gamma_ge2 * ge2 + gamma_ge3 * ge3,
                           beta_dens, sigma_dens);

  // (c) Multiplicity among carriers
  m_obs ~ categorical_logit_glm(X_pos, rep_vector(0.0, M), beta_mult);

  // (b) Serotype composition among carriers (see header)
  target += sum(beta_share .* XtWA);
  if (use_dens_weight) {
    matrix[N_pos, C] eta = X_pos * beta_share;
    vector[N_pos] lse;
    for (k in 1:N_pos) lse[k] = log_sum_exp(eta[k]);
    target += -dot_product(w_dens, lse);
  } else {
    target += -(dot_product(x_pos_colsum, beta_share[, 1])
                - categorical_logit_glm_lpmf(ones_pos | X_pos,
                                             rep_vector(0.0, C), beta_share));
  }
}

generated quantities {
  // Diagnostics carried through from transformed data
  real dens_weight_kish_ess = kish_ess_dens_weights;
  real dens_weight_sum = w_dens_sum;
  vector[T] dens_weight_sum_round = w_sum_round;
  vector[T] dens_weight_kish_ess_round = kish_ess_round;

  real dens_fold_2vs1 = exp(gamma_ge2);
  real dens_fold_3vs2 = exp(gamma_ge3);
  real dens_fold_3vs1 = exp(gamma_ge2 + gamma_ge3);

  // ---- PRIMARY: population (psi-weighted) standardisation ----
  array[T] real<lower=0, upper=1> P_overall;
  array[T, C] real<lower=0, upper=1> share;
  array[T] real<lower=0, upper=1> share_VT;
  array[T] real<lower=0, upper=1> share_NVT;
  array[T, C] real<lower=0, upper=1> P_sero;
  array[T] real<lower=0, upper=1> P_VT;
  array[T] real<lower=0, upper=1> P_NVT;

  array[T, C] real<lower=0> PR_sero;
  array[T] real<lower=0> PR_VT;
  array[T] real<lower=0> PR_NVT;
  array[T] real<lower=0> PR_overall;

  // Absolute differences versus baseline.
  array[T, C] real PD_sero;

  // ---- SECONDARY: observed-carrier standardisation ----
  array[T, C] real<lower=0, upper=1> share_oc;
  array[T, C] real<lower=0, upper=1> P_sero_oc;
  array[T] real<lower=0, upper=1> P_VT_oc;
  array[T] real<lower=0, upper=1> P_NVT_oc;

  // ---- Density summaries (marginalised over multiplicity) ----
  array[T] real<lower=0> meanD_pos;
  array[T] real<lower=0> PR_meanD_pos;
  array[T, M] real<lower=0> meanD_pos_by_m;
  array[T, M] real<lower=0> PR_meanD_pos_by_m;

  // ---- Multiplicity distribution among carriers ----
  array[T, M] real<lower=0, upper=1> P_mult;

  // ---- Optional pointwise log-likelihood ----
  // Carriage, multiplicity and density only; the share component is a
  // fractional (and, for approach iv, weighted) pseudo-likelihood.
  array[compute_log_lik ? N : 0] real log_lik_car;
  array[compute_log_lik ? N_pos : 0] real log_lik_mult;
  array[compute_log_lik ? N_pos : 0] real log_lik_dens;

  {
    real eps = 1e-12;
    real half_v = 0.5 * square(sigma_dens);

    for (t in 1:T) {
      matrix[N, C] eta_share_t = X_pred[t] * beta_share;
      matrix[N, M] eta_mult_t  = X_pred[t] * beta_mult;
      vector[N] eta_car_t      = X_pred[t] * beta_car;
      vector[N] mu_dens_base_t = X_pred[t] * beta_dens;

      real sum_psi = 0;
      real denom   = 0;     // share denominator: sum psi_i * h_i
      real numD    = 0;     // density numerator: sum psi_i * E[D_i]
      real denomD  = 0;     // sum psi_i
      vector[C] numer  = rep_vector(0.0, C);
      vector[M] numerM = rep_vector(0.0, M);
      vector[M] cf_num = rep_vector(0.0, M);

      // ---------- population standardisation over all N ----------
      for (i in 1:N) {
        real psi_i = inv_logit(eta_car_t[i]);

        // multiplicity distribution given carriage
        vector[M] pm = softmax(eta_mult_t[i]');

        // E[D | X_i, C = 1, t], marginalised over multiplicity
        real mu_base = mu_dens_base_t[i];
        real Dhat = 0;
        for (g in 1:M)
          Dhat += pm[g] * exp(mu_base
                              + multiplicity_shift(g, gamma_ge2, gamma_ge3)
                              + half_v);

        // within-stratum composition
        vector[C] q = softmax(eta_share_t[i]');

        // carrier weight: psi_i (i-iii) or psi_i * E[D_i] (iv)
        real wgt = use_dens_weight ? psi_i * Dhat : psi_i;

        sum_psi += psi_i;
        denom   += wgt;
        numer   += wgt * q;

        numD    += psi_i * Dhat;
        denomD  += psi_i;
        numerM  += psi_i * pm;

        // counterfactual mean density with all carriers at m = g
        for (g in 1:M)
          cf_num[g] += psi_i * exp(mu_base
                                   + multiplicity_shift(g, gamma_ge2, gamma_ge3)
                                   + half_v);
      }

      P_overall[t] = sum_psi / N;

      {
        vector[C] sv = numer / fmax(denom, eps);
        sv = fmax(sv, 0.0);
        sv = sv / fmax(sum(sv), eps);

        for (c in 1:C) {
          share[t, c]  = sv[c];
          P_sero[t, c] = P_overall[t] * sv[c];
        }

        real vt_sum = dot_product(sv, vt_mask_vec);
        share_VT[t]  = fmin(1.0, fmax(0.0, vt_sum));
        share_NVT[t] = 1.0 - share_VT[t];
        P_VT[t]  = P_overall[t] * share_VT[t];
        P_NVT[t] = P_overall[t] * share_NVT[t];
      }

      // mean density among carriers, marginal over multiplicity
      meanD_pos[t] = numD / fmax(denomD, eps);

      for (g in 1:M) {
        P_mult[t, g] = numerM[g] / fmax(denomD, eps);
        meanD_pos_by_m[t, g] = cf_num[g] / fmax(denomD, eps);
      }

      // ---------- secondary: observed-carrier standardisation ----------
      {
        real denom_oc = 0;
        vector[C] numer_oc = rep_vector(0.0, C);

        for (k in 1:N_pos) {
          int i  = pos_idx[k];

          // carrier weight at observed multiplicity: 1 (i-iii) or
          // predicted density (iv)
          real w = 1.0;
          if (use_dens_weight)
            w = exp(mu_dens_base_t[i]
                    + multiplicity_shift(m_obs[k], gamma_ge2, gamma_ge3)
                    + half_v);

          vector[C] q = softmax(eta_share_t[i]');

          denom_oc += w;
          numer_oc += w * q;
        }

        vector[C] sv_oc = numer_oc / fmax(denom_oc, eps);
        sv_oc = fmax(sv_oc, 0.0);
        sv_oc = sv_oc / fmax(sum(sv_oc), eps);

        real vt_oc = dot_product(sv_oc, vt_mask_vec);

        for (c in 1:C) {
          share_oc[t, c]  = sv_oc[c];
          P_sero_oc[t, c] = P_overall[t] * sv_oc[c];
        }
        P_VT_oc[t]  = P_overall[t] * fmin(1.0, fmax(0.0, vt_oc));
        P_NVT_oc[t] = P_overall[t] * (1.0 - fmin(1.0, fmax(0.0, vt_oc)));
      }
    }

    // ---------------- ratios and differences ----------------
    {
      real b_vt  = fmax(P_VT[baseline], eps);
      real b_nvt = fmax(P_NVT[baseline], eps);
      real b_ov  = fmax(P_overall[baseline], eps);
      real b_dn  = fmax(meanD_pos[baseline], eps);

      vector[C] b_sero;
      for (c in 1:C) b_sero[c] = fmax(P_sero[baseline, c], eps);

      for (t in 1:T) {
        PR_VT[t]      = P_VT[t] / b_vt;
        PR_NVT[t]     = P_NVT[t] / b_nvt;
        PR_overall[t] = P_overall[t] / b_ov;
        PR_meanD_pos[t] = meanD_pos[t] / b_dn;

        for (c in 1:C) {
          PR_sero[t, c] = P_sero[t, c] / b_sero[c];
          PD_sero[t, c] = P_sero[t, c] - P_sero[baseline, c];
        }

        for (g in 1:M)
          PR_meanD_pos_by_m[t, g] =
            meanD_pos_by_m[t, g]
            / fmax(meanD_pos_by_m[baseline, g], eps);
      }
    }

    // ---------------- pointwise log-likelihood ----------------
    if (compute_log_lik) {
      vector[N] eta_car_obs = X * beta_car;
      vector[N_pos] mu_dens_obs = X_pos * beta_dens;

      for (k in 1:N_pos)
        mu_dens_obs[k] += multiplicity_shift(m_obs[k], gamma_ge2, gamma_ge3);

      for (i in 1:N)
        log_lik_car[i] = bernoulli_logit_lpmf(carrier[i] | eta_car_obs[i]);

      for (k in 1:N_pos) {
        log_lik_mult[k] =
          categorical_logit_lpmf(m_obs[k] | (X_pos[k] * beta_mult)');
      }

      for (k in 1:N_pos)
        log_lik_dens[k] = normal_lpdf(logD_pos[k] | mu_dens_obs[k], sigma_dens);
    }
  }
}
