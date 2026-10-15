// =====================================================================
// pneumoalloc: unified Bayesian framework, approaches i-iv, all children
// inst/stan/unified.stan
//
// Identical (apart from this header) to M1to4_unified_v4_20261002.stan,
// the model of the PCV10 Nha Trang analysis. pneumoalloc always fits it
// without commune random effects (use_re = 0, J = 1).
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
// Carrier weight in (b):  h_k = 1             (i-iii)
//                         h_k = D_k / Dbar    (iv)
// and the same weight is used in the standardisation in generated
// quantities (psi_i for i-iii, psi_i * E[D] for iv).
//
// Because (a) and (c) share no parameters with (b), P_overall, P_mult
// and meanD_pos have the same posterior under all four approaches (up
// to Monte Carlo error); only the composition differs.
//
// Changes relative to M4_psiweighted_v3.stan:
//   1. use_dens_weight flag (data): h_k = 1 when 0, D_k / Dbar when 1,
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

  // (Optionally density-weighted) fractional multinomial logit for
  // serotype composition among carriers.
  //   lp += h_k * sum_s q_ks * log rho_ks
  // h_k = w_dens[k]: 1 for approaches i-iii, D_k / Dbar for iv.
  real share_partial_log(
      array[] int idx_slice, int start, int end,
      matrix X_pos, matrix beta_share,
      array[] int m_sh,
      array[,] int sh_idx,
      array[,] real sh_w,
      vector w_dens,
      array[] int commune_id_pos,
      matrix u_commune_share
  ) {
    real lp = 0;
    int L    = size(idx_slice);
    int Kx   = cols(X_pos);
    int Jcat = cols(beta_share);
    matrix[L, Kx] X_slice;

    for (n in 1:L) X_slice[n] = X_pos[idx_slice[n]];
    matrix[L, Jcat] eta_mat = X_slice * beta_share;

    for (n in 1:L) {
      int k = idx_slice[n];
      row_vector[Jcat] eta =
        eta_mat[n] + u_commune_share[commune_id_pos[k]];
      real lse = log_sum_exp(eta);
      real acc = 0;

      for (j in 1:m_sh[k]) {
        int c = sh_idx[k, j];
        acc += sh_w[k, j] * (eta[c] - lse);
      }
      lp += w_dens[k] * acc;
    }
    return lp;
  }

  // Multiplicity model: categorical logit for m in 1..M among carriers.
  real mult_partial_log(
      array[] int idx_slice, int start, int end,
      matrix X_pos, matrix beta_mult,
      array[] int m_obs,
      array[] int commune_id_pos,
      matrix u_commune_mult
  ) {
    real lp = 0;
    int L    = size(idx_slice);
    int Kx   = cols(X_pos);
    int Mcat = cols(beta_mult);
    matrix[L, Kx] X_slice;

    for (n in 1:L) X_slice[n] = X_pos[idx_slice[n]];
    matrix[L, Mcat] eta_mat = X_slice * beta_mult;

    for (n in 1:L) {
      int k = idx_slice[n];
      row_vector[Mcat] eta =
        eta_mat[n] + u_commune_mult[commune_id_pos[k]];
      lp += eta[m_obs[k]] - log_sum_exp(eta);
    }
    return lp;
  }
}

data {
  int<lower=1> T;                       // survey rounds
  int<lower=2> C;                       // serotype categories
  int<lower=1> N;                       // all sampled children
  int<lower=1> N_pos;                   // carriers
  int<lower=1> K;                       // design matrix columns
  int<lower=1> J;                       // communes
  int<lower=2, upper=5> M;              // max detected serotypes (3 here)
  int<lower=1, upper=T> baseline;

  array[C] int<lower=0, upper=1> vt_mask_sero;

  matrix[N, K] X;
  array[T] matrix[N, K] X_pred;         // counterfactual designs

  array[N] int<lower=1, upper=J> commune_id;

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

  // 0 = approaches i-iii (h_k = 1); 1 = approach iv (h_k = D_k / Dbar)
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
  real<lower=0> prior_sd_sig_car;         // [1]
  real<lower=0> prior_sd_sig_share;       // [2]
  real<lower=0> prior_sd_sig_dens;        // [1]
  real<lower=0> prior_sd_sig_mult;        // [1]

  real dens_prior_center;                 // [log(1e6)]

  // Fit WITHOUT commune random effects (0) or WITH them (1).
  // J and commune_id must still be supplied when use_re = 0.
  int<lower=0, upper=1> use_re;

  int<lower=0, upper=1> compute_log_lik;  // 0 to save output size
}

transformed data {
  matrix[N_pos, K] X_pos;
  array[N_pos] int commune_id_pos;

  vector[N_pos] w_dens;
  array[N_pos] int idx_pos;
  vector[C] vt_mask_vec;

  // Dbar is computed here rather than passed in, so that
  // sum(w_dens) = N_pos holds by construction under every approach.
  real dens_scale = mean(dens_pos);

  real w_dens_sum;
  // Kish effective sample size of the carrier weights,
  // (sum w)^2 / sum(w^2). Equals N_pos when use_dens_weight = 0.
  real kish_ess_dens_weights;

  for (k in 1:N_pos) {
    X_pos[k] = X[pos_idx[k]];
    commune_id_pos[k] = commune_id[pos_idx[k]];
  }

  for (k in 1:N_pos)
    w_dens[k] = use_dens_weight ? dens_pos[k] / dens_scale : 1.0;
  for (k in 1:N_pos) idx_pos[k] = k;
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

  // Allocation weights must sum to one for every carrier, and a
  // carrier cannot have more allocation categories than detected
  // serotypes.
  for (k in 1:N_pos) {
    real sw = 0;
    for (j in 1:m_sh[k]) sw += sh_w[k, j];
    if (abs(sw - 1.0) > 1e-8)
      reject("Active sh_w values do not sum to 1 at carrier k = ", k,
             " (sum = ", sw, ")");
    if (m_sh[k] > m_obs[k])
      reject("m_sh > m_obs at carrier k = ", k);
  }

  print("use_dens_weight = ", use_dens_weight,
        " ; Dbar = ", dens_scale,
        " ; sum(w_dens) = ", w_dens_sum, " (should equal N_pos = ", N_pos, ")",
        " ; Kish effective sample size of carrier weights = ",
        kish_ess_dens_weights, " of ", N_pos, " carriers");
}

parameters {
  // Commune random effects are declared with zero length when
  // use_re = 0, so the same file fits both models.

  // ---- Carriage (Bernoulli) ----
  vector[K] beta_car;
  vector[use_re ? J : 0] z_commune_car;
  array[use_re] real<lower=0> sigma_commune_car;

  // ---- Density among carriers ----
  vector[K] beta_dens;
  real<lower=1e-6, upper=10> sigma_dens;
  real gamma_ge2;
  real gamma_ge3;
  vector[use_re ? J : 0] z_commune_dens;
  array[use_re] real<lower=0> sigma_commune_dens;

  // ---- Multiplicity among carriers ----
  // Category 1 (single serotype) is the reference.
  matrix[K, M-1] beta_mult_raw;
  matrix[use_re ? J : 0, M-1] z_commune_mult;
  array[use_re] real<lower=0> sigma_commune_mult;

  // ---- Serotype composition among carriers ----
  // Category C is the reference.
  matrix[K, C-1] beta_share_raw;
  matrix[use_re ? J : 0, C-1] z_commune_share;
  array[use_re] real<lower=0> sigma_commune_share;
}

transformed parameters {
  matrix[K, C] beta_share =
    append_col(beta_share_raw, rep_vector(0.0, K));

  matrix[K, M] beta_mult =
    append_col(rep_vector(0.0, K), beta_mult_raw);

  vector[J] u_commune_car  = rep_vector(0.0, J);
  vector[J] u_commune_dens = rep_vector(0.0, J);
  matrix[J, C] u_commune_share = rep_matrix(0.0, J, C);
  matrix[J, M] u_commune_mult  = rep_matrix(0.0, J, M);

  if (use_re) {
    u_commune_car  = sigma_commune_car[1]  * z_commune_car;
    u_commune_dens = sigma_commune_dens[1] * z_commune_dens;

    u_commune_share =
      append_col(sigma_commune_share[1] * z_commune_share,
                 rep_vector(0.0, J));

    u_commune_mult =
      append_col(rep_vector(0.0, J),
                 sigma_commune_mult[1] * z_commune_mult);
  }

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

  if (use_re) {
    z_commune_car  ~ std_normal();
    z_commune_dens ~ std_normal();
    to_vector(z_commune_share) ~ std_normal();
    to_vector(z_commune_mult)  ~ std_normal();

    sigma_commune_car[1]   ~ normal(0, prior_sd_sig_car);
    sigma_commune_dens[1]  ~ normal(0, prior_sd_sig_dens);
    sigma_commune_share[1] ~ normal(0, prior_sd_sig_share);
    sigma_commune_mult[1]  ~ normal(0, prior_sd_sig_mult);
  }

  // ---------------- Likelihood ----------------
  {
    vector[N] eta_car_obs = X * beta_car;
    vector[N_pos] mu_dens = X_pos * beta_dens;

    for (i in 1:N) eta_car_obs[i] += u_commune_car[commune_id[i]];

    for (k in 1:N_pos)
      mu_dens[k] += u_commune_dens[commune_id_pos[k]]
                  + multiplicity_shift(m_obs[k], gamma_ge2, gamma_ge3);

    // (a) Carriage among all children
    carrier ~ bernoulli_logit(eta_car_obs);

    // (c) Density among carriers (natural-log lognormal)
    logD_pos ~ normal(mu_dens, sigma_dens);
  }

  // (c) Multiplicity among carriers
  target += reduce_sum(
    mult_partial_log, idx_pos, 2000,
    X_pos, beta_mult, m_obs, commune_id_pos, u_commune_mult
  );

  // (b) Serotype composition among carriers. w_dens is 1 for
  //     approaches i-iii and the OBSERVED D_k / Dbar for approach iv.
  //     Do NOT replace it with predicted density.
  target += reduce_sum(
    share_partial_log, idx_pos, 2000,
    X_pos, beta_share, m_sh, sh_idx, sh_w, w_dens,
    commune_id_pos, u_commune_share
  );
}

generated quantities {
  // Diagnostics carried through from transformed data
  real dens_weight_kish_ess = kish_ess_dens_weights;
  real dens_weight_sum = w_dens_sum;

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

  // ---- Between-commune heterogeneity on the prevalence scale ----
  array[T, J] real<lower=0, upper=1> P_overall_commune;

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

      array[J] real psi_sum_j = rep_array(0.0, J);
      array[J] int  n_j       = rep_array(0, J);

      // ---------- population standardisation over all N ----------
      for (i in 1:N) {
        int ci = commune_id[i];

        real psi_i = inv_logit(eta_car_t[i] + u_commune_car[ci]);

        // multiplicity distribution given carriage
        row_vector[M] em = eta_mult_t[i] + u_commune_mult[ci];
        vector[M] pm = softmax(to_vector(em));

        // E[D | X_i, C = 1, t], marginalised over multiplicity
        real mu_base = mu_dens_base_t[i] + u_commune_dens[ci];
        real Dhat = 0;
        for (g in 1:M)
          Dhat += pm[g] * exp(mu_base
                              + multiplicity_shift(g, gamma_ge2, gamma_ge3)
                              + half_v);

        // within-stratum composition
        row_vector[C] es = eta_share_t[i] + u_commune_share[ci];
        vector[C] q = softmax(to_vector(es));

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

        psi_sum_j[ci] += psi_i;
        n_j[ci] += 1;
      }

      P_overall[t] = sum_psi / N;

      for (j in 1:J)
        P_overall_commune[t, j] = psi_sum_j[j] / n_j[j];

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
          int ci = commune_id[i];

          // carrier weight at observed multiplicity: 1 (i-iii) or
          // predicted density (iv)
          real w = 1.0;
          if (use_dens_weight)
            w = exp(mu_dens_base_t[i] + u_commune_dens[ci]
                    + multiplicity_shift(m_obs[k], gamma_ge2, gamma_ge3)
                    + half_v);

          row_vector[C] es = eta_share_t[i] + u_commune_share[ci];
          vector[C] q = softmax(to_vector(es));

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

      for (i in 1:N) eta_car_obs[i] += u_commune_car[commune_id[i]];

      for (k in 1:N_pos)
        mu_dens_obs[k] += u_commune_dens[commune_id_pos[k]]
                        + multiplicity_shift(m_obs[k], gamma_ge2, gamma_ge3);

      for (i in 1:N)
        log_lik_car[i] = bernoulli_logit_lpmf(carrier[i] | eta_car_obs[i]);

      for (k in 1:N_pos) {
        row_vector[M] eta_m =
          X_pos[k] * beta_mult + u_commune_mult[commune_id_pos[k]];
        log_lik_mult[k] =
          categorical_logit_lpmf(m_obs[k] | to_vector(eta_m));
      }

      for (k in 1:N_pos)
        log_lik_dens[k] = normal_lpdf(logD_pos[k] | mu_dens_obs[k], sigma_dens);
    }
  }
}
