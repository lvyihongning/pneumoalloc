// =====================================================================
// pneumoalloc: density among carriers (total, PCV10 or non-PCV10)
// inst/stan/density.stan
//
// Identical (apart from this header) to All_density_20261002.stan of the
// PCV10 Nha Trang analysis, except that its optional commune random
// effects were removed (pneumoalloc fits no cluster effects; the analysis
// used use_re = 0, so the fitted model is unchanged).
// =====================================================================
// Lognormal regression for total pneumococcal density among carriers.
// Aligned with the density component of
// 02_approaches_1_4/M1to4_unified_v4_20261002.stan.
//
// Model:
//   log(D_k) = X_pos[k] * beta
//            + shift(m_k ; gamma[survey or pooled])
//            + eps_k,          eps_k ~ Normal(0, sigma_dens)
//
//   shift(g ; gam) = sum_{j=1}^{g-1} gam[j]
//     gam[1] : log fold-change for m = 2 versus m = 1
//     gam[2] : additional log fold-change for m = 3 versus m = 2
//
// Changes relative to the previous version:
//   1. Multiplicity effects are UNCONSTRAINED. The previous [0, 3]
//      constraint plus half-normal prior was binding: lower CrI bounds
//      of 1.01-1.03 for three-versus-two serotypes were the constraint
//      boundary, not the data.
//   2. Optional survey-by-multiplicity interaction (mult_by_survey = 1)
//      replaces the separate interaction model, so Table S5 and the
//      primary density model come from one file.
//   3. (Optional commune random effects; removed in pneumoalloc.)
//   4. Prior means and SDs for beta are supplied as data vectors, so
//      priors no longer depend on the column ORDER of the design matrix
//      and sensitivity analyses need no recompilation.
//   5. Posterior predictive quantiles use COMMON RANDOM NUMBERS across
//      survey rounds and multiplicity levels. Previously numerator and
//      denominator of every median ratio were simulated independently,
//      which added Monte Carlo noise to the reported CrIs and pinned the
//      baseline to exactly 1 while all other rounds carried noise.
//   6. mu removed from transformed parameters (it was written to the CSV
//      on every draw).
//   7. PR_shiftD added: when the survey effect is a pure additive shift
//      on the log scale, every quantile ratio equals exp(delta) exactly.
//      shift_maxdev reports departure from that condition.
//   8. Optional pointwise log-likelihood.
//
// All density summaries are on the original GE/mL scale.
// Convert to log10 GE/mL in R for plotting only.

functions {
  real vec_median(vector x) {
    int n = num_elements(x);
    vector[n] xs = sort_asc(x);
    if (n == 1) return xs[1];
    if (n % 2 == 1) return xs[(n + 1) %/% 2];
    return 0.5 * (xs[n %/% 2] + xs[n %/% 2 + 1]);
  }

  // Stan-safe type-7 quantile
  real vec_quantile(vector x, real p) {
    int n = num_elements(x);
    vector[n] xs = sort_asc(x);
    real h;
    int lo;
    int hi;
    real w;

    if (n == 1) return xs[1];
    if (p <= 0) return xs[1];
    if (p >= 1) return xs[n];

    h = (n - 1) * p + 1;

    lo = 1;
    while (lo < n && (lo + 1) <= h) lo += 1;

    hi = lo;
    if (hi < n && hi < h) hi += 1;

    w = h - lo;
    if (hi == lo) return xs[lo];
    return xs[lo] + w * (xs[hi] - xs[lo]);
  }

  // Nested multiplicity shift on the natural-log density scale.
  // g = 1 returns 0. gam has length M - 1.
  real multiplicity_shift(int g, vector gam) {
    real out = 0;
    for (j in 1:(g - 1)) out += gam[j];
    return out;
  }
}

data {
  int<lower=1> T;                          // survey rounds
  int<lower=1> N_pos;                      // carriers
  int<lower=1> K;                          // design matrix columns
  int<lower=2, upper=5> M;                 // max detected serotypes (3 here)
  int<lower=1, upper=T> baseline;

  matrix[N_pos, K] X_pos;
  array[T] matrix[N_pos, K] X_pred_pos;    // counterfactual designs

  array[N_pos] int<lower=1, upper=T> survey_pos;   // observed round
  array[N_pos] int<lower=1, upper=5> m_sh;         // detected serotypes, 1..M

  vector[N_pos] logD;                      // log(dens_pos), natural log

  // Priors as data. beta_prior_mean / beta_prior_sd are built in R from
  // colnames(X_pos), so the priors do not depend on column ORDER:
  //   is_survey        <- colnames(X_pos) %in% survey_level_names
  //   beta_prior_mean  <- ifelse(is_survey, log(1e6), 0)
  //   beta_prior_sd    <- ifelse(is_survey, 3, 1)
  vector[K] beta_prior_mean;
  vector<lower=0>[K] beta_prior_sd;

  real<lower=0> prior_sd_gamma;            // [1]
  real<lower=0> prior_rate_sigma_dens;     // [0.5]  (was exponential(1):
                                           //  prior mean 1 on the log scale
                                           //  against an observed SD near 3)

  int<lower=0, upper=1> mult_by_survey;    // survey x multiplicity interaction
  int<lower=0, upper=1> compute_log_lik;

  int<lower=1000> M_rep;                   // predictive draws per round
}

transformed data {
  array[M] int n_m = rep_array(0, M);
  int n_gam = mult_by_survey ? T : 1;

  for (k in 1:N_pos) {
    if (m_sh[k] > M)
      reject("m_sh exceeds M at carrier k = ", k);
    n_m[m_sh[k]] += 1;
  }

  for (g in 1:M)
    print("multiplicity ", g, ": n = ", n_m[g]);

  if (mult_by_survey) {
    array[T, M] int n_tm = rep_array(0, T, M);
    for (k in 1:N_pos) n_tm[survey_pos[k], m_sh[k]] += 1;
    for (t in 1:T)
      for (g in 1:M)
        if (n_tm[t, g] == 0)
          print("WARNING: survey ", t, ", multiplicity ", g,
                " has no observations; the interaction term is ",
                "prior-driven in that cell.");
  }
}

parameters {
  vector[K] beta;
  real<lower=1e-6, upper=10> sigma_dens;

  // UNCONSTRAINED multiplicity effects (was real<lower=0, upper=3>)
  array[n_gam] vector[M - 1] gamma;
}

model {
  beta ~ normal(beta_prior_mean, beta_prior_sd);
  sigma_dens ~ exponential(prior_rate_sigma_dens);

  // Symmetric normal, not half-normal: the parameters are now unbounded.
  for (g in 1:n_gam)
    gamma[g] ~ normal(0, prior_sd_gamma);

  {
    vector[N_pos] mu = X_pos * beta;

    for (k in 1:N_pos)
      mu[k] += multiplicity_shift(m_sh[k],
                                  gamma[mult_by_survey ? survey_pos[k] : 1]);

    logD ~ normal(mu, sigma_dens);
  }
}

generated quantities {
  // ---- Overall summaries, observed multiplicity distribution held fixed
  //      across rounds (so the round comparison isolates the survey
  //      effect and does not run through changes in multiplicity).
  array[T] real<lower=0> meanD;
  array[T] real<lower=0> medianD;
  array[T] real<lower=0> q025D;
  array[T] real<lower=0> q25D;
  array[T] real<lower=0> q75D;
  array[T] real<lower=0> q975D;
  array[T] real<lower=0> IQRD;

  // ---- Counterfactual: every carrier held at m = g
  array[T, M] real<lower=0> meanD_by_m;
  array[T, M] real<lower=0> medianD_by_m;
  array[T, M] real<lower=0> q025D_by_m;
  array[T, M] real<lower=0> q25D_by_m;
  array[T, M] real<lower=0> q75D_by_m;
  array[T, M] real<lower=0> q975D_by_m;
  array[T, M] real<lower=0> IQRD_by_m;

  // ---- Ratios versus baseline
  array[T] real<lower=0> PR_meanD;         // exact, analytic
  array[T] real<lower=0> PR_medianD;       // common random numbers
  array[T] real<lower=0> PR_shiftD;        // exp(mean log shift)
  array[T] real shift_maxdev;              // 0 iff a pure additive shift
  array[T, M] real<lower=0> PR_meanD_by_m;
  array[T, M] real<lower=0> PR_medianD_by_m;

  // ---- Multiplicity fold-changes (Table S5)
  array[T, M] real<lower=0> dens_fold_vs1;        // g versus 1
  array[T, M] real<lower=0> dens_fold_stepwise;   // g versus g - 1
  array[T, M] real<lower=0> RR_dens_fold_vs1;     // ratio of ratios vs baseline

  array[compute_log_lik ? N_pos : 0] real log_lik_dens;

  {
    real eps = 1e-12;
    real half_v = 0.5 * square(sigma_dens);

    // ---- Common random numbers -------------------------------------
    // Drawn ONCE per posterior iteration and reused for every survey
    // round and multiplicity level, so a ratio of two predictive
    // quantiles carries no Monte Carlo noise of its own. Declared
    // inside this block, therefore not written to the CSV.
    array[M_rep] int rep_i;
    vector[M_rep] rep_z;
    for (j in 1:M_rep) {
      rep_i[j] = discrete_range_rng(1, N_pos);
      rep_z[j] = std_normal_rng();
    }

    vector[N_pos] mu_base_ref;   // baseline linear predictor, no multiplicity

    for (t in 1:T) {
      vector[N_pos] mu_lin = X_pred_pos[t] * beta;
      vector[M - 1] gam_t = gamma[mult_by_survey ? t : 1];

      if (t == baseline) mu_base_ref = mu_lin;

      // ---- analytic mean, observed multiplicity ----
      {
        vector[N_pos] mean_ind;
        for (k in 1:N_pos)
          mean_ind[k] = exp(mu_lin[k]
                            + multiplicity_shift(m_sh[k], gam_t)
                            + half_v);
        meanD[t] = mean(mean_ind);
      }

      // ---- predictive quantiles, observed multiplicity ----
      {
        vector[M_rep] D_rep;
        for (j in 1:M_rep) {
          int k = rep_i[j];
          D_rep[j] = exp(mu_lin[k]
                         + multiplicity_shift(m_sh[k], gam_t)
                         + sigma_dens * rep_z[j]);
        }
        medianD[t] = vec_median(D_rep);
        q025D[t]   = vec_quantile(D_rep, 0.025);
        q25D[t]    = vec_quantile(D_rep, 0.25);
        q75D[t]    = vec_quantile(D_rep, 0.75);
        q975D[t]   = vec_quantile(D_rep, 0.975);
        IQRD[t]    = fmax(q75D[t] - q25D[t], 0.0);
      }

      // ---- counterfactual by multiplicity ----
      for (g in 1:M) {
        real shift_g = multiplicity_shift(g, gam_t);
        vector[N_pos] mean_ind_g;
        vector[M_rep] D_rep_g;

        for (k in 1:N_pos)
          mean_ind_g[k] = exp(mu_lin[k] + shift_g + half_v);
        meanD_by_m[t, g] = mean(mean_ind_g);

        for (j in 1:M_rep)
          D_rep_g[j] = exp(mu_lin[rep_i[j]] + shift_g
                           + sigma_dens * rep_z[j]);

        medianD_by_m[t, g] = vec_median(D_rep_g);
        q025D_by_m[t, g]   = vec_quantile(D_rep_g, 0.025);
        q25D_by_m[t, g]    = vec_quantile(D_rep_g, 0.25);
        q75D_by_m[t, g]    = vec_quantile(D_rep_g, 0.75);
        q975D_by_m[t, g]   = vec_quantile(D_rep_g, 0.975);
        IQRD_by_m[t, g]    = fmax(q75D_by_m[t, g] - q25D_by_m[t, g], 0.0);

        dens_fold_vs1[t, g] = exp(shift_g);
        dens_fold_stepwise[t, g] =
          g == 1 ? 1.0 : exp(gam_t[g - 1]);
      }
    }

    // ---- shift diagnostic ----
    // If X_pred_pos[t] differs from X_pred_pos[baseline] only in the
    // survey columns, delta is constant across carriers, the whole
    // predictive distribution shifts by delta on the log scale, and
    // every quantile ratio equals exp(delta) exactly. shift_maxdev
    // near 0 confirms this; a non-zero value means the design contains
    // survey-by-covariate interactions and only PR_meanD stays exact.
    for (t in 1:T) {
      vector[N_pos] delta =
        (X_pred_pos[t] - X_pred_pos[baseline]) * beta;
      real dbar = mean(delta);
      PR_shiftD[t] = exp(dbar);
      shift_maxdev[t] = max(abs(delta - dbar));
    }

    // ---- ratios versus baseline ----
    {
      real b_mean   = fmax(meanD[baseline], eps);
      real b_median = fmax(medianD[baseline], eps);

      for (t in 1:T) {
        PR_meanD[t]   = meanD[t] / b_mean;
        PR_medianD[t] = medianD[t] / b_median;

        for (g in 1:M) {
          PR_meanD_by_m[t, g] =
            meanD_by_m[t, g] / fmax(meanD_by_m[baseline, g], eps);
          PR_medianD_by_m[t, g] =
            medianD_by_m[t, g] / fmax(medianD_by_m[baseline, g], eps);
          RR_dens_fold_vs1[t, g] =
            dens_fold_vs1[t, g] / fmax(dens_fold_vs1[baseline, g], eps);
        }
      }
    }

    // ---- pointwise log-likelihood ----
    if (compute_log_lik) {
      vector[N_pos] mu_obs = X_pos * beta;
      for (k in 1:N_pos)
        mu_obs[k] += multiplicity_shift(m_sh[k],
                                        gamma[mult_by_survey ? survey_pos[k] : 1]);
      for (k in 1:N_pos)
        log_lik_dens[k] = normal_lpdf(logD[k] | mu_obs[k], sigma_dens);
    }
  }
}

