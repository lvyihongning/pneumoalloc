# Design matrix ~ 0 + survey + covariates for the children in `rows`, and
# one counterfactual design per survey round (every child set to that
# round). The survey indicators are the first T columns: the Stan models
# give them their own prior in the density component.
design_matrices <- function(dat, rows, covariates = NULL,
                            call = rlang::caller_env()) {
  terms_rhs <- character()
  if (!is.null(covariates)) {
    if (!inherits(covariates, "formula") || length(covariates) != 2) {
      cli::cli_abort(
        "{.arg covariates} must be a one-sided formula such as {.code ~ age + sex}.",
        call = call
      )
    }
    unknown <- setdiff(all.vars(covariates), dat$covariates)
    if (length(unknown) > 0) {
      cli::cli_abort(c(
        "{.arg covariates} uses variable{?s} not declared in {.fn prepare_carriage_data}: {.val {unknown}}.",
        "i" = "Declaring them keeps the analysis population the same for every model."
      ), call = call)
    }
    terms_rhs <- attr(stats::terms(covariates), "term.labels")
  }
  full <- stats::reformulate(c(".survey", terms_rhs), intercept = FALSE)

  df <- dat$data[rows, , drop = FALSE]
  df$.survey <- dat$survey[rows]
  # covariate levels absent from this population (e.g. among multiple
  # carriers) would give empty design columns
  single <- character()
  for (v in all.vars(covariates %||% ~ 1)) {
    if (is.factor(df[[v]])) df[[v]] <- droplevels(df[[v]])
    discrete <- is.factor(df[[v]]) || is.character(df[[v]]) || is.logical(df[[v]])
    if (discrete && length(unique(df[[v]])) < 2) single <- c(single, v)
  }
  if (length(single) > 0) {
    cli::cli_abort(c(
      "Covariate{?s} {.val {single}} take{?s/} a single value among the children of this model.",
      "i" = "Drop {cli::qty(length(single))}{?it/them} from {.arg covariates} for this model."
    ), call = call)
  }
  empty <- dat$survey_levels[table(df$.survey) == 0]
  if (length(empty) > 0) {
    cli::cli_abort("No children in this population in survey round{?s} {.val {empty}}.",
                   call = call)
  }
  X <- stats::model.matrix(full, data = df)
  if (qr(X)$rank < ncol(X)) {
    cli::cli_abort(
      "The design matrix is rank deficient (collinear covariates or a covariate level with no children).",
      call = call
    )
  }
  X_pred <- lapply(dat$survey_levels, function(lev) {
    df$.survey <- factor(lev, levels = dat$survey_levels)
    stats::model.matrix(full, data = df)
  })
  stopifnot(all(startsWith(colnames(X)[seq_along(dat$survey_levels)], ".survey")))

  strip <- function(m) {
    m <- unname(as.matrix(m))
    storage.mode(m) <- "double"
    m
  }
  list(X = strip(X), X_pred = lapply(X_pred, strip), K = ncol(X),
       colnames = sub("^\\.survey", "survey", colnames(X)))
}

# Stan data for one model, plus what the fit object keeps alongside it.
# model: "unified", "unified_multiple", "states", "multiple_types",
# "presence" or "density" (with `density` options: group, interaction,
# multiplicity, m_rep).
prepare_model_data <- function(dat, model, approach = NULL, covariates = NULL,
                               priors = carriage_priors(), density = list(),
                               call = rlang::caller_env()) {
  check_class(dat, "carriage_data", call = call)
  check_class(priors, "carriage_priors", call = call)
  lev <- dat$survey_levels
  n_rounds <- length(lev)
  if (n_rounds < 2) {
    cli::cli_abort("Models need at least two survey rounds.", call = call)
  }
  ap <- NULL
  A <- NULL

  if (model == "unified") {
    ap <- resolve_approach(approach, call = call)
    rows <- seq_len(dat$N)
    des <- design_matrices(dat, rows, covariates, call = call)
    A <- allocation_matrix(dat$sero_long, dat$N_pos, dat$m_obs, dat$sero_levels,
                           ap$name, dat$has_abundance, call = call)
    check_reference_category(A, ap$name)
    # allocation uses all detected serotypes; multiplicity is capped at M
    sp <- build_share_sparse(A, max(dat$m_obs))
    stan_data <- list(
      T = n_rounds, C = length(dat$sero_levels), N = dat$N, N_pos = dat$N_pos,
      K = des$K, M = dat$M, baseline = dat$reference,
      vt_mask_sero = dat$vt_mask,
      X = des$X, X_pred = des$X_pred,
      carrier = dat$carrier,
      m_obs = pmin(dat$m_obs, dat$M),
      m_sh_max = sp$m_sh_max, m_sh = sp$m_sh, sh_idx = sp$sh_idx, sh_w = sp$sh_w,
      pos_idx = dat$pos_idx,
      logD_pos = log(dat$dens_pos), dens_pos = dat$dens_pos,
      use_dens_weight = ap$use_dens_weight,
      prior_sd_car = priors$car, prior_sd_share = priors$share,
      prior_sd_mult = priors$mult, prior_sd_dens_survey = priors$dens_survey,
      prior_sd_dens_cov = priors$dens_cov, prior_sd_gamma = priors$gamma,
      prior_rate_sigma_dens = priors$rate_sigma_dens,
      dens_prior_center = priors$dens_center,
      compute_log_lik = 0L
    )
    labels <- list(categories = dat$sero_levels, vt_mask = dat$vt_mask)
  } else if (model == "unified_multiple") {
    ap <- resolve_approach(approach, call = call)
    mp <- multiple_population(dat)
    rows <- mp$rows
    des <- design_matrices(dat, rows, covariates, call = call)
    A <- allocation_matrix(mp$sero_long, mp$N, mp$m_obs, mp$sero_levels,
                           ap$name, dat$has_abundance, call = call)
    check_reference_category(A, ap$name)
    sp <- build_share_sparse(A, max(mp$m_obs))   # m_model: 2 or 3+
    stan_data <- list(
      T = n_rounds, C = length(mp$sero_levels), N = mp$N, K = des$K,
      baseline = dat$reference,
      vt_mask_sero = mp$vt_mask,
      X = des$X, X_pred = des$X_pred,
      m_obs = mp$m_model,
      m_sh_max = sp$m_sh_max, m_sh = sp$m_sh, sh_idx = sp$sh_idx, sh_w = sp$sh_w,
      logD = log(mp$dens), dens = mp$dens,
      use_dens_weight = ap$use_dens_weight,
      prior_sd_share = priors$share, prior_sd_mult = priors$mult,
      prior_sd_dens_survey = priors$dens_survey,
      prior_sd_dens_cov = priors$dens_cov, prior_sd_gamma = priors$gamma,
      prior_rate_sigma_dens = priors$rate_sigma_dens,
      dens_prior_center = priors$dens_center
    )
    labels <- list(categories = mp$sero_levels, vt_mask = mp$vt_mask,
                   top = mp$top, top_plus = mp$top_plus)
  } else if (model == "states") {
    rows <- seq_len(dat$N)
    des <- design_matrices(dat, rows, covariates, call = call)
    states <- state_labels(dat$M, dat$capped)
    stan_data <- list(
      N = dat$N, K = des$K, T = n_rounds, G = length(states), ref = 1L,
      prior_sd = priors$states,
      y = pmin(dat$n_sero, dat$M) + 1L, X = des$X, X_pred = des$X_pred
    )
    labels <- list(categories = states)
  } else if (model == "multiple_types") {
    mp <- multiple_population(dat)
    rows <- mp$rows
    des <- design_matrices(dat, rows, covariates, call = call)
    types <- type_labels()
    stan_data <- list(
      N = mp$N, K = des$K, T = n_rounds, G = length(types), ref = length(types),
      prior_sd = priors$multiple_types,
      y = match(dat$multiple_type[mp$k], types), X = des$X, X_pred = des$X_pred
    )
    labels <- list(categories = types)
  } else if (model == "presence") {
    ap <- list(approach = NA_integer_, name = "presence", label = "Presence-based",
               use_dens_weight = 0L)
    rows <- seq_len(dat$N)
    des <- design_matrices(dat, rows, covariates, call = call)
    y <- presence_outcomes(dat)
    stan_data <- list(
      T = n_rounds, N = dat$N, K = des$K, C = ncol(y),
      y = t(y), X = des$X, X_pred = des$X_pred,
      baseline = dat$reference, prior_sd = priors$presence
    )
    labels <- list(categories = colnames(y), vt_mask = dat$vt_mask)
  } else if (model == "density") {
    grp <- density_group(dat, density$group, density$multiplicity, call = call)
    rows <- grp$rows
    # multiplicity levels observed in this group (at least 2 for the model)
    M_grp <- max(2L, grp$m)
    if (max(grp$m) < 2) {
      cli::cli_inform(
        "No carrier has more than one {density$group} serotype; the multiplicity effect is informed by its prior only."
      )
    }
    des <- design_matrices(dat, rows, covariates, call = call)
    is_survey <- seq_len(des$K) <= n_rounds
    stan_data <- list(
      T = n_rounds, N_pos = length(rows), K = des$K, M = M_grp,
      baseline = dat$reference,
      X_pos = des$X, X_pred_pos = des$X_pred,
      survey_pos = as.integer(dat$survey[rows]),
      m_sh = grp$m,
      logD = log(grp$density),
      beta_prior_mean = ifelse(is_survey, priors$dens_center, 0),
      beta_prior_sd = ifelse(is_survey, priors$dens_survey, priors$dens_cov),
      prior_sd_gamma = priors$gamma,
      prior_rate_sigma_dens = priors$rate_sigma_dens,
      mult_by_survey = as.integer(isTRUE(density$interaction)),
      compute_log_lik = 0L,
      M_rep = as.integer(density$m_rep %||% 4000L)
    )
    labels <- list(group = density$group, multiplicity = density$multiplicity,
                   interaction = isTRUE(density$interaction),
                   density = grp$density, m = grp$m,
                   m_levels = M_grp, m_plus = any(grp$m_raw > M_grp))
  } else {
    cli::cli_abort("Unknown model {.val {model}}.", call = call)
  }

  list(stan_data = stan_data, model = model, approach = ap, rows = rows,
       design = des$colnames, allocation = A, labels = labels)
}

# The composition priors are on log ratios against the last category, so
# that category should hold a fair part of the allocation.
check_reference_category <- function(A, approach) {
  mass <- sum(A[, ncol(A)])
  if (mass < 5) {
    cli::cli_warn(c(
      "The reference category {.val {colnames(A)[ncol(A)]}} has an allocation of only {round(mass, 1)} carriers under approach {.val {approach}}.",
      "i" = "Estimates then depend more on the prior; put a common serotype or a pooled category last in {.arg categories}."
    ))
  }
}

# Child-level presence indicators (children x outcomes): overall carriage,
# VT, NVT, multiple carriage, then each serotype category.
presence_outcomes <- function(dat) {
  sl <- dat$sero_long
  any_by_child <- function(flag) {
    out <- integer(dat$N)
    out[unique(sl$child[flag])] <- 1L
    out
  }
  W <- allocation_matrix(sl, dat$N_pos, dat$m_obs, dat$sero_levels,
                         "presence", dat$has_abundance)
  sero <- matrix(0L, dat$N, ncol(W), dimnames = list(NULL, colnames(W)))
  sero[dat$pos_idx, ] <- as.integer(W)
  cbind(
    Overall = dat$carrier,
    VT = any_by_child(sl$vaccine_type),
    NVT = any_by_child(!sl$vaccine_type),
    `Multiple carriage` = as.integer(dat$n_sero >= 2),
    sero
  )
}

# Carriers of a serotype group, their group density (sum of the
# serotype-specific densities, density x relative abundance) and the
# number of serotypes used for the multiplicity effect (capped at M, and
# uncapped in m_raw).
density_group <- function(dat, group, multiplicity = "group",
                          call = rlang::caller_env()) {
  if (group == "total") {
    return(list(rows = dat$pos_idx, density = dat$dens_pos,
                m = pmin(dat$m_obs, dat$M), m_raw = dat$m_obs))
  }
  if (!dat$has_abundance) {
    cli::cli_abort(
      "Group density needs relative abundance; supply {.arg abundance} to {.fn prepare_carriage_data}.",
      call = call
    )
  }
  sl <- dat$sero_long
  in_group <- if (group == "VT") sl$vaccine_type else !sl$vaccine_type
  k <- sort(unique(sl$k[in_group]))
  sub <- sl[in_group, ]
  dens_sero <- sub$abundance * dat$dens_pos[sub$k]
  m_group <- as.vector(table(factor(sub$k, levels = k)))
  m_raw <- if (multiplicity == "group") m_group else dat$m_obs[k]
  list(
    rows = dat$pos_idx[k],
    density = as.vector(rowsum(dens_sero, factor(sub$k, levels = k))),
    m = pmin(m_raw, dat$M),
    m_raw = m_raw
  )
}

#' Data list for a Stan model
#'
#' Builds the data list passed to one of the package's Stan models. The
#' fitting functions call it internally; use it directly to inspect what
#' is passed or to fit the model with your own Stan interface.
#'
#' @param dat A [prepare_carriage_data()] object.
#' @param model `"unified"` or `"unified_multiple"` ([fit_carriage()] with
#'   `population = "all"` or `"multiple"`), `"states"`
#'   ([fit_carriage_states()]), `"multiple_types"`
#'   ([fit_multiple_types()]), `"presence"` ([fit_presence()]) or
#'   `"density"` ([fit_density()]).
#' @param approach Approach number (1-4) or name, for the unified models.
#' @param group,interaction,multiplicity,m_rep For `model = "density"`;
#'   see [fit_density()]. `group` is one of `"total"`, `"VT"` or `"NVT"`.
#' @inheritParams fit_carriage
#'
#' @return A named list matching the data block of the Stan file given by
#'   [carriage_stan_file()].
#' @export
#' @examples
#' sim <- simulate_carriage_data(seed = 1)
#' dat <- prepare_carriage_data(sim, "survey", c("sero1", "sero2", "sero3"),
#'   density = "density", abundance = c("ra1", "ra2", "ra3"),
#'   covariates = c("age_z", "sex"))
#' sd <- carriage_stan_data(dat, "unified", "density", covariates = ~ age_z + sex)
#' str(sd[c("T", "C", "N", "N_pos", "K", "M", "use_dens_weight")])
carriage_stan_data <- function(dat, model = c("unified", "unified_multiple",
                                              "states", "multiple_types",
                                              "presence", "density"),
                               approach = NULL, covariates = NULL,
                               priors = carriage_priors(), group = "total",
                               interaction = FALSE, multiplicity = "group",
                               m_rep = 4000L) {
  model <- rlang::arg_match(model)
  density <- list(group = rlang::arg_match0(group, c("total", "VT", "NVT")),
                  interaction = interaction,
                  multiplicity = rlang::arg_match0(multiplicity, c("group", "total")),
                  m_rep = m_rep)
  prepare_model_data(dat, model, approach, covariates, priors, density)$stan_data
}
