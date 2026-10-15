#' Check a fit
#'
#' Collects the checks of the PCV10 Nha Trang analysis for one fit:
#'
#' * **sampler**: divergent transitions, iterations at the maximum
#'   treedepth and E-BFMI, by chain;
#' * **convergence**: the largest R-hat and smallest bulk and tail ESS over
#'   the model parameters and over the reported estimates;
#' * **identities**: relations that hold in every posterior draw (shares
#'   and category probabilities sum to one; serotype, or VT plus NVT,
#'   prevalences sum to overall prevalence; for unadjusted unified fits,
#'   population- and carrier-standardised shares agree). A failure means an
#'   implementation error;
#' * **crude**: posterior medians against the crude estimates of
#'   [describe_prevalence()] or [describe_carriage_states()]. Unadjusted
#'   models reproduce them up to prior shrinkage. Shrinkage is visible for
#'   rare outcomes and in small populations, e.g. composition among few
#'   multiple carriers, where categories without any allocation in a survey
#'   are pulled towards the others. Adjusted models may differ.
#'
#' For the `density` approach, the summary also reports the Kish effective
#' sample size of the carrier weights \eqn{\Theta_i / \bar{\Theta}_t} as a
#' percentage of carriers. Because density is skewed, few high-density
#' carriers dominate the weighted composition likelihood; when this
#' percentage is low, the posterior intervals of serotype-specific
#' estimates understate the uncertainty.
#'
#' @param object A `pneumoalloc_fit` object, or a `pneumoalloc_fits`
#'   object (checked fit by fit).
#' @param rhat_max Largest acceptable R-hat.
#' @param ess_min_per_chain Smallest acceptable bulk and tail ESS per chain.
#' @param tol Tolerance of the identity checks.
#' @return A `pneumoalloc_check` object: a list with `summary` (one row per
#'   check, with `pass`), `sampler`, `convergence`, `identities` and
#'   `crude`. For `pneumoalloc_fits`, a list of these.
#' @export
check_fit <- function(object, rhat_max = 1.01, ess_min_per_chain = 100,
                      tol = 1e-6) {
  if (inherits(object, "pneumoalloc_fits")) {
    return(lapply(object, check_fit, rhat_max = rhat_max,
                  ess_min_per_chain = ess_min_per_chain, tol = tol))
  }
  check_class(object, "pneumoalloc_fit", arg = "object")
  nchains <- object$fit$num_chains()
  ess_min <- ess_min_per_chain * nchains

  diag <- object$fit$diagnostic_summary(quiet = TRUE)
  sampler <- tibble::tibble(
    chain = seq_along(diag$num_divergent),
    num_divergent = diag$num_divergent,
    num_max_treedepth = diag$num_max_treedepth,
    ebfmi = diag$ebfmi
  )

  params <- switch(
    object$type,
    unified = c("beta_car", "beta_dens", "sigma_dens", "gamma_ge2", "gamma_ge3",
                "beta_mult_raw", "beta_share_raw"),
    unified_multiple = c("beta_share_raw", "beta_dens", "sigma_dens",
                         "gamma_ge3", "beta_m3"),
    presence = "beta",
    density = c("beta", "sigma_dens", "gamma"),
    "beta_raw"
  )
  ps <- quiet_ess(posterior::summarise_draws(
    object$fit$draws(params), posterior::default_convergence_measures()
  ))
  est <- if (object$type == "density") marginal_density(object) else marginal_prevalence(object)
  worst <- function(rhat, bulk, tail) {
    c(max(rhat, na.rm = TRUE), min(bulk, na.rm = TRUE), min(tail, na.rm = TRUE))
  }
  convergence <- tibble::tibble(
    quantities = c("parameters", "estimates"),
    max_rhat = c(worst(ps$rhat, ps$ess_bulk, ps$ess_tail)[1],
                 worst(est$rhat, est$ess_bulk, est$ess_tail)[1]),
    min_ess_bulk = c(min(ps$ess_bulk, na.rm = TRUE), min(est$ess_bulk, na.rm = TRUE)),
    min_ess_tail = c(min(ps$ess_tail, na.rm = TRUE), min(est$ess_tail, na.rm = TRUE))
  )

  identities <- identity_checks(object, tol)
  crude <- crude_comparison(object, est)
  if (object$type == "density") {
    crude_label <- "crude: max |relative difference| of median density (%)"
    max_crude_diff <- max(abs(crude$diff_pct))
  } else {
    crude_label <- "crude: max |difference| (pp, outcomes > 0.5%)"
    big <- crude$value_crude > 0.005
    max_crude_diff <- if (any(big)) max(abs(crude$diff_pp[big])) else NA_real_
  }
  applies <- identities$applies
  identity_value <- if (any(applies)) max(identities$max_abs_deviation[applies]) else NA_real_
  identity_pass <- if (any(applies)) all(identities$pass[applies]) else NA

  # approach iv: how much information the density weights leave
  ess_rows <- NULL
  if (object$type %in% c("unified", "unified_multiple") &&
      identical(object$approach$name, "density")) {
    n_w <- if (object$type == "unified") object$data$N_pos else length(object$rows)
    ess <- stats::median(draws_array(object, "dens_weight_kish_ess"))
    ess_rows <- tibble::tibble(
      check = "density weights: effective sample size (% of carriers)",
      value = 100 * ess / n_w,
      threshold = "information; intervals assume all carriers (see ?check_fit)",
      pass = NA
    )
  }
  summary <- tibble::tibble(
    check = c("divergent transitions", "max treedepth hits", "min E-BFMI",
              "max R-hat", "min bulk ESS", "min tail ESS",
              "identities: max deviation",
              crude_label),
    value = c(sum(sampler$num_divergent), sum(sampler$num_max_treedepth),
              min(sampler$ebfmi), max(convergence$max_rhat),
              min(convergence$min_ess_bulk), min(convergence$min_ess_tail),
              identity_value,
              max_crude_diff),
    threshold = c("0", "0", ">= 0.3", paste("<", rhat_max),
                  paste(">=", ess_min), paste(">=", ess_min),
                  if (any(applies)) paste("<=", tol) else "none for this model",
                  if (is.null(object$covariates)) "small (unadjusted)" else "not expected (adjusted)"),
    pass = c(sum(sampler$num_divergent) == 0, sum(sampler$num_max_treedepth) == 0,
             min(sampler$ebfmi) >= 0.3, max(convergence$max_rhat) < rhat_max,
             min(convergence$min_ess_bulk) >= ess_min,
             min(convergence$min_ess_tail) >= ess_min,
             identity_pass, NA)
  )

  summary <- rbind(summary, ess_rows)
  structure(
    list(summary = summary, sampler = sampler, convergence = convergence,
         identities = identities, crude = crude,
         label = paste0(model_label(object),
                        if (!is.null(object$approach)) paste0(", ", object$approach$label),
                        if (is.null(object$covariates)) ", unadjusted" else ", adjusted")),
    class = "pneumoalloc_check"
  )
}

#' @export
print.pneumoalloc_check <- function(x, ...) {
  cli::cli_h2("Checks: {x$label}")
  s <- x$summary
  value <- formatC(s$value, digits = 5, format = "g")
  status <- ifelse(is.na(s$pass), "i", ifelse(s$pass, "v", "x"))
  for (i in seq_len(nrow(s))) {
    cli::cli_bullets(stats::setNames(
      paste0(s$check[i], ": ", value[i], " (", s$threshold[i], ")"), status[i]
    ))
  }
  invisible(x)
}

identity_checks <- function(object, tol) {
  sum_dev <- function(arr) max(abs(apply(arr, c(1, 2), sum) - 1))
  get <- function(v) draws_array(object, v)
  unadj <- is.null(object$covariates)
  dev <- switch(
    object$type,
    unified = {
      po <- get("P_overall")[, , 1]
      c(share_sums_to_1 = sum_dev(get("share")),
        P_sero_sums_to_P_overall = max(abs(apply(get("P_sero"), c(1, 2), sum) - po)),
        P_VT_plus_P_NVT_eq_P_overall = max(abs(get("P_VT")[, , 1] + get("P_NVT")[, , 1] - po)),
        P_mult_sums_to_1 = sum_dev(get("P_mult")),
        share_eq_share_oc = max(abs(get("share") - get("share_oc"))))
    },
    unified_multiple = c(share_sums_to_1 = sum_dev(get("share")),
                         share_VT_plus_NVT_eq_1 = max(abs(get("share_VT") + get("share_NVT") - 1))),
    states = ,
    multiple_types = c(P_class_sums_to_1 = sum_dev(get("P_class"))),
    numeric()
  )
  applies <- rep(TRUE, length(dev))
  if (object$type == "unified") applies[length(dev)] <- unadj
  tibble::tibble(check = names(dev) %||% character(), max_abs_deviation = unname(dev),
                 applies = applies, pass = unname(dev) <= tol)
}

crude_comparison <- function(object, est) {
  dat <- object$data
  if (object$type == "density") {
    obs <- describe_density(dat, object$labels$group)
    est <- est[est$outcome == "Density" & est$statistic == "median", ]
    m <- match(est$survey, obs$survey)
    return(tibble::tibble(
      statistic = "median", survey = est$survey, median = est$median,
      q2.5 = est$q2.5, q97.5 = est$q97.5, value_crude = obs$median[m],
      diff_pct = 100 * (est$median / obs$median[m] - 1)
    ))
  }
  if (object$type == "presence") {
    crude <- describe_prevalence(dat, "presence")
    st <- describe_carriage_states(dat)
    st <- st[st$outcome == "Multiple carriage", ]
    crude <- rbind(crude[, c("outcome", "survey", "value")],
                   st[, c("outcome", "survey", "value")])
  } else if (object$type %in% c("unified", "unified_multiple")) {
    pop <- if (object$type == "unified") "all" else "multiple"
    crude <- describe_prevalence(dat, object$approach$name, pop)
    est <- est[est$measure == crude$measure[1], ]
  } else {
    crude <- describe_carriage_states(dat)
    want <- if (object$type == "states") "prevalence" else "proportion among multiple carriers"
    crude <- crude[crude$measure == want, ]
  }
  key_e <- paste(est$outcome, est$survey)
  key_c <- paste(crude$outcome, crude$survey)
  m <- match(key_e, key_c)
  tibble::tibble(
    measure = est$measure, group = est$group, outcome = est$outcome,
    survey = est$survey, median = est$median,
    q2.5 = est$q2.5, q97.5 = est$q97.5,
    value_crude = crude$value[m],
    diff_pp = 100 * (est$median - crude$value[m])
  )[!is.na(m), ]
}

#' Compare allocation approaches
#'
#' In the unified framework, carriage, multiplicity and density share no
#' parameters with the composition component, so overall carriage, the
#' distribution of the number of serotypes and mean density have the same
#' posterior under every approach. Separate fits give different draws, so
#' agreement is judged against Monte Carlo error:
#' `z = (mean_a - mean_ref) / sqrt(mcse_a^2 + mcse_ref^2)`.
#' Also binds the estimates and ratios of all approaches for tables and
#' figures.
#'
#' A presence-based fit ([fit_presence()]) with the same covariates can be
#' added to the list of fits with `population = "all"`. It is then included
#' in `prevalence` and `ratios`, but not in the agreement check, which
#' applies to the allocation approaches only.
#'
#' @param fits A `pneumoalloc_fits` object from [fit_carriage()], or a list
#'   of such fits with the same data, population and covariates, optionally
#'   with a [fit_presence()] fit, e.g.
#'   `c(fit_carriage(dat), list(presence = fit_presence(dat)))`.
#' @param reference_approach Approach (number or name) used as reference.
#' @param z_flag Shared quantities with `|z|` above this are flagged.
#' @inheritParams prevalence_ratios
#' @return A list of tibbles: `shared` (the agreement check), `prevalence`
#'   ([marginal_prevalence()]) and `ratios` ([prevalence_ratios()]) of all
#'   approaches.
#' @export
compare_approaches <- function(fits, reference_approach = "density", z_flag = 4,
                               reference = "baseline") {
  ok <- is.list(fits) && length(fits) >= 2 &&
    all(vapply(fits, inherits, TRUE, "pneumoalloc_fit"))
  if (!ok) cli::cli_abort("{.arg fits} must hold at least two {.cls pneumoalloc_fit} objects.")
  all_fits <- fits
  is_presence <- vapply(fits, function(f) f$type == "presence", TRUE)
  fits <- fits[!is_presence]
  if (length(fits) < 2) cli::cli_abort("{.arg fits} must hold at least two fits of {.fn fit_carriage}.")
  f1 <- fits[[1]]
  same <- vapply(fits, function(f) {
    f$type == f1$type && f$type %in% c("unified", "unified_multiple") &&
      identical(f$rows, f1$rows) && identical(f$design, f1$design) &&
      identical(f$labels$categories, f1$labels$categories)
  }, TRUE)
  if (!all(same)) {
    cli::cli_abort("All fits must come from {.fn fit_carriage} with the same data, population and covariates.")
  }
  presence_ok <- vapply(all_fits[is_presence], function(f) {
    f1$type == "unified" && identical(f$rows, f1$rows) && identical(f$design, f1$design)
  }, TRUE)
  if (!all(presence_ok)) {
    cli::cli_abort("A presence-based fit must have the same data and covariates as the {.code population = \"all\"} fits.")
  }
  ref <- resolve_approach(reference_approach)$approach
  if (!ref %in% vapply(fits, function(f) f$approach$approach, 1L)) {
    cli::cli_abort("The reference approach {ref} is not among {.arg fits}.")
  }

  shared_vars <- if (f1$type == "unified") {
    list(P_overall = "Overall", P_mult = state_labels(f1$data$M, f1$data$capped)[-1],
         meanD_pos = "Mean density")
  } else if (f1$labels$top == 3) {
    list(P_m3 = count_label(3, f1$labels$top_plus), meanD = "Mean density")
  } else {
    list(meanD = "Mean density")
  }
  shared <- do.call(rbind, lapply(unname(fits), function(f) {
    do.call(rbind, Map(function(v, out) {
      e <- entry(draws_array(f, v), out, v, v)
      t1 <- entry_table(f, e, f$data$survey_levels)
      t1[, c("approach", "method", "group", "outcome", "survey", "mean", "mcse_mean")]
    }, names(shared_vars), shared_vars))
  }))
  names(shared)[names(shared) == "group"] <- "quantity"
  key <- paste(shared$quantity, shared$outcome, shared$survey)
  is_ref <- shared$approach == ref
  m <- match(key, key[is_ref])
  shared$mean_ref <- shared$mean[is_ref][m]
  shared$mcse_ref <- shared$mcse_mean[is_ref][m]
  shared$diff <- shared$mean - shared$mean_ref
  shared$z <- shared$diff / sqrt(shared$mcse_mean^2 + shared$mcse_ref^2)
  shared$flag <- abs(shared$z) > z_flag
  shared <- tibble::as_tibble(shared[!is_ref, ])
  if (any(shared$flag, na.rm = TRUE)) {
    cli::cli_warn(c(
      "Shared quantities differ between approaches by more than Monte Carlo error.",
      "i" = "Check that all fits used the same data, design and priors."
    ))
  }

  groups <- c(if (f1$type == "unified") "overall", "VT", "NVT", "serotype")
  all_fits <- structure(unname(all_fits), class = "pneumoalloc_fits")
  list(
    shared = shared,
    prevalence = marginal_prevalence(all_fits, groups = groups),
    ratios = prevalence_ratios(all_fits, groups = groups, reference = reference)
  )
}
