#' Prior scales
#'
#' Prior scales of all models. They are passed to Stan as data, so
#' changing them needs no recompilation. The defaults are those of the PCV10
#' Nha Trang analysis (Table S2 of its supplement). All normal priors have
#' mean 0 unless stated.
#'
#' @param car SD of the normal prior on the carriage coefficients (survey
#'   intercepts and covariates, log odds) in [fit_carriage()].
#' @param share SD of the normal prior on the composition coefficients
#'   (log share ratios).
#' @param mult SD of the normal prior on the multiplicity coefficients
#'   (log odds).
#' @param dens_survey SD of the normal prior on the survey-specific
#'   log-density levels, which is centred at `dens_center`. Used by
#'   [fit_carriage()] and [fit_density()], as are the next four.
#' @param dens_cov SD of the normal prior on the density covariate
#'   coefficients (log density ratios).
#' @param gamma SD of the normal prior on the additional log density for
#'   carriage of at least two and at least three serotypes.
#' @param rate_sigma_dens Rate of the exponential prior on the residual SD
#'   of log density.
#' @param dens_center Prior centre of the survey-specific log-density
#'   levels (natural log, e.g. `log(1e6)` genome equivalents/mL).
#' @param states SD of the normal prior on all coefficients of
#'   [fit_carriage_states()].
#' @param multiple_types SD of the normal prior on all coefficients of
#'   [fit_multiple_types()].
#' @param presence SD of the normal prior on all coefficients of
#'   [fit_presence()].
#'
#' @return A named list of class `carriage_priors`.
#' @export
#' @examples
#' carriage_priors()
#'
#' # sensitivity analysis with a wider prior on serotype composition
#' carriage_priors(share = 5)
carriage_priors <- function(car = 2.5, share = 2, mult = 2, dens_survey = 3,
                            dens_cov = 1, gamma = 1, rate_sigma_dens = 0.5,
                            dens_center = log(1e6), states = 2.5,
                            multiple_types = 2, presence = 2.5) {
  scales <- list(car = car, share = share, mult = mult,
                 dens_survey = dens_survey, dens_cov = dens_cov,
                 gamma = gamma, rate_sigma_dens = rate_sigma_dens,
                 states = states, multiple_types = multiple_types,
                 presence = presence)
  ok <- vapply(scales, function(v) is.numeric(v) && length(v) == 1 &&
                 is.finite(v) && v > 0, TRUE)
  if (!all(ok)) {
    cli::cli_abort("Prior scale{?s} must be single positive numbers: {.arg {names(scales)[!ok]}}.")
  }
  if (!is.numeric(dens_center) || length(dens_center) != 1 || !is.finite(dens_center)) {
    cli::cli_abort("{.arg dens_center} must be a single finite number.")
  }
  structure(c(scales, dens_center = dens_center), class = "carriage_priors")
}

#' @export
print.carriage_priors <- function(x, ...) {
  tab <- data.frame(
    model = c(rep("Unified framework and density models", 8), "Carriage states",
              "Types of multiple carriage", "Presence-based"),
    parameters = c(
      "Carriage: survey intercepts and covariates (log odds)",
      "Composition: survey intercepts and covariates (log share ratio)",
      "Multiplicity: survey intercepts and covariates (log odds)",
      "Density: survey-specific log density",
      "Density: covariates (log density ratio)",
      "Density: >= 2 and >= 3 serotypes (gamma)",
      "Density: residual SD",
      "(centre of the survey log-density prior)",
      "Survey intercepts and covariates (log odds)",
      "Survey intercepts and covariates (log odds)",
      "Survey intercepts and covariates, each outcome (log odds)"
    ),
    prior = c(
      sprintf("normal(0, %s)", x$car), sprintf("normal(0, %s)", x$share),
      sprintf("normal(0, %s)", x$mult),
      sprintf("normal(%s, %s)", format(signif(x$dens_center, 4)), x$dens_survey),
      sprintf("normal(0, %s)", x$dens_cov), sprintf("normal(0, %s)", x$gamma),
      sprintf("exponential(%s)", x$rate_sigma_dens), "",
      sprintf("normal(0, %s)", x$states), sprintf("normal(0, %s)", x$multiple_types),
      sprintf("normal(0, %s)", x$presence)
    )
  )
  tab <- tab[tab$prior != "", ]
  for (m in unique(tab$model)) {
    cli::cli_text("{.strong {m}}")
    rows <- tab[tab$model == m, ]
    for (i in seq_len(nrow(rows))) {
      cli::cli_text("  {rows$parameters[i]}: {.code {rows$prior[i]}}")
    }
  }
  invisible(x)
}
