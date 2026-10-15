#' Simulate carriage survey data
#'
#' Simulates repeated cross-sectional carriage surveys with multiple
#' serotype carriage, relative abundance and density, in the layout
#' expected by [prepare_carriage_data()]. Vaccine-type carriage declines
#' over the rounds while overall carriage stays roughly constant. Used in
#' the examples, tests and vignette.
#'
#' @param n_per_survey Children per survey round.
#' @param surveys Labels of the survey rounds, in order; the first is the
#'   pre-vaccination round.
#' @param vt_share Share of PCV10 serotypes among carried serotypes in each
#'   round (recycled to the number of rounds).
#' @param p_carriage Carriage prevalence at `age_z = 0`.
#' @param p_multiple Probability that a carrier carries 1, 2, 3, ...
#'   serotypes; its length (2 to 5) is the largest number of serotypes per
#'   child and the number of `sero` and `ra` columns.
#' @param seed Optional random seed; the global random number state is
#'   restored afterwards.
#'
#' @return A data frame with one row per child: `id`, `survey`, `age_z`,
#'   `sex`, the serotype columns `sero1`, `sero2`, ... (dominant serotype
#'   first), the relative abundance columns `ra1`, `ra2`, ..., `density`
#'   (`NA` for non-carriers) and four further covariates for adjusted
#'   models: `income_z` (standardised log income), `household` (household
#'   size: `"<=3"`, `"4-5"`, `">=6"`), `contact_u5` (regular contact with
#'   children under 5, 0/1) and `chronic_ill` (0/1). These four do not
#'   affect carriage in the simulation.
#' @export
#' @examples
#' sim <- simulate_carriage_data(seed = 1)
#' head(sim)
simulate_carriage_data <- function(n_per_survey = 400,
                                   surveys = c("2016", "2017", "2019"),
                                   vt_share = c(0.5, 0.3, 0.1),
                                   p_carriage = 0.4,
                                   p_multiple = c(0.75, 0.2, 0.05),
                                   seed = NULL) {
  if (!is.null(seed)) {
    return(withr::with_seed(seed, simulate_carriage_data(
      n_per_survey, surveys, vt_share, p_carriage, p_multiple
    )))
  }
  vt <- c("6B", "14", "19F", "23F")
  nvt <- c("6A", "15B/C", "19A", "11A", "35B", "23A", "34", "10A", "16F", "21")
  nvt_w <- c(0.2, 0.18, 0.15, 0.12, 0.1, 0.08, 0.07, 0.04, 0.03, 0.03)
  n_rounds <- length(surveys)
  vt_share <- rep_len(vt_share, n_rounds)
  N <- n_per_survey * n_rounds

  survey_i <- rep(seq_len(n_rounds), each = n_per_survey)
  age_z <- stats::rnorm(N)
  sex <- factor(sample(c("female", "male"), N, replace = TRUE))

  p_car <- stats::plogis(stats::qlogis(p_carriage) - 0.3 * age_z)
  carrier <- stats::runif(N) < p_car
  n_max <- length(p_multiple)
  if (n_max < 2 || n_max > 5) cli::cli_abort("{.arg p_multiple} must have 2 to 5 elements.")
  m <- ifelse(carrier, sample(seq_len(n_max), N, replace = TRUE, prob = p_multiple), 0L)

  sero <- matrix(NA_character_, N, n_max)
  ra <- matrix(NA_real_, N, n_max)
  for (i in which(carrier)) {
    w <- c(rep(vt_share[survey_i[i]] / length(vt), length(vt)),
           (1 - vt_share[survey_i[i]]) * nvt_w)
    s <- sample(c(vt, nvt), m[i], prob = w)
    a <- stats::rgamma(m[i], shape = c(4, 1.5, 0.8, 0.5, 0.4)[seq_len(m[i])])
    a <- sort(a / sum(a), decreasing = TRUE)
    sero[i, seq_len(m[i])] <- s
    ra[i, seq_len(m[i])] <- a
  }
  log_dens <- log(1e6) + 0.6 * (m - 1) - 0.2 * age_z + stats::rnorm(N, 0, 1.2)

  out <- data.frame(
    id = seq_len(N),
    survey = factor(surveys[survey_i], levels = surveys),
    age_z = age_z,
    sex = sex,
    stringsAsFactors = FALSE
  )
  colnames(sero) <- paste0("sero", seq_len(n_max))
  colnames(ra) <- paste0("ra", seq_len(n_max))
  out <- cbind(out, as.data.frame(sero, stringsAsFactors = FALSE), as.data.frame(ra))
  out$density <- ifelse(carrier, exp(log_dens), NA_real_)

  # further covariates (drawn last, so the columns above do not depend on them)
  out$income_z <- stats::rnorm(N)
  out$household <- factor(sample(c("<=3", "4-5", ">=6"), N, replace = TRUE,
                                 prob = c(0.3, 0.5, 0.2)),
                          levels = c("<=3", "4-5", ">=6"))
  out$contact_u5 <- stats::rbinom(N, 1, 0.4)
  out$chronic_ill <- stats::rbinom(N, 1, 0.05)
  out
}
