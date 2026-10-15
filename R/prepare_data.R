#' Prepare carriage survey data
#'
#' Builds the analysis data shared by every model: the analysis population
#' (complete cases for the covariates), carriers, the number of detected
#' serotypes, density, the dominant serotype, the type of multiple carriage
#' and the serotype categories.
#'
#' @param data A data frame with one row per sampled child.
#' @param survey Name of the column holding the survey round. Its factor
#'   levels (or sorted unique values) give the order of the rounds.
#' @param serotypes Names of the columns holding the detected serotypes,
#'   one column per slot (e.g. `c("sero1", "sero2", "sero3")`). Empty
#'   strings and `NA` mean "no serotype in this slot"; a child with no
#'   serotype in any slot is a non-carrier.
#' @param density Name of the column holding total pneumococcal density
#'   (e.g. genome equivalents/mL). Must be positive and finite for every
#'   carrier; ignored for non-carriers.
#' @param vaccine_types Serotypes of the vaccine-type (VT) group, e.g.
#'   [pcv10_types()] or [vaccine_types()]`("PCV13", add = "6A")`; any
#'   character vector can be used. All other detected serotypes form the
#'   non-vaccine-type (NVT) group. Names must match the spelling in `data`.
#' @param abundance Names of the columns holding the relative abundance of
#'   the serotype in the matching `serotypes` slot. Needed for the
#'   `abundance` and `density` approaches. Normalised to sum to one within
#'   each child.
#' @param categories Serotypes modelled as their own category in the
#'   unified framework; other detected serotypes are pooled into
#'   `other_label`. The default keeps every detected serotype. Every
#'   detected vaccine serotype must be in `categories`, so that the pooled
#'   category holds only non-vaccine serotypes.
#' @param other_label Label of the pooled category.
#' @param covariates Names of the covariate columns that models may adjust
#'   for. The analysis population is restricted to children with complete
#'   data on these and on `survey`, so that unadjusted and adjusted models
#'   use the same children. Prepare covariates (centring, scaling,
#'   reference levels) before calling this function.
#' @param reference Survey round used as reference for prevalence ratios,
#'   normally the pre-vaccination survey. Defaults to the first round.
#' @param dominant How the dominant serotype is chosen: `"abundance"`
#'   (highest relative abundance; ties broken by slot order) or `"first"`
#'   (the first filled slot). Defaults to `"abundance"` when `abundance` is
#'   given, otherwise `"first"`.
#' @param max_serotypes Largest number of detected serotypes distinguished
#'   in the models of carriage states, multiplicity and density (2 to 5).
#'   `NULL` (default) uses the largest number detected in the data.
#'   Children with more detected serotypes are counted in the top category,
#'   labelled e.g. `"3+ serotypes"`. The allocation approaches always use
#'   all detected serotypes (a child with four serotypes contributes 1/4 to
#'   each under `equal`). The number of `serotypes` columns can be anything.
#'
#' @return An object of class `carriage_data`.
#'
#' Serotype categories are ordered as vaccine serotypes (in the order of
#' `vaccine_types`), then non-vaccine serotypes (in the order of
#' `categories`), then the pooled category. The **last** category is the
#' reference category of the composition model. Categories that are never
#' detected are dropped.
#'
#' @export
#' @examples
#' sim <- simulate_carriage_data(seed = 1)
#' dat <- prepare_carriage_data(
#'   sim,
#'   survey = "survey",
#'   serotypes = c("sero1", "sero2", "sero3"),
#'   abundance = c("ra1", "ra2", "ra3"),
#'   density = "density",
#'   vaccine_types = pcv10_types(),
#'   categories = c(pcv10_types(), "6A", "15B/C", "19A", "11A"),
#'   covariates = c("age_z", "sex")
#' )
#' dat
prepare_carriage_data <- function(data, survey, serotypes, density,
                                  vaccine_types = pcv10_types(),
                                  abundance = NULL, categories = NULL,
                                  other_label = "Other NVT",
                                  covariates = NULL, reference = NULL,
                                  dominant = NULL, max_serotypes = NULL) {
  if (!is.data.frame(data)) cli::cli_abort("{.arg data} must be a data frame.")
  data <- as.data.frame(data)

  needed <- c(survey, serotypes, abundance, density, covariates)
  missing_cols <- setdiff(needed, names(data))
  if (length(missing_cols) > 0) {
    cli::cli_abort("Column{?s} not found in {.arg data}: {.val {missing_cols}}.")
  }
  if (length(survey) != 1 || length(density) != 1) {
    cli::cli_abort("{.arg survey} and {.arg density} must each name one column.")
  }
  if (!is.null(abundance) && length(abundance) != length(serotypes)) {
    cli::cli_abort("{.arg abundance} must name one column per {.arg serotypes} column.")
  }
  dominant <- dominant %||% if (is.null(abundance)) "first" else "abundance"
  dominant <- rlang::arg_match0(dominant, c("abundance", "first"))
  if (dominant == "abundance" && is.null(abundance)) {
    cli::cli_abort("{.code dominant = \"abundance\"} needs {.arg abundance}.")
  }

  # ---- analysis population: complete cases ----
  cc <- stats::complete.cases(data[, c(survey, covariates), drop = FALSE])
  n_dropped <- sum(!cc)
  if (n_dropped > 0) {
    cli::cli_inform(
      "Excluded {n_dropped} of {nrow(data)} children with missing survey or covariate data."
    )
  }
  d <- droplevels(data[cc, , drop = FALSE])
  rownames(d) <- NULL
  N <- nrow(d)

  survey_f <- d[[survey]]
  if (!is.factor(survey_f)) survey_f <- factor(survey_f)
  survey_f <- droplevels(survey_f)
  survey_levels <- levels(survey_f)
  if (length(survey_levels) < 1) cli::cli_abort("No survey rounds in {.arg data}.")
  reference <- reference %||% survey_levels[1]
  if (length(reference) != 1 || !reference %in% survey_levels) {
    cli::cli_abort("{.arg reference} {.val {reference}} is not a survey round.")
  }

  # ---- detected serotypes ----
  S <- as.matrix(d[, serotypes, drop = FALSE])
  S <- matrix(trimws(as.character(S)), nrow = N)
  S[!is.na(S) & S == ""] <- NA
  n_sero <- as.integer(rowSums(!is.na(S)))
  carrier <- as.integer(n_sero > 0)
  pos_idx <- which(carrier == 1L)
  N_pos <- length(pos_idx)
  if (N_pos == 0) cli::cli_abort("No carriers in the analysis population.")
  m_obs <- n_sero[pos_idx]
  M_observed <- max(m_obs)
  if (M_observed < 2) {
    cli::cli_abort(c(
      "No child carries more than one serotype.",
      "i" = "All allocation approaches then give the same result."
    ))
  }
  if (is.null(max_serotypes)) {
    if (M_observed > 5) {
      cli::cli_abort(c(
        "Up to {M_observed} serotypes are detected per child; the models distinguish at most 5.",
        "i" = "Set {.arg max_serotypes} (2 to 5); children with more serotypes are counted in the top category."
      ))
    }
    M <- M_observed
  } else {
    ok <- is.numeric(max_serotypes) && length(max_serotypes) == 1 &&
      max_serotypes %in% 2:5
    if (!ok) cli::cli_abort("{.arg max_serotypes} must be a whole number from 2 to 5.")
    M <- as.integer(max_serotypes)
    if (M > M_observed) {
      cli::cli_inform(
        "No child carries more than {M_observed} serotypes; using {.code max_serotypes = {M_observed}}."
      )
      M <- M_observed
    }
  }

  dens <- suppressWarnings(as.numeric(d[[density]]))
  dens_pos <- dens[pos_idx]
  bad_dens <- !is.finite(dens_pos) | dens_pos <= 0
  if (any(bad_dens)) {
    cli::cli_abort(c(
      "{.arg density} must be positive and finite for every carrier.",
      "x" = "{sum(bad_dens)} carrier{?s} with missing, zero or negative density."
    ))
  }

  # ---- long table: one row per carrier x detected serotype ----
  hit <- which(!is.na(S[pos_idx, , drop = FALSE]), arr.ind = TRUE)
  hit <- hit[order(hit[, 1], hit[, 2]), , drop = FALSE]
  sero_long <- data.frame(
    k = hit[, 1],
    child = pos_idx[hit[, 1]],
    slot = hit[, 2],
    serotype = S[cbind(pos_idx[hit[, 1]], hit[, 2])],
    stringsAsFactors = FALSE
  )
  sero_long$vaccine_type <- sero_long$serotype %in% vaccine_types
  if (!any(sero_long$vaccine_type)) {
    cli::cli_warn(c(
      "None of {.arg vaccine_types} is detected in the data.",
      "i" = "Check that serotype names are spelled as in the data (e.g. {.val 15B/C} versus {.val 15B})."
    ))
  }

  if (!is.null(abundance)) {
    R <- as.matrix(d[, abundance, drop = FALSE])
    ra <- suppressWarnings(as.numeric(R[cbind(sero_long$child, sero_long$slot)]))
    if (any(!is.finite(ra) | ra <= 0)) {
      cli::cli_abort(
        "{.arg abundance} must be positive and finite for every detected serotype."
      )
    }
    sero_long$abundance <- ra / stats::ave(ra, sero_long$k, FUN = sum)
  }

  first_max <- function(v) seq_along(v) == which.max(v)
  sero_long$dominant <- if (dominant == "abundance") {
    as.logical(stats::ave(sero_long$abundance, sero_long$k, FUN = first_max))
  } else {
    !duplicated(sero_long$k)
  }

  # ---- type of multiple carriage ----
  n_vt <- as.vector(rowsum(as.integer(sero_long$vaccine_type), sero_long$k))
  dominant_vt <- sero_long$vaccine_type[sero_long$dominant][order(sero_long$k[sero_long$dominant])]
  type <- ifelse(n_vt == m_obs, "VT-only",
                 ifelse(n_vt == 0, "NVT-only",
                        ifelse(dominant_vt, "Mixed, VT-dominant", "Mixed, NVT-dominant")))
  type[m_obs < 2] <- NA_character_

  # ---- serotype categories ----
  detected <- unique(sero_long$serotype)
  categories <- categories %||% detected
  vt_pooled <- setdiff(intersect(detected, vaccine_types), categories)
  if (length(vt_pooled) > 0) {
    cli::cli_abort(c(
      "All detected vaccine serotypes must be in {.arg categories}.",
      "x" = "Not in {.arg categories}: {.val {vt_pooled}}."
    ))
  }
  if (other_label %in% detected) {
    cli::cli_abort("{.arg other_label} {.val {other_label}} is also a detected serotype.")
  }
  sero_long$category <- ifelse(sero_long$serotype %in% categories,
                               sero_long$serotype, other_label)
  present <- unique(sero_long$category)
  absent <- setdiff(categories, present)
  if (length(absent) > 0) {
    cli::cli_inform("Never detected, dropped from {.arg categories}: {.val {absent}}.")
  }
  vt_cats <- intersect(vaccine_types, present)
  nvt_cats <- setdiff(intersect(categories, present), vaccine_types)
  sero_levels <- c(vt_cats, nvt_cats, if (other_label %in% present) other_label)
  if (length(sero_levels) < 2) cli::cli_abort("At least two serotype categories are needed.")

  structure(
    list(
      data = d,
      survey = survey_f,
      survey_levels = survey_levels,
      reference = match(reference, survey_levels),
      reference_level = reference,
      N = N, N_pos = N_pos, M = M, M_observed = M_observed,
      capped = M_observed > M,
      n_sero = n_sero,
      carrier = carrier, pos_idx = pos_idx, m_obs = m_obs,
      dens_pos = dens_pos,
      multiple_type = type,
      sero_long = sero_long,
      sero_levels = sero_levels,
      vt_mask = as.integer(sero_levels %in% vaccine_types),
      vaccine_types = vaccine_types,
      other_label = other_label,
      covariates = covariates,
      has_abundance = !is.null(abundance),
      dominant = dominant,
      n_dropped = n_dropped,
      columns = list(survey = survey, serotypes = serotypes,
                     abundance = abundance, density = density)
    ),
    class = "carriage_data"
  )
}

#' @export
print.carriage_data <- function(x, ...) {
  cli::cli_h1("pneumoalloc carriage data")
  cli::cli_text(
    "{x$N} children ({x$n_dropped} excluded), {x$N_pos} carriers, ",
    "{sum(x$m_obs >= 2)} with multiple serotype carriage; ",
    "{length(x$survey_levels)} survey rounds (reference {.val {x$reference_level}})."
  )
  cli::cli_text("Children by survey round and number of detected serotypes:")
  print(table(survey = x$survey,
              serotypes = factor(x$n_sero, levels = 0:x$M_observed)))
  if (x$capped) {
    cli::cli_text("Models count up to {x$M} serotypes; children with more are in the {.val {x$M}+} category.")
  }
  ref <- x$sero_levels[length(x$sero_levels)]
  cli::cli_text(
    "{length(x$sero_levels)} serotype categories ",
    "({sum(x$vt_mask)} vaccine type{?s}); reference category {.val {ref}}."
  )
  if (length(x$covariates) > 0) {
    cli::cli_text("Covariates: {.val {x$covariates}}.")
  } else {
    cli::cli_text("Covariates: none (unadjusted models only).")
  }
  approaches <- allocation_approaches()$name
  if (!x$has_abundance) approaches <- approaches[1:2]
  cli::cli_text("Available approaches: {.val {approaches}}.")
  invisible(x)
}

# The children with multiple carriage, in the layout used by the
# multiple-carriage models: their rows in `dat$data`, their carrier index,
# detected serotypes and the serotype categories present among them.
multiple_population <- function(dat) {
  k <- which(dat$m_obs >= 2)
  sl <- dat$sero_long[dat$sero_long$k %in% k, , drop = FALSE]
  sl$k <- match(sl$k, k)
  sero_levels <- dat$sero_levels[dat$sero_levels %in% unique(sl$category)]
  top <- min(dat$M, 3L)   # the multiple-carriage model distinguishes 2 and 3+
  list(
    rows = dat$pos_idx[k],
    k = k,
    N = length(k),
    m_obs = dat$m_obs[k],
    m_model = pmin(dat$m_obs[k], top),
    top = top,
    top_plus = any(dat$m_obs[k] > top),
    dens = dat$dens_pos[k],
    sero_long = sl,
    sero_levels = sero_levels,
    vt_mask = as.integer(sero_levels %in% dat$vaccine_types)
  )
}
