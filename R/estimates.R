# ---- draws ---------------------------------------------------------------

# Draws of a generated quantity as an array [draw, survey, outcome].
# Draws are ordered chain by chain.
draws_array <- function(object, variable) {
  rv <- posterior::as_draws_rvars(object$fit$draws(variable))[[variable]]
  arr <- posterior::draws_of(rv)
  d <- dim(arr)
  if (length(d) == 2) arr <- array(arr, c(d, 1)) else arr <- array(arr, d)
  arr
}

entry <- function(arr, outcome, group, measure, statistic = NULL) {
  list(arr = arr, outcome = outcome, outcome1 = outcome[1],
       group = rep_len(group, length(outcome)), measure = measure,
       statistic = statistic)
}

# The population-standardised estimates of a fit, one entry per quantity.
estimate_entries <- function(object) {
  lab <- object$labels$categories
  get <- function(v) draws_array(object, v)
  switch(
    object$type,
    unified = {
      carriers <- "proportion among carriers"
      M <- object$data$M
      list(
        entry(get("P_overall"), "Overall", "overall", "prevalence"),
        entry(get("P_VT"), "VT", "VT", "prevalence"),
        entry(get("P_NVT"), "NVT", "NVT", "prevalence"),
        entry(get("P_sero"), lab, "serotype", "prevalence"),
        entry(get("share_VT"), "VT", "VT", carriers),
        entry(get("share_NVT"), "NVT", "NVT", carriers),
        entry(get("share"), lab, "serotype", carriers),
        entry(get("P_mult"), state_labels(M, object$data$capped)[-1],
              "multiplicity", carriers)
      )
    },
    unified_multiple = {
      multi <- "proportion among multiple carriers"
      out <- list(
        entry(get("share_VT"), "VT", "VT", multi),
        entry(get("share_NVT"), "NVT", "NVT", multi),
        entry(get("share"), lab, "serotype", multi)
      )
      # with at most two serotypes there is no multiplicity to estimate
      if (object$labels$top == 3) {
        p3 <- get("P_m3")
        out <- c(out, list(entry(
          array(c(1 - p3, p3), c(dim(p3)[1:2], 2)),
          c("2 serotypes", count_label(3, object$labels$top_plus)),
          "multiplicity", multi
        )))
      }
      out
    },
    states = {
      p <- get("P_class")
      G <- dim(p)[3]
      list(
        entry(p, lab, "state", "prevalence"),
        entry(array(apply(p[, , -1, drop = FALSE], c(1, 2), sum), c(dim(p)[1:2], 1)),
              "Carriage", "carriage", "prevalence"),
        entry(array(apply(p[, , 3:G, drop = FALSE], c(1, 2), sum), c(dim(p)[1:2], 1)),
              "Multiple carriage", "multiple", "prevalence")
      )
    },
    multiple_types = list(
      entry(get("P_class"), lab, "type", "proportion among multiple carriers")
    ),
    presence = {
      groups <- c("overall", "VT", "NVT", "multiple", rep("serotype", length(lab) - 4))
      list(entry(get("P"), lab, groups, "prevalence"))
    },
    cli::cli_abort("Use {.fn marginal_density} and {.fn density_ratios} for density fits.",
                   call = rlang::caller_env(2))
  )
}

select_entries <- function(entries, groups) {
  if (is.null(groups)) return(entries)
  available <- unique(unlist(lapply(entries, `[[`, "group")))
  unknown <- setdiff(groups, available)
  if (length(unknown) > 0) {
    cli::cli_abort(c(
      "Unknown group{?s} for this model: {.val {unknown}}.",
      "i" = "Available: {.val {available}}."
    ), call = rlang::caller_env())
  }
  keep <- lapply(entries, function(e) {
    j <- which(e$group %in% groups)
    if (length(j) == 0) return(NULL)
    e$arr <- e$arr[, , j, drop = FALSE]
    e$outcome <- e$outcome[j]
    e$group <- e$group[j]
    e
  })
  Filter(Negate(is.null), keep)
}

# ---- summaries -------------------------------------------------------------

# Posterior summary of every [survey, outcome] cell of an array
summarise_array <- function(arr, nchains) {
  D <- dim(arr)[1]
  cells <- matrix(arr, nrow = D)
  q <- apply(cells, 2, stats::quantile, probs = c(0.025, 0.5, 0.975), names = FALSE)
  conv <- quiet_ess(t(apply(cells, 2, function(v) {
    m <- matrix(v, ncol = nchains)
    c(posterior::rhat(m), posterior::ess_bulk(m), posterior::ess_tail(m),
      posterior::mcse_mean(m))
  })))
  tibble::tibble(
    mean = colMeans(cells), median = q[2, ], sd = apply(cells, 2, stats::sd),
    q2.5 = q[1, ], q97.5 = q[3, ],
    rhat = conv[, 1], ess_bulk = conv[, 2], ess_tail = conv[, 3],
    mcse_mean = conv[, 4]
  )
}

# Long table for an entry whose array has survey rounds `surveys`
entry_table <- function(object, e, surveys) {
  nt <- dim(e$arr)[2]
  nj <- dim(e$arr)[3]
  vt <- object$labels$vt_mask
  info <- tibble::tibble(
    model = object$type,
    approach = object$approach$approach %||% NA_integer_,
    method = object$approach$label %||% fit_label(object),
    adjusted = !is.null(object$covariates),
    measure = e$measure,
    group = rep(e$group, each = nt),
    outcome = rep(e$outcome, each = nt),
    survey = factor(rep(surveys, times = nj), levels = object$data$survey_levels)
  )
  if (!is.null(e$statistic)) info$statistic <- e$statistic
  if (!is.null(vt)) {
    info$vaccine_type <- ifelse(
      info$group == "serotype",
      ifelse(info$outcome %in% object$data$vaccine_types, "VT", "NVT"),
      NA_character_
    )
  }
  tibble::as_tibble(cbind(info, summarise_array(e$arr, object$fit$num_chains())))
}

bind_fits <- function(object, fun, ...) {
  if (inherits(object, "pneumoalloc_fits")) {
    return(tibble::as_tibble(do.call(rbind, lapply(unname(object), fun, ...))))
  }
  check_class(object, "pneumoalloc_fit", arg = "object", call = rlang::caller_env())
  NULL
}

#' Marginal prevalence and proportions
#'
#' Population-standardised estimates by survey round: every child is set
#' to each survey round, predictions are averaged over the children within
#' each posterior draw, and the draws are summarised.
#'
#' The groups available depend on the model:
#'
#' * [fit_carriage()], all children: `"overall"`, `"VT"`, `"NVT"` and
#'   `"serotype"` prevalence, the same as proportions among carriers, and
#'   `"multiplicity"` (number of serotypes among carriers);
#' * [fit_carriage()], `population = "multiple"`: `"VT"`, `"NVT"`,
#'   `"serotype"` and `"multiplicity"` proportions among multiple carriers;
#' * [fit_carriage_states()]: `"state"` (no carriage, 1, 2, ...
#'   serotypes), `"carriage"` and `"multiple"`;
#' * [fit_multiple_types()]: `"type"`.
#'
#' @param object A `pneumoalloc_fit` or `pneumoalloc_fits` object.
#' @param groups Groups of outcomes to return (see Details); all by default.
#' @return A tibble with one row per outcome and survey round: posterior
#'   mean, median, SD, 2.5% and 97.5% quantiles (proportions, 0-1),
#'   convergence measures, and `estimate`, the text
#'   "median (2.5%, 97.5%)" in percent.
#' @export
marginal_prevalence <- function(object, groups = NULL) {
  out <- bind_fits(object, marginal_prevalence, groups = groups)
  if (!is.null(out)) return(out)
  entries <- select_entries(estimate_entries(object), groups)
  tab <- do.call(rbind, lapply(entries, entry_table, object = object,
                               surveys = object$data$survey_levels))
  tab$estimate <- fmt_ci(100 * tab$median, 100 * tab$q2.5, 100 * tab$q97.5, 2)
  tibble::as_tibble(tab)
}

#' Prevalence ratios
#'
#' Ratios of the population-standardised estimates of
#' [marginal_prevalence()] between survey rounds, formed within each
#' posterior draw. By default every round is compared with the reference
#' (pre-vaccination) round set in [prepare_carriage_data()]; with
#' `reference = "previous"` each round is compared with the round before
#' it. For proportions (among carriers or among multiple carriers) the
#' ratios are proportion ratios.
#'
#' @inheritParams marginal_prevalence
#' @param reference `"baseline"` (the reference round of the data),
#'   `"previous"` (the preceding round) or the label of a survey round.
#' @param pr_threshold Ratios are flagged as not `reportable` when the
#'   upper 97.5% limit of the reference-round estimate is below this value
#'   (default 0.1%), because they are then unstable.
#' @return A tibble with one row per outcome and compared round, with the
#'   reference round in `reference`, posterior summaries of the ratio and
#'   `estimate`, the text "median (2.5%, 97.5%)".
#' @export
prevalence_ratios <- function(object, groups = NULL, reference = "baseline",
                              pr_threshold = 0.001) {
  out <- bind_fits(object, prevalence_ratios, groups = groups,
                   reference = reference, pr_threshold = pr_threshold)
  if (!is.null(out)) return(out)
  lev <- object$data$survey_levels
  ref <- reference_index(reference, lev, object$data$reference)
  keep <- which(!is.na(ref))
  entries <- select_entries(estimate_entries(object), groups)
  tabs <- lapply(entries, function(e) {
    num <- e$arr[, keep, , drop = FALSE]
    den <- e$arr[, ref[keep], , drop = FALSE]
    ratio <- e
    ratio$arr <- num / den
    tab <- entry_table(object, ratio, lev[keep])
    tab$reference <- rep(lev[ref[keep]], times = length(e$outcome))
    den_q <- apply(den, c(2, 3), stats::quantile, probs = 0.975, names = FALSE)
    tab$reference_q97.5 <- as.vector(den_q)
    tab
  })
  tab <- do.call(rbind, tabs)
  tab$reportable <- tab$reference_q97.5 >= pr_threshold
  tab$ratio <- ifelse(tab$measure == "prevalence", "prevalence ratio", "proportion ratio")
  tab$estimate <- fmt_ci(tab$median, tab$q2.5, tab$q97.5, 2)
  first <- c("model", "approach", "method", "adjusted", "measure", "ratio",
             "group", "outcome", "survey", "reference")
  tibble::as_tibble(tab[, c(first, setdiff(names(tab), first))])
}

#' Density among carriers
#'
#' Population-standardised density among carriers by survey round.
#'
#' * [fit_density()]: the mean (analytic) and the median (posterior
#'   predictive) density, with the observed number of serotypes held fixed
#'   across rounds, and the same if every carrier carried 1, 2, ...
#'   serotypes (`outcome`);
#' * [fit_carriage()]: the mean density, marginalised over the predicted
#'   number of serotypes (and, for `population = "all"`, the mean if every
#'   carrier carried 1, 2, ... serotypes).
#'
#' The column `statistic` says which summary of the density distribution is
#' estimated.
#'
#' @param object A `pneumoalloc_fit` or `pneumoalloc_fits` object from
#'   [fit_density()] or [fit_carriage()].
#' @param statistic `"mean"` and/or `"median"`; all available by default.
#'   The median is estimated by [fit_density()] only.
#' @return A tibble of posterior summaries (density on its original scale,
#'   e.g. genome equivalents/mL).
#' @export
marginal_density <- function(object, statistic = NULL) {
  out <- bind_fits(object, marginal_density, statistic = statistic)
  if (!is.null(out)) return(out)
  entries <- select_statistic(density_entries(object), statistic)
  tab <- do.call(rbind, lapply(entries, entry_table, object = object,
                               surveys = object$data$survey_levels))
  tab$estimate <- fmt_ci(tab$median, tab$q2.5, tab$q97.5, 0)
  tibble::as_tibble(tab)
}

#' Density ratios
#'
#' Ratios of the density among carriers ([marginal_density()]) between
#' survey rounds, formed within each posterior draw, and optionally the
#' density ratios by number of detected serotypes.
#'
#' @inheritParams prevalence_ratios
#' @inheritParams marginal_density
#' @param multiplicity Also return the density ratios by number of detected
#'   serotypes (e.g. two versus one), with `outcome` naming the comparison.
#'   Without a survey x multiplicity interaction they do not depend on survey
#'   round and are reported once (`survey` is `NA`). With
#'   [fit_density()]`(interaction = TRUE)` they are reported by survey round,
#'   together with their ratio to the reference round (`measure` "ratio of
#'   density ratios").
#' @return A tibble of posterior summaries of the ratios.
#' @export
density_ratios <- function(object, reference = "baseline", statistic = NULL,
                           multiplicity = FALSE) {
  out <- bind_fits(object, density_ratios, reference = reference,
                   statistic = statistic, multiplicity = multiplicity)
  if (!is.null(out)) return(out)
  lev <- object$data$survey_levels
  ref <- reference_index(reference, lev, object$data$reference)
  keep <- which(!is.na(ref))
  entries <- density_entries(object)
  entries <- select_statistic(entries[vapply(entries, `[[`, "", "outcome1") == "Density"],
                              statistic)
  tab <- do.call(rbind, lapply(entries, function(e) {
    e$arr <- e$arr[, keep, , drop = FALSE] / e$arr[, ref[keep], , drop = FALSE]
    e$measure <- "density ratio"
    t1 <- entry_table(object, e, lev[keep])
    t1$reference <- lev[ref[keep]]
    t1
  }))
  if (multiplicity) tab <- rbind(tab, multiplicity_ratios(object, lev, ref))
  tab$estimate <- fmt_ci(tab$median, tab$q2.5, tab$q97.5, 2)
  tibble::as_tibble(tab)
}

# Number of multiplicity levels a fit distinguishes, and whether the top
# level also holds children with more serotypes ("3+")
mult_levels <- function(object) {
  if (object$type == "density" && !is.null(object$labels$m_levels)) {
    return(list(M = object$labels$m_levels, plus = object$labels$m_plus))
  }
  list(M = object$data$M, plus = object$data$capped)
}

multiplicity_ratios <- function(object, lev, ref) {
  ml <- mult_levels(object)
  M <- ml$M
  lab <- function(g) sub(" serotypes?$", "", count_label(g, ml$plus && g == M))
  vs <- function(g, h) paste(lab(g), "vs", lab(h), if (h == 1) "serotype" else "serotypes")
  pooled <- function(arr, outcome) {
    e <- entry(arr, outcome, density_group_of(object), "density ratio",
               statistic = NA_character_)
    t1 <- entry_table(object, e, lev[1])
    t1$survey <- factor(NA, levels = lev)
    t1$reference <- NA_character_
    t1
  }
  if (object$type == "unified") {
    # the density component has shifts for >= 2 and >= 3 serotypes
    top <- if (M > 3) "3+" else lab(3)
    vars <- stats::setNames("dens_fold_2vs1", paste(lab(2), "vs 1 serotype"))
    if (M >= 3) {
      vars <- c(vars, stats::setNames(c("dens_fold_3vs2", "dens_fold_3vs1"),
                                      c(paste(top, "vs 2 serotypes"),
                                        paste(top, "vs 1 serotype"))))
    }
    return(do.call(rbind, Map(function(v, nm) pooled(draws_array(object, v), nm),
                              vars, names(vars))))
  }
  if (object$type == "unified_multiple") {
    if (object$labels$top < 3) return(NULL)
    top <- sub(" serotypes$", "", count_label(3, object$labels$top_plus))
    return(pooled(draws_array(object, "dens_fold_3vs2"), paste(top, "vs 2 serotypes")))
  }
  # fit_density(): fold changes by survey, [draw, survey, serotypes]
  vs1 <- draws_array(object, "dens_fold_vs1")[, , -1, drop = FALSE]
  outcome <- vapply(2:M, function(g) vs(g, 1), "")
  arr <- vs1
  if (M >= 3) {
    step <- draws_array(object, "dens_fold_stepwise")[, , -(1:2), drop = FALSE]
    arr <- array(c(vs1, step), c(dim(vs1)[1:2], dim(vs1)[3] + dim(step)[3]))
    outcome <- c(outcome, vapply(3:M, function(g) vs(g, g - 1), ""))
  }
  if (!object$labels$interaction) {
    return(pooled(arr[, 1, , drop = FALSE], outcome))
  }
  e <- entry(arr, outcome, density_group_of(object), "density ratio",
             statistic = NA_character_)
  by_survey <- entry_table(object, e, lev)
  by_survey$reference <- NA_character_
  keep <- which(!is.na(ref))
  e$arr <- arr[, keep, , drop = FALSE] / arr[, ref[keep], , drop = FALSE]
  e$measure <- "ratio of density ratios"
  rr <- entry_table(object, e, lev[keep])
  rr$reference <- rep(lev[ref[keep]], times = length(outcome))
  rbind(by_survey, rr)
}

density_group_of <- function(object) {
  if (object$type == "density") object$labels$group else "total"
}

select_statistic <- function(entries, statistic, call = rlang::caller_env()) {
  if (is.null(statistic)) return(entries)
  available <- unique(vapply(entries, `[[`, "", "statistic"))
  bad <- setdiff(statistic, available)
  if (length(bad) > 0) {
    cli::cli_abort(c(
      "Statistic {.val {bad}} is not estimated by this model.",
      "i" = "Available: {.val {available}}. Medians come from {.fn fit_density}."
    ), call = call)
  }
  entries[vapply(entries, `[[`, "", "statistic") %in% statistic]
}

density_entries <- function(object) {
  get <- function(v) draws_array(object, v)
  ml <- mult_levels(object)
  by_m <- paste("All carriers with", state_labels(ml$M, ml$plus)[-1])
  switch(
    object$type,
    unified = {
      # the density component distinguishes 1, 2 and 3+ serotypes
      M3 <- min(ml$M, 3L)
      by_m3 <- paste("All carriers with", state_labels(M3, ml$plus || ml$M > 3)[-1])
      list(
        entry(get("meanD_pos"), "Density", "total", "density among carriers", "mean"),
        entry(get("meanD_pos_by_m")[, , seq_len(M3), drop = FALSE], by_m3, "total",
              "density among carriers", "mean")
      )
    },
    unified_multiple = list(
      entry(get("meanD"), "Density", "total", "density among multiple carriers", "mean")
    ),
    density = {
      g <- object$labels$group
      measure <- c(total = "density among carriers",
                   VT = "VT density among VT carriers",
                   NVT = "NVT density among NVT carriers")[[g]]
      list(
        entry(get("meanD"), "Density", g, measure, "mean"),
        entry(get("medianD"), "Density", g, measure, "median"),
        entry(get("meanD_by_m"), by_m, g, measure, "mean"),
        entry(get("medianD_by_m"), by_m, g, measure, "median")
      )
    },
    cli::cli_abort("Density is estimated by {.fn fit_density} and {.fn fit_carriage}.",
                   call = rlang::caller_env(2))
  )
}

#' Summarise a fit
#'
#' Collects the results of a fit in one list: [marginal_prevalence()] and
#' [prevalence_ratios()] for prevalence models, and [marginal_density()]
#' and [density_ratios()] for models that estimate density.
#'
#' @param object A `pneumoalloc_fit` or `pneumoalloc_fits` object.
#' @param reference Passed to [prevalence_ratios()] and [density_ratios()].
#' @param ... Unused.
#' @return A list of tibbles of class `pneumoalloc_summary`.
#' @export
summary.pneumoalloc_fit <- function(object, reference = "baseline", ...) {
  first <- if (inherits(object, "pneumoalloc_fits")) object[[1]] else object
  out <- list()
  if (first$type != "density") {
    out$prevalence <- marginal_prevalence(object)
    out$ratios <- prevalence_ratios(object, reference = reference)
  }
  if (first$type %in% c("unified", "unified_multiple", "density")) {
    out$density <- marginal_density(object)
    out$density_ratios <- density_ratios(object, reference = reference,
                                         multiplicity = TRUE)
  }
  structure(out, class = "pneumoalloc_summary")
}

#' @rdname summary.pneumoalloc_fit
#' @export
summary.pneumoalloc_fits <- function(object, reference = "baseline", ...) {
  summary.pneumoalloc_fit(object, reference = reference, ...)
}

#' @export
print.pneumoalloc_summary <- function(x, ...) {
  main <- x$prevalence[x$prevalence$group != "serotype", ]
  cli::cli_h2("Standardised estimates (%)")
  print(as.data.frame(main[, c("method", "measure", "outcome", "survey", "estimate")]),
        row.names = FALSE)
  main <- x$ratios[x$ratios$group != "serotype", ]
  cli::cli_h2("Ratios")
  print(as.data.frame(main[, c("method", "ratio", "outcome", "survey", "reference", "estimate")]),
        row.names = FALSE)
  cli::cli_text("")
  cli::cli_text("Tables in this list: {.field {names(x)}}.")
  invisible(x)
}

#' Write tables to an Excel workbook
#'
#' Writes each table to its own sheet. Needs the 'openxlsx' package.
#'
#' @param x A `pneumoalloc_summary`, a `pneumoalloc_fit(s)` object (which
#'   is summarised first) or a named list of data frames.
#' @param file Path of the `.xlsx` file.
#' @return `file`, invisibly.
#' @export
write_results <- function(x, file) {
  rlang::check_installed("openxlsx", reason = "to write Excel files.")
  if (inherits(x, c("pneumoalloc_fit", "pneumoalloc_fits"))) x <- summary(x)
  sheets <- Filter(Negate(is.null), unclass(x))
  openxlsx::write.xlsx(lapply(sheets, as.data.frame), file = file)
  invisible(file)
}

# Groups of serotype categories as a named list; check them against the
# categories and the reference set ("VT", "NVT" or "all").
resolve_groups <- function(groups, categories, vt_mask, among,
                           call = rlang::caller_env()) {
  if (is.character(groups)) {
    groups <- stats::setNames(list(groups), paste(groups, collapse = "+"))
  }
  if (!is.list(groups) || length(groups) == 0 ||
      !all(vapply(groups, is.character, TRUE))) {
    cli::cli_abort("{.arg groups} must be a character vector or a named list of them.",
                   call = call)
  }
  if (is.null(names(groups))) names(groups) <- rep("", length(groups))
  empty <- names(groups) == ""
  names(groups)[empty] <- vapply(groups[empty], paste, "", collapse = "+")
  pool <- switch(among, VT = categories[vt_mask == 1],
                 NVT = categories[vt_mask == 0], all = categories)
  for (nm in names(groups)) {
    unknown <- setdiff(groups[[nm]], categories)
    if (length(unknown) > 0) {
      cli::cli_abort(c(
        "Group {.val {nm}}: {.val {unknown}} {?is/are} not a serotype category of the fit.",
        "i" = "Categories: {.val {categories}}."
      ), call = call)
    }
    outside <- setdiff(groups[[nm]], pool)
    if (length(outside) > 0) {
      cli::cli_abort("Group {.val {nm}}: {.val {outside}} {?is/are} not in the {among} group.",
                     call = call)
    }
  }
  list(groups = groups, pool = pool)
}

#' Proportion of serotype groups within VT or NVT carriage
#'
#' The share of a group of serotypes in vaccine-type (`among = "VT"`),
#' non-vaccine-type (`"NVT"`) or all (`"all"`) serotype carriage, by survey
#' round, and its ratio between rounds: for example the proportion of
#' 15B/C, 15A and 23A among all non-PCV10 carriage before and after
#' vaccine introduction. Both are formed within each posterior draw as
#' `sum of the group's serotype prevalences / prevalence of the reference
#' set`, so they are standardised like the prevalences they come from.
#'
#' For fits among children with multiple carriage
#' ([fit_carriage()]`(population = "multiple")`) the proportions are within
#' multiple carriage. The presence-based model ([fit_presence()]) is not
#' supported: its serotype prevalences do not add up to VT or NVT
#' prevalence.
#'
#' @param object A `pneumoalloc_fit` or `pneumoalloc_fits` object from
#'   [fit_carriage()].
#' @param groups The serotype categories of a group, e.g.
#'   `c("15B/C", "15A", "23A")`, or a named list of several groups, e.g.
#'   `list("15B/C+15A+23A" = c("15B/C", "15A", "23A"), "19A" = "19A")`. Use
#'   the categories of the fit (pooled serotypes count as `"Other NVT"`).
#' @param among `"NVT"` (default), `"VT"` or `"all"`: the set the
#'   proportion is taken of. Every serotype of a group must belong to it.
#' @inheritParams prevalence_ratios
#' @return A tibble with, for each group, the proportion by survey round
#'   (`measure` "proportion among NVT", ...; `estimate` in percent) and the
#'   proportion ratio versus the reference round (`measure` "proportion
#'   ratio").
#' @seealso [describe_group_proportions()] for the crude proportions.
#' @export
#' @examples
#' \dontrun{
#' fit <- fit_carriage(dat, approach = "density")
#' group_proportions(fit, c("15B/C", "15A", "23A"), among = "NVT")
#' }
group_proportions <- function(object, groups, among = c("NVT", "VT", "all"),
                              reference = "baseline") {
  among <- rlang::arg_match(among)
  out <- bind_fits(object, group_proportions, groups = groups, among = among,
                   reference = reference)
  if (!is.null(out)) return(out)
  if (!object$type %in% c("unified", "unified_multiple")) {
    cli::cli_abort(c(
      "{.fn group_proportions} needs a fit of {.fn fit_carriage}.",
      "i" = if (object$type == "presence")
        "Presence-based serotype prevalences do not add up to VT or NVT prevalence."
    ))
  }
  cats <- object$labels$categories
  g <- resolve_groups(groups, cats, object$labels$vt_mask, among)
  share <- draws_array(object, "share")          # [draw, survey, category]
  denom <- apply(share[, , cats %in% g$pool, drop = FALSE], c(1, 2), sum)
  prop <- vapply(g$groups, function(s) {
    apply(share[, , cats %in% s, drop = FALSE], c(1, 2), sum) / denom
  }, matrix(0, dim(share)[1], dim(share)[2]))
  prop <- array(prop, c(dim(share)[1:2], length(g$groups)))

  lev <- object$data$survey_levels
  measure <- paste("proportion among", if (among == "all") "carriers" else among)
  e <- entry(prop, names(g$groups), "serotype group", measure)
  tab_p <- entry_table(object, e, lev)
  tab_p$reference <- NA_character_
  tab_p$estimate <- fmt_ci(100 * tab_p$median, 100 * tab_p$q2.5, 100 * tab_p$q97.5, 1)

  ref <- reference_index(reference, lev, object$data$reference)
  keep <- which(!is.na(ref))
  e$arr <- prop[, keep, , drop = FALSE] / prop[, ref[keep], , drop = FALSE]
  e$measure <- "proportion ratio"
  tab_r <- entry_table(object, e, lev[keep])
  tab_r$reference <- rep(lev[ref[keep]], times = length(g$groups))
  tab_r$estimate <- fmt_ci(tab_r$median, tab_r$q2.5, tab_r$q97.5, 2)
  tab <- rbind(tab_p, tab_r)
  tab$among <- among
  tibble::as_tibble(tab)
}
