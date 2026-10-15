#' Fit the unified framework
#'
#' Fits the unified Bayesian framework for one or more allocation
#' approaches: overall carriage (logistic, all children), serotype
#' composition among carriers (fractional multinomial logit) and density
#' among carriers (lognormal, with a sub-model for the number of
#' serotypes), fitted jointly. All fits share the population, design,
#' priors and sampler settings; only the allocation weights (and, for the
#' `density` approach, the carrier weight) differ.
#'
#' With `population = "multiple"`, the same framework is fitted among
#' children with multiple serotype carriage only. There is then no
#' carriage component: the estimand is the serotype composition of
#' multiple carriage, and the multiplicity sub-model compares three with
#' two serotypes.
#'
#' @param dat A [prepare_carriage_data()] object.
#' @param approach Which allocation approach(es) to fit, chosen by the aim
#'   of the analysis: `"all"` (default; all four, e.g. to compare them),
#'   one approach (`"dominant"`, `"equal"`, `"abundance"` or `"density"`),
#'   or several, e.g. `c("dominant", "density")`. Numbers 1-4 also work;
#'   see [allocation_approaches()].
#' @param covariates One-sided formula of covariates for an adjusted model.
#'   Any number of covariates can be used, numeric or factor, with the usual
#'   formula syntax: `~ age_z + sex + household + income_z`, transformations
#'   such as `~ age_z + I(age_z^2)`, or interactions such as `~ age_z * sex`.
#'   Survey round is always included, with one intercept per round
#'   (`~ 0 + survey + ...`). `NULL` (default) fits the unadjusted model.
#'   Variables must be declared in the `covariates` of
#'   [prepare_carriage_data()]. Avoid interactions with survey round: the
#'   prevalence ratios assume the survey effect is the same for everyone.
#' @param population `"all"` (all children) or `"multiple"` (children with
#'   multiple serotype carriage).
#' @param priors Prior scales from [carriage_priors()].
#' @param chains,parallel_chains,iter_warmup,iter_sampling,threads_per_chain,seed,max_treedepth,adapt_delta,refresh
#'   Sampler settings passed to the `$sample()` method of
#'   [cmdstanr::CmdStanModel]. The defaults (4 chains of 1000 warm-up and
#'   2000 sampling iterations) are those of the PCV10 Nha Trang analysis.
#' @param model A compiled model from [carriage_model()]. Compiled (or
#'   loaded from the cache) when `NULL`.
#' @param ... Further arguments for `$sample()`, e.g. `output_dir`.
#'
#' @return For one approach, a `pneumoalloc_fit` object: a list holding the
#'   `CmdStanMCMC` object (`fit`), the analysis data (`data`), the model
#'   type, approach, covariates and settings. For several approaches, a
#'   named list of `pneumoalloc_fit` objects of class `pneumoalloc_fits`.
#'
#' @seealso [marginal_prevalence()], [prevalence_ratios()],
#'   [density_ratios()], [check_fit()], [compare_approaches()],
#'   [save_fit()].
#' @export
#' @examples
#' \dontrun{
#' sim <- simulate_carriage_data(seed = 1)
#' dat <- prepare_carriage_data(sim, "survey", c("sero1", "sero2", "sero3"),
#'   density = "density", abundance = c("ra1", "ra2", "ra3"),
#'   covariates = c("age_z", "sex"))
#' # one approach
#' fit <- fit_carriage(dat, approach = "density", covariates = ~ age_z + sex)
#' prevalence_ratios(fit, groups = c("VT", "NVT"))
#'
#' # all four, to compare them
#' fits <- fit_carriage(dat, approach = "all", covariates = ~ age_z + sex)
#' compare_approaches(fits)
#' }
fit_carriage <- function(dat, approach = "all", covariates = NULL,
                         population = c("all", "multiple"),
                         priors = carriage_priors(), chains = 4,
                         parallel_chains = chains, iter_warmup = 1000,
                         iter_sampling = 2000, threads_per_chain = 1,
                         seed = 123, max_treedepth = 12, adapt_delta = 0.8,
                         refresh = 500, model = NULL, ...) {
  check_class(dat, "carriage_data")
  population <- rlang::arg_match(population)
  approach <- expand_approach(approach)
  type <- if (population == "all") "unified" else "unified_multiple"
  # check every approach's data before any (slow) compilation or sampling
  preps <- lapply(approach, function(a) {
    prepare_model_data(dat, type, a, covariates, priors)
  })
  model <- model %||% carriage_model(type)
  sampler <- list(chains = chains, parallel_chains = parallel_chains,
                  iter_warmup = iter_warmup, iter_sampling = iter_sampling,
                  threads_per_chain = threads_per_chain, seed = seed,
                  max_treedepth = max_treedepth, adapt_delta = adapt_delta,
                  refresh = refresh)

  fits <- lapply(preps, function(prep) {
    cli::cli_h2("Approach {prep$approach$approach}: {prep$approach$label}")
    run_model(model, prep, dat, covariates, priors, sampler, ...)
  })
  if (length(fits) == 1) return(fits[[1]])
  names(fits) <- vapply(preps, function(p) p$approach$name, "")
  structure(fits, class = "pneumoalloc_fits")
}

#' Fit overall and multiple serotype carriage
#'
#' Multinomial-logit model in which each child is classified by the number
#' of detected serotypes (no carriage, one, two, three, ...; reference: no
#' carriage). Survey-specific prevalences of each carriage state, of
#' carriage (at least one serotype) and of multiple serotype carriage (at
#' least two) are standardised over the study population.
#'
#' @inheritParams fit_carriage
#' @return A `pneumoalloc_fit` object.
#' @seealso [describe_carriage_states()] for the crude proportions.
#' @export
#' @examples
#' \dontrun{
#' sim <- simulate_carriage_data(seed = 1)
#' dat <- prepare_carriage_data(sim, "survey", c("sero1", "sero2", "sero3"),
#'   density = "density")
#' fit <- fit_carriage_states(dat)
#' prevalence_ratios(fit)
#' }
fit_carriage_states <- function(dat, covariates = NULL,
                                priors = carriage_priors(), chains = 4,
                                parallel_chains = chains, iter_warmup = 1000,
                                iter_sampling = 2000, threads_per_chain = 1,
                                seed = 123, max_treedepth = 12,
                                adapt_delta = 0.8, refresh = 500,
                                model = NULL, ...) {
  prep <- prepare_model_data(dat, "states", covariates = covariates, priors = priors)
  sampler <- list(chains = chains, parallel_chains = parallel_chains,
                  iter_warmup = iter_warmup, iter_sampling = iter_sampling,
                  threads_per_chain = threads_per_chain, seed = seed,
                  max_treedepth = max_treedepth, adapt_delta = adapt_delta,
                  refresh = refresh)
  run_model(model %||% carriage_model("categorical"), prep, dat, covariates,
            priors, sampler, ...)
}

#' Fit the types of multiple serotype carriage
#'
#' Among children with multiple serotype carriage, multinomial-logit model
#' for four types: VT-only (all detected serotypes are vaccine types),
#' NVT-only, mixed VT-dominant and mixed NVT-dominant (both kinds detected;
#' classified by the dominant serotype; reference: mixed NVT-dominant).
#' Survey-specific proportions are standardised over the children with
#' multiple carriage.
#'
#' The target acceptance rate `adapt_delta` is 0.95 by default, as in the
#' PCV10 Nha Trang analysis, because some type-by-survey cells can be
#' empty.
#'
#' @inheritParams fit_carriage
#' @return A `pneumoalloc_fit` object.
#' @seealso [describe_carriage_states()] for the crude proportions.
#' @export
fit_multiple_types <- function(dat, covariates = NULL,
                               priors = carriage_priors(), chains = 4,
                               parallel_chains = chains, iter_warmup = 1000,
                               iter_sampling = 2000, threads_per_chain = 1,
                               seed = 123, max_treedepth = 12,
                               adapt_delta = 0.95, refresh = 500,
                               model = NULL, ...) {
  prep <- prepare_model_data(dat, "multiple_types", covariates = covariates,
                             priors = priors)
  sampler <- list(chains = chains, parallel_chains = parallel_chains,
                  iter_warmup = iter_warmup, iter_sampling = iter_sampling,
                  threads_per_chain = threads_per_chain, seed = seed,
                  max_treedepth = max_treedepth, adapt_delta = adapt_delta,
                  refresh = refresh)
  run_model(model %||% carriage_model("categorical"), prep, dat, covariates,
            priors, sampler, ...)
}

#' Fit the presence-based model
#'
#' The presence-based definition used in many carriage studies, for
#' comparison with the allocation approaches. Each outcome is a child-level
#' indicator with its own logistic regression on the same design: overall
#' carriage (at least one serotype), VT and NVT carriage (at least one
#' vaccine or non-vaccine serotype), multiple carriage (at least two
#' serotypes) and each serotype category, irrespective of dominance,
#' relative abundance or density. A child carrying both a VT and an NVT
#' serotype counts for both, so VT and NVT prevalence can add up to more
#' than overall carriage. Prevalences are standardised over the study
#' population.
#'
#' @inheritParams fit_carriage
#' @return A `pneumoalloc_fit` object, with method `"Presence-based"`. It
#'   can be passed to [compare_approaches()] together with the fits of
#'   [fit_carriage()].
#' @seealso [describe_prevalence()] with `approach = "presence"` for the
#'   crude proportions.
#' @export
#' @examples
#' \dontrun{
#' sim <- simulate_carriage_data(seed = 1)
#' dat <- prepare_carriage_data(sim, "survey", c("sero1", "sero2", "sero3"),
#'   density = "density")
#' fit <- fit_presence(dat)
#' prevalence_ratios(fit, groups = c("VT", "NVT"))
#' }
fit_presence <- function(dat, covariates = NULL, priors = carriage_priors(),
                         chains = 4, parallel_chains = chains,
                         iter_warmup = 1000, iter_sampling = 2000,
                         threads_per_chain = 1, seed = 123, max_treedepth = 12,
                         adapt_delta = 0.8, refresh = 500, model = NULL, ...) {
  prep <- prepare_model_data(dat, "presence", covariates = covariates, priors = priors)
  sampler <- list(chains = chains, parallel_chains = parallel_chains,
                  iter_warmup = iter_warmup, iter_sampling = iter_sampling,
                  threads_per_chain = threads_per_chain, seed = seed,
                  max_treedepth = max_treedepth, adapt_delta = adapt_delta,
                  refresh = refresh)
  run_model(model %||% carriage_model("presence"), prep, dat, covariates,
            priors, sampler, ...)
}

#' Fit density among carriers
#'
#' Lognormal regression of density among carriers on the natural-log scale,
#' with nested effects for carriage of at least two and at least three
#' serotypes, for total pneumococcal density and for the density of
#' vaccine-type (`"VT"`) and non-vaccine-type (`"NVT"`) serotypes. Group
#' density is the sum of the serotype-specific densities (total density
#' times relative abundance) of the group's serotypes, among children
#' carrying at least one serotype of the group.
#'
#' The model reports, by survey round, the mean density (analytic) and the
#' median density (from posterior predictive draws, with the observed
#' number of serotypes held fixed across rounds); see [marginal_density()]
#' and [density_ratios()]. The PCV10 Nha Trang analysis reports median
#' density ratios.
#'
#' @inheritParams fit_carriage
#' @param group `"all"` (default; total, VT and NVT) or one or more of
#'   `"total"`, `"VT"` and `"NVT"`.
#' @param interaction Let the multiplicity effects differ by survey round
#'   (survey x multiplicity interaction)?
#' @param multiplicity For the VT and NVT models, count the serotypes of the
#'   group (`"group"`, default: a child with one VT and two NVT serotypes
#'   has one serotype in the VT model) or all detected serotypes
#'   (`"total"`). These are different estimands. The total model always
#'   uses all detected serotypes.
#' @param m_rep Posterior predictive draws per survey round for the median
#'   and quantiles.
#' @return For one group a `pneumoalloc_fit` object; for several a named
#'   list of class `pneumoalloc_fits`.
#' @seealso [describe_density()] for the observed density.
#' @export
#' @examples
#' \dontrun{
#' sim <- simulate_carriage_data(seed = 1)
#' dat <- prepare_carriage_data(sim, "survey", c("sero1", "sero2", "sero3"),
#'   density = "density", abundance = c("ra1", "ra2", "ra3"))
#' fits <- fit_density(dat)
#' density_ratios(fits, statistic = "median")
#' }
fit_density <- function(dat, group = "all", covariates = NULL,
                        interaction = FALSE,
                        multiplicity = c("group", "total"), m_rep = 4000,
                        priors = carriage_priors(), chains = 4,
                        parallel_chains = chains, iter_warmup = 1000,
                        iter_sampling = 2000, threads_per_chain = 1,
                        seed = 123, max_treedepth = 12, adapt_delta = 0.8,
                        refresh = 500, model = NULL, ...) {
  multiplicity <- rlang::arg_match(multiplicity)
  groups <- expand_density_group(group)
  preps <- lapply(groups, function(g) {
    prepare_model_data(dat, "density", covariates = covariates, priors = priors,
                       density = list(group = g, interaction = interaction,
                                      multiplicity = multiplicity, m_rep = m_rep))
  })
  model <- model %||% carriage_model("density")
  sampler <- list(chains = chains, parallel_chains = parallel_chains,
                  iter_warmup = iter_warmup, iter_sampling = iter_sampling,
                  threads_per_chain = threads_per_chain, seed = seed,
                  max_treedepth = max_treedepth, adapt_delta = adapt_delta,
                  refresh = refresh)
  fits <- lapply(preps, function(prep) {
    cli::cli_h2("Density: {density_group_label(prep$labels$group)}")
    run_model(model, prep, dat, covariates, priors, sampler, ...)
  })
  if (length(fits) == 1) return(fits[[1]])
  names(fits) <- groups
  structure(fits, class = "pneumoalloc_fits")
}

expand_density_group <- function(group, call = rlang::caller_env()) {
  all_groups <- c("total", "VT", "NVT")
  if (identical(group, "all")) return(all_groups)
  bad <- setdiff(group, all_groups)
  if (length(group) == 0 || length(bad) > 0) {
    cli::cli_abort("{.arg group} must be {.val all} or one or more of {.val {all_groups}}.",
                   call = call)
  }
  unique(group)
}

density_group_label <- function(group) {
  c(total = "all pneumococci", VT = "vaccine-type serotypes",
    NVT = "non-vaccine-type serotypes")[[group]]
}

run_model <- function(model, prep, dat, covariates, priors, sampler, ...) {
  # One BLAS/OpenMP thread per process; parallelism comes from chains and
  # reduce_sum threads.
  withr::local_envvar(c(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
                        MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1"))
  fit <- do.call(model$sample, c(list(data = prep$stan_data), sampler, list(...)))
  structure(
    list(
      fit = fit,
      data = dat,
      type = prep$model,
      approach = prep$approach,
      covariates = covariates,
      rows = prep$rows,
      design = prep$design,
      labels = prep$labels,
      allocation = prep$allocation,
      priors = priors,
      sampler = sampler
    ),
    class = "pneumoalloc_fit"
  )
}

model_label <- function(object) {
  switch(
    object$type,
    unified = "Unified framework, all children",
    unified_multiple = "Unified framework, children with multiple carriage",
    states = "Overall and multiple serotype carriage",
    multiple_types = "Types of multiple serotype carriage",
    presence = "Presence-based carriage",
    density = paste("Density among carriers,", density_group_label(object$labels$group))
  )
}

# Short label of one fit within a set of fits
fit_label <- function(object) {
  if (object$type == "density") return(density_group_label(object$labels$group))
  if (object$type == "presence") return("Presence-based")
  if (!is.null(object$approach)) {
    return(paste(object$approach$approach, object$approach$label))
  }
  model_label(object)
}

#' Save a fit with its draws
#'
#' A `CmdStanMCMC` object reads its draws lazily from CSV files that may
#' be temporary. `save_fit()` loads the draws and sampler diagnostics into
#' memory and then calls [saveRDS()], so the saved file can be summarised
#' later, or on another machine, after [readRDS()].
#'
#' @param object A `pneumoalloc_fit` or `pneumoalloc_fits` object.
#' @param file Path of the `.rds` file.
#' @return `file`, invisibly.
#' @export
save_fit <- function(object, file) {
  fits <- if (inherits(object, "pneumoalloc_fits")) object else list(object)
  for (f in fits) {
    check_class(f, "pneumoalloc_fit", arg = "object")
    f$fit$draws()
    try(f$fit$sampler_diagnostics(), silent = TRUE)
    try(f$fit$init(), silent = TRUE)
  }
  saveRDS(object, file)
  invisible(file)
}

#' @export
print.pneumoalloc_fit <- function(x, ...) {
  diag <- x$fit$diagnostic_summary(quiet = TRUE)
  cli::cli_h1("pneumoalloc fit: {model_label(x)}")
  if (x$type %in% c("unified", "unified_multiple")) {
    cli::cli_text("Approach {x$approach$approach}: {x$approach$label}.")
  }
  cli::cli_text(
    "{if (is.null(x$covariates)) 'Unadjusted' else 'Adjusted'} model ",
    "({.code {formula_text(x$covariates)}}); {length(x$rows)} ",
    "{if (x$type == 'density') 'carriers' else 'children'}."
  )
  cli::cli_text(
    "{x$sampler$chains} chain{?s} x {x$sampler$iter_sampling} draws; ",
    "{sum(diag$num_divergent)} divergent transition{?s}, ",
    "{sum(diag$num_max_treedepth)} at max treedepth."
  )
  cli::cli_text("Ratios versus {.val {x$data$reference_level}}:")
  if (x$type == "density") {
    pr <- density_ratios(x, statistic = "median")
    print(as.data.frame(pr[, c("statistic", "outcome", "survey", "estimate")]),
          row.names = FALSE)
    cli::cli_text("See {.fn marginal_density}, {.fn density_ratios} and {.fn check_fit}.")
  } else {
    groups <- switch(x$type, unified = c("VT", "NVT"),
                     unified_multiple = c("VT", "NVT"),
                     presence = c("VT", "NVT"),
                     states = c("carriage", "multiple"),
                     multiple_types = "type")
    pr <- prevalence_ratios(x, groups = groups)
    print(as.data.frame(pr[, c("measure", "outcome", "survey", "estimate")]),
          row.names = FALSE)
    cli::cli_text("See {.fn marginal_prevalence}, {.fn prevalence_ratios} and {.fn check_fit}.")
  }
  invisible(x)
}

#' @export
print.pneumoalloc_fits <- function(x, ...) {
  cli::cli_h1("pneumoalloc: {length(x)} fits")
  for (f in x) {
    diag <- f$fit$diagnostic_summary(quiet = TRUE)
    cli::cli_li(
      "{fit_label(f)}: {sum(diag$num_divergent)} divergent transition{?s}"
    )
  }
  cli::cli_text("Access one fit with {.code fits${names(x)[1]}}.")
  invisible(x)
}
