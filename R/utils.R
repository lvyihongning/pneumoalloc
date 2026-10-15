`%||%` <- function(x, y) if (is.null(x)) y else x

# "median (lower, upper)" as text
fmt_ci <- function(med, lo, hi, digits = 2) {
  f <- function(v) formatC(v, digits = digits, format = "f")
  paste0(f(med), " (", f(lo), ", ", f(hi), ")")
}

# "~ 0 + survey + age + sex" for a covariate formula (NULL = unadjusted)
formula_text <- function(covariates) {
  tl <- if (is.null(covariates)) character() else attr(stats::terms(covariates), "term.labels")
  paste(c("~ 0 + survey", tl), collapse = " + ")
}

check_class <- function(x, class, arg = rlang::caller_arg(x),
                        call = rlang::caller_env()) {
  if (!inherits(x, class)) {
    cli::cli_abort("{.arg {arg}} must be a {.cls {class}} object.", call = call)
  }
}

check_cmdstanr <- function(call = rlang::caller_env()) {
  rlang::check_installed(
    "cmdstanr",
    reason = "to compile and fit the Stan models.",
    action = function(...) {
      utils::install.packages(
        "cmdstanr",
        repos = c("https://stan-dev.r-universe.dev", getOption("repos"))
      )
    },
    call = call
  )
  version <- tryCatch(cmdstanr::cmdstan_version(error_on_NA = FALSE),
                      error = function(e) NULL)
  if (is.null(version)) {
    cli::cli_abort(c(
      "CmdStan is not installed or cannot be found.",
      "i" = "Install it once with {.code cmdstanr::install_cmdstan()}."
    ), call = call)
  }
  invisible(version)
}

# Index of the reference survey for each survey round; NA where a round
# has no reference (the reference round itself, or the first round when
# comparing with the previous round).
reference_index <- function(reference, survey_levels, default,
                            call = rlang::caller_env()) {
  n <- length(survey_levels)
  if (length(reference) != 1 || is.na(reference)) {
    cli::cli_abort("{.arg reference} must be a single value.", call = call)
  }
  if (identical(reference, "previous")) {
    return(c(NA_integer_, seq_len(n - 1)))
  }
  ref <- if (identical(reference, "baseline")) default else match(reference, survey_levels)
  if (is.na(ref)) {
    cli::cli_abort(c(
      "{.arg reference} must be {.val baseline}, {.val previous} or a survey round.",
      "i" = "Survey rounds: {.val {survey_levels}}."
    ), call = call)
  }
  idx <- rep(ref, n)
  idx[ref] <- NA_integer_
  idx
}

# posterior warns when an ESS estimate exceeds its cap (very efficient
# sampling); the capped value is still returned.
quiet_ess <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl("ESS has been capped", conditionMessage(w))) {
      invokeRestart("muffleWarning")
    }
  })
}
