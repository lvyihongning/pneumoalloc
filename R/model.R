#' Stan models
#'
#' `carriage_stan_file()` gives the path of one of the Stan programs
#' shipped with the package. `carriage_model()` compiles it with threading
#' enabled through 'cmdstanr'. Executables are cached in a user directory,
#' so each model is compiled once per package version.
#'
#' * `"unified"`: the unified framework in all children
#'   ([fit_carriage()], `population = "all"`). The model of the PCV10 Nha
#'   Trang analysis, with a faster likelihood computation (see the header
#'   of the Stan file);
#' * `"unified_multiple"`: the unified framework among children with
#'   multiple carriage ([fit_carriage()], `population = "multiple"`);
#' * `"categorical"`: the multinomial-logit model of
#'   [fit_carriage_states()] and [fit_multiple_types()];
#' * `"presence"`: the presence-based model of [fit_presence()];
#' * `"density"`: the density model of [fit_density()].
#'
#' @param model Which model.
#' @param dir Directory for the compiled executable. Defaults to a
#'   per-version folder in [tools::R_user_dir()]. Set it to a shared
#'   directory on a computing cluster.
#' @param force_recompile Recompile even if a cached executable exists.
#' @param quiet Suppress compiler output.
#' @param ... Passed to [cmdstanr::cmdstan_model()].
#'
#' @return `carriage_stan_file()`: a file path. `carriage_model()`: a
#'   `CmdStanModel` object.
#' @export
#' @examples
#' carriage_stan_file("unified")
#' \dontrun{
#' mod <- carriage_model("unified")
#' }
carriage_model <- function(model = c("unified", "unified_multiple", "categorical",
                                     "presence", "density"),
                           dir = NULL, force_recompile = FALSE, quiet = TRUE, ...) {
  model <- rlang::arg_match(model)
  check_cmdstanr()
  dir <- dir %||% file.path(
    tools::R_user_dir("pneumoalloc", which = "cache"),
    as.character(utils::packageVersion("pneumoalloc"))
  )
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  cmdstanr::cmdstan_model(
    carriage_stan_file(model),
    dir = dir,
    cpp_options = list(stan_threads = TRUE),
    force_recompile = force_recompile,
    quiet = quiet,
    ...
  )
}

#' @rdname carriage_model
#' @export
carriage_stan_file <- function(model = c("unified", "unified_multiple", "categorical",
                                     "presence", "density")) {
  model <- rlang::arg_match(model)
  system.file("stan", paste0(model, ".stan"), package = "pneumoalloc",
              mustWork = TRUE)
}
