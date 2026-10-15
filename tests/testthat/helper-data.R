example_data <- function(seed = 1, ...) {
  sim <- simulate_carriage_data(seed = seed, ...)
  suppressMessages(prepare_carriage_data(
    sim,
    survey = "survey",
    serotypes = c("sero1", "sero2", "sero3"),
    abundance = c("ra1", "ra2", "ra3"),
    density = "density",
    vaccine_types = pcv10_types(),
    categories = c(pcv10_types(), "6A", "15B/C", "19A", "11A"),
    covariates = c("age_z", "sex")
  ))
}

skip_if_no_cmdstan <- function() {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")
  version <- tryCatch(cmdstanr::cmdstan_version(error_on_NA = FALSE),
                      error = function(e) NULL)
  if (is.null(version)) skip("CmdStan is not installed")
}

# short runs: they test the plumbing, not the inference
quick <- list(chains = 2, iter_warmup = 150, iter_sampling = 150,
              refresh = 0, show_messages = FALSE)
