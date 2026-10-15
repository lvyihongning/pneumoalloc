# The package's unified model (a faster implementation) must have the same
# log density as the published Stan file (fixtures/unified_published.stan). Compiling both with model methods takes a few minutes, so
# this runs only with PNEUMOALLOC_EXTENDED_TESTS=true.
test_that("unified.stan has the log density of the published model", {
  skip_if_no_cmdstan()
  skip_if_not(identical(Sys.getenv("PNEUMOALLOC_EXTENDED_TESTS"), "true"),
              "extended tests not requested")
  skip_if_not_installed("Rcpp")
  dat <- example_data()
  dir <- withr::local_tempdir()
  compile <- function(file) {
    cmdstanr::cmdstan_model(file, dir = dir, quiet = TRUE,
                            cpp_options = list(stan_threads = TRUE),
                            compile_model_methods = TRUE)
  }
  orig <- compile(test_path("fixtures", "unified_published.stan"))
  fast <- compile(carriage_stan_file("unified"))
  for (a in c("abundance", "density")) {
    sd <- carriage_stan_data(dat, "unified", a, covariates = ~ age_z + sex)
    # the published file also takes commune data; one commune, no effects
    sd0 <- c(sd, list(J = 1L, commune_id = rep(1L, sd$N), use_re = 0L,
                      prior_sd_sig_car = 1, prior_sd_sig_share = 2,
                      prior_sd_sig_dens = 1, prior_sd_sig_mult = 1))
    args <- list(seed = 1, chains = 1, threads_per_chain = 1,
                 refresh = 0, show_messages = FALSE)
    f0 <- do.call(orig$sample, c(list(data = sd0), args,
                                 iter_warmup = 50, iter_sampling = 10))
    f1 <- do.call(fast$sample, c(list(data = sd), args, iter_warmup = 3, iter_sampling = 3))
    f0$init_model_methods()
    f1$init_model_methods()
    up <- f0$unconstrain_draws(format = "draws_matrix")
    for (i in seq_len(nrow(up))) {
      expect_equal(f1$log_prob(up[i, ]), f0$log_prob(up[i, ]), tolerance = 1e-10)
      expect_equal(f1$grad_log_prob(up[i, ]), f0$grad_log_prob(up[i, ]),
                   tolerance = 1e-8, ignore_attr = TRUE)
    }
  }
})
