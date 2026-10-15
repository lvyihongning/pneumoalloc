data_block_names <- function(model) {
  code <- readLines(carriage_stan_file(model))
  block <- code[seq(grep("^data \\{", code), grep("^transformed data", code))]
  block <- sub("//.*", "", block)
  sub(";.*", "", regmatches(block, regexpr("[A-Za-z_]+;", block)))
}

test_that("Stan data match the data blocks of the models", {
  dat <- example_data()
  files <- c(unified = "unified", unified_multiple = "unified_multiple",
             states = "categorical", multiple_types = "categorical")
  for (m in names(files)) {
    sd <- carriage_stan_data(dat, m, approach = "density", covariates = ~ age_z + sex)
    expect_setequal(names(sd), data_block_names(files[[m]]))
  }

})

test_that("designs have one column per survey round first", {
  dat <- example_data()
  sd <- carriage_stan_data(dat, "unified", "equal", covariates = ~ age_z + sex)
  n_rounds <- length(dat$survey_levels)
  expect_equal(sd$K, n_rounds + 2)
  for (t in seq_len(n_rounds)) {
    expect_true(all(sd$X_pred[[t]][, t] == 1))
    expect_equal(sd$X_pred[[t]][, -seq_len(n_rounds)], sd$X[, -seq_len(n_rounds)])
  }
  expect_equal(carriage_stan_data(dat, "states")$K, n_rounds)
  expect_equal(carriage_stan_data(dat, "states")$y, dat$n_sero + 1L)
  ty <- carriage_stan_data(dat, "multiple_types")
  expect_equal(ty$N, sum(dat$m_obs >= 2))
  expect_equal(ty$ref, 4L)
})

test_that("priors are passed through and validated", {
  dat <- example_data()
  pr <- carriage_priors(share = 5, states = 1, multiple_types = 3)
  expect_equal(carriage_stan_data(dat, "unified", 1, priors = pr)$prior_sd_share, 5)
  expect_equal(carriage_stan_data(dat, "unified_multiple", 1, priors = pr)$prior_sd_share, 5)
  expect_equal(carriage_stan_data(dat, "states", priors = pr)$prior_sd, 1)
  expect_equal(carriage_stan_data(dat, "multiple_types", priors = pr)$prior_sd, 3)
  expect_error(carriage_priors(car = -1), "positive")
  expect_error(carriage_priors(dens_center = NA), "finite")
})

test_that("covariates are validated", {
  dat <- example_data()
  expect_error(carriage_stan_data(dat, "unified", 1, covariates = ~ income), "not declared")
  expect_error(carriage_stan_data(dat, "unified", 1, covariates = y ~ age_z), "one-sided")
})

test_that("reference rounds are resolved", {
  lev <- c("a", "b", "c")
  expect_equal(reference_index("baseline", lev, 1L), c(NA, 1L, 1L))
  expect_equal(reference_index("previous", lev, 1L), c(NA, 1L, 2L))
  expect_equal(reference_index("b", lev, 1L), c(2L, NA, 2L))
  expect_error(reference_index("z", lev, 1L), "survey round")
})

test_that("simulate_carriage_data() is reproducible and leaves the RNG alone", {
  set.seed(10)
  before <- .Random.seed
  a <- simulate_carriage_data(seed = 5)
  expect_identical(.Random.seed, before)
  expect_identical(a, simulate_carriage_data(seed = 5))
})

test_that("approaches can be chosen by aim: one, several or all", {
  expect_equal(expand_approach("all"), c("dominant", "equal", "abundance", "density"))
  expect_equal(expand_approach("dominant"), "dominant")
  expect_equal(expand_approach(c("density", "equal")), c("density", "equal"))
  expect_equal(expand_approach(c(1, 4)), c("dominant", "density"))
  expect_equal(expand_approach(c("equal", 2)), "equal")
  expect_error(expand_approach(c("all", "dominant")), "on its own")
  expect_error(expand_approach("presence"), "must be one of")
  expect_error(expand_approach(character()), "empty")
})

test_that("presence and density data are built", {
  dat <- example_data()
  sd <- carriage_stan_data(dat, "presence", covariates = ~ age_z)
  expect_equal(dim(sd$y), c(4 + length(dat$sero_levels), dat$N))
  expect_equal(sd$y[1, ], dat$carrier)
  expect_true(all(sd$y[1, ] >= pmax(sd$y[2, ], sd$y[3, ])))
  expect_equal(sd$prior_sd, 2.5)

  for (g in c("total", "VT", "NVT")) {
    d <- carriage_stan_data(dat, "density", group = g, covariates = ~ age_z,
                            interaction = TRUE)
    expect_equal(d$mult_by_survey, 1L)
    expect_equal(d$beta_prior_mean, c(rep(log(1e6), 3), 0))
    expect_equal(d$beta_prior_sd, c(3, 3, 3, 1))
    expect_true(all(d$m_sh >= 1 & d$m_sh <= d$M))
  }
  expect_equal(carriage_stan_data(dat, "density", group = "total")$N_pos, dat$N_pos)
  expect_error(carriage_stan_data(dat, "density", group = "PCV13"), "must be one of")
  expect_equal(expand_density_group("all"), c("total", "VT", "NVT"))
  expect_error(expand_density_group(c("VT", "x")), "must be")

  sim <- simulate_carriage_data(seed = 1)
  no_ab <- prepare_carriage_data(sim, "survey", c("sero1", "sero2", "sero3"), density = "density")
  expect_error(carriage_stan_data(no_ab, "density", group = "VT"), "relative abundance")
})

test_that("a covariate with a single value in the population gives a clear error", {
  sim <- simulate_carriage_data(seed = 1)
  sim$grp <- factor(ifelse(rowSums(!is.na(sim[, c("sero1", "sero2", "sero3")])) >= 2, "a", "b"))
  dat <- suppressMessages(prepare_carriage_data(
    sim, "survey", c("sero1", "sero2", "sero3"), density = "density",
    categories = c(pcv10_types(), "6A", "15B/C"), covariates = "grp"
  ))
  expect_error(carriage_stan_data(dat, "unified_multiple", "equal", ~ grp), "single value")
  expect_no_error(carriage_stan_data(dat, "unified", "equal", ~ grp))
})

test_that("density models keep the multiplicity levels observed in the group", {
  dat <- example_data()
  for (g in c("total", "VT", "NVT")) {
    sd <- carriage_stan_data(dat, "density", group = g)
    expect_equal(sd$M, max(2L, sd$m_sh))
  }
})

test_that("a rare reference category gives a warning", {
  A <- cbind(a = c(1, 0.5, 1), b = c(0, 0.5, 0))
  expect_warning(check_reference_category(A, "equal"), "allocation of only 0.5")
  expect_no_warning(check_reference_category(A[rep(1:3, 10), 2:1], "equal"))
})
