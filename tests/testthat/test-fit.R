test_that("carriage states: prevalence, ratios and checks", {
  skip_if_no_cmdstan()
  dat <- example_data()
  fit <- do.call(fit_carriage_states, c(list(dat), quick))
  expect_s3_class(fit, "pneumoalloc_fit")

  est <- marginal_prevalence(fit)
  states <- est[est$group == "state", ]
  expect_equal(as.vector(tapply(states$mean, states$survey, sum)), rep(1, 3), tolerance = 1e-6)
  carriage <- est$median[est$group == "carriage"]
  crude <- describe_carriage_states(dat)
  expect_equal(carriage, crude$value[crude$outcome == "Carriage"], tolerance = 0.1)

  pr <- prevalence_ratios(fit, groups = c("carriage", "multiple"))
  expect_equal(nrow(pr), 4)
  expect_true(all(pr$reference == "2016"))
  prev <- prevalence_ratios(fit, groups = "carriage", reference = "previous")
  expect_equal(prev$reference, c("2016", "2017"))
  expect_equal(as.character(prev$survey), c("2017", "2019"))

  chk <- check_fit(fit)
  expect_s3_class(chk, "pneumoalloc_check")
  expect_true(all(chk$identities$pass))
})

test_that("types of multiple carriage", {
  skip_if_no_cmdstan()
  dat <- example_data()
  fit <- do.call(fit_multiple_types, c(list(dat, covariates = ~ age_z), quick))
  est <- marginal_prevalence(fit)
  expect_setequal(est$outcome, type_labels())
  expect_equal(as.vector(tapply(est$mean, est$survey, sum)), rep(1, 3), tolerance = 1e-6)
  expect_true(all(check_fit(fit)$identities$pass))
  expect_error(marginal_density(fit), "fit_carriage")
})

test_that("unified framework, all children, all approaches", {
  skip_if_no_cmdstan()
  dat <- example_data()
  fits <- do.call(fit_carriage, c(list(dat, approach = "all"), quick))
  expect_s3_class(fits, "pneumoalloc_fits")
  expect_named(fits, c("dominant", "equal", "abundance", "density"))

  for (f in fits) expect_true(all(check_fit(f)$identities$pass), info = f$approach$label)

  est <- marginal_prevalence(fits, groups = c("overall", "VT", "NVT"))
  prev <- est[est$measure == "prevalence", ]
  expect_equal(nrow(prev), 4 * 3 * 3)

  pr <- prevalence_ratios(fits$density, groups = "serotype")
  expect_true(all(c("reportable", "reference_q97.5", "vaccine_type") %in% names(pr)))
  expect_false(any(pr$survey == "2016"))

  # approach iv weights are normalised within each survey round
  carriers <- as.vector(table(dat$survey[dat$pos_idx]))
  w_round <- colMeans(posterior::as_draws_matrix(fits$density$fit$draws("dens_weight_sum_round")))
  expect_equal(unname(w_round), carriers, tolerance = 1e-8)
  ess_round <- colMeans(posterior::as_draws_matrix(fits$density$fit$draws("dens_weight_kish_ess_round")))
  expect_true(all(ess_round <= carriers + 1e-8))
  w_iii <- colMeans(posterior::as_draws_matrix(fits$abundance$fit$draws("dens_weight_kish_ess_round")))
  expect_equal(unname(w_iii), carriers, tolerance = 1e-8)
  ess_row <- check_fit(fits$density)$summary
  expect_true(any(grepl("effective sample size", ess_row$check)))

  cmp <- compare_approaches(fits)
  expect_false(any(cmp$shared$approach == 4))
  expect_true(all(abs(cmp$shared$z) < 6, na.rm = TRUE))

  d <- density_ratios(fits$density, multiplicity = TRUE)
  expect_equal(sum(is.na(d$survey)), 3)

  file <- withr::local_tempfile(fileext = ".rds")
  save_fit(fits, file)
  expect_equal(marginal_prevalence(readRDS(file)), marginal_prevalence(fits))
})

test_that("unified framework among multiple carriers", {
  skip_if_no_cmdstan()
  dat <- example_data()
  fits <- do.call(fit_carriage, c(list(dat, approach = c("equal", "density"),
                                       population = "multiple",
                                       covariates = ~ age_z + sex), quick))
  f <- fits$density
  expect_equal(f$type, "unified_multiple")
  est <- marginal_prevalence(f, groups = "serotype")
  expect_equal(as.vector(tapply(est$mean, est$survey, sum)), rep(1, 3), tolerance = 1e-6)
  expect_true(all(est$measure == "proportion among multiple carriers"))
  expect_true(all(check_fit(f)$identities$pass))
  expect_equal(nrow(density_ratios(f, multiplicity = TRUE)), 3)
  cmp <- compare_approaches(fits)
  expect_true(all(cmp$shared$quantity %in% c("P_m3", "meanD")))
})

test_that("presence-based model", {
  skip_if_no_cmdstan()
  dat <- example_data()
  fit <- do.call(fit_presence, c(list(dat), quick))
  expect_equal(fit$type, "presence")
  est <- marginal_prevalence(fit, groups = c("overall", "VT", "NVT", "multiple"))
  expect_setequal(est$outcome, c("Overall", "VT", "NVT", "Multiple carriage"))
  crude <- describe_prevalence(dat, "presence")
  vt <- est[est$outcome == "VT", ]
  expect_equal(vt$median, crude$value[crude$outcome == "VT"], tolerance = 0.15)
  chk <- check_fit(fit)
  expect_equal(nrow(chk$identities), 0)

  fits <- do.call(fit_carriage, c(list(dat, approach = c("dominant", "density")), quick))
  cmp <- compare_approaches(c(fits, list(presence = fit)))
  expect_true("Presence-based" %in% cmp$ratios$method)
  expect_false("Presence-based" %in% cmp$shared$method)
})

test_that("density models: groups, statistics, interaction", {
  skip_if_no_cmdstan()
  dat <- example_data()
  fits <- do.call(fit_density, c(list(dat, covariates = ~ age_z, m_rep = 1000), quick))
  expect_named(fits, c("total", "VT", "NVT"))

  md <- marginal_density(fits$total)
  expect_setequal(md$statistic, c("mean", "median"))
  med <- md[md$outcome == "Density" & md$statistic == "median", ]
  mn <- md[md$outcome == "Density" & md$statistic == "mean", ]
  expect_true(all(mn$median > med$median))   # lognormal: mean above median

  dr <- density_ratios(fits, statistic = "median", multiplicity = TRUE)
  expect_equal(sum(dr$outcome == "Density"), 3 * 2)
  expect_equal(sum(is.na(dr$survey)), 3 * 3)
  prev <- density_ratios(fits$NVT, reference = "previous", statistic = "median")
  expect_equal(prev$reference, c("2016", "2017"))
  expect_error(marginal_prevalence(fits$total), "marginal_density")

  fi <- do.call(fit_density, c(list(dat, group = "total", interaction = TRUE,
                                    m_rep = 1000), quick))
  di <- density_ratios(fi, multiplicity = TRUE, statistic = "median")
  expect_equal(sum(di$measure == "ratio of density ratios"), 3 * 2)
  expect_equal(nrow(check_fit(fi)$crude), 3)
})
