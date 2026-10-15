test_that("vaccine sets can be chosen and modified", {
  expect_length(vaccine_types("PCV10"), 10)
  expect_length(vaccine_types("PCV13"), 13)
  expect_length(vaccine_types("PCV15"), 15)
  expect_length(vaccine_types("PCV20"), 20)
  expect_length(vaccine_types("PCV21"), 21)
  expect_length(vaccine_types("PPSV23"), 23)
  expect_true(all(pcv10_types() %in% pcv13_types()))
  expect_setequal(setdiff(pcv15_types(), pcv13_types()), c("22F", "33F"))
  expect_true("6A" %in% vaccine_types("PCV10", add = "6A"))
  expect_false("15B" %in% vaccine_types("PCV20", add = "15B/C", remove = "15B"))
  expect_setequal(vaccine_types(c("PCV7", "PCV10")), pcv10_types())
  expect_error(vaccine_types("PCV99"), "Unknown vaccine")
})

test_that("the VT group can be any set of serotypes", {
  sim <- simulate_carriage_data(seed = 1)
  s3 <- c("sero1", "sero2", "sero3")
  dat <- prepare_carriage_data(sim, "survey", s3, "density",
                               vaccine_types = vaccine_types("PCV10", add = "6A"))
  expect_true("6A" %in% dat$sero_levels[dat$vt_mask == 1])
  expect_warning(prepare_carriage_data(sim, "survey", s3, "density",
                                       vaccine_types = "99X"),
                 "None of `vaccine_types`")
})

test_that("any number of covariates can be adjusted for", {
  sim <- simulate_carriage_data(seed = 1)
  cov <- c("age_z", "sex", "income_z", "household", "contact_u5", "chronic_ill")
  dat <- prepare_carriage_data(sim, "survey", c("sero1", "sero2", "sero3"), "density",
                               abundance = c("ra1", "ra2", "ra3"), covariates = cov)
  n_rounds <- length(dat$survey_levels)
  one <- carriage_stan_data(dat, "unified", "equal", covariates = ~ age_z)
  six <- carriage_stan_data(dat, "unified", "equal",
                            covariates = ~ age_z + sex + income_z + household +
                              contact_u5 + chronic_ill)
  expect_equal(one$K, n_rounds + 1)
  expect_equal(six$K, n_rounds + 7)   # household has 3 levels: 2 columns
  inter <- carriage_stan_data(dat, "states", covariates = ~ age_z * sex + I(age_z^2))
  expect_equal(inter$K, n_rounds + 4)
})

test_that("the maximum number of serotypes follows the data or the user", {
  sim5 <- simulate_carriage_data(seed = 1, p_multiple = c(0.6, 0.2, 0.1, 0.06, 0.04))
  s5 <- paste0("sero", 1:5)
  r5 <- paste0("ra", 1:5)
  d5 <- prepare_carriage_data(sim5, "survey", s5, "density", abundance = r5)
  expect_equal(d5$M, 5)
  expect_false(d5$capped)

  d3 <- prepare_carriage_data(sim5, "survey", s5, "density", abundance = r5,
                              max_serotypes = 3)
  expect_equal(d3$M, 3)
  expect_true(d3$capped)
  expect_equal(state_labels(d3$M, d3$capped),
               c("No carriage", "1 serotype", "2 serotypes", "3+ serotypes"))
  st <- describe_carriage_states(d3)
  top <- st$count[st$outcome == "3+ serotypes"]
  expect_equal(sum(top), sum(d3$n_sero >= 3))

  # allocation still uses every detected serotype
  W <- allocate_weights(d3, "equal")
  k5 <- which(d3$m_obs == 5)[1]
  expect_equal(sum(W[k5, ] > 0) >= 1, TRUE)
  expect_equal(max(d3$m_obs), 5)

  sd <- carriage_stan_data(d3, "unified", "equal")
  expect_equal(max(sd$m_obs), 3)
  expect_equal(sd$m_sh_max, 5)
  expect_equal(max(carriage_stan_data(d3, "states")$y), 4)
  expect_equal(max(carriage_stan_data(d3, "unified_multiple", "equal")$m_obs), 3)
  expect_equal(max(carriage_stan_data(d3, "density")$m_sh), 3)

  expect_message(prepare_carriage_data(sim5, "survey", s5, "density",
                                       max_serotypes = 3)$M, NA)
  sim <- simulate_carriage_data(seed = 1)
  expect_message(d <- prepare_carriage_data(sim, "survey", c("sero1", "sero2", "sero3"),
                                            "density", max_serotypes = 5),
                 "using")
  expect_equal(d$M, 3)
  expect_error(prepare_carriage_data(sim5, "survey", s5, "density", max_serotypes = 1),
               "2 to 5")

  sim2 <- simulate_carriage_data(seed = 1, p_multiple = c(0.8, 0.2))
  d2 <- prepare_carriage_data(sim2, "survey", c("sero1", "sero2"), "density",
                              abundance = c("ra1", "ra2"))
  expect_equal(d2$M, 2)
  expect_equal(carriage_stan_data(d2, "states")$G, 3)
})

test_that("models run with up to two and up to five serotypes", {
  skip_if_no_cmdstan()
  sim2 <- simulate_carriage_data(seed = 1, p_multiple = c(0.8, 0.2))
  d2 <- prepare_carriage_data(sim2, "survey", c("sero1", "sero2"), "density",
                              abundance = c("ra1", "ra2"))
  f2 <- do.call(fit_carriage, c(list(d2, approach = "equal"), quick))
  expect_equal(marginal_prevalence(f2, groups = "multiplicity")$outcome |> unique(),
               c("1 serotype", "2 serotypes"))
  expect_equal(density_ratios(f2, multiplicity = TRUE)$outcome |> unique(),
               c("Density", "2 vs 1 serotype"))
  m2 <- do.call(fit_carriage, c(list(d2, approach = "equal", population = "multiple"), quick))
  expect_false("multiplicity" %in% marginal_prevalence(m2)$group)
  expect_true(all(check_fit(f2)$identities$pass))

  sim5 <- simulate_carriage_data(seed = 1, p_multiple = c(0.6, 0.2, 0.1, 0.06, 0.04))
  s5 <- paste0("sero", 1:5)
  r5 <- paste0("ra", 1:5)
  d5 <- prepare_carriage_data(sim5, "survey", s5, "density", abundance = r5)
  f5 <- do.call(fit_carriage, c(list(d5, approach = "density"), quick))
  mult <- marginal_prevalence(f5, groups = "multiplicity")
  expect_setequal(unique(mult$outcome), paste(1:5, c("serotype", rep("serotypes", 4))))
  expect_true(all(check_fit(f5)$identities$pass))
  m5 <- do.call(fit_carriage, c(list(d5, approach = "density", population = "multiple"), quick))
  expect_true("3+ serotypes" %in% marginal_prevalence(m5, groups = "multiplicity")$outcome)

  d3 <- prepare_carriage_data(sim5, "survey", s5, "density", abundance = r5,
                              max_serotypes = 3)
  s3 <- do.call(fit_carriage_states, c(list(d3), quick))
  expect_true("3+ serotypes" %in% marginal_prevalence(s3)$outcome)
  dens <- do.call(fit_density, c(list(d3, group = "total", m_rep = 1000), quick))
  expect_setequal(density_ratios(dens, multiplicity = TRUE)$outcome,
                  c("Density", "2 vs 1 serotype", "3+ vs 1 serotype", "3+ vs 2 serotypes"))
})

test_that("custom vaccine sets can be built from scratch", {
  expect_equal(vaccine_types(NULL, add = c("6B", "14")), c("6B", "14"))
  expect_equal(vaccine_types(NULL, add = c("6B", "14"), remove = "14"), "6B")
  expect_error(vaccine_types(NULL), "add")
})

test_that("crude group proportions are consistent", {
  dat <- example_data()
  nvt <- dat$sero_levels[dat$vt_mask == 0]
  all_nvt <- describe_group_proportions(dat, nvt, among = "NVT", approach = "equal")
  expect_equal(all_nvt$value, rep(1, 3))
  two <- describe_group_proportions(dat, list(a = "6A", b = "19A"), approach = "density")
  expect_equal(unique(two$group), c("a", "b"))
  expect_true(all(two$value >= 0 & two$value <= 1))
  vt <- describe_group_proportions(dat, dat$sero_levels[dat$vt_mask == 1], among = "VT")
  expect_equal(vt$value, rep(1, 3))
  expect_error(describe_group_proportions(dat, c("6B", "6A")), "not in the NVT group")
  expect_error(describe_group_proportions(dat, "6A", approach = "presence"), "Presence-based")
})

test_that("group proportions from a fit", {
  skip_if_no_cmdstan()
  dat <- example_data()
  fits <- do.call(fit_carriage, c(list(dat, approach = c("equal", "density")), quick))
  nvt <- dat$sero_levels[dat$vt_mask == 0]
  gp <- group_proportions(fits$density, nvt)
  expect_equal(gp$median[gp$measure == "proportion among NVT"], rep(1, 3))

  grp <- list("6A+19A" = c("6A", "19A"), "11A" = "11A")
  gp <- group_proportions(fits, grp, reference = "previous")
  expect_setequal(unique(gp$outcome), names(grp))
  expect_setequal(unique(gp$method), c("Equal share", "Serotype-specific density"))
  r <- gp[gp$measure == "proportion ratio", ]
  expect_equal(unique(r$reference), c("2016", "2017"))
  crude <- describe_group_proportions(dat, grp, approach = "equal")
  p <- gp[gp$measure == "proportion among NVT" & gp$method == "Equal share", ]
  expect_equal(p$median, crude$value, tolerance = 0.25)

  pres <- do.call(fit_presence, c(list(dat), quick))
  expect_error(group_proportions(pres, "6A"), "Presence-based")
})
