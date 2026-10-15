test_that("prepare_carriage_data() defines carriers, states and categories", {
  sim <- simulate_carriage_data(seed = 2)
  dat <- example_data(seed = 2)

  expect_s3_class(dat, "carriage_data")
  expect_equal(dat$N, nrow(sim))
  expect_equal(dat$N_pos, sum(!is.na(sim$sero1)))
  expect_equal(dat$n_sero, as.integer(rowSums(!is.na(sim[, c("sero1", "sero2", "sero3")]))))
  expect_equal(dat$sero_levels[length(dat$sero_levels)], "Other NVT")
  expect_true(all(dat$sero_levels[dat$vt_mask == 1] %in% pcv10_types()))
  expect_equal(dat$reference_level, "2016")
  expect_equal(as.vector(tapply(dat$sero_long$abundance, dat$sero_long$k, sum)),
               rep(1, dat$N_pos))
})

test_that("types of multiple carriage follow the detected and dominant serotypes", {
  dat <- example_data()
  multi <- dat$m_obs >= 2
  expect_true(all(is.na(dat$multiple_type[!multi])))
  expect_true(all(dat$multiple_type[multi] %in% type_labels()))

  sl <- dat$sero_long
  for (k in which(multi)[1:20]) {
    vt <- sl$vaccine_type[sl$k == k]
    dom_vt <- sl$vaccine_type[sl$k == k & sl$dominant]
    expected <- if (all(vt)) "VT-only" else if (!any(vt)) "NVT-only" else
      if (dom_vt) "Mixed, VT-dominant" else "Mixed, NVT-dominant"
    expect_equal(dat$multiple_type[k], expected)
  }
})

test_that("the reference survey can be chosen", {
  sim <- simulate_carriage_data(seed = 1)
  dat <- prepare_carriage_data(sim, "survey", c("sero1", "sero2", "sero3"),
                               density = "density", reference = "2017")
  expect_equal(dat$reference, 2L)
  expect_error(prepare_carriage_data(sim, "survey", c("sero1", "sero2", "sero3"),
                                     density = "density", reference = "1999"),
               "not a survey round")
})

test_that("children with missing covariates are excluded", {
  sim <- simulate_carriage_data(seed = 1)
  sim$age_z[1:5] <- NA
  expect_message(
    dat <- prepare_carriage_data(sim, "survey", c("sero1", "sero2", "sero3"),
                                 density = "density", covariates = "age_z"),
    "Excluded 5"
  )
  expect_equal(dat$N, nrow(sim) - 5)
})

test_that("dominant serotype follows abundance or slot order", {
  sim <- simulate_carriage_data(seed = 3)
  i <- which(!is.na(sim$sero2))[1]
  sim[i, c("sero1", "sero2", "ra1", "ra2")] <- sim[i, c("sero2", "sero1", "ra2", "ra1")]
  args <- list(sim, "survey", c("sero1", "sero2", "sero3"), density = "density",
               abundance = c("ra1", "ra2", "ra3"))
  dom <- function(dat) dat$sero_long$serotype[dat$sero_long$child == i & dat$sero_long$dominant]
  expect_equal(dom(do.call(prepare_carriage_data, args)), sim$sero2[i])
  expect_equal(dom(do.call(prepare_carriage_data, c(args, dominant = "first"))), sim$sero1[i])
})

test_that("prepare_carriage_data() gives informative errors", {
  sim <- simulate_carriage_data(seed = 1)
  s3 <- c("sero1", "sero2", "sero3")
  expect_error(prepare_carriage_data(sim, "survey", c("sero1", "nope"), "density"), "not found")
  expect_error(prepare_carriage_data(sim, "survey", s3, "density", categories = "6A"),
               "vaccine serotypes")
  expect_error(prepare_carriage_data(sim, "survey", "sero1", "density"), "more than one serotype")
  expect_error(prepare_carriage_data(sim, "survey", s3, "density", dominant = "abundance"), "needs")
  sim$density[which(!is.na(sim$sero1))[1]] <- NA
  expect_error(prepare_carriage_data(sim, "survey", s3, "density"), "positive and finite")
})
