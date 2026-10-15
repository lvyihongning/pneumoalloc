test_that("plots are built from fits", {
  skip_if_no_cmdstan()
  dat <- example_data()
  fits <- do.call(fit_carriage, c(list(dat, approach = c("dominant", "density")), quick))
  pres <- do.call(fit_presence, c(list(dat), quick))
  st <- do.call(fit_carriage_states, c(list(dat), quick))
  ty <- do.call(fit_multiple_types, c(list(dat), quick))
  dens <- do.call(fit_density, c(list(dat, m_rep = 1000), quick))

  built <- function(p) {
    expect_s3_class(p, "ggplot")
    expect_no_error(ggplot2::ggplot_build(p))
  }
  built(plot_prevalence(c(fits, list(presence = pres)), shade_from = "2019"))
  built(plot_ratios(fits))
  built(plot_ratios(fits$density, reference = "previous"))
  built(plot_carriage_states(st))
  built(plot_carriage_states(ty))
  built(plot_density(dens))
  built(plot_density(dens$VT, statistic = "mean", observed = FALSE))
  built(plot_serotypes(fits))
  built(plot_serotype_ratios(fits, survey = "2017"))

  # stacked means add up to overall carriage
  b <- ggplot2::ggplot_build(plot_carriage_states(st, labels = FALSE))$data
  bars <- b[[length(b)]]
  carriage <- marginal_prevalence(st, groups = "carriage")$mean
  expect_equal(as.vector(tapply(bars$y, bars$x, max)), 100 * carriage, tolerance = 1e-6)

  expect_error(plot_carriage_states(fits$density), "fit_carriage_states")
  expect_error(plot_density(fits), "fit_density")
  expect_error(plot_prevalence(fits, shade_from = "1999"), "not a survey round")
})
