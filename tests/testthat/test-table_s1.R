# Totals of Table S1 (ten children): PCV10, non-PCV10 and overall carriage
# under each approach, as numbers of children.
test_that("describe_prevalence() reproduces Table S1", {
  dat <- prepare_carriage_data(table_s1, "survey", c("sero1", "sero2", "sero3"),
                               density = "density", abundance = c("ra1", "ra2", "ra3"))
  expected <- list(
    dominant  = c(VT = 5,     NVT = 1,     Overall = 6),
    equal     = c(VT = 10/3,  NVT = 8/3,   Overall = 6),
    abundance = c(VT = 3.8,   NVT = 2.2,   Overall = 6),
    density   = c(VT = 3.95,  NVT = 2.05,  Overall = 6),
    presence  = c(VT = 5,     NVT = 4,     Overall = 6)
  )
  for (a in names(expected)) {
    cr <- describe_prevalence(dat, a)
    got <- 10 * cr$value[match(c("VT", "NVT", "Overall"), cr$outcome)]
    expect_equal(got, unname(expected[[a]]), tolerance = 1e-12, info = a)
  }
})

test_that("allocation weights reproduce the rows of Table S1", {
  dat <- prepare_carriage_data(table_s1, "survey", c("sero1", "sero2", "sero3"),
                               density = "density", abundance = c("ra1", "ra2", "ra3"))
  # child 6: 6B (0.5), 15A (0.3), 23A (0.2), density 0.5 x 10^5
  k6 <- which(dat$pos_idx == 6)
  w <- function(a) allocate_weights(dat, a)[k6, c("6B", "15A", "23A")]
  expect_equal(unname(w("dominant")), c(1, 0, 0))
  expect_equal(unname(w("equal")), rep(1 / 3, 3))
  expect_equal(unname(w("abundance")), c(0.5, 0.3, 0.2))
  # density approach: theta * Theta / mean(Theta); mean density is 10^5
  expect_equal(unname(w("density")) * dat$dens_pos[k6] / mean(dat$dens_pos),
               c(0.25, 0.15, 0.10))
  expect_equal(unname(w("presence")), c(1, 1, 1))
})

test_that("models need two survey rounds", {
  dat <- prepare_carriage_data(table_s1, "survey", c("sero1", "sero2", "sero3"),
                               density = "density")
  expect_error(carriage_stan_data(dat, "states"), "two survey rounds")
})

test_that("group densities and presence outcomes follow Table S1", {
  dat <- prepare_carriage_data(table_s1, "survey", c("sero1", "sero2", "sero3"),
                               density = "density", abundance = c("ra1", "ra2", "ra3"))
  # serotype-specific density = relative abundance x total density
  vt <- density_group(dat, "VT")
  expect_equal(vt$rows, c(1, 2, 4, 5, 6))
  expect_equal(vt$density, c(1, 1, 1.4, 0.3, 0.25) * 1e5)
  expect_equal(vt$m, c(1, 1, 1, 1, 1))
  nvt <- density_group(dat, "NVT")
  expect_equal(nvt$rows, c(3, 4, 5, 6))
  expect_equal(nvt$density, c(1, 0.6, 0.2, 0.25) * 1e5)
  expect_equal(nvt$m, c(1, 1, 1, 2))
  expect_equal(density_group(dat, "NVT", multiplicity = "total")$m, c(1, 2, 2, 3))

  y <- presence_outcomes(dat)
  expect_equal(colSums(y)[c("Overall", "VT", "NVT", "Multiple carriage")],
               c(Overall = 6, VT = 5, NVT = 4, `Multiple carriage` = 3))
})
