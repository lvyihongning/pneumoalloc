test_that("allocation weights sum to one in both populations", {
  dat <- example_data()
  for (pop in c("all", "multiple")) {
    n <- if (pop == "all") dat$N_pos else sum(dat$m_obs >= 2)
    for (a in 1:4) {
      W <- allocate_weights(dat, a, pop)
      expect_equal(nrow(W), n)
      expect_equal(unname(rowSums(W)), rep(1, n))
    }
  }
})

test_that("approaches allocate as in Table S1", {
  dat <- example_data()
  W1 <- allocate_weights(dat, "dominant")
  W2 <- allocate_weights(dat, "equal")
  W3 <- allocate_weights(dat, "abundance")
  W4 <- allocate_weights(dat, "density")
  Wp <- allocate_weights(dat, "presence")

  expect_true(all(W1 %in% c(0, 1)))
  expect_equal(W3, W4)
  expect_equal(Wp, (W2 > 0) * 1)
  distinct <- tapply(dat$sero_long$category, dat$sero_long$k, function(v) !anyDuplicated(v))
  k <- which(distinct & dat$m_obs == 3)[1]
  expect_equal(sort(unname(W2[k, W2[k, ] > 0])), rep(1 / 3, 3))
  single <- dat$m_obs == 1
  expect_equal(W1[single, ], W2[single, ])
  expect_equal(W1[single, ], W3[single, ])
})

test_that("approach can be given by number or name", {
  dat <- example_data()
  expect_equal(allocate_weights(dat, 2), allocate_weights(dat, "equal"))
  expect_error(allocate_weights(dat, 5), "must be one of")
  expect_error(allocate_weights(dat, "max"), "must be one of")
})

test_that("crude estimates are consistent", {
  dat <- example_data()
  for (a in c(1:4, 5)) {
    ap <- if (a == 5) "presence" else a
    cr <- describe_prevalence(dat, ap)
    overall <- cr$value[cr$group == "overall"]
    expect_equal(overall, as.vector(tapply(dat$carrier, dat$survey, mean)))
    if (a != 5) {
      sero <- tapply(cr$value[cr$group == "serotype"], cr$survey[cr$group == "serotype"], sum)
      expect_equal(as.vector(sero), overall)
      expect_equal(cr$value[cr$group == "VT"] + cr$value[cr$group == "NVT"], overall)
    }
  }
  cm <- describe_prevalence(dat, "equal", population = "multiple")
  expect_equal(cm$value[cm$group == "overall"], rep(1, 3))

  st <- describe_carriage_states(dat)
  states <- st[st$measure == "prevalence" & st$outcome %in% state_labels(dat$M), ]
  expect_equal(as.vector(tapply(states$value, states$survey, sum)), rep(1, 3))
  ty <- st[st$measure != "prevalence", ]
  expect_equal(as.vector(tapply(ty$value, ty$survey, sum)), rep(1, 3))
})

test_that("describe_prevalence() handles a round without multiple carriers", {
  sim <- simulate_carriage_data(seed = 1)
  last <- sim$survey == levels(factor(sim$survey))[3]
  sim[last, c("sero2", "sero3")] <- NA
  dat <- suppressMessages(prepare_carriage_data(
    sim, "survey", c("sero1", "sero2", "sero3"), density = "density",
    abundance = c("ra1", "ra2", "ra3")
  ))
  lev3 <- dat$survey_levels[3]
  for (ap in c("equal", "density", "presence")) {
    cm <- describe_prevalence(dat, ap, population = "multiple")
    expect_true(all(is.na(cm$value[cm$survey == lev3])))
    expect_false(anyNA(cm$value[cm$survey != lev3]))
  }
  gp <- describe_group_proportions(dat, "15B/C", population = "multiple")
  expect_true(is.na(gp$value[gp$survey == lev3]))
})
