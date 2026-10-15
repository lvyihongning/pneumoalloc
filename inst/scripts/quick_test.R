# pneumoalloc quick test
#
# Runs every model of the package on a small simulated data set, with short
# chains, and prints the main results and checks. Run it with
#
#   source(system.file("scripts", "quick_test.R", package = "pneumoalloc"))
#
# It takes a few minutes, most of which is compiling the five Stan models
# the first time (they are cached afterwards). Short chains are for
# testing only: use the defaults (4 chains x 2000 draws) for analyses.
# With these short chains, an R-hat slightly above 1.01 (marked x in the
# checks) is expected and not a problem with the model.

library(pneumoalloc)

# ---- 1. Simulated surveys: 3 rounds, PCV10 serotypes declining ----------
sim <- simulate_carriage_data(n_per_survey = 400, seed = 2026)

dat <- prepare_carriage_data(
  sim,
  survey = "survey",
  serotypes = c("sero1", "sero2", "sero3"),
  abundance = c("ra1", "ra2", "ra3"),
  density = "density",
  vaccine_types = pcv10_types(),
  categories = c(pcv10_types(), "6A", "15B/C", "19A", "11A", "35B"),
  covariates = c("age_z", "sex"),
  reference = "2016"                       # pre-vaccination survey
)
print(dat)

short <- list(chains = 2, iter_warmup = 500, iter_sampling = 500, refresh = 0)
adj <- ~ age_z + sex

# ---- 2. Overall and multiple serotype carriage ---------------------------
states <- do.call(fit_carriage_states, c(list(dat, covariates = adj), short))
print(prevalence_ratios(states, groups = c("carriage", "multiple"))[
  , c("outcome", "survey", "estimate")])

types <- do.call(fit_multiple_types, c(list(dat, covariates = adj), short))
print(marginal_prevalence(types)[, c("outcome", "survey", "estimate")])

# ---- 3. Serotype-specific carriage: four approaches + presence-based ------
fits <- do.call(fit_carriage, c(list(dat, approach = "all", covariates = adj), short))
presence <- do.call(fit_presence, c(list(dat, covariates = adj), short))

cmp <- compare_approaches(c(fits, list(presence = presence)))
vt_nvt <- cmp$ratios[cmp$ratios$group %in% c("VT", "NVT") &
                       cmp$ratios$measure == "prevalence", ]
cat("\nPCV10 (VT) and non-PCV10 (NVT) prevalence ratios versus 2016:\n")
print(as.data.frame(vt_nvt[, c("method", "outcome", "survey", "estimate")]),
      row.names = FALSE)

# relative to the previous survey instead of the pre-vaccination survey
print(prevalence_ratios(fits$density, groups = c("VT", "NVT"),
                        reference = "previous")[
  prevalence_ratios(fits$density, groups = c("VT", "NVT"),
                    reference = "previous")$measure == "prevalence",
  c("outcome", "survey", "reference", "estimate")])

# ---- 4. Serotype composition among children with multiple carriage -------
multi <- do.call(fit_carriage, c(list(dat, approach = "density",
                                      population = "multiple", covariates = adj), short))
print(prevalence_ratios(multi, groups = c("VT", "NVT"))[
  , c("measure", "outcome", "survey", "estimate")])

# ---- 5. Density: total, PCV10 and non-PCV10 --------------------------------
dens <- do.call(fit_density, c(list(dat, covariates = adj, m_rep = 1000), short))
print(density_ratios(dens, statistic = "median", multiplicity = TRUE)[
  , c("method", "outcome", "survey", "estimate")])

# ---- 6. Checks ---------------------------------------------------------------
# With short chains some R-hat values can exceed 1.01; the identities must
# hold (deviation below 1e-6) and there should be no divergent transitions.
print(check_fit(fits$density))
print(check_fit(states))
print(check_fit(dens$total))

cat("\nQuick test finished.\n")
