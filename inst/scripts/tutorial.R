############################################################
# pneumoalloc tutorial
#
# A step-by-step introduction to the package, on simulated data.
# Open this file in RStudio and run it section by section
# (Ctrl/Cmd + Enter runs the current line or selection).
#
# Find this file after installing the package with
#   file.edit(system.file("scripts", "tutorial.R", package = "pneumoalloc"))
#
# The model fits below use SHORT chains (2 chains x 500 draws) so that
# each takes seconds. For real analyses, leave the sampler arguments out:
# the defaults are 4 chains x (1000 warm-up + 2000 sampling) iterations.
############################################################


# ==========================================================
# 0) Setup
# ==========================================================

# Once per computer (not run here):
#   install.packages("cmdstanr",
#     repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
#   cmdstanr::install_cmdstan()
#   remotes::install_local("~/Documents/pneumoalloc")

library(pneumoalloc)

# Where to write example output files (a temporary folder here; use your
# own results folder in a real analysis)
out_dir <- tempdir()


# ==========================================================
# 1) The idea: allocating a carriage episode among its serotypes
# ==========================================================

# Table S1: ten hypothetical children, six carriers. Child 6 carries
# 6B (relative abundance 0.5), 15A (0.3) and 23A (0.2).
table_s1

# prepare_carriage_data() turns a data frame into the object every other
# function uses. You tell it which columns hold what.
ex <- prepare_carriage_data(
  table_s1,
  survey    = "survey",                       # survey round
  serotypes = c("sero1", "sero2", "sero3"),   # detected serotypes, one column per slot
  abundance = c("ra1", "ra2", "ra3"),         # relative abundance of each slot
  density   = "density"                       # total density (GE/mL)
)

# How much each carrier contributes to each serotype, under each approach.
# Rows are carriers, columns are serotypes; rows sum to 1 for the four
# allocation approaches.
allocation_approaches()
allocate_weights(ex, "dominant")
allocate_weights(ex, "equal")
allocate_weights(ex, "abundance")
allocate_weights(ex, "presence")   # every detected serotype counts fully

# Crude (model-free) prevalence under one approach. `value` is a proportion.
describe_prevalence(ex, "dominant")
describe_prevalence(ex, "presence")
# PCV10 + non-PCV10 = overall carriage (0.6) for the allocation approaches,
# but more than that for the presence-based approach.


# ==========================================================
# 2) Your survey data
# ==========================================================

# The data need ONE ROW PER CHILD with:
#   - the survey round
#   - the detected serotypes, one column per slot (NA or "" if none);
#     any number of slots: 2, 3, 5, ...
#   - the relative abundance of each slot (for approaches 3 and 4 and
#     group density)
#   - total density (for carriers)
#   - covariates for adjusted models, already prepared (scaled, coded)
# simulate_carriage_data() makes such a data set: three surveys, PCV10
# serotypes declining after the first (pre-vaccination) survey, and six
# covariates (age_z, sex, income_z, household, contact_u5, chronic_ill).
sim <- simulate_carriage_data(n_per_survey = 400, seed = 2026)
head(sim)

dat <- prepare_carriage_data(
  sim,
  survey        = "survey",
  serotypes     = c("sero1", "sero2", "sero3"),
  abundance     = c("ra1", "ra2", "ra3"),
  density       = "density",
  vaccine_types = pcv10_types(),
  # serotypes modelled as their own category; all other non-PCV10
  # serotypes are pooled into "Other NVT" (the reference category)
  categories    = c(pcv10_types(), "6A", "15B/C", "19A", "11A", "35B", "23A"),
  # EVERY covariate you may adjust for (as many as you like): children with
  # missing values are excluded once, so unadjusted and adjusted models use
  # the same children
  covariates    = c("age_z", "sex", "income_z", "household",
                    "contact_u5", "chronic_ill"),
  reference     = "2016"          # pre-vaccination survey: reference for ratios
)
dat

# With your own data, replace the column names by those in your data
# (not run):
#   dat <- prepare_carriage_data(
#     my_data,
#     survey = "<survey round column>",
#     serotypes = c("<serotype 1>", "<serotype 2>", "<serotype 3>"),
#     abundance = c("<abundance 1>", "<abundance 2>", "<abundance 3>"),
#     density = "<total density column>",
#     vaccine_types = pcv10_types(),
#     categories = c(pcv10_types(), "<non-vaccine serotypes of interest>"),
#     covariates = c("<covariate 1>", "<covariate 2>", "..."),
#     reference = "<pre-vaccination survey round>"
#   )

# Crude estimates, before any model
describe_carriage_states(dat)   # no carriage / 1 / 2 / 3 serotypes; types of multiple carriage
describe_prevalence(dat, "density")
describe_density(dat)           # observed density: total, PCV10, non-PCV10 carriers


# ----------------------------------------------------------
# 2a) Covariates: choose them in the model formula
# ----------------------------------------------------------
# The adjusted models take a one-sided formula with any of the declared
# covariates, numeric or factor (factors get one column per level).
# Survey round is always in the model. We use all six here:
adj <- ~ age_z + sex + income_z + household + contact_u5 + chronic_ill

# Other formulas work as usual, for example
#   ~ age_z                               one covariate
#   ~ age_z + I(age_z^2) + sex            a quadratic age effect
#   ~ age_z * sex                         an interaction
# Leave `covariates` out for the unadjusted model.


# ----------------------------------------------------------
# 2b) Defining the vaccine-type (VT) group
# ----------------------------------------------------------
# `vaccine_types` is any set of serotypes: they form the VT group, and
# every other detected serotype the NVT group. Ready-made sets:
vaccine_types("PCV13")
setdiff(pcv20_types(), pcv13_types())          # serotypes added by PCV20
# Your own group, e.g. 6A as vaccine-associated (cross-protection):
vaccine_types("PCV10", add = "6A")
# Spell serotypes as in your data (e.g. "15B/C" from microarray):
vaccine_types("PCV20", add = "15B/C", remove = "15B")
# Or a set entirely of your own: any character vector works
vaccine_types(NULL, add = c("6B", "14", "19F", "23F"))
my_vt <- c("6B", "14", "19F", "23F", "6A")
# To use one, pass it to prepare_carriage_data(vaccine_types = ...), e.g.
#   dat_6A <- prepare_carriage_data(sim, ..., vaccine_types = vaccine_types("PCV10", add = "6A"))


# ----------------------------------------------------------
# 2c) Number of detected serotypes per child
# ----------------------------------------------------------
# By default the models distinguish as many serotypes per child as the
# data contain (here up to 3). Data with up to 5 serotypes:
sim5 <- simulate_carriage_data(seed = 1, p_multiple = c(0.6, 0.2, 0.1, 0.06, 0.04))
dat5 <- prepare_carriage_data(
  sim5, survey = "survey",
  serotypes = paste0("sero", 1:5), abundance = paste0("ra", 1:5),
  density = "density",
  max_serotypes = 3   # count 4 and 5 serotypes as "3+"; leave out to keep 1-5
)
dat5
describe_carriage_states(dat5)
# The allocation approaches always use all detected serotypes: a child with
# five serotypes contributes 1/5 to each of them under `equal`.


# ==========================================================
# 3) Step 1: overall and multiple serotype carriage
# ==========================================================

# Multinomial-logit model for the number of detected serotypes.
# covariates = NULL (the default) gives the UNADJUSTED model;
# a one-sided formula gives the ADJUSTED model.
# The first fit compiles the Stan model (a minute or two, once).
states <- fit_carriage_states(
  dat,
  covariates = adj,
  chains = 2, iter_warmup = 500, iter_sampling = 500, refresh = 0
)
states

# Standardised prevalence by survey (numbers are proportions; `estimate`
# is "median (2.5%, 97.5%)" in percent)
marginal_prevalence(states)

# Prevalence ratios. Default: each survey versus the pre-vaccination survey.
prevalence_ratios(states, groups = c("carriage", "multiple"))
# ... or versus the previous survey, or versus any survey you name:
prevalence_ratios(states, groups = c("carriage", "multiple"), reference = "previous")
prevalence_ratios(states, groups = c("carriage", "multiple"), reference = "2017")

# Among children with multiple carriage: VT-only, NVT-only, mixed.
types <- fit_multiple_types(
  dat, covariates = adj,
  chains = 2, iter_warmup = 500, iter_sampling = 500, refresh = 0
)
marginal_prevalence(types)
prevalence_ratios(types)


# ==========================================================
# 4) Step 2: serotype-specific carriage (unified framework)
# ==========================================================

# Choose the approach(es) by the aim of the analysis:
#   approach = "density"                    one approach
#   approach = c("dominant", "density")     several
#   approach = "all"                        all four (the default)

# One approach -> one fit
fit_dens <- fit_carriage(
  dat, approach = "density", covariates = adj,
  chains = 2, iter_warmup = 500, iter_sampling = 500, refresh = 0
)
fit_dens

# The results are tibbles: filter them like any data frame.
prev <- marginal_prevalence(fit_dens, groups = c("overall", "VT", "NVT"))
prev[prev$measure == "prevalence", c("outcome", "survey", "estimate")]

pr <- prevalence_ratios(fit_dens, groups = c("VT", "NVT"))
pr[pr$measure == "prevalence", c("outcome", "survey", "reference", "estimate")]

# Serotype-specific prevalence ratios. `reportable` is FALSE when the
# reference-survey prevalence is near zero (the ratio is then unstable).
pr_sero <- prevalence_ratios(fit_dens, groups = "serotype")
pr_sero[pr_sero$measure == "prevalence",
        c("outcome", "vaccine_type", "survey", "estimate", "reportable")]

# All four approaches -> a list of fits, named by approach
fits <- fit_carriage(
  dat, approach = "all", covariates = adj,
  chains = 2, iter_warmup = 500, iter_sampling = 500, refresh = 0
)
fits
names(fits)
fits$dominant


# Proportion of a serotype group within non-PCV10 (NVT) carriage, and its
# ratio versus the pre-vaccination survey, e.g. 15B/C + 23A + 35B. The
# serotypes must be categories of the fit: give them their own category in
# prepare_carriage_data(categories = ...), otherwise they are pooled into
# "Other NVT".
group_proportions(fit_dens, c("15B/C", "23A", "35B"), among = "NVT")
# several groups at once, all four approaches, versus the previous survey:
group_proportions(fits, list("15B/C+23A+35B" = c("15B/C", "23A", "35B"),
                             "19A" = "19A"),
                  among = "NVT", reference = "previous")
# among PCV10 (VT) carriage instead: among = "VT"; crude version:
describe_group_proportions(dat, c("15B/C", "23A", "35B"), approach = "density")


# ==========================================================
# 5) Comparing approaches, including the presence-based one
# ==========================================================

presence <- fit_presence(
  dat, covariates = adj,
  chains = 2, iter_warmup = 500, iter_sampling = 500, refresh = 0
)

cmp <- compare_approaches(c(fits, list(presence = presence)))

# (a) Overall carriage, multiplicity and density do not depend on the
#     allocation, so they must agree across the four approaches: |z| should
#     be small (Monte Carlo error only).
summary(cmp$shared$z)

# (b) PCV10 and non-PCV10 prevalence ratios of all approaches side by side
r <- cmp$ratios
r[r$group %in% c("VT", "NVT") & r$measure == "prevalence",
  c("method", "outcome", "survey", "estimate")]


# ==========================================================
# 6) Step 3: serotype composition among children with multiple carriage
# ==========================================================

multi <- fit_carriage(
  dat, approach = "density", population = "multiple",
  covariates = adj,
  chains = 2, iter_warmup = 500, iter_sampling = 500, refresh = 0
)
# Here the estimates are proportions among multiple carriers, and the
# ratios are proportion ratios.
prevalence_ratios(multi, groups = c("VT", "NVT"))


# ==========================================================
# 7) Step 4: density among carriers
# ==========================================================

# Total, PCV10 and non-PCV10 density (one model each)
dens <- fit_density(
  dat, covariates = adj, m_rep = 1000,
  chains = 2, iter_warmup = 500, iter_sampling = 500, refresh = 0
)
dens

# Mean and median density by survey; the column `statistic` says which
marginal_density(dens$total)

# Median density ratios between surveys, and density ratios by number of
# serotypes (2 vs 1, 3 vs 1, 3 vs 2)
density_ratios(dens, statistic = "median", multiplicity = TRUE)

# Does the multiplicity effect change between surveys? (interaction model)
dens_int <- fit_density(
  dat, group = "total", interaction = TRUE, covariates = adj,
  m_rep = 1000, chains = 2, iter_warmup = 500, iter_sampling = 500, refresh = 0
)
density_ratios(dens_int, statistic = "median", multiplicity = TRUE)


# ==========================================================
# 7b) Figures
# ==========================================================
# Each function returns a ggplot object: print it, change it with the usual
# ggplot2 syntax, combine panels (e.g. with 'patchwork') or save with
# ggplot2::ggsave().
plot_carriage_states(states)                          # Fig 1A: carriage states
plot_carriage_states(types)                           # Fig 1B: types of multiple carriage
plot_prevalence(c(fits, list(presence = presence)))   # Fig 2A-B: VT and NVT prevalence
plot_ratios(fits)                                     # Fig 2C-D: prevalence ratios
plot_density(dens)                                    # Fig 3: density
plot_serotypes(fits[c("dominant", "density")])        # Fig 4A: serotypes, 2016 vs 2019
plot_serotype_ratios(fits[c("dominant", "density")])  # Fig 4B: serotype ratios
# ggplot2::ggsave("fig2.png", plot_prevalence(fits), width = 9, height = 4.5, dpi = 300)


# ==========================================================
# 8) Checking a fit
# ==========================================================

# Sampler diagnostics, R-hat and ESS, identities that must hold in every
# draw, and the comparison with the crude estimates.
#   v = passed, x = failed, i = information
# With the short chains used here, R-hat can be slightly above 1.01;
# with the default chains it should not be.
check_fit(fit_dens)
check_fit(states)
check_fit(dens$total)

# The details are in the returned list:
chk <- check_fit(fit_dens)
names(chk)
chk$crude


# ==========================================================
# 9) Priors
# ==========================================================

# Defaults: those of the PCV10 Nha Trang analysis (Table S2)
carriage_priors()

# A sensitivity analysis with a wider prior on serotype composition
# (passed as data: no recompilation)
fit_wide <- fit_carriage(
  dat, approach = "density", covariates = adj,
  priors = carriage_priors(share = 5),
  chains = 2, iter_warmup = 500, iter_sampling = 500, refresh = 0
)
prevalence_ratios(fit_wide, groups = c("VT", "NVT"))


# ==========================================================
# 10) Saving and exporting
# ==========================================================

# Save fits with save_fit(), not saveRDS(): it loads the draws first,
# because they otherwise live in temporary files.
save_fit(fits, file.path(out_dir, "fits_adjusted.rds"))
fits_again <- readRDS(file.path(out_dir, "fits_adjusted.rds"))
prevalence_ratios(fits_again$density, groups = c("VT", "NVT"))

# All tables of a fit (or of a list of fits) in one Excel workbook,
# one sheet per table
write_results(fits, file.path(out_dir, "results_adjusted.xlsx"))

# Or collect them in R
s <- summary(fits)
names(s)


# ==========================================================
# 11) Practical tips
# ==========================================================

# * Unadjusted model: leave out `covariates`. Adjusted model: give a
#   formula with variables declared in prepare_carriage_data(covariates =).
# * Sampler settings: chains, iter_warmup, iter_sampling, seed,
#   adapt_delta, max_treedepth. Keep threads_per_chain = 1 and give each
#   chain its own core (4 chains: 4 cores).
# * On a cluster, compile each model once into a shared folder and pass it:
#     mod <- carriage_model("unified", dir = "~/stan_models")
#     fit_carriage(dat, approach = "dominant", model = mod)
# * The Stan programs: carriage_stan_file("unified"), ("unified_multiple"),
#   ("categorical"), ("presence"), ("density").
# * Help for any function: ?fit_carriage, ?prevalence_ratios, ...
#   Overview: vignette("pneumoalloc").

cat("\nTutorial finished.\n")
