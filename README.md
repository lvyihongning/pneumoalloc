# pneumoalloc

<!-- badges: start -->
<!-- badges: end -->

**pneumoalloc** estimates overall, multiple serotype and serotype-specific
pneumococcal carriage prevalence, and their prevalence ratios across
repeated cross-sectional surveys, when children can carry more than one
serotype at a time.

| Step | Model | Function |
|---|---|---|
| 1 | Carriage states: no carriage, 1, 2, 3, ... serotypes (multinomial logit) | `fit_carriage_states()` |
| 1 | Types of multiple carriage: VT-only, NVT-only, mixed VT-/NVT-dominant (multinomial logit) | `fit_multiple_types()` |
| 2 | Unified framework among all children: carriage (logistic) + serotype composition among carriers (fractional multinomial logit) + density (lognormal with multiplicity) | `fit_carriage(population = "all")` |
| 3 | Unified framework among children with multiple carriage | `fit_carriage(population = "multiple")` |
| 4 | Density among carriers: total, PCV10 and non-PCV10 (lognormal with multiplicity; optional survey x multiplicity interaction) | `fit_density()` |
| comparator | Presence-based carriage (independent logistic models) | `fit_presence()` |

In the unified framework, four approaches to allocating a carriage episode
among its serotypes differ **only** in the allocation weights, so they are
compared on the same population, design and priors:

| | approach | allocation of a carriage episode | carrier weight |
|---|---|---|---|
| i | `dominant` | all to the dominant serotype | 1 |
| ii | `equal` | 1/n to each of the n serotypes | 1 |
| iii | `abundance` | relative abundance of each serotype | 1 |
| iv | `density` | relative abundance of each serotype | density / mean density |

All estimates are standardised over the study population. Prevalence
ratios are formed within each posterior draw, relative to the
pre-vaccination survey (default), the previous survey, or any survey round.
Default priors are those of the PCV10 Nha Trang analysis and can be changed
with `carriage_priors()`.

The package was developed for the analysis of PCV10 introduction in Nha
Trang, Viet Nam (Dong et al., in preparation). Given the same data, it
builds exactly the model inputs of that analysis.

## Flexible inputs

- **Vaccine group:** `vaccine_types` takes any set of serotypes (your own,
  e.g. `c("6B", "14", "19F", "23F")`, or `vaccine_types(NULL, add = ...)`); ready-made
  sets for PCV7, PCV10, PCV10-SII, PCV13, PCV15, PCV20, PCV21 and PPSV23,
  e.g. `vaccine_types("PCV13")` or `vaccine_types("PCV10", add = "6A")`.
- **Covariates:** any number, numeric or factor, in the `covariates`
  formula, e.g. `~ age_z + sex + income_z + household`.
- **Serotypes per child:** any number of serotype columns; the models
  distinguish up to 5 per child, or fewer with `max_serotypes` (higher
  numbers pooled, e.g. "3+ serotypes").

## Data

The package contains no study data. Examples, tests and the vignette use
data simulated with `simulate_carriage_data()` and the hypothetical
ten-child example `table_s1`.

## Installation

pneumoalloc fits its models with [CmdStan](https://mc-stan.org/cmdstanr/)
through 'cmdstanr', which is not on CRAN:

``` r
install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
cmdstanr::install_cmdstan()

# install.packages("remotes")
remotes::install_local("path/to/pneumoalloc")   # or install_github() once published
```

## Tutorial

A commented, step-by-step script that teaches the whole workflow on
simulated data; open it in RStudio and run it section by section:

``` r
file.edit(system.file("scripts", "tutorial.R", package = "pneumoalloc"))
```

## Quick test

After installing, run every model on simulated data with short chains
(a few minutes, most of it compiling the Stan models the first time):

``` r
source(system.file("scripts", "quick_test.R", package = "pneumoalloc"))
```

## Example

``` r
library(pneumoalloc)

sim <- simulate_carriage_data(seed = 1)
dat <- prepare_carriage_data(
  sim,
  survey = "survey",
  serotypes = c("sero1", "sero2", "sero3"),
  abundance = c("ra1", "ra2", "ra3"),
  density = "density",
  vaccine_types = pcv10_types(),
  categories = c(pcv10_types(), "6A", "15B/C", "19A", "11A", "35B"),
  covariates = c("age_z", "sex"),
  reference = "2016"                    # pre-vaccination survey
)

# step 1: overall and multiple serotype carriage
states <- fit_carriage_states(dat, covariates = ~ age_z + sex)
prevalence_ratios(states, groups = c("carriage", "multiple"))

# step 2: serotype-specific carriage, four approaches
fits <- fit_carriage(dat, approach = "all", covariates = ~ age_z + sex)
prevalence_ratios(fits$density, groups = c("VT", "NVT"))
prevalence_ratios(fits$density, groups = c("VT", "NVT"), reference = "previous")
# share of 15B/C + 15A + 23A in non-PCV10 carriage, and its ratio over time
group_proportions(fits$density, c("15B/C", "15A", "23A"), among = "NVT")
check_fit(fits$density)
compare_approaches(fits)

# comparator: presence-based approach
presence <- fit_presence(dat, covariates = ~ age_z + sex)
compare_approaches(c(fits, list(presence = presence)))

# step 3: serotype composition of multiple carriage
multi <- fit_carriage(dat, approach = "all", population = "multiple",
                      covariates = ~ age_z + sex)

# step 4: density among carriers (total, PCV10, non-PCV10)
dens <- fit_density(dat, covariates = ~ age_z + sex)
density_ratios(dens, statistic = "median", multiplicity = TRUE)

# figures (ggplot objects)
plot_prevalence(fits)
plot_ratios(fits)
plot_serotype_ratios(fits[c("dominant", "density")])

save_fit(fits, "fits_adjusted.rds")
write_results(fits, "results_adjusted.xlsx")
```

See `vignette("pneumoalloc")` for a full walk-through, including the
worked example of Table S1 (`table_s1`).

## Development

The package follows the standard workflow of
[R Packages](https://r-pkgs.org) (usethis, devtools, roxygen2, testthat):

``` r
devtools::load_all()   # load the package
devtools::document()   # rebuild man/ and NAMESPACE from roxygen comments
devtools::test()       # run the tests (model fits need CmdStan)
devtools::check()      # R CMD check
```

The Stan models are in `inst/stan/`.
