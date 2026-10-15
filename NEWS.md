# pneumoalloc 0.1.0

* First version (previously developed as 'seroshare'). Packages the models
  of the PCV10 Nha Trang analysis:
  * `fit_carriage_states()`: overall and multiple serotype carriage
    (`1_Co-colonisation_20261002.stan`, generalised to any number of
    serotypes; same draws);
  * `fit_multiple_types()`: types of multiple carriage
    (`Proportion in co-colonisation_20261002.stan`; same draws);
  * `fit_carriage()`: the unified framework for approaches i-iv among all
    children (`M1to4_unified_v4_20261002.stan`) or among children with
    multiple carriage (`Multi_M1to4_unified_20261002.stan`).
    `approach = "all"` (default) fits all four; one or several can be
    chosen instead. The all-children model is implemented with GLM
    functions and a precomputed composition term: the published model,
    with log densities equal to rounding, about 4.6x faster on a 4-chain
    run on an Apple M4 (2.6 versus 11.9 minutes);
  * `fit_presence()`: the presence-based comparator
    (`Sensi_presence_based_20261002.stan`);
  * `fit_density()`: density among carriers for all pneumococci, PCV10 and
    non-PCV10 serotypes, with an optional survey x multiplicity
    interaction (`All_density_20261002.stan`).
* Approach iv (`density`): carrier weights are Θᵢ / Θ̄ₜ, with Θ̄ₜ the mean
  total density among carriers of the same survey round. The published
  Stan files used one mean over all rounds; normalising within each round
  keeps each round's weights summing to its number of carriers. Approaches
  i-iii are unchanged.
* No commune (cluster) random effects: they were removed from all Stan
  files; the published analyses did not use them.
* `vaccine_types()` gives the serotypes of PCV7, PCV10, PCV10-SII, PCV13,
  PCV15, PCV20, PCV21 and PPSV23, combined and with `add`/`remove`; with
  `pcv10_types()`, `pcv13_types()`, `pcv15_types()` and `pcv20_types()`;
  `vaccine_types(NULL, add = )` builds a set of your own.
  The vaccine-type group is any set of serotypes passed to
  `prepare_carriage_data(vaccine_types = )`.
* `prepare_carriage_data(max_serotypes = )`: any number of serotype
  columns; the models distinguish up to the largest number detected (2 to
  5) or a user-set maximum, with higher numbers pooled (e.g. "3+
  serotypes"). Allocation always uses all detected serotypes.
* Input checks: a serotype listed in two slots of one child is counted
  once (abundances added); `other_label` in `categories` is ignored; with
  the default `categories`, non-vaccine serotypes are in serotype order
  with the most frequently detected one last (the reference category), so
  results do not depend on row order; a warning when the reference
  category holds little allocation; a clear error for a covariate with a
  single value among a model's children.
* Multiplicity: VT and NVT density models distinguish only the numbers of
  group serotypes observed; the unified model reports density for 1, 2
  and 3+ serotypes, the levels of its density component; among multiple
  carriers, the 3-serotype term is dropped when no child has three or
  more serotypes (or `max_serotypes = 2`).
* `describe_prevalence()` and `describe_group_proportions()` handle survey
  rounds without (multiple) carriers.
* Adjusted models take any number of covariates (numeric or factor,
  transformations and interactions) through the `covariates` formula.
* `prepare_carriage_data()`, `allocate_weights()` (including
  presence-based indicators), `describe_prevalence()`,
  `describe_carriage_states()` and `describe_density()` for data and crude
  estimates.
* `marginal_prevalence()`, `prevalence_ratios()` (relative to the
  pre-vaccination survey, the previous survey or any survey),
  `marginal_density()` and `density_ratios()` (mean and/or median density,
  `statistic`) for results; `check_fit()` and `compare_approaches()` for
  checks.
* `group_proportions()` gives the proportion of a group of serotypes
  within VT, NVT or all carriage by survey round, and its ratio between
  rounds; `describe_group_proportions()` the crude counterpart.
* `carriage_priors()` sets the prior scales of all models.
* `table_s1`: the worked example of Table S1.
* `inst/scripts/quick_test.R` runs every model on simulated data.
