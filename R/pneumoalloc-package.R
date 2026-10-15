#' pneumoalloc: serotype-specific pneumococcal carriage with multiple
#' serotype carriage
#'
#' Bayesian models for repeated cross-sectional carriage surveys in which
#' a child can carry more than one serotype. Every model standardises its
#' estimates over the observed covariate distribution of the study
#' population (posterior g-computation), and prevalence ratios are formed
#' within each posterior draw.
#'
#' **Step 1: overall and multiple serotype carriage**
#'
#' * [fit_carriage_states()]: multinomial logit for the number of
#'   detected serotypes (none, one, two, three, ...), giving overall
#'   carriage, multiple serotype carriage and each carriage state;
#' * [fit_multiple_types()]: among children with multiple carriage,
#'   multinomial logit for VT-only, NVT-only, mixed VT-dominant and mixed
#'   NVT-dominant carriage.
#'
#' **Step 2: serotype-specific carriage (unified framework)**
#'
#' [fit_carriage()] fits one joint model of overall carriage (logistic),
#' serotype composition among carriers (fractional multinomial logit) and
#' density among carriers (lognormal with a multiplicity sub-model), in
#' all children (`population = "all"`) or among children with multiple
#' carriage (`population = "multiple"`). The four allocation approaches
#' differ **only** in the allocation weights of the composition
#' component:
#'
#' | approach | name        | allocation of a carriage episode      | carrier weight |
#' |----------|-------------|---------------------------------------|----------------|
#' | i        | `dominant`  | all to the dominant serotype          | 1              |
#' | ii       | `equal`     | 1 / n to each of the n serotypes      | 1              |
#' | iii      | `abundance` | relative abundance of each serotype   | 1              |
#' | iv       | `density`   | relative abundance of each serotype   | \eqn{\Theta_i / \bar{\Theta}_t} |
#'
#' Here \eqn{\Theta_i} is the total pneumococcal density of child i,
#' \eqn{\theta_{is}} the relative abundance of serotype s (so the
#' serotype-specific density is \eqn{\theta_{is} \Theta_i}), and
#' \eqn{\bar{\Theta}_t} the mean total density among the carriers of the
#' same survey round t. Normalising within each round keeps the weights of a
#' round summing to its number of carriers.
#'
#' **Results**: [marginal_prevalence()], [prevalence_ratios()],
#' [marginal_density()], [density_ratios()], [check_fit()] and
#' [compare_approaches()]. [describe_prevalence()] and
#' [describe_carriage_states()] give the crude (model-free) counterparts.
#'
#' Prior scales are set with [carriage_priors()].
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
