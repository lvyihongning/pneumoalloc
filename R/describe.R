#' Crude serotype-specific prevalence
#'
#' Model-free estimates by survey round under one approach; the crude
#' counterparts of [fit_carriage()]. For the four allocation approaches,
#' the share of serotype category s among carriers is
#' `sum_k h_k w_ks / sum_k h_k`, where `w_ks` are the allocation weights
#' ([allocate_weights()]) and `h_k` is 1, or \eqn{\Theta_k / \bar{\Theta}_t}
#' (total density relative to the mean among carriers of survey round t)
#' for the `density` approach; within a round the normalisation cancels.
#' Prevalence is the crude carriage prevalence times
#' this share, so serotype prevalences add up to carriage prevalence. For
#' the `presence` approach, prevalence is the proportion of children in
#' whom the serotype (or any vaccine or non-vaccine serotype) is detected.
#'
#' @inheritParams allocate_weights
#' @return A tibble with one row per survey round and outcome (overall,
#'   VT, NVT and each serotype category): `n` (children, or multiple
#'   carriers when `population = "multiple"`), `n_carriers`, `measure` and
#'   the crude `value` (a proportion).
#' @export
#' @examples
#' sim <- simulate_carriage_data(seed = 1)
#' dat <- prepare_carriage_data(sim, "survey", c("sero1", "sero2", "sero3"),
#'   density = "density", abundance = c("ra1", "ra2", "ra3"))
#' crude <- describe_prevalence(dat, "density")
#' crude[crude$group != "serotype", ]
describe_prevalence <- function(dat, approach = "density",
                                population = c("all", "multiple")) {
  check_class(dat, "carriage_data")
  population <- rlang::arg_match(population)
  ap <- resolve_approach(approach, allow_presence = TRUE)
  W <- allocate_weights(dat, ap$name, population)
  lev <- dat$survey_levels
  vt <- colnames(W) %in% dat$vaccine_types

  if (population == "all") {
    k_survey <- factor(dat$survey[dat$pos_idx], levels = lev)
    h <- if (ap$use_dens_weight == 1L) dat$dens_pos else rep(1, dat$N_pos)
    n <- as.vector(table(dat$survey))
    n_carriers <- as.vector(table(k_survey))
    measure <- "prevalence"
  } else {
    mp <- multiple_population(dat)
    k_survey <- factor(dat$survey[mp$rows], levels = lev)
    h <- if (ap$use_dens_weight == 1L) mp$dens else rep(1, mp$N)
    n <- n_carriers <- as.vector(table(k_survey))
    measure <- "proportion among multiple carriers"
  }
  scale <- n_carriers / n

  if (ap$name == "presence") {
    # proportion of children (or multiple carriers) with detection
    share <- rowsum(W, k_survey) / n_carriers
    share_vt <- as.vector(rowsum(as.numeric(rowSums(W[, vt, drop = FALSE]) > 0), k_survey)) / n_carriers
    share_nvt <- as.vector(rowsum(as.numeric(rowSums(W[, !vt, drop = FALSE]) > 0), k_survey)) / n_carriers
  } else {
    share <- rowsum(h * W, k_survey) / as.vector(rowsum(h, k_survey))
    share_vt <- rowSums(share[, vt, drop = FALSE])
    share_nvt <- rowSums(share[, !vt, drop = FALSE])
  }
  share <- share[lev, , drop = FALSE]
  C <- ncol(W)
  value <- unname(c(rep(1, length(lev)), share_vt, share_nvt, as.vector(share)) * scale)

  tibble::tibble(
    approach = ap$approach,
    method = ap$label,
    measure = measure,
    group = c(rep(c("overall", "VT", "NVT"), each = length(lev)),
              rep("serotype", C * length(lev))),
    outcome = c(rep(c("Overall", "VT", "NVT"), each = length(lev)),
                rep(colnames(W), each = length(lev))),
    survey = factor(rep(lev, 3 + C), levels = lev),
    n = rep(n, 3 + C),
    n_carriers = rep(n_carriers, 3 + C),
    value = value
  )
}

#' Crude carriage states and types of multiple carriage
#'
#' Model-free counterparts of [fit_carriage_states()] and
#' [fit_multiple_types()]: by survey round, the proportion of children with
#' no carriage and with 1, 2, ... detected serotypes, carriage (at least
#' one serotype), multiple carriage (at least two), and, among children with
#' multiple carriage, the proportion of each type.
#'
#' @param dat A [prepare_carriage_data()] object.
#' @return A tibble with `measure`, `outcome`, `survey`, `count`, `n` and
#'   the crude proportion `value`.
#' @export
#' @examples
#' sim <- simulate_carriage_data(seed = 1)
#' dat <- prepare_carriage_data(sim, "survey", c("sero1", "sero2", "sero3"),
#'   density = "density")
#' describe_carriage_states(dat)
describe_carriage_states <- function(dat) {
  check_class(dat, "carriage_data")
  lev <- dat$survey_levels
  states <- state_labels(dat$M, dat$capped)
  st <- table(dat$survey, factor(pmin(dat$n_sero, dat$M), levels = 0:dat$M))
  st <- cbind(st, rowSums(st[, -1, drop = FALSE]), rowSums(st[, -(1:2), drop = FALSE]))
  colnames(st) <- c(states, "Carriage", "Multiple carriage")
  n <- as.vector(table(dat$survey))

  multi <- dat$m_obs >= 2
  ty <- table(factor(dat$survey[dat$pos_idx][multi], levels = lev),
              factor(dat$multiple_type[multi], levels = type_labels()))
  n_multi <- rowSums(ty)
  value <- unname(c(as.vector(st) / rep(n, ncol(st)),
                    as.vector(ty) / rep(n_multi, ncol(ty))))

  tibble::tibble(
    measure = c(rep("prevalence", length(st)),
                rep("proportion among multiple carriers", length(ty))),
    outcome = c(rep(colnames(st), each = length(lev)),
                rep(colnames(ty), each = length(lev))),
    survey = factor(rep(lev, ncol(st) + ncol(ty)), levels = lev),
    count = c(as.vector(st), as.vector(ty)),
    n = c(rep(n, ncol(st)), rep(n_multi, ncol(ty))),
    value = value
  )
}

# "1 serotype", "2 serotypes", ...; "3+ serotypes" for a capped top category
count_label <- function(g, plus = FALSE) {
  paste0(g, if (plus) "+", if (g == 1 && !plus) " serotype" else " serotypes")
}

state_labels <- function(M, capped = FALSE) {
  c("No carriage", vapply(seq_len(M), function(g) count_label(g, capped && g == M), ""))
}

type_labels <- function() {
  c("VT-only", "NVT-only", "Mixed, VT-dominant", "Mixed, NVT-dominant")
}

#' Observed density among carriers
#'
#' Model-free summaries of density among carriers by survey round, for all
#' pneumococci and for vaccine-type and non-vaccine-type serotypes, the
#' crude counterparts of [fit_density()]. Group density is the sum of the
#' serotype-specific densities (total density times relative abundance) of
#' the group's serotypes, among children carrying at least one of them.
#'
#' @param dat A [prepare_carriage_data()] object.
#' @inheritParams fit_density
#' @return A tibble with `group`, `survey`, the number of carriers `n`, and
#'   the mean, geometric mean, median and quartiles of density.
#' @export
#' @examples
#' sim <- simulate_carriage_data(seed = 1)
#' dat <- prepare_carriage_data(sim, "survey", c("sero1", "sero2", "sero3"),
#'   density = "density", abundance = c("ra1", "ra2", "ra3"))
#' describe_density(dat)
describe_density <- function(dat, group = "all") {
  check_class(dat, "carriage_data")
  groups <- expand_density_group(group)
  tabs <- lapply(groups, function(g) {
    grp <- density_group(dat, g)
    s <- factor(dat$survey[grp$rows], levels = dat$survey_levels)
    stat <- function(f) as.vector(tapply(grp$density, s, f))
    q <- function(p) function(x) stats::quantile(x, p, names = FALSE)
    tibble::tibble(
      group = g,
      survey = factor(dat$survey_levels, levels = dat$survey_levels),
      n = as.vector(table(s)),
      mean = stat(mean),
      geometric_mean = stat(function(x) exp(mean(log(x)))),
      median = stat(stats::median),
      q25 = stat(q(0.25)),
      q75 = stat(q(0.75))
    )
  })
  tibble::as_tibble(do.call(rbind, tabs))
}

#' Crude proportion of serotype groups within VT or NVT carriage
#'
#' Model-free counterpart of [group_proportions()]: the share of a group of
#' serotypes in VT, NVT or all carriage, by survey round, from the crude
#' allocation-based prevalences of [describe_prevalence()].
#'
#' @inheritParams describe_prevalence
#' @inheritParams group_proportions
#' @return A tibble with `group`, `survey` and the crude proportion `value`.
#' @export
#' @examples
#' sim <- simulate_carriage_data(seed = 1)
#' dat <- prepare_carriage_data(sim, "survey", c("sero1", "sero2", "sero3"),
#'   density = "density", abundance = c("ra1", "ra2", "ra3"))
#' describe_group_proportions(dat, c("15B/C", "23A"), among = "NVT",
#'                            approach = "density")
describe_group_proportions <- function(dat, groups, among = c("NVT", "VT", "all"),
                                       approach = "density",
                                       population = c("all", "multiple")) {
  check_class(dat, "carriage_data")
  among <- rlang::arg_match(among)
  population <- rlang::arg_match(population)
  if (identical(approach, "presence")) {
    cli::cli_abort("Presence-based serotype prevalences do not add up to VT or NVT prevalence.")
  }
  cr <- describe_prevalence(dat, approach, population)
  cr <- cr[cr$group == "serotype", ]
  cats <- unique(cr$outcome)
  g <- resolve_groups(groups, cats, as.integer(cats %in% dat$vaccine_types), among)
  lev <- dat$survey_levels
  by_survey <- function(s) {
    as.vector(tapply(cr$value[cr$outcome %in% s],
                     factor(cr$survey[cr$outcome %in% s], levels = lev), sum))
  }
  denom <- by_survey(g$pool)
  tibble::tibble(
    approach = unique(cr$approach),
    method = unique(cr$method),
    among = among,
    group = rep(names(g$groups), each = length(lev)),
    survey = factor(rep(lev, length(g$groups)), levels = lev),
    value = unlist(lapply(g$groups, function(s) by_survey(s) / denom), use.names = FALSE)
  )
}
