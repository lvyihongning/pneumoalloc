#' Allocation weights
#'
#' The contribution of each carrier to each serotype category under one
#' approach:
#'
#' * `dominant`: 1 for the category of the dominant serotype;
#' * `equal`: 1/n for each of the n detected serotypes;
#' * `abundance` and `density`: the relative abundance of each detected
#'   serotype. The `density` approach also weights carriers by their
#'   density inside the model (and in [describe_prevalence()]), so its
#'   allocation equals that of `abundance`;
#' * `presence`: 1 for every category with at least one detected serotype
#'   (rows then sum to more than one for multiple carriers).
#'
#' Serotypes pooled into the same category have their weights added, so
#' the rows of the four allocation approaches sum to one.
#'
#' @param dat A [prepare_carriage_data()] object.
#' @param approach Approach number (1-4) or name, see
#'   [allocation_approaches()], or `"presence"`.
#' @param population `"all"` for all carriers or `"multiple"` for children
#'   with multiple serotype carriage.
#'
#' @return A numeric matrix with one row per carrier (in the order of the
#'   carriers in `dat`) and one column per serotype category.
#' @export
#' @examples
#' sim <- simulate_carriage_data(seed = 1)
#' dat <- prepare_carriage_data(sim, "survey", c("sero1", "sero2", "sero3"),
#'   density = "density", abundance = c("ra1", "ra2", "ra3"))
#' W <- allocate_weights(dat, "equal", population = "multiple")
#' head(W[, colSums(W) > 0])
allocate_weights <- function(dat, approach, population = c("all", "multiple")) {
  check_class(dat, "carriage_data")
  population <- rlang::arg_match(population)
  ap <- resolve_approach(approach, allow_presence = TRUE)
  if (population == "all") {
    allocation_matrix(dat$sero_long, dat$N_pos, dat$m_obs, dat$sero_levels,
                      ap$name, dat$has_abundance)
  } else {
    mp <- multiple_population(dat)
    allocation_matrix(mp$sero_long, mp$N, mp$m_obs, mp$sero_levels,
                      ap$name, dat$has_abundance)
  }
}

allocation_matrix <- function(sl, n, m_obs, sero_levels, name, has_abundance,
                              call = rlang::caller_env()) {
  if (name %in% c("abundance", "density") && !has_abundance) {
    cli::cli_abort(
      "Approach {.val {name}} needs relative abundance; supply {.arg abundance} to {.fn prepare_carriage_data}.",
      call = call
    )
  }
  a <- switch(
    name,
    dominant = as.numeric(sl$dominant),
    equal = 1 / m_obs[sl$k],
    abundance = ,
    density = sl$abundance,
    presence = rep(1, nrow(sl))
  )
  cell <- (match(sl$category, sero_levels) - 1L) * n + sl$k
  agg <- rowsum(a, cell)
  A <- matrix(0, nrow = n, ncol = length(sero_levels),
              dimnames = list(NULL, sero_levels))
  A[as.integer(rownames(agg))] <- agg[, 1]
  if (name == "presence") {
    A[] <- as.numeric(A > 0)
  } else if (any(abs(rowSums(A) - 1) > 1e-8)) {
    cli::cli_abort("Internal error: allocation weights do not sum to one.", call = call)
  }
  A
}

# Dense carriers x categories matrix -> sparse form used by the Stan models
build_share_sparse <- function(A, m_sh_max) {
  n <- nrow(A)
  m_sh <- integer(n)
  sh_idx <- matrix(1L, nrow = n, ncol = m_sh_max)
  sh_w <- matrix(0, nrow = n, ncol = m_sh_max)
  for (k in seq_len(n)) {
    nz <- which(A[k, ] > 0)
    mi <- length(nz)
    if (mi < 1L || mi > m_sh_max) {
      cli::cli_abort("Internal error: carrier {k} has {mi} allocation categories.")
    }
    m_sh[k] <- mi
    sh_idx[k, seq_len(mi)] <- nz
    sh_w[k, seq_len(mi)] <- A[k, nz]
  }
  list(m_sh_max = as.integer(m_sh_max), m_sh = m_sh,
       sh_idx = sh_idx, sh_w = sh_w)
}
