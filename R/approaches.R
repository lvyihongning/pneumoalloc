#' Allocation approaches
#'
#' The four approaches to allocating a carriage episode among its detected
#' serotypes. They differ only in the allocation weights passed to the
#' composition component of the unified model and, for `density`, in the
#' carrier weight. The presence-based definition (`presence`), in which
#' every detected serotype counts fully, is available in
#' [allocate_weights()] and [describe_prevalence()] for comparison; it is
#' not an allocation, so it cannot be fitted with [fit_carriage()].
#'
#' @return A tibble with one row per approach: `approach` (1-4), `name`
#'   (used as the `approach` argument elsewhere), `label` and
#'   `use_dens_weight` (1 when carriers are weighted by their density).
#' @export
#' @examples
#' allocation_approaches()
allocation_approaches <- function() {
  tibble::tibble(
    approach = 1:4,
    name = c("dominant", "equal", "abundance", "density"),
    label = c("Dominant serotype", "Equal share", "Relative abundance",
              "Serotype-specific density"),
    use_dens_weight = c(0L, 0L, 0L, 1L)
  )
}

# Approaches requested from fit_carriage(): "all", or one or more names or
# numbers. Returns approach names without duplicates.
expand_approach <- function(approach, call = rlang::caller_env()) {
  if (length(approach) == 0) cli::cli_abort("{.arg approach} is empty.", call = call)
  if ("all" %in% approach) {
    if (length(approach) > 1) {
      cli::cli_abort("Use {.val all} on its own, not together with other approaches.",
                     call = call)
    }
    return(allocation_approaches()$name)
  }
  unique(vapply(as.list(approach), function(a) resolve_approach(a, call = call)$name, ""))
}

# One approach, given by number (1-4) or name, as a list
resolve_approach <- function(approach, allow_presence = FALSE,
                             call = rlang::caller_env()) {
  tab <- allocation_approaches()
  if (allow_presence) {
    tab <- rbind(tab, tibble::tibble(approach = 5L, name = "presence",
                                     label = "Presence-based",
                                     use_dens_weight = 0L))
  }
  idx <- NA_integer_
  if (is.character(approach) && length(approach) == 1 && grepl("^[0-9]+$", approach)) {
    approach <- as.integer(approach)
  }
  if (length(approach) == 1 && !is.na(approach)) {
    idx <- if (is.character(approach)) {
      match(approach, tab$name)
    } else {
      match(as.integer(approach), tab$approach[tab$name != "presence"])
    }
  }
  if (is.na(idx)) {
    cli::cli_abort(
      "{.arg approach} must be one of 1-4 or {.or {.val {tab$name}}}.",
      call = call
    )
  }
  as.list(tab[idx, ])
}
