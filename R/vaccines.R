#' Vaccine serotype sets
#'
#' Serotypes covered by pneumococcal vaccines, for the `vaccine_types`
#' argument of [prepare_carriage_data()]. The serotypes given there form the
#' vaccine-type (VT) group; every other detected serotype, including
#' non-encapsulated pneumococci, forms the non-vaccine-type (NVT) group. Any
#' character vector can be used, so the VT group can be defined freely, for
#' example to reclassify a cross-protected serotype as vaccine-associated.
#'
#' Serotype names must be spelled as in your data. Serotyping methods
#' differ: microarray reports e.g. `"15B/C"` where a vaccine list has
#' `"15B"`, and `"6C"` or `"20A"`/`"20B"` may be reported separately. Use
#' `add` and `remove` to match your nomenclature, and check the result
#' against the vaccine's product information.
#'
#' | vaccine | serotypes |
#' |---|---|
#' | PCV7 | 4, 6B, 9V, 14, 18C, 19F, 23F |
#' | PCV10 | PCV7 + 1, 5, 7F (Synflorix) |
#' | PCV10-SII | 1, 5, 6A, 6B, 7F, 9V, 14, 19A, 19F, 23F (Pneumosil) |
#' | PCV13 | PCV10 + 3, 6A, 19A |
#' | PCV15 | PCV13 + 22F, 33F |
#' | PCV20 | PCV15 + 8, 10A, 11A, 12F, 15B |
#' | PCV21 | 3, 6A, 7F, 8, 9N, 10A, 11A, 12F, 15A, 15C, 16F, 17F, 19A, 20A, 22F, 23A, 23B, 24F, 31, 33F, 35B (V116) |
#' | PPSV23 | 1, 2, 3, 4, 5, 6B, 7F, 8, 9N, 9V, 10A, 11A, 12F, 14, 15B, 17F, 18C, 19A, 19F, 20, 22F, 23F, 33F |
#'
#' @param vaccine One or more of the vaccines in the table; with several,
#'   their serotypes are combined. `NULL` starts from an empty set, so that
#'   `add` defines a set of your own.
#' @param add Serotypes to add to the set.
#' @param remove Serotypes to remove from the set.
#' @return A character vector of serotypes.
#' @export
#' @examples
#' vaccine_types("PCV13")
#'
#' # PCV10 plus 6A as vaccine-associated (cross-protection)
#' vaccine_types("PCV10", add = "6A")
#'
#' # PCV20 with microarray nomenclature for 15B
#' vaccine_types("PCV20", add = "15B/C", remove = "15B")
#'
#' # serotypes in PCV20 but not in PCV13
#' setdiff(pcv20_types(), pcv13_types())
#'
#' # a set of your own (any character vector also works directly in
#' # prepare_carriage_data(vaccine_types = ...))
#' vaccine_types(NULL, add = c("6B", "14", "19F", "23F"))
vaccine_types <- function(vaccine = "PCV10", add = NULL, remove = NULL) {
  sets <- vaccine_sets()
  if (is.null(vaccine)) {
    if (length(add) == 0) cli::cli_abort("With {.code vaccine = NULL}, give the serotypes in {.arg add}.")
    return(setdiff(unique(as.character(add)), remove))
  }
  unknown <- setdiff(vaccine, names(sets))
  if (length(vaccine) == 0 || length(unknown) > 0) {
    cli::cli_abort(c(
      "Unknown vaccine{?s}: {.val {unknown}}.",
      "i" = "Available: {.val {names(sets)}}."
    ))
  }
  out <- unique(c(unlist(sets[vaccine], use.names = FALSE), add))
  setdiff(out, remove)
}

#' @rdname vaccine_types
#' @export
pcv10_types <- function() vaccine_types("PCV10")

#' @rdname vaccine_types
#' @export
pcv13_types <- function() vaccine_types("PCV13")

#' @rdname vaccine_types
#' @export
pcv15_types <- function() vaccine_types("PCV15")

#' @rdname vaccine_types
#' @export
pcv20_types <- function() vaccine_types("PCV20")

vaccine_sets <- function() {
  pcv7 <- c("4", "6B", "9V", "14", "18C", "19F", "23F")
  pcv10 <- c("1", "4", "5", "6B", "7F", "9V", "14", "18C", "19F", "23F")
  pcv13 <- c(pcv10, "3", "6A", "19A")
  pcv15 <- c(pcv13, "22F", "33F")
  pcv20 <- c(pcv15, "8", "10A", "11A", "12F", "15B")
  list(
    PCV7 = pcv7,
    PCV10 = pcv10,
    `PCV10-SII` = c("1", "5", "6A", "6B", "7F", "9V", "14", "19A", "19F", "23F"),
    PCV13 = pcv13,
    PCV15 = pcv15,
    PCV20 = pcv20,
    PCV21 = c("3", "6A", "7F", "8", "9N", "10A", "11A", "12F", "15A", "15C",
              "16F", "17F", "19A", "20A", "22F", "23A", "23B", "24F", "31",
              "33F", "35B"),
    PPSV23 = c("1", "2", "3", "4", "5", "6B", "7F", "8", "9N", "9V", "10A",
               "11A", "12F", "14", "15B", "17F", "18C", "19A", "19F", "20",
               "22F", "23F", "33F")
  )
}
