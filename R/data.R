#' Worked example of the allocation approaches (Table S1)
#'
#' A hypothetical population of ten children, six of them carriers, from
#' Table S1 of the supplement of the PCV10 Nha Trang analysis. It shows how
#' one carriage episode contributes to serotype-specific prevalence under
#' the four allocation approaches and the presence-based approach.
#'
#' @format A data frame with 10 rows and 9 columns:
#' \describe{
#'   \item{id}{Child.}
#'   \item{survey}{Survey round (one round, `"Example"`).}
#'   \item{sero1, sero2, sero3}{Detected serotypes, dominant first.}
#'   \item{ra1, ra2, ra3}{Relative abundance of each detected serotype.}
#'   \item{density}{Total pneumococcal density, genome equivalents/mL.}
#' }
#' @examples
#' dat <- prepare_carriage_data(table_s1, "survey", c("sero1", "sero2", "sero3"),
#'   density = "density", abundance = c("ra1", "ra2", "ra3"))
#' crude <- lapply(c("dominant", "equal", "abundance", "density", "presence"),
#'   function(a) describe_prevalence(dat, a))
#' # number of children: PCV10 and non-PCV10 carriage under each approach
#' sapply(crude, function(x) 10 * x$value[x$group %in% c("VT", "NVT")])
"table_s1"
