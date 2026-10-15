# ---- shared helpers ----------------------------------------------------------

# A fit, a pneumoalloc_fits object or a plain list of fits -> list of fits
as_fit_list <- function(object, arg = rlang::caller_arg(object),
                        call = rlang::caller_env()) {
  if (inherits(object, "pneumoalloc_fit")) return(list(object))
  ok <- is.list(object) && length(object) > 0 &&
    all(vapply(object, inherits, TRUE, "pneumoalloc_fit"))
  if (!ok) {
    cli::cli_abort("{.arg {arg}} must be a {.cls pneumoalloc_fit} object or a list of them.",
                   call = call)
  }
  unname(object)
}

bind_tables <- function(fits, fun, ...) {
  tibble::as_tibble(do.call(rbind, lapply(fits, fun, ...)))
}

method_colours <- function() {
  c(`Dominant serotype` = "#5AB4AC", `Equal share` = "#4393C3",
    `Relative abundance` = "#9970AB", `Serotype-specific density` = "#D6604D",
    `Presence-based` = "grey35")
}

method_fills <- function() {
  c(`Dominant serotype` = "#B2DFDB", `Equal share` = "#9ECAE1",
    `Relative abundance` = "#CAB2D6", `Serotype-specific density` = "#EF8A91",
    `Presence-based` = "grey75")
}

vt_colours <- function() c(VT = "#CC6677", NVT = "#4477AA")

method_levels <- function(methods) {
  known <- names(method_colours())
  c(intersect(known, methods), setdiff(methods, known))
}

wrap_label <- function(x) gsub(" ", "\n", x)

#' ggplot2 theme of pneumoalloc figures
#'
#' A clean theme used by the plotting functions; add it to your own plots
#' for a consistent look.
#'
#' @param base_size Base font size.
#' @return A ggplot2 theme.
#' @export
theme_pneumoalloc <- function(base_size = 11) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey95", colour = NA),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "top",
      legend.title = ggplot2::element_blank(),
      plot.title.position = "plot"
    )
}

# Dashed line after the reference survey and an optional shaded period
time_layers <- function(levels, reference_level, vaccine_intro, shade_from) {
  out <- list()
  if (!is.null(shade_from)) {
    i <- match(shade_from, levels)
    if (is.na(i)) cli::cli_abort("{.arg shade_from} {.val {shade_from}} is not a survey round.",
                                 call = rlang::caller_env(2))
    out <- c(out, list(ggplot2::annotate("rect", xmin = i - 0.5, xmax = Inf,
                                         ymin = -Inf, ymax = Inf,
                                         fill = "grey80", alpha = 0.3)))
  }
  if (isTRUE(vaccine_intro)) {
    out <- c(out, list(ggplot2::geom_vline(xintercept = match(reference_level, levels) + 0.5,
                                           linetype = "dashed", linewidth = 0.35)))
  }
  out
}

#' Plot standardised prevalence by survey round
#'
#' Bars with 95% credible intervals of [marginal_prevalence()], one colour
#' per approach, one panel per outcome (by default VT and NVT carriage, as
#' in Fig 2A-B of the PCV10 Nha Trang analysis).
#'
#' @param object A `pneumoalloc_fit`, a `pneumoalloc_fits` object or a list
#'   of fits, e.g. `c(fits, list(presence = fit_presence(dat)))`.
#' @param groups Groups of outcomes to show, see [marginal_prevalence()].
#' @param measure Which measure to show when a model gives several (e.g.
#'   `"proportion among carriers"`); the first by default.
#' @param vaccine_intro Draw a dashed line after the reference
#'   (pre-vaccination) survey round?
#' @param shade_from Optional survey round from which the background is
#'   shaded, e.g. to mark surveys during the COVID-19 pandemic.
#' @return A ggplot object.
#' @export
#' @examples
#' \dontrun{
#' plot_prevalence(fits)
#' plot_prevalence(c(fits, list(presence = fit_presence(dat))), shade_from = "2019")
#' }
plot_prevalence <- function(object, groups = c("VT", "NVT"), measure = NULL,
                            vaccine_intro = TRUE, shade_from = NULL) {
  fits <- as_fit_list(object)
  tab <- bind_tables(fits, marginal_prevalence, groups = groups)
  measure <- measure %||% tab$measure[1]
  tab <- tab[tab$measure == measure, ]
  dat <- fits[[1]]$data
  tab$method <- factor(tab$method, levels = method_levels(unique(tab$method)))
  tab$outcome <- factor(tab$outcome, levels = unique(tab$outcome))
  dodge <- ggplot2::position_dodge(width = 0.8)
  ggplot2::ggplot(tab, ggplot2::aes(x = .data$survey, y = 100 * .data$median,
                                    fill = .data$method)) +
    time_layers(dat$survey_levels, dat$reference_level, vaccine_intro, shade_from) +
    ggplot2::geom_col(position = dodge, width = 0.75, colour = "black",
                      linewidth = 0.25, alpha = 0.85) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = 100 * .data$q2.5, ymax = 100 * .data$q97.5),
                           position = dodge, width = 0.25, linewidth = 0.3) +
    ggplot2::facet_wrap(~ outcome, scales = "free_y") +
    ggplot2::scale_fill_manual(values = method_fills(), drop = TRUE) +
    ggplot2::scale_x_discrete(labels = wrap_label) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::labs(x = NULL, y = paste0(tools::toTitleCase(measure), " (%)")) +
    theme_pneumoalloc()
}

#' Plot prevalence ratios by survey round
#'
#' Prevalence ratios of [prevalence_ratios()] with 95% credible intervals,
#' one colour per approach, on a log scale with a reference line at 1 (as
#' in Fig 2C-D of the PCV10 Nha Trang analysis).
#'
#' @inheritParams plot_prevalence
#' @inheritParams prevalence_ratios
#' @param measure Which measure to show (`"prevalence"` by default).
#' @return A ggplot object.
#' @export
plot_ratios <- function(object, groups = c("VT", "NVT"), reference = "baseline",
                        measure = NULL, vaccine_intro = FALSE, shade_from = NULL) {
  fits <- as_fit_list(object)
  tab <- bind_tables(fits, prevalence_ratios, groups = groups, reference = reference)
  measure <- measure %||% tab$measure[1]
  tab <- tab[tab$measure == measure, ]
  dat <- fits[[1]]$data
  tab$method <- factor(tab$method, levels = method_levels(unique(tab$method)))
  tab$outcome <- factor(tab$outcome, levels = unique(tab$outcome))
  ref_lab <- if (identical(reference, "previous")) "the previous survey" else unique(tab$reference)[1]
  dodge <- ggplot2::position_dodge(width = 0.6)
  ggplot2::ggplot(tab, ggplot2::aes(x = .data$survey, y = .data$median,
                                    ymin = .data$q2.5, ymax = .data$q97.5,
                                    colour = .data$method)) +
    time_layers(dat$survey_levels, dat$reference_level, vaccine_intro, shade_from) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40",
                        linewidth = 0.35) +
    ggplot2::geom_linerange(position = dodge, linewidth = 0.6) +
    ggplot2::geom_point(position = dodge, size = 1.6) +
    ggplot2::facet_wrap(~ outcome, scales = "free_y") +
    ggplot2::scale_colour_manual(values = method_colours(), drop = TRUE) +
    ggplot2::scale_x_discrete(labels = wrap_label) +
    ggplot2::scale_y_log10(labels = function(x) format(x, scientific = FALSE,
                                                       drop0trailing = TRUE)) +
    ggplot2::labs(x = NULL, y = paste("Ratio versus", ref_lab)) +
    theme_pneumoalloc()
}

#' Plot carriage states or types of multiple carriage
#'
#' Stacked bars of standardised estimates by survey round, as in Fig 1 of the
#' PCV10 Nha Trang analysis: for [fit_carriage_states()], the prevalence of
#' carrying one, two, three, ... serotypes (stacked to overall carriage);
#' for [fit_multiple_types()], the composition of multiple carriage (stacked
#' to 100%). Posterior means are stacked, so the bars add up exactly to
#' overall carriage (or 100%).
#'
#' @param object A fit of [fit_carriage_states()] or [fit_multiple_types()].
#' @inheritParams plot_prevalence
#' @param labels Print the values on the bars?
#' @return A ggplot object.
#' @export
plot_carriage_states <- function(object, vaccine_intro = TRUE, shade_from = NULL,
                                 labels = TRUE) {
  check_class(object, "pneumoalloc_fit")
  if (!object$type %in% c("states", "multiple_types")) {
    cli::cli_abort("{.arg object} must come from {.fn fit_carriage_states} or {.fn fit_multiple_types}.")
  }
  dat <- object$data
  if (object$type == "states") {
    tab <- marginal_prevalence(object, groups = "state")
    tab <- tab[tab$outcome != "No carriage", ]
    ylab <- "Prevalence (%)"
    fills <- c("#A6CEE3", "#CAB2D6", "#FDBF6F", "#B2DF8A", "#FB9A99")
  } else {
    tab <- marginal_prevalence(object)
    ylab <- "Proportion of multiple carriage (%)"
    fills <- c("#F4A3A8", "#9ECAE1", "#FDBB84", "#B2DFDB")
  }
  tab$outcome <- factor(tab$outcome, levels = rev(unique(tab$outcome)))
  tab$value <- 100 * tab$mean
  fills <- stats::setNames(rev(fills[seq_along(levels(tab$outcome))]), levels(tab$outcome))
  p <- ggplot2::ggplot(tab, ggplot2::aes(x = .data$survey, y = .data$value,
                                         fill = .data$outcome)) +
    time_layers(dat$survey_levels, dat$reference_level, vaccine_intro, shade_from) +
    ggplot2::geom_col(width = 0.6, colour = "white", linewidth = 0.2) +
    ggplot2::scale_fill_manual(values = fills, breaks = rev(levels(tab$outcome))) +
    ggplot2::scale_x_discrete(labels = wrap_label) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::labs(x = NULL, y = ylab) +
    theme_pneumoalloc()
  if (labels) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = ifelse(.data$value >= 1, sprintf("%.1f", .data$value), "")),
      position = ggplot2::position_stack(vjust = 0.5), size = 2.8
    )
  }
  p
}

#' Plot serotype-specific prevalence in two survey rounds
#'
#' Prevalence of each serotype category in one survey round against another
#' (log scales, with the line of no change), coloured by vaccine type, as in
#' Fig 4A of the PCV10 Nha Trang analysis. With several approaches, each has
#' its own marker and the points of a serotype are joined.
#'
#' @inheritParams plot_prevalence
#' @param x_survey,y_survey Survey rounds on the horizontal and vertical
#'   axes. Default: the reference round and the last round.
#' @param intervals Draw 95% credible intervals (for the last approach
#'   only, to keep the plot readable)?
#' @param label Label the points with the serotype names? Labels are
#'   placed with 'ggrepel' when it is installed.
#' @return A ggplot object.
#' @export
plot_serotypes <- function(object, x_survey = NULL, y_survey = NULL,
                           intervals = TRUE, label = TRUE) {
  fits <- as_fit_list(object)
  dat <- fits[[1]]$data
  x_survey <- x_survey %||% dat$reference_level
  y_survey <- y_survey %||% utils::tail(dat$survey_levels, 1)
  tab <- bind_tables(fits, marginal_prevalence, groups = "serotype")
  tab <- tab[tab$measure == tab$measure[1], ]
  get <- function(s) tab[tab$survey == s, ]
  a <- get(x_survey)
  b <- get(y_survey)
  key <- function(t) paste(t$method, t$outcome)
  m <- match(key(a), key(b))
  wide <- data.frame(method = a$method, outcome = a$outcome, vt = a$vaccine_type,
                     x = 100 * a$median, x_lo = 100 * a$q2.5, x_hi = 100 * a$q97.5,
                     y = 100 * b$median[m], y_lo = 100 * b$q2.5[m], y_hi = 100 * b$q97.5[m])
  wide$method <- factor(wide$method, levels = method_levels(unique(wide$method)))
  last <- levels(wide$method)[length(levels(wide$method))]
  sub <- wide[wide$method == last, ]
  shapes <- c(16, 1, 2, 5, 4)[seq_along(levels(wide$method))]
  shapes <- stats::setNames(rev(shapes), levels(wide$method))
  p <- ggplot2::ggplot(wide, ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey30",
                         linewidth = 0.35)
  if (length(levels(wide$method)) > 1) {
    p <- p + ggplot2::geom_line(ggplot2::aes(group = .data$outcome), colour = "grey65",
                                linewidth = 0.3, linetype = "dashed")
  }
  if (intervals) {
    p <- p +
      ggplot2::geom_linerange(data = sub, ggplot2::aes(xmin = .data$x_lo, xmax = .data$x_hi,
                                                       colour = .data$vt), alpha = 0.6) +
      ggplot2::geom_linerange(data = sub, ggplot2::aes(ymin = .data$y_lo, ymax = .data$y_hi,
                                                       colour = .data$vt), alpha = 0.6)
  }
  p <- p +
    ggplot2::geom_point(ggplot2::aes(colour = .data$vt, shape = .data$method), size = 2.4,
                        stroke = 0.6) +
    ggplot2::scale_colour_manual(values = vt_colours(),
                                 labels = c(VT = "Vaccine types", NVT = "Non-vaccine types")) +
    ggplot2::scale_shape_manual(values = shapes) +
    ggplot2::scale_x_log10(labels = function(x) format(x, scientific = FALSE, drop0trailing = TRUE)) +
    ggplot2::scale_y_log10(labels = function(x) format(x, scientific = FALSE, drop0trailing = TRUE)) +
    ggplot2::labs(x = paste0("Prevalence in ", x_survey, " (%)"),
                  y = paste0("Prevalence in ", y_survey, " (%)")) +
    theme_pneumoalloc() +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_line(colour = "grey92"),
                   legend.position = "right")
  if (label) {
    lab_aes <- ggplot2::aes(label = .data$outcome, colour = .data$vt)
    p <- p + if (rlang::is_installed("ggrepel")) {
      ggrepel::geom_text_repel(data = sub, lab_aes, size = 3, show.legend = FALSE,
                               box.padding = 0.3, max.overlaps = Inf)
    } else {
      ggplot2::geom_text(data = sub, lab_aes, size = 3, vjust = -0.8, show.legend = FALSE)
    }
  }
  p
}

#' Plot serotype-specific prevalence ratios
#'
#' Prevalence ratios of each serotype category in one survey round versus
#' the reference round, with 95% credible intervals, on a log scale, coloured
#' by vaccine type and dodged by approach (as in Fig 4B of the PCV10 Nha
#' Trang analysis). Ratios flagged as not `reportable` (reference
#' prevalence near zero) are drawn hollow.
#'
#' @inheritParams plot_prevalence
#' @param survey Survey round compared with the reference. Default: the
#'   last round.
#' @inheritParams prevalence_ratios
#' @return A ggplot object.
#' @export
plot_serotype_ratios <- function(object, survey = NULL, reference = "baseline") {
  fits <- as_fit_list(object)
  dat <- fits[[1]]$data
  survey <- survey %||% utils::tail(dat$survey_levels, 1)
  tab <- bind_tables(fits, prevalence_ratios, groups = "serotype", reference = reference)
  tab <- tab[tab$measure == tab$measure[1] & tab$survey == survey, ]
  if (nrow(tab) == 0) cli::cli_abort("No ratios for survey round {.val {survey}}.")
  tab$method <- factor(tab$method, levels = method_levels(unique(tab$method)))
  tab$outcome <- factor(tab$outcome, levels = unique(tab$outcome))
  shapes <- c(16, 15, 17, 18, 4)[seq_along(levels(tab$method))]
  dodge <- ggplot2::position_dodge(width = 0.7)
  ggplot2::ggplot(tab, ggplot2::aes(x = .data$outcome, y = .data$median,
                                    ymin = .data$q2.5, ymax = .data$q97.5,
                                    colour = .data$vaccine_type, shape = .data$method,
                                    group = .data$method)) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40",
                        linewidth = 0.35) +
    ggplot2::geom_linerange(position = dodge, linewidth = 0.5) +
    ggplot2::geom_point(position = dodge, size = 2,
                        ggplot2::aes(alpha = ifelse(.data$reportable, "yes", "no"))) +
    ggplot2::scale_alpha_manual(values = c(yes = 1, no = 0.35), guide = "none") +
    ggplot2::scale_colour_manual(values = vt_colours(),
                                 labels = c(VT = "Vaccine types", NVT = "Non-vaccine types")) +
    ggplot2::scale_shape_manual(values = stats::setNames(shapes, levels(tab$method))) +
    ggplot2::scale_x_discrete(labels = wrap_label) +
    ggplot2::scale_y_log10(labels = function(x) format(x, scientific = FALSE, drop0trailing = TRUE)) +
    ggplot2::labs(x = "Serotype", y = paste0("Prevalence ratio, ", survey, " versus ",
                                               unique(tab$reference)[1])) +
    theme_pneumoalloc()
}

#' Plot density among carriers
#'
#' Observed density among carriers by survey round (violins, log10 scale)
#' with the model-based median (or mean) and its 95% credible interval, one
#' panel per serotype group, as in Fig 3 of the PCV10 Nha Trang analysis.
#'
#' @param object A fit or fits of [fit_density()].
#' @param statistic `"median"` (default) or `"mean"`.
#' @inheritParams plot_prevalence
#' @param observed Draw the observed densities?
#' @return A ggplot object.
#' @export
plot_density <- function(object, statistic = "median", observed = TRUE,
                         vaccine_intro = TRUE, shade_from = NULL) {
  fits <- as_fit_list(object)
  if (!all(vapply(fits, function(f) f$type == "density", TRUE))) {
    cli::cli_abort("{.arg object} must come from {.fn fit_density}.")
  }
  dat <- fits[[1]]$data
  labs_g <- c(total = "All pneumococci", VT = "Vaccine types", NVT = "Non-vaccine types")
  est <- bind_tables(fits, marginal_density, statistic = statistic)
  est <- est[est$outcome == "Density", ]
  est$panel <- factor(labs_g[est$group], levels = labs_g)
  p <- ggplot2::ggplot() +
    time_layers(dat$survey_levels, dat$reference_level, vaccine_intro, shade_from)
  if (observed) {
    obs <- do.call(rbind, lapply(fits, function(f) {
      g <- density_group(dat, f$labels$group)
      data.frame(panel = labs_g[[f$labels$group]],
                 survey = factor(dat$survey[g$rows], levels = dat$survey_levels),
                 log10_density = log10(g$density))
    }))
    obs$panel <- factor(obs$panel, levels = labs_g)
    p <- p + ggplot2::geom_violin(data = obs, ggplot2::aes(x = .data$survey, y = .data$log10_density),
                                  fill = "#F4A3A8", colour = NA, alpha = 0.6, scale = "width",
                                  width = 0.7)
  }
  p +
    ggplot2::geom_linerange(data = est, ggplot2::aes(x = .data$survey, ymin = log10(.data$q2.5),
                                                     ymax = log10(.data$q97.5)), linewidth = 0.7) +
    ggplot2::geom_point(data = est, ggplot2::aes(x = .data$survey, y = log10(.data$median)),
                        size = 1.8) +
    ggplot2::facet_wrap(~ panel) +
    ggplot2::scale_x_discrete(labels = wrap_label, drop = FALSE) +
    ggplot2::labs(x = NULL, y = paste0("Density (log10 GE/mL); model ", statistic, " (95% CrI)")) +
    theme_pneumoalloc()
}
