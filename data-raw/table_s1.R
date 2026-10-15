# Table S1 of the supplement of the PCV10 Nha Trang analysis: a
# hypothetical population of ten children, six of them carriers,
# illustrating the four allocation approaches and the presence-based
# approach. Density in genome equivalents/mL (Table S1 gives 10^5 GE/mL).
table_s1 <- data.frame(
  id = 1:10,
  survey = factor("Example"),
  sero1 = c("6B", "19F", "15A", "6B", "19F", "6B", NA, NA, NA, NA),
  sero2 = c(NA, NA, NA, "15A", "23A", "15A", NA, NA, NA, NA),
  sero3 = c(NA, NA, NA, NA, NA, "23A", NA, NA, NA, NA),
  ra1 = c(1, 1, 1, 0.7, 0.6, 0.5, NA, NA, NA, NA),
  ra2 = c(NA, NA, NA, 0.3, 0.4, 0.3, NA, NA, NA, NA),
  ra3 = c(NA, NA, NA, NA, NA, 0.2, NA, NA, NA, NA),
  density = c(1, 1, 1, 2, 0.5, 0.5, NA, NA, NA, NA) * 1e5,
  stringsAsFactors = FALSE
)
usethis::use_data(table_s1, overwrite = TRUE)
