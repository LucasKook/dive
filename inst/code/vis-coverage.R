### Results coverage simulation DIVE
### LK 2025

set.seed(12)
save <- TRUE

# DEPs --------------------------------------------------------------------

devtools::load_all()
library("tidyverse")
library("ggpubr")
library("tram")
library("dare")

# Read --------------------------------------------------------------------

odir <- file.path("inst/results/coverage")
files <- list.files(odir, "*.csv", full.names = TRUE)
dat <- tibble(path = files) |>
  mutate(data = map(path, ~ read_csv(.x, show_col_types = FALSE))) |>
  unnest(data) |>
  mutate(path = case_when(
    str_detect(path, "res-1") ~ "Scenario 1",
    str_detect(path, "res-2") ~ "Scenario 2",
    str_detect(path, "res-3") ~ "Scenario 3",
    str_detect(path, "res-4") ~ "Scenario 4"
  ), method = factor(method, levels = c("DIVE", "TRAM"), labels = c("DIVE", "CCDF")))

pd <- dat |>
  group_by(path, n, method, scenario, D, run) |>
  summarize(
    const = sum((ORACLE * (1 - ORACLE))[-1] * diff(Y)),
    cover = sum(((lwr <= ORACLE & upr >= ORACLE) * ORACLE * (1 - ORACLE) / const)[-1] * diff(Y)),
    absdev = sum((abs((lwr <= ORACLE & upr >= ORACLE) - 0.8) * ORACLE * (1 - ORACLE) / const)[-1] * diff(Y))
  )

p1 <- ggplot(pd, aes(x = ordered(n), y = cover, color = method)) +
  ggbeeswarm::geom_quasirandom(aes(shape = factor(D)), alpha = 0.1) +
  stat_summary(size = rel(0.5), fun.data = mean_se, fun.args = list(mult = 0)) +
  geom_hline(yintercept = 0.8, color = "darkred", linetype = 2) +
  theme_bw() +
  labs(x = "Sample size", y = "Average coverage", color = "Method", shape = "Treatment", tag = "A") +
  lims(y = c(0, 1)) +
  theme(legend.position = "top", text = element_text(size = 13.5)) +
  scale_color_brewer(palette = "Dark2") +
  scale_shape_manual(values = c("0" = 4, "1" = 5))

p2 <- ggplot(pd, aes(x = ordered(n), y = absdev, color = method)) +
  ggbeeswarm::geom_quasirandom(aes(shape = factor(D)), alpha = 0.1) +
  stat_summary(size = rel(0.5), fun.data = mean_se, fun.args = list(mult = 0)) +
  theme_bw() +
  labs(x = "Sample size", y = "Average absolute deviation", color = "Method", shape = "Treatment", tag = "B") +
  lims(y = c(0, 1)) +
  theme(legend.position = "top", text = element_text(size = 13.5)) +
  scale_color_brewer(palette = "Dark2") +
  scale_shape_manual(values = c("0" = 4, "1" = 5))

ggarrange(p1, p2, align = "hv", common.legend = TRUE)

if (save) {
  ggsave("inst/figures/sim-coverage.pdf", height = 4.5, width = 11)
}
