### Simulation DIVE
### LK 2024

set.seed(12)
save <- TRUE

# DEPs --------------------------------------------------------------------

devtools::load_all()
library("tidyverse")
library("tram")
library("dare")

# Read --------------------------------------------------------------------

odir <- file.path("inst/results/simulations", "2024-05-13")
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

dat |> pivot_longer(CmV:KS, names_to = "metric", values_to = "value") |>
  ggplot(aes(x = n, y = value, color = method, group = interaction(n, method))) +
  facet_grid(metric ~ path, scales = "free") +
  geom_boxplot(outlier.shape = NA, position = position_dodge(0.08)) +
  ggbeeswarm::geom_quasirandom(alpha = 0.1, dodge.width = 0.08, size = 0.7,
                               width = 0.04) +
  theme_bw() +
  scale_color_brewer(palette = "Dark2") +
  scale_x_log10(breaks = sort(unique(dat$n))) +
  scale_y_log10() +
  theme(legend.position = "top", text = element_text(size = 13.5),
        axis.text.x.bottom = element_text(angle = 30, hjust = 1, vjust = 1)) +
  labs(x = "Sample size", y = "Estimation error")

if (save)
  ggsave("inst/figures/sim-results.pdf", height = 5, width = 8)
