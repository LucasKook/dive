### Simulation DIVE
### LK 2024

set.seed(12)

### File names
save <- TRUE

# DEPs --------------------------------------------------------------------

devtools::load_all()
library("tidyverse")
library("tram")
library("dare")

# Settings ----------------------------------------------------------------

dgp <- \(scenario) switch(
  scenario,
  "1" = function(n = 1e2, do = FALSE) {
    Z <- rlogis(n)
    H <- rlogis(n)
    D <- as.numeric(4 * Z + (1 - do) * 4 * H > rlogis(n))
    Y <- -10 + 8 * D + 6 * H
    data.frame(Y = Y, H = H, D = D, Z = Z)
  },
  "2" = function(n = 1e2, do = FALSE) {
    Z <- rlogis(n)
    H <- rlogis(n)
    D <- as.numeric(4 * Z + (1 - do) * 4 * H > rlogis(n))
    Y <- 16 * D + 6 * H + H * rlogis(n)
    data.frame(Y = Y, H = H, D = D, Z = Z)
  },
  "3" = function(n = 1e2, do = FALSE) {
    Z <- rlogis(n)
    H <- rlogis(n)
    D <- as.numeric(4 * Z + (1 - do) * 4 * H > rlogis(n))
    Y <- log(1 + exp(18 + 8 * D + 6 * H))
    data.frame(Y = Y, H = H, D = D, Z = Z)
  },
  "4" = function(n = 1e2, do = FALSE) {
    Z <- rlogis(n)
    H <- rlogis(n)
    D <- as.numeric(4 * Z + (1 - do) * 4 * H > rlogis(n))
    Y <- -6 * D + (6 + 3 * D) * H
    data.frame(Y = Y, H = H, D = D, Z = Z)
  }
)

# ORACLE ------------------------------------------------------------------

pd <- do.call("rbind", lapply(1:4, \(scen) {
  ### Interventional
  dint <- dgp(scen)(1e4, do = TRUE)
  dint$mode <- "interventional"
  ### Observational
  dobs <- dgp(scen)(1e4)
  dobs$mode <- "observational"
  bind_rows(dobs, dint) |> mutate(scenario = scen)
}))

ggplot(pd, aes(x = Y, color = factor(D), linetype = mode)) +
  facet_grid(~ scenario, scales = "free",
             labeller = label_bquote(cols = Scenario~.(scenario))) +
  stat_ecdf() +
  theme_bw() +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Y", y = "CDF", linetype = element_blank(), color = "D") +
  theme(legend.position = "top", text = element_text(size = 13.5)) +
  scale_linetype_manual(values = c(4, 1))

ggsave("inst/figures/sim-scenarios.pdf", height = 3, width = 9)
