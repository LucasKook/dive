# Distributional random forest with control function
# LK March 2023

set.seed(241068)

# Dependencies ------------------------------------------------------------

library("ranger")
library("randomForest")
library("coin")
library("tidyverse")
devtools::load_all()

# Run ---------------------------------------------------------------------

n <- 1e4
prs <- rnorm(6, 1, 0.05)

res <- replicate(5e1, {
  ### Data under intervention on D (d0) and observational (d)
  d0 <- dgp_foster(n, doD = TRUE, prs = prs)
  d1 <- dgp_foster(n, doD = FALSE, prs = prs)

  ### Formulae and data to loop over
  fm_list <- list("DXH" = Y ~ D + X + H, "DX" = Y ~ D + X, "D" = Y ~ D)
  d_list <- list("interventional" = d0, "observational" = d1)

    out <- lapply(fm_list, \(fm) {
      ret <- lapply(d_list, \(d) {
        ### GLM
        p <- glm_marginal_predictions(fm, data = d)
        GLM <- OR(p[, "p1"], p[, "p0"], log)
        CFX <- unname(p[, "cfx"])
        # GLM <- ATE(p[, "p1"], p[, "p0"])

        ### COR
        parCOR <- indep_iv(fm, ~ Z, data = d, "COR")
        pCOR <- indep_marginal_predictions(parCOR, d)
        COR <- OR(pCOR[, "p1"], pCOR[, "p0"], log)
        # COR <- ATE(pCOR[, "p1"], pCOR[, "p0"])

        ### IND
        parIND <- indep_iv(fm, ~ 0 + Z, data = d, "IND")
        pIND <- indep_marginal_predictions(parIND, d)
        IND <- OR(pIND[, "p1"], pIND[, "p0"], log)
        # IND <- ATE(pIND[, "p1"], pIND[, "p0"])

        ### NCTL
        S1 <- glm(D ~ Z, data = d, family = "binomial")
        d$R <- d$D - predict(S1, type = "response")
        pS2 <- glm_marginal_predictions(update(fm, . ~ . + R), data = d)
        NCTL <- OR(pS2[, "p1"], pS2[, "p0"], log)
        # NCTL <- ATE(pS2[, "p1"], pS2[, "p0"])

        ### RF CTRL
        cf <- ranger(D ~ Z, data = d, probability = TRUE)
        d$ps <- d$D - predict(cf, data = d)$predictions[, 2]
        pRF <- ranger_marginal_predictions(update(fm, factor(Y) ~ . + ps), data = d)
        RF <- OR(pRF[, "p1"], pRF[, "p0"], log)
        # RF <- ATE(pRF[, "p1"], pRF[, "p0"])

        c(GLM = GLM, COR = COR, IND = IND, NCTL = NCTL, RF = RF, CFX = CFX)
      })
      names(ret) <- names(d_list)
      bind_rows(ret, .id = "dataset")
    })

    bind_rows(out, .id = "formula") |>
      pivot_longer(GLM:CFX, names_to = "method", values_to = "estimate")

}, simplify = FALSE) %>% bind_rows()

write_csv(res, "inst/results/break-foster.csv")

oracle <- res %>%
  filter(method == "GLM", dataset == "interventional", formula == "DXH") %>%
  pull(estimate) %>% mean()

# Plot and save -----------------------------------------------------------

ggplot(res, aes(x = method, y = estimate)) +
  geom_boxplot() +
  ggbeeswarm::geom_quasirandom(width = 0.1, alpha = 0.3) +
  facet_grid(dataset ~ formula) +
  geom_hline(yintercept = oracle, color = "darkred", linetype = 2) +
  labs(y = "log OR") +
  theme_bw()

ggsave("inst/results/break-foster.pdf", height = 5, width = 7)
