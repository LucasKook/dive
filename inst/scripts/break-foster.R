setting <- commandArgs(TRUE)[1]

# Distributional random forest with control function
# LK March 2023

set.seed(241068)

# Dependencies ------------------------------------------------------------

library("ranger")
library("coin")
library("tidyverse")
devtools::load_all()

### Params
n <- 3e3
nsim <- 1e1
if (is.na(setting))
  setting <- c("cond", "cond-noX", "marg", "marg-noX")[2]
metric <- c("log-OR", "ATE")[2]

if (setting == "cond") {
  prs <- rep(1, 6)
  DGP <- \(n, doD) dgp_foster(n = n, doD = doD, prs = prs)
  fm_list <- list("DXH" = Y ~ X + H, "DX" = Y ~ X, "D" = Y ~ 1)
  coracle <- "DXH"
} else if (setting == "cond-noX") {
  DGP <- dgp_ex1_binary
  fm_list <- list("DH" = Y ~ H, "D" = Y ~ 1)
  coracle <- "DH"
} else if (setting == "marg") {
  prs <- rep(1, 4)
  DGP <- \(n, doD) marginal_dgp_foster(n = n, doD = doD, prs = prs)
  fm_list <- list("DX" = Y ~ X, "D" = Y ~ 1)
  coracle <- "D"
} else if (setting == "marg-noX") {
  DGP <- marginal_dgp_ex1_binary
  fm_list <- list("D" = Y ~ 1)
  coracle <- "D"
} else {
  stop("Setting not implemented.")
}

if (metric == "log-OR") {
  EVAL <- \(p1, p2) OR(p1, p2, log)
} else if (metric == "ATE") {
  EVAL <- \(p1, p2) ATE(p1, p2)
} else {
  stop("Metric not implemented.")
}

### For saving
odir <- file.path("inst", "results", Sys.Date(), setting, metric)
if (!dir.exists(odir))
  dir.create(odir, recursive = TRUE)

# Run ---------------------------------------------------------------------

res <- replicate(nsim, {
  ### Data under intervention on D (d0) and observational (d)
  d0 <- DGP(n, doD = TRUE)
  d1 <- DGP(n, doD = FALSE)
  d_list <- list("interventional" = d0, "observational" = d1)

  ### Formulae and data to loop over

    out <- lapply(fm_list, \(fm) {
      ret <- lapply(d_list, \(d) {
        ### GLM
        p <- glm_marginal_predictions(fm, data = d)
        GLM <- EVAL(p[, "p1"], p[, "p0"])
        CFX <- NULL # unname(p[, "cfx"])
        # if (metric == "ATE")
        #   CFX <- NULL

        ### COR
        parCOR <- indep_iv(update(fm, . ~ . + D), ~ Z, data = d, "COR")
        pCOR <- indep_marginal_predictions(parCOR, d)
        COR <- EVAL(pCOR[, "p1"], pCOR[, "p0"])

        ### IND
        IND <- max(replicate(10, {
          parIND <- indep_iv(update(fm, . ~ + D), ~ 0 + Z, data = d, "IND", ytrafo = rank)
          pIND <- indep_marginal_predictions(parIND, d)
          EVAL(pIND[, "p1"], pIND[, "p0"])
        }))

        ### NCTL
        S1 <- glm(D ~ Z, data = d, family = "binomial")
        d$R <- randomized_pit(1 - predict(S1, type = "response"), d$D)
        pS2 <- glm_marginal_predictions(update(fm, . ~ . + R), data = d)
        NCTL <- EVAL(pS2[, "p1"], pS2[, "p0"])

        ### RF CTRL
        cf <- ranger(factor(D) ~ Z, data = d, probability = TRUE)
        ps <- predict(cf, data = d)$predictions[, 1]
        RF <- mean(replicate(10, {
          d$ps <- randomized_pit(ps, d$D)
          ufm <- update(fm, factor(Y) ~ . + ps)
          pRF <- ranger_marginal_predictions(ufm, data = d)
          EVAL(pRF[, "p1"], pRF[, "p0"])
        }))

        c(GLM = GLM, COR = COR, IND = IND, NCTL = NCTL, CFX = CFX, RF = RF)
      })
      names(ret) <- names(d_list)
      bind_rows(ret, .id = "dataset")
    })

    bind_rows(out, .id = "formula") |>
      pivot_longer(GLM:RF, names_to = "method", values_to = "estimate")

}, simplify = FALSE) %>% bind_rows()

write_csv(res, file.path(odir, "break-foster.csv"))

oracle <- res %>%
  filter(method == "GLM", dataset == "interventional", formula == coracle) %>%
  pull(estimate) %>% mean()

# Plot and save -----------------------------------------------------------

ggplot(res, aes(x = method, y = estimate)) +
  geom_boxplot(outlier.shape = NA) +
  ggbeeswarm::geom_quasirandom(width = 0.1, alpha = 0.3) +
  facet_grid(dataset ~ formula) +
  geom_hline(yintercept = oracle, color = "darkred", linetype = 2) +
  labs(y = metric) +
  theme_bw() +
  labs(caption = paste0("Setting: ", setting))

ggsave(file.path(odir, "break-foster.pdf"))
