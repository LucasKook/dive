# Distributional random forest with control function
# LK March 2023

set.seed(241068)

# Dependencies ------------------------------------------------------------

library("ranger")
library("randomForest")
library("coin")
devtools::load_all()

# Data --------------------------------------------------------------------

n <- 1e5

### Data under intervention on D (d0) and observational (d)
d0 <- dgp_ex1_binary(n, doD = TRUE)
d1 <- dgp_ex1_binary(n, doD = FALSE)

# Run ---------------------------------------------------------------------

fm_list <- list("DH" = Y ~ D + H, "D" = Y ~ D)
# d_list <- list("interventional" = d0, "observational" = d1)
d_list <- list("observational" = d1)

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
    d$ps <- d$D - predict(cf, data = d)$predictions[, 1]
    pRF <- ranger_marginal_predictions(update(fm, . ~ . + ps), data = d)
    RF <- OR(pRF[, "p1"], pRF[, "p0"], log)
    # RF <- ATE(pRF[, "p1"], pRF[, "p0"])

    c(GLM = GLM, COR = COR, IND = IND, NCTL = NCTL, RF = RF, CFX = CFX)
  })
  names(ret) <- names(d_list)
  dplyr::bind_rows(ret, .id = "dataset")
})

(res <- dplyr::bind_rows(out, .id = "formula"))
