
library("tram")
library("survival")
library("trtf")
library("ranger")
library("randomForest")

set.seed(6)

# Data --------------------------------------------------------------------

data("GBSG2", package = "TH.data")
GBSG2$surv <- with(GBSG2, Surv(time, cens))
GBSG2$time <- as.double(GBSG2$time)

# Model -------------------------------------------------------------------

rf <- ranger(time ~ horTh, data = GBSG2, quantreg = TRUE, num.trees = 5e3)
m <- Coxph(time | horTh ~ 1, data = GBSG2, prob = c(0.001, 0.999), order = 10)

ps <- seq(1e-2, 1 - 1e-2, length.out = 1e2)
p0 <- predict(rf, type = "quantiles", data = nd <- data.frame(horTh = unique(GBSG2$horTh)), quantiles = ps)

pm0 <- as.double(predict(m, newdata = nd[1,,drop=FALSE], type = "quantile", prob = ps))
pm1 <- as.double(predict(m, newdata = nd[2,,drop=FALSE], type = "quantile", prob = ps))

plot(ps, p0$predictions[1, ], type = "l")
lines(ps, p0$predictions[2, ], col = 2)

plot(ps, apply(p0$predictions, 2, diff), type = "l")
lines(ps, pm1 - pm0)

tn <- predict(rf, data = nd[1,,drop=FALSE], type = "terminalNodes")
summary(tn$predictions[1,])

# CDF forest --------------------------------------------------------------

learn <- sample.int(nrow(GBSG2), size = floor(0.8 * nrow(GBSG2)))


crf <- randomForest(time ~ horTh + age + menostat + tsize + tgrade + pnodes +
                      progrec + estrec, data = GBSG2[learn,])

rfw <- predict(crf, newdata = GBSG2, proximity = TRUE)$proximity[learn,][,-learn]
rfw <- rfw / matrix(pmax(colSums(rfw), .Machine$double.eps),
                    ncol = nrow(GBSG2) - length(learn),
                    nrow = length(learn), byrow = TRUE)

idx0 <- which(GBSG2$horTh[-learn] == "no")
idx1 <- which(GBSG2$horTh[-learn] == "yes")

pcdf <- Vectorize(\(y, idx) c(t(rfw[, idx]) %*% as.numeric(GBSG2$time[learn] <= y)), "y")
ts <- seq(min(GBSG2$time), max(GBSG2$time), length.out = 1e2)
plot(ts, colMeans(pcdf(ts, idx0)), type = "l")
lines(ts, colMeans(pcdf(ts, idx1)), type = "l", col = 2)
plot(ecdf(GBSG2$time[GBSG2$horTh == "no"]), add = TRUE, col = rgb(.1, .1, .1, .1))
plot(ecdf(GBSG2$time[GBSG2$horTh == "yes"]), add = TRUE, col = rgb(.9, .1, .1, .1))

# m0 <- Coxph(surv ~ 1, data = GBSG2, prob = c(0.001, 0.999))
# ctrl <- ctree_control(minsplit = 50, minbucket = 20, mincriterion = 0)
# tf <- traforest(m0, formula = surv ~ horTh, control = ctrl, ntree = 50,
#                 mtry = 1, trace = TRUE, data = GBSG2)
#
# predict(tf, newdata = data.frame(surv = 1:10, horTh = unique(GBSG2$horTh)),
#         type = "distribution")
