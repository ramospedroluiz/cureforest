# cureforest versus an ordinary random survival forest

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
if (!requireNamespace("randomForestSRC", quietly = TRUE)) {
  install.packages("randomForestSRC")
}

remotes::install_github("ramospedroluiz/cureforest", upgrade = "never")
library(cureforest)
library(survival)

data(cure_latency_demo)
train <- subset(cure_latency_demo, cohort == "train")
test <- subset(cure_latency_demo, cohort == "test")
model_formula <- Surv(time, event) ~
  cure_marker + latency_marker + age + stage + biomarker + noise

fit.cure <- cureforest(
  model_formula, data = train, ntree = 200, mtry = 6,
  nodesize = 80, maxdepth = 1,
  tail.start = 7, tail.end = 9.5,
  n.cores = 6, seed = 20260725
)

fit.rsf <- randomForestSRC::rfsrc(
  model_formula, data = train, ntree = 200, mtry = 6,
  nodesize = 80, nodedepth = 1, splitrule = "logrank",
  nthread = 6, seed = 20260725
)

pred.cure <- predict(fit.cure, newdata = test, type = "cure")
pred.rsf <- predict(fit.rsf, newdata = test)
head(pred.cure)
head(pred.rsf$predicted)

comparison <- compare_forest_fits(
  fit.cure, fit.rsf, newdata = test, n.cores = 6
)
print(comparison)
