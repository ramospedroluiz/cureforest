# Two-minute cureforest versus random survival forest example

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
if (!requireNamespace("randomForestSRC", quietly = TRUE)) {
  install.packages("randomForestSRC")
}

remotes::install_github("ramospedroluiz/cureforest", upgrade = "never")
library(cureforest)

data(cure_latency_demo)
table(cure_latency_demo$cohort)

result <- compare_cureforest(
  formula = survival::Surv(time, event) ~
    cure_marker + latency_marker + age + stage + biomarker + noise,
  data = cure_latency_demo,
  split = "cohort",
  train = "train",
  test = "test",
  ntree = 200,
  maxdepth = 1,
  n.cores = 6,
  seed = 20260725
)

print(result)
