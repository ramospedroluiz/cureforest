simulate_cure_latency <- function(n, seed, cohort) {
  set.seed(seed)
  cure_marker <- rbinom(n, 1L, 0.5)
  latency_marker <- rnorm(n)
  age <- pmin(85, pmax(30, round(rnorm(n, 61, 11))))
  stage <- factor(
    sample(c("I", "II", "III"), n, TRUE, c(0.30, 0.45, 0.25)),
    levels = c("I", "II", "III")
  )
  biomarker <- rnorm(n)
  noise <- rnorm(n)

  cure_probability <- plogis(
    -0.4 + 3.5 * (cure_marker - 0.5) - 0.25 * (stage == "III")
  )
  cured <- rbinom(n, 1L, cure_probability) == 1L
  susceptible_time <- 1.6 * (
    -log(runif(n)) /
      exp(4.0 * latency_marker + 0.15 * (stage == "III"))
  )^(1 / 1.5)
  event_time <- ifelse(cured, Inf, susceptible_time)
  censoring_time <- rep(10, n)

  data.frame(
    time = pmin(event_time, censoring_time),
    event = as.integer(event_time <= censoring_time),
    cure_marker = factor(cure_marker, labels = c("low", "high")),
    latency_marker = latency_marker,
    age = age,
    stage = stage,
    biomarker = biomarker,
    noise = noise,
    cohort = factor(cohort, levels = c("train", "test"))
  )
}

cure_latency_demo <- rbind(
  simulate_cure_latency(2200L, 20260725L, "train"),
  simulate_cure_latency(1800L, 20260726L, "test")
)
rownames(cure_latency_demo) <- NULL

dir.create("data", showWarnings = FALSE)
save(
  cure_latency_demo,
  file = file.path("data", "cure_latency_demo.rda"),
  compress = "xz"
)
