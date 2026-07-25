# Contributing

Issues and focused pull requests are welcome. For changes to the statistical
method, please describe the estimand, censoring assumptions, follow-up
maturity requirements, and the tests used to validate the change.

Before opening a pull request, run:

```r
devtools::test()
devtools::check()
```

Please do not include registry data, patient-level data, credentials, or
generated check directories in commits.
