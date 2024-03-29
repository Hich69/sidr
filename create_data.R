
dat = data.frame(USUBJID = 1:1000,
                 VALUE = runif(n = 1000, min = 0.5, max = 0.75),
                 GLD = rbinom(n = 1000, size = 1, prob = 0.5))
data_test <- dat

#save(data_test, file = "data/data_test.RData")

devtools::use_data(data_test, overwrite = TRUE)
