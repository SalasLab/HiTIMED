test_that("errors if bad parameters", {
  library(HiTIMED)
  data("Example_Beta")
  expect_error(HiTIMED_deconvolution(Example_Beta,
                                   h=7
  ))
})
