test_that("errors if bad parameters", {
  library(ExperimentHub)
  Example_Beta<-query(ExperimentHub(), "HiTIMED")[["EH8092"]]
  expect_error(HiTIMED_deconvolution(Example_Beta,
                                   h=7
  ))
  expect_error(HiBED_deconvolution(Example_Beta,"COAD",
                                   h="one"
  ))
  HiTIMED_Library<-query(ExperimentHub(), "HiTIMED")[["EH8093"]]
  expect_error(HiBED_deconvolution(HiTIMED_Library,"COAD",
                                   h=1,"tumor"
  ))
})
