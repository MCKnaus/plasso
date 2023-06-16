test_that("norm_w_to_n handles non-numeric weights", {
  
  w = c("a", "b", "c")
  
  expect_error(norm_w_to_n(w), "Weights do not have the required format (vector or n x 1 matrix)", fixed=TRUE)
})


test_that("norm_w_to_n handles valid weights (no treatment indicators)", {
  w = c(0.2, 0.3, 0.5)
  result = norm_w_to_n(w)
  
  expect_equal(sum(result), length(w))
})

test_that("norm_w_to_n handles valid weights (with treatment indicators)", {
  w = c(0.2, 0.3, 0.5)
  d = c(1, 0, 1)
  result = norm_w_to_n(w, d)
  
  expect_equal(sum(result[d == 1]), sum(result[d == 0]))
  expect_equal(sum(result), length(w)*2)
})

test_that("norm_w_to_n handles matrix weights (with treatment indicators)", {
  w = matrix(c(1, 1, 1), nrow = 3, ncol = 1)
  d = c(1, 0, 1)
  result = norm_w_to_n(w, d)
  
  expect_equal(sum(result[d == 1]), sum(result[d == 0]))
  expect_equal(sum(result), nrow(w)*2)
  expect_identical(result, as.matrix(c(1.5, 3, 1.5)))
})