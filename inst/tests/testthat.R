library("testthat")

context("sim_style_data()")

style.all <- c("ERS1", "ERS2", "ARS", "ADRS")
d1 <- sim_style_data(ndimc=1, style=style.all, var.s=.01, reversed=2/3)

test_that("Items are positively correlated and in admittable range", {  
  expect_that(all(apply(d1$dat, 3, cor) > 0), is_true())
  expect_that(all(d1$dat >= 0), is_true())
  expect_that(all(d1$dat <= d1$categories), is_true())
})

test_that("Threshold parameters sum to zero", {
  expect_that(sum(
    d1$item.parameters[which(names(d1$item.parameters) == "categ1"):
                         length(d1$item.parameters)]), equals(0))
})
