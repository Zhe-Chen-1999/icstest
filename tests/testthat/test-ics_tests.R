# Test file for ICS tests
# tests/testthat/test-ics_tests.R

# Helper function to generate test data
generate_test_data <- function(n_clusters = 20, min_size = 10, max_size = 30, 
                               informative = FALSE, seed = 123) {
  set.seed(seed)
  
  cluster_sizes <- sample(min_size:max_size, n_clusters, replace = TRUE)
  n <- sum(cluster_sizes)
  cluster_id <- rep(1:n_clusters, cluster_sizes)
  Z <- rep(rbinom(n_clusters, 1, 0.5), cluster_sizes)
  
  # Generate outcomes
  Y0 <- rnorm(n) + rep(rnorm(n_clusters), cluster_sizes)
  
  if (informative) {
    # Make treatment effect depend on cluster size
    tau <- rep(0.5 * (cluster_sizes > median(cluster_sizes)), cluster_sizes)
  } else {
    # Constant treatment effect
    tau <- rep(0.5, n)
  }
  
  Y <- Y0 + tau * Z
  
  list(Y = Y, Z = Z, cluster_id = cluster_id, cluster_sizes = cluster_sizes, n = n)
}

# Tests for model_assisted_test
test_that("model_assisted_test works with basic input", {
  data <- generate_test_data()
  
  result <- model_assisted_test(data$Y, data$Z, data$cluster_id)
  
  # Check output structure
  expect_type(result, "list")
  expect_named(result, c("method", "t_stat", "p_value"))
  expect_equal(result$method, "Model-Assisted Test")
  expect_type(result$t_stat, "double")
  expect_type(result$p_value, "double")
  expect_true(result$p_value >= 0 && result$p_value <= 1)
})

test_that("model_assisted_test handles covariates", {
  data <- generate_test_data()
  
  # Create cluster-level covariates
  covariates <- data.frame(
    x1 = rnorm(length(unique(data$cluster_id))),
    x2 = runif(length(unique(data$cluster_id)))
  )
  
  result <- model_assisted_test(data$Y, data$Z, data$cluster_id, 
                                covariates = covariates)
  
  expect_type(result, "list")
  expect_named(result, c("method", "t_stat", "p_value"))
  expect_true(result$p_value >= 0 && result$p_value <= 1)
})


# Tests for randomization_test
test_that("randomization_test works with basic input", {
  data <- generate_test_data()
  
  result <- randomization_test(data$Y, data$Z, data$cluster_id, n_perms = 100)
  
  # Check output structure
  expect_type(result, "list")
  expect_named(result, c("method", "statistic", "p_value"))
  expect_equal(result$method, "Randomization-Based Test")
  expect_type(result$statistic, "double")
  expect_type(result$p_value, "double")
  expect_true(result$p_value >= 0 && result$p_value <= 1)
})


# Tests for model_based_test
test_that("model_based_test works with basic input", {
  skip_if_not_installed("geepack")
  
  data <- generate_test_data()
  
  # Create data frame for model
  model_data <- data.frame(
    Y = data$Y,
    Z = data$Z,
    cluster_size = data$cluster_sizes[data$cluster_id],
    cluster_id = data$cluster_id
  )
  
  formula <- as.formula("Y ~ Z * cluster_size")
  
  result <- model_based_test(formula, "Z:cluster_size", model_data, 
                             model_data$cluster_id)
  
  expect_type(result, "double")
  expect_true(result >= 0 && result <= 1)
  expect_true(is.finite(result))
})

test_that("model_based_test error handling works", {
  skip_if_not_installed("geepack")
  
  data <- generate_test_data()
  
  model_data <- data.frame(
    Y = data$Y,
    Z = data$Z,
    cluster_size = data$cluster_sizes[data$cluster_id],
    cluster_id = data$cluster_id
  )
  
  formula <- as.formula("Y ~ Z * cluster_size")
  
  # Test with wrong interaction term name
  expect_error(
    model_based_test(formula, "Z:wrong_name", model_data, model_data$cluster_id),
    "Interaction term 'Z:wrong_name' not found in the model"
  )
})

# Integration tests
test_that("all three methods give reasonable results on same data", {
  skip_if_not_installed("geepack")
  
  # Generate data with informative cluster sizes
  data <- generate_test_data(informative = TRUE)
  
  # Model-assisted test
  ma_result <- model_assisted_test(data$Y, data$Z, data$cluster_id)
  
  # Randomization test (small number of perms for speed)
  rand_result <- randomization_test(data$Y, data$Z, data$cluster_id, n_perms = 100)
  
  # Model-based test
  model_data <- data.frame(
    Y = data$Y,
    Z = data$Z,
    cluster_size = data$cluster_sizes[data$cluster_id],
    cluster_id = data$cluster_id
  )
  formula <- as.formula("Y ~ Z * cluster_size")
  mb_result <- model_based_test(formula, "Z:cluster_size", model_data, 
                                model_data$cluster_id)
  
  # All should give valid p-values
  expect_true(ma_result$p_value >= 0 && ma_result$p_value <= 1)
  expect_true(rand_result$p_value >= 0 && rand_result$p_value <= 1)
  expect_true(mb_result >= 0 && mb_result <= 1)
  
})

test_that("calculate_test_statistic helper function works", {
  data <- generate_test_data()
  
  # Calculate inputs for helper function
  Yvec <- tapply(data$Y, data$cluster_id, mean)
  Zvec <- tapply(data$Z, data$cluster_id, mean)
  cluster_sizes <- as.numeric(table(data$cluster_id))
  num_clusters <- length(unique(data$cluster_id))
  n <- length(data$Y)
  w <- cluster_sizes / n - 1 / num_clusters
  wYvec <- w * Yvec * num_clusters
  
  # Test helper function
  stat <- calculate_test_statistic(wYvec, Zvec, cluster_sizes, NULL, TRUE)
  
  expect_type(stat, "double")
  expect_true(is.finite(stat))
})
