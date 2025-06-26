#' Model-Assisted Test for Informative Cluster Size
#'
#' Performs a model-assisted test for informative cluster sizes in cluster
#' randomized trials.
#'
#' @param Y Numeric vector of individual-level outcomes
#' @param Z Numeric vector of individual-level treatment assignments (0/1)
#' @param cluster_id Vector identifying cluster membership
#' @param covariates Optional matrix/data.frame of cluster-level covariates
#' @param adjust_size Logical indicating whether to adjust for cluster size (default TRUE)
#'
#' @return List containing test statistic and p-value
#' @export
#' @import lmtest sandwich dplyr
#' @examples
#' # Generate example data
#' n_clusters <- 200
#' cluster_sizes <- sample(20:80, n_clusters, replace = TRUE)
#' n <- sum(cluster_sizes)
#' cluster_id <- rep(1:n_clusters, cluster_sizes)
#' Z <- rep(rbinom(n_clusters, 1, 0.5), cluster_sizes)
#' Y0 <- rnorm(n) + rep(rnorm(n_clusters), cluster_sizes)
#' tau <- rnorm(n) + rep(0.5 * (cluster_sizes > 50), cluster_sizes)
#' Y <- Y0 + tau * Z
#'
#' # Run test
#' result <- model_assisted_test(Y, Z, cluster_id)
#' print(result)
model_assisted_test <- function(Y, Z, cluster_id, covariates = NULL,
                                adjust_size = TRUE) {

  # Input validation
  if (length(Y) != length(Z) || length(Y) != length(cluster_id)) {
    stop("Y, Z, and cluster_id must have the same length")
  }

  # Calculate cluster-level summaries
  Yvec <- tapply(Y, cluster_id, mean)
  Zvec <- tapply(Z, cluster_id, mean)
  cluster_sizes <- as.numeric(table(cluster_id))
  num_clusters <- length(unique(cluster_id))
  n <- length(Y)

  # Calculate weights for model-assisted approach
  w <- cluster_sizes / n - 1 / num_clusters
  wYvec <- w * Yvec * num_clusters

  # Build model formula
  if (is.null(covariates)) {
    model_formula <- wYvec ~ Zvec
    model_data <- data.frame(wYvec = wYvec, Zvec = Zvec)
  } else {
    # Adjust for cluster-level covariates
    if (nrow(covariates) != num_clusters) {
      stop("Number of rows in covariates must equal number of clusters")
    }

    # Center covariate columns
    covariates <- covariates %>%
      mutate(across(everything(), ~ . - mean(., na.rm = TRUE)))

    model_data <- data.frame(wYvec = wYvec, Zvec = Zvec, covariates)

    # Add cluster size adjustment if requested
    if (adjust_size) {
      cs_cent <- cluster_sizes - mean(cluster_sizes)
      model_data$cs_cent <- cs_cent

      # Create interaction terms with treatment
      covar_names <- paste0("Zvec:", names(covariates))
      size_interaction <- "Zvec:cs_cent"

      model_formula <- as.formula(paste("wYvec ~ Zvec +",
                                        paste(c(names(covariates), "cs_cent",
                                                covar_names, size_interaction),
                                              collapse = " + ")))
    } else {
       covar_names <- paste0("Zvec:", names(covariates))
       model_formula <- as.formula(paste("wYvec ~ Zvec +",
                                         paste(c(names(covariates), covar_names),
                                               collapse = " + ")))
      #covar_names <- paste0("Zvec*", names(covariates))
      #model_formula <- as.formula(paste("wYvec ~ ", paste(covar_names, collapse = " + ")))
    }
  }

  # Fit model and get robust standard errors
  model <- lm(model_formula, data = model_data)
  robust_vcov <- vcovHC(model, type = "HC0")
  test_result <- coeftest(model, vcov = robust_vcov)

  # Return t-statistic and p-value
  t_stat <- test_result["Zvec", "t value"]
  p_value <- test_result["Zvec", "Pr(>|t|)"]

  return(list(method = "Model-Assisted Test",
              t_stat = t_stat,
              p_value = p_value))
}

#' Randomization-Based Test for Informative Cluster Size
#'
#' Performs a randomization-based test for informative cluster sizes in cluster
#' randomized trials.
#'
#' @param Y Numeric vector of individual-level outcomes
#' @param Z Numeric vector of individual-level treatment assignments (0/1)
#' @param cluster_id Vector identifying cluster membership
#' @param covariates Optional matrix/data.frame of cluster-level covariates
#' @param adjust_size Logical indicating whether to adjust for cluster size
#' @param n_perms Number of permutations for randomization test (default 5000)
#'
#' @return List containing test statistic and permutation-based p-value
#' @export
#' @import lmtest sandwich dplyr
#' @examples
#' # Generate example data
#' n_clusters <- 200
#' cluster_sizes <- sample(20:80, n_clusters, replace = TRUE)
#' n <- sum(cluster_sizes)
#' cluster_id <- rep(1:n_clusters, cluster_sizes)
#' Z <- rep(rbinom(n_clusters, 1, 0.5), cluster_sizes)
#' Y0 <- rnorm(n) + rep(rnorm(n_clusters), cluster_sizes)
#' tau <- rnorm(n) + rep(0.5 * (cluster_sizes > 50), cluster_sizes)
#' Y <- Y0 + tau * Z
#'
#' # Run test
#' result <- randomization_test(Y, Z, cluster_id, n_perms = 5000)
#' print(result)
randomization_test <- function(Y, Z, cluster_id, covariates = NULL,
                               adjust_size = TRUE, n_perms = 5000) {

  # Input validation
  if (length(Y) != length(Z) || length(Y) != length(cluster_id)) {
    stop("Y, Z, and cluster_id must have the same length")
  }

  # Calculate cluster-level summaries
  Yvec <- tapply(Y, cluster_id, mean)
  Zvec <- tapply(Z, cluster_id, mean)
  cluster_sizes <- as.numeric(table(cluster_id))
  num_clusters <- length(unique(cluster_id))
  n <- length(Y)

  # Calculate weights
  w <- cluster_sizes / n - 1 / num_clusters
  wYvec <- w * Yvec * num_clusters

  # Number of treated clusters
  n_trt_clusters <- sum(Zvec)

  # Calculate observed test statistic
  obs_stat <- calculate_test_statistic(wYvec, Zvec, cluster_sizes, covariates, adjust_size)

  # Permutation test
  perm_stats <- numeric(n_perms)
  for (i in 1:n_perms) {
    # Permute cluster-level treatment assignment
    trt_perm <- sample(1:num_clusters, n_trt_clusters)
    Zvec_perm <- numeric(num_clusters)
    Zvec_perm[trt_perm] <- 1

    # Calculate permuted test statistic
    perm_stats[i] <- calculate_test_statistic(wYvec, Zvec_perm, cluster_sizes,
                                              covariates, adjust_size)
  }

  # Calculate p-value
  p_value <- mean(abs(perm_stats) >= abs(obs_stat))

  list(method = "Randomization-Based Test",
       statistic = obs_stat,
       p_value = p_value
  )
}

#' Helper function to calculate test statistic
#' @keywords internal
calculate_test_statistic <- function(wYvec, Zvec, cluster_sizes, covariates, adjust_size) {

  num_clusters <- length(Zvec)

  if (is.null(covariates)) {
    model_data <- data.frame(wYvec = wYvec, Zvec = Zvec)
    model_formula <- wYvec ~ Zvec
  } else {

    # Center covariate columns
    covariates <- covariates %>%
      mutate(across(everything(), ~ . - mean(., na.rm = TRUE)))

    model_data <- data.frame(wYvec = wYvec, Zvec = Zvec, covariates)

    if (adjust_size) {
      cs_cent <- cluster_sizes - mean(cluster_sizes)
      model_data$cs_cent <- cs_cent

      covar_names <- paste0("Zvec:", names(covariates))
      size_interaction <- "Zvec:cs_cent"

      model_formula <- as.formula(paste("wYvec ~ Zvec +",
                                        paste(c(names(covariates), "cs_cent",
                                                covar_names, size_interaction),
                                              collapse = " + ")))
    } else {
      covar_names <- paste0("Zvec:", names(covariates))
      model_formula <- as.formula(paste("wYvec ~ Zvec +",
                                        paste(c(names(covariates), covar_names),
                                              collapse = " + ")))
    }
  }

  model <- lm(model_formula, data = model_data)
  robust_vcov <- vcovHC(model, type = "HC0")
  test_result <- coeftest(model, vcov = robust_vcov)

  # Return t-statistic
  test_result["Zvec", "t value"]

}

#' Model-Based Test for Informative Cluster Size
#'
#' Performs a model-based test using generalized estimating equations for
#' informative cluster sizes in cluster randomized trials.
#'
#' @param model_formula An R formula (class `formula`) specifying the model,
#'   including the interaction term between treatment and cluster size.
#' @param interaction_term_name A character string specifying the name of the interaction term between treatment and cluster size
#'   (e.g., `"Z:cluster_size"`) for which the p-value will be returned.
#' @param model_data A data frame containing all the variables required by the `model_formula`:
#'   the outcome, treatment assignment, cluster size, cluster index, and any covariates to be adjusted for.
#' @param cluster_id A vector or column in `model_data` indicating the cluster membership of each unit.
#'
#' @return A numeric value, the p-value for the interaction term between treatment and cluster size.
#' @export
#' @import geepack
#' @examples
#' library(geepack)
#'
#' # Generate example data
#' n_clusters <- 20
#' cluster_sizes <- sample(20:80, n_clusters, replace = TRUE)
#' n <- sum(cluster_sizes)
#' cluster_id <- rep(1:n_clusters, cluster_sizes)
#' cluster_size <- cluster_sizes[cluster_id]
#' Z <- rep(rbinom(n_clusters, 1, 0.5), cluster_sizes)
#' Y0 <- rnorm(n) + rep(rnorm(n_clusters), cluster_sizes)
#' tau <- 0.5 * log(cluster_size)
#' Y <- Y0 + tau * Z
#'
#' data <- data.frame(Y = Y, Z = Z, cluster_size = cluster_size,
#'                          cluster_id = cluster_id)
#'
#' formula <- as.formula("Y ~ Z * log(cluster_size)")
#'
#' # Run the test
#' model_based_test(formula, "Z:log(cluster_size)", data, cluster_id)
model_based_test <- function(model_formula, interaction_term_name, model_data, cluster_id) {

  # Fit the GEE model
  fit <- geepack::geeglm(model_formula, id = cluster_id, data = model_data,
                         corstr = "independence")

  # Get coefficient names
  coef_names <- rownames(summary(fit)$coefficients)

  # Check if interaction_term_name is in the model
  if (!interaction_term_name %in% coef_names) {
    stop(paste0("Interaction term '", interaction_term_name, "' not found in the model.\n",
                "Available terms are:\n", paste(coef_names, collapse = ", ")))
  }

  # Extract p-value
  p_value <- summary(fit)$coefficients[interaction_term_name, "Pr(>|W|)"]

  return(p_value)
}
