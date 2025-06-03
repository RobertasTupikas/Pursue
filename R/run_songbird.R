#' Run Songbird-style multinomial regression in R using a phyloseq object
#'
#' @param physeq A `phyloseq` object (optional). Use this if your data is already in a phyloseq format.
#' @param otu_table A numeric matrix or data.frame of OTU counts (samples as rows, features as columns). Required if `physeq` is not provided.
#' @param metadata A data.frame of sample metadata (samples as rows, features as columns). Required if `physeq` is not provided.
#' @param formula A model formula as a string (e.g., "~ Depth")
#' @param differential_prior Standard deviation of the normal prior on beta coefficients (default = 1.0)
#' @param learning_rate Learning rate for Adam optimizer (default = 0.001)
#' @param clipnorm Gradient clipping value (default = 10.0)
#' @param batch_size Size of training batches (default = 5)
#' @param seed Random seed for R and TensorFlow (default = 0)
#' @param num_test_samples Number of random test samples (default = 5)
#' @param epochs Number of training epochs (default = 1000)
#'
#' @return A list with beta coefficients
#'
run_songbird <- function(physeq = NULL,
                         otu_table = NULL,
                         metadata = NULL,
                         formula,
                         differential_prior = 1.0,
                         learning_rate = 0.001,
                         clipnorm = 10.0,
                         batch_size = 5,
                         seed = 0,
                         num_test_samples = 5,
                         epochs = 1000) {

  library(tensorflow)
  library(reticulate)
  library(phyloseq)
  library(dplyr)

  # --- INPUT VALIDATION ---

  # Check formula
  if (missing(formula) || !is.character(formula)) {
    stop("Argument 'formula' must be provided as a character string (e.g., \"~ treatment\").")
  }

  # Check that learning_rate is positive
  if (!is.numeric(learning_rate) || length(learning_rate) != 1 || learning_rate <= 0) {
    stop("learning_rate must be a single positive numeric value.")
  }

  # Check that batch_size is a positive integer
  if (!is.numeric(batch_size) || batch_size <= 0 || batch_size != as.integer(batch_size)) {
    stop("batch_size must be a positive integer.")
  }

  # Check that num_test_samples is a non-negative integer
  if (!is.numeric(num_test_samples) || num_test_samples < 0 || num_test_samples != as.integer(num_test_samples)) {
    stop("num_test_samples must be a non-negative integer.")
  }

  # Check that epochs is a positive integer
  if (!is.numeric(epochs) || epochs <= 0 || epochs != as.integer(epochs)) {
    stop("epochs must be a positive integer.")
  }

  # Check prior
  if (!is.numeric(differential_prior) || differential_prior <= 0) {
    stop("differential_prior must be a positive number.")
  }

  # If physeq not provided, check otu_table and metadata are non-null and matching
  if (is.null(physeq)) {
    if (is.null(otu_table) || is.null(metadata)) {
      stop("Either 'physeq' or both 'otu_table' and 'metadata' must be provided.")
    }

    if (!is.matrix(otu_table) && !is.data.frame(otu_table)) {
      stop("otu_table must be a matrix or data.frame.")
    }

    if (!is.data.frame(metadata)) {
      stop("metadata must be a data.frame.")
    }

    if (nrow(otu_table) != nrow(metadata)) {
      stop("Number of rows (samples) must match between 'otu_table' and 'metadata'.")
    }
  }

  # Check num_test_samples is less than number of samples
  n_samples <- if (!is.null(physeq)) nsamples(physeq) else nrow(otu_table)
  if (num_test_samples >= n_samples) {
    stop("num_test_samples must be less than the number of samples (", n_samples, ").")
  }

  # User must run: use_condaenv() externally
  tfp <- reticulate::import("tensorflow_probability", delay_load = TRUE)
  tfd <- tfp$distributions

  tf$compat$v1$disable_eager_execution()
  tf$compat$v1$reset_default_graph()
  tf$random$set_seed(as.integer(seed))
  set.seed(seed)

  # Input: phyloseq or separate tables
  if (!is.null(physeq)) {
    otu_mat <- as(otu_table(physeq), "matrix")
    if (taxa_are_rows(otu_table(physeq))) {
      otu_mat <- t(otu_mat)
    }
    metadata <- as(sample_data(physeq), "data.frame")
  } else if (!is.null(otu_table) && !is.null(metadata)) {
    otu_mat <- as.matrix(otu_table)
    metadata <- as.data.frame(metadata)
  } else {
    stop("Please provide either a 'physeq' object OR both 'otu_table' and 'metadata'.")
  }

  if (nrow(otu_mat) != nrow(metadata)) {
    stop("The number of samples (rows) must match between OTU table and metadata.")
  }

  # Design matrix
  design <- model.matrix(as.formula(formula), metadata)

  # Tensor shapes
  trainX <- design
  trainY <- otu_mat
  D <- as.integer(ncol(trainY))
  p <- as.integer(ncol(trainX))

  # Placeholders and model
  X <- tf$compat$v1$placeholder(tf$float32, shape(NULL, p))
  Y <- tf$compat$v1$placeholder(tf$float32, shape(NULL, D))
  qbeta <- tf$Variable(tf$random$normal(shape = c(p, D - 1L), mean = 0.0, stddev = differential_prior))
  eta <- tf$matmul(X, qbeta)
  zero_col <- tf$zeros(shape = c(tf$shape(eta)[1], 1L), dtype = tf$float32)
  eta_full <- tf$concat(list(zero_col, eta), axis = 1L)
  phi <- tf$nn$log_softmax(eta_full)

  # Loss function
  prior <- tfd$Normal(loc = 0.0, scale = differential_prior)
  log_prior <- tf$reduce_sum(prior$log_prob(qbeta))
  total_count <- tf$reduce_sum(Y, axis = 1L)
  likelihood <- tfd$Multinomial(total_count = total_count, logits = phi)
  log_likelihood <- tf$reduce_sum(likelihood$log_prob(Y))
  loss <- -(log_prior + log_likelihood * (nrow(trainX) / 5))

  # Optimizer
  optimizer <- tf$compat$v1$train$AdamOptimizer(learning_rate, beta1 = 0.9, beta2 = 0.99)
  grads_vars <- optimizer$compute_gradients(loss)
  grads <- lapply(grads_vars, `[[`, 1)
  vars <- lapply(grads_vars, `[[`, 2)
  clipped <- tf$clip_by_global_norm(grads, clipnorm)[[1]]
  train_step <- optimizer$apply_gradients(Map(list, clipped, vars))

  # Split data
  all_idx <- 1:nrow(trainX)
  test_idx <- sample(all_idx, num_test_samples)
  train_idx <- setdiff(all_idx, test_idx)

  trainX_split <- trainX[train_idx, ]
  trainY_split <- trainY[train_idx, ]
  testX_split <- trainX[test_idx, ]
  testY_split <- trainY[test_idx, ]

  # CLR functions
  clr_inv <- function(x) {
    expx <- exp(x)
    expx / rowSums(expx)
  }
  clr <- function(x) {
    logx <- log(x)
    sweep(logx, 1, rowMeans(logx))
  }

  # Run session
  sess <- tf$compat$v1$Session()
  sess$run(tf$compat$v1$global_variables_initializer())
  for (step in 1:epochs) {
    idx <- sample(1:nrow(trainX_split), batch_size, replace = TRUE)
    batchX <- trainX_split[idx, , drop = FALSE]
    batchY <- trainY_split[idx, , drop = FALSE]
    sess$run(train_step, feed_dict = dict(X = batchX, Y = batchY))
    if (step %% 1000 == 0) {
      test_loss <- sess$run(loss,
                               feed_dict = dict(X = testX_split, Y = testY_split))
      cat("Epoch", step, "Test loss:", test_loss, "\n")
    }
  }

  # Extract coefficients
  final_beta <- sess$run(qbeta)
  beta_clr <- clr(clr_inv(cbind(0, final_beta)))
  beta_clr <- t(beta_clr)

  colnames(beta_clr) <- colnames(trainX)
  rownames(beta_clr) <- colnames(trainY)

  return(beta = beta_clr)
}
