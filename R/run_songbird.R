#' Inverse CLR transform
#' @keywords internal
#' @export
clr_inv <- function(x) {
  expx <- exp(x)
  expx / rowSums(expx)
}

#' CLR transform
#' @keywords internal
#' @export
clr <- function(x) {
  logx <- log(x)
  sweep(logx, 1, rowMeans(logx))
}

#' Fit one Songbird model replicate
#' @keywords internal
#' @export
run_single_model <- function(trainX, trainY,
                             testX,  testY,
                             differential_prior = 1,
                             learning_rate      = 1e-3,
                             clipnorm           = 10,
                             batch_size         = 5,
                             epochs             = 1000) {

  suppressPackageStartupMessages(library(tensorflow))
  tf$compat$v1$disable_eager_execution()
  tfp <- reticulate::import("tensorflow_probability", delay_load = TRUE)
  tfd <- tfp$distributions

  tf$compat$v1$reset_default_graph()

  X <- tf$compat$v1$placeholder(tf$float32, shape(NULL, ncol(trainX)))
  Y <- tf$compat$v1$placeholder(tf$float32, shape(NULL, ncol(trainY)))

  qbeta <- tf$Variable(
    tf$random$normal(
      shape = c(ncol(trainX), ncol(trainY) - 1L),
      mean  = 0,
      stddev = differential_prior
    )
  )

  eta       <- tf$matmul(X, qbeta)
  zero_col  <- tf$zeros(shape = c(tf$shape(eta)[1], 1L), dtype = tf$float32)
  eta_full  <- tf$concat(list(zero_col, eta), axis = 1L)
  phi       <- tf$nn$log_softmax(eta_full)

  prior          <- tfd$Normal(loc = 0, scale = differential_prior)
  log_prior      <- tf$reduce_sum(prior$log_prob(qbeta))
  total_count    <- tf$reduce_sum(Y, axis = 1L)
  likelihood     <- tfd$Multinomial(total_count = total_count, logits = phi)
  log_likelihood <- tf$reduce_sum(likelihood$log_prob(Y))

  loss <- -(log_prior + log_likelihood * (nrow(trainX) / 5))

  optimizer  <- tf$compat$v1$train$AdamOptimizer(learning_rate, 0.9, 0.99)
  grads_vars <- optimizer$compute_gradients(loss)
  grads      <- lapply(grads_vars, `[[`, 1)
  vars       <- lapply(grads_vars, `[[`, 2)
  clipped    <- tf$clip_by_global_norm(grads, clipnorm)[[1]]
  train_step <- optimizer$apply_gradients(Map(list, clipped, vars))

  sess <- tf$compat$v1$Session()
  sess$run(tf$compat$v1$global_variables_initializer())

  for (step in seq_len(epochs)) {
    idx    <- sample(seq_len(nrow(trainX)), batch_size, replace = TRUE)
    batchX <- trainX[idx, , drop = FALSE]
    batchY <- trainY[idx, , drop = FALSE]
    sess$run(train_step, feed_dict = dict(X = batchX, Y = batchY))
  }

  final_beta <- sess$run(qbeta)
  beta_clr   <- clr(clr_inv(cbind(0, final_beta)))
  t(beta_clr)
}

#' Run Songbird-style multinomial regression with bootstrapped p-values
#'
#' @param physeq  A phyloseq object (optional).
#' @param otu_table OTU table matrix (required if physeq not provided).
#' @param metadata Sample-metadata data.frame (required if physeq not provided).
#' @param formula  Model formula as string, e.g. "~ Depth".
#' @param reference_levels Named list to specify reference levels for categorical variables.
#' @param differential_prior  Std. deviation of normal prior on β.
#' @param learning_rate       Learning rate for Adam.
#' @param clipnorm            Gradient-clipping value.
#' @param batch_size          Size of mini-batches.
#' @param seed                RNG seed.
#' @param num_test_samples    Number of test samples per replicate.
#' @param epochs              Training epochs.
#' @param n_bootstrap         Number of bootstrap replicates.
#' @param n_cores             Parallel workers.
#' @param use_gpu             Allow GPU (default TRUE).
#'
#' @return list(beta_mean, beta_sd, beta_pval)
#' @export
run_songbird <- function(physeq      = NULL,
                         otu_table   = NULL,
                         metadata    = NULL,
                         formula,
                         reference_levels = NULL,
                         differential_prior = 1,
                         learning_rate      = 1e-3,
                         clipnorm           = 10,
                         batch_size         = 5,
                         seed               = 0,
                         num_test_samples   = 5,
                         epochs             = 1000,
                         n_bootstrap        = 100,
                         n_cores            = 1,
                         use_gpu            = TRUE) {

  suppressPackageStartupMessages(library(tensorflow))
  library(reticulate)
  library(phyloseq)
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(future.apply))

  tf$compat$v1$disable_eager_execution()
  tf$config$threading$set_intra_op_parallelism_threads(as.integer(n_cores))
  tf$config$threading$set_inter_op_parallelism_threads(as.integer(n_cores))
  if (!use_gpu) {
    tf$config$set_visible_devices(list(), device_type = "GPU")
    message("GPU disabled – running on CPU.")
  }

  if (!is.null(physeq)) {
    otu_mat <- as(phyloseq::otu_table(physeq), "matrix")
    if (phyloseq::taxa_are_rows(phyloseq::otu_table(physeq))) {
      otu_mat <- t(otu_mat)
    }
    metadata <- as(phyloseq::sample_data(physeq), "data.frame")
  } else {
    stopifnot(!is.null(otu_table), !is.null(metadata))
    otu_mat  <- as.matrix(otu_table)
    metadata <- as.data.frame(metadata)
  }

  if (!is.null(reference_levels)) {
    for (var in names(reference_levels)) {
      ref <- reference_levels[[var]]
      metadata[[var]] <- factor(metadata[[var]], levels = c(ref, setdiff(unique(metadata[[var]]), ref)))
    }
  }

  design <- model.matrix(as.formula(formula), metadata)

  set.seed(seed)
  future::plan(multisession, workers = n_cores)

  results <- future.apply::future_lapply(
    seq_len(n_bootstrap),
    function(i) {
      cat(sprintf("Running bootstrap %d of %d\n", i, n_bootstrap)); flush.console()
      idx      <- sample(seq_len(nrow(otu_mat)), replace = TRUE)
      trainX   <- design[idx, , drop = FALSE]
      trainY   <- otu_mat[idx, , drop = FALSE]
      test_idx <- sample(seq_len(nrow(trainX)), num_test_samples)

      run_single_model(
        trainX, trainY,
        trainX[test_idx, , drop = FALSE],
        trainY[test_idx, , drop = FALSE],
        differential_prior, learning_rate, clipnorm, batch_size, epochs
      )
    },
    future.seed = TRUE,
    future.globals = list(
      run_single_model     = songbirdR::run_single_model,
      clr                  = clr,
      clr_inv              = clr_inv,
      differential_prior   = differential_prior,
      learning_rate        = learning_rate,
      clipnorm             = clipnorm,
      batch_size           = batch_size,
      epochs               = epochs,
      num_test_samples     = num_test_samples,
      design               = design,
      otu_mat              = otu_mat
    ),
    future.packages = c("tensorflow", "reticulate")
  )

  beta_array <- simplify2array(results)
  beta_mean  <- apply(beta_array, c(1, 2), mean)
  beta_sd    <- apply(beta_array, c(1, 2), sd)
  beta_pval <- matrix(NA, nrow = dim(beta_array)[1], ncol = dim(beta_array)[2],
                      dimnames = dimnames(beta_mean))

  for (j in 1:dim(beta_array)[2]) {
    for (i in 1:dim(beta_array)[1]) {
      boots <- beta_array[i, j, ]
      obs   <- beta_mean[i, j]
      beta_pval[i, j] <- mean(abs(boots - mean(boots)) >= abs(obs))
    }
  }

  rownames(beta_mean) <- rownames(beta_sd) <- rownames(beta_pval) <- colnames(otu_mat)
  colnames(beta_mean) <- colnames(beta_sd) <- colnames(beta_pval) <- colnames(design)

  list(beta_mean = beta_mean,
       beta_sd   = beta_sd,
       beta_pval = beta_pval,
       beta_array = beta_array)
}
