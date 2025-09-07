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


#' Fit one regression model
#' @keywords internal
#' @export
run_single_model <- function(trainX, trainY,
                             testX,  testY,
                             differential_prior = 1,
                             learning_rate      = 1e-3,
                             clipnorm           = 10,
                             batch_size         = 5,
                             epochs             = 1000,
                             re_info            = NULL,
                             re_lambda          = 1e-3,
                             early_stopping     = TRUE,
                             patience           = 2,
                             min_epochs         = 200,
                             eval_every         = 1000,
                             tol                = 1e-4,
                             lr_decay_on_plateau= TRUE,
                             lr_decay_factor    = 0.5,
                             lr_min             = 1e-6,
                             verbose            = FALSE,
                             tb_logdir          = NULL,
                             print_prefix       = "") {

  suppressPackageStartupMessages(library(tensorflow))
  tf$compat$v1$disable_eager_execution()
  tfp <- reticulate::import("tensorflow_probability", delay_load = TRUE)
  tfd <- tfp$distributions

  tf$compat$v1$reset_default_graph()

  # TF Placeholders
  X   <- tf$compat$v1$placeholder(tf$float32, shape(NULL, ncol(trainX)))
  Y   <- tf$compat$v1$placeholder(tf$float32, shape(NULL, ncol(trainY)))
  IDX <- tf$compat$v1$placeholder(tf$int32,   shape(NULL))

  # Coefficients (K-1) relative to implicit reference class
  qbeta <- tf$Variable(
    tf$random$normal(
      shape = c(ncol(trainX), ncol(trainY) - 1L),
      mean  = 0,
      stddev = differential_prior
    ), name = "qbeta"
  )

  # Linear predictor (K-1), then expand to K with a zero column
  eta       <- tf$matmul(X, qbeta)
  zero_col  <- tf$zeros(shape = c(tf$shape(eta)[1], 1L), dtype = tf$float32)
  eta_full  <- tf$concat(list(zero_col, eta), axis = 1L)

  #Random effects
  re_contribs <- list()
  re_penalty  <- tf$zeros(shape(), dtype = tf$float32)

  if (!is.null(re_info) && length(re_info) > 0) {
    K <- as.integer(ncol(trainY))
    for (i in seq_along(re_info)) {
      n_levels <- as.integer(re_info[[i]]$n_levels)
      idx_full <- tf$constant(as.integer(re_info[[i]]$idx0), dtype = tf$int32)

      RE_raw <- tf$Variable(tf$zeros(shape = c(n_levels, K), dtype = tf$float32),
                            name = paste0("RE_", re_info[[i]]$name))
      RE_ctr <- RE_raw - tf$reduce_mean(RE_raw, axis = as.integer(1), keepdims = TRUE)

      lvl_mb <- tf$gather(params = idx_full, indices = IDX)
      re_mb  <- tf$gather(params = RE_ctr,  indices = lvl_mb)

      re_contribs[[length(re_contribs) + 1L]] <- re_mb
      re_penalty <- re_penalty + tf$constant(re_lambda, tf$float32) * tf$nn$l2_loss(RE_ctr)
    }
  }

  logits <- if (length(re_contribs) > 0) {
    Reduce(function(a, b) a + b, re_contribs) + eta_full
  } else {
    eta_full
  }

  # Softmax
  phi <- tf$nn$log_softmax(logits)

  # Priors & likelihood
  prior          <- tfd$Normal(loc = 0, scale = differential_prior)
  log_prior      <- tf$reduce_sum(prior$log_prob(qbeta))
  total_count    <- tf$reduce_sum(Y, axis = 1L)
  likelihood     <- tfd$Multinomial(total_count = total_count, logits = phi)
  log_likelihood <- tf$reduce_sum(likelihood$log_prob(Y))

  # Loss formulation
  loss <- -(log_prior + log_likelihood * (nrow(trainX) / 5))
  if (length(re_contribs) > 0) loss <- loss + re_penalty

  # Optimizer
  lr_var     <- tf$Variable(learning_rate, trainable = FALSE, dtype = tf$float32, name = "lr")
  optimizer  <- tf$compat$v1$train$AdamOptimizer(lr_var, 0.9, 0.99)
  gv         <- optimizer$compute_gradients(loss)
  grads      <- lapply(gv, `[[`, 1)
  vars       <- lapply(gv, `[[`, 2)
  clipped    <- tf$clip_by_global_norm(grads, clipnorm)[[1]]
  train_step <- optimizer$apply_gradients(Map(list, clipped, vars))

  #evaluate full-batch loss given row index vector
  eval_loss <- function(sess, Xmat, Ymat, idx0) {
    sess$run(loss, feed_dict = dict(X = Xmat, Y = Ymat, IDX = as.integer(idx0)))
  }

  # Session initialization
  sess <- tf$compat$v1$Session()
  sess$run(tf$compat$v1$global_variables_initializer())

  # TensorBoard setup
  use_tb <- !is.null(tb_logdir) && nzchar(tb_logdir)
  if (use_tb) {
    dir.create(tb_logdir, recursive = TRUE, showWarnings = FALSE)
    writer  <- tf$compat$v1$summary$FileWriter(tb_logdir, sess$graph)
    val_ph  <- tf$compat$v1$placeholder(tf$float32, shape())
    s_train <- tf$compat$v1$summary$scalar("train_loss", loss)
    s_val   <- tf$compat$v1$summary$scalar("val_loss",  val_ph)
    s_lr    <- tf$compat$v1$summary$scalar("lr", lr_var)
    merged  <- tf$compat$v1$summary$merge(list(s_train, s_lr))
  }

  n_train <- nrow(trainX)

  # Internal validation split
  has_val <- n_train >= 10
  if (has_val) {
    vsize   <- max(5L, ceiling(0.1 * n_train))
    val_idx <- sample(seq_len(n_train), vsize)
    tr_idx  <- setdiff(seq_len(n_train), val_idx)
    X_tr    <- trainX[tr_idx, , drop = FALSE]
    Y_tr    <- trainY[tr_idx, , drop = FALSE]
    X_va    <- trainX[val_idx, , drop = FALSE]
    Y_va    <- trainY[val_idx, , drop = FALSE]
    idx_tr0 <- as.integer(tr_idx - 1L)
    idx_va0 <- as.integer(val_idx - 1L)
  } else {
    X_tr <- trainX; Y_tr <- trainY
    idx_tr0 <- as.integer(seq_len(n_train) - 1L)
    X_va <- NULL; Y_va <- NULL; idx_va0 <- NULL
  }

  # Early-stopping bookkeeping
  best_val <- Inf
  best_weights <- NULL
  no_improve <- 0L
  steps_used <- 0L

  # Training loop
  for (step in seq_len(epochs)) {
    idx_rows <- sample(seq_len(n_train), batch_size, replace = TRUE)
    batchX   <- trainX[idx_rows, , drop = FALSE]
    batchY   <- trainY[idx_rows, , drop = FALSE]
    sess$run(train_step, feed_dict = dict(
      X = batchX, Y = batchY, IDX = as.integer(idx_rows - 1L)
    ))

    if (step %% eval_every == 0L) {
      tr_loss <- eval_loss(sess, X_tr, Y_tr, idx_tr0)
      va_loss <- if (has_val) eval_loss(sess, X_va, Y_va, idx_va0) else NA_real_

      # TensorBoard logging
      if (use_tb) {
        writer$add_summary(
          sess$run(merged, feed_dict = dict(X = X_tr, Y = Y_tr, IDX = as.integer(idx_tr0))),
          step
        )
        if (has_val && is.finite(va_loss)) {
          writer$add_summary(sess$run(s_val, feed_dict = dict(val_ph = va_loss)), step)
        }
      }

      # Early-stopping check
      improved <- FALSE
      if (early_stopping && has_val) {
        rel_gain <- (best_val - va_loss) / (abs(best_val) + 1e-12)
        if (is.finite(va_loss) && (va_loss + tol < best_val || rel_gain > tol)) {
          best_val <- va_loss
          best_weights <- tf$compat$v1$trainable_variables() |> lapply(sess$run)
          no_improve <- 0L
          improved <- TRUE
        } else {
          no_improve <- no_improve + 1L
          if (lr_decay_on_plateau && no_improve %% patience == 0L) {
            new_lr <- max(lr_min, sess$run(lr_var) * lr_decay_factor)
            sess$run(lr_var$assign(new_lr))
          }
        }
      }

      if (verbose) {
        cur_lr <- sess$run(lr_var)
        cat(sprintf(
          "%s%s epoch %d | train=%.6f | val=%s | best=%s | lr=%.3g%s\n",
          if (nchar(print_prefix)) paste0(print_prefix, " | ") else "",
          "regression",
          step,
          tr_loss,
          if (is.na(va_loss)) "NA" else sprintf("%.6f", va_loss),
          if (is.infinite(best_val)) "Inf" else sprintf("%.6f", best_val),
          cur_lr,
          if (improved) " *" else ""
        ))
        flush.console()
      }

      if (early_stopping && has_val && step >= min_epochs && no_improve >= patience) {
        steps_used <- step
        if (verbose) {
          cat(sprintf("%sEarly stopping at epoch %d (no improvement in %d evaluations)\n",
                      if (nchar(print_prefix)) paste0(print_prefix, " | ") else "",
                      step, no_improve))
          flush.console()
        }
        break
      }
    }

    steps_used <- step
  }

  # Final TensorBoard close
  if (use_tb) writer$close()

  # Restore best weights if ever improved
  if (!is.null(best_weights)) {
    trainables <- tf$compat$v1$trainable_variables()
    for (i in seq_along(trainables)) {
      sess$run(trainables[[i]]$assign(best_weights[[i]]))
    }
  }

  # Extract coefficients and convert
  final_beta <- sess$run(qbeta)
  clr_inv <- function(x) { expx <- exp(x); expx / rowSums(expx) }
  clr     <- function(x) { logx <- log(x); sweep(logx, 1, rowMeans(logx)) }
  beta_clr <- clr(clr_inv(cbind(0, final_beta)))

  # Return taxa x coefs matrix
  t(beta_clr)
}



#' Data handling and passing regression models for bootstrapping
#' @keywords internal
#' @export

run_regression <- function(physeq = NULL,
                                       otu_table = NULL,
                                       metadata  = NULL,
                                       formula,
                                       random_effects    = NULL,
                                       reference_levels  = NULL,
                                       resample = c("cluster","row"),
                                       cluster_var = NULL,
                                       differential_prior = 1,
                                       learning_rate      = 1e-3,
                                       clipnorm           = 10,
                                       batch_size         = 5,
                                       seed               = 0,
                                       num_test_samples   = 5,
                                       epochs             = 1000,
                                       n_boot             = 100,
                                       n_cores            = 1,
                                       parallel           = c("auto","no","multicore","snow"),
                                       re_lambda          = 1e-3,
                                       early_stopping     = TRUE,
                                       patience           = 20,
                                       min_epochs         = 200,
                                       eval_every         = 10,
                                       tol                = 1e-4,
                                       lr_decay_on_plateau= TRUE,
                                       lr_decay_factor    = 0.5,
                                       lr_min             = 1e-6,
                                       verbose            = FALSE,
                                       tb_logdir          = NULL,
                                       print_prefix       = "") {

  suppressPackageStartupMessages(library(boot))

  if (!is.null(physeq)) {
    suppressPackageStartupMessages(library(phyloseq))
    metadata  <- as(phyloseq::sample_data(physeq), "data.frame")
    otu_table <- as(phyloseq::otu_table(physeq), "matrix")
    if (phyloseq::taxa_are_rows(phyloseq::otu_table(physeq))) otu_table <- t(otu_table)
  }
  stopifnot(!is.null(otu_table), !is.null(metadata))
  otu_mat  <- as.matrix(otu_table)
  metadata <- as.data.frame(metadata)

  # Reference levels for fixed effects
  if (!is.null(reference_levels)) {
    for (var in names(reference_levels)) {
      ref <- reference_levels[[var]]
      metadata[[var]] <- factor(metadata[[var]],
                                levels = c(ref, setdiff(unique(metadata[[var]]), ref)))
    }
  }

  # Fixed-effects design and dims
  design <- model.matrix(as.formula(formula), metadata)
  K <- ncol(otu_mat)
  P <- ncol(design)
  taxa_names <- colnames(otu_mat)
  coef_names <- colnames(design)

  # Random-effects variables
  re_vars <- character(0)
  if (!is.null(random_effects)) {
    if (inherits(random_effects, "formula")) {
      re_vars <- all.vars(update(random_effects, . ~ .))
    } else if (is.character(random_effects)) {
      re_vars <- random_effects
    } else {
      stop("random_effects must be NULL, a formula, or a character vector.")
    }
    re_vars <- unique(re_vars)
  }


  # Choose boot parallel mode
  parallel <- match.arg(parallel)

  if (parallel == "auto") {
    parallel <- if (.Platform$OS.type == "windows") "snow" else if (n_cores > 1) "multicore" else "no"
  }

  resample <- match.arg(resample)

  if (resample == "cluster") {
    if (is.null(cluster_var)) stop("cluster_var must be provided when resample='cluster'.")
    clusters    <- as.factor(metadata[[cluster_var]])
    cluster_ids <- levels(clusters)
    n_clusters  <- length(cluster_ids)
    if (n_clusters < 2L) stop("Need ≥2 clusters for clustered resampling.")
  }

  boot_counter   <- 0L
  local_verbose  <- isTRUE(verbose) && parallel == "no"
  base_prefix    <- if (nzchar(print_prefix)) print_prefix else "boot"

  # Statistic function for bootstraps
  stat_fun <- function(data, indices) {

    # follow bootstrap counter for diagnostics (only in serial mode)
    if (parallel == "no") { boot_counter <<- boot_counter + 1L; boot_id <- boot_counter } else { boot_id <- NA_integer_ }

    # map resample indices -> row indices
    keep_idx <- if (resample == "cluster") {
      picked <- cluster_ids[indices]
      unlist(lapply(picked, function(cl) which(clusters == cl)), use.names = FALSE)
    } else {
      indices
    }
    if (length(keep_idx) == 0L) stop("Empty bootstrap resample.")

    trainX <- design[keep_idx, , drop = FALSE]
    trainY <- otu_mat[keep_idx, , drop = FALSE]

    # RE info (built on full metadata, indexed by keep_idx)
    re_info <- if (length(re_vars) > 0) lapply(re_vars, function(v) {
      fac_all <- base::addNA(factor(metadata[[v]]), ifany = TRUE)
      fac_sub <- factor(fac_all[keep_idx], levels = levels(fac_all))
      lvl_idx <- as.integer(fac_sub); if (anyNA(lvl_idx)) lvl_idx[is.na(lvl_idx)] <- nlevels(fac_sub)
      list(name = v, n_levels = as.integer(nlevels(fac_sub)), idx0 = as.integer(lvl_idx - 1L))
    }) else NULL
    if (!is.null(re_info)) names(re_info) <- re_vars

    # validation
    test_idx <- if (nrow(trainX) > 0) sample(seq_len(nrow(trainX)), min(num_test_samples, nrow(trainX))) else integer(0)

    # per-replicate TensorBoard dir (serial only; requires run_songbird has tb_logdir in scope)
    this_tb_dir <- if (parallel == "no" && !is.null(tb_logdir) && nzchar(tb_logdir))
      file.path(tb_logdir, sprintf("rep_%03d", boot_id)) else NULL

    # fit one replicate
    if (parallel == "no") {
      fit_b <- run_single_model(
        trainX, trainY,
        if (length(test_idx)) trainX[test_idx, , drop = FALSE] else trainX[0, , drop = FALSE],
        if (length(test_idx)) trainY[test_idx, , drop = FALSE] else trainY[0, , drop = FALSE],
        differential_prior, learning_rate, clipnorm, batch_size, epochs,
        re_info = re_info, re_lambda = re_lambda,
        early_stopping      = early_stopping,
        patience            = patience,
        min_epochs          = min_epochs,
        eval_every          = eval_every,
        tol                 = tol,
        lr_decay_on_plateau = lr_decay_on_plateau,
        lr_decay_factor     = lr_decay_factor,
        lr_min              = lr_min,
        verbose             = local_verbose,
        print_prefix        = if (local_verbose) sprintf("%s %d", base_prefix, boot_id) else "",
        tb_logdir           = this_tb_dir
      )
    } else {
      fit_b <- run_single_model(
        trainX, trainY,
        if (length(test_idx)) trainX[test_idx, , drop = FALSE] else trainX[0, , drop = FALSE],
        if (length(test_idx)) trainY[test_idx, , drop = FALSE] else trainY[0, , drop = FALSE],
        differential_prior, learning_rate, clipnorm, batch_size, epochs,
        re_info = re_info, re_lambda = re_lambda,
        early_stopping      = early_stopping,
        patience            = patience,
        min_epochs          = min_epochs,
        eval_every          = eval_every,
        tol                 = tol,
        lr_decay_on_plateau = lr_decay_on_plateau,
        lr_decay_factor     = lr_decay_factor,
        lr_min              = lr_min,
        verbose             = FALSE,
        print_prefix        = "",
        tb_logdir           = NULL
      )
    }

    beta_b <- if (is.matrix(fit_b)) fit_b else if (!is.null(fit_b$beta_KP)) fit_b$beta_KP else fit_b$beta_mean
    as.vector(beta_b)
  }

  # Data for boot function
  boot_data <- if (resample == "cluster") seq_len(n_clusters) else seq_len(nrow(design))

  # Run bootstraps
  set.seed(seed)
  boot_out <- boot::boot(data = boot_data, statistic = stat_fun, R = n_boot,
                         sim = "ordinary", parallel = parallel, ncpus = n_cores)

  # Reassemble coefficient arrays/matrices from boot output
  t0_vec <- boot_out$t0
  t_mat  <- boot_out$t

  # If t0 wasn't computed (rare), fall back to mean of boots (but normally present)
  if (length(t0_vec) != K*P) {
    t0_vec <- colMeans(t_mat, na.rm = TRUE)
  }

  # Final output
  B <- nrow(t_mat)
  beta_array <- array(NA_real_, dim = c(K, P, B),
                      dimnames = list(taxa_names, coef_names, paste0("boot", seq_len(B))))
  for (b in seq_len(B)) {
    beta_array[ , , b] <- matrix(t_mat[b, ], nrow = K, ncol = P, byrow = FALSE,
                                 dimnames = list(taxa_names, coef_names))
  }

  beta_mean <- matrix(t0_vec, nrow = K, ncol = P, byrow = FALSE,
                      dimnames = list(taxa_names, coef_names))
  beta_sd   <- apply(beta_array, c(1, 2), sd, na.rm = TRUE)

  beta_pval <- matrix(NA_real_, nrow = K, ncol = P,
                      dimnames = list(taxa_names, coef_names))
  for (i in seq_len(K)) {
    for (j in seq_len(P)) {
      boots <- beta_array[i, j, ]
      t0    <- beta_mean[i, j]
      mboot <- mean(boots, na.rm=TRUE)
      obs_c <- t0 - mboot
      B     <- sum(!is.na(boots))
      beta_pval[i, j] <- (sum(abs(boots - mboot) >= abs(obs_c), na.rm=TRUE) + 1) / (B + 1)
    }
  }


  list(
    boot       = boot_out,
    beta_mean  = beta_mean,
    beta_sd    = beta_sd,
    beta_pval  = beta_pval,
    beta_array = beta_array,
    dims       = list(K = K, P = P, taxa_names = taxa_names, coef_names = coef_names)
  )
}
