# Packages: MASS, gbm, glmnet, Matching, CBPS, sbw, quadprog.
# Set repetitions, case_ids and action below. Parallel workers: 100.
# action: "run" or "summarize". Initial seed: 123.
# Case   n     p     Outcome
# 1      200   100   Y1
# 2      200   100   Y2
# 3      500   200   Y1
# 4      500   200   Y2
# 5      500   100   Y1
# 6      500   100   Y2
# 7      1000  100   Y1
# 8      1000  100   Y2
# Covariate correlation: 0.8.
repetitions <- 500L

case_ids <- 1:8

action <- "run"

settings <- list(penalties = 10^seq(1, -6, length.out = 8), folds = 5L, trees = 300L, learning_rate = 0.1, depth = 3L, 
    maxit = 5000L, tolerance = 1e-04, ps_ridge = c(0.001, 0.01, 0.1, 1))

script_location <- function() {
    a <- grep("^--file=", commandArgs(), value = TRUE)
    if (length(a)) 
        return(normalizePath(sub("^--file=", "", a[1]), winslash = "/"))
    for (f in rev(sys.frames())) if (!is.null(f$ofile)) 
        return(normalizePath(f$ofile, winslash = "/"))
    normalizePath("01_simulation_and_summary.R", winslash = "/", mustWork = TRUE)
}

bind_rows <- function(x) {
    x <- Filter(function(z) !is.null(z) && nrow(z) > 0, x)
    if (!length(x)) 
        return(data.frame())
    do.call(rbind, x)
}

softmax <- function(x) {
    z <- exp(x - max(x))
    z/sum(z)
}

soft_threshold <- function(x, a) sign(x) * pmax(abs(x) - a, 0)

make_basis <- function(X, D) {
    X <- as.matrix(X)
    mu <- colMeans(X[D == 1, , drop = FALSE])
    B <- cbind(X, sweep(X, 2, mu, "-")^2)
    colnames(B) <- c(colnames(X), paste0(colnames(X), "_sq"))
    sc <- apply(X, 2, sd)
    sc[sc < 1e-10] <- 1
    Z <- sweep(sweep(X, 2, colMeans(X), "-"), 2, sc, "/")
    list(B = B, Z = Z, center = mu, scale = rep(1, ncol(B)))
}

entropy_fit <- function(B, D, l1, l2, cfg, init = NULL) {
    C <- sweep(B[D == 0, , drop = FALSE], 2, colMeans(B[D == 1, , drop = FALSE]), "-")
    theta <- if (is.null(init)) 
        rep(0, ncol(C))
    else init
    z <- theta
    accel <- 1
    L <- 1
    converged <- FALSE
    smooth <- function(v) {
        eta <- drop(C %*% v)
        mm <- max(eta)
        w <- softmax(eta)
        list(f = mm + log(mean(exp(eta - mm))), g = drop(crossprod(C, w)), w = w)
    }
    for (iter in seq_len(cfg$maxit)) {
        s <- smooth(z)
        for (bt in seq_len(60L)) {
            candidate <- soft_threshold(z - s$g/L, l1/L)/(1 + l2/L)
            delta <- candidate - z
            sn <- smooth(candidate)
            if (is.finite(sn$f) && sn$f <= s$f + sum(s$g * delta) + L * sum(delta^2)/2 + 1e-12) 
                break
            L <- L * 2
        }
        g <- sn$g + l2 * candidate
        kkt <- max(ifelse(abs(candidate) > 1e-09, abs(g + l1 * sign(candidate)), pmax(abs(g) - l1, 0)))
        if (kkt < cfg$tolerance) {
            theta <- candidate
            converged <- TRUE
            break
        }
        an <- (1 + sqrt(1 + 4 * accel^2))/2
        z <- candidate + ((accel - 1)/an) * (candidate - theta)
        theta <- candidate
        accel <- an
        if (iter%%100L == 0L) {
            z <- theta
            accel <- 1
        }
        L <- max(L/1.2, 1e-08)
    }
    final <- smooth(theta)
    list(w = final$w, theta = theta, converged = converged, iterations = iter, kkt = kkt, objective = final$f + 
        l1 * sum(abs(theta)) + l2 * sum(theta^2)/2)
}

weight_diagnostics <- function(B, D, w) {
    gap <- colMeans(B[D == 1, , drop = FALSE]) - drop(crossprod(w, B[D == 0, , drop = FALSE]))
    scale <- apply(B, 2, sd)
    scale[!is.finite(scale) | scale < 1e-10] <- 1
    gap <- gap/scale
    list(max_smd = max(abs(gap)), mean_smd = mean(abs(gap)), ess = 1/sum(w^2), max_weight = max(w), gap = gap)
}

tune_weights <- function(B, D, method, cfg) {
    g <- switch(method, EBL1 = data.frame(l1 = cfg$penalties, l2 = 0), EBL2 = data.frame(l1 = 0, l2 = cfg$penalties), 
        EBenet = expand.grid(l1 = cfg$penalties, l2 = cfg$penalties))
    ctrl <- which(D == 0)
    trt <- which(D == 1)
    p <- ncol(B)/2
    target <- colMeans(B[trt, , drop = FALSE])
    sd_t <- apply(B[trt, seq_len(p), drop = FALSE], 2, sd)
    sd_t[!is.finite(sd_t) | sd_t < 1e-10] <- 1
    fold <- cfg$balance_fold
    if (is.null(fold)) 
        fold <- sample(rep(seq_len(cfg$folds), length.out = length(ctrl)))
    tab <- g
    tab$cv_smd <- Inf
    tab$converged <- FALSE
    details <- vector("list", nrow(g))
    for (j in seq_len(nrow(g))) {
        losses <- rep(NA_real_, cfg$folds)
        for (k in seq_len(cfg$folds)) {
            train <- c(ctrl[fold != k], trt)
            val <- ctrl[fold == k]
            f <- entropy_fit(B[train, , drop = FALSE], D[train], g$l1[j], g$l2[j], cfg)
            if (!isTRUE(f$converged)) 
                next
            w <- softmax(drop(B[val, , drop = FALSE] %*% f$theta))
            gap <- (target[seq_len(p)] - drop(crossprod(w, B[val, seq_len(p), drop = FALSE])))/sd_t
            losses[k] <- mean(abs(gap))
        }
        details[[j]] <- losses
        if (all(is.finite(losses))) {
            tab$cv_smd[j] <- mean(losses)
            tab$converged[j] <- TRUE
        }
    }
    good <- which(tab$converged)
    if (!length(good)) 
        return(list(w = NULL, converged = FALSE, tuning = tab, message = "No candidate converged in all validation folds"))
    j <- good[which.min(tab$cv_smd[good])]
    f <- entropy_fit(B, D, g$l1[j], g$l2[j], cfg)
    tab$selected <- seq_len(nrow(tab)) == j
    f$tuning <- tab
    f$fold_smd <- details
    f
}

stratified_folds <- function(D, K) {
    fold <- integer(length(D))
    for (d in 0:1) {
        ii <- which(D == d)
        fold[ii] <- sample(rep(seq_len(K), length.out = length(ii)))
    }
    fold
}

crossfit_outcome <- function(X, Y, D, cfg) {
    fold <- stratified_folds(D, cfg$folds)
    pred <- rep(NA_real_, length(Y))
    for (k in seq_len(cfg$folds)) {
        tr <- which(D == 0 & fold != k)
        va <- which(fold == k)
        if (length(tr) < 5) 
            stop("Too few training controls")
        fit <- gbm::gbm(Y ~ ., data = data.frame(Y = Y[tr], X[tr, , drop = FALSE]), distribution = "gaussian", 
            n.trees = cfg$trees, interaction.depth = cfg$depth, shrinkage = cfg$learning_rate, verbose = FALSE)
        pred[va] <- predict(fit, newdata = as.data.frame(X[va, , drop = FALSE]), n.trees = cfg$trees)
    }
    list(pred = pred, fold = fold)
}

ridge_ps <- function(B, D, cfg, fold) {
    X <- cbind(1, B)
    fit_one <- function(ii, lam) {
        xx <- X[ii, , drop = FALSE]
        yy <- D[ii]
        n <- length(ii)
        loss <- function(b) {
            eta <- drop(xx %*% b)
            mean(pmax(eta, 0) + log1p(exp(-abs(eta))) - yy * eta) + lam * sum(b[-1]^2)/2
        }
        grad <- function(b) drop(crossprod(xx, plogis(xx %*% b) - yy))/n + c(0, lam * b[-1])
        f <- optim(rep(0, ncol(X)), loss, grad, method = "BFGS", control = list(maxit = 1000))
        if (f$convergence != 0) 
            stop("Ridge PS optimizer did not converge")
        f$par
    }
    score <- numeric(length(cfg$ps_ridge))
    for (j in seq_along(score)) for (k in seq_len(cfg$folds)) {
        b <- fit_one(which(fold != k), cfg$ps_ridge[j])
        ii <- which(fold == k)
        pr <- pmin(pmax(plogis(X[ii, , drop = FALSE] %*% b), 1e-08), 1 - 1e-08)
        score[j] <- score[j] - sum(D[ii] * log(pr) + (1 - D[ii]) * log1p(-pr))/length(D)
    }
    j <- which.min(score)
    b <- fit_one(seq_along(D), cfg$ps_ridge[j])
    ps <- pmin(pmax(plogis(drop(X %*% b)), 1e-08), 1 - 1e-08)
    list(w = {
        v <- ps[D == 0]/(1 - ps[D == 0])
        v/sum(v)
    }, ps = ps, theta = b, converged = TRUE, tuning = data.frame(lambda = cfg$ps_ridge, cv_loss = score))
}

original_ps <- function(X, D, design = 1, alpha = 1, lambda_choice = "lambda.min", nfolds = 5) {
    if (!is.matrix(X) && !is.data.frame(X)) {
        stop("X must be a matrix or data frame")
    }
    if (length(D) != nrow(X)) {
        stop("Length of D must equal nrow(X)")
    }
    if (!all(D %in% c(0, 1))) {
        stop("D must contain only 0 and 1")
    }
    X_df <- as.data.frame(X)
    n <- nrow(X_df)
    p <- ncol(X_df)
    if (p == 0) 
        stop("X must have at least one column")
    all_covs <- colnames(X_df)
    if (is.null(all_covs)) {
        all_covs <- paste0("X", 1:p)
        colnames(X_df) <- all_covs
    }
    stopifnot(design == 1)
    covs_to_use <- all_covs
    formula_str <- paste("~", paste(covs_to_use, collapse = " + "))
    model_mat <- model.matrix(as.formula(formula_str), data = X_df)[, -1, drop = FALSE]
    if (ncol(model_mat) == 0) {
        stop("Could not construct the covariate model matrix; check column names")
    }
    y <- as.factor(D)
    cv_fit <- glmnet::cv.glmnet(x = model_mat, y = y, family = "binomial", alpha = alpha, nfolds = nfolds, type.measure = "deviance", 
        standardize = TRUE, trace.it = 0)
    best_lambda <- ifelse(lambda_choice == "lambda.1se", cv_fit$lambda.1se, cv_fit$lambda.min)
    ps <- predict(cv_fit, newx = model_mat, s = best_lambda, type = "response")[, 1]
    list(model = cv_fit, ps = as.vector(ps), covs_used = covs_to_use, lambda_used = best_lambda, design = design, 
        alpha = alpha)
}

original_cbps <- function(X, Y, D, method = "over", alpha = 0.5, nfolds = 5, lambda_choice = "lambda.min", verbose = FALSE) {
    df <- data.frame(D = D, X)
    n <- nrow(df)
    p <- ncol(X)
    if (p == 0) 
        stop("X must have at least one column")
    if (verbose) 
        cat("Number of variables:", p, "\n")
    X_mat <- as.matrix(X)
    y_vec <- as.factor(D)
    cv_fit <- glmnet::cv.glmnet(x = X_mat, y = y_vec, family = "binomial", alpha = alpha, nfolds = nfolds, type.measure = "deviance", 
        standardize = TRUE)
    best_lambda <- ifelse(lambda_choice == "lambda.1se", cv_fit$lambda.1se, cv_fit$lambda.min)
    beta <- as.matrix(coef(cv_fit, s = best_lambda))
    sel_idx <- which(beta[-1, , drop = FALSE] != 0)
    sel_vars <- colnames(X_mat)[sel_idx]
    if (length(sel_vars) == 0) {
        warning("Elastic net selected no variables; using all variables.")
        sel_vars <- colnames(X_mat)
    }
    if (verbose) 
        cat("Number of selected variables:", length(sel_vars), "\n")
    rhs <- paste(sel_vars, collapse = " + ")
    form <- as.formula(paste("D ~ ", rhs))
    fit <- CBPS::CBPS(form, data = df, ATT = TRUE, method = method, twostep = TRUE)
    w <- fit$weights
    tau_hat <- with(df, weighted.mean(Y[D == 1], w[D == 1]) - weighted.mean(Y[D == 0], w[D == 0]))
    return(list(w = w[D == 0]/sum(w[D == 0]), converged = TRUE, estimate = tau_hat))
}

original_match <- function(Y, D, X, weight) {
    fit <- Matching::Match(Y = Y, Tr = D, X = X, M = 1, replace = TRUE, Weight = weight)
    v <- numeric(length(D))
    ww <- fit$weights
    if (is.null(ww)) 
        ww <- rep(1, length(fit$index.control))
    for (j in seq_along(fit$index.control)) v[fit$index.control[j]] <- v[fit$index.control[j]] + ww[j]
    list(w = v[D == 0]/sum(v[D == 0]), estimate = as.numeric(fit$est), converged = TRUE)
}

method_names <- function() c("Unadjusted", "IPW", "PSM", "MD", "CBPS_e", "CBPS_o", "EBL1", "EBL2", "EBenet", "RidgePS")

augmented_methods <- function() c("EBL1", "EBL2", "EBenet")

result_names <- function() c(method_names(), paste0(augmented_methods(), "_AUG"))

analyze_data <- function(dat, cfg, keep = TRUE) {
    X <- dat$X
    D <- dat$D
    Y <- dat$Y
    if (min(table(factor(D, levels = 0:1))) < 5) 
        stop("Fewer than 5 subjects in a group")
    basis <- make_basis(X, D)
    B <- basis$B
    cfg$balance_fold <- sample(rep(seq_len(cfg$folds), length.out = sum(D == 0)))
    cf <- crossfit_outcome(X, Y, D, cfg)
    n0 <- sum(D == 0)
    fits <- list()
    rows <- list()
    all_warnings <- character()
    for (method in method_names()) {
        started <- proc.time()[3]
        f <- tryCatch(withCallingHandlers({
            if (method == "Unadjusted") 
                list(w = rep(1/n0, n0), converged = TRUE)
            else if (method == "IPW") {
                psf <- original_ps(X, D, design = 1)
                ps <- pmin(pmax(psf$ps, 1e-08), 1 - 1e-08)
                v <- ps[D == 0]/(1 - ps[D == 0])
                list(w = v/sum(v), ps = ps, converged = TRUE, tuning = data.frame(lambda = psf$lambda_used))
            }
            else if (method == "PSM") {
                if (is.null(fits$IPW$ps)) 
                  stop("Original PS fit unavailable")
                original_match(Y, D, fits$IPW$ps, 1)
            }
            else if (method == "MD") 
                original_match(Y, D, X, 2)
            else if (method %in% c("CBPS_e", "CBPS_o")) 
                original_cbps(X, Y, D, method = if (method == "CBPS_e") 
                  "exact"
                else "over")
            else if (method == "RidgePS") 
                ridge_ps(B, D, cfg, cf$fold)
            else tune_weights(B, D, method, cfg)
        }, warning = function(w) {
            all_warnings <<- c(all_warnings, paste(method, conditionMessage(w)))
            invokeRestart("muffleWarning")
        }), error = function(e) list(w = NULL, converged = FALSE, message = conditionMessage(e)))
        elapsed <- unname(proc.time()[3] - started)
        valid <- isTRUE(f$converged) && length(f$w) == n0 && all(is.finite(f$w))
        if (!valid) 
            f$message <- if (is.null(f$message)) 
                "Optimization did not meet tolerance"
            else f$message
        diagnostic <- if (valid) 
            weight_diagnostics(B, D, f$w)
        else list(max_smd = NA_real_, mean_smd = NA_real_, ess = NA_real_, max_weight = NA_real_)
        variants <- if (method %in% augmented_methods()) 
            c(FALSE, TRUE)
        else FALSE
        for (aug in variants) {
            est <- if (!valid) 
                NA_real_
            else if (aug) 
                mean(Y[D == 1] - cf$pred[D == 1]) - sum(f$w * (Y[D == 0] - cf$pred[D == 0]))
            else if (!is.null(f$estimate)) 
                f$estimate
            else mean(Y[D == 1]) - sum(f$w * Y[D == 0])
            rows[[length(rows) + 1L]] <- data.frame(method = paste0(method, if (aug) 
                "_AUG"
            else ""), estimate = est, success = valid, max_smd = diagnostic$max_smd, mean_smd = diagnostic$mean_smd, 
                ess = diagnostic$ess, max_weight = diagnostic$max_weight, seconds = elapsed, message = if (valid) 
                  ""
                else f$message, stringsAsFactors = FALSE)
        }
        fits[[method]] <- f
    }
    list(results = bind_rows(rows), fits = if (keep) fits else NULL, outcome = if (keep) cf else NULL, basis = if (keep) basis[c("center", 
        "scale")] else NULL, warnings = unique(all_warnings))
}

gen_one_population <- function(N, p, rho) {
    X1 <- rbinom(N, 1, 0.3)
    X2 <- rbinom(N, 1, 0.2)
    n_cont <- p - 2
    Sigma <- matrix(rho, nrow = n_cont, ncol = n_cont)
    diag(Sigma) <- 1
    X_cont <- MASS::mvrnorm(n = N, mu = rep(0, n_cont), Sigma = Sigma)
    colnames(X_cont) <- paste0("X", 3:p)
    X <- cbind(X1, X2, X_cont)
    base_terms <- X
    X1 <- X[, "X1"]
    X2 <- X[, "X2"]
    X3 <- X[, "X3"]
    X4 <- X[, "X4"]
    X5 <- X[, "X5"]
    X6 <- X[, "X6"]
    higher_terms <- cbind(X3^2, X4^2, X5^2, X6^2, X3 * X4, X3 * X1, X4 * X5, X4 * X6, X4 * X2, X5 * X2)
    colnames(higher_terms) <- c("X3^2", "X4^2", "X5^2", "X6^2", "X3X4", "X3X1", "X4X5", "X4X6", "X4X2", "X5X2")
    selected_higher <- higher_terms
    beta_base <- rnorm(ncol(base_terms), mean = 0, sd = 1)
    beta_higher <- rnorm(ncol(selected_higher), mean = 0, sd = 1)
    f_X <- as.numeric(base_terms %*% beta_base) + as.numeric(selected_higher %*% beta_higher)
    eps <- rnorm(N, mean = 0, sd = sqrt(8))
    lin <- f_X - 2 + eps
    D <- as.integer(lin > 0)
    u <- rnorm(N, mean = 0, sd = 1)
    y10 <- X3 + X4 + X5 - X6 + X1 + X2 + u
    y11 <- y10 + 2
    Y1 <- ifelse(D == 1, y11, y10)
    y20 <- X3 + X3^2 - X6 * X2 + u
    y21 <- y20 + 2
    Y2 <- ifelse(D == 1, y21, y20)
    list(X = X, D = D, Y = cbind(Y1 = Y1, Y2 = Y2))
}

draw_by_ratio_adaptive <- function(N, ratio_ct, p, rho, max_iter = 50) {
    n1_needed <- round(N/(1 + ratio_ct))
    n0_needed <- N - n1_needed
    batch_size <- N
    iter <- 0
    while (iter < max_iter) {
        iter <- iter + 1
        X_acc <- NULL
        D_acc <- integer(0)
        Y_acc <- NULL
        pop <- gen_one_population(N = batch_size, p = p, rho = rho)
        Xb <- pop$X
        Db <- pop$D
        Yb <- pop$Y
        X_acc <- Xb
        D_acc <- Db
        Y_acc <- Yb
        if (sum(D_acc == 1) >= n1_needed && sum(D_acc == 0) >= n0_needed) {
            idx1 <- which(D_acc == 1)
            idx0 <- which(D_acc == 0)
            pick1 <- sample(idx1, n1_needed)
            pick0 <- sample(idx0, n0_needed)
            idx <- c(pick1, pick0)
            return(list(X = X_acc[idx, , drop = FALSE], D = D_acc[idx], Y = Y_acc[idx, , drop = FALSE]))
        }
        batch_size <- ceiling(batch_size * 1.1)
    }
    stop("Insufficient treated or control observations after maximum sampling attempts")
}

generate_data <- function(s) {
    samp <- draw_by_ratio_adaptive(N = s$n, ratio_ct = s$ratio_ct, p = s$p, rho = s$rho, 
        max_iter = 50)
    list(X = samp$X, D = samp$D, Y = as.numeric(samp$Y[, s$outcome_id]), true_att = 2)
}

simulation_job <- function(task, cfg) {
    assign(".Random.seed", task$rng, envir = .GlobalEnv)
    dat <- generate_data(task$s)
    ans <- tryCatch(analyze_data(dat, cfg), error = function(e) list(results = data.frame(method = result_names(), 
        estimate = NA_real_, success = FALSE, max_smd = NA_real_, mean_smd = NA_real_, ess = NA_real_, max_weight = NA_real_, 
        seconds = NA_real_, message = conditionMessage(e))))
    ans$results$truth <- dat$true_att
    ans$data <- dat[c("D", "Y", "true_att")]
    if (task$rep == 1L) 
        ans$data$X <- dat$X
    ans$scenario <- task$s$scenario
    ans$rep <- task$rep
    ans$rng <- task$rng
    ans
}

save_checkpoint <- function(x, path) {
    tmp <- paste0(path, ".tmp")
    saveRDS(x, tmp, compress = FALSE)
    if (!file.copy(tmp, path, overwrite = TRUE)) 
        stop("Could not save checkpoint")
    unlink(tmp)
}

sbw_config <- list(grid = c(1e-04, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1), tuning_resamples = 1000L, 
    solver = "quadprog", scale = "group")

fit_sbw <- function(X, D) {
    B <- make_basis(X, D)$B
    ctrl <- which(D == 0)
    trt <- which(D == 1)
    target <- colMeans(B[trt, , drop = FALSE])
    sd0 <- apply(B[ctrl, , drop = FALSE], 2, sd)
    constant <- !is.finite(sd0) | sd0 < 1e-12
    gap0 <- colMeans(B[ctrl, , drop = FALSE]) - target
    if (any(constant & abs(gap0) > 1e-08)) 
        stop("SBW infeasible: a control-constant feature differs from the treated target")
    keep <- !constant
    if (!any(keep)) 
        stop("No nonconstant balancing features")
    dd <- data.frame(row_id = seq_len(nrow(X)), D = D, B[, keep, drop = FALSE], check.names = FALSE)
    covariates <- colnames(B)[keep]
    feasible <- logical(length(sbw_config$grid))
    reason <- rep("", length(feasible))
    for (j in seq_along(feasible)) {
        trial <- tryCatch(sbw::sbw(dat = dd, ind = "D", out = NULL, bal = list(bal_cov = covariates, bal_alg = FALSE, 
            bal_std = sbw_config$scale, bal_tol = sbw_config$grid[j]), wei = list(wei_sum = TRUE, wei_pos = TRUE), 
            sol = list(sol_nam = sbw_config$solver, sol_dis = FALSE), par = list(par_est = "att", par_tar = NULL), 
            mes = FALSE), error = function(e) e)
        feasible[j] <- !inherits(trial, "error")
        if (!feasible[j]) 
            reason[j] <- conditionMessage(trial)
    }
    feasibility <- data.frame(tolerance = sbw_config$grid, feasible = feasible, message = reason)
    if (!any(feasible)) 
        return(list(error = "All prespecified SBW tolerances are infeasible", feasibility = feasibility))
    f <- sbw::sbw(dat = dd, ind = "D", out = NULL, bal = list(bal_cov = covariates, bal_alg = TRUE, bal_std = sbw_config$scale, 
        bal_gri = sbw_config$grid[feasible], bal_sam = sbw_config$tuning_resamples), wei = list(wei_sum = TRUE, 
        wei_pos = TRUE), sol = list(sol_nam = sbw_config$solver, sol_dis = FALSE), par = list(par_est = "att", 
        par_tar = NULL), mes = FALSE)
    output <- f$dat_weights
    stopifnot(all(c("row_id", "sbw_weights") %in% names(output)))
    w <- as.numeric(output$sbw_weights[match(ctrl, output$row_id)])
    if (length(w) != length(ctrl) || any(!is.finite(w)) || any(w < -1e-07) || sum(w) <= 0) 
        stop("SBW returned invalid control weights")
    w <- pmax(w, 0)
    w <- w/sum(w)
    delta <- as.numeric(f$bal$bal_tol)
    if (length(delta) != 1L || !is.finite(delta)) 
        stop("Could not identify selected SBW tolerance")
    gap <- drop(crossprod(w, B[ctrl, , drop = FALSE])) - target
    if (any(abs(gap) > delta * sd0 + 1e-06 * (1 + sd0 + abs(target)))) 
        stop("SBW solution violates the selected balance constraints")
    list(w = w, selected_tolerance = delta, cstat = f$cstat, feasibility = feasibility, balance_parameters = f$balance_parameters, 
        status = f$status, diagnostics = weight_diagnostics(B, D, w))
}

sbw_job <- function(task) {
    assign(".Random.seed", task$rng, envir = .GlobalEnv)
    started <- proc.time()[3]
    warnings <- character()
    dat <- NULL
    f <- tryCatch(withCallingHandlers({
        dat <- generate_data(task$s)
        fit_sbw(dat$X, dat$D)
    }, warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
    }), error = function(e) list(error = conditionMessage(e)))
    ok <- is.null(f$error)
    diagnostic <- if (ok) 
        f$diagnostics
    else list(max_smd = NA_real_, mean_smd = NA_real_, ess = NA_real_, max_weight = NA_real_)
    est <- if (ok) 
        mean(dat$Y[dat$D == 1]) - sum(f$w * dat$Y[dat$D == 0])
    else NA_real_
    row <- data.frame(method = "SBW", estimate = est, success = ok, max_smd = diagnostic$max_smd, mean_smd = diagnostic$mean_smd, 
        ess = diagnostic$ess, max_weight = diagnostic$max_weight, seconds = unname(proc.time()[3] - started), message = if (ok) 
            ""
        else f$error, truth = 2, selected_tolerance = if (ok) 
            f$selected_tolerance
        else NA_real_, seed = 123L, scenario = task$s$scenario, rep = task$rep, stringsAsFactors = FALSE)
    saved_data <- if (is.null(dat)) 
        NULL
    else dat[c("D", "Y", "true_att")]
    if (!is.null(dat) && task$rep == 1L) 
        saved_data$X <- dat$X
    list(results = row, w = if (ok) f$w else NULL, tuning = if (ok) f[c("selected_tolerance", "cstat", "balance_parameters", 
        "status", "feasibility")] else f["feasibility"], data = saved_data, warnings = unique(warnings), rng = task$rng, 
        scenario = task$s$scenario, rep = task$rep)
}

script_file <- script_location()

scenarios <- data.frame(scenario = 1:8, n = c(200L, 200L, 500L, 500L, 500L, 500L, 1000L, 1000L), p = c(100L, 100L, 
    200L, 200L, 100L, 100L, 100L, 100L), ratio_ct = 1, outcome_id = rep(1:2, 4), rho = 0.8, 
    design = rep(c("Y1", "Y2"), 4))

stream_slots <- c(1L, 2L, 13L, 14L, 10L, 11L, 22L, 23L)

output_dir <- file.path(dirname(script_file), "simulation_results")

if (dir.exists(file.path(dirname(script_file), ".Rlib"))) .libPaths(c(file.path(dirname(script_file), ".Rlib"), 
    .libPaths()))

case_tasks <- function(id) {
    set.seed(123, kind = "L'Ecuyer-CMRG")
    stream <- .Random.seed
    for (j in seq_len(stream_slots[id] - 1L)) stream <- parallel::nextRNGStream(stream)
    out <- vector("list", repetitions)
    for (r in seq_len(repetitions)) {
        out[[r]] <- list(key = paste(id, r, sep = "_"), rep = r, s = scenarios[id, , drop = FALSE], rng = stream)
        for (j in seq_len(24L)) stream <- parallel::nextRNGStream(stream)
    }
    out
}

combined_job <- function(task) {
    primary <- simulation_job(task, settings)
    separate <- sbw_job(task)
    stopifnot(identical(primary$data$D, separate$data$D), identical(primary$data$Y, separate$data$Y))
    if (!is.null(primary$data$X)) 
        stopifnot(identical(primary$data$X, separate$data$X))
    primary$sbw <- separate
    primary
}

summarize_results <- function() {
    files <- list.files(output_dir, "^case_[0-9]+\\.rds$", full.names = TRUE)
    if (!length(files))
        stop("No case files in ", output_dir)
    states <- lapply(files, readRDS)
    ids <- vapply(states, function(st) as.integer(st$case_id), integer(1))
    if (anyDuplicated(ids)) stop("Duplicate case files")
    if (!all(vapply(states, function(st) identical(st$config, states[[1]]$config), logical(1))))
        stop("Case files have different settings")
    if (any(!vapply(states, function(st) isTRUE(st$completed), logical(1))))
        warning("Some cases are incomplete; summary uses available repetitions")
    rows <- list()
    for (st in states) {
        if (!isTRUE(all.equal(st$scenario, scenarios[st$case_id, , drop = FALSE], check.attributes = FALSE))) 
            stop("Scenario metadata differ; do not merge by case number alone")
        for (z in st$simulation) {
            a <- z$results
            b <- z$sbw$results
            a <- a[, intersect(names(a), names(b)), drop = FALSE]
            b <- b[, names(a), drop = FALSE]
            v <- rbind(a, b)
            v$rep <- z$rep
            v$scenario <- st$case_id
            rows[[length(rows) + 1L]] <- v
        }
    }
    raw <- do.call(rbind, rows)
    if (anyDuplicated(paste(raw$scenario, raw$rep, raw$method))) 
        stop("Duplicate repetitions")
    avg <- function(x) if (any(is.finite(x))) 
        mean(x[is.finite(x)])
    else NA_real_
    sd_finite <- function(x) {
        x <- x[is.finite(x)]
        if (length(x) > 1L) sd(x) else NA_real_
    }
    mcse <- function(x) sd_finite(x) / sqrt(sum(is.finite(x)))
    out <- lapply(split(raw, interaction(raw$scenario, raw$method, drop = TRUE)), function(v) {
        ok <- v$success & is.finite(v$estimate)
        e <- v$estimate[ok] - v$truth[ok]
        info <- scenarios[scenarios$scenario == v$scenario[1],
            c("scenario", "n", "p", "ratio_ct", "design"), drop = FALSE]
        cbind(info, data.frame(method = v$method[1], attempted = nrow(v),
            successful = sum(ok), failure_rate = mean(!ok),
            mcse_failure = sqrt(mean(!ok) * mean(ok) / nrow(v)),
            bias = avg(e), mcse_bias = mcse(e), mse = avg(e^2),
            rmse = sqrt(avg(e^2)),
            mcse_rmse = if (is.finite(avg(e^2)) && avg(e^2) > 0)
                mcse(e^2) / (2 * sqrt(avg(e^2))) else NA_real_,
            mae = avg(abs(e)), mcse_mae = mcse(abs(e)),
            empirical_sd = sd_finite(v$estimate[ok]),
            mean_ess = avg(v$ess[ok]), mean_max_smd = avg(v$max_smd[ok]),
            mean_max_weight = avg(v$max_weight[ok])))
    })
    out <- do.call(rbind, out)
    out <- out[order(out$scenario, out$method), ]
    write.csv(out, file.path(output_dir, "summary.csv"), row.names = FALSE, na = "NA")
    message("Saved summary.csv")
    invisible(out)
}

run_simulations <- function() {
    stopifnot(length(repetitions) == 1L, is.finite(repetitions), repetitions >= 1, repetitions == as.integer(repetitions), 
        length(case_ids) > 0, all(case_ids %in% 1:8))
    needed <- c("MASS", "gbm", "glmnet", "Matching", "CBPS", "sbw")
    missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing)) 
        stop("Install missing packages: ", paste(missing, collapse = ", "))
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    ids <- sort(unique(as.integer(case_ids)))
    print(scenarios[ids, ], row.names = FALSE)
    config <- list(scenarios = scenarios, settings = settings,
                   sbw = sbw_config, streams = stream_slots)
    cl <- parallel::makePSOCKcluster(100L)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterCall(cl, function(paths) .libPaths(paths), .libPaths())
    parallel::clusterExport(cl, c("script_location", "bind_rows", "softmax", "soft_threshold", "make_basis", "entropy_fit", "weight_diagnostics", "tune_weights", "stratified_folds", "crossfit_outcome", "ridge_ps", "original_ps", "original_cbps", "original_match", "method_names", "augmented_methods", "result_names", "analyze_data", "gen_one_population", "draw_by_ratio_adaptive", "generate_data", "simulation_job", "save_checkpoint", "fit_sbw", "sbw_job", "case_tasks", "combined_job", "summarize_results", "run_simulations", "settings", "sbw_config"), envir = environment(run_simulations))
    for (id in ids) {
        path <- file.path(output_dir, sprintf("case_%02d.rds", id))
        if (file.exists(path)) {
            state <- readRDS(path)
            if (!identical(state$config, config)) 
                stop("Existing file differs: ", path)
        }
        else state <- list(case_id = id, scenario = scenarios[id, , drop = FALSE], seed = 123L, requested_sim = repetitions, 
            simulation = list(), completed = FALSE, config = config, session = sessionInfo())
        tasks <- Filter(function(t) is.null(state$simulation[[t$key]]), case_tasks(id))
        if (length(tasks)) 
            for (wave in split(tasks, ceiling(seq_along(tasks)/100L))) {
                ans <- parallel::parLapplyLB(cl, wave, combined_job)
                for (j in seq_along(wave)) state$simulation[[wave[[j]]$key]] <- ans[[j]]
                save_checkpoint(state, path)
                message("Case ", id, ": ", length(state$simulation), "/", repetitions)
            }
        state$completed <- length(state$simulation) == repetitions
        save_checkpoint(state, path)
    }
    summarize_results()
}

if (action == "run") {
    run_simulations()
} else if (action == "summarize") {
    summarize_results()
} else stop('action must be "run" or "summarize"')
