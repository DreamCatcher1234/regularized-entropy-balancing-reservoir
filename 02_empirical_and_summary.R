# Packages: gbm, glmnet, Matching, CBPS, sbw, quadprog.
# Set repetitions and action below. Parallel workers: 100.
# action: "bootstrap" or "summarize". Bootstrap seeds: 123, 124, ...
# Place camels_br_att_data.rds in this script's directory.
repetitions <- 500L
action <- "bootstrap"

script_location <- function() {
    args <- grep("^--file=", commandArgs(), value = TRUE)
    if (length(args)) 
        return(normalizePath(sub("^--file=", "", args[1]), winslash = "/"))
    for (f in rev(sys.frames())) {
        if (!is.null(f$ofile)) 
            return(normalizePath(f$ofile, winslash = "/"))
    }
    normalizePath("02_empirical_and_summary.R", winslash = "/", mustWork = TRUE)
}

script_dir <- dirname(script_location())
data_file <- file.path(script_dir, "camels_br_att_data.rds")
bootstrap_file <- file.path(script_dir, "camels_att_bootstrap_recomputed.rds")
settings <- list(penalties = 10^seq(1, -6, length.out = 8), folds = 5L, trees = 300L, learning_rate = 0.1, depth = 3L, 
    maxit = 5000L, tolerance = 1e-04, ps_ridge = c(0.001, 0.01, 0.1, 1))

sbw_config <- list(grid = c(1e-04, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1), tuning_resamples = 1000L, 
    solver = "quadprog", scale = "manual")

att_methods <- c("Unadjusted", "IPW", "PSM", "MD", "CBPS_e", "CBPS_o", "RidgePS", "SBW", "EBL1", "EBL1_AUG", "EBL2", 
    "EBL2_AUG", "EBenet", "EBenet_AUG")

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
    cv_fit <- glmnet::cv.glmnet(x = model_mat, y = y, family = "binomial", alpha = alpha, nfolds = nfolds, foldid = .att_fold, 
        type.measure = "deviance", standardize = TRUE, trace.it = 0)
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
    cv_fit <- glmnet::cv.glmnet(x = X_mat, y = y_vec, family = "binomial", alpha = alpha, nfolds = nfolds, foldid = .att_fold, 
        type.measure = "deviance", standardize = TRUE)
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

fit_sbw <- function(X, D) {
    B <- make_basis(X, D)$B
    ctrl <- which(D == 0)
    trt <- which(D == 1)
    target <- colMeans(B[trt, , drop = FALSE])
    sd0 <- apply(B, 2, sd)
    constant <- !is.finite(sd0) | sd0 < 1e-12
    keep <- !constant
    if (!any(keep)) 
        stop("No nonconstant balancing features")
    scaled_B <- sweep(B[, keep, drop = FALSE], 2, sd0[keep], "/")
    dd <- data.frame(row_id = seq_len(nrow(X)), D = D, scaled_B, check.names = FALSE)
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
        status = f$status, tolerance_scale = "Full-sample SD of each balancing feature; official manual mode", 
        diagnostics = weight_diagnostics(B, D, w))
}

grouped_folds <- function(D, id, K) {
    fold <- integer(length(D))
    for (d in 0:1) {
        ids <- unique(id[D == d])
        if (length(ids) < K) 
            stop("Too few unique catchments for validation folds")
        assigned <- sample(rep(seq_len(K), length.out = length(ids)))
        fold[D == d] <- assigned[match(id[D == d], ids)]
    }
    fold
}

stratified_folds <- function(D, K) .att_fold

fit_att_sample <- function(index, job_seed, keep_details = FALSE) {
    set.seed(job_seed)
    dat <- att_data
    X <- dat$X[index, , drop = FALSE]
    D <- dat$D[index]
    Y <- dat$Y[index]
    id <- dat$id[index]
    center <- colMeans(X)
    scale <- apply(X, 2, sd)
    scale[!is.finite(scale) | scale < 1e-12] <- 1
    X <- sweep(sweep(X, 2, center, "-"), 2, scale, "/")
    .att_fold <<- grouped_folds(D, id, settings$folds)
    cfg <- settings
    cfg$balance_fold <- .att_fold[D == 0]
    basis <- make_basis(X, D)
    B <- basis$B
    n0 <- sum(D == 0)
    outcome_warnings <- character()
    cf <- tryCatch(withCallingHandlers(crossfit_outcome(X, Y, D, cfg), warning = function(w) {
        outcome_warnings <<- c(outcome_warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
    }), error = function(e) list(error = conditionMessage(e)))
    outcome_ok <- is.null(cf$error) && length(cf$pred) == length(Y) && all(is.finite(cf$pred))
    fits <- list()
    rows <- list()
    for (method in c("Unadjusted", "IPW", "PSM", "MD", "CBPS_e", "CBPS_o", "RidgePS", "SBW", "EBL1", "EBL2", "EBenet")) {
        if (keep_details) 
            message("Original sample: ", method)
        started <- proc.time()[3]
        warnings <- character()
        f <- tryCatch(withCallingHandlers({
            if (method == "Unadjusted") {
                list(w = rep(1/n0, n0), converged = TRUE)
            }
            else if (method == "IPW") {
                psf <- original_ps(X, D, design = 1)
                ps <- pmin(pmax(psf$ps, 1e-08), 1 - 1e-08)
                v <- ps[D == 0]/(1 - ps[D == 0])
                list(w = v/sum(v), ps = ps, converged = TRUE, tuning = data.frame(lambda = psf$lambda_used))
            }
            else if (method == "PSM") {
                if (is.null(fits$IPW$ps)) 
                  stop("IPW propensity fit unavailable")
                original_match(Y, D, fits$IPW$ps, 1)
            }
            else if (method == "MD") {
                original_match(Y, D, X, 2)
            }
            else if (method %in% c("CBPS_e", "CBPS_o")) {
                original_cbps(X, Y, D, method = if (method == "CBPS_e") 
                  "exact"
                else "over")
            }
            else if (method == "RidgePS") {
                ridge_ps(B, D, cfg, .att_fold)
            }
            else if (method == "SBW") {
                z <- fit_sbw(X, D)
                if (!is.null(z$error)) 
                  stop(z$error)
                z$converged <- TRUE
                z
            }
            else {
                tune_weights(B, D, method, cfg)
            }
        }, warning = function(w) {
            warnings <<- c(warnings, conditionMessage(w))
            invokeRestart("muffleWarning")
        }), error = function(e) list(converged = FALSE, message = conditionMessage(e)))
        valid <- isTRUE(f$converged) && length(f$w) == n0 && all(is.finite(f$w)) && all(f$w >= -1e-08) && abs(sum(f$w) - 
            1) < 1e-06
        message_text <- if (valid) 
            ""
        else if (!is.null(f$message)) 
            f$message
        else "Optimization did not converge or returned invalid weights"
        diagnostics <- if (valid) 
            weight_diagnostics(B, D, f$w)
        else list(ess = NA_real_, max_smd = NA_real_, max_weight = NA_real_)
        variants <- if (method %in% c("EBL1", "EBL2", "EBenet")) 
            c(FALSE, TRUE)
        else FALSE
        for (aug in variants) {
            ok <- valid && (!aug || outcome_ok)
            estimate <- if (!ok) 
                NA_real_
            else if (aug) {
                mean(Y[D == 1] - cf$pred[D == 1]) - sum(f$w * (Y[D == 0] - cf$pred[D == 0]))
            }
            else if (!is.null(f$estimate)) 
                f$estimate
            else {
                mean(Y[D == 1]) - sum(f$w * Y[D == 0])
            }
            ok <- ok && is.finite(estimate)
            if (!ok) 
                estimate <- NA_real_
            msg <- if (valid && aug && !outcome_ok) 
                paste("Outcome model:", cf$error)
            else message_text
            rows[[length(rows) + 1L]] <- data.frame(method = paste0(method, if (aug) 
                "_AUG"
            else ""), estimate = estimate, success = ok, ess = diagnostics$ess, max_smd = diagnostics$max_smd, 
                max_weight = diagnostics$max_weight, seconds = unname(proc.time()[3] - started), message = msg, 
                warnings = paste(unique(c(warnings, if (aug) outcome_warnings)), collapse = " | "), stringsAsFactors = FALSE)
        }
        fits[[method]] <- f
    }
    rows <- do.call(rbind, rows)
    rows <- rows[match(att_methods, rows$method), , drop = FALSE]
    rownames(rows) <- NULL
    list(results = rows, details = if (keep_details) list(fits = fits, outcome = cf, X = X, D = D, Y = Y, center = center, 
        scale = scale, folds = .att_fold, basis_center = basis$center) else NULL)
}

bootstrap_job <- function(b) {
    tryCatch({
        z <- fit_att_sample(bootstrap_index[[b]], bootstrap_seeds[b])
        z$results$replicate <- b
        z$results
    }, error = function(e) {
        data.frame(method = att_methods, estimate = NA_real_, success = FALSE, ess = NA_real_, max_smd = NA_real_, 
            max_weight = NA_real_, seconds = NA_real_, message = conditionMessage(e), warnings = "", replicate = b, 
            stringsAsFactors = FALSE)
    })
}

summarize_att <- function(original, boot) {
    rows <- lapply(att_methods, function(method) {
        point <- original$results[original$results$method == method, ]
        v <- boot$estimate[boot$method == method & boot$success]
        v <- v[is.finite(v)]
        n_ok <- length(v)
        interval <- if (n_ok >= 2L) 
            unname(quantile(v, c(0.025, 0.975), type = 7))
        else c(NA_real_, NA_real_)
        note <- c("Individual-catchment iid bootstrap; spatial dependence not modeled", "Fixed completed data; imputation uncertainty not included", 
            "Causal interpretation requires timing and confounding assumptions")
        if (method %in% c("PSM", "MD")) {
            note <- c(note, "Nearest-neighbor ordinary-bootstrap interval is descriptive; 95% coverage not established")
        }
        if (n_ok < repetitions) 
            note <- c(note, "Summary conditional on successful fits; inspect failures")
        data.frame(method = method, att_original = point$estimate, att_boot_mean = if (n_ok) 
            mean(v)
        else NA_real_, att_boot_sd = if (n_ok >= 2L) 
            sd(v)
        else NA_real_, ci95_lower = interval[1], ci95_upper = interval[2], bootstrap_requested = repetitions, bootstrap_success = n_ok, 
            bootstrap_failed = repetitions - n_ok, failure_rate = (repetitions - n_ok)/repetitions, original_success = point$success, 
            original_message = point$message, ci_method = "Bootstrap percentile (2.5%, 97.5%)", notes = paste(note, 
                collapse = "; "), stringsAsFactors = FALSE)
    })
    do.call(rbind, rows)
}

main_att <- function() {
    stopifnot(length(repetitions) == 1L, is.finite(repetitions),
              repetitions >= 1L, repetitions == as.integer(repetitions))
    if (file.exists(file.path(script_dir, "camels_att_bootstrap_recomputed.rds")))
        stop("Output RDS already exists. Preserve it elsewhere before a new bootstrap run, or use action = summarize.")
    needed <- c("gbm", "glmnet", "Matching", "CBPS", "sbw", "quadprog")
    missing <- needed[!vapply(needed, requireNamespace, quietly = TRUE, FUN.VALUE = TRUE)]
    if (length(missing)) {
        stop("Run install.packages(c(\"", paste(missing, collapse = "\", \""), "\")) first.")
    }
    input_file <- data_file
    if (!file.exists(input_file)) 
        stop("Missing data file: ", input_file)
    z <- readRDS(input_file)
    att_data <<- list(X = as.matrix(z$X), D = z$D, Y = z$Y, id = z$id)
    n <- nrow(att_data$X)
    stopifnot(length(z$D) == n, length(z$Y) == n, length(z$id) == n, !anyDuplicated(z$id), all(is.finite(att_data$X)), 
        all(is.finite(z$Y)), all(z$D %in% 0:1), length(unique(z$D)) == 2L)
    bootstrap_seeds <<- 123L + seq_len(repetitions) - 1L
    bootstrap_index <<- lapply(bootstrap_seeds, function(seed) {
        set.seed(seed)
        sample.int(n, n, replace = TRUE)
    })
    message("Fitting original sample...")
    original <- fit_att_sample(seq_len(n), 123L, keep_details = TRUE)
    state <- list(seed = 123L, bootstrap_repetitions = repetitions, data = att_data, source = z$source, definitions = z$definitions, 
        causal_design = z$causal_design, preprocessing = z$preprocessing, settings = settings, sbw_config = sbw_config, 
        empirical_adaptations = c("X standardized again within each resample; all 66 input columns retained", "SBW tolerance in full-sample SD units using official manual mode", 
            "Validation folds group all copies of the same original catchment together"), bootstrap_scheme = "Ordinary iid resampling of 897 catchments with replacement", 
        bootstrap_index = bootstrap_index, bootstrap_seeds = bootstrap_seeds, original = original, bootstrap = NULL, 
        summary = NULL, complete = FALSE, package_versions = vapply(needed, function(p) as.character(packageVersion(p)), 
            character(1)), session_info = sessionInfo())
    rds_file <- file.path(script_dir, "camels_att_bootstrap_recomputed.rds")
    saveRDS(state, rds_file)
    workers <- 100L
    message("Starting ", repetitions, " bootstrap repetitions on ", workers, " workers...")
    thread_vars <- c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS")
    old_threads <- Sys.getenv(thread_vars, unset = NA_character_)
    do.call(Sys.setenv, as.list(setNames(rep("1", length(thread_vars)), thread_vars)))
    on.exit({
        for (nm in thread_vars) {
            if (is.na(old_threads[nm])) Sys.unsetenv(nm) else do.call(Sys.setenv, setNames(list(old_threads[nm]), 
                nm))
        }
    }, add = TRUE)
    cl <- parallel::makeCluster(workers, outfile = if (.Platform$OS.type == "windows") 
        "nul:"
    else "/dev/null")
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterCall(cl, function(paths) .libPaths(paths), .libPaths())
    functions <- c("softmax", "soft_threshold", "make_basis", "entropy_fit", "weight_diagnostics", "tune_weights", 
        "stratified_folds", "crossfit_outcome", "ridge_ps", "original_ps", "original_cbps", "original_match", "fit_sbw", 
        "grouped_folds", "fit_att_sample", "bootstrap_job")
    parallel::clusterExport(cl, c(functions, "att_data", "att_methods", "settings", "sbw_config", "bootstrap_index", 
        "bootstrap_seeds", "repetitions"), envir = .GlobalEnv)
    answers <- list()
    for (start in seq.int(1L, repetitions, by = workers)) {
        ids <- seq.int(start, min(repetitions, start + workers - 1L))
        results <- parallel::parLapplyLB(cl, ids, bootstrap_job)
        answers <- c(answers, results)
        state$bootstrap <- do.call(rbind, answers)
        state$completed_repetitions <- length(answers)
        saveRDS(state, rds_file)
        message("Bootstrap completed: ", length(answers), "/", repetitions)
    }
    state$summary <- summarize_att(original, state$bootstrap)
    state$complete <- TRUE
    saveRDS(state, rds_file)
    csv_file <- file.path(script_dir, "camels_att_summary.csv")
    write.csv(state$summary, csv_file, row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
    print(state$summary[, c("method", "att_original", "att_boot_mean", "att_boot_sd", "ci95_lower", "ci95_upper", 
        "bootstrap_success")], row.names = FALSE)
    message("Saved: ", csv_file, "\nSaved: ", rds_file)
    invisible(state$summary)
}

summarize_saved <- function() {
    if (!file.exists(bootstrap_file)) stop("Missing saved bootstrap file: ", bootstrap_file)
    state <- readRDS(bootstrap_file)
    repetitions <<- state$bootstrap_repetitions
    if (is.null(repetitions) || repetitions < 1L) stop("Missing bootstrap repetition metadata")
    original <- state$original
    original$results <- original$results[original$results$method %in% att_methods, , drop = FALSE]
    boot <- state$bootstrap[state$bootstrap$method %in% att_methods, , drop = FALSE]
    out <- summarize_att(original, boot)
    write.csv(out, file.path(script_dir, "camels_att_summary.csv"), row.names = FALSE, na = "NA")
    invisible(out)
}

if (action == "bootstrap") {
    main_att()
} else if (action == "summarize") {
    summarize_saved()
} else stop('action must be "bootstrap" or "summarize"')
