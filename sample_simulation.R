# Simulation sample code
# Required packages (install once if needed):
# install.packages(c("MASS", "ncvreg", "L0Learn"))

# ----- Settings -----
scenarios <- c("Scenario I")
n_list <- c(250)
n_reps <- 2L
seed <- 20262027L
save_dir <- "simulation_output"

# Full simulation: use the following settings instead.
# scenarios <- c("Scenario I", "Scenario II", "Scenario III")
# n_list <- c(250, 500, 1000, 2500)
# n_reps <- 1000L

packages <- c("MASS", "ncvreg", "L0Learn")
missing <- packages[!vapply(packages, requireNamespace, logical(1L), quietly = TRUE)]
if (length(missing)) stop("Install missing packages first: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(ncvreg))

# ----- Simulation functions -----
get_scenario_config <- function(scenario = c("Scenario I", "Scenario II", "Scenario III"), rho = 0.5) {
    scenario <- match.arg(scenario)
    if (scenario == "Scenario I") {
        p <- q <- k <- m <- 10
        beta_true <- c(2, 1.5, 1, rep(0, p - 3))
        gamma_true <- c(1.2, 1.4, -1, rep(0, q - 3))
        eta_true <- c(0.8, 1, -0.5, rep(0, m - 3))
        tau_true <- c(3, -1.5, 2.7, rep(0, k - 3))
        alpha0_true <- 1
        beta0_true <- 1.5
        gamma0_true <- 0.9
        tau0_true <- -1
    }
    if (scenario == "Scenario II") {
        p <- q <- k <- m <- 30
        beta_true <- c(2, 1.5, 1, rep(0, p - 3))
        gamma_true <- c(1.2, 1.4, -1, rep(0, q - 3))
        eta_true <- c(0.8, 1, -0.5, rep(0, m - 3))
        tau_true <- c(3, -1.5, 2.7, rep(0, k - 3))
        alpha0_true <- 1
        beta0_true <- 1.5
        gamma0_true <- 0.9
        tau0_true <- -1
    }
    if (scenario == "Scenario III") {
        p <- q <- k <- m <- 10
        beta_true <- c(2, 1.5, 1, rep(0, p - 3))
        gamma_true <- c(0.5, 0.5, -0.5, rep(0, q - 3))
        eta_true <- c(0.8, 1, -0.5, rep(0, m - 3))
        tau_true <- c(3, -1.5, 2.7, rep(0, k - 3))
        alpha0_true <- 1
        beta0_true <- 1.5
        gamma0_true <- 0.5
        tau0_true <- -1
    }
    list(scenario = scenario, K1 = p, K2 = q, K3 = k, K4 = m, p = p, q = q, k = k, m = m, rho = rho, 
        beta_true = beta_true, gamma_true = gamma_true, eta_true = eta_true, tau_true = tau_true, alpha0_true = alpha0_true, 
        beta0_true = beta0_true, gamma0_true = gamma0_true, tau0_true = tau0_true)
}


make_sigma <- function(dim, rho) {
    outer(1:dim, 1:dim, function(i, j) rho^abs(i - j))
}


generate_data_scenario <- function(n, scenario = c("Scenario I", "Scenario II", "Scenario III"), rho = 0.5) {
    cfg <- get_scenario_config(scenario, rho = rho)
    p <- cfg$p
    q <- cfg$q
    k <- cfg$k
    m <- cfg$m
    U <- mvrnorm(n, mu = rep(0, p), Sigma = make_sigma(p, rho))
    V <- mvrnorm(n, mu = rep(0, q), Sigma = make_sigma(q, rho))
    W <- mvrnorm(n, mu = rep(0, k), Sigma = make_sigma(k, rho))
    Z <- mvrnorm(n, mu = rep(0, m), Sigma = make_sigma(m, rho))
    x <- runif(n, -20, 20)
    knot_true <- cfg$tau0_true + as.vector(W %*% cfg$tau_true)
    x_pos <- pmax(x - knot_true, 0)
    y_pure <- cfg$alpha0_true + (cfg$beta0_true + as.vector(U %*% cfg$beta_true)) * x + (cfg$gamma0_true + 
        as.vector(V %*% cfg$gamma_true)) * x_pos + as.vector(Z %*% cfg$eta_true)
    y <- y_pure + rnorm(n, mean = 0, sd = 1)
    true_params <- c(alpha0 = cfg$alpha0_true, beta0 = cfg$beta0_true, gamma0 = cfg$gamma0_true, tau0 = cfg$tau0_true, 
        setNames(cfg$beta_true, paste0("beta_u_", 1:p)), setNames(cfg$gamma_true, paste0("gamma_v_", 
            1:q)), setNames(cfg$eta_true, paste0("eta_z_", 1:m)), setNames(cfg$tau_true, paste0("tau_w_", 
            1:k)))
    list(scenario = cfg$scenario, y = y, x = x, U = U, V = V, W = W, Z = Z, true_params = true_params, 
        K1 = p, K2 = q, K3 = k, K4 = m, p = p, q = q, k = k, m = m, rho = rho)
}


compute_nr_init <- function(y, x, U, V, W, Z, tol = 0.01, max_iter = 500) {
    x <- as.numeric(x)
    y <- as.numeric(y)
    U <- as.matrix(U)
    V <- as.matrix(V)
    W <- as.matrix(W)
    Z <- as.matrix(Z)
    p_u <- ncol(U)
    p_v <- ncol(V)
    p_w <- ncol(W)
    p_z <- ncol(Z)
    if (is.null(colnames(U)) && p_u > 0) 
        colnames(U) <- paste0("U", 1:p_u)
    if (is.null(colnames(V)) && p_v > 0) 
        colnames(V) <- paste0("V", 1:p_v)
    if (is.null(colnames(W)) && p_w > 0) 
        colnames(W) <- paste0("W", 1:p_w)
    if (is.null(colnames(Z)) && p_z > 0) 
        colnames(Z) <- paste0("Z", 1:p_z)
    tau <- c(median(x, na.rm = TRUE), rep(0, p_w))
    for (iter in 1:max_iter) {
        knot <- tau[1] + as.vector(W %*% tau[-1])
        x_pos <- pmax(x - knot, 0)
        x_u <- if (p_u > 0) 
            sweep(U, 1, x, "*")
        else matrix(nrow = length(x), ncol = 0)
        x_v <- if (p_v > 0) 
            sweep(V, 1, x_pos, "*")
        else matrix(nrow = length(x), ncol = 0)
        Xlin <- cbind(x, x_u, x_pos, x_v, Z)
        fit <- lm.fit(cbind(1, Xlin), y)
        b <- fit$coefficients
        b[is.na(b)] <- 0
        idx_gamma0 <- 3 + p_u
        idx_gamma_v <- if (p_v > 0) 
            (4 + p_u):(3 + p_u + p_v)
        else integer(0)
        gamma0_hat <- b[idx_gamma0]
        gamma_v_hat <- if (p_v > 0) 
            b[idx_gamma_v]
        else numeric(0)
        slope_v <- gamma0_hat + if (p_v > 0) 
            as.vector(V %*% gamma_v_hat)
        else 0
        I_gt <- as.numeric(x > knot)
        W_mat <- cbind(1, W)
        d_mu_dtau <- -as.vector(slope_v * I_gt) * W_mat
        resid <- y - as.vector(cbind(1, Xlin) %*% b)
        Un <- colMeans(resid * d_mu_dtau)
        Jn <- (t(d_mu_dtau) %*% d_mu_dtau)/length(y)
        step <- tryCatch(solve(Jn, Un), error = function(e) rep(NA_real_, length(Un)))
        if (any(is.na(step))) 
            break
        tau_new <- tau + as.numeric(step)
        if (max(abs(step)) < tol) {
            tau <- tau_new
            break
        }
        tau <- tau_new
    }
    knot <- tau[1] + as.vector(W %*% tau[-1])
    x_pos <- pmax(x - knot, 0)
    x_u <- if (p_u > 0) 
        sweep(U, 1, x, "*")
    else matrix(nrow = length(x), ncol = 0)
    x_v <- if (p_v > 0) 
        sweep(V, 1, x_pos, "*")
    else matrix(nrow = length(x), ncol = 0)
    Xlin <- cbind(x, x_u, x_pos, x_v, Z)
    fit <- lm.fit(cbind(1, Xlin), y)
    b <- fit$coefficients
    b[is.na(b)] <- 0
    beta_u <- if (p_u > 0) 
        b[3:(2 + p_u)]
    else numeric(0)
    gamma_v <- if (p_v > 0) 
        b[(4 + p_u):(3 + p_u + p_v)]
    else numeric(0)
    eta <- if (p_z > 0) 
        b[(4 + p_u + p_v):(3 + p_u + p_v + p_z)]
    else numeric(0)
    out <- c(alpha0 = b[1], beta0 = b[2], gamma0 = b[3 + p_u], tau0 = tau[1], beta_u, gamma_v, eta, tau[-1])
    names(out) <- c("alpha0", "beta0", "gamma0", "tau0", if (p_u > 0) paste0("beta_u_", 1:p_u) else NULL, 
        if (p_v > 0) paste0("gamma_v_", 1:p_v) else NULL, if (p_z > 0) paste0("eta_z_", 1:p_z) else NULL, 
        if (p_w > 0) paste0("tau_w_", 1:p_w) else NULL)
    out
}


ensure_named_full_estimate <- function(est_vec, true_params) {
    out <- setNames(rep(NA_real_, length(true_params)), names(true_params))
    common_names <- intersect(names(est_vec), names(true_params))
    out[common_names] <- est_vec[common_names]
    out
}


init_coef_store <- function() {
    data.frame(scenario = character(), n = integer(), method = character(), rep = integer(), param = character(), 
        truth = numeric(), estimate = numeric(), error = numeric(), stringsAsFactors = FALSE)
}


append_coef_store <- function(store, scenario, n, method, rep, est_vec, true_vec, only_true_nonzero = TRUE) {
    est_vec <- ensure_named_full_estimate(est_vec, true_vec)
    true_vec <- true_vec[names(est_vec)]
    keep_idx <- rep(TRUE, length(true_vec))
    if (only_true_nonzero) {
        keep_idx <- abs(true_vec) > 0
    }
    tmp <- data.frame(scenario = scenario, n = n, method = method, rep = rep, param = names(true_vec)[keep_idx], 
        truth = as.numeric(true_vec[keep_idx]), estimate = as.numeric(est_vec[keep_idx]), error = as.numeric(est_vec[keep_idx] - 
            true_vec[keep_idx]), stringsAsFactors = FALSE)
    rbind(store, tmp)
}


save_checkpoint <- function(coef_store, summary_store, save_dir, file_tag, inference_store = NULL, inference_summary_store = NULL, 
    inference_avg_store = NULL) {
    if (!dir.exists(save_dir)) 
        dir.create(save_dir, recursive = TRUE)
    saveRDS(coef_store, file = file.path(save_dir, paste0("coef_store_", file_tag, ".rds")))
    write.csv(coef_store, file = file.path(save_dir, paste0("coef_store_", file_tag, ".csv")), row.names = FALSE)
    write.csv(summary_store, file = file.path(save_dir, paste0("summary_", file_tag, ".csv")), row.names = FALSE)
    if (!is.null(inference_store)) {
        saveRDS(inference_store, file = file.path(save_dir, paste0("inference_store_", file_tag, ".rds")))
        write.csv(inference_store, file = file.path(save_dir, paste0("inference_store_", file_tag, ".csv")), 
            row.names = FALSE)
    }
    if (!is.null(inference_summary_store)) {
        write.csv(inference_summary_store, file = file.path(save_dir, paste0("inference_param_summary_", 
            file_tag, ".csv")), row.names = FALSE)
    }
    if (!is.null(inference_avg_store)) {
        write.csv(inference_avg_store, file = file.path(save_dir, paste0("inference_avg_summary_", file_tag, 
            ".csv")), row.names = FALSE)
    }
}


init_inference_store <- function() {
    data.frame(scenario = character(), n = integer(), method = character(), rep = integer(), param = character(), 
        truth = numeric(), estimate = numeric(), se = numeric(), lower = numeric(), upper = numeric(), 
        covered = logical(), selected = logical(), stringsAsFactors = FALSE)
}


append_inference_store <- function(store, scenario, n, method, rep, est_vec, se_vec, true_vec, active_params, 
    z_crit = 1.96) {
    est_vec <- ensure_named_full_estimate(est_vec, true_vec)
    se_full <- setNames(rep(NA_real_, length(true_vec)), names(true_vec))
    if (!is.null(se_vec)) {
        common_names <- intersect(names(se_vec), names(true_vec))
        se_full[common_names] <- se_vec[common_names]
    }
    params <- intersect(active_params, names(true_vec))
    est <- as.numeric(est_vec[params])
    se <- as.numeric(se_full[params])
    truth <- as.numeric(true_vec[params])
    lower <- est - z_crit * se
    upper <- est + z_crit * se
    valid_ci <- is.finite(lower) & is.finite(upper)
    covered <- valid_ci & (truth >= lower) & (truth <= upper)
    selected <- is.finite(est) & abs(est) > 0
    tmp <- data.frame(scenario = scenario, n = n, method = method, rep = rep, param = params, truth = truth, 
        estimate = est, se = se, lower = lower, upper = upper, covered = covered, selected = selected, 
        stringsAsFactors = FALSE)
    rbind(store, tmp)
}


compute_inference_param_summary <- function(inference_store) {
    if (nrow(inference_store) == 0) {
        return(data.frame())
    }
    spl <- split(inference_store, list(inference_store$scenario, inference_store$n, inference_store$method, 
        inference_store$param), drop = TRUE)
    out <- do.call(rbind, lapply(spl, function(d) {
        err <- d$estimate - d$truth
        data.frame(Scenario = unique(d$scenario), n = unique(d$n), Method = unique(d$method), Parameter = unique(d$param), 
            Bias = mean(err, na.rm = TRUE), AbsBias = mean(abs(err), na.rm = TRUE), MCSD = sd(d$estimate, 
                na.rm = TRUE), RMSE = sqrt(mean(err^2, na.rm = TRUE)), AVESE = mean(d$se, na.rm = TRUE), 
            CP = mean(d$covered, na.rm = TRUE) * 100, SelectedRate = mean(d$selected, na.rm = TRUE) * 
                100, Nvalid = sum(is.finite(d$estimate)), stringsAsFactors = FALSE)
    }))
    rownames(out) <- NULL
    out <- out[order(out$Scenario, out$n, out$Method, out$Parameter), ]
    out
}


compute_inference_avg_summary <- function(param_summary) {
    if (nrow(param_summary) == 0) {
        return(data.frame())
    }
    spl <- split(param_summary, list(param_summary$Scenario, param_summary$n, param_summary$Method), 
        drop = TRUE)
    out <- do.call(rbind, lapply(spl, function(d) {
        data.frame(Scenario = unique(d$Scenario), n = unique(d$n), Method = unique(d$Method), Bias = mean(d$Bias, 
            na.rm = TRUE), AbsBias = mean(d$AbsBias, na.rm = TRUE), MCSD = mean(d$MCSD, na.rm = TRUE), 
            RMSE = mean(d$RMSE, na.rm = TRUE), AVESE = mean(d$AVESE, na.rm = TRUE), CP = mean(d$CP, na.rm = TRUE), 
            SelectedRate = mean(d$SelectedRate, na.rm = TRUE), stringsAsFactors = FALSE)
    }))
    rownames(out) <- NULL
    out <- out[order(out$Scenario, out$n, out$Method), ]
    out
}


fit_global_onestep_penalty <- function(y, x, U, V, W, Z, penalty_type = c("SCAD", "MCP", "lasso")) {
    penalty_type <- match.arg(penalty_type)
    x <- as.numeric(x)
    y <- as.numeric(y)
    U <- as.matrix(U)
    V <- as.matrix(V)
    W <- as.matrix(W)
    Z <- as.matrix(Z)
    p_u <- ncol(U)
    p_v <- ncol(V)
    p_w <- ncol(W)
    p_z <- ncol(Z)
    if (is.null(colnames(U)) && p_u > 0) 
        colnames(U) <- paste0("U", 1:p_u)
    if (is.null(colnames(V)) && p_v > 0) 
        colnames(V) <- paste0("V", 1:p_v)
    if (is.null(colnames(W)) && p_w > 0) 
        colnames(W) <- paste0("W", 1:p_w)
    if (is.null(colnames(Z)) && p_z > 0) 
        colnames(Z) <- paste0("Z", 1:p_z)
    X_u_names <- if (p_u > 0) 
        paste0("U_", colnames(U))
    else character(0)
    X_v_names <- if (p_v > 0) 
        paste0("V_", colnames(V))
    else character(0)
    Z_names <- if (p_z > 0) 
        paste0("Z_", colnames(Z))
    else character(0)
    W_names <- if (p_w > 0) 
        paste0("W_", colnames(W))
    else character(0)
    init <- compute_nr_init(y, x, U, V, W, Z)
    gamma0_init <- init["gamma0"]
    tau0_init <- init["tau0"]
    gamma_v_init <- if (p_v > 0) 
        init[paste0("gamma_v_", 1:p_v)]
    else numeric(0)
    tau_w_init <- if (p_w > 0) 
        init[paste0("tau_w_", 1:p_w)]
    else numeric(0)
    knot_curr <- tau0_init + as.vector(W %*% tau_w_init)
    x_pos <- pmax(x - knot_curr, 0)
    x_u <- if (p_u > 0) 
        sweep(U, 1, x, "*")
    else matrix(nrow = length(x), ncol = 0)
    x_v <- if (p_v > 0) 
        sweep(V, 1, x_pos, "*")
    else matrix(nrow = length(x), ncol = 0)
    if (p_u > 0) 
        colnames(x_u) <- X_u_names
    if (p_v > 0) 
        colnames(x_v) <- X_v_names
    if (p_z > 0) 
        colnames(Z) <- Z_names
    X_main <- cbind(X_main = x, x_u, X_pos_main = x_pos, x_v, Z)
    slope_v <- gamma0_init + if (p_v > 0) 
        as.vector(V %*% gamma_v_init)
    else 0
    I_gt <- as.numeric(x > knot_curr)
    d_mu_d_tau0 <- -slope_v * I_gt
    d_mu_d_tauw <- if (p_w > 0) 
        sweep(W, 1, d_mu_d_tau0, "*")
    else NULL
    if (p_w > 0) 
        colnames(d_mu_d_tauw) <- paste0("d_tau_", W_names)
    X_tau <- cbind(d_tau0 = d_mu_d_tau0, d_mu_d_tauw)
    tau_init_vec <- c(tau0_init, tau_w_init)
    y_pseudo <- y + as.vector(X_tau %*% tau_init_vec)
    X_global <- cbind(X_main, X_tau)
    if (any(is.na(X_global)) || any(is.na(y_pseudo)) || any(is.infinite(X_global)) || any(is.infinite(y_pseudo))) {
        stop("X_global or y_pseudo contains NA/Inf.")
    }
    penalty.factor <- c(1, rep(1, p_u), 1, rep(1, p_v), rep(1, p_z), 1, rep(1, p_w))
    fit_pen <- ncvreg(X = X_global, y = y_pseudo, penalty = penalty_type, penalty.factor = penalty.factor, 
        eps = 0.001, max.iter = 1000)
    best_id <- which.min(BIC(fit_pen))
    raw_co <- as.numeric(coef(fit_pen, lambda = fit_pen$lambda[best_id]))
    names(raw_co) <- c("(Intercept)", colnames(X_global))
    final <- c(alpha0 = raw_co["(Intercept)"], beta0 = raw_co["X_main"], gamma0 = raw_co["X_pos_main"], 
        tau0 = raw_co["d_tau0"], if (p_u > 0) raw_co[X_u_names] else NULL, if (p_v > 0) raw_co[X_v_names] else NULL, 
        if (p_z > 0) raw_co[Z_names] else NULL, if (p_w > 0) raw_co[paste0("d_tau_", W_names)] else NULL)
    names(final) <- c("alpha0", "beta0", "gamma0", "tau0", if (p_u > 0) paste0("beta_u_", 1:p_u) else NULL, 
        if (p_v > 0) paste0("gamma_v_", 1:p_v) else NULL, if (p_z > 0) paste0("eta_z_", 1:p_z) else NULL, 
        if (p_w > 0) paste0("tau_w_", 1:p_w) else NULL)
    param_to_col <- c(beta0 = "X_main", gamma0 = "X_pos_main", tau0 = "d_tau0", setNames(X_u_names, paste0("beta_u_", 
        1:p_u)), setNames(X_v_names, paste0("gamma_v_", 1:p_v)), setNames(Z_names, paste0("eta_z_", 1:p_z)), 
        setNames(paste0("d_tau_", W_names), paste0("tau_w_", 1:p_w)))
    col_to_param <- setNames(names(param_to_col), unname(param_to_col))
    selected_cols <- names(raw_co)[-1][abs(raw_co[-1]) > 0]
    selected_cols <- intersect(selected_cols, colnames(X_global))
    refit_coef_full <- setNames(rep(0, length(final)), names(final))
    refit_se_full <- setNames(rep(NA_real_, length(final)), names(final))
    refit_coef_full["alpha0"] <- raw_co["(Intercept)"]
    inference_ok <- FALSE
    if (length(selected_cols) > 0) {
        X_refit <- X_global[, selected_cols, drop = FALSE]
        fit_refit <- tryCatch(lm.fit(cbind(`(Intercept)` = 1, X_refit), y_pseudo), error = function(e) NULL)
        if (!is.null(fit_refit)) {
            coef_refit <- fit_refit$coefficients
            coef_refit[is.na(coef_refit)] <- 0
            names(coef_refit) <- c("(Intercept)", selected_cols)
            resid_refit <- y_pseudo - as.vector(cbind(`(Intercept)` = 1, X_refit) %*% coef_refit)
            df_resid <- length(y_pseudo) - length(coef_refit)
            XtX_inv <- tryCatch(solve(crossprod(cbind(`(Intercept)` = 1, X_refit))), error = function(e) NULL)
            if (!is.null(XtX_inv) && df_resid > 0) {
                sigma2_hat <- sum(resid_refit^2)/df_resid
                se_refit <- sqrt(pmax(diag(sigma2_hat * XtX_inv), 0))
                names(se_refit) <- names(coef_refit)
                refit_coef_full["alpha0"] <- coef_refit["(Intercept)"]
                refit_se_full["alpha0"] <- se_refit["(Intercept)"]
                for (cc in selected_cols) {
                  pp <- col_to_param[cc]
                  if (!is.na(pp) && pp %in% names(refit_coef_full)) {
                    refit_coef_full[pp] <- coef_refit[cc]
                    refit_se_full[pp] <- se_refit[cc]
                  }
                }
                inference_ok <- TRUE
            }
        }
    }
    list(method = penalty_type, best_lambda = fit_pen$lambda[best_id], final_coef = final, refit_coef = refit_coef_full, 
        refit_se = refit_se_full, selected_cols = selected_cols, inference_ok = inference_ok, init_coef = init)
}


fit_onestep_ci_selection <- function(y, x, U, V, W, Z, z_crit = 1.96) {
    x <- as.numeric(x)
    y <- as.numeric(y)
    U <- as.matrix(U)
    V <- as.matrix(V)
    W <- as.matrix(W)
    Z <- as.matrix(Z)
    p_u <- ncol(U)
    p_v <- ncol(V)
    p_w <- ncol(W)
    p_z <- ncol(Z)
    X_u_names <- if (p_u > 0) 
        paste0("U_", 1:p_u)
    else character(0)
    X_v_names <- if (p_v > 0) 
        paste0("V_", 1:p_v)
    else character(0)
    Z_names <- if (p_z > 0) 
        paste0("Z_", 1:p_z)
    else character(0)
    W_names <- if (p_w > 0) 
        paste0("W_", 1:p_w)
    else character(0)
    init <- compute_nr_init(y, x, U, V, W, Z)
    gamma0_init <- unname(init["gamma0"])
    tau0_init <- unname(init["tau0"])
    gamma_v_init <- if (p_v > 0) 
        unname(init[paste0("gamma_v_", 1:p_v)])
    else numeric(0)
    tau_w_init <- if (p_w > 0) 
        unname(init[paste0("tau_w_", 1:p_w)])
    else numeric(0)
    knot_curr <- tau0_init + as.vector(W %*% tau_w_init)
    x_pos <- pmax(x - knot_curr, 0)
    x_u <- if (p_u > 0) 
        sweep(U, 1, x, "*")
    else matrix(nrow = length(x), ncol = 0)
    x_v <- if (p_v > 0) 
        sweep(V, 1, x_pos, "*")
    else matrix(nrow = length(x), ncol = 0)
    if (p_u > 0) 
        colnames(x_u) <- X_u_names
    if (p_v > 0) 
        colnames(x_v) <- X_v_names
    if (p_z > 0) 
        colnames(Z) <- Z_names
    X_main <- cbind(X_main = x, x_u, X_pos_main = x_pos, x_v, Z)
    slope_v <- gamma0_init + if (p_v > 0) 
        as.vector(V %*% gamma_v_init)
    else 0
    I_gt <- as.numeric(x > knot_curr)
    d_tau0 <- -slope_v * I_gt
    d_tauw <- if (p_w > 0) 
        sweep(W, 1, d_tau0, "*")
    else NULL
    if (p_w > 0) 
        colnames(d_tauw) <- paste0("d_tau_W_", 1:p_w)
    X_tau <- cbind(d_tau0 = d_tau0, d_tauw)
    y_pseudo <- y + as.vector(X_tau %*% c(tau0_init, tau_w_init))
    X_design <- cbind(`(Intercept)` = 1, X_main, X_tau)
    fit <- lm.fit(X_design, y_pseudo)
    coef_hat <- fit$coefficients
    coef_hat[is.na(coef_hat)] <- 0
    names(coef_hat) <- colnames(X_design)
    resid <- y_pseudo - as.vector(X_design %*% coef_hat)
    df_resid <- nrow(X_design) - ncol(X_design)
    if (df_resid <= 0) 
        stop("Residual df <= 0.")
    XtX_inv <- tryCatch(solve(crossprod(X_design)), error = function(e) NULL)
    if (is.null(XtX_inv)) 
        stop("X'X singular.")
    sigma2_hat <- sum(resid^2)/df_resid
    se_hat <- sqrt(pmax(diag(sigma2_hat * XtX_inv), 0))
    names(se_hat) <- colnames(X_design)
    select_cols <- setdiff(colnames(X_design), "(Intercept)")
    lower <- coef_hat[select_cols] - z_crit * se_hat[select_cols]
    upper <- coef_hat[select_cols] + z_crit * se_hat[select_cols]
    keep <- !(lower <= 0 & upper >= 0)
    keep[!is.finite(keep)] <- FALSE
    coef_sel <- coef_hat
    coef_sel[select_cols[!keep]] <- 0
    final <- c(alpha0 = coef_sel["(Intercept)"], beta0 = coef_sel["X_main"], gamma0 = coef_sel["X_pos_main"], 
        tau0 = coef_sel["d_tau0"], if (p_u > 0) coef_sel[X_u_names] else NULL, if (p_v > 0) coef_sel[X_v_names] else NULL, 
        if (p_z > 0) coef_sel[Z_names] else NULL, if (p_w > 0) coef_sel[paste0("d_tau_W_", 1:p_w)] else NULL)
    names(final) <- c("alpha0", "beta0", "gamma0", "tau0", if (p_u > 0) paste0("beta_u_", 1:p_u) else NULL, 
        if (p_v > 0) paste0("gamma_v_", 1:p_v) else NULL, if (p_z > 0) paste0("eta_z_", 1:p_z) else NULL, 
        if (p_w > 0) paste0("tau_w_", 1:p_w) else NULL)
    list(final_coef = final, se = se_hat, init_coef = init)
}


fit_oracle_onestep <- function(y, x, U, V, W, Z, true_params) {
    x <- as.numeric(x)
    y <- as.numeric(y)
    U <- as.matrix(U)
    V <- as.matrix(V)
    W <- as.matrix(W)
    Z <- as.matrix(Z)
    p_u <- ncol(U)
    p_v <- ncol(V)
    p_w <- ncol(W)
    p_z <- ncol(Z)
    X_u_names <- if (p_u > 0) 
        paste0("U_", 1:p_u)
    else character(0)
    X_v_names <- if (p_v > 0) 
        paste0("V_", 1:p_v)
    else character(0)
    Z_names <- if (p_z > 0) 
        paste0("Z_", 1:p_z)
    else character(0)
    W_names <- if (p_w > 0) 
        paste0("W_", 1:p_w)
    else character(0)
    init <- compute_nr_init(y, x, U, V, W, Z)
    gamma0_init <- unname(init["gamma0"])
    tau0_init <- unname(init["tau0"])
    gamma_v_init <- if (p_v > 0) 
        unname(init[paste0("gamma_v_", 1:p_v)])
    else numeric(0)
    tau_w_init <- if (p_w > 0) 
        unname(init[paste0("tau_w_", 1:p_w)])
    else numeric(0)
    knot_curr <- tau0_init + as.vector(W %*% tau_w_init)
    x_pos <- pmax(x - knot_curr, 0)
    x_u <- if (p_u > 0) 
        sweep(U, 1, x, "*")
    else matrix(nrow = length(x), ncol = 0)
    x_v <- if (p_v > 0) 
        sweep(V, 1, x_pos, "*")
    else matrix(nrow = length(x), ncol = 0)
    if (p_u > 0) 
        colnames(x_u) <- X_u_names
    if (p_v > 0) 
        colnames(x_v) <- X_v_names
    if (p_z > 0) 
        colnames(Z) <- Z_names
    X_main <- cbind(X_main = x, x_u, X_pos_main = x_pos, x_v, Z)
    slope_v <- gamma0_init + if (p_v > 0) 
        as.vector(V %*% gamma_v_init)
    else 0
    I_gt <- as.numeric(x > knot_curr)
    d_mu_d_tau0 <- -slope_v * I_gt
    d_mu_d_tauw <- if (p_w > 0) 
        sweep(W, 1, d_mu_d_tau0, "*")
    else NULL
    if (p_w > 0) 
        colnames(d_mu_d_tauw) <- paste0("d_tau_", W_names)
    X_tau <- cbind(d_tau0 = d_mu_d_tau0, d_mu_d_tauw)
    tau_init_vec <- c(tau0_init, tau_w_init)
    y_pseudo <- y + as.vector(X_tau %*% tau_init_vec)
    X_global <- cbind(X_main, X_tau)
    active_param_names <- setdiff(names(true_params)[abs(true_params) > 0], "alpha0")
    param_to_col <- c(beta0 = "X_main", gamma0 = "X_pos_main", tau0 = "d_tau0", setNames(X_u_names, paste0("beta_u_", 
        1:p_u)), setNames(X_v_names, paste0("gamma_v_", 1:p_v)), setNames(Z_names, paste0("eta_z_", 1:p_z)), 
        setNames(paste0("d_tau_", W_names), paste0("tau_w_", 1:p_w)))
    active_cols <- unname(param_to_col[active_param_names])
    active_cols <- active_cols[!is.na(active_cols)]
    X_oracle <- if (length(active_cols) > 0) 
        X_global[, active_cols, drop = FALSE]
    else NULL
    fit <- if (is.null(X_oracle)) {
        lm.fit(cbind(1), y_pseudo)
    }
    else {
        lm.fit(cbind(1, X_oracle), y_pseudo)
    }
    coef_red <- fit$coefficients
    coef_red[is.na(coef_red)] <- 0
    names(coef_red) <- if (is.null(X_oracle)) 
        "(Intercept)"
    else c("(Intercept)", colnames(X_oracle))
    final <- setNames(rep(0, length(true_params)), names(true_params))
    final["alpha0"] <- coef_red["(Intercept)"]
    if ("X_main" %in% names(coef_red)) 
        final["beta0"] <- coef_red["X_main"]
    if ("X_pos_main" %in% names(coef_red)) 
        final["gamma0"] <- coef_red["X_pos_main"]
    if ("d_tau0" %in% names(coef_red)) 
        final["tau0"] <- coef_red["d_tau0"]
    for (j in 1:p_u) {
        nm <- paste0("beta_u_", j)
        colnm <- paste0("U_", j)
        if (colnm %in% names(coef_red)) 
            final[nm] <- coef_red[colnm]
    }
    for (j in 1:p_v) {
        nm <- paste0("gamma_v_", j)
        colnm <- paste0("V_", j)
        if (colnm %in% names(coef_red)) 
            final[nm] <- coef_red[colnm]
    }
    for (j in 1:p_z) {
        nm <- paste0("eta_z_", j)
        colnm <- paste0("Z_", j)
        if (colnm %in% names(coef_red)) 
            final[nm] <- coef_red[colnm]
    }
    for (j in 1:p_w) {
        nm <- paste0("tau_w_", j)
        colnm <- paste0("d_tau_W_", j)
        if (colnm %in% names(coef_red)) 
            final[nm] <- coef_red[colnm]
    }
    list(final_coef = final, active_params = c("alpha0", active_param_names), init_coef = init)
}


build_prwd_working <- function(dat, helpers = list(compute_nr_init = compute_nr_init)) {
    x <- as.numeric(dat$x)
    y <- as.numeric(dat$y)
    U <- as.matrix(dat$U)
    V <- as.matrix(dat$V)
    W <- as.matrix(dat$W)
    Z <- as.matrix(dat$Z)
    p_u <- ncol(U)
    p_v <- ncol(V)
    p_w <- ncol(W)
    p_z <- ncol(Z)
    if (is.null(colnames(U)) && p_u > 0) 
        colnames(U) <- paste0("U", 1:p_u)
    if (is.null(colnames(V)) && p_v > 0) 
        colnames(V) <- paste0("V", 1:p_v)
    if (is.null(colnames(W)) && p_w > 0) 
        colnames(W) <- paste0("W", 1:p_w)
    if (is.null(colnames(Z)) && p_z > 0) 
        colnames(Z) <- paste0("Z", 1:p_z)
    U_names <- if (p_u > 0) 
        paste0("U_", colnames(U))
    else character(0)
    V_names <- if (p_v > 0) 
        paste0("V_", colnames(V))
    else character(0)
    Z_names <- if (p_z > 0) 
        paste0("Z_", colnames(Z))
    else character(0)
    W_names <- if (p_w > 0) 
        paste0("W_", colnames(W))
    else character(0)
    init <- helpers$compute_nr_init(y, x, U, V, W, Z)
    gamma0_init <- init["gamma0"]
    tau0_init <- init["tau0"]
    gamma_v_init <- if (p_v > 0) 
        init[paste0("gamma_v_", 1:p_v)]
    else numeric(0)
    tau_w_init <- if (p_w > 0) 
        init[paste0("tau_w_", 1:p_w)]
    else numeric(0)
    knot_curr <- tau0_init + as.vector(W %*% tau_w_init)
    x_pos <- pmax(x - knot_curr, 0)
    x_u <- if (p_u > 0) 
        sweep(U, 1, x, "*")
    else matrix(nrow = length(x), ncol = 0)
    x_v <- if (p_v > 0) 
        sweep(V, 1, x_pos, "*")
    else matrix(nrow = length(x), ncol = 0)
    if (p_u > 0) 
        colnames(x_u) <- U_names
    if (p_v > 0) 
        colnames(x_v) <- V_names
    if (p_z > 0) 
        colnames(Z) <- Z_names
    X_main <- cbind(X_main = x, x_u, X_pos_main = x_pos, x_v, Z)
    slope_v <- gamma0_init + if (p_v > 0) 
        as.vector(V %*% gamma_v_init)
    else 0
    d_mu_d_tau0 <- -slope_v * as.numeric(x > knot_curr)
    d_mu_d_tauw <- if (p_w > 0) 
        sweep(W, 1, d_mu_d_tau0, "*")
    else NULL
    if (p_w > 0) 
        colnames(d_mu_d_tauw) <- paste0("d_tau_", W_names)
    X_tau <- cbind(d_tau0 = d_mu_d_tau0, d_mu_d_tauw)
    y_pseudo <- y + as.vector(X_tau %*% c(tau0_init, tau_w_init))
    X_global <- cbind(X_main, X_tau)
    if (any(!is.finite(X_global)) || any(!is.finite(y_pseudo))) 
        stop("Non-finite working design/response.")
    param_to_col <- c(beta0 = "X_main", gamma0 = "X_pos_main", tau0 = "d_tau0", setNames(U_names, paste0("beta_u_", 
        seq_len(p_u))), setNames(V_names, paste0("gamma_v_", seq_len(p_v))), setNames(Z_names, paste0("eta_z_", 
        seq_len(p_z))), setNames(paste0("d_tau_", W_names), paste0("tau_w_", seq_len(p_w))))
    list(X = X_global, y = y_pseudo, init = init, param_to_col = param_to_col)
}


default_solver_controls <- function() {
    list(algorithm = "CDPSI", nLambda = 100L, maxIters = 2000L, retryMaxIters = 10000L, rtol = 1e-07, 
        atol = 1e-10, maxSwaps = 100L, scaleDownFactor = 0.8)
}


gaussian_bic <- function(rss, k, n) {
    n * log(pmax(rss/n, .Machine$double.xmin)) + k * log(n)
}


fit_l0_working <- function(X, y, controls = default_solver_controls()) {
    if (!requireNamespace("L0Learn", quietly = TRUE)) 
        stop("L0Learn is not installed in the configured R library.")
    X <- as.matrix(X)
    y <- as.numeric(y)
    if (any(!is.finite(X)) || any(!is.finite(y))) 
        stop("Non-finite L0 inputs.")
    if (!ncol(X) || nrow(X) != length(y)) 
        stop("Invalid working design dimensions.")
    had_rng <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_rng) 
        rng_before <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({
        if (had_rng) assign(".Random.seed", rng_before, envir = .GlobalEnv) else if (exists(".Random.seed", 
            envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    run_path <- function(max_iters) {
        warnings <- character()
        fit <- withCallingHandlers(L0Learn::L0Learn.fit(x = X, y = y, loss = "SquaredError", penalty = "L0", 
            algorithm = controls$algorithm, maxSuppSize = min(ncol(X), nrow(X) - 2L), nLambda = controls$nLambda, 
            nGamma = 1L, maxIters = max_iters, rtol = controls$rtol, atol = controls$atol, maxSwaps = controls$maxSwaps, 
            scaleDownFactor = controls$scaleDownFactor, screenSize = ncol(X), intercept = TRUE, excludeFirstK = 0L), 
            warning = function(w) {
                warnings <<- c(warnings, conditionMessage(w))
                invokeRestart("muffleWarning")
            })
        if (length(fit$beta) != 1L) 
            stop("Pure L0 should return one gamma path.")
        beta <- as.matrix(fit$beta[[1L]])
        a0 <- as.numeric(fit$a0[[1L]])
        converged <- as.logical(fit$converged[[1L]])
        if (nrow(beta) != ncol(X) || ncol(beta) != length(a0) || length(converged) != length(a0)) 
            stop("Malformed L0 path.")
        rss <- colSums((y - sweep(X %*% beta, 2L, a0, "+"))^2)
        k <- colSums(beta != 0)
        bic <- gaussian_bic(rss, k, length(y))
        finite <- is.finite(rss) & is.finite(bic) & is.finite(a0) & colSums(!is.finite(beta)) == 0L
        if (!any(finite)) 
            stop("No finite candidate in L0 regularization path.")
        list(fit = fit, beta = beta, a0 = a0, converged = converged, rss = rss, k = k, bic = bic, finite = finite, 
            warnings = warnings)
    }
    path <- run_path(controls$maxIters)
    best_all <- which.min(ifelse(path$finite, path$bic, Inf))
    retried <- FALSE
    if (!isTRUE(path$converged[best_all])) {
        path <- run_path(controls$retryMaxIters)
        retried <- TRUE
        best_all <- which.min(ifelse(path$finite, path$bic, Inf))
    }
    if (!isTRUE(path$converged[best_all])) 
        stop("BIC-best L0 candidate did not converge, including the high-iteration retry.")
    eligible <- path$finite & !is.na(path$converged) & path$converged
    best <- which.min(ifelse(eligible, path$bic, Inf))
    coef <- c(`(Intercept)` = path$a0[best], setNames(path$beta[, best], colnames(X)))
    list(coefficients = coef, diagnostics = list(lambda = as.numeric(path$fit$lambda[[1L]][best]), bic = path$bic[best], 
        rss = path$rss[best], support_size = path$k[best], path_candidates = length(path$bic), selected_path_index = best, 
        selected_path_endpoint = best == length(path$bic), nonconverged_candidates = sum(!path$converged, 
            na.rm = TRUE), retried = retried, selected_converged = TRUE, warnings = path$warnings, controls = controls))
}


fit_l0_onestep <- function(dat, helpers = list(compute_nr_init = compute_nr_init), controls = default_solver_controls()) {
    work <- build_prwd_working(dat, helpers)
    ans <- fit_l0_working(work$X, work$y, controls)
    final <- setNames(rep(NA_real_, length(dat$true_params)), names(dat$true_params))
    final["alpha0"] <- ans$coefficients["(Intercept)"]
    for (param in setdiff(names(final), "alpha0")) final[param] <- ans$coefficients[work$param_to_col[param]]
    if (any(!is.finite(final))) 
        stop("L0 coefficient mapping produced missing/non-finite original parameters.")
    list(method = "L0", best_lambda = ans$diagnostics$lambda, final_coef = final, init_coef = work$init, 
        diagnostics = ans$diagnostics)
}


run_all_methods_on_dataset <- function(dat, z_crit = 1.96) {
    out <- list()
    out[["SCAD"]] <- tryCatch(fit_global_onestep_penalty(dat$y, dat$x, dat$U, dat$V, dat$W, dat$Z, penalty_type = "SCAD"), 
        error = function(e) NULL)
    out[["MCP"]] <- tryCatch(fit_global_onestep_penalty(dat$y, dat$x, dat$U, dat$V, dat$W, dat$Z, penalty_type = "MCP"), 
        error = function(e) NULL)
    out[["lasso"]] <- tryCatch(fit_global_onestep_penalty(dat$y, dat$x, dat$U, dat$V, dat$W, dat$Z, penalty_type = "lasso"), 
        error = function(e) NULL)
    out[["L0"]] <- tryCatch(fit_l0_onestep(dat), error = function(e) {
        warning(paste("L0 fit failed:", conditionMessage(e)), call. = FALSE)
        NULL
    })
    out[["CI-based"]] <- tryCatch(fit_onestep_ci_selection(dat$y, dat$x, dat$U, dat$V, dat$W, dat$Z, 
        z_crit = z_crit), error = function(e) NULL)
    out[["Oracle"]] <- tryCatch(fit_oracle_onestep(dat$y, dat$x, dat$U, dat$V, dat$W, dat$Z, true_params = dat$true_params), 
        error = function(e) NULL)
    out
}


compute_one_method_summary <- function(results_matrix, true_vals, knot_rmse_vec) {
    keep <- complete.cases(results_matrix)
    results_matrix <- results_matrix[keep, , drop = FALSE]
    knot_rmse_vec <- knot_rmse_vec[keep]
    if (nrow(results_matrix) == 0) {
        return(data.frame(TN = NA_real_, FN = NA_real_, OP = NA_real_, L2E = NA_real_, CPE = NA_real_, 
            Tau_TN = NA_real_, Tau_FN = NA_real_, Tau_OP = NA_real_, stringsAsFactors = FALSE))
    }
    true_vec <- as.numeric(true_vals)
    names(true_vec) <- names(true_vals)
    check_idx <- setdiff(seq_along(true_vec), match("alpha0", names(true_vec)))
    sel_est <- abs(results_matrix[, check_idx, drop = FALSE]) > 0
    sel_true <- abs(true_vec[check_idx]) > 0
    TN <- mean(rowSums(!sel_est[, !sel_true, drop = FALSE]))
    FN <- mean(rowSums(!sel_est[, sel_true, drop = FALSE]))
    OP <- mean(apply(sel_est, 1, function(r) all(r == sel_true))) * 100
    L2E <- mean(apply(results_matrix, 1, function(r) sqrt(sum((r - true_vec)^2))))
    CPE <- mean(knot_rmse_vec)
    tau_idx <- grep("^tau_w_", names(true_vec))
    tau_sel_est <- abs(results_matrix[, tau_idx, drop = FALSE]) > 0
    tau_sel_true <- abs(true_vec[tau_idx]) > 0
    Tau_TN <- mean(rowSums(!tau_sel_est[, !tau_sel_true, drop = FALSE]))
    Tau_FN <- mean(rowSums(!tau_sel_est[, tau_sel_true, drop = FALSE]))
    Tau_OP <- mean(apply(tau_sel_est, 1, function(r) all(r == tau_sel_true))) * 100
    data.frame(TN = TN, FN = FN, OP = OP, L2E = L2E, CPE = CPE, Tau_TN = Tau_TN, Tau_FN = Tau_FN, Tau_OP = Tau_OP, 
        stringsAsFactors = FALSE)
}


run_all_simulations <- function(scenarios = c("Scenario I", "Scenario II", "Scenario III"), n_list = c(250, 
    500, 1000, 2500), n_reps = 1000, rho = 0.5, z_crit = 1.96, seed = 20262027, checkpoint_every = 20, 
    save_dir = "sim_results") {
    set.seed(seed)
    if (!dir.exists(save_dir)) 
        dir.create(save_dir, recursive = TRUE)
    all_coef_store <- init_coef_store()
    all_summary <- data.frame()
    all_l0_summary <- data.frame()
    inference_method_names <- c("SCAD", "MCP")
    active_tau_params <- c("tau0", "tau_w_1", "tau_w_2", "tau_w_3")
    all_inference_store <- init_inference_store()
    all_inference_param_summary <- data.frame()
    all_inference_avg_summary <- data.frame()
    method_names <- c("SCAD", "MCP", "lasso", "L0", "CI-based", "Oracle")
    for (scenario_name in scenarios) {
        cat("\n============================================================\n")
        cat("Running", scenario_name, "\n")
        for (n_sample in n_list) {
            cat("------------------------------------------------------------\n")
            cat(sprintf("Scenario = %s, n = %d\n", scenario_name, n_sample))
            template_dat <- generate_data_scenario(n = 10, scenario = scenario_name, rho = rho)
            true_vals <- template_dat$true_params
            p_total <- length(true_vals)
            tau_names <- grep("^tau_w_", names(true_vals), value = TRUE)
            method_results <- lapply(method_names, function(m) {
                matrix(NA_real_, nrow = n_reps, ncol = p_total, dimnames = list(NULL, names(true_vals)))
            })
            names(method_results) <- method_names
            method_knot_rmse <- lapply(method_names, function(m) rep(NA_real_, n_reps))
            names(method_knot_rmse) <- method_names
            for (i in 1:n_reps) {
                dat <- generate_data_scenario(n = n_sample, scenario = scenario_name, rho = rho)
                fits <- run_all_methods_on_dataset(dat, z_crit = z_crit)
                for (method_name in method_names) {
                  ans <- fits[[method_name]]
                  if (is.null(ans)) 
                    next
                  est <- ensure_named_full_estimate(ans$final_coef, dat$true_params)
                  method_results[[method_name]][i, ] <- est
                  knot_true <- dat$true_params["tau0"] + as.vector(dat$W %*% dat$true_params[tau_names])
                  knot_hat <- est["tau0"] + as.vector(dat$W %*% est[tau_names])
                  method_knot_rmse[[method_name]][i] <- sqrt(mean((knot_hat - knot_true)^2))
                  all_coef_store <- append_coef_store(store = all_coef_store, scenario = scenario_name, 
                    n = n_sample, method = method_name, rep = i, est_vec = est, true_vec = dat$true_params, 
                    only_true_nonzero = TRUE)
                  if (method_name %in% inference_method_names && !is.null(ans$refit_coef) && !is.null(ans$refit_se)) {
                    all_inference_store <- append_inference_store(store = all_inference_store, scenario = scenario_name, 
                      n = n_sample, method = method_name, rep = i, est_vec = ans$refit_coef, se_vec = ans$refit_se, 
                      true_vec = dat$true_params, active_params = active_tau_params, z_crit = z_crit)
                  }
                }
                if (i%%checkpoint_every == 0) {
                  cat(sprintf("  finished %d / %d\n", i, n_reps))
                  all_inference_param_summary <- compute_inference_param_summary(all_inference_store)
                  all_inference_avg_summary <- compute_inference_avg_summary(all_inference_param_summary)
                  save_checkpoint(coef_store = all_coef_store, summary_store = all_summary, save_dir = save_dir, 
                    file_tag = "checkpoint", inference_store = all_inference_store, inference_summary_store = all_inference_param_summary, 
                    inference_avg_store = all_inference_avg_summary)
                }
            }
            for (method_name in method_names) {
                summ <- compute_one_method_summary(results_matrix = method_results[[method_name]], true_vals = true_vals, 
                  knot_rmse_vec = method_knot_rmse[[method_name]])
                summ$Scenario <- scenario_name
                summ$n <- n_sample
                summ$Method <- method_name
                all_summary <- rbind(all_summary, summ[, c("Scenario", "n", "Method", "TN", "FN", "OP", 
                  "L2E", "CPE", "Tau_TN", "Tau_FN", "Tau_OP")])
                if (method_name == "L0") {
                  l0_summ <- summ[, c("Scenario", "n", "Method", "TN", "FN", "OP", "L2E", "CPE", "Tau_TN", 
                    "Tau_FN", "Tau_OP")]
                  l0_summ$NRequested <- n_reps
                  l0_summ$NProcessed <- n_reps
                  l0_summ$NSuccess <- sum(complete.cases(method_results[["L0"]]))
                  l0_summ$NFailed <- n_reps - l0_summ$NSuccess
                  all_l0_summary <- rbind(all_l0_summary, l0_summ)
                }
            }
            all_inference_param_summary <- compute_inference_param_summary(all_inference_store)
            all_inference_avg_summary <- compute_inference_avg_summary(all_inference_param_summary)
            write.csv(all_summary, file = file.path(save_dir, paste0("summary_up_to_", gsub(" ", "", 
                scenario_name), "_n", n_sample, ".csv")), row.names = FALSE)
            write.csv(all_inference_param_summary, file = file.path(save_dir, paste0("inference_param_summary_up_to_", 
                gsub(" ", "", scenario_name), "_n", n_sample, ".csv")), row.names = FALSE)
            write.csv(all_inference_avg_summary, file = file.path(save_dir, paste0("inference_avg_summary_up_to_", 
                gsub(" ", "", scenario_name), "_n", n_sample, ".csv")), row.names = FALSE)
            saveRDS(all_coef_store, file = file.path(save_dir, paste0("coef_store_up_to_", gsub(" ", 
                "", scenario_name), "_n", n_sample, ".rds")))
            saveRDS(all_inference_store, file = file.path(save_dir, paste0("inference_store_up_to_", 
                gsub(" ", "", scenario_name), "_n", n_sample, ".rds")))
        }
    }
    if (!dir.exists(save_dir)) 
        dir.create(save_dir, recursive = TRUE)
    all_inference_param_summary <- compute_inference_param_summary(all_inference_store)
    all_inference_avg_summary <- compute_inference_avg_summary(all_inference_param_summary)
    scad_inference_avg <- subset(all_inference_avg_summary, Method == "SCAD")
    mcp_inference_avg <- subset(all_inference_avg_summary, Method == "MCP")
    scad_inference_param <- subset(all_inference_param_summary, Method == "SCAD")
    mcp_inference_param <- subset(all_inference_param_summary, Method == "MCP")
    saveRDS(all_coef_store, file = file.path(save_dir, "all_coef_store_nonzero.rds"))
    write.csv(all_coef_store, file = file.path(save_dir, "all_coef_store_nonzero.csv"), row.names = FALSE)
    write.csv(all_summary, file = file.path(save_dir, "all_summary_results.csv"), row.names = FALSE)
    write.csv(all_l0_summary, file = file.path(save_dir, "L0_summary_all.csv"), row.names = FALSE)
    tau_selection_summary <- all_summary[, c("Scenario", "n", "Method", "Tau_TN", "Tau_FN", "Tau_OP")]
    write.csv(tau_selection_summary, file = file.path(save_dir, "all_tau_selection_summary.csv"), row.names = FALSE)
    saveRDS(all_inference_store, file = file.path(save_dir, "all_inference_store_active_tau.rds"))
    write.csv(all_inference_store, file = file.path(save_dir, "all_inference_store_active_tau.csv"), 
        row.names = FALSE)
    write.csv(all_inference_param_summary, file = file.path(save_dir, "all_inference_param_summary_active_tau.csv"), 
        row.names = FALSE)
    write.csv(all_inference_avg_summary, file = file.path(save_dir, "all_inference_avg_summary_active_tau.csv"), 
        row.names = FALSE)
    write.csv(scad_inference_avg, file = file.path(save_dir, "SCAD_inference_avg_summary_active_tau.csv"), 
        row.names = FALSE)
    write.csv(mcp_inference_avg, file = file.path(save_dir, "MCP_inference_avg_summary_active_tau.csv"), 
        row.names = FALSE)
    write.csv(scad_inference_param, file = file.path(save_dir, "SCAD_inference_param_summary_active_tau.csv"), 
        row.names = FALSE)
    write.csv(mcp_inference_param, file = file.path(save_dir, "MCP_inference_param_summary_active_tau.csv"), 
        row.names = FALSE)
    invisible(list(coef_store = all_coef_store, summary = all_summary, L0_summary = all_l0_summary, inference_store = all_inference_store, 
        inference_param_summary = all_inference_param_summary, inference_avg_summary = all_inference_avg_summary, 
        SCAD_inference_avg = scad_inference_avg, MCP_inference_avg = mcp_inference_avg, tau_selection_summary = tau_selection_summary))
}



# ----- Run -----
results <- run_all_simulations(
  scenarios = scenarios, n_list = n_list, n_reps = n_reps,
  rho = 0.5, z_crit = 1.96, seed = seed,
  checkpoint_every = min(20L, n_reps), save_dir = save_dir
)
cat("Results saved to:", normalizePath(save_dir, winslash = "/"), "\n")
