library(MASS)
library(policytree)
library(Matching)
library(caret)
library(parallel)
library(snowfall)
library(policytree)
library(grf)
library(DiagrammeR)
library(glmnet)
library(spatialreg)
library(spdep)
library(mgcv)
library(data.table)
library(corrplot)
library(MatchIt)
library(Matching)
library(fields)
library(ggplot2)
library(DAPSm)
library(stringr)
library(arepa)
library(optmatch)
library(lpSolve)
library(data.table)

set.seed(1)

inv_logit <- function(x) {
  1 / (1 + exp(-x))
}

matern_cor <- function(d, range = 0.25, nu = 1.46) {
  x <- d / range
  out <- numeric(length(x))
  zero_idx <- (d == 0)
  out[zero_idx] <- 1
  nz <- !zero_idx
  x_nz <- x[nz]
  out[nz] <- (1 / (2^(nu - 1) * gamma(nu))) *
    (x_nz^nu) * besselK(x_nz, nu)
  out
}




generate_dataset <- function(
    n = 800, 
    nu = 1.46,
    range = 0.25,
    sigma_eps = 1,
    seed=1,
    scenario="tree_X"
) {
  set.seed(seed)
  lon <- runif(n, 0, 1)
  lat <- runif(n, 0, 1)
  coords <- cbind(lon = lon, lat = lat)
  
  dmat <- fields::rdist(coords)

  R <- matrix(
    matern_cor(as.vector(dmat), range = range, nu = nu),
    nrow = n, ncol = n
  )

  diag(R) <- diag(R) + 1e-8
  U <- as.numeric(MASS::mvrnorm(n = 1, mu = rep(0, n), Sigma = R))
  U <- scale(U)[, 1] 
  X1 <- rnorm(n, 0, 1)
  X2 <- rnorm(n, 0, 1)
  X3 <- rnorm(n, 0, 1)
  X4 <- rnorm(n, 0, 1)
  

  lp_w <- -0.2 + 0.1 * X1 + 0.2 * X2 - 0.3 * X3 - 0.2 * X4 + 0.8 * U
  p_true <- inv_logit(lp_w)
  W <- rbinom(n, size = 1, prob = p_true)
  mu0 <- 2 + 0.50 * X1 + 0.40 * X2 + 0.50 * X3 - 0.70 * X4 + 2.00 * U + 0.75 * X1 * U

  if (scenario == "tree_X") {
    rule_tree <- (X1 > 0 & X2 > 0)
    tau <- 0.75 * (2 * as.numeric(rule_tree) - 1) 
  } 
  else if (scenario == "nontree_X") {
    rule_nontree <- (2 * X2 - exp(1 + X1) + 2 > 0)
    tau <- 0.75 * (2 * as.numeric(rule_nontree) - 1) 
  }
  else if (scenario == "tree_XU"){
    rule_tree <- (X1>0 & U<0)
    tau <- 0.75 * (2 * as.numeric(rule_tree) - 1)
  }
  else if (scenario == "nontree_XU"){
    rule_nontree <- sqrt(2) * (X2 - U) - exp(1 + X1) + 2 > 0
    tau <- 0.75 * (2 * as.numeric(rule_nontree) - 1)
  }
  
  
  
  Y0_mean <- mu0
  Y1_mean <- mu0 + tau
  eps <- rnorm(n, mean = 0, sd = sigma_eps)
  Y0_true <- Y0_mean + eps
  Y1_true <- Y1_mean + eps
  Y <- ifelse(W == 1, Y1_true, Y0_true)
  
  dat <- data.frame(
    id = 1:n,
    Longitude = lon,
    Latitude = lat,
    U = U,
    X1 = X1,
    X2 = X2,
    X3 = X3,
    X4 = X4,
    W = W,
    Y = Y,
    p_true = p_true,
    tau=tau,
    Y0_mean = Y0_mean,
    Y1_mean = Y1_mean,
    Y0_true = Y0_true,
    Y1_true = Y1_true
  )
  

  return(dat)
  
}






#################################################
generate_train_test_data <- function(
    n_train = 800,
    n_test = 2000,
    nu = 1.46,
    range = 0.25,
    sigma_eps = 1,
    seed = 1,
    scenario = "tree_X"
) {

  full.dat <- generate_dataset(
    n = n_train + n_test,
    nu = nu,
    range = range,
    sigma_eps = sigma_eps,
    seed = seed,
    scenario = scenario
  )

  train_idx <- seq_len(n_train)
  test_idx <- n_train + seq_len(n_test)
  
  train.dat <- full.dat[train_idx,,drop = FALSE]
  test.dat <- full.dat[test_idx,,drop = FALSE]

  train.dat$id <- seq_len(n_train)
  test.dat$id <- seq_len(n_test)
  

  naive_ps <- glm(
    W ~ X1 + X2 + X3 + X4,
    data = train.dat,
    family = "binomial"
  )
  
  train.dat$prop_score <- as.numeric(fitted(naive_ps))
  
  opt_rule <- as.integer(test.dat$tau > 0)
  
  opt_value <- mean(ifelse(opt_rule == 1,test.dat$Y1_true,test.dat$Y0_true))
  
  return(
    list(
      train.dat = train.dat,
      test.dat = test.dat,
      opt_rule = opt_rule,
      opt_value = opt_value
    )
  )
}





StandDist <- function(M) {
  rng <- range(M, finite = TRUE)
  if (!is.finite(rng[1]) ||
      !is.finite(rng[2]) ||
      rng[1] == rng[2]) {
    return(matrix(0, nrow = nrow(M), ncol = ncol(M)))
  }
  (M - rng[1]) / (rng[2] - rng[1])
}



calculate_asdm_for_weight <- function(weight,
                                      W,
                                      Y,
                                      X,
                                      S,
                                      prop_score = NULL,
                                      m = 2,
                                      caliper = 3,
                                      caliper_type = "DAPS",
                                      coord_dist = FALSE) {
  W <- as.integer(W)
  X <- as.matrix(X)
  S <- as.matrix(S)
  prop_score <- as.numeric(prop_score)
  
  treated_pos <- which(W == 1)
  control_pos <- which(W == 0)
  n_treated <- length(treated_pos)
  n_control <- length(control_pos)
  n_total <- length(W)
  
  if (n_treated == 0 || n_control == 0) {
    stop("Both treated and control units are required.")
  }
  
  treated <- data.frame(
    Longitude = S[treated_pos, 1],
    Latitude = S[treated_pos, 2],
    prop.scores = prop_score[treated_pos]
  )
  control <- data.frame(
    Longitude = S[control_pos, 1],
    Latitude = S[control_pos, 2],
    prop.scores = prop_score[control_pos]
  )
  
  if (coord_dist) {
    dist.mat <- fields::rdist.earth(
      cbind(treated$Longitude, treated$Latitude),
      cbind(control$Longitude, control$Latitude)
    )
  } else {
    dist.mat <- fields::rdist(
      cbind(treated$Longitude, treated$Latitude),
      cbind(control$Longitude, control$Latitude)
    )
  }
  
  stand.dist.mat <- StandDist(dist.mat)
  ps.diff <- abs(outer(treated$prop.scores, control$prop.scores, "-"))
  dapscore <- (1 - weight) * stand.dist.mat + weight * ps.diff
  
  M <- dapscore
  if (is.finite(caliper) && caliper_type == "DAPS") {
    threshold <- caliper * sd(as.vector(dapscore), na.rm = TRUE)
    if (is.finite(threshold) && threshold > 0) {
      M[dapscore > threshold] <- Inf
    }
  } else if (is.finite(caliper) && caliper_type == "PS") {
    threshold <- caliper * sd(prop_score, na.rm = TRUE)
    if (is.finite(threshold) && threshold > 0) {
      M[ps.diff > threshold] <- Inf
    }
  }
  
  treated_matches <- matrix(NA_integer_, nrow = n_treated, ncol = m)
  dropped_treated <- logical(n_treated)
  for (i in seq_len(n_treated)) {
    dist_i <- M[i, ]
    ok <- which(is.finite(dist_i))
    if (length(ok) < m) {
      dropped_treated[i] <- TRUE
      next
    }
    sorted_control <- ok[order(dist_i[ok])][seq_len(m)]
    treated_matches[i, ] <- control_pos[sorted_control]
  }
  
  control_matches <- matrix(NA_integer_, nrow = n_control, ncol = m)
  dropped_control <- logical(n_control)
  M_t <- t(M)
  for (j in seq_len(n_control)) {
    dist_j <- M_t[j, ]
    ok <- which(is.finite(dist_j))
    if (length(ok) < m) {
      dropped_control[j] <- TRUE
      next
    }
    sorted_treated <- ok[order(dist_j[ok])][seq_len(m)]
    control_matches[j, ] <- treated_pos[sorted_treated]
  }
  
  match_index_control <- rep(NA_integer_, n_total * m)
  match_index_treated <- rep(NA_integer_, n_total * m)
  
  for (i in 1:n_total) {
    start_idx <- (i - 1) * m + 1
    end_idx <- i * m
    if (W[i] == 0) {
      match_index_control[start_idx:end_idx] <- rep(i, m)  
    } else {
      treated_idx_in_group <- which(treated_pos == i)  
      match_index_control[start_idx:end_idx] <- treated_matches[treated_idx_in_group, ]
    }
  }
  
  for (i in 1:n_total) {
    start_idx <- (i - 1) * m + 1
    end_idx <- i * m
    if (W[i] == 1) {
      match_index_treated[start_idx:end_idx] <- rep(i, m)  
    } else {
      control_idx_in_group <- which(control_pos == i)  
      match_index_treated[start_idx:end_idx] <- control_matches[control_idx_in_group, ]
    }
  }
  
  K_M <- rep(0, n_total)
  for (i in 1:n_total) {
    if (W[i] == 1) {
      matched_index <- match_index_control[(m * i - m + 1):(m * i)]
    } else {
      matched_index <- match_index_treated[(m * i - m + 1):(m * i)]
    }
    
    matched_index <- matched_index[!is.na(matched_index)]
    if (length(matched_index) == 0) next
    
    K_M[matched_index] <- K_M[matched_index] + 1
  }
  

  w_asdm <- 1 + K_M / m

  valid_treated_pos <- treated_pos[!dropped_treated]
  valid_control_pos <- control_pos[!dropped_control]

  
  orig_var_trt <- apply(X[treated_pos, , drop = FALSE], 2, stats::var, na.rm = TRUE)
  orig_var_cnt <- apply(X[control_pos, , drop = FALSE], 2, stats::var, na.rm = TRUE)
  asdm.vec <- rep(NA_real_, ncol(X))
  names(asdm.vec) <- colnames(X)
  
  for (j in seq_len(ncol(X))) {
    sum_wt_trt <- sum(w_asdm[valid_treated_pos], na.rm = TRUE)
    sum_wt_cnt <- sum(w_asdm[valid_control_pos], na.rm = TRUE)
    
    mean_trt <- sum(X[valid_treated_pos, j] * w_asdm[valid_treated_pos], na.rm = TRUE) / sum_wt_trt
    mean_cnt <- sum(X[valid_control_pos, j] * w_asdm[valid_control_pos], na.rm = TRUE) / sum_wt_cnt
    denom <- sqrt((orig_var_trt[j] + orig_var_cnt[j]) / 2)
    asdm.vec[j] <- if (is.finite(denom) && denom > 0) {
      abs(mean_trt - mean_cnt) / denom
    } else {
      NA_real_
    }
  }
  
  match.res <- list(
    index.control = match_index_control,
    index.treated = match_index_treated,
    K_M = K_M
  )
  
  list(
    weight = weight,
    asdm.vec = asdm.vec,
    K_M = K_M,
    match.res = match.res,
    treated_matches = treated_matches,
    control_matches = control_matches
  )
}




find_optimal_weight <- function(candidate_weights,
                                asdm_caliper = 0.15,
                                W,
                                X,
                                Y,
                                S,
                                prop_score = NULL,
                                m = 2,
                                caliper = 3,
                                caliper_type = "DAPS",
                                coord_dist = FALSE) {
  candidate_weights <- sort(candidate_weights)
  

  all_results <- vector("list", 0)
  
  for (w in candidate_weights) {
    cat("Processing weight:", w, "\n")
    
    result <- calculate_asdm_for_weight(
      weight = w, W = W, X = X, Y = Y, S = S, prop_score = prop_score,
      m = m, caliper = caliper,
      caliper_type = caliper_type,
      coord_dist = coord_dist)
    
    all_asdm_below_caliper <- all(result$asdm.vec < asdm_caliper, na.rm = FALSE)
    all_asdm_below_caliper <- isTRUE(all_asdm_below_caliper)
    

    if (all_asdm_below_caliper) {
      finite_asdm <- result$asdm.vec[is.finite(result$asdm.vec)]
      
      return(list(
        optimal_weight = w,
        asdm.vec = result$asdm.vec,
        K_M = result$K_M,
        max_asdm = if (length(finite_asdm) > 0) max(finite_asdm) else NA_real_,
        mean_asdm = if (length(finite_asdm) > 0) mean(finite_asdm) else NA_real_,
        selection_rule = "first_weight_satisfying_asdm_caliper"
      ))
    }
    

    finite_asdm <- result$asdm.vec[is.finite(result$asdm.vec)]
    
    if (length(finite_asdm) > 0) {
      all_results[[length(all_results) + 1]] <- list(
        weight = w,
        max_asdm = max(finite_asdm),
        mean_asdm = mean(finite_asdm),
        asdm.vec = result$asdm.vec,
        K_M = result$K_M
      )
    }
  }
  

  if (length(all_results) == 0) {
    return(list(
      optimal_weight = NA_real_,
      asdm.vec = NA,
      K_M = NA,
      max_asdm = NA_real_,
      mean_asdm = NA_real_,
      selection_rule = "no_valid_weight"
    ))
  }
  
  max_vals <- sapply(all_results, function(x) x$max_asdm)
  min_max_asdm <- min(max_vals)
  idx_stage1 <- which(max_vals == min_max_asdm)
  
  if (length(idx_stage1) > 1) {
    mean_vals <- sapply(all_results[idx_stage1], function(x) x$mean_asdm)
    min_mean_asdm <- min(mean_vals)
    idx_stage2 <- idx_stage1[which(mean_vals == min_mean_asdm)]
  } else {
    idx_stage2 <- idx_stage1
  }
  
  if (length(idx_stage2) > 1) {
    weight_vals <- sapply(all_results[idx_stage2], function(x) x$weight)
    best_idx <- idx_stage2[which.min(weight_vals)]
  } else {
    best_idx <- idx_stage2
  }
  cat("Using the compromise weight:", all_results[[best_idx]]$weight, "\n")
  
  return(list(
    optimal_weight = all_results[[best_idx]]$weight,
    asdm.vec = all_results[[best_idx]]$asdm.vec,
    K_M = all_results[[best_idx]]$K_M,
    max_asdm = all_results[[best_idx]]$max_asdm,
    mean_asdm = all_results[[best_idx]]$mean_asdm,
    selection_rule = "min_max_asdm_then_mean_asdm_then_smallest_weight"
  ))
}




MBlearning <- function(X,
                       W,
                       Y,
                       S,
                       prop_score = NULL,
                       m = 1,
                       biascorrect = TRUE,
                       biascorrect.lasso = FALSE,
                       tree_depth = 2,
                       coord_dist = FALSE,
                       caliper_type = "DAPS",
                       weight,
                       caliper = 3) {
  W <- as.integer(W)
  X <- as.matrix(X)
  S <- as.matrix(S)
  prop_score <- as.numeric(prop_score)
  
  treated_pos <- which(W == 1)
  control_pos <- which(W == 0)
  n_treated <- length(treated_pos)
  n_control <- length(control_pos)
  n_total <- length(W)
  
  treated <- data.frame(
    Longitude = S[treated_pos, 1],
    Latitude = S[treated_pos, 2],
    prop.scores = prop_score[treated_pos]
  )
  
  control <- data.frame(
    Longitude = S[control_pos, 1],
    Latitude = S[control_pos, 2],
    prop.scores = prop_score[control_pos]
  )
  

  if (coord_dist) {
    dist.mat <- fields::rdist.earth(cbind(treated$Longitude, 
                                          treated$Latitude), cbind(control$Longitude, control$Latitude))
  }else {
    dist.mat <- fields::rdist(cbind(treated$Longitude, treated$Latitude), 
                              cbind(control$Longitude, control$Latitude))
  }
  
  stand.dist.mat <- StandDist(dist.mat)
  
  treatment_indicator <- c(rep(1, nrow(treated)), rep(0, nrow(control)))
  propensity_scores <- c(treated$prop.scores, control$prop.scores)
  ps.diff <- as.matrix(optmatch::match_on(treatment_indicator ~ 
                                            propensity_scores, method = "euclidean"))
  
  dapscore <- (1 - weight) * stand.dist.mat + weight * ps.diff
  
  if (caliper_type == "DAPS") {
    M <- dapscore + optmatch::caliper(dapscore, caliper * sd(dapscore))
  } else if (caliper_type == "PS") {
    M <- dapscore + optmatch::caliper(ps.diff, caliper * sd(propensity_scores))
  }
  
  
  M <- as.matrix(M)
  
  treated_matches <- matrix(NA, nrow = n_treated, ncol = m)  
  dropped_treated <- logical(n_treated) 
  for (i in 1:n_treated) {
    dist_i <- M[i, ]
    ok <- which(is.finite(dist_i))                           
    if (length(ok) < m) {                                    
      dropped_treated[i] <- TRUE                              
      next  
    }
    sorted_control <- ok[order(dist_i[ok])][1:m]             
    treated_matches[i, ] <- control_pos[sorted_control]      
  }
  
  
  control_matches <- matrix(NA, nrow = n_control, ncol = m)  
  dropped_control <- logical(n_control)                      
  M_t <- t(M)  
  for (j in 1:n_control) {
    dist_j <- M_t[j, ]
    ok <- which(is.finite(dist_j))                           
    if (length(ok) < m) {                                    
      dropped_control[j] <- TRUE                              
      next
    }
    sorted_treated <- ok[order(dist_j[ok])][1:m]             
    control_matches[j, ] <- treated_pos[sorted_treated]      
  }
 
  match_index_control <- rep(NA_integer_, n_total * m)
  for (i in 1:n_total) {
    start_idx <- (i - 1) * m + 1
    end_idx <- i * m
    if (W[i] == 0) {
      match_index_control[start_idx:end_idx] <- rep(i, m)  
    } else {
      treated_idx_in_group <- which(treated_pos == i)  
      match_index_control[start_idx:end_idx] <- treated_matches[treated_idx_in_group, ]
    }
  }
  
  
  
  match_index_treated <- rep(NA_integer_, n_total * m)
  for (i in 1:n_total) {
    start_idx <- (i - 1) * m + 1
    end_idx <- i * m
    if (W[i] == 1) {
      match_index_treated[start_idx:end_idx] <- rep(i, m)  
    } else {
      control_idx_in_group <- which(control_pos == i)  
      match_index_treated[start_idx:end_idx] <- control_matches[control_idx_in_group, ]
    }
  }
  
  match.res=list()
  match.res[["index.control"]]=match_index_control
  match.res[["index.treated"]]=match_index_treated
  
  K_M <- rep(0, n_total)  
  for (i in 1:n_total) {
    if (W[i] == 1) {
      matched_index <- match_index_control[(m * i - m + 1):(m * i)]
    } else {
      matched_index <- match_index_treated[(m * i - m + 1):(m * i)]
    }
    
    matched_index <- matched_index[!is.na(matched_index)]
    if (length(matched_index) == 0) next
    
    K_M[matched_index] <- K_M[matched_index] + 1
  }
 
  
  if (!biascorrect) {
   
    Y0_simple = Y0_bc = Y1_simple = Y1_bc = rep(NA_real_, n_total)
    for (i in 1:dim(X)[1]) {
      if (W[i] == 1) {
        Y1_simple[i] = Y1_bc[i] = Y[i]
        matched_index = match.res[["index.control"]][(m*i-m+1):(m*i)]
        matched_index <- matched_index[!is.na(matched_index)]
        if (length(matched_index) == 0) {
          Y0_simple[i] <- NA
          Y0_bc[i] <- NA
          next
        }
        Y0_simple[i] = mean(Y[matched_index], na.rm = TRUE)
        Y0_bc[i] = Y0_simple[i]
      } else {
        Y0_simple[i] = Y0_bc[i] = Y[i]
        matched_index = match.res[["index.treated"]][(m*i-m+1):(m*i)]
        matched_index <- matched_index[!is.na(matched_index)]
        if (length(matched_index) == 0) {
          Y1_simple[i] <- NA
          Y1_bc[i] <- NA
          next
        }
        Y1_simple[i] = mean(Y[matched_index], na.rm = TRUE)
        Y1_bc[i] = Y1_simple[i]
      }
    }
  } else if (!biascorrect.lasso) {
    
    X <- as.matrix(X)
    storage.mode(X) <- "double"
    x_names <- colnames(X)
    
    data_gam <- data.frame(Y = Y, X, lon = S[, 1], lat = S[, 2])
    colnames(data_gam) <- c("Y", x_names, "lon", "lat")
    
    idx1 <- which(W==1 & K_M>0)
    idx0 <- which(W==0 & K_M>0)
    
    
    gam_formula <- as.formula(
      paste0(" Y ~ ", paste(colnames(X), collapse = " + "), " + s(lon, lat, bs='tp',k=8)")
    )
    
    
    
    reg_Trt1_gam <- mgcv::gam(gam_formula, data = data_gam[idx1, , drop=FALSE],
                              weights = K_M[idx1], method = "REML")
    reg_Trt0_gam <- mgcv::gam(gam_formula, data = data_gam[idx0, , drop=FALSE],
                              weights = K_M[idx0], method = "REML")
    
    reg_beta1 <- coef(reg_Trt1_gam)[colnames(X)]
    reg_beta0 <- coef(reg_Trt0_gam)[colnames(X)]
    
    
    terms1 <- predict(reg_Trt1_gam, newdata = data_gam, type = "terms")
    terms0 <- predict(reg_Trt0_gam, newdata = data_gam, type = "terms")
    
    
    sc1 <- grep("^s\\(lon,\\s*lat\\)$", colnames(terms1))
    sc0 <- grep("^s\\(lon,\\s*lat\\)$", colnames(terms0))
    
    if (length(sc1) != 1 || length(sc0) != 1) stop("GAM smooth term not found in predict(type='terms').")
    
    g1_all <- terms1[, sc1]
    g0_all <- terms0[, sc0]
    
    
    Y0_simple = Y0_bc = Y1_simple = Y1_bc = rep(NA_real_, n_total)
    for(i in 1:dim(X)[1]){
      if(W[i]==1){
        Y1_simple[i] = Y1_bc[i] = Y[i]
        matched_index = match.res[["index.control"]][(m*i-m+1):(m*i)]
        
        matched_index <- matched_index[!is.na(matched_index)]     
        if (length(matched_index) == 0) {                         
          Y0_simple[i] <- NA                                      
          Y0_bc[i] <- NA
          next
        }
        
        Y0_simple[i] = mean(Y[matched_index], na.rm = TRUE)
        Y0_bc[i] = Y0_simple[i] + X[i,]%*%reg_beta0 - apply(matrix(X[matched_index,],ncol = dim(X)[2]), 2, mean)%*%reg_beta0+
          g0_all[i] - mean(g0_all[matched_index], na.rm = TRUE)  
        ## --------------------------------
      }else{
        Y0_simple[i] = Y0_bc[i] = Y[i]
        matched_index = match.res[["index.treated"]][(m*i-m+1):(m*i)]
        matched_index <- matched_index[!is.na(matched_index)]
        if (length(matched_index) == 0) {
          Y1_simple[i] <- NA
          Y1_bc[i] <- NA
          next
        }
        Y1_simple[i] = mean(Y[matched_index], na.rm = TRUE)
        Y1_bc[i] = Y1_simple[i] + X[i,]%*%reg_beta1 - apply(matrix(X[matched_index,],ncol = dim(X)[2]), 2, mean)%*%reg_beta1+
          g1_all[i] - mean(g1_all[matched_index], na.rm = TRUE)
      }  
    }
  }
  
  else {
    
    
    X <- as.matrix(X)
    storage.mode(X) <- "double"
    
    idx1 <- which(W == 1 & K_M > 0)
    idx0 <- which(W == 0 & K_M > 0)
    
    X_interaction <- matrix(
      NA,
      nrow = nrow(X),
      ncol = choose(ncol(X), 2)
    )
    
    ind <- 1
    for (a in 1:(ncol(X) - 1)) {
      for (b in ((a + 1):ncol(X))) {
        X_interaction[, ind] <- X[, a] * X[, b]
        ind <- ind + 1
      }
    }
    
    Phi_X <- cbind(X, X^2, X_interaction)
    storage.mode(Phi_X) <- "double"
    
    
    fit_phi_spatial_gamlasso <- function(idx, max_iter = 10, tol = 1e-4, k_space = 4) {
      
      Phi_train <- Phi_X[idx, , drop = FALSE]
      y_train <- Y[idx]
      w_train <- K_M[idx]
      
      data_train <- data.frame(
        Y = y_train,
        lon = S[idx, 1],
        lat = S[idx, 2]
      )
      
      
      nfolds_use <- min(10, length(idx))
      foldid_use <- sample(rep(seq_len(nfolds_use), length.out = length(idx)))
      
      init_glm <- glmnet::cv.glmnet(
        x = Phi_train,
        y = y_train,
        weights = w_train,
        alpha = 1,
        intercept = TRUE,
        nfolds = nfolds_use,
        foldid = foldid_use,
        lambda.min.ratio = 1e-3
      )
      
      beta_init <- as.numeric(
        coef(init_glm, s = "lambda.min")
      )[-1]
      
      eta_phi <- as.numeric(Phi_train %*% beta_init)
      
      old_full_fit <- rep(Inf, length(idx))
      
      gam_formula_spatial <- as.formula(
        paste0(
          "Y ~ offset(.offset_phi) + ",
          "s(lon, lat, bs = 'tp', k = ", k_space, ")"
        )
      )
      
      
      for (iter in seq_len(max_iter)) {
        
        
        data_train$.offset_phi <- eta_phi
        
        gam_fit <- mgcv::gam(
          gam_formula_spatial,
          data = data_train,
          weights = w_train,
          method = "REML"
        )
        
        
        eta_spatial <- as.numeric(
          predict(gam_fit, newdata = data_train, type = "link")
        ) - eta_phi
        
        
        glm_fit <- glmnet::cv.glmnet(
          x = Phi_train,
          y = y_train,
          weights = w_train,
          offset = eta_spatial,
          alpha = 1,
          intercept = FALSE,
          nfolds = nfolds_use,
          foldid = foldid_use,
          lambda.min.ratio = 1e-3
        )
        
        beta_hat <- as.numeric(
          coef(glm_fit, s = "lambda.min")
        )[-1]
        
        eta_phi_new <- as.numeric(Phi_train %*% beta_hat)
        full_fit_new <- eta_spatial + eta_phi_new
        
        
        if (
          iter > 1 &&
          mean((full_fit_new - old_full_fit)^2, na.rm = TRUE) < tol
        ) {
          eta_phi <- eta_phi_new
          break
        }
        
        eta_phi <- eta_phi_new
        old_full_fit <- full_fit_new
      }
      
      return(list(
        gam_fit = gam_fit,
        glm_fit = glm_fit,
        beta_hat = beta_hat
      ))
    }
    
    
    
    predict_phi_spatial_gamlasso <- function(fit_obj) {
      
      eta_phi_all <- as.numeric(Phi_X %*% fit_obj$beta_hat)
      
      data_predict <- data.frame(
        lon = S[, 1],
        lat = S[, 2]
      )
      data_predict$.offset_phi <- eta_phi_all
      
      mu_all <- as.numeric(
        predict(fit_obj$gam_fit, newdata = data_predict, type = "link")
      )
      
      return(as.numeric(mu_all))
    }
    
    
    
    reg_Trt1_gamlasso <- fit_phi_spatial_gamlasso(idx1)
    reg_Trt0_gamlasso <- fit_phi_spatial_gamlasso(idx0)
    
    mu1_all <- predict_phi_spatial_gamlasso(reg_Trt1_gamlasso)
    mu0_all <- predict_phi_spatial_gamlasso(reg_Trt0_gamlasso)
    
    
    Y0_simple <- Y0_bc <- Y1_simple <- Y1_bc <- rep(NA_real_, n_total)
    
    for (i in 1:nrow(X)) {
      
      if (W[i] == 1) {
        
        Y1_simple[i] <- Y1_bc[i] <- Y[i]
        
        matched_index <- match.res[["index.control"]][
          (m * i - m + 1):(m * i)
        ]
        
        matched_index <- matched_index[!is.na(matched_index)]
        
        if (length(matched_index) == 0) {
          Y0_simple[i] <- NA
          Y0_bc[i] <- NA
          next
        }
        
        Y0_simple[i] <- mean(Y[matched_index], na.rm = TRUE)
        Y0_bc[i] <- Y0_simple[i] +
          mu0_all[i] -
          mean(mu0_all[matched_index], na.rm = TRUE)
        
      } else {
        
        Y0_simple[i] <- Y0_bc[i] <- Y[i]
        
        matched_index <- match.res[["index.treated"]][
          (m * i - m + 1):(m * i)
        ]
        
        matched_index <- matched_index[!is.na(matched_index)]
        
        if (length(matched_index) == 0) {
          Y1_simple[i] <- NA
          Y1_bc[i] <- NA
          next
        }
        
        Y1_simple[i] <- mean(Y[matched_index], na.rm = TRUE)
        Y1_bc[i] <- Y1_simple[i] +
          mu1_all[i] -
          mean(mu1_all[matched_index], na.rm = TRUE)
      }
    }
  }
  
  
  contrastY_simple = Y1_simple-Y0_simple
  contrastY_bc = Y1_bc-Y0_bc
  
  
  if(biascorrect){
    Gammai=contrastY_bc
  }else{
    Gammai=contrastY_simple
  }
  Gammai <- as.numeric(Gammai)
  if (!any(is.finite(Gammai))) {
    Gammai <- as.numeric(contrastY_simple)
  }
  
  
  valid_idx <- is.finite(Gammai)  
  
  X_valid <- X[valid_idx, , drop = FALSE]  
  Gamma_valid <- Gammai[valid_idx]         
  
  X_policy_valid <- cbind(X_valid, S[valid_idx, , drop = FALSE])
  
  
  if (tree_depth < 3) {
    matching_tr <- policy_tree(X_policy_valid, Gamma = cbind(-Gamma_valid, Gamma_valid), depth = tree_depth)
  } else {
    matching_tr <- hybrid_policy_tree(X = X_policy_valid, Gamma = cbind(-Gamma_valid, Gamma_valid), depth = tree_depth)
  }
  
  return(list(obj=matching_tr,contrastY_simple=contrastY_simple,contrastY_bc=contrastY_bc))
}





MBlearning_non<-function(X,W,Y,m = 1,biascorrect=TRUE,tree_depth=2,biascorrect.lasso=FALSE){
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  W <- as.integer(W)
  Y <- as.numeric(Y)
  
  match.res  <- Match(Y=Y, Tr=W, X=X, M=m, ties = F,BiasAdjust = F, Weight = 2, estimand='ATE')
  
  
  K_M=rep(0,dim(X)[1]) 
  for(i in 1:dim(X)[1]){
    if(W[i]==1){
      matched_index = match.res[["index.control"]][(m*i-m+1):(m*i)]
      K_M[matched_index]=K_M[matched_index]+1 
    } else {
      matched_index = match.res[["index.treated"]][(m*i-m+1):(m*i)]
      K_M[matched_index]=K_M[matched_index]+1 
    }
  }
  
  if (!biascorrect) {
    
    Y0_simple = Y0_bc = Y1_simple = Y1_bc = rep(NA_real_, dim(X)[1])
    for(i in 1:dim(X)[1]){
      if(W[i]==1){
        Y1_simple[i] = Y1_bc[i] = Y[i]
        matched_index = match.res[["index.control"]][(m*i-m+1):(m*i)]
        matched_index <- matched_index[!is.na(matched_index)]
        if (length(matched_index) == 0) {
          Y0_simple[i] <- NA
          Y0_bc[i] <- NA
          next
        }
        Y0_simple[i] = mean(Y[matched_index], na.rm = TRUE)
        Y0_bc[i] = Y0_simple[i]
      }else{
        Y0_simple[i] = Y0_bc[i] = Y[i]
        matched_index = match.res[["index.treated"]][(m*i-m+1):(m*i)]
        matched_index <- matched_index[!is.na(matched_index)]
        if (length(matched_index) == 0) {
          Y1_simple[i] <- NA
          Y1_bc[i] <- NA
          next
        }
        Y1_simple[i] = mean(Y[matched_index], na.rm = TRUE)
        Y1_bc[i] = Y1_simple[i]
      }
    }
  } else if(!biascorrect.lasso){
    
    reg_Trt1=lm(Y[which(W==1)]~X[which(W==1),],weights = K_M[which(W==1)])
    reg_Trt0=lm(Y[which(W==0)]~X[which(W==0),],weights = K_M[which(W==0)])
    
    reg_beta1=coefficients(reg_Trt1)[-1]
    reg_beta0=coefficients(reg_Trt0)[-1]
    
    
    Y0_simple = Y0_bc = Y1_simple = Y1_bc = rep(0,dim(X)[1])
    for(i in 1:dim(X)[1]){
      if(W[i]==1){
        Y1_simple[i] = Y1_bc[i] = Y[i]
        matched_index = match.res[["index.control"]][(m*i-m+1):(m*i)]
        Y0_simple[i] = mean(Y[matched_index])
        Y0_bc[i] = Y0_simple[i] + X[i,]%*%reg_beta0 - apply(matrix(X[matched_index,],ncol = dim(X)[2]), 2, mean)%*%reg_beta0
      }else{
        Y0_simple[i] = Y0_bc[i] = Y[i]
        matched_index = match.res[["index.treated"]][(m*i-m+1):(m*i)]
        Y1_simple[i] = mean(Y[matched_index])
        Y1_bc[i] = Y1_simple[i] + X[i,]%*%reg_beta1 - apply(matrix(X[matched_index,],ncol = dim(X)[2]), 2, mean)%*%reg_beta1
      }
    }
  } else {
    
    
    X_interaction=matrix(NA,nrow=nrow(X),ncol = choose(ncol(X),2))
    ind=1
    for(i in 1:(ncol(X)-1)){
      for(j in ((i+1):ncol(X))){
        X_interaction[,ind]=X[,i]*X[,j]
        ind=ind+1
      }}
    
    DataX=cbind(X,X^2,X_interaction)
    
    reg_Trt1=cv.glmnet(x=DataX[which(W==1),], y = Y[which(W==1)],weights=K_M[which(W==1)],lambda.min.ratio=1e-3)
    reg_Trt0=cv.glmnet(x=DataX[which(W==0),], y = Y[which(W==0)],weights=K_M[which(W==0)],lambda.min.ratio=1e-3)
    
    reg_beta1=as.numeric(coef(reg_Trt1,s='lambda.min'))[-1]
    reg_beta0=as.numeric(coef(reg_Trt0,s='lambda.min'))[-1]
    
    
    Y0_simple = Y0_bc = Y1_simple = Y1_bc = rep(0,dim(X)[1])
    for(i in 1:dim(X)[1]){
      if(W[i]==1){
        Y1_simple[i] = Y1_bc[i] = Y[i]
        matched_index = match.res[["index.control"]][(m*i-m+1):(m*i)]
        Y0_simple[i] = mean(Y[matched_index])
        Y0_bc[i] = Y0_simple[i] + DataX[i,]%*%reg_beta0 - apply(matrix(DataX[matched_index,],ncol = dim(DataX)[2]), 2, mean)%*%reg_beta0
      }else{
        Y0_simple[i] = Y0_bc[i] = Y[i]
        matched_index = match.res[["index.treated"]][(m*i-m+1):(m*i)]
        Y1_simple[i] = mean(Y[matched_index])
        Y1_bc[i] = Y1_simple[i] + DataX[i,]%*%reg_beta1 - apply(matrix(DataX[matched_index,],ncol = dim(DataX)[2]), 2, mean)%*%reg_beta1
      }
    }
    
  }
  
  contrastY_simple = Y1_simple-Y0_simple
  contrastY_bc = Y1_bc-Y0_bc
  
  
  if(biascorrect){
    Gammai=contrastY_bc
  }else{
    Gammai=contrastY_simple
  }
  Gammai <- as.numeric(Gammai)
  if (!any(is.finite(Gammai))) {
    Gammai <- as.numeric(contrastY_simple)
  }
  
  valid_idx <- is.finite(Gammai)
  X_valid <- X[valid_idx, , drop = FALSE]
  Gamma_valid <- Gammai[valid_idx]
  
  
  if(tree_depth<3){
    matching_tr <- policy_tree(X_valid, Gamma = cbind(-Gamma_valid,Gamma_valid), depth = tree_depth)
  } else {
    matching_tr <- hybrid_policy_tree(X = X_valid,Gamma = cbind(-Gamma_valid,Gamma_valid),depth = tree_depth)
  }
  
  return(list(obj=matching_tr,contrastY_simple=contrastY_simple,contrastY_bc=contrastY_bc))
}











get_true_value <- function(policy, dat) {
  mean(ifelse(policy == 1, dat$Y1_true, dat$Y0_true))
}

evaluate_policy_sim <- function(policy, dat, opt_rule, opt_value) {
  value <- get_true_value(policy, dat)
  regret <- opt_value - value
  accuracy <- mean(policy == opt_rule)
  
  return(list(
    value = value,
    regret = regret,
    accuracy = accuracy
  ))
}





compareAllmethod_sim <- function(
    seed = 1,
    n_train = 800,
    n_test = 2000,
    nu = 1.46,
    range = 0.25,
    sigma_eps = 1,
    candidate_weights = seq(0.1, 1, by = 0.1),
    asdm_caliper = 0.1,
    caliper = 3,
    caliper_type = "DAPS",
    coord_dist = FALSE,
    tree_depth = 2,
    scenario="tree_X"
) {
  
  result <- matrix(0, 1, 15)
  colnames(result) <- c("MB-M1", "MB-M2", "MB-LR-M1",
                        "MB-LR-M2", "MB-LASSO-M1", "MB-LASSO-M2", "Policy Tree",
                        "MB-M1-non", "MB-M2-non", "MB-LR-M1-non",
                        "MB-LR-M2-non", "MB-LASSO-M1-non", "MB-LASSO-M2-non","Qlearning-Para","Qlearning-RF")
  rownames(result) <- "value"
  
  acc_result <- matrix(0, 1, 15)
  colnames(acc_result) <- colnames(result)
  rownames(acc_result) <- "accuracy"
  
  regret_result <- matrix(0, 1, 15)
  colnames(regret_result) <- colnames(result)
  rownames(regret_result) <- "regret"
  
  

  data.obj <- generate_train_test_data(
    n_train = n_train,
    n_test = n_test,
    nu = nu,
    range = range,
    sigma_eps = sigma_eps,
    seed = seed,
    scenario = scenario
  )
  
  train.dat <- data.obj$train.dat
  test.dat <- data.obj$test.dat
  opt_rule <- data.obj$opt_rule
  oracle_value <- data.obj$opt_value
  #############################
  
  X_train <- as.matrix(train.dat[, c("X1", "X2", "X3", "X4")])
  W_train <- train.dat$W
  Y_train <- train.dat$Y
  S_train <- as.matrix(train.dat[, c("Longitude", "Latitude")])
  prop_score_train <- train.dat$prop_score
  
  X_test <- as.matrix(test.dat[, c("X1", "X2", "X3", "X4")])
  S_test <- as.matrix(test.dat[, c("Longitude", "Latitude")])
  X_policy_test_spatial <- cbind(X_test, S_test)
  
  # -------
  optimal_result1 <- find_optimal_weight(
    candidate_weights = candidate_weights,
    asdm_caliper = asdm_caliper,
    W = W_train,
    X = X_train,
    Y = Y_train,
    S = S_train,
    prop_score = prop_score_train,
    m = 1,
    caliper = caliper,
    caliper_type = caliper_type,
    coord_dist = coord_dist
  )
  
  optimal_result2 <- find_optimal_weight(
    candidate_weights = candidate_weights,
    asdm_caliper = asdm_caliper,
    W = W_train,
    X = X_train,
    Y = Y_train,
    S = S_train,
    prop_score = prop_score_train,
    m = 2,
    caliper = caliper,
    caliper_type = caliper_type,
    coord_dist = coord_dist
  )
  
  # -------- MB-M1 --------
  fit1 <- MBlearning(
    X = X_train, W = W_train, Y = Y_train, S = S_train, prop_score = prop_score_train,
    m = 1, biascorrect = FALSE, tree_depth = tree_depth,
    coord_dist = coord_dist, caliper_type = caliper_type,
    weight = optimal_result1$optimal_weight, caliper = caliper
  )
  policy1 <- predict(fit1$obj, X_policy_test_spatial) - 1
  eval1 <- evaluate_policy_sim(policy1, test.dat, opt_rule, oracle_value)
  
  result[1, 1] <- eval1$value
  acc_result[1, 1] <- eval1$accuracy
  regret_result[1, 1] <- eval1$regret
  
  # -------- MB-M2 --------
  fit2 <- MBlearning(
    X = X_train, W = W_train, Y = Y_train, S = S_train, prop_score = prop_score_train,
    m = 2, biascorrect = FALSE, tree_depth = tree_depth,
    coord_dist = coord_dist, caliper_type = caliper_type,
    weight = optimal_result2$optimal_weight, caliper = caliper
  )
  policy2 <- predict(fit2$obj, X_policy_test_spatial) - 1
  eval2 <- evaluate_policy_sim(policy2, test.dat, opt_rule, oracle_value)
  
  result[1, 2] <- eval2$value
  acc_result[1, 2] <- eval2$accuracy
  regret_result[1, 2] <- eval2$regret
  
  # -------- MB-LR-M1 --------
  fit3 <- MBlearning(
    X = X_train, W = W_train, Y = Y_train, S = S_train, prop_score = prop_score_train,
    m = 1, biascorrect = TRUE, biascorrect.lasso = FALSE,
    tree_depth = tree_depth,
    coord_dist = coord_dist, caliper_type = caliper_type,
    weight = optimal_result1$optimal_weight, caliper = caliper
  )
  policy3 <- predict(fit3$obj, X_policy_test_spatial) - 1
  eval3 <- evaluate_policy_sim(policy3, test.dat, opt_rule, oracle_value)
  
  result[1, 3] <- eval3$value
  acc_result[1, 3] <- eval3$accuracy
  regret_result[1, 3] <- eval3$regret
  
  # -------- MB-LR-M2 --------
  fit4 <- MBlearning(
    X = X_train, W = W_train, Y = Y_train, S = S_train, prop_score = prop_score_train,
    m = 2, biascorrect = TRUE, biascorrect.lasso = FALSE,
    tree_depth = tree_depth,
    coord_dist = coord_dist, caliper_type = caliper_type,
    weight = optimal_result2$optimal_weight, caliper = caliper
  )
  policy4 <- predict(fit4$obj, X_policy_test_spatial) - 1
  eval4 <- evaluate_policy_sim(policy4, test.dat, opt_rule, oracle_value)
  
  result[1, 4] <- eval4$value
  acc_result[1, 4] <- eval4$accuracy
  regret_result[1, 4] <- eval4$regret
  
  # -------- MB-LASSO-M1 --------
  fit5 <- MBlearning(
    X = X_train, W = W_train, Y = Y_train, S = S_train, prop_score = prop_score_train,
    m = 1, biascorrect = TRUE, biascorrect.lasso = TRUE,
    tree_depth = tree_depth,
    coord_dist = coord_dist, caliper_type = caliper_type,
    weight = optimal_result1$optimal_weight, caliper = caliper
  )
  policy5 <- predict(fit5$obj, X_policy_test_spatial) - 1
  eval5 <- evaluate_policy_sim(policy5, test.dat, opt_rule, oracle_value)
  
  result[1, 5] <- eval5$value
  acc_result[1, 5] <- eval5$accuracy
  regret_result[1, 5] <- eval5$regret
  
  # -------- MB-LASSO-M2 --------
  fit6 <- MBlearning(
    X = X_train, W = W_train, Y = Y_train, S = S_train, prop_score = prop_score_train,
    m = 2, biascorrect = TRUE, biascorrect.lasso = TRUE,
    tree_depth = tree_depth,
    coord_dist = coord_dist, caliper_type = caliper_type,
    weight = optimal_result2$optimal_weight, caliper = caliper
  )
  policy6 <- predict(fit6$obj, X_policy_test_spatial) - 1
  eval6 <- evaluate_policy_sim(policy6, test.dat, opt_rule, oracle_value)
  
  result[1, 6] <- eval6$value
  acc_result[1, 6] <- eval6$accuracy
  regret_result[1, 6] <- eval6$regret
  
  # -------- Policy Tree --------
  c.forest <- causal_forest(X = X_train, Y = Y_train, W = W_train)
  dr.scores <- double_robust_scores(c.forest)
  
  if (tree_depth < 3) {
    fitted_tree <- policy_tree(X = X_train, Gamma = dr.scores, depth = tree_depth)
  } else {
    fitted_tree <- hybrid_policy_tree(X = X_train, Gamma = dr.scores, depth = tree_depth)
  }
  
  policy7 <- predict(fitted_tree, X_test) - 1
  eval7 <- evaluate_policy_sim(policy7, test.dat, opt_rule, oracle_value)
  
  result[1, 7] <- eval7$value
  acc_result[1, 7] <- eval7$accuracy
  regret_result[1, 7] <- eval7$regret
  
  # -------- MB-M1-non --------
  fit8 <- MBlearning_non(
    X = X_train, W = W_train, Y = Y_train,
    m = 1, biascorrect = FALSE, tree_depth = tree_depth)
  policy8 <- predict(fit8$obj, X_test) - 1
  eval8 <- evaluate_policy_sim(policy8, test.dat, opt_rule, oracle_value)
  
  result[1, 8] <- eval8$value
  acc_result[1, 8] <- eval8$accuracy
  regret_result[1, 8] <- eval8$regret
  
  # -------- MB-M2-non --------
  fit9 <- MBlearning_non(
    X = X_train, W = W_train, Y = Y_train,
    m = 2, biascorrect = FALSE, tree_depth = tree_depth)
  policy9 <- predict(fit9$obj, X_test) - 1
  eval9 <- evaluate_policy_sim(policy9, test.dat, opt_rule, oracle_value)
  
  result[1, 9] <- eval9$value
  acc_result[1, 9] <- eval9$accuracy
  regret_result[1, 9] <- eval9$regret
  
  # -------- MB-LR-M1-non --------
  fit10 <- MBlearning_non(
    X = X_train, W = W_train, Y = Y_train,
    m = 1, biascorrect = TRUE, biascorrect.lasso = FALSE,
    tree_depth = tree_depth)
  policy10 <- predict(fit10$obj, X_test) - 1
  eval10 <- evaluate_policy_sim(policy10, test.dat, opt_rule, oracle_value)
  
  result[1, 10] <- eval10$value
  acc_result[1, 10] <- eval10$accuracy
  regret_result[1, 10] <- eval10$regret
  
  # -------- MB-LR-M2-non --------
  fit11 <- MBlearning_non(
    X = X_train, W = W_train, Y = Y_train,
    m = 2, biascorrect = TRUE, biascorrect.lasso = FALSE,
    tree_depth = tree_depth)
  policy11 <- predict(fit11$obj, X_test) - 1
  eval11 <- evaluate_policy_sim(policy11, test.dat, opt_rule, oracle_value)
  
  result[1, 11] <- eval11$value
  acc_result[1, 11] <- eval11$accuracy
  regret_result[1, 11] <- eval11$regret
  
  # -------- MB-LASSO-M1-non --------
  fit12 <- MBlearning_non(
    X = X_train, W = W_train, Y = Y_train,
    m = 1, biascorrect = TRUE, biascorrect.lasso = TRUE,
    tree_depth = tree_depth)
  policy12 <- predict(fit12$obj, X_test) - 1
  eval12 <- evaluate_policy_sim(policy12, test.dat, opt_rule, oracle_value)
  
  result[1, 12] <- eval12$value
  acc_result[1, 12] <- eval12$accuracy
  regret_result[1, 12] <- eval12$regret
  
  # -------- MB-LASSO-M2-non --------
  fit13 <- MBlearning_non(
    X = X_train, W = W_train, Y = Y_train,
    m = 2, biascorrect = TRUE, biascorrect.lasso = TRUE,
    tree_depth = tree_depth)
  policy13 <- predict(fit13$obj, X_test) - 1
  eval13 <- evaluate_policy_sim(policy13, test.dat, opt_rule, oracle_value)
  
  result[1, 13] <- eval13$value
  acc_result[1, 13] <- eval13$accuracy
  regret_result[1, 13] <- eval13$regret

# -------- Q-learning (Para) --------
q_lasso_0 <- glmnet::cv.glmnet(
  x = model.matrix(
    ~ 0 + (X1 + X2 + X3 + X4)^2 +
      I(X1^2) + I(X2^2) + I(X3^2) + I(X4^2),
    train.dat[train.dat$W == 0, ]
  ),
  y = train.dat$Y[train.dat$W == 0],
  alpha = 1,
  lambda.min.ratio = 1e-3
)

q_lasso_1 <- glmnet::cv.glmnet(
  x = model.matrix(
    ~ 0 + (X1 + X2 + X3 + X4)^2 +
      I(X1^2) + I(X2^2) + I(X3^2) + I(X4^2),
    train.dat[train.dat$W == 1, ]
  ),
  y = train.dat$Y[train.dat$W == 1],
  alpha = 1,
  lambda.min.ratio = 1e-3
)

q0_hat_lasso <- as.numeric(predict(q_lasso_0,
  newx = model.matrix(~ 0 + (X1 + X2 + X3 + X4)^2 +I(X1^2) + I(X2^2) + I(X3^2) + I(X4^2),
    test.dat),s = "lambda.min"))

q1_hat_lasso <- as.numeric(predict(q_lasso_1,
  newx = model.matrix(~ 0 + (X1 + X2 + X3 + X4)^2 +I(X1^2) + I(X2^2) + I(X3^2) + I(X4^2),
    test.dat),s = "lambda.min"))

policy14 <- as.integer(q1_hat_lasso > q0_hat_lasso)
eval14 <- evaluate_policy_sim(policy14, test.dat, opt_rule, oracle_value)

result[1, 14] <- eval14$value
acc_result[1, 14] <- eval14$accuracy
regret_result[1, 14] <- eval14$regret
  
  # -------- Q-learning (RF) --------
  q_rf_variables <- c(
    "X1", "X2", "X3", "X4"
  )
  
  q_rf_0 <- grf::regression_forest(
    X = as.matrix(
      train.dat[train.dat$W == 0, q_rf_variables, drop = FALSE]
    ),
    Y = train.dat$Y[train.dat$W == 0],
    num.trees = 500,
    seed = 10000 + seed
  )
  
  q_rf_1 <- grf::regression_forest(
    X = as.matrix(
      train.dat[train.dat$W == 1, q_rf_variables, drop = FALSE]
    ),
    Y = train.dat$Y[train.dat$W == 1],
    num.trees = 500,
    seed = 20000 + seed
  )
  
  q0_hat_rf <- as.numeric(
    predict(
      q_rf_0,
      newdata = as.matrix(test.dat[, q_rf_variables, drop = FALSE])
    )$predictions
  )
  
  q1_hat_rf <- as.numeric(
    predict(
      q_rf_1,
      newdata = as.matrix(test.dat[, q_rf_variables, drop = FALSE])
    )$predictions
  )
  
  policy15 <- as.integer(q1_hat_rf > q0_hat_rf)
  eval15 <- evaluate_policy_sim(policy15, test.dat, opt_rule, oracle_value)
  
  result[1, 15] <- eval15$value
  acc_result[1, 15] <- eval15$accuracy
  regret_result[1, 15] <- eval15$regret
  
  return(list(
    value = result,
    accuracy = acc_result,
    regret = regret_result,
    optimal_weight_m1 = optimal_result1$optimal_weight,
    optimal_weight_m2 = optimal_result2$optimal_weight,
    oracle_value = oracle_value
  ))
}




run_simulation_reps <- function(
    replication = 100,
    n_train = 800,
    n_test = 2000,
    nu = 1.46,
    range = 0.25,
    sigma_eps = 1,
    asdm_caliper=0.1,
    tree_depth = 2,
    scenario="tree_X"
) {
  all_value <- list()
  all_acc <- list()
  all_regret <- list()
  all_oracle_value <- numeric(replication)
  
  for (i in 1:replication) {
    cat("Replication:", i, "\n")
    
    res_i <- compareAllmethod_sim(
      seed = i,
      n_train = n_train,
      n_test = n_test,
      nu = nu,
      range = range,
      sigma_eps = sigma_eps,
      asdm_caliper = asdm_caliper,
      tree_depth = tree_depth,scenario=scenario
    )
    
    all_value[[i]] <- res_i$value
    all_acc[[i]] <- res_i$accuracy
    all_regret[[i]] <- res_i$regret
    all_oracle_value[i] <- res_i$oracle_value
  }
  
  return(list(
    value = do.call(rbind, all_value),
    accuracy = do.call(rbind, all_acc),
    regret = do.call(rbind, all_regret),
    oracle_value = all_oracle_value
  ))
}



res_all_tree_X <- run_simulation_reps(replication = 100, n_train = 800, range=0.4, n_test = 2000,asdm_caliper=0.1,tree_depth = 2,scenario="tree_X")

apply(res_all_tree_X$value, 2, mean)
apply(res_all_tree_X$accuracy, 2, mean)
apply(res_all_tree_X$regret, 2, mean)
mean(res_all_tree_X$oracle_value)


res_all_nontree_X <- run_simulation_reps(replication = 100, n_train = 800,range=0.4, n_test = 2000,asdm_caliper=0.1,tree_depth = 2,scenario="nontree_X")

apply(res_all_nontree_X$value, 2, mean)
apply(res_all_nontree_X$accuracy, 2, mean)
apply(res_all_nontree_X$regret, 2, mean)
mean(res_all_nontree_X$oracle_value)


res_all_tree_XU <- run_simulation_reps(replication = 100, n_train = 800,range = 0.4, n_test = 2000,asdm_caliper = 0.1,tree_depth = 2,scenario = "tree_XU")

apply(res_all_tree_XU$value, 2, mean)
apply(res_all_tree_XU$accuracy, 2, mean)
apply(res_all_tree_XU$regret, 2, mean)
mean(res_all_tree_XU$oracle_value)


res_all_nontree_XU <- run_simulation_reps(replication = 100, n_train = 800,range = 0.4,n_test = 2000,asdm_caliper = 0.1,tree_depth = 2,scenario = "nontree_XU")

apply(res_all_nontree_XU$value, 2, mean)
apply(res_all_nontree_XU$accuracy, 2, mean)
apply(res_all_nontree_XU$regret, 2, mean)
mean(res_all_nontree_XU$oracle_value)


#################################
method_axis_labels <- c(
  "MB-M1" = "Spatial-MB-M1",
  "MB-M2" = "Spatial-MB-M2",
  "MB-LR-M1" = "Spatial-MB-Regression-M1",
  "MB-LR-M2" = "Spatial-MB-Regression-M2",
  "MB-LASSO-M1" = "Spatial-MB-Shrinkage-M1",
  "MB-LASSO-M2" = "Spatial-MB-Shrinkage-M2",
  "Policy Tree" = "Policy Tree",
  "MB-M1-non" = "Non-Spatial-MB-M1",
  "MB-M2-non" = "Non-Spatial-MB-M2",
  "MB-LR-M1-non" = "Non-Spatial-MB-Regression-M1",
  "MB-LR-M2-non" = "Non-Spatial-MB-Regression-M2",
  "MB-LASSO-M1-non" = "Non-Spatial-MB-Shrinkage-M1",
  "MB-LASSO-M2-non" = "Non-Spatial-MB-Shrinkage-M2"
)

# tree-X 
dt_tree <- as.data.table(res_all_tree_X$value)
dt_tree[, rep := seq_len(.N)]


dt_tree_long <- data.table::melt(
  dt_tree,
  id.vars = "rep",
  variable.name = "method",
  value.name = "value"
)

setDT(dt_tree_long)
dt_tree_long[, scenario := "Tree-X"]

# nontree-X
dt_nontree <- as.data.table(res_all_nontree_X$value)
dt_nontree[, rep := seq_len(.N)]

dt_nontree_long <- data.table::melt(
  dt_nontree,
  id.vars = "rep",
  variable.name = "method",
  value.name = "value"
)

setDT(dt_nontree_long)
dt_nontree_long[, scenario := "Non-tree-X"]


dt_all <- rbind(dt_tree_long, dt_nontree_long)
dt_all[, method := factor(method, levels = colnames(res_all_tree_X$value))]
dt_all[, scenario := factor(scenario, levels = c("Tree-X", "Non-tree-X"))]


p <- ggplot(dt_all, aes(x = method, y = value, fill = method)) +
  geom_boxplot(outlier.size = 0.6) +
  facet_wrap(~ scenario, nrow = 1) +
  scale_x_discrete(labels = method_axis_labels) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    strip.text = element_text(size = 14),
    axis.title = element_text(size = 12)
  ) +
  labs(
    x = "Method",
    y = "Value Function"
  )

print(p)

#################################XU
# tree-XU 
dt_tree_XU <- as.data.table(res_all_tree_XU$value)
dt_tree_XU[, rep := seq_len(.N)]


dt_tree_XU_long <- data.table::melt(
  dt_tree_XU,
  id.vars = "rep",
  variable.name = "method",
  value.name = "value"
)

setDT(dt_tree_XU_long)
dt_tree_XU_long[, scenario := "Tree-XU"]

# nontree-XU
dt_nontree_XU <- as.data.table(res_all_nontree_XU$value)
dt_nontree_XU[, rep := seq_len(.N)]

dt_nontree_XU_long <- data.table::melt(
  dt_nontree_XU,
  id.vars = "rep",
  variable.name = "method",
  value.name = "value"
)

setDT(dt_nontree_XU_long)
dt_nontree_XU_long[, scenario := "Non-tree-XU"]


dt_all_XU <- rbind(dt_tree_XU_long, dt_nontree_XU_long)
dt_all_XU[, method := factor(method, levels = colnames(res_all_tree_XU$value))]
dt_all_XU[, scenario := factor(scenario, levels = c("Tree-XU", "Non-tree-XU"))]


p_XU <- ggplot(dt_all_XU, aes(x = method, y = value, fill = method)) +
  geom_boxplot(outlier.size = 0.6) +
  facet_wrap(~ scenario, nrow = 1) +
  scale_x_discrete(labels = method_axis_labels) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    strip.text = element_text(size = 14),
    axis.title = element_text(size = 12)
  ) +
  labs(
    x = "Method",
    y = "Value Function"
  )

print(p_XU)



