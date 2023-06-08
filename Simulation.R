library(DiceKriging)
library(mined)

inv_logit <- function(x) {exp(x)/(1 + exp(x))}

Simulation_random <- function(B, n_train,
                              Covariate_model, Propensity_model = NULL, Outcome_model,
                              X_test, target = "ATE") {
  target_est <- numeric(B)

  pb <- txtProgressBar(min = 0, max = B, style = 3)
  for (b in 1:B) {
    X_train <- Covariate_model(n_train)
    if (is.null(Propensity_model)) {
      A_train <- rbinom(n_train, 1, 0.5)
    } else {
      Propensity_train <- inv_logit(apply(X_train, 1, Propensity_model))
      A_train <- as.numeric(runif(length(Propensity_train)) < Propensity_train)
    }
    Y_train <- mapply(Outcome_model, asplit(X_train, 1), A_train)

    idx_t <- A_train == 1
    idx_c <- A_train == 0
    fit_t <- km(~ 1,
                design = X_train[idx_t,],
                response = Y_train[idx_t],
                covtype = "matern5_2",
                control = list(pop.size = 50, trace = F),
                nugget.estim = T)
    fit_c <- km(~ 1,
                design = X_train[idx_c,],
                response = Y_train[idx_c],
                covtype = "matern5_2",
                control = list(pop.size = 50, trace = F),
                nugget.estim = T)

    pred_t <- predict(fit_t, X_test, type = "SK", se.compute = T)
    pred_c <- predict(fit_c, X_test, type = "SK", se.compute = T)

    if (target == "ATE") {
      target_est[b] <- mean(pred_t$mean - pred_c$mean)
    } else if (target == "ATTE") {
      # fit_propensity <- glm(A_train ~ .^2, family = "binomial", data = as_tibble(cbind(X_train, A_train)))
      # pred_propensity <- predict(fit_propensity, as_tibble(X_test), type = "response")
      pred_propensity <- inv_logit(apply(X_test, 1, Propensity_model))

      target_est[b] <- sum(pred_propensity*(pred_t$mean - pred_c$mean)/sum(pred_propensity))
    } else if (target == "ATO") {
      # fit_propensity <- glm(A_train ~ .^2, family = "binomial", data = as_tibble(cbind(X_train, A_train)))
      # pred_propensity <- predict(fit_propensity, as_tibble(X_test), type = "response")
      pred_propensity <- inv_logit(apply(X_test, 1, Propensity_model))

      target_est[b] <- sum(pred_propensity*(1 - pred_propensity)*(pred_t$mean - pred_c$mean))/sum(pred_propensity*(1 - pred_propensity))
    }

    setTxtProgressBar(pb, b)
  }
  close(pb)

  return(target_est)
}

Simulation_space_filling <- function(B, n_train, n_cand,
                                     Covariate_model, Propensity_model = NULL, Outcome_model,
                                     X_test, target = "ATE") {
  target_est <- numeric(B)

  pb <- txtProgressBar(min = 0, max = B, style = 3)
  for (b in 1:B) {
    X_cand <- Covariate_model(n_cand)
    X_train <- SelectMinED(X_cand, log(rep(1, n_cand)), n_train)$points

    if (is.null(Propensity_model)) {
      A_train <- rbinom(n_train, 1, 0.5)
    } else {
      Propensity_train <- inv_logit(apply(X_train, 1, Propensity_model))
      A_train <- as.numeric(runif(length(Propensity_train)) < Propensity_train)
    }
    Y_train <- mapply(Outcome_model, asplit(X_train, 1), A_train)

    idx_t <- A_train == 1
    idx_c <- A_train == 0
    fit_t <- km(~ 1,
                design = X_train[idx_t,],
                response = Y_train[idx_t],
                covtype = "matern5_2",
                control = list(pop.size = 50, trace = F),
                nugget.estim = T)
    fit_c <- km(~ 1,
                design = X_train[idx_c,],
                response = Y_train[idx_c],
                covtype = "matern5_2",
                control = list(pop.size = 50, trace = F),
                nugget.estim = T)

    pred_t <- predict(fit_t, X_test, type = "SK", se.compute = T)
    pred_c <- predict(fit_c, X_test, type = "SK", se.compute = T)

    if (target == "ATE") {
      target_est[b] <- mean(pred_t$mean - pred_c$mean)
    } else if (target == "ATTE") {
      fit_propensity <- glm(A_train ~ .^2, family = "binomial", data = as_tibble(cbind(X_train, A_train)))
      pred_propensity <- predict(fit_propensity, as_tibble(X_test), type = "response")
      # pred_propensity <- inv_logit(apply(X_test, 1, Propensity_model))

      target_est[b] <- sum(pred_propensity*(pred_t$mean - pred_c$mean)/sum(pred_propensity))
    } else if (target == "ATO") {
      fit_propensity <- glm(A_train ~ .^2, family = "binomial", data = as_tibble(cbind(X_train, A_train)))
      pred_propensity <- predict(fit_propensity, as_tibble(X_test), type = "response")
      # pred_propensity <- inv_logit(apply(X_test, 1, Propensity_model))

      target_est[b] <- sum(pred_propensity*(1 - pred_propensity)*(pred_t$mean - pred_c$mean))/sum(pred_propensity*(1 - pred_propensity))
    }

    setTxtProgressBar(pb, b)
  }
  close(pb)

  return(target_est)
}

Simulation1_VR <- function(B, n_train, n_init, n_batch = 10,
                           Covariate_model, Outcome_model,
                           X_test) {
  n_test <- nrow(X_test)
  target_est <- numeric(B)

  b <- 1
  pb <- txtProgressBar(min = 0, max = B, style = 3)
  while (b <= B) {
    try({
      X_cand <- Covariate_model(n_train)
      Y_cand <- cbind(mapply(Outcome_model, asplit(X_cand, 1), 0),
                      mapply(Outcome_model, asplit(X_cand, 1), 1))

      X_train <- X_cand[1:n_init, ]
      A_train <- as.numeric(1:n_init %in% sample(n_init, n_init %/% 2))
      Y_train <- Y_cand[cbind(1:n_init, A_train + 1)]

      idx_t <- A_train == 1
      idx_c <- A_train == 0
      fit_t <- km(~ 1,
                  design = X_train[idx_t,],
                  response = Y_train[idx_t],
                  covtype = "matern5_2",
                  control = list(pop.size = 50, trace = F),
                  nugget.estim = T)
      fit_c <- km(~ 1,
                  design = X_train[idx_c,],
                  response = Y_train[idx_c],
                  covtype = "matern5_2",
                  control = list(pop.size = 50, trace = F),
                  nugget.estim = T)

      while (length(A_train) < n_train) {
        pred_t <- predict(fit_t, rbind(X_test, X_cand[length(A_train) + 1,]), type = "SK", cov.compute = T)
        pred_c <- predict(fit_c, rbind(X_test, X_cand[length(A_train) + 1,]), type = "SK", cov.compute = T)

        Reduc_t <- sum(pred_t$cov[1:n_test, n_test + 1]^2)/pred_t$cov[n_test + 1, n_test + 1]
        Reduc_c <- sum(pred_c$cov[1:n_test, n_test + 1]^2)/pred_c$cov[n_test + 1, n_test + 1]

        A_train <- c(A_train, if_else(Reduc_t >= Reduc_c, 1, 0))

        X_train <- X_cand[1:length(A_train),]
        Y_train <- Y_cand[cbind(1:length(A_train), A_train + 1)]

        idx_t <- A_train == 1
        idx_c <- A_train == 0
        if (length(A_train) %% n_batch) {
          fit_t <- km(~ 1,
                      design = X_train[idx_t,],
                      response = Y_train[idx_t],
                      covtype = "matern5_2",
                      coef.trend = fit_t@trend.coef,
                      coef.cov = fit_t@covariance@range.val,
                      coef.var = fit_t@covariance@sd2,
                      nugget = fit_t@covariance@nugget.estim,
                      nugget.estim = F)

          fit_c <- km(~ 1,
                      design = X_train[idx_c,],
                      response = Y_train[idx_c],
                      covtype = "matern5_2",
                      coef.trend = fit_c@trend.coef,
                      coef.cov = fit_c@covariance@range.val,
                      coef.var = fit_c@covariance@sd2,
                      nugget = fit_c@covariance@nugget.estim,
                      nugget.estim = F)
        } else {
          fit_t <- km(~ 1,
                      design = X_train[idx_t,],
                      response = Y_train[idx_t],
                      covtype = "matern5_2",
                      control = list(pop.size = 50, trace = F),
                      nugget.estim = T)
          fit_c <- km(~ 1,
                      design = X_train[idx_c,],
                      response = Y_train[idx_c],
                      covtype = "matern5_2",
                      control = list(pop.size = 50, trace = F),
                      nugget.estim = T)
        }
      }
      pred_t <- predict(fit_t, X_test, type = "SK", se.compute = T)
      pred_c <- predict(fit_c, X_test, type = "SK", se.compute = T)
      target_est[b] <- mean(pred_t$mean - pred_c$mean)
      setTxtProgressBar(pb, b)
      b <- b + 1
    }, silent = T)
  }
  close(pb)
  return(target_est)
}

Simulation2_VR <- function(B, n_train, n_cand, n_init,
                           Covariate_model, Propensity_model = NULL, Outcome_model,
                           X_test, target = "ATE") {
  n_test <- nrow(X_test)
  target_est <- numeric(B)

  pb <- txtProgressBar(min = 0, max = B, style = 3)
  for (b in 1:B) {
    X_cand <- Covariate_model(n_cand)
    if (!is.null(Propensity_model)) {
      Propensity_cand <- inv_logit(apply(X_cand, 1, Propensity_model))
      A_cand <- as.numeric(runif(length(Propensity_cand)) < Propensity_cand)
    }
    Y_cand <- cbind(mapply(Outcome_model, asplit(X_cand, 1), 0),
                    mapply(Outcome_model, asplit(X_cand, 1), 1))

    idx <- sample(n_cand, n_init)
    if (is.null(Propensity_model)) {
      A_train <- as.numeric(1:n_init %in% sample(n_init, n_init %/% 2))
    } else {
      while (sum(A_cand[idx]) <= 2) {
        idx <- sample(n_cand, n_init)
      }
    }

    while (length(idx) <= n_train) {
      X_train <- X_cand[idx,]
      if (!is.null(Propensity_model)) {
        A_train <- A_cand[idx]
      }
      Y_train <- Y_cand[cbind(idx, A_train + 1)]

      idx_t <- A_train == 1
      idx_c <- A_train == 0
      fit_t <- km(~ 1,
                  design = X_train[idx_t,],
                  response = Y_train[idx_t],
                  covtype = "matern5_2",
                  control = list(pop.size = 50, trace = F),
                  nugget.estim = T)
      fit_c <- km(~ 1,
                  design = X_train[idx_c,],
                  response = Y_train[idx_c],
                  covtype = "matern5_2",
                  control = list(pop.size = 50, trace = F),
                  nugget.estim = T)

      pred_t <- predict(fit_t, rbind(X_test, X_cand), type = "SK", cov.compute = T)
      pred_c <- predict(fit_c, rbind(X_test, X_cand), type = "SK", cov.compute = T)

      # fit_propensity <- glm(A_train ~ .^2, family = "binomial", data = as_tibble(cbind(X_train, A_train)))
      # pred_propensity <- predict(fit_propensity, as_tibble(rbind(X_test, X_cand)), type = "response")
      pred_propensity <- inv_logit(apply(rbind(X_test, X_cand), 1, Propensity_model))

      if (target == "ATE") {
        w <- rep(1, n_test)
      } else if (target == "ATTE") {
        w <- pred_propensity[1:n_test]
      } else if (target == "ATO") {
        w <- pred_propensity[1:n_test]*(1 - pred_propensity[1:n_test])
      }
      Reduc_t <- colSums(w*pred_t$cov[1:n_test, (n_test + 1):(n_test + n_cand)])^2/pred_t$cov[cbind((n_test + 1):(n_test + n_cand), (n_test + 1):(n_test + n_cand))]
      Reduc_c <- colSums(w*pred_c$cov[1:n_test, (n_test + 1):(n_test + n_cand)])^2/pred_c$cov[cbind((n_test + 1):(n_test + n_cand), (n_test + 1):(n_test + n_cand))]

      if (is.null(Propensity_model)) {
        Reduc <- cbind(Reduc_c, Reduc_t)
        Reduc[idx, ] <- 0

        tmp <- which(Reduc == max(Reduc), arr.ind = T)[1, ]
        idx <- c(idx, tmp[1])
        A_train <- c(A_train, tmp[2] - 1)
      } else {
        Reduc <- pred_propensity[(n_test + 1):(n_test + n_cand)]*Reduc_t + (1 - pred_propensity[(n_test + 1):(n_test + n_cand)])*Reduc_c
        Reduc[idx] <- 0

        idx <- c(idx, which.max(Reduc))
      }
    }

    if (target == "ATE") {
      target_est[b] <- mean(pred_t$mean[1:n_test] - pred_c$mean[1:n_test])
    } else if (target == "ATTE") {
      target_est[b] <- sum(pred_propensity[1:n_test]*(pred_t$mean[1:n_test] - pred_c$mean[1:n_test])/sum(pred_propensity[1:n_test]))
    } else if (target == "ATO") {
      target_est[b] <- sum(pred_propensity[1:n_test]*(1 - pred_propensity[1:n_test])*(pred_t$mean[1:n_test] - pred_c$mean[1:n_test]))/sum(pred_propensity[1:n_test]*(1 - pred_propensity[1:n_test]))
    }

    setTxtProgressBar(pb, b)
  }
  close(pb)

  return(target_est)
}

Simulation1_ALC <- function(B, n_train, n_init, n_batch = 10,
                            Covariate_model, Outcome_model,
                            X_test) {
  n_test <- nrow(X_test)
  target_est <- numeric(B)
  
  b <- 1
  pb <- txtProgressBar(min = 0, max = B, style = 3)
  while (b <= B) {
    try({
      X_cand <- Covariate_model(n_train)
      Y_cand <- cbind(mapply(Outcome_model, asplit(X_cand, 1), 0),
                      mapply(Outcome_model, asplit(X_cand, 1), 1))
      
      X_train <- X_cand[1:n_init, ]
      A_train <- as.numeric(1:n_init %in% sample(n_init, n_init %/% 2))
      Y_train <- Y_cand[cbind(1:n_init, A_train + 1)]
      
      idx_t <- A_train == 1
      idx_c <- A_train == 0
      fit_t <- km(~ 1,
                  design = X_train[idx_t,],
                  response = Y_train[idx_t],
                  covtype = "matern5_2",
                  control = list(pop.size = 50, trace = F),
                  nugget.estim = T)
      fit_c <- km(~ 1,
                  design = X_train[idx_c,],
                  response = Y_train[idx_c],
                  covtype = "matern5_2",
                  control = list(pop.size = 50, trace = F),
                  nugget.estim = T)
      
      while (length(A_train) < n_train) {
        pred_t <- predict(fit_t, rbind(X_test, X_cand[length(A_train) + 1,]), type = "SK", se.compute = T)
        pred_c <- predict(fit_c, rbind(X_test, X_cand[length(A_train) + 1,]), type = "SK", se.compute = T)
        
        Reduc_t <- pred_t$sd[n_test + 1]
        Reduc_c <- pred_c$sd[n_test + 1]
        
        A_train <- c(A_train, if_else(Reduc_t >= Reduc_c, 1, 0))
        
        X_train <- X_cand[1:length(A_train),]
        Y_train <- Y_cand[cbind(1:length(A_train), A_train + 1)]
        
        idx_t <- A_train == 1
        idx_c <- A_train == 0
        if (length(A_train) %% n_batch) {
          fit_t <- km(~ 1,
                      design = X_train[idx_t,],
                      response = Y_train[idx_t],
                      covtype = "matern5_2",
                      coef.trend = fit_t@trend.coef,
                      coef.cov = fit_t@covariance@range.val,
                      coef.var = fit_t@covariance@sd2,
                      nugget = fit_t@covariance@nugget.estim,
                      nugget.estim = F)
          
          fit_c <- km(~ 1,
                      design = X_train[idx_c,],
                      response = Y_train[idx_c],
                      covtype = "matern5_2",
                      coef.trend = fit_c@trend.coef,
                      coef.cov = fit_c@covariance@range.val,
                      coef.var = fit_c@covariance@sd2,
                      nugget = fit_c@covariance@nugget.estim,
                      nugget.estim = F)
        } else {
          fit_t <- km(~ 1,
                      design = X_train[idx_t,],
                      response = Y_train[idx_t],
                      covtype = "matern5_2",
                      control = list(pop.size = 50, trace = F),
                      nugget.estim = T)
          fit_c <- km(~ 1,
                      design = X_train[idx_c,],
                      response = Y_train[idx_c],
                      covtype = "matern5_2",
                      control = list(pop.size = 50, trace = F),
                      nugget.estim = T)
        }
      }
      pred_t <- predict(fit_t, X_test, type = "SK", se.compute = T)
      pred_c <- predict(fit_c, X_test, type = "SK", se.compute = T)
      target_est[b] <- mean(pred_t$mean - pred_c$mean)
      setTxtProgressBar(pb, b)
      b <- b + 1
    }, silent = T)
  }
  close(pb)
  return(target_est)
}

Simulation2_ALC <- function(B, n_train, n_cand, n_init,
                            Covariate_model, Propensity_model = NULL, Outcome_model,
                            X_test, target = "ATE") {
  n_test <- nrow(X_test)
  target_est <- numeric(B)

  pb <- txtProgressBar(min = 0, max = B, style = 3)
  for (b in 1:B) {
    X_cand <- Covariate_model(n_cand)
    if (!is.null(Propensity_model)) {
      Propensity_cand <- inv_logit(apply(X_cand, 1, Propensity_model))
      A_cand <- as.numeric(runif(length(Propensity_cand)) < Propensity_cand)
    }
    Y_cand <- cbind(mapply(Outcome_model, asplit(X_cand, 1), 0),
                    mapply(Outcome_model, asplit(X_cand, 1), 1))

    idx <- sample(n_cand, n_init)
    if (is.null(Propensity_model)) {
      A_train <- as.numeric(1:n_init %in% sample(n_init, n_init %/% 2))
    } else {
      while (sum(A_cand[idx]) <= 2) {
        idx <- sample(n_cand, n_init)
      }
    }
    while (length(idx) <= n_train) {
      X_train <- X_cand[idx,]
      if (!is.null(Propensity_model)) {
        A_train <- A_cand[idx]
      }
      Y_train <- Y_cand[cbind(idx, A_train + 1)]

      idx_t <- A_train == 1
      idx_c <- A_train == 0
      fit_t <- km(~ 1,
                  design = X_train[idx_t,],
                  response = Y_train[idx_t],
                  covtype = "matern5_2",
                  control = list(pop.size = 50, trace = F),
                  nugget.estim = T)
      fit_c <- km(~ 1,
                  design = X_train[idx_c,],
                  response = Y_train[idx_c],
                  covtype = "matern5_2",
                  control = list(pop.size = 50, trace = F),
                  nugget.estim = T)

      pred_t <- predict(fit_t, rbind(X_test, X_cand), type = "SK", se.compute = T)
      pred_c <- predict(fit_c, rbind(X_test, X_cand), type = "SK", se.compute = T)

      Reduc_t <- pred_t$sd[(n_test + 1):(n_test + n_cand)]^2
      Reduc_c <- pred_c$sd[(n_test + 1):(n_test + n_cand)]^2

      if (is.null(Propensity_model)) {
        Reduc <- cbind(Reduc_c, Reduc_t)
        Reduc[idx, ] <- 0

        tmp <- which(Reduc == max(Reduc), arr.ind = T)[1, ]
        idx <- c(idx, tmp[1])
        A_train <- c(A_train, tmp[2] - 1)
      } else {
        # fit_propensity <- glm(A_train ~ .^2, family = "binomial", data = as_tibble(cbind(X_train, A_train)))
        # pred_propensity <- predict(fit_propensity, as_tibble(rbind(X_test, X_cand)), type = "response")
        pred_propensity <- inv_logit(apply(rbind(X_test, X_cand), 1, Propensity_model))

        Reduc <- pred_propensity[(n_test + 1):(n_test + n_cand)]*Reduc_t + (1 - pred_propensity[(n_test + 1):(n_test + n_cand)])*Reduc_c
        Reduc[idx] <- 0

        idx <- c(idx, which.max(Reduc))
      }
    }

    if (target == "ATE") {
      target_est[b] <- mean(pred_t$mean[1:n_test] - pred_c$mean[1:n_test])
    } else if (target == "ATTE") {
      target_est[b] <- sum(pred_propensity[1:n_test]*(pred_t$mean[1:n_test] - pred_c$mean[1:n_test]))/sum(pred_propensity[1:n_test])
    } else if (target == "ATO") {
      target_est[b] <- sum(pred_propensity[1:n_test]*(1 - pred_propensity[1:n_test])*(pred_t$mean[1:n_test] - pred_c$mean[1:n_test]))/sum(pred_propensity[1:n_test]*(1 - pred_propensity[1:n_test]))
    }
    setTxtProgressBar(pb, b)
  }
  close(pb)

  return(target_est)
}

Simulation3_random <- function(B, n_train,
                               Covariate_model, Propensity_model, Outcome_model) {
  Effect <- numeric(B)
  
  pb <- txtProgressBar(min = 0, max = B, style = 3)
  for (b in 1:B) {
    X_train <- Covariate_model(n_train)
    Propensity_train <- inv_logit(apply(X_train, 1, Propensity_model))
    A_train <- as.numeric(runif(length(Propensity_train)) < Propensity_train)
    Y_train <- mapply(Outcome_model, asplit(X_train, 1), A_train)
    
    Effect[b] <- sum(A_train*(mapply(Outcome_model, asplit(X_train, 1), 1) - mapply(Outcome_model, asplit(X_train, 1), 0)))
    setTxtProgressBar(pb, b)
  }
  close(pb)
  
  return(Effect)
}

Simulation3_greedy <- function(B, n_train, n_cand, n_init,
                               Covariate_model, Propensity_model, Outcome_model) {
  Effect <- numeric(B)
  
  pb <- txtProgressBar(min = 0, max = B, style = 3)
  for (b in 1:B) {
    X_cand <- Covariate_model(n_cand)
    Propensity_cand <- inv_logit(apply(X_cand, 1, Propensity_model))
    A_cand <- as.numeric(runif(length(Propensity_cand)) < Propensity_cand)
    Y_cand <- cbind(mapply(Outcome_model, asplit(X_cand, 1), 0),
                    mapply(Outcome_model, asplit(X_cand, 1), 1))
    
    idx <- sample(n_cand, n_init)
    while (sum(A_cand[idx]) <= 2) {
      idx <- sample(n_cand, n_init)
    }
    
    while (length(idx) <= n_train) {
      X_train <- X_cand[idx,]
      A_train <- A_cand[idx]
      Y_train <- Y_cand[cbind(idx, A_train + 1)]
      
      idx_t <- A_train == 1
      idx_c <- A_train == 0
      fit_t <- km(~ 1,
                  design = X_train[idx_t,],
                  response = Y_train[idx_t],
                  covtype = "matern5_2",
                  control = list(pop.size = 50, trace = F),
                  nugget.estim = T)
      fit_c <- km(~ 1,
                  design = X_train[idx_c,],
                  response = Y_train[idx_c],
                  covtype = "matern5_2",
                  control = list(pop.size = 50, trace = F),
                  nugget.estim = T)
      
      pred_t <- predict(fit_t, X_cand, type = "SK", se.compute = T)
      pred_c <- predict(fit_c, X_cand, type = "SK", se.compute = T)
      
      # fit_propensity <- glm(A_train ~ .^2, family = "binomial", data = as_tibble(cbind(X_train, A_train)))
      # pred_propensity <- predict(fit_propensity, as_tibble(rbind(X_test, X_cand)), type = "response")
      pred_propensity <- inv_logit(apply(X_cand, 1, Propensity_model))
      
      Expected_effect <- pred_propensity*(pred_t$mean - pred_c$mean)
      Expected_effect[idx] <- -Inf
      idx <- c(idx, which.max(Expected_effect))
    }
    
    Effect[b] <- sum(A_train*(mapply(Outcome_model, asplit(X_train, 1), 1) - mapply(Outcome_model, asplit(X_train, 1), 0)))
    
    setTxtProgressBar(pb, b)
  }
  close(pb)
  
  return(Effect)
}


Simulation3_UCB <- function(B, n_train, n_cand, n_init,
                            Covariate_model, Propensity_model, Outcome_model, c = 1) {
  Effect <- numeric(B)
  
  pb <- txtProgressBar(min = 0, max = B, style = 3)
  for (b in 1:B) {
    X_cand <- Covariate_model(n_cand)
    Propensity_cand <- inv_logit(apply(X_cand, 1, Propensity_model))
    A_cand <- as.numeric(runif(length(Propensity_cand)) < Propensity_cand)
    Y_cand <- cbind(mapply(Outcome_model, asplit(X_cand, 1), 0),
                    mapply(Outcome_model, asplit(X_cand, 1), 1))
    
    idx <- sample(n_cand, n_init)
    while (sum(A_cand[idx]) <= 2) {
      idx <- sample(n_cand, n_init)
    }
    
    while (length(idx) <= n_train) {
      X_train <- X_cand[idx,]
      A_train <- A_cand[idx]
      Y_train <- Y_cand[cbind(idx, A_train + 1)]
      
      idx_t <- A_train == 1
      idx_c <- A_train == 0
      fit_t <- km(~ 1,
                  design = X_train[idx_t,],
                  response = Y_train[idx_t],
                  covtype = "matern5_2",
                  control = list(pop.size = 50, trace = F),
                  nugget.estim = T)
      fit_c <- km(~ 1,
                  design = X_train[idx_c,],
                  response = Y_train[idx_c],
                  covtype = "matern5_2",
                  control = list(pop.size = 50, trace = F),
                  nugget.estim = T)
      
      pred_t <- predict(fit_t, X_cand, type = "SK", se.compute = T)
      pred_c <- predict(fit_c, X_cand, type = "SK", se.compute = T)
      
      # fit_propensity <- glm(A_train ~ .^2, family = "binomial", data = as_tibble(cbind(X_train, A_train)))
      # pred_propensity <- predict(fit_propensity, as_tibble(rbind(X_test, X_cand)), type = "response")
      pred_propensity <- inv_logit(apply(X_cand, 1, Propensity_model))
      
      Expected_effect <- pred_propensity*(pred_t$mean - pred_c$mean)
      beta <- c*log(length(idx))
      UCB <- Expected_effect + sqrt(beta*(pred_propensity*(pred_t$sd^2 + pred_c$sd^2) +
                                            pred_propensity*(1 - pred_propensity)*(pred_t$mean - pred_c$mean)^2))
      UCB[idx] <- -Inf
      idx <- c(idx, which.max(UCB))
    }
    
    Effect[b] <- sum(A_train*(mapply(Outcome_model, asplit(X_train, 1), 1) - mapply(Outcome_model, asplit(X_train, 1), 0)))
    
    setTxtProgressBar(pb, b)
  }
  close(pb)
  
  return(Effect)
}

Eval <- function(est, true) {
  list(bias = mean(est) - true, RMSE = sqrt(mean((est - true)^2)))
}
