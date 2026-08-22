# IWLS penalizzato (Kosmidis-Firth/Jeffreys) con criterio di ARMIJO halving
# La funzione individua il valore del vettore dei parametri beta che massimizzano la massima verosimiglianza penalizzata per dati binari con diversi legami. 
# L'Articolo su cui si basa la struttura della funzione implementata è Kosmidis-Firth 2020. Per maggiori informazione rispetto ai fodnamenti teorici si consiglia di leggere il lavoro originale, mentre oppure la tesi dell'autore.  

IWLS_AJKF <- function(formula, data = NULL, # funziona anche senza dataset
                      a = 1 / 2, # valore della penalità, 1/2 per eliminare la distorsione del proimo ordine con legame canonico
                      link = "logit",
                      beta_start = NULL,
                      tol = 1e-12,
                      maxit = 1000,
                      max_half = 30, # backtracking max steps
                      c1 = 1e-4, # costante di Armijo (aumento sufficiente)
                      rho = 0.5,
                      na.action = na.omit) {
  # Salva la chiamata per riproducibilità
  call <- match.call()

  # Estrai il dataframe in input e i termini
  mf <- model.frame(formula, data = data, na.action = na.action)
  terms_obj <- attr(mf, "terms")

  # Matrice del disegno (intercetta + covariate) e risposte
  X <- model.matrix(formula, data = data)
  y <- mf[, 1][, 1] # estrae i successi dalla formula
  m <- mf[, 1][, 1] + mf[, 1][, 2]

  valid_links <- c("logit", "probit", "cloglog", "loglog", "cauchit")
  if (!link %in% valid_links) {
    stop("Gli unici link ammessi sono: ", paste(valid_links, collapse = ", "), ".")
  }

  n <- nrow(X)
  p <- ncol(X)
  df <- n - p
  if (is.null(n)) stop("'X' deve essere una matrice.")
  if (!is.null(beta_start) && length(beta_start) != p) {
    stop("`beta_start` deve essere NULL oppure un vettore di lunghezza p.")
  }

  source("AlgoritmiImplementati/.links_GLM.R")
  lnk <- .links_GLM()[[link]] # estrae funzione legame e derivate del legame scelto

  # Soglie di clipping per stabilità numerica
  eps_pi <- 1e-10
  w_max <- 1e8
  half <- 0
  tot_armijo <- numeric(maxit)

  conv <- matrix(NA, nrow = maxit, ncol = p)
  colnames(conv) <- colnames(X)
  converged <- TRUE

  #### Log-verosimiglianza penalizzata: l*(b) = l(b) + a*log|I(b)|
  # Serve per calcolare AIC, BIC e criterio di Armijo
  Penalized_likelihood <- function(beta, a) {
    eta <- as.vector(X %*% beta)
    Pi <- lnk$G(eta)
    Pi <- pmin(pmax(Pi, eps_pi), 1 - eps_pi)
    d <- lnk$g(eta)
    ll <- sum(y * log(Pi) + (m - y) * log(1 - Pi))
    w_m <- m * d^2 / (Pi * (1 - Pi))
    R <- qr.R(qr(sqrt(w_m) * X))
    log_det_I <- 2 * sum(log(abs(diag(R))))
    as.numeric(ll + a * log_det_I)
  }

  Deviance <- function(obs, fitt, m) {
    # Protezione 0*log(0) = 0
    t1 <- ifelse(obs > 0, obs * log(obs / fitt), 0)
    t2 <- ifelse(obs < 1, (1 - obs) * log((1 - obs) / (1 - fitt)), 0)
    D <- 2 * m * (t1 + t2)
    residuals <- sign(obs - fitt) * sqrt(pmax(D, 0)) # pmax per evitare problemi numerici
    return(list(
      dev.residuals = residuals,
      deviance = sum(D)
    ))
  }

  phi_hat <- function(obs, fitt, m, df) {
    pearson_sq <- m * (obs - fitt)^2 / (fitt * (1 - fitt))
    X2 <- sum(pearson_sq)
    phi_hat <- X2 / df
    Presiduals <- sign(obs - fitt) * sqrt(pmax(pearson_sq, 0))
    return(list(phi_hat = phi_hat, pearson.res = Presiduals, X2 = X2))
  }

  ##### Estrazione in modo compatto le quantità finali per evitare duplicazioni
  finalize <- function(beta_new, iter, delta_trace, converged) {
    eta_f <- as.vector(X %*% beta_new)
    Pi_f <- lnk$G(eta_f)
    Pi_f <- pmin(pmax(Pi_f, eps_pi), 1 - eps_pi)
    d_f <- lnk$g(eta_f)
    dd_f <- lnk$gdash(eta_f)
    w_f <- (d_f^2) / (Pi_f * (1 - Pi_f))
    w_m_f <- m * w_f
    WX_f <- qr(sqrt(w_m_f) * X)
    h_f <- rowSums(qr.Q(WX_f)^2)
    R_f <- qr.R(WX_f)
    q_f <- dd_f / w_f + Pi_f
    Cov.mat <- chol2inv(R_f)

    dev <- Deviance(y / m, Pi_f, m)
    x2 <- phi_hat(y / m, Pi_f, m, df)

    # Calcolo del modello nullo per devianza e test
    pi0 <- sum(y) / sum(m)
    beta_null <- lnk$invG(pi0)
    Pi0_null <- rep(lnk$G(beta_null), n)
    Pi0_null <- pmin(pmax(Pi0_null, eps_pi), 1 - eps_pi)
    null_deviance <- Deviance(y / m, Pi0_null, m)
    null_dev <- null_deviance$deviance
    df_null <- n - 1
    pen_lik <- Penalized_likelihood(beta_new, a)

    # Calcolo dei test
    sd <- sqrt(diag(Cov.mat))
    z_stat <- beta_new / sd
    p_values <- format((2 * (1 - pnorm(abs(z_stat)))), scientific = TRUE, digits = 3)

    coef.test <- data.frame(
      coef = beta_new,
      s.d = sd,
      z.stat = z_stat,
      p.values = p_values
    )
    null.X2 <- null_deviance$deviance - dev$deviance

    list(
      call = call,
      formula = formula,
      terms = terms_obj,
      Convergence = data.frame(
        Iteration = seq_len(iter),
        Delta_beta = delta_trace[seq_len(iter)]
      ),
      converged = converged,
      Armijo_halving = tot_armijo[1:iter],
      beta_log = conv[1:iter, , drop = FALSE],
      Penalized_likelihood = pen_lik,
      Penalized_Score = as.vector(
        crossprod(X, (d_f / (Pi_f * (1 - Pi_f))) *
          (y - m * Pi_f + 2 * a * h_f * (q_f - 0.5)))
      ),
      probs = Pi_f,
      lin.pred = eta_f,
      weights = w_f,
      Cov_mat = Cov.mat,
      df_null = df_null,
      null_deviance = null_dev,
      df_model = df,
      deviance = dev$deviance,
      residuals.dev = dev$dev.residuals,
      residuals.X2 = x2$pearson.res,
      phi_hat = x2$phi_hat,
      leverage = h_f,
      AIC = -2 * pen_lik + 2 * p,
      BIC = -2 * pen_lik + p * log(n),
      test = coef.test,
      Null.test = format(pchisq(null.X2, p - 1, lower.tail = FALSE), scientific = TRUE, digits = 3),
      beta = beta_new
    )
  }

  # Inizializzazione
  if (is.null(beta_start)) {
    # costante aumentata per dataset piccoli (Cordeiro-McCullagh) proporzionale alla distorsione
    c_start <- 0.5 * (p / n)
    Pi0 <- (y + c_start) / (m + 2 * c_start)
    eta0 <- lnk$invG(Pi0)
    d0 <- lnk$g(eta0)
    w_m0 <- m * (d0^2) / (Pi0 * (1 - Pi0))
    qr0 <- qr(sqrt(w_m0) * X)
    beta_old <- solve(qr0, sqrt(w_m0) * eta0)
    conv[1, ] <- beta_old
  } else {
    beta_old <- beta_start
    conv[1, ] <- beta_old
  }

  delta_trace <- numeric(maxit)
  delta_trace[1] <- NA

  # Ciclo IWLS
  for (iter in 2:maxit) {
    eta <- as.vector(X %*% beta_old)
    Pi <- lnk$G(eta)
    Pi <- pmin(pmax(Pi, eps_pi), 1 - eps_pi) # clip prima del denominatore
    d <- lnk$g(eta)
    d <- ifelse(abs(d) < 1e-12, sign(d) * 1e-12, d) # clip a denominatore tra valori estremi
    dd <- lnk$gdash(eta)

    w <- (d^2) / (Pi * (1 - Pi))
    w <- pmin(pmax(w, 1e-14), w_max) # tetto finito per evitare pesi esplosivi
    w_m <- m * w
    WX <- qr(sqrt(w_m) * X)
    h <- rowSums(qr.Q(WX)^2)
    q <- dd / w + Pi

    # Update IWLS
    z_star <- eta + (y / m - Pi + 2 * a * h * (q - 0.5) / m) / d
    beta_full_step <- solve(WX, sqrt(w_m) * z_star) # classico passo IWLS
    step_direction <- beta_full_step - beta_old # variazione (Δbeta = I^{-1} U)

    # Bypass convergenza se il gradiente è già infinitesimo
    delta_pre <- sqrt(sum(step_direction^2))
    if (delta_pre < tol) {
      conv[iter, ] <- beta_full_step
      delta_trace[iter] <- delta_pre
      tot_armijo[iter] <- 0
      return(finalize(beta_full_step, iter, delta_trace, converged = TRUE))
    }

    # LINE SEARCH: BACKTRACKING (N&W - ottimizzazione numerica 3.1)
    score_old <- as.vector(
      crossprod(X, (d / (Pi * (1 - Pi))) * (y - m * Pi + 2 * a * h * (q - 0.5)))
    )
    dir_deriv <- sum(score_old * step_direction) # derivata direzionale
    f_old <- Penalized_likelihood(beta_old, a) # controllo sulla funzione da massimizzare

    alpha <- 1.0 # partenza
    armijo_ok <- FALSE
    half <- 0 # reset

    while (half < max_half) {
      half <- half + 1 # tiene conto di quanto si riduce il passo
      beta_candidate <- beta_old + alpha * step_direction
      f_cand <- Penalized_likelihood(beta_candidate, a) # -Inf se inammissibile

      # Criterio Armijo per massimizzazione
      # Condizione curvatura di Wolfe non necessaria per proprietà del Backtracking (N&W, pag. 42)
      if (is.finite(f_cand) && f_cand >= (f_old + c1 * alpha * dir_deriv - 1e-10)) {
        beta_new <- beta_candidate
        armijo_ok <- TRUE
        break
      }
      alpha <- alpha * rho # dimezzamento in caso non converga
    }

    tot_armijo[iter] <- half - 1 # storico armijo

    if (armijo_ok == FALSE) {
      # Nessun passo soddisfa Armijo, non accettare punti peggiori
      conv[iter, ] <- beta_old
      delta_trace[iter] <- 0
      message(sprintf("Iterazione %d: nessun valore di alpha individuato, provare con Rho diverso.", iter))
      return(finalize(beta_old, iter, delta_trace, converged = FALSE))
    }

    conv[iter, ] <- beta_new
    delta <- sqrt(sum((beta_new - beta_old)^2)) # Criterio differenze iterate (Russo cap. 9)
    delta_trace[iter] <- delta

    if (delta < tol) {
      return(finalize(beta_new, iter, delta_trace, converged = TRUE))
    }

    beta_old <- beta_new
  }

  warning("Nessuna convergenza in ", maxit, " iterazioni.")
  finalize(beta_new, maxit, delta_trace, converged = FALSE)
}
