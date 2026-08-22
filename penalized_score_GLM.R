# Funzione che valuta lo score penalizzato in un dato valore delvettore dei parametri beta.
Penalized_score_GLM <- function(X, y, m, beta,
                                a           = 1 / 2,
                                link        = "Logit") {
  
  valid_links <- c("Logit", "Probit", "cLogLog", "LogLog", "Cauchit")
  if (!link %in% valid_links)
    stop("The only Links allowed are: Logit, Probit, cLogLog, LogLog and Cauchit.")
  
  if (length(beta) != ncol(X))
    stop("`beta` deve avere lunghezza pari al numero di colonne di X.")
  
  # importiamo i legami con le relative funzioni di interesse
  source("AlgoritmiImplementati/.links_GLM.R")
  
  lnk <- .links_GLM()[[link]]
  
  # Predittore lineare e quantità del legame
  eta <- as.vector(X %*% beta) # grazie a beta individuiamo il predittore lineare
  Pi  <- lnk$G(eta)
  d   <- lnk$g(eta)        # G'(η)
  dd  <- lnk$gdash(eta)    # G''(η)
  
  # Pesi 
  w   <- (d^2) / (Pi * (1 - Pi))
  w_m <- m * w
  
  # Decomposizione QR
  V  <- sqrt(w_m) * X
  WX <- qr(V)
  h  <- rowSums(qr.Q(WX)^2)    # h_i = diag(H)
  
  q <- dd / w + Pi
  
  # Score penalizzato:
  score <- as.vector(
    crossprod(X, (d / (Pi * (1 - Pi))) * (y - m * Pi + 2 * a * h * (q - 0.5)))
  )
  
  nomi <- if (!is.null(colnames(X))) colnames(X) else paste0("beta_", seq_along(score))
  names(score) <- nomi
  
  norm_inf <- max(abs(score))
  norm_2   <- sqrt(sum(score^2))
  
    return(list(
      score    = score,
      norm_inf = norm_inf,
      norm_2   = norm_2
    ))
}