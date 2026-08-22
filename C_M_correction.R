# Algoritmo per l'individuazione della distorsione del primo ordine per le stime di massima verosimiglianza a posteriori. La struttura ella funzione segue l'Articolo di Cordeiro McCullagh 1991 ed è stata implememntata per i dati binari con diversi legami.

C_M_correction <- function(X, m, beta_ml, link = 'logit'){
  
  source("AlgoritmiImplementati/.links_GLM.R")
  lnk <- .links_GLM()[[link]]
  
  # Predittore lineare e quantità del legame
  eta <- as.vector(X %*% beta_ml)
  Pi  <- lnk$G(eta)
  d   <- lnk$g(eta)        # G'(eta)
  dd  <- lnk$gdash(eta)    # G''(eta)
  
  # Pesi w_i = m * G'(eta)^2 / V(mu)
  V_mu <- Pi * (1 - Pi)
  w    <- (d^2) / V_mu
  w_m  <- m * w
  
  # Decomposizione QR di W_m^{1/2} X  -> leverages h_i
  V   <- sqrt(w_m) * X
  QR  <- qr(V)
  h   <- rowSums(qr.Q(QR)^2)
  
  # Vettore zeta_i = (1/2) * h_i * G''(eta_i) / G'(eta_i)
  zeta <- 0.5 * h * dd / d
  
  # Il bias del primo ordine è b = -(X' W X)^{-1} X' zeta.
  # Risolviamo il sistema lineare: (X' W X) B = - X' zeta
  Info_matrix <- crossprod(V)         # Equivale a X' W X
  Term_noto   <- crossprod(X, -zeta)  # Equivale a - X' zeta
  
  # Calcolo del bias
  B <- solve(Info_matrix, Term_noto)
  
  return(B)
}