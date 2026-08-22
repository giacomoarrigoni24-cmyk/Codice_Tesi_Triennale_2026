#   GRAFICI PER LA TESI
#   Codici utilizzati per le figure, organizzati per capitoli.
#

#   - I file sorgente ("AlgoritmiImplementati/IWLS_AJKF.R", ".links_GLM.R",
#     "C_M_correction.R") devono essere presenti nella directory di lavoro.
#   - I dataset "Heart", "Beetles", "endometrial" provengono dal pacchetto MLGdata.
#   - Il pacchetto brglm2 è richiesto solo per i confronti con la funzione brglmFit.

# Librerie
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(palmerpenguins)
library(MLGdata)

# Funzioni
source("AlgoritmiImplementati/.links_GLM.R")
source("AlgoritmiImplementati/IWLS_AJKF.R")
source("AlgoritmiImplementati/C_M_correction.R")

# Tema per i grafici
tema_base <- theme_gray(base_size = 14)

###############################################################################
#                          CAPITOLO 1
#
# FIGURA 1.1 — Log-verosimiglianza per un modello lineare senza intercetta.
# Dati penguins, covariata lunghezza del becco (standardizzata)
###############################################################################

peng <- na.omit(penguins)
peng$stdmas <- (peng$body_mass_g - mean(peng$body_mass_g)) / sd(peng$body_mass_g)

mod_lineare <- glm(stdmas ~ bill_length_mm - 1, data = peng, family = gaussian())
disp_mod <- summary(mod_lineare)$dispersion

y_lin <- peng$stdmas
x_lin <- peng$bill_length_mm
beta_hat_lm <- coef(mod_lineare)
n_lin <- nrow(peng)

beta_seq_lm <- seq(beta_hat_lm - 0.002, beta_hat_lm + 0.002, length.out = 200)

loglik_valori <- sapply(beta_seq_lm, function(b) {
  -(n_lin / 2) * log(2 * pi * disp_mod) -
    (1 / (2 * disp_mod)) * sum((y_lin - x_lin * b)^2)
})

loglik_hat <- -(n_lin / 2) * log(2 * pi * disp_mod) -
  (1 / (2 * disp_mod)) * sum((y_lin - x_lin * beta_hat_lm)^2)

df_loglik <- data.frame(beta = beta_seq_lm, loglik = loglik_valori)

ggplot(df_loglik, aes(x = beta, y = loglik)) +
  geom_vline(xintercept = 0, col = "black", linetype = "dashed") +
  geom_line(linewidth = 1.2, color = "steelblue") +
  geom_hline(yintercept = loglik_hat, col = "red", linetype = "dotdash") +
  geom_point(x = beta_hat_lm, y = loglik_hat, color = "red", size = 3) +
  labs(
    x = expression(beta),
    y = expression(U(beta)),
    subtitle = paste0(
      "Stima OLS \u03b2: ", round(beta_hat_lm, 6),
      "        Log-verosimiglianza: ", round(loglik_hat, 2)
    )
  ) +
  tema_base


###############################################################################
# FIGURA 1.2 — Consistenza asintotica dello stimatore di massima verosimiglianza
# Simulazione da un modello binomiale con legame canonico, beta = (-1, 0.8).
###############################################################################

set.seed(123)
true_beta_cons <- c(-1.0, 0.8)
legame_cons <- "logit"
lnk_cons <- .links_GLM()[[legame_cons]]
m_cons <- 1

esperimenti_n <- c(
  15, 25, seq(30, 200, 10), seq(250, 1000, 50),
  seq(1000, 10000, 100), seq(10000, 100000, 1000)
)

risultati_cons <- matrix(NA, ncol = 2, nrow = length(esperimenti_n))
colnames(risultati_cons) <- c("intercetta", "var_1")

for (i in seq_along(esperimenti_n)) {
  n_g <- esperimenti_n[i]
  X_sim <- cbind(rep(1, n_g), rnorm(n_g))

  eta_sim <- X_sim %*% true_beta_cons
  Pi_sim <- lnk_cons$G(eta_sim)
  y_sim <- as.vector(rbinom(n_g, m_cons, Pi_sim))

  X_data <- as.data.frame(X_sim)
  X_data$y <- y_sim
  X_data$m_g <- m_cons

  modello_ML <- glm(cbind(y, m_g - y) ~ . - 1,
    data = X_data, family = binomial(link = "logit")
  )

  mu_hat <- fitted(modello_ML)
  separazione <- any(mu_hat >= 1 - 1e-12) | any(mu_hat <= 1e-12)

  if (!separazione && modello_ML$converged) {
    risultati_cons[i, 1] <- coef(modello_ML)[1] - true_beta_cons[1]
    risultati_cons[i, 2] <- coef(modello_ML)[2] - true_beta_cons[2]
  }
}


df_cons <- as.data.frame(risultati_cons)
df_cons$n <- esperimenti_n
# saveRDS(df_cons, 'consist_beta_canonico')

# Grafico della convergenza (a partire da n = 450 per leggibilità)
ggplot(df_cons[25:nrow(df_cons), ], aes(x = n, y = var_1)) +
  geom_line(col = "steelblue", linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", col = "red", linewidth = 0.8) +
  labs(
    x        = "Numero di osservazioni (n)",
    y        = expression(hat(beta)[1] - beta[1]),
    subtitle = "Convergenza empirica dello SMV"
  ) +
  scale_x_continuous(breaks = c(
    450, seq(3000, 10000, 4000),
    seq(10000, 100000, 20000)
  )) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


###############################################################################
# FIGURA 1.3 — Distribuzione asintotica dello SMV per una regressione binaria
# 1000 campioni di dimensione 100 da una Bernoulli, beta = (-1, 0.8).
###############################################################################

set.seed(123)
n_dist <- 100
N_dist <- 500
true_beta_dist <- c(-1.0, 0.8)
lnk_dist <- .links_GLM()[["logit"]]

# otteniamo varianza teorica come info di fisher
X_dist <- cbind(rep(1, n_dist), rnorm(n_dist))
eta_dist <- X_dist %*% true_beta_dist
Pi_dist <- lnk_dist$G(eta_dist)
g_dist <- lnk_dist$g(eta_dist)
w_dist <- as.vector(g_dist^2 / (Pi_dist * (1 - Pi_dist)))
R_dist <- qr.R(qr(sqrt(w_dist) * X_dist))
Cov_asint <- chol2inv(R_dist)
sd_asint <- sqrt(diag(Cov_asint))

stime_beta2 <- rep(NA, N_dist)

for (i in 1:N_dist) {
  y_d <- rbinom(n_dist, 1, Pi_dist)
  X_data_d <- data.frame(V2 = X_dist[, 2], y = y_d)

  mod_d <- glm(y ~ V2, data = X_data_d, family = binomial())
  mu_d <- fitted(mod_d)
  sep_d <- any(mu_d >= 1 - 1e-12) | any(mu_d <= 1e-12)

  if (!sep_d) stime_beta2[i] <- coef(mod_d)[2]
}

df_dist <- data.frame(stima = stime_beta2)

ggplot(df_dist, aes(x = stima)) +
  geom_histogram(aes(y = after_stat(density)),
    bins = 40,
    fill = "steelblue", alpha = 0.6
  ) +
  stat_function(
    fun = dnorm,
    args = list(mean = true_beta_dist[2], sd = sd_asint[2]),
    color = "darkred", linewidth = 1
  ) +
  geom_vline(xintercept = mean(stime_beta2, na.rm = TRUE)) +
  labs(
    subtitle = "Distribuzione empirica dello stimatore di massima verosimiglianze rispetto a quella teorica",
    x = expression(beta[2]),
    y = "Densità"
  ) +
  tema_base

# Test di normalità
shapiro.test(stime_beta2[!is.na(stime_beta2)])


###############################################################################
# FIGURA 1.4 — Ordine di convergenza dello SMV (scala log-log)
# Stima della velocità di convergenza: pendenza attesa circa -0.5.
###############################################################################

df_loglog <- df_cons[25:nrow(df_cons), ]

fit_loglog <- lm(log(abs(var_1)) ~ log(n), data = df_loglog)

ggplot(df_loglog, aes(x = log(n), y = log(abs(var_1)))) +
  geom_line(col = "steelblue", linewidth = 0.4, alpha = 0.6) +
  geom_smooth(
    method = "lm", col = "darkred", linetype = "dashed",
    se = FALSE, linewidth = 0.7
  ) +
  labs(
    x        = "log(n)",
    y        = expression("log(" ~ hat(beta)[1] - beta[1] ~ ")"),
    subtitle = "Ordine di convergenza dello SMV"
  ) +
  tema_base

fit_loglog$coefficients[2]

###############################################################################
# FIGURA 1.5 — GLM binomiale (cLogLog) vs Regressione lineare classica
# Dataset Beetles (Bliss, 1935).
###############################################################################

Beetles$logdose_m <- Beetles$logdose - mean(Beetles$logdose)
Beetles$prop <- Beetles$uccisi / Beetles$num

mod_glm_beet <- glm(cbind(uccisi, num - uccisi) ~ logdose_m - 1,
  data = Beetles, family = binomial(link = "cloglog")
)
mod_lm_beet <- lm(prop ~ logdose, data = Beetles)

griglia_beet <- seq(min(Beetles$logdose) - 0.02,
  max(Beetles$logdose) + 0.02,
  length.out = 20000
)
X_grid_beet <- data.frame(
  logdose = griglia_beet,
  logdose_m = griglia_beet - mean(Beetles$logdose)
)

pred_lm <- predict(mod_lm_beet, newdata = X_grid_beet)
pred_glm <- predict(mod_glm_beet, newdata = X_grid_beet, type = "response")

df_pred <- data.frame(
  logdose = rep(griglia_beet, 2),
  Probabilita = c(pred_lm, pred_glm),
  Modello = rep(c("Regressione lineare", "GLM cloglog"),
    each = length(griglia_beet)
  )
)

ggplot(df_pred, aes(x = logdose, y = Probabilita)) +
  geom_line(col = "steelblue", linewidth = 1) +
  geom_point(
    data = Beetles, aes(x = logdose, y = prop),
    color = "blue", size = 2
  ) +
  facet_grid(. ~ Modello) +
  tema_base


###############################################################################
# FIGURA 1.6 — Linearità dello score per un modello lineare
# Dati penguins, covariata lunghezza del becco (standardizzata).
###############################################################################

score_valori_lm <- sapply(beta_seq_lm, function(b) {
  sum(x_lin * (y_lin - x_lin * b))
})

df_score_lm <- data.frame(beta = beta_seq_lm, score = score_valori_lm)

ggplot(df_score_lm, aes(x = beta, y = score)) +
  geom_hline(yintercept = 0, col = "black") +
  geom_vline(xintercept = 0, col = "black") +
  geom_line(linewidth = 1.2, color = "steelblue") +
  geom_point(aes(x = beta_hat_lm, y = 0), color = "black", size = 3) +
  labs(
    x        = expression(beta),
    y        = expression(U(beta)),
    subtitle = paste("Stima OLS \u03b2:", round(beta_hat_lm, 6))
  ) +
  scale_x_continuous(breaks = seq(-0.001, 0.004, by = 0.0005)) +
  tema_base


###############################################################################
# FIGURA 1.7 — Non linearità dello score per un GLM binomiale (Beetles)
# Legame canonico (logit), modello senza intercetta.
###############################################################################

mod_glm_score <- glm(cbind(uccisi, num - uccisi) ~ logdose_m - 1,
  data = Beetles, family = binomial(link = "logit")
)
beta_hat_beet <- coef(mod_glm_score)

Score_bin <- function(xi, y, m_i, beta) {
  lnk <- .links_GLM()[["logit"]]
  eta <- as.numeric(xi * beta)
  pi_i <- lnk$G(eta)
  varianza <- pi_i * (1 - pi_i)
  grad <- (y - m_i * pi_i) / varianza * lnk$g(eta) * xi
  return(sum(grad))
}

beta_seq_beet <- seq(beta_hat_beet - 8, beta_hat_beet + 8, length.out = 2000)

score_valori_beet <- mapply(
  FUN      = Score_bin,
  beta     = beta_seq_beet,
  MoreArgs = list(xi = Beetles$logdose_m, y = Beetles$uccisi, m_i = Beetles$num)
)

df_score_beet <- data.frame(beta = beta_seq_beet, score = score_valori_beet)

ggplot(df_score_beet, aes(x = beta, y = score)) +
  geom_hline(yintercept = 0) +
  geom_line(col = "steelblue", linewidth = 1) +
  annotate("point", x = beta_hat_beet, y = 0, size = 2.5, col = "black") +
  labs(
    x        = expression(beta),
    y        = expression(U(beta)),
    subtitle = paste("Stima IRLS \u03b2:", round(beta_hat_beet, 3))
  ) +
  scale_x_continuous(breaks = seq(24, 40, by = 2)) +
  scale_y_continuous(breaks = seq(-2, 3, by = 0.5)) +
  tema_base


###############################################################################
# FIGURA 1.8 — Sensibilità delle previsioni al parametro
# Modello lineare vs modello logistico (modelli senza intercetta).
###############################################################################

x_sens <- seq(-10, 10, length.out = 300)

df_lineare <- data.frame(
  x    = rep(x_sens, 2),
  y    = c(1 * x_sens, 2 * x_sens),
  Beta = factor(rep(c("1", "2"), each = length(x_sens)))
)

p_lin <- ggplot(df_lineare, aes(x = x, y = y, color = Beta)) +
  geom_line(linewidth = 1.2) +
  labs(
    title = "Modello Lineare",
    y = "Valore atteso",
    x = expression("Predittore" ~ eta)
  ) +
  scale_color_manual(values = c("steelblue", "darkblue")) +
  theme(legend.position = "bottom") +
  coord_cartesian(xlim = c(-10, 10), ylim = c(-15, 15))

df_logistico <- data.frame(
  x    = rep(x_sens, 2),
  y    = c(plogis(1 * x_sens), plogis(3 * x_sens)),
  Beta = factor(rep(c("1", "3"), each = length(x_sens)))
)

p_log <- ggplot(df_logistico, aes(x = x, y = y, color = Beta)) +
  geom_line(linewidth = 1.2) +
  labs(
    title = "Modello Logistico",
    y = "Probabilità",
    x = expression("Predittore" ~ eta)
  ) +
  scale_color_manual(values = c("red", "darkred")) +
  theme(legend.position = "bottom")

grid.arrange(p_lin, p_log, nrow = 1, ncol = 2)


###############################################################################
#
#                          CAPITOLO 2
###############################################################################


###############################################################################
# FIGURA 2.1 — Distorsione dello SMV per un modello binomiale
# Confronto ML, Cordeiro-McCullagh (CM) e Kosmidis-Firth (KF).
# 5000 repliche, modello senza intercetta, beta = 0.8.
###############################################################################

set.seed(123)
true_beta_bias <- 0.8
m_bias <- 1
legame_bias <- "logit"
lnk_bias <- .links_GLM()[[legame_bias]]
ns_bias <- c(10, 15, 20, 30, 50, 75, 100)
N_bias <- 5000

risultati_bias <- data.frame()

for (n_g in ns_bias) {
  beta_ML <- rep(NA, N_bias)
  beta_CM <- rep(NA, N_bias)
  beta_KF <- rep(NA, N_bias)
  n_fail_ML <- 0

  for (i in 1:N_bias) {
    X_b <- matrix(rnorm(n_g), ncol = 1)
    eta_b <- X_b %*% true_beta_bias
    Pi_b <- lnk_bias$G(eta_b)
    y_b <- as.vector(rbinom(n_g, m_bias, Pi_b))

    X_data_b <- as.data.frame(X_b)
    X_data_b$y <- y_b

    # Stima KF (Firth)
    mod_KF <- IWLS_AJKF(cbind(y, m_bias - y) ~ . - 1,
      data = X_data_b, link = legame_bias,
      a = 1 / 2, maxit = 5000
    )
    if (!is.null(mod_KF) && isTRUE(mod_KF$converged)) {
      beta_KF[i] <- mod_KF$beta[1]
    }

    # Stima ML classica
    mod_ML <- glm(cbind(y, m_bias - y) ~ . - 1,
      data = X_data_b, family = binomial(link = legame_bias),
      maxit = 5000
    )

    mu_hat_b <- fitted(mod_ML)
    if (!mod_ML$converged || any(mu_hat_b >= 1 - 1e-12) || any(mu_hat_b <= 1e-12)) {
      n_fail_ML <- n_fail_ML + 1
      next
    }

    beta_ML[i] <- coef(mod_ML)[1]
    beta_CM[i] <- beta_ML[i] - C_M_correction(X_b, m_bias, beta_ML[i])
  }

  risultati_bias <- rbind(risultati_bias, data.frame(
    n = n_g,
    stimatore = c("ML", "CM", "KF"),
    bias = c(
      mean(beta_ML, na.rm = TRUE) - true_beta_bias,
      mean(beta_CM, na.rm = TRUE) - true_beta_bias,
      mean(beta_KF, na.rm = TRUE) - true_beta_bias
    )
  ))
}

risultati_bias$abs_bias <- abs(risultati_bias$bias)

# Stima pendenze in scala log-log
conv_ML <- lm(log(abs_bias) ~ log(n),
  data = filter(risultati_bias, stimatore == "ML")
)
conv_CM <- lm(log(abs_bias) ~ log(n),
  data = filter(risultati_bias, stimatore == "CM")
)
conv_KF <- lm(log(abs_bias) ~ log(n),
  data = filter(risultati_bias, stimatore == "KF")
)

ggplot(risultati_bias, aes(x = n, y = abs_bias, color = stimatore)) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  scale_x_log10(breaks = ns_bias) +
  scale_y_log10() +
  scale_color_manual(values = c("ML" = "red", "CM" = "steelblue", "KF" = "#0B2857")) +
  labs(
    x = "log(n)",
    y = expression("|" ~ E(hat(beta)[1] - beta[1]) ~ "|"),
    subtitle = paste0(
      "Distorsione dello stimatore di massima verosimiglianza\n",
      "Tasso stimato  n^(", round(coef(conv_ML)[2], 3), ")"
    )
  ) +
  tema_base


###############################################################################
# FIGURA 2.2a — Separazione completa nel dataset iris
# La specie setosa è perfettamente discriminata dalla lunghezza del petalo.
###############################################################################


taglio_iris <- 1 / 2 * (max(subset(iris$Petal.Length, iris$Species == "setosa")) + min(subset(iris$Petal.Length, iris$Species != "setosa")))

ggplot(iris, aes(
  x = Petal.Length,
  y = as.numeric(Species == "setosa"),
  col = Species
)) +
  geom_point(size = 2.3, alpha = 0.5) +
  geom_vline(
    xintercept = taglio_iris, col = "darkred",
    linetype = "dashed", linewidth = 0.7
  ) +
  scale_y_continuous(breaks = c(0, 1)) +
  labs(
    x        = "Lunghezza del petalo (cm)",
    y        = "1 (specie = setosa)",
    subtitle = "Separazione completa nel dataset iris"
  ) +
  tema_base


###############################################################################
# FIGURA 2.2b — Log-verosimiglianza divergente in presenza di separazione
# Modello logistico a un parametro (covariata centrata, senza intercetta).
###############################################################################

lnk_iris <- .links_GLM()[["logit"]]
X_iris <- matrix(iris$Petal.Length - taglio_iris, ncol = 1)
y_iris <- as.numeric(iris$Species != "setosa")
m_iris <- rep(1, nrow(X_iris))

pen_ll_iris <- function(beta, a) {
  eta <- as.vector(X_iris %*% beta)
  Pi <- pmin(pmax(lnk_iris$G(eta), 1e-10), 1 - 1e-10)
  w_m <- m_iris * lnk_iris$g(eta)^2 / (Pi * (1 - Pi))
  R <- qr.R(qr(sqrt(w_m) * X_iris))
  sum(y_iris * log(Pi) + (m_iris - y_iris) * log(1 - Pi)) +
    a * 2 * sum(log(abs(diag(R))))
}

griglia_iris <- seq(0.1, 15, length.out = 1000)

df_iris_ll <- data.frame(
  beta   = griglia_iris,
  loglik = sapply(griglia_iris, pen_ll_iris, a = 0)
)

ggplot(df_iris_ll, aes(x = beta, y = loglik)) +
  geom_line(linewidth = 1.2, col = "steelblue") +
  labs(
    x        = expression(beta),
    y        = expression(l^"*"(beta)),
    subtitle = "Divergenza log-verosimiglianza in presenza di separazione completa (iris)"
  ) +
  ylim(-20, 0) +
  xlim(1, 15) +
  tema_base


###############################################################################
# FIGURA 2.4 — Distorsione del primo ordine per la logistica
# b(theta)/n = -(1 - 2*pi) / (2*n*pi*(1 - pi)),  con n = 20.
###############################################################################

n_bias_plot <- 20
pi_seq_bias <- seq(0.02, 0.98, length.out = 400)

df_bias_plot <- data.frame(
  pi   = pi_seq_bias,
  bias = -(1 - 2 * pi_seq_bias) / (2 * n_bias_plot * pi_seq_bias * (1 - pi_seq_bias))
)

ggplot(df_bias_plot, aes(x = pi, y = bias)) +
  geom_hline(yintercept = 0, color = "black") +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey") +
  geom_line(color = "steelblue", linewidth = 1.2) +
  geom_point(aes(x = 0.5, y = 0), color = "brown", size = 2.5) +
  labs(
    x        = expression(pi),
    y        = expression(b(theta) / n),
    subtitle = paste0("Distorsione del primo ordine (n = ", n_bias_plot, ")")
  ) +
  tema_base


###############################################################################
# FIGURA 2.5 — Confronto ML vs KF vs CM: rapporto MSE
###############################################################################

set.seed(123)
true_beta_mse <- c(-1.0, 0.8)
p_mse <- length(true_beta_mse)
legame_mse <- "logit"
lnk_mse <- .links_GLM()[[legame_mse]]
N_mse <- 10000
ns_mse <- c(15, 25, 35, 50, 75, 100, 150, 200)
ms_mse <- c(1, 5, 20)

risultati_mse <- data.frame()
tempo_inizio_mse <- Sys.time()

for (m_g in ms_mse) {
  for (n_g in ns_mse) {
    # Matrici per raccogliere le stime
    stime_ML <- matrix(NA, nrow = N_mse, ncol = p_mse)
    stime_CM <- matrix(NA, nrow = N_mse, ncol = p_mse)
    stime_KF <- matrix(NA, nrow = N_mse, ncol = p_mse)
    n_fail_ML <- 0

    for (i in 1:N_mse) {
      X_mse <- cbind(rep(1, n_g), rnorm(n_g))
      eta_mse <- X_mse %*% true_beta_mse
      Pi_mse <- lnk_mse$G(eta_mse)
      y_mse <- as.vector(rbinom(n_g, m_g, Pi_mse))

      X_data_mse <- as.data.frame(X_mse)
      X_data_mse$y <- y_mse
      X_data_mse$m_g <- m_g

      # Stima KF (Firth)
      mod_KF_mse <- tryCatch(
        IWLS_AJKF(cbind(y, m_g - y) ~ . - 1,
          data = X_data_mse, link = legame_mse,
          a = 1 / 2, maxit = 5000
        ),
        error = function(e) NULL
      )
      if (!is.null(mod_KF_mse) && isTRUE(mod_KF_mse$converged)) {
        stime_KF[i, ] <- mod_KF_mse$beta
      }

      # Stima ML classica
      mod_ML_mse <- suppressWarnings(
        glm(cbind(y, m_g - y) ~ . - 1,
          data = X_data_mse, family = binomial(link = legame_mse),
          maxit = 5000
        )
      )

      mu_hat_mse <- fitted(mod_ML_mse)
      if (!mod_ML_mse$converged ||
        any(mu_hat_mse >= 1 - 1e-12) || any(mu_hat_mse <= 1e-12)) {
        n_fail_ML <- n_fail_ML + 1
        next
      }

      stime_ML[i, ] <- coef(mod_ML_mse)
      stime_CM[i, ] <- coef(mod_ML_mse) -
        C_M_correction(X_mse, m_g, coef(mod_ML_mse))
    }

    # Calcolo MSE e bias per ciascun parametro e stimatore
    nomi_param <- c("Intercetta", "Var_1")

    for (j in 1:p_mse) {
      mse_ml <- mean((stime_ML[, j] - true_beta_mse[j])^2, na.rm = TRUE)
      mse_cm <- mean((stime_CM[, j] - true_beta_mse[j])^2, na.rm = TRUE)
      mse_kf <- mean((stime_KF[, j] - true_beta_mse[j])^2, na.rm = TRUE)

      bias_ml <- mean(stime_ML[, j], na.rm = TRUE) - true_beta_mse[j]
      bias_cm <- mean(stime_CM[, j], na.rm = TRUE) - true_beta_mse[j]
      bias_kf <- mean(stime_KF[, j], na.rm = TRUE) - true_beta_mse[j]

      risultati_mse <- rbind(risultati_mse, data.frame(
        n = n_g, m = m_g, parametro = nomi_param[j],
        stimatore = c("ML", "CM", "K_Firth"),
        MSE = c(mse_ml, mse_cm, mse_kf),
        bias = c(bias_ml, bias_cm, bias_kf)
      ))
    }
  }
}

# Distorsione in scala log-log
ris_bias_mse <- risultati_mse %>%
  filter(parametro == "Var_1") %>%
  select(n, m, stimatore, bias) %>%
  mutate(
    abs_bias = abs(bias),
    N_eff = n * m
  )

# Stima delle pendenze per m = 1
conv_ML_mse <- lm(log(abs_bias) ~ log(n),
  data = filter(ris_bias_mse, stimatore == "ML", m == 1)
)
conv_CM_mse <- lm(log(abs_bias) ~ log(n),
  data = filter(ris_bias_mse, stimatore == "CM", m == 1)
)
conv_KF_mse <- lm(log(abs_bias) ~ log(n),
  data = filter(ris_bias_mse, stimatore == "K_Firth", m == 1)
)

ggplot(
  filter(ris_bias_mse, m == 1),
  aes(x = n, y = abs_bias, color = stimatore)
) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_manual(values = c(
    "ML" = "red", "CM" = "steelblue",
    "K_Firth" = "darkblue"
  )) +
  labs(
    x = "log(n)",
    y = expression(log("|" ~ E(hat(beta)[1] - beta[1]) ~ "|")),
    subtitle = paste0(
      "Distorsione dello stimatore di massima verosimiglianza\n",
      "Tasso stimato  n^(", round(coef(conv_ML_mse)[2], 3), ")"
    )
  ) +
  tema_base

# Rapporto MSE
df_ratio_mse <- risultati_mse %>%
  select(n, m, parametro, stimatore, MSE) %>%
  filter(!(m == 1 & n == 25)) %>%
  pivot_wider(names_from = stimatore, values_from = MSE) %>%
  mutate(
    Ratio_ML = ML / K_Firth,
    Ratio_CM = CM / K_Firth
  ) %>%
  pivot_longer(
    cols = c(Ratio_ML, Ratio_CM),
    names_to = "Confronto", values_to = "Efficienza_Relativa"
  ) %>%
  mutate(Confronto = ifelse(Confronto == "Ratio_ML",
    "ML vs Firth", "CM vs Firth"
  ))

ggplot(df_ratio_mse, aes(x = n, y = Efficienza_Relativa, color = Confronto)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.8) +
  facet_grid(parametro ~ m,
    scales = "free_y",
    labeller = labeller(m = label_both)
  ) +
  labs(x = "Numerosità campionaria", y = "Rapporto MSE") +
  tema_base

###############################################################################
#
#                          CAPITOLO 3
#      Penalizzazione di Jeffreys nei GLM per dati binari
#
###############################################################################


###############################################################################
# FIGURA 3.1 — Inverse delle funzioni legame per la regressione binaria
# Logit, Probit, cLogLog, LogLog, Cauchit.
###############################################################################

x_legami <- seq(-7, 7, length.out = 300)
links_fn <- .links_GLM()

df_legami <- data.frame(
  x = rep(x_legami, 5),
  y = c(
    links_fn[["logit"]]$G(x_legami),
    links_fn[["probit"]]$G(x_legami),
    links_fn[["cloglog"]]$G(x_legami),
    links_fn[["loglog"]]$G(x_legami),
    links_fn[["cauchit"]]$G(x_legami)
  ),
  Legame = factor(rep(c("Logit", "Probit", "Cloglog", "Loglog", "Cauchit"),
    each = length(x_legami)
  ))
)

ggplot(df_legami, aes(x = x, y = y, color = Legame)) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5) +
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.5) +
  scale_color_manual(values = c(
    "Logit" = "steelblue", "Probit" = "darkred", "Cloglog" = "darkgreen",
    "Loglog" = "darkorange", "Cauchit" = "purple"
  )) +
  labs(
    x     = expression("Predittore lineare" ~ eta),
    y     = expression("Probabilità prevista" ~ pi),
    color = "Funzione Legame"
  ) +
  coord_cartesian(xlim = c(-7, 7), ylim = c(0, 1)) +
  theme(legend.position = "bottom") +
  tema_base


###############################################################################
# FIGURA 3.2 — Log-verosimiglianza penalizzata (iris, separazione completa)
# Per a = 0 la stima diverge; per a > 0 il massimo esiste finito.
###############################################################################

valori_a <- c(0, 0.25, 0.5, 1)

df_pen_iris <- do.call(rbind, lapply(valori_a, function(a) {
  data.frame(
    a      = factor(a),
    beta   = griglia_iris,
    loglik = sapply(griglia_iris, pen_ll_iris, a = a)
  )
}))

ggplot(df_pen_iris, aes(x = beta, y = loglik, colour = a)) +
  geom_line(linewidth = 1) +
  labs(
    x        = expression(beta),
    y        = expression(l^"*"(beta)),
    colour   = "a",
    subtitle = "Effetto della penalizzazione in presenza di separazione completa (iris)"
  ) +
  ylim(-20, 0) +
  xlim(1, 15) +
  tema_base


###############################################################################
# FIGURA 3.3 — Effetto di contrazione delle probabilità penalizzate
# Curve di livello di log|X'W(pi)X| e frecce ML -> KF.
###############################################################################

# Funzione per calcolare il log-determinante dell'informazione
# su un modello con matrice disegno X = cbind(1, x) fissato (n = 2 osservazioni,
# m = 9, ispirato a.


x_shrink <- c(-1, 1)
X_shrink <- cbind(1, x_shrink)
m_shrink <- 9

log_det_info <- function(pi1, pi2) {
  w <- c(pi1 * (1 - pi1), pi2 * (1 - pi2)) * m_shrink
  XtWX <- t(X_shrink) %*% diag(w) %*% X_shrink
  det_val <- det(XtWX)
  if (det_val <= 0) {
    return(NA)
  } # separazione
  log(det_val)
}

grid_shrink <- expand.grid(
  pi1 = seq(0.01, 0.99, length.out = 150),
  pi2 = seq(0.01, 0.99, length.out = 150)
)
grid_shrink$logdet <- mapply(log_det_info, grid_shrink$pi1, grid_shrink$pi2)


combinazioni <- expand.grid(y1 = 0:m_shrink, y2 = 0:m_shrink) # tutti e 100 i dataset

frecce <- data.frame(
  pi1_ML = NA, pi2_ML = NA,
  pi1_KF = NA, pi2_KF = NA
)

for (k in 1:nrow(combinazioni)) {
  y1_k <- combinazioni$y1[k]
  y2_k <- combinazioni$y2[k]

  y_vec <- c(y1_k, y2_k)
  dat_k <- data.frame(
    y   = y_vec,
    m_g = rep(m_shrink, 2),
    V1  = x_shrink
  )

  # Stima ML
  mod_ml_k <- tryCatch(
    suppressWarnings(
      glm(cbind(y, m_g - y) ~ V1,
        data = dat_k,
        family = binomial(link = "logit")
      )
    ),
    error = function(e) NULL
  )

  if (!is.null(mod_ml_k) && mod_ml_k$converged) {
    pi_ml_k <- fitted(mod_ml_k)
    frecce[k, "pi1_ML"] <- pi_ml_k[1]
    frecce[k, "pi2_ML"] <- pi_ml_k[2]
  }

  # Stima KF (Firth, a = 1/2)
  mod_kf_k <- tryCatch(
    IWLS_AJKF(cbind(y, m_g - y) ~ V1,
      data = dat_k,
      link = "logit", a = 1 / 2, maxit = 5000
    ),
    error = function(e) NULL # non dovrebbe mai succedere, prcauzione
  )

  # non deve mai divergere
  if (!is.null(mod_kf_k) && isTRUE(mod_kf_k$converged)) {
    frecce[k, "pi1_KF"] <- mod_kf_k$probs[1]
    frecce[k, "pi2_KF"] <- mod_kf_k$probs[2]
  }
}

# Rimuoviamo le righe in cui ML o KF non sono finite (separazione)
frecce <- frecce[complete.cases(frecce), ]


# Livelli delle curve basati sui quantili della distribuzione dei valori,
livelli_contour <- quantile(grid_shrink$logdet,
  probs = seq(0.05, 0.95, length.out = 8),
  na.rm = TRUE
)

ggplot() +
  # Curve di livello con breaks ai quantili
  geom_contour(
    data = grid_shrink,
    aes(x = pi1, y = pi2, z = logdet),
    breaks = livelli_contour,
    colour = "darkgrey", linewidth = 0.4
  ) +
  # Frecce da ML a KF
  geom_segment(
    data = frecce,
    aes(
      x = pi1_ML, y = pi2_ML,
      xend = pi1_KF, yend = pi2_KF
    ),
    arrow = arrow(length = unit(0.12, "cm"), type = "closed"),
    colour = "black", linewidth = 0.4, alpha = 0.7
  ) +
  # Punto di massima informazione (1/2, 1/2)
  annotate("point",
    x = 0.5, y = 0.5,
    color = "red", size = 3, shape = 3, stroke = 1.5
  ) +
  labs(
    x        = expression(pi[1]),
    y        = expression(pi[2]),
    subtitle = "Effetto contrazione delle probabilità penalizzate\nLink: Logit"
  ) +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
  tema_base

###############################################################################
# FIGURA 3.4 — Distribuzione delle stime ML vs KF per diversi legami
# N = 1000 regressioni binomiali, n = 10, beta = 0.8.
###############################################################################
# raw_shrink <- readRDS("esperimenti_shrinkage_10")

set.seed(123)
true_beta_shrink <- 0.8
m_shrink_sim <- 1
N_shrink <- 1000
n_focus <- 10
a <- 1 / 2
links_shrink <- c("logit", "probit", "cloglog", "cauchit")

# Disegno fissato (generato una sola volta così da poter fissare il massimo teorico)
X_fixed_shrink <- rnorm(n_focus)

# Funzione dei pesi omega(eta) per il cLogLog (per il centro teorico)
omega_cloglog <- function(eta) exp(2 * eta) / (exp(exp(eta)) - 1)

# Centro teorico beta* = argmax log|X'W(beta)X| calcolato per il cLogLog numericamente
beta1_max_cloglog_val <- {
  log_det_info_clog <- function(b) {
    log(sum(omega_cloglog(b * X_fixed_shrink) * X_fixed_shrink^2))
  }
  optimize(log_det_info_clog, interval = c(-3, 3), maximum = TRUE)$maximum # prendiamo il max del log determinante per questa sizuazione senza intercetta
}

centro_teorico_shrink <- data.frame(
  link   = links_shrink,
  centro = c(0, 0, beta1_max_cloglog_val, 0)
)

# Funzione di simulazione per una cella (link, n)
simula_cella_shrink <- function(link, n_g, N, true_beta, m_g, lnk, X_fix) {
  beta_ML <- rep(NA_real_, N)
  beta_KF <- rep(NA_real_, N)
  sep_ML <- logical(N)

  eta_fix <- X_fix * true_beta
  Pi_fix <- lnk$G(eta_fix)

  for (i in seq_len(N)) {
    y_s <- rbinom(n_g, m_g, Pi_fix)
    dat_s <- data.frame(y = y_s, m_g = m_g, V1 = X_fix)

    # Stima KF (Firth-Kosmidis)
    mk_s <- tryCatch(
      IWLS_AJKF(cbind(y, m_g - y) ~ V1 - 1,
        data = dat_s,
        link = link, a, maxit = 5000, rho = 0.5
      ),
      error = function(e) NULL
    )
    if (!is.null(mk_s) && isTRUE(mk_s$converged)) {
      beta_KF[i] <- mk_s$beta[1]
    }

    # Stima ML classica
    mml_s <- tryCatch(
      suppressWarnings(glm(cbind(y, m_g - y) ~ V1 - 1,
        data = dat_s,
        family = binomial(link = link),
        control = glm.control(maxit = 5000)
      )),
      error = function(e) NULL
    )

    if (is.null(mml_s) || !mml_s$converged) {
      sep_ML[i] <- TRUE
      next
    }
    mu_s <- fitted(mml_s)
    if (any(mu_s >= 1 - 1e-8) || any(mu_s <= 1e-8) || abs(coef(mml_s)[1]) > 30) {
      sep_ML[i] <- TRUE
      next
    }
    beta_ML[i] <- coef(mml_s)[1]
  }

  data.frame(
    link = link, n = n_g,
    beta_ML = beta_ML, beta_KF = beta_KF, sep_ML = sep_ML
  )
}

# simulazione
raw_shrink <- list()
for (link_s in links_shrink) {
  lnk_s <- .links_GLM()[[link_s]]
  raw_shrink[[link_s]] <- simula_cella_shrink(
    link_s, n_focus, N_shrink, true_beta_shrink,
    m_shrink_sim, lnk_s, X_fixed_shrink
  )
}
raw_shrink <- bind_rows(raw_shrink)

# Tabella riassuntiva: distorsione e fattore di shrinkage
summary_shrink <- raw_shrink %>%
  left_join(centro_teorico_shrink, by = "link") %>%
  group_by(link) %>%
  summarise(
    bias_ML = round(mean(beta_ML, na.rm = TRUE) - true_beta_shrink, 3),
    bias_KF = round(mean(beta_KF, na.rm = TRUE) - true_beta_shrink, 3),
    shrinkage = {
      ok <- is.finite(beta_ML) & is.finite(beta_KF)
      ct <- first(centro)
      round(mean(abs(beta_KF[ok] - ct)) / mean(abs(beta_ML[ok] - ct)), 3)
    },
    .groups = "drop"
  )

# Densità delle stime ML vs KF per ciascun legame
plot_df_shrink <- raw_shrink %>%
  pivot_longer(c(beta_ML, beta_KF),
    names_to = "stimatore", values_to = "stima"
  ) %>%
  filter(is.finite(stima)) %>%
  left_join(centro_teorico_shrink, by = "link")

ggplot(plot_df_shrink, aes(x = stima, colour = stimatore, fill = stimatore)) +
  geom_density(alpha = 0.15) +
  geom_vline(xintercept = true_beta_shrink, linetype = 2) +
  geom_vline(aes(xintercept = centro), linetype = 2, colour = "red") +
  facet_wrap(~link, scales = "fixed") +
  scale_x_continuous(limits = c(-3, 3), breaks = seq(-3, 3, 1)) +
  labs(
    subtitle = sprintf("Distribuzione delle stime (n = %d)", n_focus),
    x = expression(hat(beta)),
    y = "Densità"
  ) +
  tema_base


###############################################################################
# FIGURA 3.5 — Effetto di contrazione delle probabilità stimate (mtcars)
# Regressione binaria con legame canonico: tipo di cambio ~ peso.
###############################################################################

mod_mtcars_ML <- glm(am ~ wt, data = mtcars, family = binomial(link = "logit"))
mod_mtcars_KF <- IWLS_AJKF(cbind(am, 1 - am) ~ wt,
  data = mtcars, link = "logit", a = 1 / 2
)

pi_ML <- fitted(mod_mtcars_ML)
pi_KF <- mod_mtcars_KF$probs
h_ii <- mod_mtcars_KF$leverage

df_mtcars_prob <- data.frame(
  wt           = mtcars$wt,
  am           = mtcars$am,
  ML           = pi_ML,
  KF           = pi_KF,
  classe_reale = factor(mtcars$am, labels = c("Automatico (0)", "Manuale (1)"))
)

df_mtcars_long <- df_mtcars_prob %>%
  pivot_longer(cols = c(ML, KF), names_to = "Metodo", values_to = "prob")

# Pannello sinistro: probabilità stimate
p_prob <- ggplot(df_mtcars_long, aes(x = wt, y = prob, color = Metodo)) +
  geom_line(linewidth = 1) +
  geom_point(
    data = df_mtcars_prob,
    aes(x = wt, y = am, col = classe_reale),
    inherit.aes = FALSE, size = 2, alpha = 0.7
  ) +
  scale_color_manual(values = c("black", "orange", "steelblue", "red")) +
  labs(
    x        = "Peso dell'auto (wt)",
    y        = "Probabilità Cambio Manuale",
    subtitle = "Effetto contrazione delle probabilità stimate"
  ) +
  tema_base

# Pannello destro: penalità additiva
df_penalty <- data.frame(
  pi_ML = pi_ML,
  penalty = h_ii * (1 / 2 - pi_ML),
  Metodo = as.numeric(mtcars$am)
)

p_penalty <- ggplot(df_penalty, aes(x = pi_ML, y = penalty, col = as.factor(Metodo))) +
  geom_point(size = 2, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0.5, linetype = "dotted", color = "grey") +
  scale_color_manual(values = c("black", "orange")) +
  labs(
    x        = "Probabilità stimata ML",
    y        = expression(h[ii](1 / 2 - pi[i])),
    subtitle = expression("Penalità additiva: " * h[ii](1 / 2 - pi[i]))
  ) +
  tema_base +
  theme(legend.position = "none")


grid.arrange(p_prob, p_penalty, nrow = 1)


###############################################################################
# FIGURA 3.6 — Ordine di convergenza empirico dell'algoritmo IWLS
###############################################################################
set.seed(123)
stima_ordine_conv <- function(err, tol = 1e-10) {
  k <- length(err)
  if (k < 4) { # scarta prime iterazioni
    return(c(q = NA, r = NA))
  }
  et <- err[-k]
  et1 <- err[-1]
  ok <- et > tol & et1 > tol
  if (sum(ok) < 3) { # solo se convergenza non troppo rapida
    return(c(q = NA, r = NA))
  }
  fit <- lm(log(et1[ok]) ~ log(et[ok])) # stima dell'ordine di covnergenza
  c(q = unname(coef(fit)[2]), r = unname(exp(coef(fit)[1])))
}

mc_ordine_conv <- function(M = 300, n = 50, a = 1 / 2,
                           legami = c("logit", "probit", "cloglog", "cauchit"),
                           beta_true = c(-1, 0.8)) {
  lnk_ref <- .links_GLM()[["logit"]]
  out <- vector("list", 0)

  for (m_rep in seq_len(M)) {
    
    x_mc <- rnorm(n)
    Pi_true_mc <- lnk_ref$G(beta_true[1] + beta_true[2] * x_mc)
    y_mc <- rbinom(n, 1, Pi_true_mc)
    dat_mc <- data.frame(succ = y_mc, fail = 1 - y_mc, x = x_mc)

    for (l in legami) {
      fit_mc <- tryCatch(
        IWLS_AJKF(cbind(succ, fail) ~ x, data = dat_mc, link = l, a = a),
        error = function(e) NULL
      )
      if (is.null(fit_mc) || !isTRUE(fit_mc$converged)) next

      err_mc <- apply(
        fit_mc$beta_log, 1,
        function(b) sqrt(sum((b - fit_mc$beta)^2))
      )
      s_mc <- stima_ordine_conv(err_mc)

      out[[length(out) + 1]] <- data.frame(
        rep = m_rep, link = l,
        q = s_mc["q"], r = s_mc["r"],
        iter = nrow(fit_mc$beta_log)
      )
    }
  }
  do.call(rbind, out) # unisci un un solo dataset
}

MC_conv <- mc_ordine_conv(M = 300, n = 50)

ggplot(MC_conv, aes(x = link, y = q)) +
  geom_boxplot(
    alpha = 0.5, outlier.size = 0.6,
    fill = "lightblue", width = 0.55
  ) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "black") +
  labs(
    x        = "Legame",
    y        = "Ordine di convergenza stimato",
    subtitle = "Distribuzione dell'ordine di convergenza empirico IWLS (n=50)"
  ) +
  tema_base


###############################################################################
#
#              FIGURE AGGIUNTIVE

###############################################################################
# Funzione dei pesi omega(pi) per i cinque legami
###############################################################################

legami_pesi <- c("logit", "probit", "cloglog", "loglog", "cauchit")
z_pesi <- seq(0.001, 0.999, length.out = 5000)

df_pesi <- do.call(rbind, lapply(legami_pesi, function(l) {
  lnk_p <- .links_GLM()[[l]]
  eta_p <- lnk_p$invG(z_pesi)
  w_p <- lnk_p$g(eta_p)^2 / (z_pesi * (1 - z_pesi))
  data.frame(link = l, z = z_pesi, w = w_p)
}))

df_massimi <- df_pesi %>%
  group_by(link) %>%
  slice_max(w, n = 1) %>%
  ungroup()

ggplot(df_pesi, aes(x = z, y = w, colour = link)) +
  geom_line(linewidth = 0.6) +
  geom_point(data = df_massimi, size = 2.5) +
  labs(
    x        = expression(pi),
    y        = expression(bar(omega)(pi)),
    colour   = "Legame",
    subtitle = "Funzione dei pesi rispetto alla media e relativi punti di massimo"
  ) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
  tema_base


###############################################################################
# Distribuzione a priori di Jeffreys per la regressione logistica univariata
###############################################################################

beta_jeff <- seq(-8, 8, length.out = 500)
pi_jeff <- plogis(beta_jeff)
jeffreys <- sqrt(pi_jeff * (1 - pi_jeff))

df_jeff <- data.frame(beta = beta_jeff, prior = jeffreys)

ggplot(df_jeff, aes(x = beta, y = prior)) +
  geom_area(fill = "darkblue", alpha = 0.15) +
  geom_line(color = "blue", linewidth = 1.2) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  annotate("text",
    x = 0.4, y = max(jeffreys) * 0.95,
    label = expression(beta == 0), color = "red", hjust = 0
  ) +
  labs(
    x = expression(beta),
    y = expression(sqrt(pi(beta) * (1 - pi(beta))))
  ) +
  tema_base


###############################################################################
# Penalità di Jeffreys: contour su dati mtcars
###############################################################################

x_heat <- as.numeric(scale(mtcars$wt))
X_heat <- cbind(1, x_heat)

compute_sqrt_det <- function(b0, b1) {
  eta_h <- b0 + b1 * x_heat
  pi_h <- plogis(eta_h)
  w_h <- pi_h * (1 - pi_h)
  XtWX_h <- t(X_heat) %*% (w_h * X_heat) / 2
  sqrt(det(XtWX_h))
}

grid_heat <- expand.grid(
  beta0 = seq(-4, 4, length.out = 150),
  beta1 = seq(-4, 4, length.out = 150)
)
grid_heat$value <- mapply(compute_sqrt_det, grid_heat$beta0, grid_heat$beta1)

ggplot(grid_heat, aes(x = beta0, y = beta1, fill = value)) +
  geom_raster() +
  scale_fill_gradient(low = "white", high = "darkgreen") +
  geom_contour(aes(z = value), color = "white", alpha = 0.7, linewidth = 0.3) +
  geom_point(aes(x = 0, y = 0), color = "red", size = 2.5, inherit.aes = FALSE) +
  labs(
    subtitle = "Penalità di Jeffreys",
    x = expression(beta[0]),
    y = expression(beta[1])
  ) +
  theme(legend.position = "none") +
  tema_base
