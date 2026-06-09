# =============================================================================
# 07_blrt_multivariado_v3.R
#
# Bootstrap LRT (BLRT) sobre a mistura NB multivariada V3
# (5 variáveis: consultas, ps, exames, internacoes, terapias).
#
# Testa apenas os dois pares críticos:
#   - H0: g = 3  vs  H1: g = 4   (validação da escolha BIC do pôster)
#   - H0: g = 4  vs  H1: g = 5   (confirmação de que g = 4 é teto)
#
# Reusa o EM customizado (independência condicional NB) do script
# comparar_variaveis.R, replicado aqui para autocontentação.
#
# Caching agressivo por par:
#   - fits_v3.rds        (ajustes em g = 3, 4, 5 nos dados reais)
#   - blrt_3vs4.rds      (49 réplicas bootstrap)
#   - blrt_4vs5.rds      (49 réplicas bootstrap)
#
# Saídas:
#   docs/assets/img/blrt_multivariado_v3/
#   src/03_analises/rds/blrt_multivariado_v3/
# =============================================================================


# ── 0. Pacotes ────────────────────────────────────────────────────────────────
pkgs <- c("dplyr", "ggplot2", "gt", "webshot2", "patchwork")
novos <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(novos) > 0) {
  install.packages(novos, repos = "https://cloud.r-project.org", quiet = TRUE)
}
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(gt); library(patchwork)
})


# ── 1. Raiz do projeto + diretórios ──────────────────────────────────────────
encontrar_raiz <- function(marcador = "poster_figures.R",
                            inicio = getwd(), max_subir = 6L) {
  dir_atual <- normalizePath(inicio, winslash = "/", mustWork = FALSE)
  for (i in seq_len(max_subir)) {
    if (file.exists(file.path(dir_atual, marcador))) return(dir_atual)
    pai <- dirname(dir_atual)
    if (pai == dir_atual) break
    dir_atual <- pai
  }
  stop("Nao encontrei a raiz do projeto. Rode setwd() para a raiz.")
}

PROJ_ROOT <- encontrar_raiz()
setwd(PROJ_ROOT)
cat("Raiz do projeto: ", PROJ_ROOT, "\n", sep = "")

IMG_DIR <- file.path(PROJ_ROOT, "docs/assets/img/blrt_multivariado_v3")
RDS_DIR <- file.path(PROJ_ROOT, "src/03_analises/rds/blrt_multivariado_v3")
dir.create(IMG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RDS_DIR, recursive = TRUE, showWarnings = FALSE)


# ── 2. Parâmetros do experimento ─────────────────────────────────────────────
set.seed(42)
CSV_PATH    <- "dados/dados_simulados.csv"
VARS_V3     <- c("consultas", "ps", "exames", "internacoes", "terapias")
K_GRID      <- c(3L, 4L, 5L)              # apenas o necessário p/ os pares
N_INIT_NB   <- 8L                          # reinícios no ajuste principal
N_INIT_BOOT <- 3L                          # reinícios em cada bootstrap fit
MAX_ITER    <- 200L
MAX_ITER_BT <- 100L                        # iter.max nos fits do bootstrap
TOL_NB      <- 1e-5
BLRT_B      <- 49L                         # bootstrap replicates por par


# ── 3. Funções do EM NB multivariado (replicadas de comparar_variaveis.R) ────

# Atualiza theta_hj por MLE univariado dado mu_hj e pesos tau_h
.update_theta <- function(y_j, tau_h, mu_hj) {
  neg_ll <- function(lt) {
    th <- exp(lt)
    -sum(tau_h * dnbinom(y_j, size = th, mu = mu_hj, log = TRUE))
  }
  exp(tryCatch(optimize(neg_ll, interval = c(-4, 8))$minimum,
               error = function(e) 0))
}

# EM para mistura NB multivariada (independência condicional)
# Init 1 = K-means em log1p(Y); demais = aleatório
nb_mix_em <- function(Y, g,
                      max_iter = MAX_ITER, tol = TOL_NB,
                      n_init = N_INIT_NB, seed = 42L) {
  n <- nrow(Y); p <- ncol(Y)
  Y <- as.matrix(Y)
  best_loglik <- -Inf; best_fit <- NULL

  for (init in seq_len(n_init)) {
    if (init == 1) {
      set.seed(seed)
      km  <- tryCatch(kmeans(log1p(Y), centers = g,
                              nstart = 10, iter.max = 100),
                       error = function(e) NULL)
      if (is.null(km)) next
      tau <- matrix(0.02 / max(g - 1, 1), n, g)
      for (i in seq_len(n)) tau[i, km$cluster[i]] <- 0.98
    } else {
      set.seed(seed + init)
      tau <- matrix(runif(n * g) + 0.01, n, g)
    }
    tau   <- tau / rowSums(tau)
    theta <- matrix(2, g, p)
    loglik_old <- -Inf
    converged  <- FALSE

    for (iter in seq_len(max_iter)) {
      pi_h <- colMeans(tau); n_h <- colSums(tau)
      mu   <- pmax(t(tau) %*% Y / n_h, 1e-6)
      for (h in seq_len(g))
        for (j in seq_len(p))
          theta[h, j] <- .update_theta(Y[, j], tau[, h], mu[h, j])
      theta <- pmax(theta, 0.05)

      log_f <- matrix(0, n, g)
      for (h in seq_len(g)) {
        log_f[, h] <- log(pi_h[h]) +
          rowSums(dnbinom(Y,
                          size = matrix(theta[h, ], n, p, byrow = TRUE),
                          mu   = matrix(mu[h, ],    n, p, byrow = TRUE),
                          log  = TRUE))
      }
      log_max  <- apply(log_f, 1, max)
      log_norm <- log(rowSums(exp(log_f - log_max))) + log_max
      loglik   <- sum(log_norm)
      tau_new  <- exp(log_f - log_norm)
      tau_new  <- tau_new / rowSums(tau_new)

      if (iter > 1 && abs(loglik - loglik_old) < tol) { converged <- TRUE; break }
      loglik_old <- loglik; tau <- tau_new
    }

    if (loglik > best_loglik) {
      best_loglik <- loglik
      k_params    <- 2L * g * p + (g - 1L)
      best_fit    <- list(
        g = g, pi = pi_h, mu = mu, theta = theta, tau = tau_new,
        loglik = loglik, n_iter = iter, converged = converged,
        classification = apply(tau_new, 1, which.max),
        BIC = -2 * loglik + k_params * log(n),
        AIC = -2 * loglik + 2 * k_params,
        np  = k_params
      )
    }
  }
  best_fit
}


# ── 4. Simulação paramétrica sob H0 (multivariada NB) ────────────────────────
# Dado um fit com (pi, mu [g x p], theta [g x p]):
#   1. Para cada i = 1..n_sim: sorteia componente h ~ Multinomial(pi)
#   2. Para j = 1..p: amostra y_ij ~ NB(mu[h,j], theta[h,j])  independentemente
simular_nb_v3 <- function(fit, n_sim) {
  g <- fit$g; p <- ncol(fit$mu)
  comp <- sample.int(g, size = n_sim, replace = TRUE, prob = fit$pi)
  Y    <- matrix(0L, n_sim, p)
  for (h in seq_len(g)) {
    idx <- which(comp == h)
    if (!length(idx)) next
    for (j in seq_len(p)) {
      Y[idx, j] <- rnbinom(length(idx),
                           size = fit$theta[h, j],
                           mu   = fit$mu[h, j])
    }
  }
  colnames(Y) <- colnames(fit$mu)
  Y
}


# ── 5. BLRT para um par (k0, k1) ─────────────────────────────────────────────
blrt_passo_v3 <- function(fit_k0, fit_k1, B = BLRT_B, nrep = N_INIT_BOOT,
                          seed = 1000L) {
  set.seed(seed)
  k0 <- fit_k0$g; k1 <- fit_k1$g
  lrt_obs <- 2 * (fit_k1$loglik - fit_k0$loglik)
  n_sim  <- nrow(fit_k0$tau)

  ll_best_k <- function(Y_b, k_val) {
    fit_b <- tryCatch(nb_mix_em(Y_b, g = k_val,
                                 max_iter = MAX_ITER_BT,
                                 n_init = nrep,
                                 seed = sample.int(.Machine$integer.max, 1)),
                       error = function(e) NULL)
    if (is.null(fit_b)) NA_real_ else fit_b$loglik
  }

  cat("\n    progresso: ")
  step <- max(1L, B %/% 10L)
  lrt_b <- vapply(seq_len(B), function(b) {
    Y_b <- simular_nb_v3(fit_k0, n_sim)
    ll0 <- ll_best_k(Y_b, k0)
    ll1 <- ll_best_k(Y_b, k1)
    out <- if (anyNA(c(ll0, ll1))) NA_real_ else 2 * (ll1 - ll0)
    if (is.na(out)) cat("x") else cat(".")
    if (b %% step == 0L) cat(sprintf("[%d/%d]", b, B))
    flush.console()
    out
  }, numeric(1))
  cat("\n")

  lrt_v <- lrt_b[!is.na(lrt_b)]
  pval  <- (sum(lrt_v >= lrt_obs) + 1) / (length(lrt_v) + 1)
  list(k0 = k0, k1 = k1,
       lrt_obs = round(lrt_obs, 2),
       pval = round(pval, 4),
       B_valido = length(lrt_v),
       lrt_boot = lrt_v)
}


# ── 6. Carregar dados V3 ─────────────────────────────────────────────────────
dados <- read.csv(CSV_PATH, stringsAsFactors = FALSE)
cat(sprintf("Dados carregados: %d obs, %d colunas\n",
            nrow(dados), ncol(dados)))

Y_v3      <- as.matrix(dados[, VARS_V3])
n         <- nrow(Y_v3)
grupo_int <- as.integer(factor(dados$grupo_latente,
                                levels = c("baixo_uso",
                                            "ambulatorial_coordenado",
                                            "agudo_hospitalar", "atipico")))

cat(sprintf("V3 multivariado: n = %d · p = %d (%s)\n",
            n, ncol(Y_v3), paste(VARS_V3, collapse = ", ")))


# ── 7. Ajustar modelos g = 3, 4, 5 nos dados reais (com cache) ───────────────
fits_cache <- file.path(RDS_DIR, "fits_v3.rds")

if (file.exists(fits_cache)) {
  cat("\n[cache] carregando fits_v3 de RDS\n")
  fits <- readRDS(fits_cache)
} else {
  cat("\n", strrep("=", 70), "\n", sep = "")
  cat("Ajustando modelos NB multivariados em g = 3, 4, 5\n")
  cat(strrep("=", 70), "\n\n", sep = "")

  fits <- list()
  for (g_val in K_GRID) {
    cat(sprintf("  g = %d ... ", g_val))
    t0 <- proc.time()["elapsed"]
    fit_g <- nb_mix_em(Y_v3, g = g_val, n_init = N_INIT_NB,
                        max_iter = MAX_ITER, seed = 42L + g_val)
    dt <- proc.time()["elapsed"] - t0
    cat(sprintf("logLik=%.1f | BIC=%.1f | iter=%d | %.0fs\n",
                fit_g$loglik, fit_g$BIC, fit_g$n_iter, dt))
    fits[[as.character(g_val)]] <- fit_g
  }
  saveRDS(fits, fits_cache)
}


# ── 8. BLRT g = 3 vs g = 4 (com cache) ───────────────────────────────────────
blrt34_cache <- file.path(RDS_DIR, "blrt_3vs4.rds")

cat("\n", strrep("=", 70), "\n", sep = "")
cat("BLRT g=3 vs g=4 (B=", BLRT_B, ")\n", sep = "")
cat(strrep("=", 70), sep = "")

if (file.exists(blrt34_cache)) {
  cat("\n[cache] carregando blrt_3vs4 de RDS\n")
  res_34 <- readRDS(blrt34_cache)
} else {
  cat("\n  Iniciando bootstrap paramétrico...")
  t0 <- proc.time()["elapsed"]
  res_34 <- blrt_passo_v3(fits[["3"]], fits[["4"]],
                          B = BLRT_B, nrep = N_INIT_BOOT, seed = 1003L)
  dt <- proc.time()["elapsed"] - t0
  cat(sprintf("  LRT_obs = %.1f | p = %.4f | B_validos = %d | %.0fs\n",
              res_34$lrt_obs, res_34$pval, res_34$B_valido, dt))
  saveRDS(res_34, blrt34_cache)
}


# ── 9. BLRT g = 4 vs g = 5 (com cache) ───────────────────────────────────────
blrt45_cache <- file.path(RDS_DIR, "blrt_4vs5.rds")

cat("\n", strrep("=", 70), "\n", sep = "")
cat("BLRT g=4 vs g=5 (B=", BLRT_B, ")\n", sep = "")
cat(strrep("=", 70), sep = "")

if (file.exists(blrt45_cache)) {
  cat("\n[cache] carregando blrt_4vs5 de RDS\n")
  res_45 <- readRDS(blrt45_cache)
} else {
  cat("\n  Iniciando bootstrap paramétrico...")
  t0 <- proc.time()["elapsed"]
  res_45 <- blrt_passo_v3(fits[["4"]], fits[["5"]],
                          B = BLRT_B, nrep = N_INIT_BOOT, seed = 1004L)
  dt <- proc.time()["elapsed"] - t0
  cat(sprintf("  LRT_obs = %.1f | p = %.4f | B_validos = %d | %.0fs\n",
              res_45$lrt_obs, res_45$pval, res_45$B_valido, dt))
  saveRDS(res_45, blrt45_cache)
}


# ── 10. Tabela de resultados ─────────────────────────────────────────────────
tab_blrt <- data.frame(
  g0         = c(res_34$k0, res_45$k0),
  g1         = c(res_34$k1, res_45$k1),
  LRT_obs    = c(res_34$lrt_obs, res_45$lrt_obs),
  p_valor    = c(res_34$pval, res_45$pval),
  B_validos  = c(res_34$B_valido, res_45$B_valido),
  rejeita_H0 = c(
    ifelse(res_34$pval <= 0.05, "Sim", "Nao"),
    ifelse(res_45$pval <= 0.05, "Sim", "Nao")
  )
)

write.csv(tab_blrt, file.path(RDS_DIR, "tab01_blrt_v3.csv"), row.names = FALSE)

gt01 <- gt(tab_blrt) %>%
  tab_header(
    title    = md("**Tabela 1.** Bootstrap LRT — V3 multivariado"),
    subtitle = md(sprintf("B = %d réplicas paramétricas · NB multivariada (5 variáveis) · n = %d",
                          BLRT_B, n))
  ) %>%
  cols_label(g0 = md("*g*"), g1 = md("*g*+1"),
             LRT_obs = "LRT observado", p_valor = "p-valor",
             B_validos = md("*B* válidos"), rejeita_H0 = "Rejeita H₀") %>%
  fmt_number(columns = LRT_obs, decimals = 1) %>%
  fmt_number(columns = p_valor, decimals = 4) %>%
  tab_style(style = cell_fill(color = "#f0f9ee"),
            locations = cells_body(rows = rejeita_H0 == "Sim")) %>%
  tab_style(style = cell_text(color = "white", weight = "bold"),
            locations = cells_column_labels(everything())) %>%
  tab_options(table.font.size = px(12),
              column_labels.background.color = "#2c3e50",
              column_labels.border.top.style = "hidden",
              heading.align = "left")

tryCatch({
  gtsave(gt01, file.path(IMG_DIR, "tab01_blrt_v3.png"),
         vwidth = 800, vheight = 240, zoom = 2, expand = 5)
  message("  -> salvo: tab01_blrt_v3.png")
}, error = function(e) message("  FALHA tab01: ", conditionMessage(e)))


# ── 11. Figura · distribuições bootstrap ─────────────────────────────────────
plot_blrt_um <- function(r) {
  dfh <- data.frame(lrt = r$lrt_boot)
  ggplot(dfh, aes(x = lrt)) +
    geom_histogram(bins = 25, fill = "grey70", colour = "white") +
    geom_vline(xintercept = r$lrt_obs, colour = "#e74c3c",
               linewidth = 1, linetype = "dashed") +
    annotate("text", x = r$lrt_obs, y = Inf, vjust = 2, hjust = -0.05,
             label = sprintf("LRT_obs = %.1f\np = %.4f", r$lrt_obs, r$pval),
             colour = "#c0392b", size = 3.2) +
    labs(title = sprintf("g = %d  vs  g = %d", r$k0, r$k1),
         x = "LRT bootstrap", y = "Frequência") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))
}

fig01 <- (plot_blrt_um(res_34) | plot_blrt_um(res_45)) +
  plot_annotation(
    title    = "Bootstrap LRT em V3 multivariado",
    subtitle = sprintf("B = %d réplicas paramétricas por par · NB multivariada (5 variáveis)",
                       BLRT_B),
    theme    = theme_minimal()
  )

ggsave(filename = file.path(IMG_DIR, "fig01_blrt_distribuicoes.png"),
       plot = fig01, width = 11, height = 4.5, dpi = 150, bg = "white")
message("  -> salvo: fig01_blrt_distribuicoes.png")


# ── 12. Resumo ────────────────────────────────────────────────────────────────
cat("\n", strrep("=", 70), "\n", sep = "")
cat("RESUMO — BLRT multivariado V3\n")
cat(strrep("=", 70), "\n", sep = "")
cat(sprintf("Variáveis: %s\n", paste(VARS_V3, collapse = ", ")))
cat(sprintf("n = %d · p = %d · B = %d\n", n, length(VARS_V3), BLRT_B))
cat("\n")
cat(sprintf("g = 3 vs g = 4: LRT = %.1f · p = %.4f · %s\n",
            res_34$lrt_obs, res_34$pval,
            ifelse(res_34$pval <= 0.05, "REJEITA H0 (g=4 melhor)",
                   "Nao rejeita")))
cat(sprintf("g = 4 vs g = 5: LRT = %.1f · p = %.4f · %s\n",
            res_45$lrt_obs, res_45$pval,
            ifelse(res_45$pval <= 0.05, "REJEITA H0",
                   "Nao rejeita (g=4 e suficiente)")))
cat("\nSaidas em: ", IMG_DIR, "\n", sep = "")
cat("  fig01_blrt_distribuicoes.png\n")
cat("  tab01_blrt_v3.png\n")
cat(strrep("=", 70), "\n", sep = "")
