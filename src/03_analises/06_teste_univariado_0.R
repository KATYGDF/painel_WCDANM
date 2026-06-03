# =============================================================================
# 06_teste_univariado_0.R
#
# TESTE UNIVARIADO 0 — BASELINE METODOLÓGICO
#
# Mistura Gaussiana univariada simulada em condições próximas ao ideal:
#   - Distribuição Normal (não NB, não sobredispersa)
#   - 4 componentes com médias bem espaçadas (mu = 10, 30, 50, 70)
#   - σ constante = 5 (separação de 4σ entre vizinhos)
#   - Pesos desbalanceados como nos dados reais: π = (0.50, 0.25, 0.20, 0.05)
#
# Mesmo protocolo dos Testes 1 e 2:
#   - Inicialização k-means + EM
#   - 5 critérios: AIC, BIC, ICL, AIC3, Bootstrap LRT
#   - Bootstrap de estabilidade via ARI
#
# Hipótese: em condições próximas das ideais (distribuição Gauss, μ bem
# espaçadas, σ moderado, pesos realistas), todos os critérios recuperam
# g = 4 com ARI > 0,95. O teste isola por exclusão a causa do fracasso
# do Teste 2: nem pesos desbalanceados, nem sobreposição leve são
# suficientes para quebrar a recuperação univariada.
#
# Saidas:
#   docs/assets/img/teste_univariado_0/
#   src/03_analises/rds/teste_univariado_0/
# =============================================================================


# ── 0. Pacotes ────────────────────────────────────────────────────────────────
pkgs <- c("dplyr", "tidyr", "ggplot2", "flexmix", "gt", "webshot2", "patchwork")
novos <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(novos) > 0) {
  install.packages(novos, repos = "https://cloud.r-project.org", quiet = TRUE)
}

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
  library(flexmix); library(gt); library(patchwork)
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

IMG_DIR <- file.path(PROJ_ROOT, "docs/assets/img/teste_univariado_0")
RDS_DIR <- file.path(PROJ_ROOT, "src/03_analises/rds/teste_univariado_0")
dir.create(IMG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RDS_DIR, recursive = TRUE, showWarnings = FALSE)


# ── 2. Estetica ───────────────────────────────────────────────────────────────
COR_GRUPO <- c(
  G1 = "#4a7c59",   # verde
  G2 = "#2980b9",   # azul
  G3 = "#e8a33c",   # laranja
  G4 = "#c0392b"    # vermelho
)

TEMA <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor  = element_blank(),
    strip.text        = element_text(face = "bold", size = 10),
    strip.background  = element_rect(fill = "#f0f0f0", colour = NA),
    plot.title        = element_text(face = "bold", size = 12),
    plot.subtitle     = element_text(colour = "grey40", size = 9),
    plot.caption      = element_text(colour = "grey55", size = 8, hjust = 0),
    legend.position   = "bottom"
  )

salvar_fig <- function(plot_obj, nome, largura = 8, altura = 5, dpi = 150) {
  ggsave(filename = file.path(IMG_DIR, paste0(nome, ".png")),
         plot = plot_obj, width = largura, height = altura,
         dpi = dpi, bg = "white")
  message("  -> salvo: ", nome, ".png")
}

salvar_gt <- function(gt_obj, nome, w_in = 8, h_in = 4) {
  tryCatch({
    gtsave(gt_obj, file.path(IMG_DIR, paste0(nome, ".png")),
           vwidth  = round(w_in * 100),
           vheight = round(h_in * 100),
           zoom = 2, expand = 5)
    message("  -> salvo: ", nome, ".png")
  }, error = function(e) message("  FALHA ", nome, ": ", conditionMessage(e)))
}


# ── 3. Parâmetros do experimento ──────────────────────────────────────────────
set.seed(1234)
N          <- 5000
PI_TRUE    <- c(0.50, 0.25, 0.20, 0.05)   # mesma proporcao real
MU_TRUE    <- c(10, 30, 50, 70)            # 4σ entre vizinhos
SIGMA_TRUE <- 5                            # σ constante
K_TRUE     <- 4
K_GRID     <- 2:6
BLRT_B     <- 29
BLRT_NREP  <- 2
BOOT_B     <- 50
G_LABEL    <- c("G1 · pop. base", "G2 · grupo médio",
                "G3 · grupo intensivo", "G4 · grupo raro")


# ── 4. Simulacao Gaussiana ────────────────────────────────────────────────────
grupo_real <- sample(1:4, size = N, replace = TRUE, prob = PI_TRUE)
valor      <- rnorm(N, mean = MU_TRUE[grupo_real], sd = SIGMA_TRUE)

df <- data.frame(
  grupo_real = factor(grupo_real, levels = 1:4,
                       labels = c("G1", "G2", "G3", "G4")),
  valor      = valor
)
grupo_int <- as.integer(df$grupo_real)

cat(sprintf("Dados simulados: n=%d, k_true=%d, sigma=%.1f\n",
            N, K_TRUE, SIGMA_TRUE))
cat(sprintf("  mu verdadeiros : %s\n", paste(MU_TRUE, collapse = ", ")))
cat(sprintf("  pi verdadeiros : %s\n", paste(PI_TRUE, collapse = ", ")))
cat("Distribuicao por grupo verdadeiro:\n")
print(table(df$grupo_real))


# ── 5. FIG 1 · Distribuição empírica ─────────────────────────────────────────
fig01 <- ggplot(df, aes(x = valor)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 60, fill = "grey80", colour = "white", alpha = 0.7) +
  geom_density(aes(fill = grupo_real, colour = grupo_real),
               alpha = 0.30, linewidth = 0.8) +
  scale_fill_manual(values = COR_GRUPO, name = "Grupo verdadeiro",
                    labels = G_LABEL) +
  scale_colour_manual(values = COR_GRUPO, name = "Grupo verdadeiro",
                      labels = G_LABEL) +
  labs(
    x        = "Valor (simulado: Gauss, μ = 10/30/50/70, σ = 5)",
    y        = "Densidade",
    title    = "Teste Univariado 0 — Mistura Gaussiana de 4 componentes",
    subtitle = sprintf("n = %d · pesos = (0,50; 0,25; 0,20; 0,05) · σ = 5 (constante)",
                       N),
    caption  = "Baseline metodológico — distribuição ideal, separação clara, rarefação realista (G4 = 5%)"
  ) +
  TEMA

salvar_fig(fig01, "fig01_dados_simulados", largura = 9, altura = 4.5)


# =============================================================================
# 6. Funções auxiliares (Gaussian mixture)
# =============================================================================

# Ajuste via flexmix Gaussiano (default - sem FLXMRnegbin)
fit_kmeans_em <- function(k_val, seed = 200) {
  set.seed(seed)
  if (k_val == 1L) {
    return(tryCatch(suppressWarnings(flexmix(
      valor ~ 1, data = df, k = 1L,
      control = list(iter.max = 500, tolerance = 1e-7, minprior = 0.005)
    )), error = function(e) NULL))
  }
  km <- kmeans(df$valor, centers = k_val, nstart = 25, iter.max = 100)
  tryCatch(suppressWarnings(flexmix(
    valor ~ 1, data = df, k = k_val,
    cluster = km$cluster,
    control = list(iter.max = 500, tolerance = 1e-7, minprior = 0.005)
  )), error = function(e) NULL)
}

extrair_metricas <- function(fit) {
  if (is.null(fit)) return(NULL)
  ll <- as.numeric(logLik(fit))
  np <- attr(logLik(fit), "df")
  list(
    k_efetivo = fit@k,
    logLik    = round(ll, 2),
    np        = np,
    AIC       = round(-2 * ll + 2 * np, 2),
    BIC       = round(-2 * ll + np * log(N), 2),
    ICL       = round(ICL(fit), 2),
    AIC3      = round(-2 * ll + 3 * np, 2)
  )
}

adj_rand_index <- function(x, y) {
  x <- as.integer(factor(x)); y <- as.integer(factor(y))
  tab <- table(x, y)
  a  <- sum(choose(tab, 2))
  rs <- rowSums(tab); cs <- colSums(tab); n <- sum(tab)
  b  <- sum(choose(rs, 2)) - a
  cv <- sum(choose(cs, 2)) - a
  exp_a <- (a + b) * (a + cv) / choose(n, 2)
  den   <- 0.5 * ((a + b) + (a + cv)) - exp_a
  if (den <= 0) return(1)
  (a - exp_a) / den
}

# Simulação paramétrica sob H0 (Gaussiana).
# parameters() do flexmix Gaussian retorna matriz 2xk: linha 1 = mean, linha 2 = sigma
simular_normal_mistura <- function(fit, n_sim) {
  pesos <- prior(fit)
  pars  <- parameters(fit)
  k     <- fit@k
  pars_m <- if (is.matrix(pars)) pars else matrix(pars, nrow = 2)
  mus    <- pars_m[1L, ]
  sigmas <- pars_m[2L, ]
  comp <- sample.int(k, size = n_sim, replace = TRUE, prob = pesos)
  y <- numeric(n_sim)
  for (h in seq_len(k)) {
    idx <- which(comp == h)
    if (length(idx))
      y[idx] <- rnorm(length(idx), mean = mus[h], sd = sigmas[h])
  }
  y
}

# BLRT paramétrico — H0: k vs H1: k+1.
blrt_passo <- function(fit_k0, fit_k1, B = BLRT_B, nrep = BLRT_NREP, seed = 1000) {
  set.seed(seed)
  k0 <- fit_k0@k; k1 <- fit_k1@k
  lrt_obs <- 2 * (as.numeric(logLik(fit_k1)) - as.numeric(logLik(fit_k0)))

  ll_best_k <- function(df_b, k_val) {
    best <- -Inf
    for (i in seq_len(nrep)) {
      m <- tryCatch(suppressWarnings(flexmix(
        valor ~ 1, data = df_b, k = k_val,
        control = list(iter.max = 200, minprior = 0.005, tolerance = 1e-6)
      )), error = function(e) NULL)
      if (!is.null(m) && m@k == k_val) {
        ll <- as.numeric(logLik(m))
        if (!is.na(ll) && ll > best) best <- ll
      }
    }
    if (is.infinite(best)) NA_real_ else best
  }

  cat("\n    progresso: ")
  step <- max(1L, B %/% 10L)
  lrt_b <- vapply(seq_len(B), function(b) {
    y_b  <- simular_normal_mistura(fit_k0, N)
    df_b <- data.frame(valor = y_b)
    ll0  <- ll_best_k(df_b, k0); ll1 <- ll_best_k(df_b, k1)
    out <- if (anyNA(c(ll0, ll1))) NA_real_ else 2 * (ll1 - ll0)
    if (is.na(out)) cat("x") else cat(".")
    if (b %% step == 0L) cat(sprintf("[%d/%d]", b, B))
    flush.console()
    out
  }, numeric(1))
  cat("\n  ")

  lrt_v <- lrt_b[!is.na(lrt_b)]
  pval  <- (sum(lrt_v >= lrt_obs) + 1) / (length(lrt_v) + 1)
  list(k0 = k0, k1 = k1, lrt_obs = round(lrt_obs, 2),
       pval = round(pval, 4), B_valido = length(lrt_v), lrt_boot = lrt_v)
}

# Pareamento húngaro
perms_all <- function(n) {
  if (n == 1L) return(matrix(1L, 1L, 1L))
  p <- perms_all(n - 1L)
  do.call(rbind, lapply(seq_len(n), function(i) cbind(i, p + (p >= i))))
}
hungarian_match <- function(conf_mat) {
  nr <- nrow(conf_mat); nc <- ncol(conf_mat); k <- max(nr, nc)
  sq <- matrix(0L, k, k); sq[seq_len(nr), seq_len(nc)] <- conf_mat
  perms <- perms_all(k)
  scores <- apply(perms, 1L, function(p) sum(diag(sq[, p, drop = FALSE])))
  a <- as.integer(perms[which.max(scores), ])
  ifelse(a[seq_len(nr)] <= nc, a[seq_len(nr)], NA_integer_)
}

remap_hungarian <- function(real_int, est_int) {
  est_levels <- sort(unique(est_int))
  conf <- table(real = real_int,
                est  = factor(est_int, levels = est_levels))
  assign_h <- hungarian_match(conf)
  remap <- setNames(rep(NA_integer_, length(est_levels)), as.character(est_levels))
  for (g in seq_along(assign_h)) {
    if (!is.na(assign_h[g]) && assign_h[g] <= length(est_levels))
      remap[as.character(est_levels[assign_h[g]])] <- g
  }
  as.integer(remap[as.character(est_int)])
}


# =============================================================================
# 7. MODULO A · Ajuste em K_GRID
# =============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("MODULO A · k-means + EM (Gaussian), k = ", paste(K_GRID, collapse = ", "), "\n", sep = "")
cat(strrep("=", 70), "\n\n", sep = "")

fits_cache <- file.path(RDS_DIR, "fits.rds")
tab_cache  <- file.path(RDS_DIR, "tab_init.rds")

if (file.exists(fits_cache) && file.exists(tab_cache)) {
  cat("  [cache] carregando fits + tab_init de RDS\n")
  fits     <- readRDS(fits_cache)
  tab_init <- readRDS(tab_cache)
} else {
  fits     <- setNames(vector("list", length(K_GRID)), as.character(K_GRID))
  tab_init <- list()
  for (k_val in K_GRID) {
    cat(sprintf("  kmeans+EM (Gauss) k=%d ... ", k_val))
    t0 <- proc.time()["elapsed"]
    fit_kv <- fit_kmeans_em(k_val)
    dt <- proc.time()["elapsed"] - t0
    fits[[as.character(k_val)]] <- fit_kv
    m <- extrair_metricas(fit_kv)
    if (is.null(m)) {
      cat(sprintf("FALHOU (%.1fs)\n", dt))
      tab_init[[length(tab_init) + 1]] <- data.frame(
        k = k_val, k_efetivo = NA, logLik = NA, np = NA,
        AIC = NA, BIC = NA, ICL = NA, AIC3 = NA, tempo_s = round(dt, 1)
      )
    } else {
      cat(sprintf("logLik=%.1f | BIC=%.1f | np=%d | %.1fs\n",
                  m$logLik, m$BIC, m$np, dt))
      tab_init[[length(tab_init) + 1]] <- data.frame(
        k = k_val, k_efetivo = m$k_efetivo, logLik = m$logLik, np = m$np,
        AIC = m$AIC, BIC = m$BIC, ICL = m$ICL, AIC3 = m$AIC3,
        tempo_s = round(dt, 1)
      )
    }
  }
  tab_init <- do.call(rbind, tab_init)
  saveRDS(fits, fits_cache)
  saveRDS(tab_init, tab_cache)
}

write.csv(tab_init, file.path(RDS_DIR, "tab01_criterios.csv"),
          row.names = FALSE)


# ── FIG 2 · Critérios normalizados ───────────────────────────────────────────
norm01 <- function(v) {
  r <- range(v, na.rm = TRUE)
  if (diff(r) == 0) return(rep(0, length(v)))
  (v - r[1]) / diff(r)
}

valid <- !is.na(tab_init$logLik)
k_aic  <- tab_init$k[valid][which.min(tab_init$AIC[valid])]
k_bic  <- tab_init$k[valid][which.min(tab_init$BIC[valid])]
k_icl  <- tab_init$k[valid][which.min(tab_init$ICL[valid])]
k_aic3 <- tab_init$k[valid][which.min(tab_init$AIC3[valid])]

cat(sprintf("\nk* por criterio: AIC=%d | BIC=%d | ICL=%d | AIC3=%d\n",
            k_aic, k_bic, k_icl, k_aic3))

tab_long <- data.frame(
  k    = rep(tab_init$k[valid], 4),
  crit = rep(c("AIC", "BIC", "ICL", "AIC3"), each = sum(valid)),
  vn   = c(norm01(tab_init$AIC[valid]),
            norm01(tab_init$BIC[valid]),
            norm01(tab_init$ICL[valid]),
            norm01(tab_init$AIC3[valid]))
)

fig02 <- ggplot(tab_long, aes(x = k, y = vn, colour = crit, shape = crit)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 3) +
  geom_vline(xintercept = K_TRUE, linetype = "dashed",
             colour = "grey55", linewidth = 0.5) +
  annotate("text", x = K_TRUE, y = 1.05, label = "k verdadeiro = 4",
           colour = "grey45", size = 3, vjust = 0) +
  scale_x_continuous(breaks = K_GRID) +
  scale_colour_manual(values = c(AIC = "#27ae60", BIC = "#e74c3c",
                                  ICL = "#2980b9", AIC3 = "#8e44ad"),
                      name = "Critério") +
  scale_shape_manual(values = c(AIC = 15, BIC = 16, ICL = 17, AIC3 = 18),
                     name = "Critério") +
  labs(
    title    = "Critérios de seleção normalizados [0 = melhor, 1 = pior]",
    subtitle = sprintf("k* — AIC: %d | BIC: %d | ICL: %d | AIC3: %d",
                       k_aic, k_bic, k_icl, k_aic3),
    x = "Número de componentes (k)",
    y = "Valor normalizado"
  ) +
  TEMA

salvar_fig(fig02, "fig02_criterios_por_k", largura = 8, altura = 5)


# =============================================================================
# 8. MODULO B · Bootstrap LRT
# =============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("MODULO B · Bootstrap LRT (H0: k vs H1: k+1)\n")
cat(strrep("=", 70), "\n\n", sep = "")

k_blrt_range <- K_GRID[-length(K_GRID)]
blrt_cache <- file.path(RDS_DIR, "resultados_blrt.rds")

if (file.exists(blrt_cache)) {
  cat("  [cache] carregando resultados_blrt de RDS\n")
  resultados_blrt <- readRDS(blrt_cache)
  k_blrt_escolhido <- NA_integer_
  for (k0 in sort(as.integer(names(resultados_blrt)))) {
    if (is.na(k_blrt_escolhido) && resultados_blrt[[as.character(k0)]]$pval > 0.05)
      k_blrt_escolhido <- k0
  }
} else {
  resultados_blrt  <- list()
  k_blrt_escolhido <- NA_integer_
  cat(sprintf("Configuracao: B=%d bootstrap | nrep=%d reinicios por ajuste\n",
              BLRT_B, BLRT_NREP))
  for (k0 in k_blrt_range) {
    k1 <- k0 + 1L
    fit_k0 <- fits[[as.character(k0)]]
    fit_k1 <- fits[[as.character(k1)]]
    if (is.null(fit_k0) || is.null(fit_k1)) next
    cat(sprintf("  BLRT k=%d vs k=%d ... ", k0, k1))
    t0 <- proc.time()["elapsed"]
    res <- blrt_passo(fit_k0, fit_k1, B = BLRT_B, nrep = BLRT_NREP,
                      seed = 1000 + k0)
    dt <- proc.time()["elapsed"] - t0
    cat(sprintf("LRT=%.1f | p=%.4f | B_validos=%d | %.0fs\n",
                res$lrt_obs, res$pval, res$B_valido, dt))
    resultados_blrt[[as.character(k0)]] <- res
    if (is.na(k_blrt_escolhido) && res$pval > 0.05)
      k_blrt_escolhido <- k0
  }
  saveRDS(resultados_blrt, blrt_cache)
}
if (is.na(k_blrt_escolhido))
  k_blrt_escolhido <- max(k_blrt_range) + 1L

cat(sprintf("\nk* por BLRT: %d\n", k_blrt_escolhido))

tab_blrt <- do.call(rbind, lapply(resultados_blrt, function(r) {
  data.frame(k0 = r$k0, k1 = r$k1,
             LRT_obs = r$lrt_obs, p_valor = r$pval,
             B_validos = r$B_valido,
             rejeita_H0 = ifelse(r$pval <= 0.05, "Sim", "Nao"))
}))
write.csv(tab_blrt, file.path(RDS_DIR, "tab02_blrt.csv"), row.names = FALSE)

plot_blrt_um <- function(r) {
  dfh <- data.frame(lrt = r$lrt_boot)
  ggplot(dfh, aes(x = lrt)) +
    geom_histogram(bins = 25, fill = "grey70", colour = "white") +
    geom_vline(xintercept = r$lrt_obs, colour = "#e74c3c",
               linewidth = 1, linetype = "dashed") +
    annotate("text", x = r$lrt_obs, y = Inf, vjust = 2, hjust = -0.05,
             label = sprintf("LRT_obs\np = %.4f", r$pval),
             colour = "#c0392b", size = 3) +
    labs(title = sprintf("k=%d vs k=%d", r$k0, r$k1),
         x = "LRT bootstrap", y = "Frequência") +
    TEMA + theme(plot.title = element_text(size = 11))
}
paineis_blrt <- lapply(resultados_blrt, plot_blrt_um)
if (length(paineis_blrt) >= 1) {
  fig03 <- wrap_plots(paineis_blrt, ncol = 2) +
    plot_annotation(
      title = "Bootstrap LRT — distribuição sob H0",
      subtitle = sprintf("B=%d réplicas paramétricas (Gauss)", BLRT_B),
      theme = TEMA
    )
  salvar_fig(fig03, "fig03_blrt_distribuicoes",
             largura = 10,
             altura = 3.5 * ceiling(length(paineis_blrt) / 2))
}


# =============================================================================
# 9. MODULO C · Estabilidade
# =============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("MODULO C · Selecao final + estabilidade\n")
cat(strrep("=", 70), "\n\n", sep = "")

votos <- c(AIC = k_aic, BIC = k_bic, ICL = k_icl, AIC3 = k_aic3,
           BLRT = k_blrt_escolhido)
k_final <- as.integer(names(sort(table(votos), decreasing = TRUE))[1])
cat(sprintf("Votos por k:  AIC=%d | BIC=%d | ICL=%d | AIC3=%d | BLRT=%d\n",
            k_aic, k_bic, k_icl, k_aic3, k_blrt_escolhido))
cat(sprintf("k_final (consenso): %d\n", k_final))

fit_final <- fits[[as.character(k_final)]]

ari_cache <- file.path(RDS_DIR, "ari_boot.rds")
if (file.exists(ari_cache)) {
  cat(sprintf("\n[cache] carregando ari_vec (B=%d) de RDS\n", BOOT_B))
  ari_vec <- readRDS(ari_cache)
} else {
  cat(sprintf("\nBootstrap de estabilidade (B=%d) ...\n  progresso: ", BOOT_B))
  t0_boot <- proc.time()["elapsed"]
  set.seed(7777)
  cl_ref <- clusters(fit_final)
  step_b <- max(1L, BOOT_B %/% 10L)
  ari_vec <- vapply(seq_len(BOOT_B), function(b) {
    idx <- sample.int(N, replace = TRUE)
    df_b <- data.frame(valor = valor[idx])
    m_best <- NULL; ll_best <- -Inf
    for (i in seq_len(BLRT_NREP)) {
      m_i <- tryCatch(suppressWarnings(flexmix(
        valor ~ 1, data = df_b, k = k_final,
        control = list(iter.max = 200, minprior = 0.005, tolerance = 1e-6)
      )), error = function(e) NULL)
      if (!is.null(m_i) && m_i@k == k_final) {
        ll <- as.numeric(logLik(m_i))
        if (!is.na(ll) && ll > ll_best) { ll_best <- ll; m_best <- m_i }
      }
    }
    out <- if (is.null(m_best)) NA_real_
           else adj_rand_index(cl_ref[idx], clusters(m_best))
    if (is.na(out)) cat("x") else cat(".")
    if (b %% step_b == 0L) cat(sprintf("[%d/%d]", b, BOOT_B))
    flush.console()
    out
  }, numeric(1))
  dt_boot <- proc.time()["elapsed"] - t0_boot
  cat(sprintf("\n  Concluido em %.0fs\n", dt_boot))
  saveRDS(ari_vec, ari_cache)
}

ari_valid <- ari_vec[!is.na(ari_vec)]
ari_med <- if (length(ari_valid)) mean(ari_valid) else NA
ari_ic  <- if (length(ari_valid)) {
  quantile(ari_valid, c(0.025, 0.975))
} else {
  c(NA, NA)
}
cat(sprintf("ARI medio = %.3f | IC95%% = [%.3f, %.3f]\n",
            ari_med, ari_ic[1], ari_ic[2]))

if (length(ari_valid)) {
  fig04 <- ggplot(data.frame(ari = ari_valid), aes(x = ari)) +
    geom_histogram(bins = 25, fill = "#3498db", colour = "white", alpha = 0.85) +
    geom_vline(xintercept = ari_med, colour = "#e74c3c",
               linewidth = 1, linetype = "dashed") +
    geom_vline(xintercept = 0.80, colour = "#f39c12",
               linewidth = 0.7, linetype = "dotted") +
    annotate("text", x = 0.80, y = Inf, vjust = 1.5, hjust = 1.1,
             label = "limiar 0,80", colour = "#d68910", size = 3) +
    labs(
      title = sprintf("Estabilidade do modelo k=%d (bootstrap)", k_final),
      subtitle = sprintf("ARI medio=%.3f | IC95%% [%.3f, %.3f] | B=%d",
                         ari_med, ari_ic[1], ari_ic[2], length(ari_valid)),
      x = "Adjusted Rand Index", y = "Frequência"
    ) + TEMA
  salvar_fig(fig04, "fig04_estabilidade_ari", largura = 8, altura = 4.5)
}


# =============================================================================
# 10. MODULO D · Avaliacao em k = K_TRUE
# =============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat(sprintf("MODULO D · Avaliacao em k = K_TRUE = %d\n", K_TRUE))
cat(strrep("=", 70), "\n\n", sep = "")

fit_K   <- fits[[as.character(K_TRUE)]]
est_K   <- as.integer(clusters(fit_K))
est_par <- remap_hungarian(grupo_int, est_K)
ari_K   <- adj_rand_index(grupo_int, est_K)
acu_K   <- mean(est_par == grupo_int, na.rm = TRUE)
conf_K  <- table(real = grupo_int, est = est_K)
pur_K   <- sum(apply(conf_K, 2L, max)) / sum(conf_K)
f1_grp  <- vapply(seq_len(4), function(g) {
  tp <- sum(grupo_int == g & est_par == g, na.rm = TRUE)
  fp <- sum(grupo_int != g & est_par == g, na.rm = TRUE)
  fn <- sum(grupo_int == g & est_par != g, na.rm = TRUE)
  if (tp + fp == 0 || tp + fn == 0) return(0)
  prec <- tp / (tp + fp); rec <- tp / (tp + fn)
  if (prec + rec == 0) 0 else 2 * prec * rec / (prec + rec)
}, numeric(1))
sens_grp <- vapply(seq_len(4), function(g) {
  n_g <- sum(grupo_int == g)
  if (n_g == 0) NA_real_ else sum(est_par == g & grupo_int == g) / n_g
}, numeric(1))

cat(sprintf("ARI = %.3f | Acuracia = %.3f | Pureza = %.3f | F1 macro = %.3f\n",
            ari_K, acu_K, pur_K, mean(f1_grp)))
cat("Sensibilidade por grupo:\n")
for (g in 1:4) cat(sprintf("  G%d: %.3f\n", g, sens_grp[g]))


# Tabela 3 · métricas globais
tab03 <- data.frame(
  Métrica = c("ARI", "Acurácia", "Pureza", "F1 macro"),
  Valor   = round(c(ari_K, acu_K, pur_K, mean(f1_grp)), 3)
)
write.csv(tab03, file.path(RDS_DIR, "tab03_metricas_globais.csv"),
          row.names = FALSE)

gt03 <- gt(tab03) %>%
  tab_header(
    title    = md(sprintf("**Tabela 3.** Métricas globais no modelo k = %d (Gauss)", K_TRUE)),
    subtitle = md("Pareamento via algoritmo húngaro.")
  ) %>%
  fmt_number(columns = Valor, decimals = 3) %>%
  tab_style(style = cell_text(color = "white", weight = "bold"),
            locations = cells_column_labels(everything())) %>%
  tab_options(table.font.size = px(12),
              column_labels.background.color = "#2c3e50",
              column_labels.border.top.style = "hidden",
              heading.align = "left")
salvar_gt(gt03, "tab03_metricas_globais", w_in = 6, h_in = 2.5)


# Tabela 4 · parâmetros estimados (mu, sigma, pi) por componente
pars_K <- parameters(fit_K)
pars_m <- if (is.matrix(pars_K)) pars_K else matrix(pars_K, nrow = 2)
mus_est    <- pars_m[1L, ]
sigmas_est <- pars_m[2L, ]
prior_est  <- as.numeric(prior(fit_K))

est_levels <- sort(unique(est_K))
conf_full <- table(real = grupo_int,
                   est = factor(est_K, levels = est_levels))
assign_h <- hungarian_match(conf_full)

mu_aligned    <- vapply(assign_h, function(j) if (is.na(j)) NA_real_ else mus_est[j],    numeric(1))
sigma_aligned <- vapply(assign_h, function(j) if (is.na(j)) NA_real_ else sigmas_est[j], numeric(1))
pi_aligned    <- vapply(assign_h, function(j) if (is.na(j)) NA_real_ else prior_est[j],  numeric(1))

tab04 <- data.frame(
  Grupo         = G_LABEL,
  n_real        = as.integer(table(grupo_int)),
  mu_real       = MU_TRUE,
  mu_est        = round(mu_aligned, 2),
  sigma_real    = SIGMA_TRUE,
  sigma_est     = round(sigma_aligned, 2),
  pi_real       = PI_TRUE,
  pi_est        = round(pi_aligned, 3),
  sensibilidade = round(sens_grp, 3),
  F1            = round(f1_grp, 3)
)
write.csv(tab04, file.path(RDS_DIR, "tab04_parametros_diagnostico.csv"),
          row.names = FALSE)

gt04 <- gt(tab04) %>%
  tab_header(
    title    = md(sprintf("**Tabela 4.** Recuperação paramétrica + diagnóstico por grupo (k = %d)",
                          K_TRUE)),
    subtitle = md("Alinhamento via algoritmo húngaro. Reais vs. estimados.")
  ) %>%
  cols_label(
    Grupo = "Grupo", n_real = md("*n* real"),
    mu_real = md("μ real"), mu_est = md("μ est."),
    sigma_real = md("σ real"), sigma_est = md("σ est."),
    pi_real = md("π real"), pi_est = md("π est."),
    sensibilidade = "Sensib.", F1 = "F1"
  ) %>%
  fmt_number(columns = c(mu_real, mu_est, sigma_real, sigma_est), decimals = 2) %>%
  fmt_number(columns = c(pi_real, pi_est, sensibilidade, F1), decimals = 3) %>%
  tab_style(style = cell_text(color = "white", weight = "bold"),
            locations = cells_column_labels(everything())) %>%
  tab_options(table.font.size = px(11),
              column_labels.background.color = "#2c3e50",
              column_labels.border.top.style = "hidden",
              heading.align = "left")
salvar_gt(gt04, "tab04_parametros_diagnostico", w_in = 11, h_in = 3.5)


# Fig 5 · densidade ajustada
ord_idx <- order(mus_est)
mus_sorted    <- mus_est[ord_idx]
sigmas_sorted <- sigmas_est[ord_idx]
prior_sorted  <- prior_est[ord_idx]
grade_x <- seq(min(valor), max(valor), length.out = 400)
dens_comp <- do.call(rbind, lapply(seq_along(mus_sorted), function(h) {
  data.frame(
    x = grade_x,
    densidade = prior_sorted[h] * dnorm(grade_x, mean = mus_sorted[h],
                                         sd = sigmas_sorted[h]),
    componente = sprintf("C%d (μ=%.1f, σ=%.2f)",
                         h, mus_sorted[h], sigmas_sorted[h])
  )
}))
dens_total <- data.frame(
  x = grade_x,
  densidade = sapply(grade_x, function(xv) {
    sum(prior_sorted * dnorm(xv, mean = mus_sorted, sd = sigmas_sorted))
  })
)

fig05 <- ggplot() +
  geom_histogram(data = df, aes(x = valor, y = after_stat(density)),
                 bins = 60, fill = "grey80", colour = "white", alpha = 0.7) +
  geom_line(data = dens_comp,
            aes(x = x, y = densidade, colour = componente), linewidth = 0.9) +
  geom_line(data = dens_total, aes(x = x, y = densidade),
            colour = "black", linewidth = 1, linetype = "dashed") +
  scale_colour_brewer(palette = "Set1") +
  labs(
    title    = sprintf("Modelo ajustado (k = %d, Gauss, k-means + EM)", K_TRUE),
    subtitle = sprintf("ARI = %.3f · Acurácia = %.3f · F1 macro = %.3f",
                       ari_K, acu_K, mean(f1_grp)),
    x = "Valor", y = "Densidade",
    colour = "Componente estimado"
  ) + TEMA

salvar_fig(fig05, "fig05_modelo_ajustado", largura = 10, altura = 5)


# Fig 6 · matriz de confusão
conf_match <- conf_K[, assign_h]
colnames(conf_match) <- G_LABEL
rownames(conf_match) <- G_LABEL

conf_df <- as.data.frame.table(conf_match, responseName = "n")
names(conf_df)[1:2] <- c("real", "est")
conf_df$pct <- 100 * conf_df$n / N

fig06 <- ggplot(conf_df, aes(x = est, y = real, fill = n)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", n, pct)),
            size = 3.2, colour = "grey15") +
  scale_fill_gradient(low = "#eef5fb", high = "#2980b9", guide = "none") +
  scale_x_discrete(position = "top") +
  labs(
    x = "Componente estimado",
    y = "Grupo real",
    title = sprintf("Matriz de confusão — k = %d (Gauss)", K_TRUE),
    subtitle = sprintf("Acurácia = %.1f%%", 100 * acu_K)
  ) +
  TEMA + theme(panel.grid = element_blank(),
               axis.text.x = element_text(face = "bold"),
               axis.text.y = element_text(face = "bold"))

salvar_fig(fig06, "fig06_matriz_confusao", largura = 9, altura = 5.5)


# =============================================================================
# 11. RESUMO
# =============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("RESUMO — Teste Univariado 0\n")
cat(strrep("=", 70), "\n", sep = "")
cat(sprintf("Saidas em: %s\n", IMG_DIR))
cat("Figuras: fig01..fig06 (distribuicoes/criterios/blrt/estabilidade/modelo/confusao)\n")
cat("Tabelas: tab01..tab04 (criterios/blrt/metricas/parametros)\n")
cat(sprintf("\nVotos por k: AIC=%d | BIC=%d | ICL=%d | AIC3=%d | BLRT=%d\n",
            k_aic, k_bic, k_icl, k_aic3, k_blrt_escolhido))
cat(sprintf("k_final (consenso) = %d\n", k_final))
cat(sprintf("Em k=%d: ARI=%.3f | Acuracia=%.3f | Pureza=%.3f | F1 macro=%.3f\n",
            K_TRUE, ari_K, acu_K, pur_K, mean(f1_grp)))
cat(sprintf("Sensibilidade G4 (raro, n=%d): %.3f\n",
            sum(grupo_int == 4), sens_grp[4]))
cat(strrep("=", 70), "\n", sep = "")
