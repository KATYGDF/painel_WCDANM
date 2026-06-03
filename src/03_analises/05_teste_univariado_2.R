# =============================================================================
# 05_teste_univariado_2.R
#
# TESTE UNIVARIADO 2
# Recuperação de mistura NB sobre a variável `consultas` da simulação do
# projeto, com:
#   - Dados: dados/dados_simulados.csv (4 grupos latentes calibrados a ANS)
#   - Variável: apenas `consultas`
#   - Inicialização: k-means + EM (única estratégia)
#   - theta: LIVRE (estimado por componente)
#   - Critérios: AIC, BIC, ICL, AIC3, Bootstrap LRT
#
# Pergunta: o flexmix com k-means init e theta livre consegue recuperar
# os 4 grupos verdadeiros sobre a variável consultas isoladamente?
# Com qual qualidade de classificação?
#
# Saidas:
#   docs/assets/img/teste_univariado_2/ ............ figuras + tabelas PNG
#   src/03_analises/rds/teste_univariado_2/ ........ RDS de cache + CSV
# =============================================================================


# ── 0. Pacotes ────────────────────────────────────────────────────────────────
pkgs <- c("dplyr", "tidyr", "ggplot2", "flexmix", "gt", "webshot2", "patchwork")
novos <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(novos) > 0) {
  install.packages(novos, repos = "https://cloud.r-project.org", quiet = TRUE)
}
if (!requireNamespace("countreg", quietly = TRUE))
  stop("Instale: install.packages('countreg', repos='http://R-Forge.R-project.org')")

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
  library(flexmix); library(gt); library(patchwork)
})
FLXMRnegbin <- countreg::FLXMRnegbin


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

IMG_DIR <- file.path(PROJ_ROOT, "docs/assets/img/teste_univariado_2")
RDS_DIR <- file.path(PROJ_ROOT, "src/03_analises/rds/teste_univariado_2")
dir.create(IMG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RDS_DIR, recursive = TRUE, showWarnings = FALSE)


# ── 2. Estética ───────────────────────────────────────────────────────────────
COR_GRUPO <- c(
  baixo_uso               = "#4a7c59",
  ambulatorial_coordenado = "#2980b9",
  agudo_hospitalar        = "#e8a33c",
  atipico                 = "#c0392b"
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
CSV_PATH   <- "dados/dados_simulados.csv"
K_GRID     <- 2:6
K_TRUE     <- 4
BLRT_B     <- 99
BLRT_NREP  <- 3
BOOT_B     <- 100
ROTULOS    <- c("baixo_uso", "ambulatorial_coordenado",
                "agudo_hospitalar", "atipico")
G_LABEL    <- c("G1 · baixo uso", "G2 · ambulatorial",
                "G3 · agudo / hospitalar", "G4 · padrão atípico")


# ── 4. Carregar dados ─────────────────────────────────────────────────────────
dados_raw <- read.csv(CSV_PATH, stringsAsFactors = FALSE)
cat(sprintf("Dados carregados: %d obs, %d colunas\n",
            nrow(dados_raw), ncol(dados_raw)))

df <- data.frame(
  grupo_real = factor(dados_raw$grupo_latente, levels = ROTULOS),
  consultas  = as.integer(dados_raw$consultas)
)
N         <- nrow(df)
consultas <- df$consultas
grupo_int <- as.integer(df$grupo_real)

cat(sprintf("Univariado em `consultas`: n = %d\n", N))
cat("Distribuicao por grupo verdadeiro:\n")
print(table(df$grupo_real))


# ── 5. FIG 1 · Distribuição empírica por grupo verdadeiro ───────────────────
df_long <- df %>%
  filter(consultas <= quantile(consultas, 0.97)) %>%
  transmute(grupo = grupo_real, valor = consultas)

fig01 <- ggplot(df_long, aes(x = valor)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 50, fill = "grey80", colour = "white", alpha = 0.7) +
  geom_density(aes(fill = grupo, colour = grupo),
               alpha = 0.25, linewidth = 0.8) +
  scale_fill_manual(values = COR_GRUPO, name = "Grupo verdadeiro",
                    labels = G_LABEL) +
  scale_colour_manual(values = COR_GRUPO, name = "Grupo verdadeiro",
                      labels = G_LABEL) +
  labs(
    x        = "Número de consultas",
    y        = "Densidade",
    title    = "Teste Univariado 2 — consultas por grupo latente",
    subtitle = sprintf("Dados da simulação calibrada (n = %d) · eixo truncado em P97",
                       N),
    caption  = "θ varia entre grupos (calibração ANS); modelo será ajustado com θ livre"
  ) +
  TEMA

salvar_fig(fig01, "fig01_dados_simulados", largura = 8, altura = 4.5)


# =============================================================================
# 6. Funções auxiliares (theta LIVRE)
# =============================================================================

# Ajuste único via flexmix com k-means init + EM, theta livre.
fit_kmeans_em <- function(k_val, seed = 200) {
  set.seed(seed)
  if (k_val == 1L) {
    return(tryCatch(suppressWarnings(flexmix(
      consultas ~ 1, data = df, k = 1L,
      model = FLXMRnegbin(),   # theta LIVRE (NULL implícito)
      control = list(iter.max = 500, tolerance = 1e-7, minprior = 0.01)
    )), error = function(e) NULL))
  }
  km <- kmeans(log1p(df$consultas), centers = k_val,
               nstart = 25, iter.max = 100)
  tryCatch(suppressWarnings(flexmix(
    consultas ~ 1, data = df, k = k_val,
    cluster = km$cluster,
    model   = FLXMRnegbin(),  # theta LIVRE — sem argumento theta
    control = list(iter.max = 500, tolerance = 1e-7, minprior = 0.01)
  )), error = function(e) NULL)
}

# Extrai métricas penalizadas (theta livre afeta `df` automaticamente).
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

# Adjusted Rand Index (base R).
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

# Simulação paramétrica sob H0 com theta LIVRE: lê mu e theta de parameters(fit).
# Em FLXMRnegbin(theta = NULL), parameters() retorna matriz 2xk:
#   linha 1 = log(mu) (intercepto)
#   linha 2 = theta estimado por componente
simular_nb_mistura_livre <- function(fit, n_sim) {
  pesos <- prior(fit)
  pars  <- parameters(fit)
  k     <- fit@k
  pars_m <- if (is.matrix(pars)) pars else matrix(pars, nrow = 2)
  log_mus <- pars_m[1L, ]
  thetas  <- pars_m[2L, ]
  comp <- sample.int(k, size = n_sim, replace = TRUE, prob = pesos)
  y <- integer(n_sim)
  for (h in seq_len(k)) {
    idx <- which(comp == h)
    if (length(idx))
      y[idx] <- rnbinom(length(idx), size = thetas[h], mu = exp(log_mus[h]))
  }
  y
}

# BLRT paramétrico (theta livre) — H0: k vs H1: k+1.
blrt_passo <- function(fit_k0, fit_k1, B = BLRT_B, nrep = BLRT_NREP, seed = 1000) {
  set.seed(seed)
  k0 <- fit_k0@k; k1 <- fit_k1@k
  lrt_obs <- 2 * (as.numeric(logLik(fit_k1)) - as.numeric(logLik(fit_k0)))

  ll_best_k <- function(df_b, k_val) {
    best <- -Inf
    for (i in seq_len(nrep)) {
      m <- tryCatch(suppressWarnings(flexmix(
        consultas ~ 1, data = df_b, k = k_val,
        model = FLXMRnegbin(),  # theta livre
        control = list(iter.max = 200, minprior = 0.01, tolerance = 1e-6)
      )), error = function(e) NULL)
      if (!is.null(m) && m@k == k_val) {
        ll <- as.numeric(logLik(m))
        if (!is.na(ll) && ll > best) best <- ll
      }
    }
    if (is.infinite(best)) NA_real_ else best
  }

  lrt_b <- vapply(seq_len(B), function(b) {
    y_b  <- simular_nb_mistura_livre(fit_k0, N)
    df_b <- data.frame(consultas = y_b)
    ll0  <- ll_best_k(df_b, k0); ll1 <- ll_best_k(df_b, k1)
    if (anyNA(c(ll0, ll1))) NA_real_ else 2 * (ll1 - ll0)
  }, numeric(1))

  lrt_v <- lrt_b[!is.na(lrt_b)]
  pval  <- (sum(lrt_v >= lrt_obs) + 1) / (length(lrt_v) + 1)
  list(k0 = k0, k1 = k1, lrt_obs = round(lrt_obs, 2),
       pval = round(pval, 4), B_valido = length(lrt_v), lrt_boot = lrt_v)
}

# Pareamento húngaro (permutação exaustiva, k <= 8 OK).
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
# 7. MODULO A · Ajuste em K_GRID com k-means + EM (theta livre)
# =============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("MODULO A · k-means + EM com theta livre, k = ", paste(K_GRID, collapse = ", "), "\n", sep = "")
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
    cat(sprintf("  kmeans+EM (theta livre) k=%d ... ", k_val))
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


# ── FIG 2 · Critérios normalizados vs k ──────────────────────────────────────
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

cat(sprintf("\nk* por critério: AIC=%d | BIC=%d | ICL=%d | AIC3=%d\n",
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
  scale_colour_manual(values = c(AIC  = "#27ae60", BIC  = "#e74c3c",
                                  ICL  = "#2980b9", AIC3 = "#8e44ad"),
                      name = "Critério") +
  scale_shape_manual(values = c(AIC  = 15, BIC  = 16, ICL  = 17, AIC3 = 18),
                     name = "Critério") +
  labs(
    title    = "Critérios de seleção normalizados [0 = melhor, 1 = pior]",
    subtitle = sprintf("k* — AIC: %d | BIC: %d | ICL: %d | AIC3: %d  (theta livre)",
                       k_aic, k_bic, k_icl, k_aic3),
    x = "Número de componentes (k)",
    y = "Valor normalizado"
  ) +
  TEMA

salvar_fig(fig02, "fig02_criterios_por_k", largura = 8, altura = 5)


# =============================================================================
# 8. MODULO B · Bootstrap LRT sequencial
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
    if (is.null(fit_k0) || is.null(fit_k1)) {
      cat(sprintf("  k=%d vs k=%d : fit indisponivel — pulando\n", k0, k1))
      next
    }
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

# Fig 3 · BLRT distribuicoes
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
      subtitle = sprintf("B=%d réplicas paramétricas (theta livre)", BLRT_B),
      theme = TEMA
    )
  salvar_fig(fig03, "fig03_blrt_distribuicoes",
             largura = 10,
             altura = 3.5 * ceiling(length(paineis_blrt) / 2))
}


# =============================================================================
# 9. MODULO C · Selecao final, estabilidade, classificacao
# =============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("MODULO C · Selecao final + estabilidade + classificacao\n")
cat(strrep("=", 70), "\n\n", sep = "")

votos <- c(AIC = k_aic, BIC = k_bic, ICL = k_icl, AIC3 = k_aic3,
           BLRT = k_blrt_escolhido)
k_final <- as.integer(names(sort(table(votos), decreasing = TRUE))[1])
cat(sprintf("Votos por k:  AIC=%d | BIC=%d | ICL=%d | AIC3=%d | BLRT=%d\n",
            k_aic, k_bic, k_icl, k_aic3, k_blrt_escolhido))
cat(sprintf("k_final (consenso): %d\n", k_final))

fit_final <- fits[[as.character(k_final)]]

# Bootstrap de estabilidade
ari_cache <- file.path(RDS_DIR, "ari_boot.rds")
if (file.exists(ari_cache)) {
  cat(sprintf("\n[cache] carregando ari_vec (B=%d) de RDS\n", BOOT_B))
  ari_vec <- readRDS(ari_cache)
} else {
  cat(sprintf("\nBootstrap de estabilidade (B=%d) ... ", BOOT_B))
  t0_boot <- proc.time()["elapsed"]
  set.seed(7777)
  cl_ref <- clusters(fit_final)
  ari_vec <- vapply(seq_len(BOOT_B), function(b) {
    idx <- sample.int(N, replace = TRUE)
    df_b <- data.frame(consultas = consultas[idx])
    m_best <- NULL; ll_best <- -Inf
    for (i in seq_len(BLRT_NREP)) {
      m_i <- tryCatch(suppressWarnings(flexmix(
        consultas ~ 1, data = df_b, k = k_final,
        model = FLXMRnegbin(),
        control = list(iter.max = 200, minprior = 0.01, tolerance = 1e-6)
      )), error = function(e) NULL)
      if (!is.null(m_i) && m_i@k == k_final) {
        ll <- as.numeric(logLik(m_i))
        if (!is.na(ll) && ll > ll_best) { ll_best <- ll; m_best <- m_i }
      }
    }
    if (is.null(m_best)) NA_real_
    else adj_rand_index(cl_ref[idx], clusters(m_best))
  }, numeric(1))
  dt_boot <- proc.time()["elapsed"] - t0_boot
  cat(sprintf("%.0fs\n", dt_boot))
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

# Fig 4 · estabilidade
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
# 10. MODULO D · Avaliacao no k = K_TRUE (recuperacao verdadeira)
# =============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat(sprintf("MODULO D · Avaliacao em k = K_TRUE = %d\n", K_TRUE))
cat(strrep("=", 70), "\n\n", sep = "")

fit_K  <- fits[[as.character(K_TRUE)]]
if (is.null(fit_K)) {
  stop("fit em k=K_TRUE indisponivel; nao e possivel avaliar recuperacao verdadeira.")
}

# Métricas globais
est_K     <- as.integer(clusters(fit_K))
est_par   <- remap_hungarian(grupo_int, est_K)
ari_K     <- adj_rand_index(grupo_int, est_K)
acu_K     <- mean(est_par == grupo_int, na.rm = TRUE)
conf_K    <- table(real = grupo_int, est = est_K)
pur_K     <- sum(apply(conf_K, 2L, max)) / sum(conf_K)
f1_grp    <- vapply(seq_len(4), function(g) {
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

# Tabela 3 · métricas globais
tab03 <- data.frame(
  Métrica = c("ARI", "Acurácia", "Pureza", "F1 macro"),
  Valor   = round(c(ari_K, acu_K, pur_K, mean(f1_grp)), 3)
)
write.csv(tab03, file.path(RDS_DIR, "tab03_metricas_globais.csv"),
          row.names = FALSE)

gt03 <- gt(tab03) %>%
  tab_header(
    title    = md(sprintf("**Tabela 3.** Métricas globais no modelo k = %d (theta livre)",
                          K_TRUE)),
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


# Tabela 4 · parâmetros estimados (incluindo theta por componente!)
pars_K   <- parameters(fit_K)
pars_m   <- if (is.matrix(pars_K)) pars_K else matrix(pars_K, nrow = 2)
mus_est  <- exp(pars_m[1L, ])
thetas_est <- pars_m[2L, ]
prior_est  <- as.numeric(prior(fit_K))

# Ordena por mu para alinhar com grupos crescentes
ord_idx <- order(mus_est)
# Hungarian: assign_h[g] = posicao em est_K_levels para grupo real g
est_levels <- sort(unique(est_K))
conf_full <- table(real = grupo_int,
                   est = factor(est_K, levels = est_levels))
assign_h <- hungarian_match(conf_full)

# mu/theta/pi alinhados ao grupo real
mu_aligned    <- vapply(assign_h, function(j) if (is.na(j)) NA_real_ else mus_est[j], numeric(1))
theta_aligned <- vapply(assign_h, function(j) if (is.na(j)) NA_real_ else thetas_est[j], numeric(1))
pi_aligned    <- vapply(assign_h, function(j) if (is.na(j)) NA_real_ else prior_est[j], numeric(1))

tab04 <- data.frame(
  Grupo       = G_LABEL,
  n_real      = as.integer(table(grupo_int)),
  mu_est      = round(mu_aligned, 2),
  theta_est   = round(theta_aligned, 2),
  pi_est      = round(pi_aligned, 3),
  sensibilidade = round(sens_grp, 3),
  F1          = round(f1_grp, 3)
)
write.csv(tab04, file.path(RDS_DIR, "tab04_parametros_diagnostico.csv"),
          row.names = FALSE)

gt04 <- gt(tab04) %>%
  tab_header(
    title    = md(sprintf("**Tabela 4.** Parâmetros estimados + diagnóstico por grupo (k = %d)",
                          K_TRUE)),
    subtitle = md("μ, θ e π estimados pelo flexmix com θ livre. Alinhamento via algoritmo húngaro.")
  ) %>%
  cols_label(
    Grupo = "Grupo verdadeiro",
    n_real = md("*n* real"),
    mu_est = md("μ estimado"),
    theta_est = md("θ estimado"),
    pi_est = md("π estimado"),
    sensibilidade = "Sensibilidade",
    F1 = "F1"
  ) %>%
  fmt_number(columns = c(mu_est, theta_est), decimals = 2) %>%
  fmt_number(columns = c(pi_est, sensibilidade, F1), decimals = 3) %>%
  tab_style(style = cell_text(color = "white", weight = "bold"),
            locations = cells_column_labels(everything())) %>%
  tab_options(table.font.size = px(11),
              column_labels.background.color = "#2c3e50",
              column_labels.border.top.style = "hidden",
              heading.align = "left")
salvar_gt(gt04, "tab04_parametros_diagnostico", w_in = 10, h_in = 3.5)


# ── Fig 5 · Densidade ajustada vs dados (theta livre por componente) ─────────
mus_sorted    <- mus_est[ord_idx]
thetas_sorted <- thetas_est[ord_idx]
prior_sorted  <- prior_est[ord_idx]
grade_x <- seq(0, quantile(consultas, 0.97), length.out = 400)
dens_comp <- do.call(rbind, lapply(seq_along(mus_sorted), function(h) {
  data.frame(
    x = grade_x,
    densidade = prior_sorted[h] * dnbinom(round(grade_x),
                                           size = thetas_sorted[h],
                                           mu   = mus_sorted[h]),
    componente = sprintf("C%d (μ=%.1f, θ=%.2f)",
                         h, mus_sorted[h], thetas_sorted[h])
  )
}))
dens_total <- data.frame(
  x = grade_x,
  densidade = sapply(grade_x, function(xv) {
    sum(prior_sorted * dnbinom(round(xv),
                                size = thetas_sorted,
                                mu   = mus_sorted))
  })
)

df_hist <- df %>% filter(consultas <= max(grade_x))

fig05 <- ggplot() +
  geom_histogram(data = df_hist, aes(x = consultas, y = after_stat(density)),
                 bins = 50, fill = "grey80", colour = "white", alpha = 0.7) +
  geom_line(data = dens_comp,
            aes(x = x, y = densidade, colour = componente),
            linewidth = 0.9) +
  geom_line(data = dens_total,
            aes(x = x, y = densidade),
            colour = "black", linewidth = 1, linetype = "dashed") +
  scale_colour_brewer(palette = "Set1") +
  labs(
    title    = sprintf("Modelo ajustado: k = %d, k-means + EM, θ livre", K_TRUE),
    subtitle = sprintf("ARI = %.3f · Acurácia = %.3f · F1 macro = %.3f",
                       ari_K, acu_K, mean(f1_grp)),
    x = "Número de consultas", y = "Densidade",
    colour = "Componente estimado"
  ) + TEMA

salvar_fig(fig05, "fig05_modelo_ajustado", largura = 10, altura = 5)


# ── Fig 6 · Matriz de confusão ───────────────────────────────────────────────
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
    title = sprintf("Matriz de confusão — k = %d (k-means + EM, θ livre)", K_TRUE),
    subtitle = sprintf("Pareamento via algoritmo húngaro · Acurácia = %.1f%%",
                       100 * acu_K)
  ) +
  TEMA + theme(panel.grid = element_blank(),
               axis.text.x = element_text(face = "bold"),
               axis.text.y = element_text(face = "bold"))

salvar_fig(fig06, "fig06_matriz_confusao", largura = 9, altura = 5.5)


# ── Fig 7 · Strip plot 4x2 — classificacao por grupo ─────────────────────────
LABEL_PANEL <- c(verdadeiro = "verdadeiro", kmeans = "k-means + EM (θ livre)")
mapped_assign <- list(verdadeiro = grupo_int, kmeans = est_par)

grid_data <- do.call(rbind, lapply(seq_len(4), function(g_focus) {
  do.call(rbind, lapply(names(mapped_assign), function(meth) {
    data.frame(
      g_focus_lbl = G_LABEL[g_focus],
      method_lbl  = LABEL_PANEL[meth],
      consultas   = consultas,
      in_focus    = !is.na(mapped_assign[[meth]]) &
                    mapped_assign[[meth]] == g_focus
    )
  }))
}))
grid_data <- grid_data[grid_data$consultas <= quantile(consultas, 0.97), ]
grid_data$g_focus_lbl <- factor(grid_data$g_focus_lbl, levels = G_LABEL)
grid_data$method_lbl  <- factor(grid_data$method_lbl, levels = LABEL_PANEL)

COR_FOCUS <- setNames(unname(COR_GRUPO[ROTULOS]), G_LABEL)

fig07 <- ggplot(grid_data, aes(x = consultas, y = 0)) +
  geom_jitter(data = subset(grid_data, !in_focus),
              colour = "grey85", alpha = 0.20, size = 0.35, height = 0.4) +
  geom_jitter(data = subset(grid_data, in_focus),
              aes(colour = g_focus_lbl),
              alpha = 0.65, size = 0.5, height = 0.4) +
  facet_grid(g_focus_lbl ~ method_lbl, switch = "y") +
  scale_colour_manual(values = COR_FOCUS, guide = "none") +
  scale_y_continuous(limits = c(-1, 1), breaks = NULL) +
  labs(
    x = "Número de consultas", y = NULL,
    title = sprintf("Classificação por grupo — k = %d (univariado)", K_TRUE),
    subtitle = "Cor viva = grupo em foco · cinza = demais"
  ) +
  TEMA + theme(legend.position = "none",
               strip.text.y.left = element_text(angle = 0, hjust = 1,
                                                  face = "bold"),
               strip.text.x = element_text(face = "bold"))

salvar_fig(fig07, "fig07_classificacao_grid", largura = 10, altura = 7)


# ── Fig 8 · Taxa de erro por faixa do número de consultas ───────────────────
faixas    <- c(0, 5, 15, 35, Inf)
faixa_lbl <- c("0–4", "5–14", "15–34", "≥ 35")
faixa_i   <- cut(consultas, breaks = faixas, right = FALSE, labels = faixa_lbl)

erro_data <- do.call(rbind, lapply(faixa_lbl, function(fx) {
  idx <- faixa_i == fx & !is.na(faixa_i)
  if (sum(idx) == 0) return(NULL)
  data.frame(
    faixa = fx,
    n = sum(idx),
    err = mean(est_par[idx] != grupo_int[idx], na.rm = TRUE)
  )
}))
erro_data$faixa <- factor(erro_data$faixa, levels = faixa_lbl)

fig08 <- ggplot(erro_data, aes(x = faixa, y = err)) +
  geom_col(width = 0.6, fill = "#16a085") +
  geom_text(aes(label = sprintf("%.1f%% (n=%d)", 100 * err, n)),
            vjust = -0.3, size = 3) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.12))) +
  labs(
    x = "Faixa do número de consultas",
    y = "Taxa de erro",
    title = sprintf("Erro de classificação por faixa (k = %d)", K_TRUE),
    subtitle = "k-means + EM com θ livre"
  ) + TEMA

salvar_fig(fig08, "fig08_erro_por_faixa", largura = 9, altura = 5)


# =============================================================================
# 11. RESUMO
# =============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("RESUMO — Teste Univariado 2\n")
cat(strrep("=", 70), "\n", sep = "")
cat(sprintf("Saidas em: %s\n", IMG_DIR))
cat("Figuras:\n")
cat("  fig01_dados_simulados.png       (distribuicoes verdadeiras)\n")
cat("  fig02_criterios_por_k.png       (AIC/BIC/ICL/AIC3 vs k)\n")
cat("  fig03_blrt_distribuicoes.png    (Bootstrap LRT)\n")
cat("  fig04_estabilidade_ari.png      (estabilidade ARI)\n")
cat("  fig05_modelo_ajustado.png       (densidades ajustadas)\n")
cat("  fig06_matriz_confusao.png       (real x estimado em k=4)\n")
cat("  fig07_classificacao_grid.png    (strip 4x2 ground truth vs kmeans)\n")
cat("  fig08_erro_por_faixa.png        (taxa de erro por faixa de consultas)\n")
cat("Tabelas (PNG + CSV):\n")
cat("  tab01_criterios                 (k vs AIC/BIC/ICL/AIC3/tempo)\n")
cat("  tab02_blrt                      (BLRT sequencial)\n")
cat("  tab03_metricas_globais          (ARI/Acuracia/Pureza/F1 em k=4)\n")
cat("  tab04_parametros_diagnostico    (mu/theta/pi/sens/F1 por grupo)\n")
cat(sprintf("\nVotos por k: AIC=%d | BIC=%d | ICL=%d | AIC3=%d | BLRT=%d\n",
            k_aic, k_bic, k_icl, k_aic3, k_blrt_escolhido))
cat(sprintf("k_final (consenso) = %d\n", k_final))
cat(sprintf("Em k=%d: ARI=%.3f | Acuracia=%.3f | Pureza=%.3f | F1 macro=%.3f\n",
            K_TRUE, ari_K, acu_K, pur_K, mean(f1_grp)))
cat(strrep("=", 70), "\n", sep = "")
