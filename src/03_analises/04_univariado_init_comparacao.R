# =============================================================================
# 04_univariado_init_comparacao.R
#
# Comparacao de estrategias de inicializacao do EM em mistura NB univariada
# (variavel: consultas) com avaliacao por AIC, BIC, ICL e Bootstrap LRT.
#
# Objetivo: identificar a melhor estrategia para separar grupos sobrepostos
# (cenario com theta = 2.5, mu = 2, 12, 25, 50; G3 e G4 muito sobrepostos).
#
# Referencias metodologicas:
#   - Scharl, Grun & Leisch (2010), Bioinformatics 26(3):370-377
#   - Biernacki, Celeux & Govaert (2003), CSDA 41:561-575
#   - Celeux & Govaert (1992), CSDA 14:315-332     [CEM]
#   - Celeux & Diebolt (1985), Comp. Stat. Quart.   [SEM]
#   - McLachlan (1987) bootstrap LRT para misturas
#
# Saidas:
#   docs/assets/img/univariado_init/ ......... figuras PNG
#   src/03_analises/rds/univariado_init/ ..... resultados + tabelas CSV
# =============================================================================


# ── 0. Pacotes ────────────────────────────────────────────────────────────────
pkgs <- c("dplyr", "tidyr", "ggplot2", "flexmix", "gt", "webshot2", "patchwork")
novos <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(novos) > 0) {
  message("Instalando: ", paste(novos, collapse = ", "))
  install.packages(novos, repos = "https://cloud.r-project.org", quiet = TRUE)
}
if (!requireNamespace("countreg", quietly = TRUE))
  stop("Instale: install.packages('countreg', repos='http://R-Forge.R-project.org')")

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
  library(flexmix); library(gt); library(patchwork)
})
FLXMRnegbin <- countreg::FLXMRnegbin


# ── 1. Raiz do projeto + diretorios ──────────────────────────────────────────
# Detecta a raiz do projeto procurando o arquivo-marcador (poster_figures.R).
# Garante que IMG_DIR e RDS_DIR caiam sempre na pasta correta,
# independentemente do diretorio de trabalho onde o script foi invocado.
encontrar_raiz <- function(marcador = "poster_figures.R",
                            inicio = getwd(), max_subir = 6L) {
  dir_atual <- normalizePath(inicio, winslash = "/", mustWork = FALSE)
  for (i in seq_len(max_subir)) {
    if (file.exists(file.path(dir_atual, marcador))) return(dir_atual)
    pai <- dirname(dir_atual)
    if (pai == dir_atual) break  # chegou na raiz do sistema de arquivos
    dir_atual <- pai
  }
  stop("Nao encontrei a raiz do projeto (procurando por '", marcador,
       "'). Rode setwd() para a raiz e tente novamente.")
}

PROJ_ROOT <- encontrar_raiz()
setwd(PROJ_ROOT)
cat("Raiz do projeto: ", PROJ_ROOT, "\n", sep = "")

IMG_DIR <- file.path(PROJ_ROOT, "docs/assets/img/univariado_init")
RDS_DIR <- file.path(PROJ_ROOT, "src/03_analises/rds/univariado_init")
dir.create(IMG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RDS_DIR, recursive = TRUE, showWarnings = FALSE)


# ── 2. Estetica ───────────────────────────────────────────────────────────────
COR_GRUPO <- c(
  baixo_uso               = "#4a7c59",   # G1 verde
  ambulatorial_coordenado = "#2980b9",   # G2 azul
  agudo_hospitalar        = "#e8a33c",   # G3 laranja/dourado
  atipico                 = "#c0392b"    # G4 vermelho
)

COR_INIT <- c(
  random  = "#7f8c8d",   # cinza
  kmeans  = "#16a085",   # verde-azulado
  cem_em  = "#8e44ad",   # roxo
  sem_em  = "#d35400"    # laranja
)

LABEL_INIT <- c(
  random  = "Random multistart",
  kmeans  = "k-means + EM",
  cem_em  = "CEM + EM",
  sem_em  = "SEM + EM"
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


# ── 3. Parametros da experiencia ──────────────────────────────────────────────
set.seed(1234)
N         <- 5000
PI_TRUE   <- c(0.25, 0.25, 0.25, 0.25)
MU_TRUE   <- c(2, 12, 25, 50)
THETA_NB  <- 2.5
K_TRUE    <- 4
K_GRID    <- 2:6                                     # k para varrer
NREP      <- 20                                      # reinicios por estrategia
BLRT_B    <- 99                                      # bootstrap BLRT
BLRT_NREP <- 3                                       # reinicios por boot
BOOT_B    <- 100                                     # estabilidade
ROTULOS   <- c("baixo_uso", "ambulatorial_coordenado",
               "agudo_hospitalar", "atipico")


# ── 4. Simulacao dos dados ────────────────────────────────────────────────────
grupo_real <- sample(1:4, size = N, replace = TRUE, prob = PI_TRUE)
consultas  <- rnbinom(n = N, size = THETA_NB, mu = MU_TRUE[grupo_real])

df <- data.frame(
  grupo_real = factor(grupo_real, levels = 1:4, labels = ROTULOS),
  consultas  = consultas
)

cat(sprintf("Dados simulados: n=%d, k_true=%d, theta=%.1f\n",
            N, K_TRUE, THETA_NB))
cat(sprintf("  mu verdadeiros : %s\n", paste(MU_TRUE, collapse = ", ")))
cat(sprintf("  pi verdadeiros : %s\n", paste(PI_TRUE, collapse = ", ")))


# ── 5. FIG 1 · Distribuicao dos grupos verdadeiros ───────────────────────────
df_long <- df %>%
  filter(consultas <= quantile(consultas, 0.97)) %>%
  transmute(grupo = grupo_real, valor = consultas)

fig01 <- ggplot(df_long, aes(x = valor)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 50, fill = "grey80", colour = "white", alpha = 0.7) +
  geom_density(aes(fill = grupo, colour = grupo),
               alpha = 0.25, linewidth = 0.8) +
  scale_fill_manual(values = COR_GRUPO, name = "Grupo verdadeiro") +
  scale_colour_manual(values = COR_GRUPO, name = "Grupo verdadeiro") +
  labs(
    x        = "Numero de consultas",
    y        = "Densidade",
    title    = "Distribuicao de consultas por grupo latente (gabarito)",
    subtitle = sprintf("n = %d  ·  theta = %.1f  ·  mu = %s  ·  eixo truncado em P97",
                       N, THETA_NB, paste(MU_TRUE, collapse = ", ")),
    caption  = "G3 (agudo) e G4 (atipico) muito sobrepostos sob theta baixo"
  ) +
  TEMA

salvar_fig(fig01, "fig01_dados_simulados", largura = 8, altura = 4.5)


# =============================================================================
# 6. Funcoes auxiliares
# =============================================================================

# Ajuste para um k, dado um cluster inicial (NULL = aleatorio).
fit_um <- function(cluster_init, k_val, classify_mode = "weighted", seed = 42,
                   iter_max = 500, tolerance = 1e-7, minprior = 0.01) {
  set.seed(seed)
  ctrl <- list(iter.max = iter_max, tolerance = tolerance,
               minprior = minprior, classify = classify_mode)
  args <- list(
    formula = consultas ~ 1, data = df, k = k_val,
    model = FLXMRnegbin(theta = THETA_NB), control = ctrl
  )
  if (!is.null(cluster_init)) args$cluster <- cluster_init
  tryCatch(suppressWarnings(do.call(flexmix, args)),
           error = function(e) NULL)
}

# Random multistart: nrep ajustes EM com posteriors aleatorios.
init_random <- function(k_val, nrep = NREP, seed = 100) {
  best <- NULL; best_ll <- -Inf
  for (i in seq_len(nrep)) {
    fit_i <- fit_um(cluster_init = sample.int(k_val, N, replace = TRUE),
                    k_val = k_val, classify_mode = "weighted",
                    seed = seed + i)
    if (!is.null(fit_i) && fit_i@k == k_val) {
      ll <- as.numeric(logLik(fit_i))
      if (!is.na(ll) && ll > best_ll) { best_ll <- ll; best <- fit_i }
    }
  }
  best
}

# Inicializacao por k-means em log(1+x).
init_kmeans <- function(k_val, seed = 200) {
  set.seed(seed)
  if (k_val == 1L) return(fit_um(NULL, 1L, seed = seed))
  km <- kmeans(log1p(df$consultas), centers = k_val, nstart = 25, iter.max = 100)
  fit_um(cluster_init = km$cluster, k_val = k_val, seed = seed)
}

# CEM + EM: nrep arranques aleatorios com classify="hard" (CEM),
# pega o melhor por logLik e refina com EM completo (classify="weighted").
init_cem_em <- function(k_val, nrep = NREP, seed = 300) {
  best_cem <- NULL; best_ll <- -Inf
  for (i in seq_len(nrep)) {
    fit_i <- fit_um(cluster_init = sample.int(k_val, N, replace = TRUE),
                    k_val = k_val, classify_mode = "hard",
                    seed = seed + i, iter_max = 100, tolerance = 1e-5)
    if (!is.null(fit_i) && fit_i@k == k_val) {
      ll <- as.numeric(logLik(fit_i))
      if (!is.na(ll) && ll > best_ll) { best_ll <- ll; best_cem <- fit_i }
    }
  }
  if (is.null(best_cem)) return(NULL)
  fit_um(cluster_init = clusters(best_cem), k_val = k_val,
         classify_mode = "weighted", seed = seed)
}

# SEM + EM: identico ao CEM-EM mas com classify="random" (estocastico).
init_sem_em <- function(k_val, nrep = NREP, seed = 400) {
  best_sem <- NULL; best_ll <- -Inf
  for (i in seq_len(nrep)) {
    fit_i <- fit_um(cluster_init = sample.int(k_val, N, replace = TRUE),
                    k_val = k_val, classify_mode = "random",
                    seed = seed + i, iter_max = 100, tolerance = 1e-5)
    if (!is.null(fit_i) && fit_i@k == k_val) {
      ll <- as.numeric(logLik(fit_i))
      if (!is.na(ll) && ll > best_ll) { best_ll <- ll; best_sem <- fit_i }
    }
  }
  if (is.null(best_sem)) return(NULL)
  fit_um(cluster_init = clusters(best_sem), k_val = k_val,
         classify_mode = "weighted", seed = seed)
}

# Extrai metricas de um fit flexmix.
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
  a <- sum(choose(tab, 2))
  rs <- rowSums(tab); cs <- colSums(tab); n <- sum(tab)
  b <- sum(choose(rs, 2)) - a
  cv <- sum(choose(cs, 2)) - a
  exp_a <- (a + b) * (a + cv) / choose(n, 2)
  den <- 0.5 * ((a + b) + (a + cv)) - exp_a
  if (den <= 0) return(1)
  (a - exp_a) / den
}

# Simula bootstrap parametrico sob H0 (k0 componentes).
simular_nb_mistura <- function(fit, n_sim, theta_fixo = THETA_NB) {
  pesos <- prior(fit); pars <- as.numeric(parameters(fit)); k <- fit@k
  log_mus <- pars; thetas <- rep(theta_fixo, k)
  comp <- sample.int(k, size = n_sim, replace = TRUE, prob = pesos)
  y <- integer(n_sim)
  for (h in seq_len(k)) {
    idx <- which(comp == h)
    if (length(idx)) y[idx] <- rnbinom(length(idx), size = thetas[h],
                                        mu = exp(log_mus[h]))
  }
  y
}

# BLRT: H0: k vs H1: k+1, paramétrico.
blrt_passo <- function(fit_k0, fit_k1, B = BLRT_B, nrep = BLRT_NREP, seed = 1000) {
  set.seed(seed)
  k0 <- fit_k0@k; k1 <- fit_k1@k
  lrt_obs <- 2 * (as.numeric(logLik(fit_k1)) - as.numeric(logLik(fit_k0)))

  ll_best_k <- function(df_b, k_val) {
    best <- -Inf
    for (i in seq_len(nrep)) {
      m <- tryCatch(suppressWarnings(flexmix(
        consultas ~ 1, data = df_b, k = k_val,
        model = FLXMRnegbin(theta = THETA_NB),
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
    y_b <- simular_nb_mistura(fit_k0, N)
    df_b <- data.frame(consultas = y_b)
    ll0 <- ll_best_k(df_b, k0); ll1 <- ll_best_k(df_b, k1)
    if (anyNA(c(ll0, ll1))) NA_real_ else 2 * (ll1 - ll0)
  }, numeric(1))

  lrt_v <- lrt_b[!is.na(lrt_b)]
  pval <- (sum(lrt_v >= lrt_obs) + 1) / (length(lrt_v) + 1)
  list(k0 = k0, k1 = k1, lrt_obs = round(lrt_obs, 2),
       pval = round(pval, 4), B_valido = length(lrt_v), lrt_boot = lrt_v)
}

# Pareamento hungaro via permutacao exaustiva (k <= 8 OK).
perms_all <- function(n) {
  if (n == 1L) return(matrix(1L, 1L, 1L))
  p <- perms_all(n - 1L)
  do.call(rbind, lapply(seq_len(n), function(i) {
    cbind(i, p + (p >= i))
  }))
}
hungarian_match <- function(conf_mat) {
  nr <- nrow(conf_mat); nc <- ncol(conf_mat); k <- max(nr, nc)
  sq <- matrix(0L, k, k); sq[seq_len(nr), seq_len(nc)] <- conf_mat
  perms <- perms_all(k)
  scores <- apply(perms, 1L, function(p) sum(diag(sq[, p, drop = FALSE])))
  a <- as.integer(perms[which.max(scores), ])
  ifelse(a[seq_len(nr)] <= nc, a[seq_len(nr)], NA_integer_)
}


# =============================================================================
# 7. MODULO A · Comparacao de inicializacoes × k
# =============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("MODULO A · Comparacao de inicializacoes\n")
cat(strrep("=", 70), "\n\n", sep = "")

INITS <- list(
  random = init_random,
  kmeans = init_kmeans,
  cem_em = init_cem_em,
  sem_em = init_sem_em
)

# Cache: se fits.rds + tab_init.rds existem, pula o loop pesado.
fits_cache_path <- file.path(RDS_DIR, "fits.rds")
tab_init_cache  <- file.path(RDS_DIR, "tab_init.rds")
USE_CACHE_A <- file.exists(fits_cache_path) && file.exists(tab_init_cache)

if (USE_CACHE_A) {
  cat("  [cache] carregando fits e tab_init de RDS — pulando ajuste\n")
  fits     <- readRDS(fits_cache_path)
  tab_init <- readRDS(tab_init_cache)
} else {
  fits     <- setNames(vector("list", length(INITS)), names(INITS))
  tab_init <- list()

  for (init_nm in names(INITS)) {
    fits[[init_nm]] <- setNames(vector("list", length(K_GRID)),
                                as.character(K_GRID))
    for (k_val in K_GRID) {
      cat(sprintf("  %-8s k=%d ... ", init_nm, k_val))
      t0 <- proc.time()["elapsed"]
      fit_kv <- INITS[[init_nm]](k_val)
      dt <- proc.time()["elapsed"] - t0
      fits[[init_nm]][[as.character(k_val)]] <- fit_kv
      m <- extrair_metricas(fit_kv)
      if (is.null(m)) {
        cat(sprintf("FALHOU (%.1fs)\n", dt))
        tab_init[[length(tab_init) + 1]] <- data.frame(
          init = init_nm, k = k_val, k_efetivo = NA, logLik = NA, np = NA,
          AIC = NA, BIC = NA, ICL = NA, AIC3 = NA, tempo_s = round(dt, 1)
        )
      } else {
        cat(sprintf("logLik=%.1f | BIC=%.1f | %.1fs\n",
                    m$logLik, m$BIC, dt))
        tab_init[[length(tab_init) + 1]] <- data.frame(
          init = init_nm, k = k_val,
          k_efetivo = m$k_efetivo, logLik = m$logLik, np = m$np,
          AIC = m$AIC, BIC = m$BIC, ICL = m$ICL, AIC3 = m$AIC3,
          tempo_s = round(dt, 1)
        )
      }
    }
  }
  tab_init <- do.call(rbind, tab_init)
  saveRDS(fits, fits_cache_path)
  saveRDS(tab_init, tab_init_cache)
}

write.csv(tab_init, file.path(RDS_DIR, "tab01_inicializacoes.csv"),
          row.names = FALSE)


# ── FIG 2 · logLik vs k por inicializacao ─────────────────────────────────────
fig02 <- ggplot(tab_init, aes(x = k, y = logLik, colour = init, shape = init)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_x_continuous(breaks = K_GRID) +
  scale_colour_manual(values = COR_INIT, labels = LABEL_INIT, name = "Inicializacao") +
  scale_shape_manual(values = c(random = 16, kmeans = 17, cem_em = 15, sem_em = 18),
                     labels = LABEL_INIT, name = "Inicializacao") +
  labs(
    title    = "log-verossimilhanca por estrategia de inicializacao",
    subtitle = sprintf("Mistura NB univariada (theta=%.1f, n=%d) — %d reinicios por estrategia",
                       THETA_NB, N, NREP),
    x = "Numero de componentes (k)",
    y = "log-verossimilhanca (maior = melhor)"
  ) +
  TEMA

salvar_fig(fig02, "fig02_loglik_por_init", largura = 8, altura = 5)


# =============================================================================
# 8. MODULO B · Selecao do k por criterios (na melhor inicializacao)
# =============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("MODULO B · Selecao de k por AIC, BIC, ICL, AIC3\n")
cat(strrep("=", 70), "\n\n", sep = "")

# Identifica a melhor inicializacao por k = melhor logLik
melhor_por_k <- tab_init %>%
  filter(!is.na(logLik)) %>%
  group_by(k) %>%
  slice_max(logLik, n = 1, with_ties = FALSE) %>%
  ungroup()

cat("Melhor inicializacao por k (criterio: maior logLik):\n")
print(melhor_por_k[, c("k", "init", "logLik", "BIC", "ICL")], row.names = FALSE)

# k* por criterio
k_aic  <- melhor_por_k$k[which.min(melhor_por_k$AIC)]
k_bic  <- melhor_por_k$k[which.min(melhor_por_k$BIC)]
k_icl  <- melhor_por_k$k[which.min(melhor_por_k$ICL)]
k_aic3 <- melhor_por_k$k[which.min(melhor_por_k$AIC3)]

cat(sprintf("\nk* por criterio:  AIC=%d | BIC=%d | ICL=%d | AIC3=%d\n",
            k_aic, k_bic, k_icl, k_aic3))

# Normaliza criterios para um grafico comparativo
norm01 <- function(v) {
  r <- range(v, na.rm = TRUE)
  if (diff(r) == 0) return(rep(0, length(v)))
  (v - r[1]) / diff(r)
}

tab_crit_long <- data.frame(
  k    = rep(melhor_por_k$k, 4),
  crit = rep(c("AIC", "BIC", "ICL", "AIC3"), each = nrow(melhor_por_k)),
  vn   = c(norm01(melhor_por_k$AIC),
            norm01(melhor_por_k$BIC),
            norm01(melhor_por_k$ICL),
            norm01(melhor_por_k$AIC3))
)

fig03 <- ggplot(tab_crit_long, aes(x = k, y = vn, colour = crit, shape = crit)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 3) +
  geom_vline(xintercept = K_TRUE, linetype = "dashed",
             colour = "grey55", linewidth = 0.5) +
  annotate("text", x = K_TRUE, y = 1.05, label = "k verdadeiro = 4",
           colour = "grey45", size = 3, vjust = 0) +
  scale_x_continuous(breaks = K_GRID) +
  scale_colour_manual(values = c(AIC  = "#27ae60", BIC  = "#e74c3c",
                                  ICL  = "#2980b9", AIC3 = "#8e44ad"),
                      name = "Criterio") +
  scale_shape_manual(values = c(AIC  = 15, BIC  = 16, ICL  = 17, AIC3 = 18),
                     name = "Criterio") +
  labs(
    title    = "Criterios de selecao normalizados [0 = melhor, 1 = pior]",
    subtitle = sprintf("k* — AIC: %d | BIC: %d | ICL: %d | AIC3: %d",
                       k_aic, k_bic, k_icl, k_aic3),
    x = "Numero de componentes (k)",
    y = "Valor normalizado"
  ) +
  TEMA

salvar_fig(fig03, "fig03_criterios_por_k", largura = 8, altura = 5)


# =============================================================================
# 9. MODULO C · Bootstrap LRT sequencial
# =============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("MODULO C · Bootstrap LRT (H0: k vs H1: k+1)\n")
cat(strrep("=", 70), "\n\n", sep = "")

# Usa fits da melhor inicializacao por k
melhor_fits <- setNames(vector("list", nrow(melhor_por_k)),
                        as.character(melhor_por_k$k))
for (i in seq_len(nrow(melhor_por_k))) {
  k_i <- melhor_por_k$k[i]; init_i <- melhor_por_k$init[i]
  melhor_fits[[as.character(k_i)]] <- fits[[init_i]][[as.character(k_i)]]
}

# BLRT testa k = 2, 3, 4, 5 contra k+1
k_blrt_range <- intersect(K_GRID[-length(K_GRID)],
                          as.integer(names(melhor_fits)))

blrt_cache_path <- file.path(RDS_DIR, "resultados_blrt.rds")

if (file.exists(blrt_cache_path)) {
  cat("  [cache] carregando resultados_blrt de RDS — pulando ", BLRT_B,
      " replicacoes bootstrap por par\n", sep = "")
  resultados_blrt  <- readRDS(blrt_cache_path)
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
    fit_k0 <- melhor_fits[[as.character(k0)]]
    fit_k1 <- melhor_fits[[as.character(k1)]]
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
  saveRDS(resultados_blrt, blrt_cache_path)
}

if (is.na(k_blrt_escolhido))
  k_blrt_escolhido <- max(k_blrt_range) + 1L

cat(sprintf("\nk* por BLRT: %d\n", k_blrt_escolhido))

# Tabela BLRT
tab_blrt <- do.call(rbind, lapply(resultados_blrt, function(r) {
  data.frame(k0 = r$k0, k1 = r$k1,
             LRT_obs = r$lrt_obs, p_valor = r$pval,
             B_validos = r$B_valido,
             rejeita_H0 = ifelse(r$pval <= 0.05, "Sim", "Nao"))
}))
write.csv(tab_blrt, file.path(RDS_DIR, "tab03_blrt.csv"), row.names = FALSE)


# ── FIG 4 · BLRT — distribuicoes bootstrap (multi-painel) ─────────────────────
plot_blrt_um <- function(r) {
  dfh <- data.frame(lrt = r$lrt_boot)
  ggplot(dfh, aes(x = lrt)) +
    geom_histogram(bins = 25, fill = "grey70", colour = "white") +
    geom_vline(xintercept = r$lrt_obs, colour = "#e74c3c",
               linewidth = 1, linetype = "dashed") +
    annotate("text", x = r$lrt_obs,
             y = Inf, vjust = 2, hjust = -0.05,
             label = sprintf("LRT_obs\np = %.4f", r$pval),
             colour = "#c0392b", size = 3) +
    labs(
      title = sprintf("k=%d vs k=%d", r$k0, r$k1),
      x = "LRT bootstrap", y = "Frequencia"
    ) +
    TEMA + theme(plot.title = element_text(size = 11))
}

paineis_blrt <- lapply(resultados_blrt, plot_blrt_um)
if (length(paineis_blrt) >= 1) {
  fig04 <- wrap_plots(paineis_blrt, ncol = 2) +
    plot_annotation(
      title = "Bootstrap LRT — distribuicao sob H0 para cada par (k, k+1)",
      subtitle = sprintf("B=%d replicacoes paramétricas | linha vermelha = LRT observado",
                         BLRT_B),
      theme = TEMA
    )
  salvar_fig(fig04, "fig04_blrt_distribuicoes",
             largura = 10, altura = 3.5 * ceiling(length(paineis_blrt) / 2))
}


# =============================================================================
# 10. MODULO D · Selecao final, estabilidade e classificacao
# =============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("MODULO D · Selecao final + estabilidade + classificacao\n")
cat(strrep("=", 70), "\n\n", sep = "")

# Consenso entre criterios
votos <- c(AIC = k_aic, BIC = k_bic, ICL = k_icl, AIC3 = k_aic3,
           BLRT = k_blrt_escolhido)
k_final <- as.integer(names(sort(table(votos), decreasing = TRUE))[1])
cat(sprintf("Votos por k:  AIC=%d | BIC=%d | ICL=%d | AIC3=%d | BLRT=%d\n",
            k_aic, k_bic, k_icl, k_aic3, k_blrt_escolhido))
cat(sprintf("k_final (consenso): %d\n", k_final))

fit_final <- melhor_fits[[as.character(k_final)]]
init_final <- melhor_por_k$init[melhor_por_k$k == k_final]
cat(sprintf("Inicializacao do modelo final: %s\n", init_final))

# Bootstrap de estabilidade (nao-paramétrico) — com cache
ari_boot_cache <- file.path(RDS_DIR, "ari_boot.rds")

if (file.exists(ari_boot_cache)) {
  cat(sprintf("\n[cache] carregando ari_vec (B=%d) de RDS\n", BOOT_B))
  ari_vec <- readRDS(ari_boot_cache)
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
        model = FLXMRnegbin(theta = THETA_NB),
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
  saveRDS(ari_vec, ari_boot_cache)
}

ari_valid <- ari_vec[!is.na(ari_vec)]
ari_med <- if (length(ari_valid)) mean(ari_valid) else NA
ari_ic  <- if (length(ari_valid)) quantile(ari_valid, c(0.025, 0.975)) else c(NA, NA)
cat(sprintf("ARI medio = %.3f | IC95%% = [%.3f, %.3f] | B_validos = %d\n",
            ari_med, ari_ic[1], ari_ic[2], length(ari_valid)))


# ── FIG 5 · Bootstrap de estabilidade ──────────────────────────────────────
if (length(ari_valid)) {
  fig05 <- ggplot(data.frame(ari = ari_valid), aes(x = ari)) +
    geom_histogram(bins = 25, fill = "#3498db", colour = "white", alpha = 0.85) +
    geom_vline(xintercept = ari_med, colour = "#e74c3c",
               linewidth = 1, linetype = "dashed") +
    geom_vline(xintercept = 0.80, colour = "#f39c12",
               linewidth = 0.7, linetype = "dotted") +
    annotate("text", x = 0.80, y = Inf, vjust = 1.5, hjust = 1.1,
             label = "limiar 0.80", colour = "#d68910", size = 3) +
    labs(
      title = sprintf("Estabilidade do modelo k=%d (bootstrap nao-paramétrico)",
                      k_final),
      subtitle = sprintf("ARI medio=%.3f | IC95%% [%.3f, %.3f] | B=%d",
                         ari_med, ari_ic[1], ari_ic[2], length(ari_valid)),
      x = "Adjusted Rand Index (bootstrap vs original)",
      y = "Frequencia"
    ) + TEMA
  salvar_fig(fig05, "fig05_estabilidade_ari", largura = 8, altura = 4.5)
}


# ── Avaliacao de classificacao do modelo final ───────────────────────────────
cl_final <- clusters(fit_final)
conf <- table(real = df$grupo_real,
              est  = factor(cl_final, levels = seq_len(k_final)))

if (k_final >= 4 && nrow(conf) == ncol(conf)) {
  assign_h <- hungarian_match(conf)
  conf_match <- conf[, assign_h]
  acc_global <- sum(diag(conf_match)) / sum(conf_match)
  sens <- diag(conf_match) / rowSums(conf_match)
  cat(sprintf("\nAcuracia global apos pareamento hungaro: %.1f%%\n",
              100 * acc_global))
  cat("Sensibilidade por grupo:\n")
  for (g in seq_along(sens))
    cat(sprintf("  G%d %-26s: %.1f%%\n", g,
                paste0("(", rownames(conf)[g], ")"), 100 * sens[g]))
} else {
  cat(sprintf("\nk_final=%d != 4 grupos verdadeiros — pareamento limitado.\n",
              k_final))
  acc_global <- NA; sens <- rep(NA, 4)
}


# ── FIG 6 · Densidade ajustada vs grupos verdadeiros ──────────────────────────
pars_final <- as.numeric(parameters(fit_final))
mus_final <- sort(exp(pars_final))
prior_final <- prior(fit_final)
prior_ord <- prior_final[order(exp(pars_final))]

grade_x <- seq(0, quantile(consultas, 0.97), length.out = 400)
dens_comp <- do.call(rbind, lapply(seq_along(mus_final), function(h) {
  data.frame(
    x = grade_x,
    densidade = prior_ord[h] * dnbinom(round(grade_x), size = THETA_NB,
                                        mu = mus_final[h]),
    componente = factor(sprintf("Comp. %d (mu=%.1f)", h, mus_final[h]))
  )
}))

dens_total <- data.frame(
  x = grade_x,
  densidade = sapply(grade_x, function(xv) {
    sum(prior_ord * dnbinom(round(xv), size = THETA_NB, mu = mus_final))
  })
)

fig06 <- ggplot() +
  geom_histogram(data = df %>% filter(consultas <= max(grade_x)),
                 aes(x = consultas, y = after_stat(density)),
                 bins = 50, fill = "grey80", colour = "white", alpha = 0.7) +
  geom_line(data = dens_comp, aes(x = x, y = densidade, colour = componente),
            linewidth = 0.9) +
  geom_line(data = dens_total, aes(x = x, y = densidade),
            colour = "black", linewidth = 1, linetype = "dashed") +
  labs(
    title    = sprintf("Modelo final ajustado: k = %d componentes (init: %s)",
                       k_final, LABEL_INIT[init_final]),
    subtitle = sprintf("Acuracia global apos pareamento hungaro: %s",
                       ifelse(is.na(acc_global), "n/a",
                              sprintf("%.1f%%", 100 * acc_global))),
    x = "Numero de consultas",
    y = "Densidade",
    colour = "Componente estimado"
  ) +
  TEMA

salvar_fig(fig06, "fig06_modelo_final", largura = 9, altura = 5)


# =============================================================================
# 10b. MODULO E · Avaliacao comparativa por metodo (k = K_TRUE)
# =============================================================================
# Avalia os 4 metodos no mesmo k = K_TRUE = 4 (cardinalidade verdadeira),
# gerando: metricas (ARI/Acuracia/Pureza/F1), recuperacao parametrica,
# diagnostico por grupo, ARI cruzado, e visualizacoes 1D analogas ao
# "PCA-by-group" do caso multivariado.
# =============================================================================

cat("\n", strrep("=", 70), "\n", sep = "")
cat(sprintf("MODULO E · Avaliacao comparativa por metodo (k = %d)\n", K_TRUE))
cat(strrep("=", 70), "\n\n", sep = "")

K_AVALIA   <- K_TRUE
grupo_int  <- as.integer(df$grupo_real)
ROTULOS_PT <- c("baixo uso", "ambulatorial", "agudo / hospitalar",
                "padrão atípico")
G_LABEL    <- sprintf("G%d · %s", 1:4, ROTULOS_PT)


# ── Helper: re-rotular estimativas via pareamento hungaro ─────────────────────
remap_hungarian <- function(real_int, est_int) {
  est_levels <- sort(unique(est_int))
  conf <- table(real = real_int,
                est  = factor(est_int, levels = est_levels))
  assign_h <- hungarian_match(conf)  # assign_h[g] = coluna de conf p/ grupo real g
  remap <- setNames(rep(NA_integer_, length(est_levels)), as.character(est_levels))
  for (g in seq_along(assign_h)) {
    if (!is.na(assign_h[g]) && assign_h[g] <= length(est_levels)) {
      remap[as.character(est_levels[assign_h[g]])] <- g
    }
  }
  as.integer(remap[as.character(est_int)])
}

# ── Helper: metricas para um fit no k = K_TRUE ───────────────────────────────
metricas_metodo <- function(fit, real_int) {
  if (is.null(fit)) return(NULL)
  est     <- as.integer(clusters(fit))
  est_par <- remap_hungarian(real_int, est)

  ari      <- adj_rand_index(real_int, est)
  acuracia <- mean(est_par == real_int, na.rm = TRUE)
  tab      <- table(real_int, est)
  pureza   <- sum(apply(tab, 2L, max)) / sum(tab)

  k_real <- length(unique(real_int))
  f1_grp <- vapply(seq_len(k_real), function(g) {
    tp <- sum(real_int == g & est_par == g, na.rm = TRUE)
    fp <- sum(real_int != g & est_par == g, na.rm = TRUE)
    fn <- sum(real_int == g & est_par != g, na.rm = TRUE)
    if (tp + fp == 0 || tp + fn == 0) return(0)
    prec <- tp / (tp + fp); rec <- tp / (tp + fn)
    if (prec + rec == 0) 0 else 2 * prec * rec / (prec + rec)
  }, numeric(1))

  sens <- vapply(seq_len(k_real), function(g) {
    n_g <- sum(real_int == g)
    if (n_g == 0) NA_real_ else sum(est_par == g & real_int == g) / n_g
  }, numeric(1))

  list(
    ARI = ari, Acuracia = acuracia, Pureza = pureza,
    F1_macro = mean(f1_grp), F1_grp = f1_grp, Sens = sens,
    est = est, est_par = est_par
  )
}

# ── Coleta dos fits e metricas no k = K_TRUE para os 4 metodos ───────────────
fits_E       <- list()
metricas_E   <- list()
mapped_assign <- list(verdadeiro = grupo_int)  # adiciona referencia "true"

for (init_nm in names(INITS)) {
  fit_m <- fits[[init_nm]][[as.character(K_AVALIA)]]
  fits_E[[init_nm]] <- fit_m
  if (is.null(fit_m)) {
    cat(sprintf("  %-8s : FIT INDISPONIVEL (colapsou em k=%d)\n",
                init_nm, K_AVALIA))
    metricas_E[[init_nm]] <- NULL
    mapped_assign[[init_nm]] <- rep(NA_integer_, N)
    next
  }
  m <- metricas_metodo(fit_m, grupo_int)
  metricas_E[[init_nm]] <- m
  mapped_assign[[init_nm]] <- m$est_par
  cat(sprintf("  %-8s : ARI=%.3f | Acur=%.3f | Pureza=%.3f | F1=%.3f\n",
              init_nm, m$ARI, m$Acuracia, m$Pureza, m$F1_macro))
}


# ── Tabela 5 · Metricas globais por metodo ───────────────────────────────────
tab05 <- do.call(rbind, lapply(names(INITS), function(nm) {
  m <- metricas_E[[nm]]
  if (is.null(m)) {
    data.frame(metodo = LABEL_INIT[nm],
               ARI = NA, Acuracia = NA, Pureza = NA, F1_macro = NA,
               status = "colapsou")
  } else {
    data.frame(metodo = LABEL_INIT[nm],
               ARI = round(m$ARI, 3),
               Acuracia = round(m$Acuracia, 3),
               Pureza = round(m$Pureza, 3),
               F1_macro = round(m$F1_macro, 3),
               status = "ok")
  }
}))
write.csv(tab05, file.path(RDS_DIR, "tab05_metricas_por_metodo.csv"),
          row.names = FALSE)

gt05 <- gt(tab05) %>%
  tab_header(
    title    = md(sprintf("**Tabela 5.** Métricas de classificação por método (*k* = %d)",
                          K_AVALIA)),
    subtitle = md("Pareamento via algoritmo húngaro. F1 macro = média não ponderada das F1 por grupo.")
  ) %>%
  cols_label(metodo = "Método", ARI = "ARI", Acuracia = "Acurácia",
             Pureza = "Pureza", F1_macro = "F1 (macro)",
             status = "Status") %>%
  fmt_number(columns = c(ARI, Acuracia, Pureza, F1_macro), decimals = 3) %>%
  tab_style(style = cell_text(color = "white", weight = "bold"),
            locations = cells_column_labels(everything())) %>%
  tab_options(table.font.size = px(12),
              column_labels.background.color = "#2c3e50",
              column_labels.border.top.style = "hidden",
              heading.align = "left")

# Salva via gtsave (mesmo padrao das outras tabelas)
tryCatch({
  gtsave(gt05, file.path(IMG_DIR, "tab05_metricas_por_metodo.png"),
         vwidth = 900, vheight = 320, zoom = 2, expand = 5)
  message("  -> salvo: tab05_metricas_por_metodo.png")
}, error = function(e) message("  FALHA tab05: ", conditionMessage(e)))


# ── Tabela 6 · Recuperacao parametrica (mu, pi) ──────────────────────────────
tab06 <- list()
for (nm in names(INITS)) {
  m <- metricas_E[[nm]]; fit_m <- fits_E[[nm]]
  if (is.null(m) || is.null(fit_m)) next
  # Ordena componentes pelo mu estimado, depois aplica hungarian para o grupo verdadeiro
  pars_m   <- as.numeric(parameters(fit_m))
  mus_est  <- exp(pars_m)
  pesos_est <- as.numeric(prior(fit_m))

  est_levels <- sort(unique(m$est))
  conf <- table(real = grupo_int,
                est  = factor(m$est, levels = est_levels))
  assign_h <- hungarian_match(conf)
  # assign_h[g] = posicao em est_levels que mapeia para grupo real g
  mu_alinhado  <- vapply(seq_along(assign_h), function(g) {
    if (is.na(assign_h[g]) || assign_h[g] > length(mus_est)) NA_real_
    else mus_est[assign_h[g]]
  }, numeric(1))
  pi_alinhado <- vapply(seq_along(assign_h), function(g) {
    if (is.na(assign_h[g]) || assign_h[g] > length(pesos_est)) NA_real_
    else pesos_est[assign_h[g]]
  }, numeric(1))

  for (g in 1:4) {
    tab06[[length(tab06) + 1]] <- data.frame(
      metodo  = LABEL_INIT[nm],
      grupo   = G_LABEL[g],
      mu_real = MU_TRUE[g],
      mu_est  = round(mu_alinhado[g], 2),
      mu_err  = round(mu_alinhado[g] - MU_TRUE[g], 2),
      pi_real = PI_TRUE[g],
      pi_est  = round(pi_alinhado[g], 3),
      pi_err  = round(pi_alinhado[g] - PI_TRUE[g], 3)
    )
  }
}
tab06 <- do.call(rbind, tab06)
write.csv(tab06, file.path(RDS_DIR, "tab06_recuperacao_parametros.csv"),
          row.names = FALSE)

gt06 <- gt(tab06, groupname_col = "metodo") %>%
  tab_header(
    title    = md(sprintf("**Tabela 6.** Recuperação paramétrica por método (*k* = %d)",
                          K_AVALIA)),
    subtitle = md("μ e π estimados, alinhados aos grupos verdadeiros via algoritmo húngaro. Erro = estimado − verdadeiro.")
  ) %>%
  cols_label(grupo = "Grupo",
             mu_real = md("μ real"), mu_est = md("μ estimado"), mu_err = md("Erro μ"),
             pi_real = md("π real"), pi_est = md("π estimado"), pi_err = md("Erro π")) %>%
  fmt_number(columns = c(mu_est, mu_err), decimals = 2) %>%
  fmt_number(columns = c(pi_real, pi_est, pi_err), decimals = 3) %>%
  tab_style(style = cell_text(color = "white", weight = "bold"),
            locations = cells_column_labels(everything())) %>%
  tab_options(table.font.size = px(11),
              column_labels.background.color = "#2c3e50",
              column_labels.border.top.style = "hidden",
              heading.align = "left",
              row_group.background.color = "#ecf0f1",
              row_group.font.weight = "bold")

tryCatch({
  gtsave(gt06, file.path(IMG_DIR, "tab06_recuperacao_parametros.png"),
         vwidth = 950, vheight = 700, zoom = 2, expand = 5)
  message("  -> salvo: tab06_recuperacao_parametros.png")
}, error = function(e) message("  FALHA tab06: ", conditionMessage(e)))


# ── Tabela 7 · Diagnostico por grupo (sensibilidade, F1, n) ──────────────────
tab07 <- list()
for (nm in names(INITS)) {
  m <- metricas_E[[nm]]
  if (is.null(m)) next
  for (g in 1:4) {
    n_real    <- sum(grupo_int == g)
    n_correto <- sum(m$est_par == g & grupo_int == g, na.rm = TRUE)
    tab07[[length(tab07) + 1]] <- data.frame(
      metodo = LABEL_INIT[nm],
      grupo  = G_LABEL[g],
      n_real = n_real,
      n_correto = n_correto,
      sensibilidade = round(m$Sens[g], 3),
      F1 = round(m$F1_grp[g], 3)
    )
  }
}
tab07 <- do.call(rbind, tab07)
write.csv(tab07, file.path(RDS_DIR, "tab07_diagnostico_por_grupo.csv"),
          row.names = FALSE)

gt07 <- gt(tab07, groupname_col = "metodo") %>%
  tab_header(
    title    = md(sprintf("**Tabela 7.** Diagnóstico por grupo (*k* = %d)", K_AVALIA)),
    subtitle = md("Sensibilidade = proporção do grupo verdadeiro corretamente classificada.")
  ) %>%
  cols_label(grupo = "Grupo", n_real = md("*n* verdadeiro"),
             n_correto = md("*n* correto"),
             sensibilidade = "Sensibilidade", F1 = "F1") %>%
  fmt_number(columns = c(sensibilidade, F1), decimals = 3) %>%
  tab_style(style = cell_text(color = "white", weight = "bold"),
            locations = cells_column_labels(everything())) %>%
  tab_options(table.font.size = px(11),
              column_labels.background.color = "#2c3e50",
              column_labels.border.top.style = "hidden",
              heading.align = "left",
              row_group.background.color = "#ecf0f1",
              row_group.font.weight = "bold")

tryCatch({
  gtsave(gt07, file.path(IMG_DIR, "tab07_diagnostico_por_grupo.png"),
         vwidth = 800, vheight = 650, zoom = 2, expand = 5)
  message("  -> salvo: tab07_diagnostico_por_grupo.png")
}, error = function(e) message("  FALHA tab07: ", conditionMessage(e)))


# ── Tabela 8 · ARI cruzado entre metodos ─────────────────────────────────────
# IMPORTANTE: metricas_E[[init_nm]] <- NULL APAGA o elemento da lista em R,
# entao names(metricas_E) ja contem apenas os metodos que convergiram.
metodos_validos <- names(metricas_E)
ari_mat <- matrix(NA_real_,
                  nrow = length(metodos_validos),
                  ncol = length(metodos_validos),
                  dimnames = list(LABEL_INIT[metodos_validos],
                                  LABEL_INIT[metodos_validos]))
for (i in seq_along(metodos_validos)) {
  for (j in seq_along(metodos_validos)) {
    ari_mat[i, j] <- adj_rand_index(
      metricas_E[[metodos_validos[i]]]$est,
      metricas_E[[metodos_validos[j]]]$est
    )
  }
}
tab08 <- data.frame(metodo_i = rownames(ari_mat), round(ari_mat, 3),
                    check.names = FALSE)
write.csv(tab08, file.path(RDS_DIR, "tab08_ari_cruzado.csv"), row.names = FALSE)

gt08 <- gt(tab08) %>%
  tab_header(
    title    = md(sprintf("**Tabela 8.** ARI cruzado entre métodos (*k* = %d)", K_AVALIA)),
    subtitle = md("ARI = 1 indica partição idêntica. Valores próximos de 1 indicam métodos com soluções equivalentes.")
  ) %>%
  cols_label(metodo_i = "Método") %>%
  fmt_number(columns = -metodo_i, decimals = 3) %>%
  tab_style(style = cell_text(color = "white", weight = "bold"),
            locations = cells_column_labels(everything())) %>%
  tab_options(table.font.size = px(12),
              column_labels.background.color = "#2c3e50",
              column_labels.border.top.style = "hidden",
              heading.align = "left")

tryCatch({
  gtsave(gt08, file.path(IMG_DIR, "tab08_ari_cruzado.png"),
         vwidth = 700, vheight = 280, zoom = 2, expand = 5)
  message("  -> salvo: tab08_ari_cruzado.png")
}, error = function(e) message("  FALHA tab08: ", conditionMessage(e)))


# ── FIG 7 · Densidades ajustadas por metodo (4 paineis lado a lado) ─────────
LABEL_PANEL <- c(
  random  = "Random multistart",
  kmeans  = "k-means + EM",
  cem_em  = "CEM + EM",
  sem_em  = "SEM + EM"
)

dens_metodos <- list()
for (nm in names(INITS)) {
  fit_m <- fits_E[[nm]]
  if (is.null(fit_m)) next
  mus_est   <- sort(exp(as.numeric(parameters(fit_m))))
  ord       <- order(exp(as.numeric(parameters(fit_m))))
  pesos_est <- as.numeric(prior(fit_m))[ord]
  grade_x   <- seq(0, quantile(consultas, 0.97), length.out = 300)

  for (h in seq_along(mus_est)) {
    dens_metodos[[length(dens_metodos) + 1]] <- data.frame(
      metodo     = LABEL_PANEL[nm],
      componente = sprintf("C%d (μ=%.1f)", h, mus_est[h]),
      x          = grade_x,
      densidade  = pesos_est[h] * dnbinom(round(grade_x), size = THETA_NB,
                                          mu = mus_est[h])
    )
  }
}
dens_metodos <- do.call(rbind, dens_metodos)
dens_metodos$metodo <- factor(dens_metodos$metodo,
                              levels = LABEL_PANEL[names(LABEL_PANEL) %in%
                                                   names(fits_E)[!sapply(fits_E, is.null)]])

df_hist <- df[df$consultas <= quantile(consultas, 0.97), ]

fig07 <- ggplot() +
  geom_histogram(data = df_hist,
                 aes(x = consultas, y = after_stat(density)),
                 bins = 50, fill = "grey85", colour = "white", alpha = 0.7) +
  geom_line(data = dens_metodos,
            aes(x = x, y = densidade, colour = componente),
            linewidth = 0.85) +
  facet_wrap(~ metodo, ncol = 2) +
  scale_colour_brewer(palette = "Set1") +
  labs(
    x = "Número de consultas",
    y = "Densidade",
    colour = "Componente estimado",
    title = sprintf("Densidades ajustadas por método (*k* = %d)", K_AVALIA),
    subtitle = "Histograma cinza = dados; linhas coloridas = componentes NB estimados. CEM colapsou neste cenário."
  ) +
  TEMA + theme(legend.position = "right", plot.title = element_text(face = "bold"))

salvar_fig(fig07, "fig07_densidades_por_metodo", largura = 11, altura = 6)


# ── FIG 8 · Grid 4x5 — classificacao por grupo (analogo ao PCA-by-group) ─────
LABEL_PANEL_FULL <- c(verdadeiro = "verdadeiro", LABEL_PANEL)

grid_data <- do.call(rbind, lapply(seq_len(4), function(g_focus) {
  do.call(rbind, lapply(names(mapped_assign), function(meth) {
    data.frame(
      g_focus_lbl  = G_LABEL[g_focus],
      method_lbl   = LABEL_PANEL_FULL[meth],
      consultas    = consultas,
      in_focus     = !is.na(mapped_assign[[meth]]) &
                     mapped_assign[[meth]] == g_focus
    )
  }))
}))
grid_data <- grid_data[grid_data$consultas <= quantile(consultas, 0.97), ]
grid_data$g_focus_lbl <- factor(grid_data$g_focus_lbl, levels = G_LABEL)
grid_data$method_lbl  <- factor(grid_data$method_lbl,
                                levels = unname(LABEL_PANEL_FULL))

# Cor por linha (grupo em foco)
COR_FOCUS <- setNames(unname(COR_GRUPO[c("baixo_uso",
                                          "ambulatorial_coordenado",
                                          "agudo_hospitalar",
                                          "atipico")]),
                      G_LABEL)

fig08 <- ggplot(grid_data, aes(x = consultas, y = 0)) +
  geom_jitter(data = subset(grid_data, !in_focus),
              colour = "grey85", alpha = 0.20, size = 0.35, height = 0.4) +
  geom_jitter(data = subset(grid_data, in_focus),
              aes(colour = g_focus_lbl),
              alpha = 0.65, size = 0.5, height = 0.4) +
  facet_grid(g_focus_lbl ~ method_lbl, switch = "y") +
  scale_colour_manual(values = COR_FOCUS, guide = "none") +
  scale_y_continuous(limits = c(-1, 1), breaks = NULL) +
  labs(
    x = "Número de consultas",
    y = NULL,
    title = sprintf("Classificação por grupo — k = %d (univariado)", K_AVALIA),
    subtitle = "Cada linha destaca um grupo verdadeiro · cor viva = em foco · cinza = demais"
  ) +
  TEMA + theme(
    legend.position = "none",
    strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold"),
    strip.text.x      = element_text(face = "bold", size = 9),
    panel.spacing.x   = unit(0.3, "lines")
  )

salvar_fig(fig08, "fig08_classificacao_grid", largura = 14, altura = 7)


# ── FIG 9 · Certeza posterior (max tau) por metodo e grupo verdadeiro ────────
cert_data <- list()
for (nm in names(INITS)) {
  fit_m <- fits_E[[nm]]
  if (is.null(fit_m)) next
  post <- posterior(fit_m)
  max_tau <- apply(post, 1L, max)
  cert_data[[length(cert_data) + 1]] <- data.frame(
    metodo       = LABEL_PANEL[nm],
    grupo_real   = G_LABEL[grupo_int],
    max_tau      = max_tau
  )
}
cert_data <- do.call(rbind, cert_data)
cert_data$metodo     <- factor(cert_data$metodo,
                               levels = LABEL_PANEL[names(LABEL_PANEL) %in%
                                                    names(fits_E)[!sapply(fits_E, is.null)]])
cert_data$grupo_real <- factor(cert_data$grupo_real, levels = G_LABEL)

fig09 <- ggplot(cert_data, aes(x = max_tau, fill = grupo_real, colour = grupo_real)) +
  geom_density(alpha = 0.30, linewidth = 0.7) +
  geom_vline(xintercept = 0.9, linetype = "dotted", colour = "grey40") +
  facet_wrap(~ metodo, ncol = 2) +
  scale_fill_manual(values = setNames(unname(COR_GRUPO[c("baixo_uso",
                                                           "ambulatorial_coordenado",
                                                           "agudo_hospitalar",
                                                           "atipico")]),
                                        G_LABEL),
                    name = "Grupo verdadeiro") +
  scale_colour_manual(values = setNames(unname(COR_GRUPO[c("baixo_uso",
                                                             "ambulatorial_coordenado",
                                                             "agudo_hospitalar",
                                                             "atipico")]),
                                          G_LABEL),
                      name = "Grupo verdadeiro") +
  scale_x_continuous(limits = c(0.25, 1), breaks = seq(0.3, 1, 0.2)) +
  labs(
    x = expression(max[k] ~ hat(tau)[ik]),
    y = "Densidade",
    title = "Certeza máxima de classificação por método e grupo verdadeiro",
    subtitle = "Linha pontilhada = limiar 0,90. G3 e G4 tendem a ter certeza menor (zona de sobreposição)."
  ) +
  TEMA + theme(plot.title = element_text(face = "bold"))

salvar_fig(fig09, "fig09_certeza_por_grupo", largura = 11, altura = 6)


# ── FIG 10 · Taxa de erro por faixa do numero de consultas ──────────────────
faixas    <- c(0, 5, 15, 35, Inf)
faixa_lbl <- c("0–4", "5–14", "15–34", "≥ 35")
faixa_i   <- cut(consultas, breaks = faixas, right = FALSE, labels = faixa_lbl)

erro_data <- list()
for (nm in names(INITS)) {
  m <- metricas_E[[nm]]
  if (is.null(m)) next
  err <- m$est_par != grupo_int
  for (fx in faixa_lbl) {
    idx <- faixa_i == fx & !is.na(faixa_i)
    if (sum(idx) == 0) next
    erro_data[[length(erro_data) + 1]] <- data.frame(
      metodo = LABEL_PANEL[nm],
      faixa  = fx,
      n      = sum(idx),
      err    = mean(err[idx], na.rm = TRUE)
    )
  }
}
erro_data <- do.call(rbind, erro_data)
erro_data$faixa  <- factor(erro_data$faixa, levels = faixa_lbl)
erro_data$metodo <- factor(erro_data$metodo,
                           levels = LABEL_PANEL[names(LABEL_PANEL) %in%
                                                 names(fits_E)[!sapply(fits_E, is.null)]])

fig10 <- ggplot(erro_data, aes(x = faixa, y = err, fill = metodo)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", 100 * err)),
            position = position_dodge(width = 0.75),
            vjust = -0.3, size = 3) +
  scale_fill_manual(values = setNames(unname(COR_INIT[names(LABEL_PANEL)]),
                                       LABEL_PANEL),
                    name = "Método") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.12))) +
  labs(
    x = "Faixa do número de consultas",
    y = "Taxa de erro",
    title = "Erro de classificação por faixa do número de consultas",
    subtitle = "Zona de sobreposição G3↔G4 (faixa 15–34) concentra os maiores erros."
  ) +
  TEMA + theme(plot.title = element_text(face = "bold"))

salvar_fig(fig10, "fig10_erro_por_faixa", largura = 11, altura = 5.5)


cat("\nModulo E concluido.\n")
cat(sprintf("  Tabelas:  tab05 (metricas) · tab06 (parametros) · tab07 (por grupo) · tab08 (ARI cruzado)\n"))
cat(sprintf("  Figuras:  fig07 (densidades) · fig08 (grid 4x5) · fig09 (certeza) · fig10 (erro/faixa)\n"))


# =============================================================================
# 11. TABELAS FORMATADAS (gt → PNG + CSV)
# =============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("Gerando tabelas formatadas (gt → PNG + CSV)\n")
cat(strrep("=", 70), "\n\n", sep = "")

salvar_gt <- function(gt_obj, nome, w_in = 8, h_in = 5) {
  png_path <- file.path(IMG_DIR, paste0(nome, ".png"))
  tryCatch({
    gtsave(gt_obj, png_path,
           vwidth = round(w_in * 100), vheight = round(h_in * 100),
           zoom = 2, expand = 5)
    message("  -> salvo: ", nome, ".png")
  }, error = function(e) message("  FALHA em ", nome, ": ", conditionMessage(e)))
}


# ── Tabela 1: resultados por inicializacao × k ────────────────────────────────
tab01 <- tab_init %>%
  mutate(init = LABEL_INIT[init]) %>%
  arrange(init, k)

gt01 <- gt(tab01) %>%
  tab_header(
    title    = md("**Tabela 1.** Comparacao de inicializacoes × numero de componentes"),
    subtitle = md(sprintf("Mistura NB univariada (theta=%.1f, n=%d) — %d reinicios por estrategia",
                          THETA_NB, N, NREP))
  ) %>%
  cols_label(
    init = "Inicializacao", k = md("*k*"), k_efetivo = md("*k* efetivo"),
    logLik = "log-Lik", np = md("*p*"),
    AIC = "AIC", BIC = "BIC", ICL = "ICL", AIC3 = "AIC3",
    tempo_s = "Tempo (s)"
  ) %>%
  fmt_number(columns = c(logLik, AIC, BIC, ICL, AIC3), decimals = 1) %>%
  fmt_number(columns = tempo_s, decimals = 1) %>%
  data_color(columns = BIC,
             colors = scales::col_numeric(
               palette = c("#27ae60", "#f1c40f", "#e74c3c"),
               domain = NULL)) %>%
  tab_style(style = cell_text(color = "white", weight = "bold"),
            locations = cells_column_labels(everything())) %>%
  tab_options(
    table.font.size = px(11),
    column_labels.background.color = "#2c3e50",
    column_labels.border.top.style = "hidden",
    heading.align = "left",
    heading.title.font.size = px(13)
  )
salvar_gt(gt01, "tab01_inicializacoes", w_in = 9, h_in = 5)


# ── Tabela 2: criterios por k (melhor inicializacao) ──────────────────────────
tab02 <- melhor_por_k %>%
  mutate(init = LABEL_INIT[init],
         escolhido_por = sapply(k, function(kk) {
           crits <- c()
           if (kk == k_aic)  crits <- c(crits, "AIC")
           if (kk == k_bic)  crits <- c(crits, "BIC")
           if (kk == k_icl)  crits <- c(crits, "ICL")
           if (kk == k_aic3) crits <- c(crits, "AIC3")
           if (kk == k_blrt_escolhido) crits <- c(crits, "BLRT")
           if (!length(crits)) "" else paste(crits, collapse = ", ")
         })) %>%
  select(k, init, logLik, AIC, BIC, ICL, AIC3, tempo_s, escolhido_por)

write.csv(tab02, file.path(RDS_DIR, "tab02_criterios_melhor_init.csv"),
          row.names = FALSE)

gt02 <- gt(tab02) %>%
  tab_header(
    title    = md("**Tabela 2.** Criterios de selecao na melhor inicializacao por *k*"),
    subtitle = md(sprintf("Melhor inicializacao por *k* = maior log-verossimilhanca | *k* verdadeiro = %d",
                          K_TRUE))
  ) %>%
  cols_label(k = md("*k*"), init = "Melhor inicializacao",
             logLik = "log-Lik", tempo_s = "Tempo (s)",
             escolhido_por = "Selecionado por") %>%
  fmt_number(columns = c(logLik, AIC, BIC, ICL, AIC3), decimals = 1) %>%
  fmt_number(columns = tempo_s, decimals = 1) %>%
  tab_style(style = list(cell_fill(color = "#fdf6e3"), cell_text(weight = "bold")),
            locations = cells_body(rows = k == k_final)) %>%
  tab_style(style = cell_text(color = "white", weight = "bold"),
            locations = cells_column_labels(everything())) %>%
  tab_options(
    table.font.size = px(11),
    column_labels.background.color = "#2c3e50",
    column_labels.border.top.style = "hidden",
    heading.align = "left"
  )
salvar_gt(gt02, "tab02_criterios_por_k", w_in = 9, h_in = 4)


# ── Tabela 3: BLRT sequencial ─────────────────────────────────────────────────
if (length(resultados_blrt)) {
  gt03 <- gt(tab_blrt) %>%
    tab_header(
      title = md("**Tabela 3.** Bootstrap LRT sequencial (H0: *k* vs H1: *k*+1)"),
      subtitle = md(sprintf("B = %d replicacoes paramétricas | k* (BLRT) = %d",
                            BLRT_B, k_blrt_escolhido))
    ) %>%
    cols_label(k0 = md("*k*"), k1 = md("*k*+1"),
               LRT_obs = "LRT observado", p_valor = "p-valor",
               B_validos = md("*B* validos"), rejeita_H0 = "Rejeita H0") %>%
    fmt_number(columns = LRT_obs, decimals = 1) %>%
    fmt_number(columns = p_valor, decimals = 4) %>%
    tab_style(style = cell_fill(color = "#fdf6e3"),
              locations = cells_body(rows = rejeita_H0 == "Nao")) %>%
    tab_style(style = cell_text(color = "white", weight = "bold"),
              locations = cells_column_labels(everything())) %>%
    tab_options(
      table.font.size = px(11),
      column_labels.background.color = "#2c3e50",
      column_labels.border.top.style = "hidden",
      heading.align = "left"
    )
  salvar_gt(gt03, "tab03_blrt", w_in = 8, h_in = 3.5)
}


# ── Tabela 4: classificacao do modelo final ───────────────────────────────────
if (!is.na(acc_global)) {
  tab04 <- data.frame(
    grupo_real = rownames(conf),
    mu_verdadeiro = MU_TRUE,
    n_verdadeiro = as.integer(rowSums(conf)),
    componente_atribuido = assign_h,
    n_correto = as.integer(diag(conf_match)),
    sensibilidade_pct = round(100 * sens, 1)
  )
  write.csv(tab04, file.path(RDS_DIR, "tab04_modelo_final.csv"),
            row.names = FALSE)

  gt04 <- gt(tab04) %>%
    tab_header(
      title    = md(sprintf("**Tabela 4.** Avaliacao do modelo final (*k* = %d)",
                            k_final)),
      subtitle = md(sprintf("Inicializacao: %s | Acuracia global = %.1f%% | ARI medio (bootstrap) = %.3f",
                            LABEL_INIT[init_final], 100 * acc_global, ari_med))
    ) %>%
    cols_label(
      grupo_real = "Grupo verdadeiro",
      mu_verdadeiro = md("*mu* verdadeiro"),
      n_verdadeiro = md("*n* verdadeiro"),
      componente_atribuido = "Componente",
      n_correto = md("*n* correto"),
      sensibilidade_pct = "Sensibilidade (%)"
    ) %>%
    data_color(columns = sensibilidade_pct,
               colors = scales::col_numeric(
                 palette = c("#e74c3c", "#f1c40f", "#27ae60"),
                 domain = c(0, 100))) %>%
    tab_style(style = cell_text(color = "white", weight = "bold"),
              locations = cells_column_labels(everything())) %>%
    tab_options(
      table.font.size = px(11),
      column_labels.background.color = "#2c3e50",
      column_labels.border.top.style = "hidden",
      heading.align = "left"
    )
  salvar_gt(gt04, "tab04_modelo_final", w_in = 9, h_in = 3.5)
}


# =============================================================================
# 12. RESUMO FINAL
# =============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("RESUMO\n")
cat(strrep("=", 70), "\n", sep = "")
cat(sprintf("Saidas em: %s\n", IMG_DIR))
cat("Figuras:\n")
cat("  fig01_dados_simulados.png          (distribuicoes verdadeiras)\n")
cat("  fig02_loglik_por_init.png          (logLik vs k por estrategia)\n")
cat("  fig03_criterios_por_k.png          (AIC/BIC/ICL/AIC3 vs k)\n")
cat("  fig04_blrt_distribuicoes.png       (bootstrap LRT, multi-painel)\n")
cat("  fig05_estabilidade_ari.png         (bootstrap ARI)\n")
cat("  fig06_modelo_final.png             (modelo ajustado vs dados)\n")
cat("Tabelas (PNG + CSV em RDS_DIR):\n")
cat("  tab01_inicializacoes               (init × k: logLik/BIC/ICL/AIC/tempo)\n")
cat("  tab02_criterios_por_k              (k vs todos os criterios)\n")
cat("  tab03_blrt                         (BLRT sequencial)\n")
cat("  tab04_modelo_final                 (sensibilidade por grupo)\n")
cat(sprintf("\nk_final = %d  |  init = %s  |  acuracia = %s\n",
            k_final, init_final,
            ifelse(is.na(acc_global), "n/a", sprintf("%.1f%%", 100 * acc_global))))
cat(strrep("=", 70), "\n", sep = "")
