# =============================================================================
# 08_metricas_multivariado_v3.R
#
# Análise de métricas substantivas dos fits NB multivariados em V3.
# Reusa os fits já gerados pelo script 07 (cache fits_v3.rds), sem refit.
#
# Calcula, para cada g em {3, 4, 5}:
#   - Métricas globais: ARI, Acurácia, Pureza, F1 macro
#   - Diagnóstico por grupo: n real, n correto, sensibilidade, F1,
#     AvePP (probabilidade posterior média no componente correto)
#
# Saídas:
#   docs/assets/img/blrt_multivariado_v3/
#     tab02_metricas_globais.png       (g x métrica)
#     tab03_metricas_por_grupo.png     (g x grupo x sens/F1/AvePP)
#     fig02_sensibilidade_por_grupo.png (barras facetadas por g)
#     fig03_matrizes_confusao.png       (3 matrizes lado a lado)
#   src/03_analises/rds/blrt_multivariado_v3/
#     tab02_metricas_globais.csv
#     tab03_metricas_por_grupo.csv
# =============================================================================


# ── 0. Pacotes ────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(gt); library(patchwork); library(tidyr)
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
cat("Raiz: ", PROJ_ROOT, "\n", sep = "")

IMG_DIR <- file.path(PROJ_ROOT, "docs/assets/img/blrt_multivariado_v3")
RDS_DIR <- file.path(PROJ_ROOT, "src/03_analises/rds/blrt_multivariado_v3")


# ── 2. Constantes ─────────────────────────────────────────────────────────────
ROTULOS <- c("baixo_uso", "ambulatorial_coordenado",
             "agudo_hospitalar", "atipico")
G_LABEL <- c("G1 · baixo uso", "G2 · ambulatorial",
             "G3 · agudo / hospitalar", "G4 · atípico")
COR_GRUPO <- c(
  "G1 · baixo uso"           = "#4a7c59",
  "G2 · ambulatorial"        = "#2980b9",
  "G3 · agudo / hospitalar"  = "#e8a33c",
  "G4 · atípico"             = "#c0392b"
)


# ── 3. Carregar fits + dados ─────────────────────────────────────────────────
fits  <- readRDS(file.path(RDS_DIR, "fits_v3.rds"))
dados <- read.csv("dados/dados_simulados.csv", stringsAsFactors = FALSE)

grupo_int <- as.integer(factor(dados$grupo_latente, levels = ROTULOS))
N         <- length(grupo_int)
G_VALS    <- as.integer(names(fits))   # 3, 4, 5

cat(sprintf("n = %d · fits disponíveis: g = %s\n",
            N, paste(G_VALS, collapse = ", ")))
cat("Tamanhos reais por grupo:\n")
print(table(factor(dados$grupo_latente, levels = ROTULOS)))


# ── 4. Helpers ────────────────────────────────────────────────────────────────
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

# Algoritmo húngaro via permutação exaustiva (k <= 8)
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

# Re-rotula `est` para o espaço dos grupos verdadeiros via pareamento ótimo.
# Quando há mais componentes que grupos (ex: g=5 vs 4 reais), os componentes
# "extras" recebem NA — observações ali ficam fora da diagonal e contam erro.
remap_hungarian <- function(real_int, est_int) {
  est_levels <- sort(unique(est_int))
  conf <- table(real = real_int,
                est  = factor(est_int, levels = est_levels))
  assign_h <- hungarian_match(conf)
  remap <- setNames(rep(NA_integer_, length(est_levels)),
                    as.character(est_levels))
  for (g in seq_along(assign_h)) {
    if (!is.na(assign_h[g]) && assign_h[g] <= length(est_levels))
      remap[as.character(est_levels[assign_h[g]])] <- g
  }
  as.integer(remap[as.character(est_int)])
}

# F1 por grupo (binário: g vs resto)
f1_grupo <- function(real_int, est_par, g) {
  tp <- sum(real_int == g & est_par == g, na.rm = TRUE)
  fp <- sum(real_int != g & est_par == g, na.rm = TRUE)
  fn <- sum(real_int == g & est_par != g, na.rm = TRUE)
  if (tp + fp == 0 || tp + fn == 0) return(0)
  prec <- tp / (tp + fp); rec <- tp / (tp + fn)
  if (prec + rec == 0) 0 else 2 * prec * rec / (prec + rec)
}


# ── 5. Loop sobre g = 3, 4, 5 ────────────────────────────────────────────────
tab_global  <- list()
tab_grupo   <- list()
matrizes    <- list()
mappings    <- list()

for (g_val in G_VALS) {
  fit_g     <- fits[[as.character(g_val)]]
  est       <- as.integer(fit_g$classification)
  est_par   <- remap_hungarian(grupo_int, est)
  conf      <- table(real = grupo_int, est = est)
  est_lvls  <- sort(unique(est))
  assign_h  <- hungarian_match(conf)

  # Métricas globais
  ari <- adj_rand_index(grupo_int, est)
  acu <- mean(est_par == grupo_int, na.rm = TRUE)
  pur <- sum(apply(conf, 2L, max)) / sum(conf)
  f1g <- vapply(1:4, function(g) f1_grupo(grupo_int, est_par, g), numeric(1))
  f1m <- mean(f1g)

  tab_global[[length(tab_global) + 1]] <- data.frame(
    g            = g_val,
    ARI          = round(ari, 4),
    Acuracia     = round(acu, 4),
    Pureza       = round(pur, 4),
    F1_macro     = round(f1m, 4),
    logLik       = round(fit_g$loglik, 1),
    BIC          = round(fit_g$BIC, 1)
  )

  # Métricas por grupo verdadeiro
  for (g_real in 1:4) {
    n_real    <- sum(grupo_int == g_real)
    n_correto <- sum(est_par == g_real & grupo_int == g_real, na.rm = TRUE)
    sens      <- n_correto / n_real
    # AvePP: posterior média no componente Hungarian-matched a g_real
    if (!is.na(assign_h[g_real]) && assign_h[g_real] <= length(est_lvls)) {
      h_match  <- assign_h[g_real]
      idx_real <- which(grupo_int == g_real)
      avepp    <- mean(fit_g$tau[idx_real, h_match])
    } else {
      avepp <- NA_real_
    }
    tab_grupo[[length(tab_grupo) + 1]] <- data.frame(
      g         = g_val,
      grupo     = G_LABEL[g_real],
      n_real    = n_real,
      n_correto = n_correto,
      sensibilidade = round(sens, 4),
      F1            = round(f1g[g_real], 4),
      AvePP         = round(avepp, 4)
    )
  }

  # Salva matriz de confusão para visualização
  matrizes[[as.character(g_val)]] <- list(
    conf = conf, assign_h = assign_h,
    est_par = est_par, acu = acu
  )
  mappings[[as.character(g_val)]] <- est_par

  cat(sprintf("g=%d  ARI=%.3f Acu=%.3f Pur=%.3f F1=%.3f\n",
              g_val, ari, acu, pur, f1m))
}

tab_global <- do.call(rbind, tab_global)
tab_grupo  <- do.call(rbind, tab_grupo)

write.csv(tab_global, file.path(RDS_DIR, "tab02_metricas_globais.csv"),
          row.names = FALSE)
write.csv(tab_grupo,  file.path(RDS_DIR, "tab03_metricas_por_grupo.csv"),
          row.names = FALSE)


# ── 6. Tabela 2 (gt) — métricas globais ──────────────────────────────────────
gt02 <- gt(tab_global) %>%
  tab_header(
    title    = md("**Tabela 2.** Métricas globais por *g* — V3 multivariado"),
    subtitle = md(sprintf("Pareamento via algoritmo húngaro · n = %d", N))
  ) %>%
  cols_label(
    g = md("*g*"),
    ARI = "ARI", Acuracia = "Acurácia", Pureza = "Pureza",
    F1_macro = "F1 macro", logLik = "log-Lik", BIC = "BIC"
  ) %>%
  fmt_number(columns = c(ARI, Acuracia, Pureza, F1_macro), decimals = 3) %>%
  fmt_number(columns = c(logLik, BIC), decimals = 1) %>%
  tab_style(style = list(cell_fill(color = "#fdf6e3"),
                          cell_text(weight = "bold")),
            locations = cells_body(rows = g == 4)) %>%
  tab_style(style = cell_text(color = "white", weight = "bold"),
            locations = cells_column_labels(everything())) %>%
  tab_options(table.font.size = px(12),
              column_labels.background.color = "#2c3e50",
              column_labels.border.top.style = "hidden",
              heading.align = "left")

tryCatch({
  gtsave(gt02, file.path(IMG_DIR, "tab02_metricas_globais.png"),
         vwidth = 900, vheight = 220, zoom = 2, expand = 5)
  message("  -> tab02_metricas_globais.png")
}, error = function(e) message("  FALHA tab02: ", conditionMessage(e)))


# ── 7. Tabela 3 (gt) — diagnóstico por grupo, agrupada por g ─────────────────
gt03 <- gt(tab_grupo, groupname_col = "g") %>%
  tab_header(
    title    = md("**Tabela 3.** Diagnóstico por grupo verdadeiro × *g* — V3"),
    subtitle = md("Sensibilidade = proporção do grupo real corretamente classificada. AvePP = posterior média no componente correto.")
  ) %>%
  cols_label(
    grupo = "Grupo verdadeiro",
    n_real = md("*n* real"), n_correto = md("*n* correto"),
    sensibilidade = "Sensibilidade", F1 = "F1", AvePP = "AvePP"
  ) %>%
  fmt_number(columns = c(sensibilidade, F1, AvePP), decimals = 3) %>%
  data_color(
    columns = sensibilidade,
    fn = scales::col_numeric(
      palette = c("#e74c3c", "#f1c40f", "#27ae60"),
      domain  = c(0, 1)
    )
  ) %>%
  tab_style(style = cell_text(color = "white", weight = "bold"),
            locations = cells_column_labels(everything())) %>%
  tab_options(table.font.size = px(11),
              column_labels.background.color = "#2c3e50",
              column_labels.border.top.style = "hidden",
              heading.align = "left",
              row_group.background.color = "#ecf0f1",
              row_group.font.weight = "bold")

tryCatch({
  gtsave(gt03, file.path(IMG_DIR, "tab03_metricas_por_grupo.png"),
         vwidth = 950, vheight = 650, zoom = 2, expand = 5)
  message("  -> tab03_metricas_por_grupo.png")
}, error = function(e) message("  FALHA tab03: ", conditionMessage(e)))


# ── 8. Figura 2 — sensibilidade por grupo, facetado por g ────────────────────
df_sens <- tab_grupo %>%
  mutate(g_lbl = sprintf("g = %d", g),
         grupo = factor(grupo, levels = G_LABEL))

fig02 <- ggplot(df_sens, aes(x = grupo, y = sensibilidade, fill = grupo)) +
  geom_col(width = 0.75) +
  geom_text(aes(label = sprintf("%.1f%%", 100 * sensibilidade)),
            vjust = -0.4, size = 3) +
  geom_hline(yintercept = 0.80, linetype = "dashed",
             colour = "grey55", linewidth = 0.4) +
  facet_wrap(~ g_lbl, ncol = 3) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1.10),
                     breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
  scale_fill_manual(values = COR_GRUPO, guide = "none") +
  labs(
    title    = "Sensibilidade por grupo verdadeiro · V3 multivariado",
    subtitle = "Linha tracejada = limiar 80 %  ·  pareamento via algoritmo húngaro",
    x = NULL, y = "Sensibilidade"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 35, hjust = 1, size = 8),
        strip.text  = element_text(face = "bold"),
        panel.grid.minor = element_blank())

ggsave(filename = file.path(IMG_DIR, "fig02_sensibilidade_por_grupo.png"),
       plot = fig02, width = 12, height = 5, dpi = 150, bg = "white")
message("  -> fig02_sensibilidade_por_grupo.png")


# ── 9. Figura 3 — matrizes de confusão (g=3, 4, 5) ───────────────────────────
plot_conf <- function(g_val) {
  m       <- matrizes[[as.character(g_val)]]
  conf    <- m$conf
  est_par <- m$est_par
  # Recria matriz "alinhada": linhas = grupos reais, colunas = grupo mapeado
  est_lvls <- sort(unique(as.integer(colnames(conf))))
  conf_mat <- matrix(0L, 4, g_val)
  rownames(conf_mat) <- G_LABEL
  colnames(conf_mat) <- sprintf("C%d", seq_len(g_val))
  for (i in seq_len(nrow(conf)))
    for (j in seq_len(ncol(conf)))
      conf_mat[i, j] <- conf[i, j]

  df <- as.data.frame.table(conf_mat, responseName = "n")
  names(df)[1:2] <- c("real", "est")
  df$pct <- 100 * df$n / sum(df$n)

  ggplot(df, aes(x = est, y = real, fill = n)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = ifelse(n == 0, "",
                                  sprintf("%d\n(%.1f%%)", n, pct))),
              size = 2.9, colour = "grey15") +
    scale_fill_gradient(low = "#eef5fb", high = "#2980b9", guide = "none") +
    labs(title = sprintf("g = %d  (Acurácia = %.1f%%)",
                         g_val, 100 * m$acu),
         x = "Componente estimado", y = "Grupo real") +
    theme_minimal(base_size = 10) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(face = "bold", size = 11),
          axis.text.y = element_text(face = "bold"),
          axis.text.x = element_text(face = "bold"))
}

fig03 <- plot_conf(3) | plot_conf(4) | plot_conf(5)
fig03 <- fig03 + plot_annotation(
  title    = "Matrizes de confusão — V3 multivariado",
  subtitle = "Pareamento húngaro entre componentes estimados e grupos verdadeiros",
  theme    = theme_minimal()
)

ggsave(filename = file.path(IMG_DIR, "fig03_matrizes_confusao.png"),
       plot = fig03, width = 14, height = 5.5, dpi = 150, bg = "white")
message("  -> fig03_matrizes_confusao.png")


# ── 10. Resumo ────────────────────────────────────────────────────────────────
cat("\n", strrep("=", 70), "\n", sep = "")
cat("RESUMO — Métricas multivariadas V3\n")
cat(strrep("=", 70), "\n", sep = "")

cat("\nGlobais (por g):\n")
print(tab_global, row.names = FALSE)

cat("\nSensibilidade por grupo (por g):\n")
sens_wide <- tab_grupo %>%
  select(g, grupo, sensibilidade) %>%
  pivot_wider(names_from = g, values_from = sensibilidade,
              names_prefix = "g=")
print(as.data.frame(sens_wide), row.names = FALSE)

cat("\nSaídas em:", IMG_DIR, "\n")
cat("  tab02_metricas_globais.png\n")
cat("  tab03_metricas_por_grupo.png\n")
cat("  fig02_sensibilidade_por_grupo.png\n")
cat("  fig03_matrizes_confusao.png\n")
cat(strrep("=", 70), "\n", sep = "")
