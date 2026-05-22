# =============================================================================
# poster_figures.R
# Poster-ready figures — all labels in English · G1–G4 notation
#
# Generates (in poster/figs/):
#   fig1_eda_distributions.png  — histogram + density for 5 count variables
#   fig2_results_table.png      — comparative results table (V1–V3 × 2 methods)
#   fig3_confusion_matrix.png   — best model confusion matrix (V3 · NB Mixture)
#   fig4_pca_by_group.png       — PCA projection by group (4 rows × 3 columns)
#
# HOW TO ADJUST FIGURE SIZE:
#   Edit the W/H variables in the DIMENSIONS block below, then re-source.
#   windows() opens a preview — text and elements scale correctly because
#   ggsave() re-renders at the exact W × H you specify (300 dpi).
#
# Run from project root (C:/Katy/Doutorado/painel_WCDANM):
#   source("poster_figures.R")
#
# One-time dependency for fig2:
#   install.packages("webshot2")
# =============================================================================


# ── DIMENSIONS (inches) — edit here, then re-source ──────────────────────────
#   Rule of thumb for poster: width ≈ column width in inches, height as needed.

W1 <- 20;   H1 <- 5.5   # fig1 · EDA distributions  (wide strip)
W3 <- 12;   H3 <- 6     # fig3 · Confusion matrices  (side-by-side)
W4 <- 15;   H4 <- 17    # fig4 · PCA by group        (tall 4-row panel)
DPI <- 300              # resolution for all saved files

# ── 0. Packages & directories ─────────────────────────────────────────────────

pkgs  <- c("dplyr", "tidyr", "ggplot2", "patchwork", "scales",
           "clue", "aricode", "gt", "webshot2")
novos <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(novos) > 0) {
  message("Installing: ", paste(novos, collapse = ", "))
  install.packages(setdiff(novos, "webshot2"),
                   repos = "https://cloud.r-project.org", quiet = TRUE)
  if ("webshot2" %in% novos)
    install.packages("webshot2", repos = "https://cloud.r-project.org", quiet = TRUE)
}

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
  library(patchwork); library(scales)
  library(clue); library(aricode)
  library(gt)
})

DIR_OUT  <- "poster/figs"
RDS_DIR  <- "src/03_analises/rds"
CSV_DATA <- "dados/dados_simulados.csv"

if (!dir.exists(DIR_OUT)) dir.create(DIR_OUT, recursive = TRUE)

cat(strrep("─", 60), "\n")
cat("poster_figures.R — generating", DIR_OUT, "\n")
cat(strrep("─", 60), "\n\n")


# ── 1. Load data & saved model results ───────────────────────────────────────

df_raw <- read.csv(CSV_DATA, stringsAsFactors = FALSE)
col_g  <- c("grupo_latente", "grupo", "grupo_real", "grupo_idx")
col_g  <- col_g[col_g %in% names(df_raw)][1]
if (is.na(col_g)) stop("Group column not found in CSV.")

df <- df_raw
names(df)[names(df) == col_g] <- "grupo"
df$grupo <- factor(df$grupo,
                   levels = c("baixo_uso", "ambulatorial_coordenado",
                              "agudo_hospitalar", "atipico"))
grupo_int <- as.integer(df$grupo)
n         <- nrow(df)

fits_nb    <- readRDS(file.path(RDS_DIR, "nb_fits_por_conjunto.rds"))
km_results <- readRDS(file.path(RDS_DIR, "km_results_por_conjunto.rds"))

cat(sprintf("Data loaded: %d obs · 4 true groups\n\n", n))


# ── 2. Aesthetic constants ────────────────────────────────────────────────────

# Colour palette: G1=blue · G2=green · G3=orange · G4=pink
# Order in factor: 1=baixo_uso 2=ambulatorial 3=agudo 4=atipico
COR <- c("G1" = "#7BAFD4",   # Low Utilization
         "G2" = "#6EBF8B",   # Coordinated Outpatient
         "G3" = "#E8B96B",   # Acute / Hospital
         "G4" = "#D98585")   # Atypical

G_LABELS <- c("G1 · Low Utilization",
              "G2 · Coordinated Outpatient",
              "G3 · Acute / Hospital",
              "G4 · Atypical")

# English names for the 5 count variables
VARS_V3 <- c("consultas", "ps", "exames", "internacoes", "terapias")
VAR_EN  <- c(consultas   = "Outpatient Visits",
             ps          = "ER Visits",
             exames      = "Diagnostic Exams",
             internacoes = "Hospitalizations",
             terapias    = "Therapy Sessions")

# Poster theme — slightly larger than default for print legibility
TEMA_P <- theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor  = element_blank(),
    strip.text        = element_text(face = "bold", size = 11),
    strip.background  = element_rect(fill = "#f2f2f2", colour = NA),
    plot.title        = element_text(face = "bold", size = 14),
    plot.subtitle     = element_text(colour = "grey40", size = 10),
    plot.caption      = element_text(colour = "grey55", size = 9, hjust = 0),
    legend.position   = "bottom",
    legend.title      = element_text(face = "bold", size = 10),
    legend.text       = element_text(size = 10)
  )


# ── 3. Helper functions ───────────────────────────────────────────────────────

# Hungarian optimal pairing: maps estimated labels → true group numbering
parear_otimo <- function(real, estimado) {
  g     <- max(max(real), max(estimado))
  custo <- matrix(0, g, g)
  for (i in seq_len(max(real)))
    for (j in seq_len(max(estimado)))
      custo[i, j] <- -sum(real == i & estimado == j)
  sol <- clue::solve_LSAP(custo + abs(min(custo)))
  order(as.integer(sol))[estimado]
}

# Confusion matrix data frame with Gx labels
conf_df_en <- function(real, estimado) {
  est_par <- parear_otimo(real, estimado)
  g_lbl   <- paste0("G", 1:4)
  ct <- as.data.frame(table(
    True     = factor(paste0("G", real),     levels = g_lbl),
    Estimate = factor(paste0("G", est_par),  levels = g_lbl)
  ))
  ct$pct      <- ct$Freq / n * 100
  ct$True_rev <- factor(ct$True, levels = rev(g_lbl))
  ct
}

# Single confusion matrix ggplot
plot_conf_en <- function(ct, titulo, subtitulo) {
  ggplot(ct, aes(Estimate, True_rev)) +
    geom_tile(aes(fill = Freq), colour = "white", linewidth = 0.9) +
    geom_text(
      aes(label = ifelse(Freq > 0,
                         paste0(Freq, "\n(", sprintf("%.1f%%", pct), ")"), "")),
      size = 3.4, lineheight = 1.1, colour = "grey10"
    ) +
    scale_fill_gradient(low = "#eef5fb", high = "#2980b9", name = "n") +
    scale_x_discrete(position = "top") +
    labs(x = "Estimated group", y = "True group",
         title = titulo, subtitle = subtitulo) +
    TEMA_P +
    theme(panel.grid    = element_blank(),
          legend.position = "none",
          axis.text     = element_text(size = 11, face = "bold"))
}


# =============================================================================
# FIG 1 · EDA — Count distributions by latent group (5 variables)
# =============================================================================

cat("Generating fig1_eda_distributions.png ...\n")

df_long <- df %>%
  select(grupo, all_of(VARS_V3)) %>%
  pivot_longer(-grupo, names_to = "variable", values_to = "count") %>%
  mutate(
    grupo_lbl = paste0("G", as.integer(grupo)),
    variable  = factor(VAR_EN[variable], levels = VAR_EN)
  )

# Truncate at P97 per variable (same criterion as analysis)
df_long <- df_long %>%
  group_by(variable) %>%
  filter(count <= quantile(count, 0.97)) %>%
  ungroup()

fig1 <- ggplot(df_long, aes(x = count)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 40, fill = "grey82", colour = "white", alpha = 0.8) +
  geom_density(aes(fill = grupo_lbl, colour = grupo_lbl),
               alpha = 0.28, linewidth = 0.8) +
  scale_fill_manual(values = COR, name = "True group", labels = G_LABELS) +
  scale_colour_manual(values = COR, name = "True group", labels = G_LABELS) +
  facet_wrap(~ variable, scales = "free", ncol = 5) +
  labs(
    x        = "Count",
    y        = "Density",
    title    = "Count Variable Distributions by Latent Group",
    subtitle = sprintf(
      "n = %d  ·  x-axis truncated at P97 per variable  ·  density curves per true group",
      n),
    caption  = "True groups known from simulation (ground truth). Mixture structure visible in all variables."
  ) +
  guides(
    fill   = guide_legend(nrow = 1, override.aes = list(alpha = 0.55, linewidth = 0)),
    colour = "none"
  ) +
  TEMA_P

windows(width = W1, height = H1)   # preview — resize to taste, then adjust W1/H1 above
print(fig1)
ggsave(file.path(DIR_OUT, "fig1_eda_distributions.png"),
       fig1, width = W1, height = H1, dpi = DPI, bg = "white")
cat("  ✓ fig1_eda_distributions.png  (", W1, "×", H1, "in ·", DPI, "dpi)\n\n")


# =============================================================================
# FIG 2 · Results table (saved via gt + webshot2)
# =============================================================================

cat("Generating fig2_results_table.png ...\n")

tab_data <- data.frame(
  Set = c("V1", "", "V2", "", "V3", ""),
  Variables = c(
    "Outpatient Visits · Diagnostic Exams · ER Visits", "",
    "V1  +  Hospitalizations",                          "",
    "V2  +  Therapy Sessions",                          ""
  ),
  Method   = c("NB Mixture", "K-means",
               "NB Mixture", "K-means",
               "NB Mixture ★", "K-means"),   # ★
  g        = c(3L, 3L, 3L, 4L, 4L, 4L),
  ARI      = c(0.630, 0.570, 0.658, 0.580, 0.729, 0.673),
  Accuracy = c(0.820, 0.797, 0.836, 0.777, 0.877, 0.815),
  Purity   = c(0.796, 0.768, 0.814, 0.803, 0.861, 0.843),
  stringsAsFactors = FALSE
)

gt_tbl <- gt(tab_data) %>%
  tab_header(
    title    = md("**Latent Group Recovery — Variable Sets × Methods**"),
    subtitle = md("*g* selected by BIC (NB Mixture) and Silhouette (K-means) &nbsp;·&nbsp; n = 5 000")
  ) %>%
  cols_label(
    Set = "Set", Variables = "Variables included",
    Method = "Method", g = md("*g*"),
    ARI = "ARI", Accuracy = "Accuracy", Purity = "Purity"
  ) %>%
  fmt_number(columns = c(ARI, Accuracy, Purity), decimals = 3) %>%
  # Highlight best row
  tab_style(
    style     = list(cell_fill(color = "#fdf6e3"),
                     cell_text(weight = "bold")),
    locations = cells_body(rows = 5)
  ) %>%
  # Highlight V3 set label
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_body(columns = Set, rows = 5)
  ) %>%
  # Center numeric columns
  tab_style(
    style     = cell_text(align = "center"),
    locations = list(cells_body(columns = c(g, ARI, Accuracy, Purity)),
                     cells_column_labels(columns = c(g, ARI, Accuracy, Purity)))
  ) %>%
  # Dim empty Variable/Set cells
  tab_style(
    style     = cell_text(color = "grey70"),
    locations = cells_body(columns = c(Set, Variables), rows = c(2, 4, 6))
  ) %>%
  tab_footnote(
    footnote  = md("★ Best overall — BIC selected *g* = 4, matching the true number of latent groups."),
    locations = cells_body(columns = Method, rows = 5)
  ) %>%
  tab_footnote(
    footnote  = "Groups: G1 = Low Utilization · G2 = Coordinated Outpatient · G3 = Acute/Hospital · G4 = Atypical",
    locations = cells_column_labels(columns = Variables)
  ) %>%
  tab_style(
    style     = cell_text(color = "white", weight = "bold"),
    locations = cells_column_labels(everything())
  ) %>%
  tab_options(
    table.font.size                = px(14),
    table.width                    = pct(100),
    heading.align                  = "left",
    heading.title.font.size        = px(16),
    heading.subtitle.font.size     = px(13),
    column_labels.background.color = "#2c3e50",
    column_labels.border.top.style = "hidden",
    table_body.hlines.color        = "#e0e0e0",
    stub.background.color          = "white",
    source_notes.font.size         = px(11)
  )

tryCatch({
  gtsave(gt_tbl, file.path(DIR_OUT, "fig2_results_table.png"),
         zoom = 2, expand = 15)
  cat("  ✓ fig2_results_table.png\n\n")
}, error = function(e) {
  message("  ✗ fig2 failed: ", conditionMessage(e))
  message("    Run once: install.packages('webshot2')")
  message("    Then rerun this script.\n")
})


# =============================================================================
# FIG 3 · Confusion matrix — V3, NB Mixture (best) vs K-means
# =============================================================================

cat("Generating fig3_confusion_matrix.png ...\n")

nb_v3 <- fits_nb[["V3"]]
km_v3 <- km_results[["V3"]]

ct_nb <- conf_df_en(grupo_int, nb_v3$classification)
ct_km <- conf_df_en(grupo_int, km_v3$cluster)

ari_nb  <- round(aricode::ARI(grupo_int, nb_v3$classification), 3)
ari_km  <- round(aricode::ARI(grupo_int, km_v3$cluster),        3)
acc_nb  <- round(mean(parear_otimo(grupo_int, nb_v3$classification) == grupo_int), 3)
acc_km  <- round(mean(parear_otimo(grupo_int, km_v3$cluster)        == grupo_int), 3)

p_nb <- plot_conf_en(
  ct_nb,
  titulo    = "NB Mixture  (V3)",
  subtitulo = sprintf("ARI = %.3f  ·  Accuracy = %.3f", ari_nb, acc_nb)
)

p_km <- plot_conf_en(
  ct_km,
  titulo    = "K-means  (V3)",
  subtitulo = sprintf("ARI = %.3f  ·  Accuracy = %.3f", ari_km, acc_km)
)

fig3 <- (p_nb | p_km) +
  plot_annotation(
    title   = "Confusion Matrices — Best Variable Set (V3: 5 variables)",
    caption = paste0(
      "n = ", n, "  ·  4 true latent groups  ·  ",
      "Labels paired via Hungarian (Kuhn-Munkres) algorithm\n",
      "G1 = Low Utilization  ·  G2 = Coordinated Outpatient  ·  ",
      "G3 = Acute/Hospital  ·  G4 = Atypical"
    ),
    theme = TEMA_P + theme(plot.title = element_text(face = "bold", size = 15))
  )

windows(width = W3, height = H3)   # preview — resize to taste, then adjust W3/H3 above
print(fig3)
ggsave(file.path(DIR_OUT, "fig3_confusion_matrix.png"),
       fig3, width = W3, height = H3, dpi = DPI, bg = "white")
cat("  ✓ fig3_confusion_matrix.png  (", W3, "×", H3, "in ·", DPI, "dpi)\n\n")


# =============================================================================
# FIG 4 · PCA projection by group — 4 rows (G1–G4) × 3 columns (True | NB | KM)
# =============================================================================

cat("Generating fig4_pca_by_group.png ...\n")

# PCA on V3 variables: log1p + scale
Y_v3 <- as.matrix(df[, VARS_V3])
pca  <- prcomp(scale(log1p(Y_v3)))
pct1 <- round(summary(pca)$importance[2, 1] * 100, 1)
pct2 <- round(summary(pca)$importance[2, 2] * 100, 1)

nb_paired <- parear_otimo(grupo_int, nb_v3$classification)
km_paired <- parear_otimo(grupo_int, km_v3$cluster)

df_pca <- data.frame(
  PC1  = pca$x[, 1],
  PC2  = pca$x[, 2],
  True = paste0("G", grupo_int),
  NB   = paste0("G", nb_paired),
  KM   = paste0("G", km_paired)
)

# Long: one row per observation × panel type
df_pca_long <- bind_rows(
  df_pca %>% mutate(panel = "True Group",  shown = True),
  df_pca %>% mutate(panel = "NB Mixture",  shown = NB),
  df_pca %>% mutate(panel = "K-means",     shown = KM)
) %>%
  mutate(panel = factor(panel, levels = c("True Group", "NB Mixture", "K-means")))

# Build one row (3-panel strip) for a given group g ∈ {1,2,3,4}
plot_pca_row <- function(g_idx) {
  g_lbl   <- paste0("G", g_idx)
  g_color <- COR[[g_lbl]]

  df_pca_long %>%
    mutate(
      highlight = (shown == g_lbl),
      col_val   = ifelse(highlight, g_lbl, "other")
    ) %>%
    arrange(highlight) %>%                     # grey points rendered first
    ggplot(aes(PC1, PC2, colour = col_val,
               alpha = col_val, size = col_val)) +
    geom_point(shape = 16) +
    scale_colour_manual(
      values = c(setNames(g_color, g_lbl), other = "grey84"),
      guide  = "none"
    ) +
    scale_alpha_manual(
      values = c(setNames(0.80, g_lbl), other = 0.12),
      guide  = "none"
    ) +
    scale_size_manual(
      values = c(setNames(1.0, g_lbl), other = 0.25),
      guide  = "none"
    ) +
    facet_wrap(~ panel, nrow = 1) +
    labs(
      x        = sprintf("PC1 (%.1f%% var.)", pct1),
      y        = sprintf("PC2 (%.1f%% var.)", pct2),
      subtitle = G_LABELS[g_idx]
    ) +
    TEMA_P +
    theme(
      strip.text    = element_text(face = "bold", size = 10),
      plot.subtitle = element_text(face = "bold", colour = g_color, size = 11),
      axis.title    = element_text(size = 8),
      axis.text     = element_text(size = 7),
      plot.margin   = margin(3, 6, 3, 6)
    )
}

fig4 <- wrap_plots(lapply(1:4, plot_pca_row), ncol = 1) +
  plot_annotation(
    title   = "PCA Classification by Group — V3 (log1p + scaled)",
    caption = paste0(
      "Each row highlights one true group (coloured); remaining points in grey.\n",
      "Columns: ground truth | NB Mixture estimate | K-means estimate.\n",
      "V3: Outpatient Visits · ER Visits · Diagnostic Exams · ",
      "Hospitalizations · Therapy Sessions  ·  n = ", n
    ),
    theme = TEMA_P + theme(plot.title = element_text(face = "bold", size = 15))
  )

windows(width = W4, height = H4)   # preview — resize to taste, then adjust W4/H4 above
print(fig4)
ggsave(file.path(DIR_OUT, "fig4_pca_by_group.png"),
       fig4, width = W4, height = H4, dpi = DPI, bg = "white")
cat("  ✓ fig4_pca_by_group.png  (", W4, "×", H4, "in ·", DPI, "dpi)\n\n")


# ── Summary ───────────────────────────────────────────────────────────────────

cat(strrep("─", 60), "\n")
cat("Done. Figures saved to:", DIR_OUT, "\n\n")
cat("To change size: edit W1/H1, W3/H3, W4/H4 at the top of the\n")
cat("script and re-source — text will scale proportionally.\n")
cat(strrep("─", 60), "\n")
