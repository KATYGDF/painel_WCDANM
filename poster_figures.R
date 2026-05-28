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
#   ggsave() renders at the exact W × H specified — text scales correctly.
#
# Run from project root (C:/Katy/Doutorado/painel_WCDANM):
#   source("poster_figures.R")
#
# One-time dependency for fig2:
#   install.packages("webshot2")
# =============================================================================


# ── DIMENSIONS — edit here, then re-source ───────────────────────────────────
#
# HOW TO GET SIZES FROM CANVA:
#   Click the white box → Avançados → Largura / Altura (in mm)
#   Paste the values below. The script converts mm → inches automatically.
#
# White content box (from Canva "Avançados"):
CANVA_W_MM <- 353.06    # Largura
CANVA_H_MM <- 200.68    # Altura

DPI <- 300              # 300 dpi for crisp text at poster print size

# Fraction of box height for each figure (adjust until preview looks right;
# fig1 + fig2 should sum to ~1.0 if they share the same box):
FRAC_FIG1  <- 0.42      # EDA distributions   (upper strip)   — height fraction
FRAC_FIG2  <- 0.55      # Results table       (lower portion) — height fraction
FRAC_FIG3  <- 0.50      # Confusion matrices  (half height)   — height fraction
FRAC_FIG4  <- 1.00      # PCA by group        (full height)   — height fraction

W_FIG2     <- 0.50      # ← table width as fraction of box width (reduce to shrink)

# Methodology panel (fig0) — left column dimensions in mm
W0_MM <- 160            # ← width of the methodology column in Canva (measure it)
H0_MM <- 200            # ← height of the methodology column in Canva

# Derived sizes in inches (do not edit below this line)
BOX_W <- CANVA_W_MM / 25.4
BOX_H <- CANVA_H_MM / 25.4

W1 <- BOX_W;              H1 <- BOX_H * FRAC_FIG1
W2 <- BOX_W * W_FIG2;     H2 <- BOX_H * FRAC_FIG2
W3 <- BOX_W;              H3 <- BOX_H * FRAC_FIG3
W4 <- BOX_W*0.4;          H4 <- BOX_H * FRAC_FIG4

cat(sprintf(paste0(
  "Box: %.2f × %.2f in (%.0f × %.0f mm)\n",
  "  fig1 (EDA):     %.2f × %.2f in\n",
  "  fig2 (table):   %.2f × %.2f in\n",
  "  fig3 (conf.):   %.2f × %.2f in\n",
  "  fig4 (PCA):     %.2f × %.2f in\n\n"),
  BOX_W, BOX_H, CANVA_W_MM, CANVA_H_MM,
  W1, H1, W2, H2, W3, H3, W4, H4))

# ── 0. Packages & directories ─────────────────────────────────────────────────

pkgs  <- c("dplyr", "tidyr", "ggplot2", "patchwork", "scales",
           "clue", "aricode", "gt", "webshot2", "ggtext")
novos <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(novos) > 0) {
  message("Installing: ", paste(novos, collapse = ", "))
  install.packages(setdiff(novos, "webshot2"),
                   repos = "https://cloud.r-project.org", quiet = TRUE)
  if ("webshot2" %in% novos)
    install.packages("webshot2", repos = "https://cloud.r-project.org", quiet = TRUE)
}

# Confirm critical packages loaded — stop early with a clear message if missing
for (pkg in c("ggtext", "patchwork", "gt")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop(sprintf("Package '%s' could not be installed. Run: install.packages('%s')", pkg, pkg))
}

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
  library(patchwork); library(scales)
  library(clue); library(aricode)
  library(gt); library(ggtext)
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
         "G3" = "#E8B96B",   # Acute
         "G4" = "#D98585")   # Atypical

G_LABELS <- c("G1 · Low utilization",
              "G2 · Coordinated outpatient",
              "G3 · Acute",
              "G4 · Atypical")

# Shared legend string — used in captions / footnotes of all figures
LEGENDA_EN <- paste0(
  "G1 = Low utilization  ·  G2 = Coordinated outpatient  ·  ",
  "G3 = Acute  ·  G4 = Atypical"
)

# English names for the 5 count variables
VARS_V3 <- c("consultas", "ps", "exames", "internacoes", "terapias")
VAR_EN  <- c(consultas   = "medical consultations",
             ps          = "emergency visits",
             exames      = "diagnostic exams",
             internacoes = "hospitalizations",
             terapias    = "therapies")

# Figure / table captions — numbered, bold label + normal text
# Format: "**Figure X.** Description. Legend."
CAP1 <- paste0(
  "**Figure 1.** Distribution of count variables by latent group (n = 5,000). ",
  "Grey bars: overall histogram (density scale). Coloured curves: kernel density per true group; ",
  "x-axis truncated at the 97th percentile.\n"
  
)
CAP2 <- paste0(
  "**Table 1.** Latent group recovery metrics for three variable sets and two methods. ",
  "*g* selected by BIC (NB Mixture) or Silhouette (K-means). ",
  "★ Best overall result."
)
CAP3 <- paste0(
  "**Figure 2.** Confusion matrices for the best variable set (V3: 5 variables), ",
  "NB Mixture (left) vs. K-means (right). ",
  "Diagonal = correct assignments. Labels paired via Hungarian algorithm.\n"
  
)
CAP4 <- paste0(
  "**Figure 3.** PCA projection (log1p + scaled) of V3 variables. ",
  "Each row highlights one true group (coloured); remaining points in grey. ",
  "Columns: ground truth | NB Mixture estimate | K-means estimate.\n"
  
)

# Caption element: element_textbox_simple wraps text automatically AND
# renders **bold** markdown — solves captions being cut off at figure edge
cap_elem <- ggtext::element_textbox_simple(
  colour    = "grey20",
  size      = 10,
  hjust     = 0.5,        # centre-align the box itself
  halign    = 0.5,        # centre-align text inside the box
  margin    = margin(t = 8),
  padding   = margin(0, 4, 0, 4)
)

# Poster theme
TEMA_P <- theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor  = element_blank(),
    strip.text        = element_text(face = "bold", size = 11),
    strip.background  = element_rect(fill = "#f2f2f2", colour = NA),
    plot.title        = element_text(face = "bold", size = 14),
    plot.subtitle     = element_text(colour = "grey40", size = 10),
    plot.caption      = cap_elem,
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
    caption  = CAP1
  ) +
  guides(
    fill   = guide_legend(nrow = 1, override.aes = list(alpha = 0.55, linewidth = 0)),
    colour = "none"
  ) +
  TEMA_P

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
    "medical consultations · diagnostic exams · emergency visits", "",
    "V1  +  hospitalizations",                          "",
    "V2  +  therapies",                          ""
  ),
  Method   = c("NB Mixture", "K-means",
               "NB Mixture", "K-means",
               "NB Mixture ★", "K-means"),   # ★
  g        = c(3L, 3L, 3L, 4L, 4L, 4L),
  ARI      = c(0.630, 0.570, 0.658, 0.580, 0.729, 0.673),
  Accuracy = c(0.820, 0.797, 0.836, 0.777, 0.877, 0.815),
  stringsAsFactors = FALSE
)

gt_tbl <- gt(tab_data) %>%
  tab_header(
    title    = md("**Table 1.** Latent Group Recovery — Variable Sets × Methods"),
    subtitle = md("*g* selected by BIC (NB Mixture) and Silhouette (K-means) &nbsp;·&nbsp; n = 5 000")
  ) %>%
  cols_label(
    Set = "Set", Variables = "Variables included",
    Method = "Method", g = md("*g*"),
    ARI = "ARI", Accuracy = "Accuracy"
  ) %>%
  fmt_number(columns = c(ARI, Accuracy), decimals = 3) %>%
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
    locations = list(cells_body(columns = c(g, ARI, Accuracy)),
                     cells_column_labels(columns = c(g, ARI, Accuracy)))
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
  tab_source_note(
    source_note = md(paste0("**Legend:** ", LEGENDA_EN))
  ) %>%
  tab_style(
    style     = cell_text(color = "white", weight = "bold"),
    locations = cells_column_labels(everything())
  ) %>%
  # Proportional column widths — fills viewport without leaving dead whitespace
  cols_width(
    Set       ~ pct(6),
    Variables ~ pct(40),
    Method    ~ pct(20),
    g         ~ pct(7),
    ARI       ~ pct(13),
    Accuracy  ~ pct(14)
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

# vwidth/vheight: viewport in px (zoom=2 doubles it, giving W2*DPI output pixels)
tryCatch({
  gtsave(gt_tbl, file.path(DIR_OUT, "fig2_results_table.png"),
         vwidth  = round(W2 * DPI / 2),
         vheight = round(H2 * DPI / 2),
         zoom = 2, expand = 0)
  cat("  ✓ fig2_results_table.png  (", round(W2, 2), "×", round(H2, 2), "in ·", DPI, "dpi)\n\n")
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
    caption = CAP3,
    theme   = TEMA_P + theme(plot.title = element_text(face = "bold", size = 15))
  )

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
      values = c(setNames(g_color, g_lbl), other = "grey60"),
      guide  = "none"
    ) +
    scale_alpha_manual(
      values = c(setNames(0.75, g_lbl), other = 0.25),
      guide  = "none"
    ) +
    scale_size_manual(
      values = c(setNames(0.6, g_lbl), other = 0.15),
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
    title   = "PCA Classification",
    caption = CAP4,
    theme   = TEMA_P + theme(plot.title = element_text(face = "bold", size = 12))
  )

ggsave(file.path(DIR_OUT, "fig4_pca_by_group.png"),
       fig4, width = W4, height = H4, dpi = DPI, bg = "white")
cat("  ✓ fig4_pca_by_group.png  (", W4, "×", H4, "in ·", DPI, "dpi)\n\n")


# =============================================================================
# FIG 0 · Methodology panel (HTML → PNG via webshot2)
# =============================================================================

cat("Generating fig0_methodology.png ...\n")

meth_html <- sprintf('<!DOCTYPE html>
<html>
<head>
<meta charset="UTF-8">
<style>
  * { box-sizing: border-box; margin: 0; padding: 0; }
  body {
    font-family: "Segoe UI", Arial, sans-serif;
    font-size: 15px; line-height: 1.5;
    width: %dpx; background: white;
  }
  .header {
    background: #1a4a72; color: white;
    font-weight: bold; font-size: 18px;
    padding: 10px 16px; letter-spacing: 0.05em;
    text-transform: uppercase;
  }
  .body { padding: 14px 16px; }
  .block { margin-bottom: 11px; }
  .label { font-weight: bold; }
  .method-title { font-weight: bold; margin-bottom: 4px; }
  .formula-wrap {
    text-align: center; margin: 8px 0 4px;
    font-family: "Cambria Math", "Georgia", serif;
    font-size: 15px;
  }
  .where {
    font-size: 13px; color: #333; margin: 4px 0 0 8px;
    line-height: 1.65;
  }
  sub, sup { font-size: 0.72em; }
  i { font-style: italic; }
  hr { border: none; border-top: 1px solid #ddd; margin: 10px 0; }
</style>
</head>
<body>
<div class="header">Methodology</div>
<div class="body">

  <div class="block">
    <span class="label">Data:</span>
    A synthetic dataset of 5,000 beneficiaries structured into four latent
    groups and described by five count variables: <i>medical consultations,
    emergency visits, diagnostic exams, hospitalizations</i> and <i>therapies.</i>
    Marginal distributions calibrated to ANS panel indicators.
  </div>

  <div class="block">
    <span class="label">Latent profiles:</span>
    G1 &middot; Low utilization &nbsp;&nbsp;
    G2 &middot; Coordinated outpatient &nbsp;&nbsp;
    G3 &middot; Acute &nbsp;&nbsp;
    G4 &middot; Atypical pattern
  </div>

  <hr>

  <div class="block">
    <div class="method-title">M1 &mdash; Negative Binomial Mixture.</div>
    Each observation <b>x</b><sub><i>i</i></sub> = (<i>x</i><sub><i>i</i>1</sub>, &hellip;,
    <i>x</i><sub><i>ip</i></sub>) &isin; &#8469;<sup><i>p</i></sup>
    is a multivariate count vector with <i>p</i> &isin; {3, 4, 5} variables.
    Under conditional independence of the <i>p</i> variables given the latent component,
    the component-specific density factors as a product of <i>p</i> independent NB marginals.
    Parameters estimated by maximum likelihood via the EM algorithm;
    number of components <i>g</i> &isin; {1, &hellip;, 4} selected by BIC.
    <div class="formula-wrap">
      <i>f</i>(<b>x</b><sub><i>i</i></sub> | &Theta;) =
      &nbsp;<span style="font-size:1.3em">&sum;</span><sub style="font-size:0.65em"><i>k</i>=1</sub><sup style="font-size:0.65em"><i>g</i></sup>&nbsp;
      &pi;<sub><i>k</i></sub>
      &nbsp;<span style="font-size:1.3em">&prod;</span><sub style="font-size:0.65em"><i>j</i>=1</sub><sup style="font-size:0.65em"><i>p</i></sup>&nbsp;
      NB(<i>x</i><sub><i>ij</i></sub> | &mu;<sub><i>kj</i></sub>, &phi;<sub><i>kj</i></sub>)
    </div>
    <div class="where">
      &pi;<sub><i>k</i></sub> is the weight of component <i>k</i>, with
        &sum;<sub><i>k</i></sub>&pi;<sub><i>k</i></sub> = 1 and &pi;<sub><i>k</i></sub> &ge; 0<br>
      &mu;<sub><i>kj</i></sub> is the NB mean of variable <i>j</i> in component <i>k</i><br>
      &phi;<sub><i>kj</i></sub> is the corresponding dispersion parameter<br>
      &Theta; = {&pi;<sub><i>k</i></sub>, &mu;<sub><i>kj</i></sub>, &phi;<sub><i>kj</i></sub>}
        denotes the full set of parameters
    </div>
  </div>

  <hr>

  <div class="block">
    <div class="method-title">M2 &mdash; K-means (baseline).</div>
    Applied to log-transformed and standardised counts.
    Number of clusters selected by the Silhouette criterion.
    <div class="formula-wrap">
      arg&nbsp;min<sub>{<i>C</i><sub>1</sub>,&hellip;,<i>C</i><sub><i>k</i></sub>}</sub>
      &nbsp;<span style="font-size:1.3em">&sum;</span><sub style="font-size:0.65em"><i>j</i>=1</sub><sup style="font-size:0.65em"><i>k</i></sup>&nbsp;
      <span style="font-size:1.3em">&sum;</span><sub style="font-size:0.65em"><i>x</i><sub><i>i</i></sub>&isin;<i>C</i><sub><i>j</i></sub></sub>
      &nbsp;&#8214;<i>x</i><sub><i>i</i></sub> &minus; &mu;<sub><i>j</i></sub>&#8214;<sup>2</sup>
    </div>
    <div class="where">
      <i>C</i><sub>1</sub>, &hellip;, <i>C</i><sub><i>k</i></sub>
        are the <i>k</i> partitions (clusters)<br>
      &mu;<sub><i>j</i></sub> = <sup>1</sup>&frasl;<sub>|<i>C</i><sub><i>j</i></sub>|</sub>
        &sum;<sub><i>x</i><sub><i>i</i></sub>&isin;<i>C</i><sub><i>j</i></sub></sub>
        <i>x</i><sub><i>i</i></sub> is the centroid of cluster <i>j</i><br>
      &nbsp;&#8214;&nbsp;&middot;&nbsp;&#8214;<sup>2</sup>
        is the squared Euclidean distance
    </div>
  </div>

</div>
</body>
</html>', round(W0_MM / 25.4 * DPI / 2))

tmp_meth <- tempfile(fileext = ".html")
writeLines(meth_html, tmp_meth, useBytes = FALSE)

tryCatch({
  webshot2::webshot(
    tmp_meth,
    file.path(DIR_OUT, "fig0_methodology.png"),
    vwidth  = round(W0_MM / 25.4 * DPI / 2),
    vheight = round(H0_MM / 25.4 * DPI / 2),
    zoom    = 2,
    delay   = 0.2
  )
  cat("  ✓ fig0_methodology.png  (", W0_MM, "×", H0_MM, "mm)\n\n")
}, error = function(e) {
  message("  ✗ fig0 failed: ", conditionMessage(e))
})

# ── Summary ───────────────────────────────────────────────────────────────────

cat(strrep("─", 60), "\n")
cat("Done. Figures saved to:", DIR_OUT, "\n\n")
cat("To change size: edit W1/H1, W3/H3, W4/H4 at the top of the\n")
cat("script and re-source — text will scale proportionally.\n")
cat(strrep("─", 60), "\n")
