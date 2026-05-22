# ============================================================
# RECUPERAR OS 4 GRUPOS COM NB - 3 ESTRATEGIAS
#   E1: NB padrao com inicializacao por quantis (theta livre)
#   E2: NB com theta fixo no valor verdadeiro
#   E3: NB com theta fixo + inicializacao por quantis
# Comparados ao modelo aleatorio padrao (baseline) e ao oraculo.
# ============================================================

library(dplyr)
library(ggplot2)
library(flexmix)
if (!requireNamespace("countreg", quietly = TRUE)) {
  stop("Instale: install.packages('countreg', repos='http://R-Forge.R-project.org')")
}
FLXMRnegbin <- countreg::FLXMRnegbin

set.seed(1234)

# ------------------------------------------------------------
# 1) SIMULACAO
# ------------------------------------------------------------
n         <- 5000
pi_grupos <- c(0.25, 0.25, 0.25, 0.25)
mu_grupos <- c(2, 12, 25, 50)
theta_nb  <- 50

rotulos_grupo <- c("baixo_uso", "ambulatorial_coordenado",
                   "agudo_hospitalar", "atipico")

grupo_real <- sample(1:4, size = n, replace = TRUE, prob = pi_grupos)
consultas  <- rnbinom(n = n, size = theta_nb, mu = mu_grupos[grupo_real])

df <- data.frame(
  grupo_real = factor(grupo_real, levels = 1:4, labels = rotulos_grupo),
  consultas  = consultas
)

# ------------------------------------------------------------
# 2) Formato longo para o ggplot (mesmo padrao do script original)
# ------------------------------------------------------------
df_long <- df %>%
  mutate(grupo = grupo_real) %>%
  transmute(grupo,
            variavel = "consultas",
            valor    = consultas)

# Truncar no P97 (igual ao seu codigo)
p97 <- quantile(df_long$valor, 0.97)
df_long <- df_long %>% filter(valor <= p97)

# ------------------------------------------------------------
# 3) Estetica (suas cores e tema)
# ------------------------------------------------------------
COR_GRUPO <- c(
  agudo_hospitalar        = "#c0392b",
  ambulatorial_coordenado = "#2980b9",
  atipico                 = "#8e44ad",
  baixo_uso               = "#4a7c59"
)

TEMA <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor  = element_blank(),
    strip.text        = element_text(face = "bold", size = 9),
    strip.background  = element_rect(fill = "#f0f0f0", colour = NA),
    plot.title        = element_text(face = "bold", size = 12),
    plot.subtitle     = element_text(colour = "grey40", size = 9),
    plot.caption      = element_text(colour = "grey55", size = 8, hjust = 0),
    legend.position   = "bottom"
  )

# ------------------------------------------------------------
# 4) Grafico
# ------------------------------------------------------------
fig08 <- ggplot(df_long, aes(x = valor)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 40, fill = "grey78", colour = "white", alpha = 0.7) +
  geom_density(aes(fill = grupo, colour = grupo),
               alpha = 0.25, linewidth = 0.7) +
  scale_fill_manual(values = COR_GRUPO, name = "Grupo") +
  scale_colour_manual(values = COR_GRUPO, name = "Grupo") +
  facet_wrap(~ variavel, scales = "free", ncol = 3) +
  labs(
    x        = "Contagem",
    y        = "Densidade",
    title    = "Distribuicao de consultas por grupo latente",
    subtitle = sprintf("n = %d  ·  eixo x truncado no P97  ·  densidades sobrepostas por grupo", n),
    caption  = "Grupos verdadeiros (gabarito simulado)"
  ) +
  TEMA

print(fig08)
# ------------------------------------------------------------
# 2) LOGLIK-ORACULO (teto teorico para qualquer modelo k=4)
# ------------------------------------------------------------
# log p(y) = log sum_h pi_h * NB(y | mu_h, theta)
ll_obs <- function(y, mus, thetas, pis) {
  # mus, thetas, pis sao vetores de tamanho k
  k <- length(mus)
  lp <- sapply(1:k, function(h) {
    log(pis[h]) + dnbinom(y, size = thetas[h], mu = mus[h], log = TRUE)
  })
  # log-sum-exp por linha
  m <- apply(lp, 1, max)
  m + log(rowSums(exp(lp - m)))
}

ll_oraculo <- sum(ll_obs(consultas,
                         mus    = mu_grupos,
                         thetas = rep(theta_nb, 4),
                         pis    = pi_grupos))
cat(sprintf("\nLogLik-oraculo (parametros verdadeiros): %.1f\n", ll_oraculo))

# ------------------------------------------------------------
# 3) INICIALIZACAO POR QUANTIS
# ------------------------------------------------------------
# 4 cortes que aproximam os 4 grupos verdadeiros (50/25/20/5)
brks <- quantile(consultas, probs = c(0, 0.50, 0.75, 0.95, 1))
brks <- unique(brks)   # garante quebras unicas
init_cluster <- as.integer(cut(consultas, breaks = brks,
                               include.lowest = TRUE))
cat("\nDistribuicao da inicializacao por quantis:\n")
print(table(init_cluster))

# Diagnostico: pureza da inicializacao vs grupo real
cat("\nMatriz de confusao da INICIALIZACAO vs grupo real:\n")
print(table(real = df$grupo_real, init = init_cluster))

# ------------------------------------------------------------
# 4) BASELINE: NB padrao com inicializacao aleatoria (k=2..5)
# ------------------------------------------------------------
set.seed(42)
fx_random <- suppressWarnings(stepFlexmix(
  consultas ~ 1, data = df, k = 2:5, nrep = 10,
  model = FLXMRnegbin(),
  control = list(iter.max = 500, minprior = 0.01, tolerance = 1e-7)
))
cat("\n--- Baseline (init aleatoria) ---\n"); print(fx_random)

# ------------------------------------------------------------
# 5) E1: NB com inicializacao por quantis, k=4 (theta livre)
# ------------------------------------------------------------
cat("\n--- E1: NB k=4 + init por quantis (theta livre) ---\n")
fit_E1 <- suppressWarnings(flexmix(
  consultas ~ 1, data = df, k = 4,
  cluster = init_cluster,
  model   = FLXMRnegbin(),
  control = list(iter.max = 1000, minprior = 0.01, tolerance = 1e-8)
))
cat(sprintf("k_efetivo=%d  logLik=%.1f  BIC=%.1f\n",
            fit_E1@k, as.numeric(logLik(fit_E1)), BIC(fit_E1)))

# ------------------------------------------------------------
# 6) E2: NB com theta fixo = 2.5 (init aleatoria)
# ------------------------------------------------------------
cat("\n--- E2: NB k=2..5 com theta=2.5 fixo (init aleatoria) ---\n")
set.seed(42)
fx_fixed <- suppressWarnings(stepFlexmix(
  consultas ~ 1, data = df, k = 2:5, nrep = 10,
  model   = FLXMRnegbin(theta = theta_nb),
  control = list(iter.max = 500, minprior = 0.01, tolerance = 1e-7)
))
print(fx_fixed)

# ------------------------------------------------------------
# 7) E3: NB com theta fixo + init por quantis, k=4
# ------------------------------------------------------------
cat("\n--- E3: NB k=4 + theta=2.5 fixo + init por quantis ---\n")
fit_E3 <- suppressWarnings(flexmix(
  consultas ~ 1, data = df, k = 4,
  cluster = init_cluster,
  model   = FLXMRnegbin(theta = theta_nb),
  control = list(iter.max = 1000, minprior = 0.01, tolerance = 1e-8)
))
cat(sprintf("k_efetivo=%d  logLik=%.1f  BIC=%.1f\n",
            fit_E3@k, as.numeric(logLik(fit_E3)), BIC(fit_E3)))

# ------------------------------------------------------------
# 8) RESUMO COMPARATIVO
# ------------------------------------------------------------
resumo <- function(nome, fit) {
  if (is.null(fit)) return(data.frame(modelo=nome, k_ef=NA, logLik=NA, BIC=NA))
  data.frame(
    modelo = nome,
    k_ef   = fit@k,
    logLik = round(as.numeric(logLik(fit)), 1),
    BIC    = round(BIC(fit), 1)
  )
}

best_random_k4 <- tryCatch(getModel(fx_random, which = "4"), error = function(e) NULL)
best_fixed_k4  <- tryCatch(getModel(fx_fixed,  which = "4"), error = function(e) NULL)
best_random_BIC <- getModel(fx_random, which = "BIC")
best_fixed_BIC  <- getModel(fx_fixed,  which = "BIC")

tab_resumo <- do.call(rbind, list(
  data.frame(modelo = "Oraculo (k=4 verdadeiro)",
             k_ef = 4, logLik = round(ll_oraculo, 1), BIC = NA),
  resumo("Baseline aleatorio (melhor BIC)",      best_random_BIC),
  resumo("Baseline aleatorio k=4 forcado",       best_random_k4),
  resumo("E1: init quantis, theta livre, k=4",   fit_E1),
  resumo("E2: theta fixo (melhor BIC)",          best_fixed_BIC),
  resumo("E2: theta fixo k=4 forcado",           best_fixed_k4),
  resumo("E3: theta fixo + init quantis, k=4",   fit_E3)
))
cat("\n=========================================================\n")
cat("RESUMO COMPARATIVO\n")
cat("=========================================================\n")
print(tab_resumo, row.names = FALSE)
cat(sprintf("\nDelta logLik p/ favorecer k=4 sobre k=2 no BIC: ~%.1f (com theta livre)\n",
            4 * log(n) / 2))

# ------------------------------------------------------------
# 9) Para o melhor candidato a 4 grupos, ver parametros + confusao
# ------------------------------------------------------------
avaliar <- function(nome, fit) {
  if (is.null(fit) || fit@k != 4) {
    cat(sprintf("\n[%s] k=%s, pulando.\n", nome, ifelse(is.null(fit),"NA",fit@k)))
    return(invisible(NULL))
  }
  cat(sprintf("\n--- [%s] parametros estimados ---\n", nome))
  pars <- parameters(fit)
  print(round(pars, 3))
  # Com theta livre, pars e matriz (2 x k); com theta fixo, pars e vetor (k)
  if (is.matrix(pars)) {
    intercepts <- pars[1, ]
  } else {
    intercepts <- as.numeric(pars)
  }
  mus <- exp(intercepts)
  ord <- order(mus)
  cat("mu estimados (ordenados):  ", paste(round(mus[ord], 2), collapse = ", "), "\n")
  cat("mu verdadeiros (ordenados):", paste(mu_grupos, collapse = ", "), "\n")
  cat("pesos estimados (ord):     ", paste(round(prior(fit)[ord], 3), collapse = ", "), "\n")
  cat("pesos verdadeiros (ord):   ", paste(sort(pi_grupos), collapse = ", "), "\n")

  cl <- factor(clusters(fit), levels = 1:fit@k)   # garante todas as k colunas
  conf <- table(real = df$grupo_real, est = cl)
  cat("\nMatriz de confusao:\n"); print(conf)

  k <- fit@k
  if (k <= 6) {
    perms   <- gtools::permutations(k, k)
    acertos <- apply(perms, 1, function(p) sum(diag(conf[, p, drop = FALSE])))
    cat(sprintf("Pureza apos pareamento otimo: %.1f%%\n",
                100 * max(acertos) / sum(conf)))
  }

  post     <- posterior(fit)
  max_post <- apply(post, 1, max)
  cat(sprintf("Certeza media: %.3f  |  %%>0.9: %.1f%%\n",
              mean(max_post), 100 * mean(max_post > 0.9)))
}

avaliar("E1: init quantis, theta livre",      fit_E1)
avaliar("Baseline aleatorio k=4 forcado",     best_random_k4)
avaliar("E2: theta fixo k=4 forcado",         best_fixed_k4)
avaliar("E3: theta fixo + init quantis",      fit_E3)
