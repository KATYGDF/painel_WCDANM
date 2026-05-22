# ============================================================
# GRAFICO COMPARATIVO: densidade ajustada por estrategia
# Assume que o script de recuperacao ja foi executado e que
# os objetos fit_E1, fit_E3, best_random_k4, best_fixed_k4,
# df, mu_grupos, pi_grupos, theta_nb existem no ambiente.
# ============================================================

library(ggplot2)
library(dplyr)
library(tidyr)

# ------------------------------------------------------------
# 1) Funcao: densidade da mistura NB a partir de um fit flexmix
# ------------------------------------------------------------
densidade_mistura <- function(y_grid, fit, theta_fixo = NULL) {
  pars <- parameters(fit)
  if (is.matrix(pars)) {
    mus    <- exp(pars[1, ])
    thetas <- pars[2, ]
  } else {
    mus    <- exp(as.numeric(pars))
    thetas <- rep(theta_fixo, length(mus))
  }
  pis <- prior(fit)
  sapply(y_grid, function(y) {
    sum(pis * dnbinom(y, size = thetas, mu = mus))
  })
}

# Densidade teorica do modelo verdadeiro (oraculo)
densidade_oraculo <- function(y_grid, mus, thetas, pis) {
  sapply(y_grid, function(y) {
    sum(pis * dnbinom(y, size = thetas, mu = mus))
  })
}

# ------------------------------------------------------------
# 2) Monta data.frame com densidades de cada estrategia
# ------------------------------------------------------------
y_max  <- quantile(df$consultas, 0.995)  # corta cauda extrema p/ visualizar
y_grid <- 0:ceiling(y_max)

df_dens <- bind_rows(
  data.frame(y = y_grid, dens = densidade_oraculo(
    y_grid, mu_grupos, rep(theta_nb, 4), pi_grupos),
    modelo = "Oraculo (verdadeiro, k=4)"),
  data.frame(y = y_grid, dens = densidade_mistura(y_grid, best_random_k4),
    modelo = "Baseline aleatorio (k=4 forcado)"),
  data.frame(y = y_grid, dens = densidade_mistura(y_grid, fit_E1),
    modelo = "E1: init quantis, theta livre"),
  data.frame(y = y_grid, dens = densidade_mistura(y_grid, best_fixed_k4,
                                                  theta_fixo = theta_nb),
    modelo = "E2: theta fixo (k=4 forcado)"),
  data.frame(y = y_grid, dens = densidade_mistura(y_grid, fit_E3,
                                                  theta_fixo = theta_nb),
    modelo = "E3: theta fixo + init quantis")
)
df_dens$modelo <- factor(df_dens$modelo,
                         levels = c("Oraculo (verdadeiro, k=4)",
                                    "Baseline aleatorio (k=4 forcado)",
                                    "E1: init quantis, theta livre",
                                    "E2: theta fixo (k=4 forcado)",
                                    "E3: theta fixo + init quantis"))

# ------------------------------------------------------------
# 3) Linhas verticais nas medias verdadeiras
# ------------------------------------------------------------
df_mus <- data.frame(mu = mu_grupos,
                     grupo = paste0("G", 1:4),
                     y_pos = max(df_dens$dens) * 0.95)

# ------------------------------------------------------------
# 4) GRAFICO PRINCIPAL: histograma + densidade ajustada
# ------------------------------------------------------------
df_hist <- df %>% filter(consultas <= y_max)

p_main <- ggplot() +
  geom_histogram(data = df_hist,
                 aes(x = consultas, y = after_stat(density)),
                 bins = 60, fill = "grey80", color = "white", linewidth = 0.1) +
  geom_line(data = df_dens,
            aes(x = y, y = dens, color = modelo),
            linewidth = 1) +
  geom_vline(data = df_mus, aes(xintercept = mu),
             linetype = "dashed", color = "red", alpha = 0.5) +
  geom_text(data = df_mus,
            aes(x = mu, y = y_pos, label = grupo),
            color = "red", size = 3, vjust = 0, hjust = -0.2) +
  facet_wrap(~ modelo, ncol = 1) +
  scale_color_brewer(palette = "Dark2", guide = "none") +
  labs(title = "Densidade ajustada por estrategia vs. histograma observado",
       subtitle = "Linhas vermelhas tracejadas = medias verdadeiras (2, 12, 25, 50)",
       x = "Consultas", y = "Densidade") +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold"))

print(p_main)

# ------------------------------------------------------------
# 5) GRAFICO SECUNDARIO: dumbbell de mu estimados vs verdadeiros
# ------------------------------------------------------------
extrai_mus <- function(fit, nome) {
  pars <- parameters(fit)
  mus  <- if (is.matrix(pars)) exp(pars[1, ]) else exp(as.numeric(pars))
  pis  <- prior(fit)
  data.frame(modelo = nome, mu_est = sort(mus),
             peso = pis[order(mus)],
             grupo_alvo = paste0("G", 1:length(mus)))
}

df_mu_comp <- bind_rows(
  data.frame(modelo = "Verdadeiro", mu_est = mu_grupos,
             peso = pi_grupos,
             grupo_alvo = paste0("G", 1:4)),
  extrai_mus(best_random_k4, "Baseline aleatorio"),
  extrai_mus(fit_E1,         "E1"),
  extrai_mus(best_fixed_k4,  "E2"),
  extrai_mus(fit_E3,         "E3")
)
df_mu_comp$modelo <- factor(df_mu_comp$modelo,
                            levels = c("Verdadeiro","Baseline aleatorio",
                                       "E1","E2","E3"))

p_mu <- ggplot(df_mu_comp,
               aes(x = mu_est, y = modelo, color = grupo_alvo, size = peso)) +
  geom_point(alpha = 0.8) +
  geom_vline(xintercept = mu_grupos,
             linetype = "dashed", color = "grey60") +
  scale_size_continuous(range = c(3, 10), name = "Peso π") +
  scale_color_brewer(palette = "Set1", name = "Grupo (ord)") +
  labs(title = "Medias estimadas por estrategia",
       subtitle = "Linhas tracejadas = medias verdadeiras  |  tamanho do ponto = peso do componente",
       x = "μ estimada", y = NULL) +
  theme_minimal(base_size = 11)

print(p_mu)
