# ============================================================
# fig08 adaptado para o experimento univariado (so consultas)
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
pi_grupos <- c(0.50, 0.25, 0.20, 0.05)
mu_grupos <- c(2, 12, 25, 50)
theta_nb  <- 2.5

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
