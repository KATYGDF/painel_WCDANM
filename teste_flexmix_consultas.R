# ============================================================
# EXPERIMENTO CONTROLADO - FLEXMIX UNIVARIADO EM CONSULTAS
# Objetivo: verificar se, com medias BEM SEPARADAS,
# o stepFlexmix recupera os 4 grupos latentes.
# Espelha a analise da Parte 4 (todas as funcoes do flexmix).
# ============================================================

library(dplyr)
library(ggplot2)
library(flexmix)

# FLXMRnegbin vem do countreg (R-Forge).
# Se nao estiver instalado, rode UMA vez:
#   install.packages("countreg", repos = "http://R-Forge.R-project.org")
if (!requireNamespace("countreg", quietly = TRUE)) {
  stop("Pacote 'countreg' nao encontrado. Instale com:\n",
       "  install.packages(\"countreg\", repos = \"http://R-Forge.R-project.org\")")
}
FLXMRnegbin <- countreg::FLXMRnegbin

set.seed(1234)

# ------------------------------------------------------------
# 1) SIMULACAO: apenas consultas, NB, medias bem separadas
# ------------------------------------------------------------
n         <- 5000
pi_grupos <- c(0.50, 0.25, 0.20, 0.05)
mu_grupos <- c(2, 12, 25, 50)          # bem separadas
theta_nb  <- 2.5

grupo_real <- sample(1:4, size = n, replace = TRUE, prob = pi_grupos)
consultas  <- rnbinom(n = n, size = theta_nb, mu = mu_grupos[grupo_real])

df <- data.frame(
  id         = 1:n,
  grupo_real = factor(grupo_real, levels = 1:4,
                      labels = c("G1_baixo","G2_ambul","G3_agudo","G4_atipico")),
  consultas  = consultas
)

cat("\n--- Estatisticas por grupo verdadeiro ---\n")
print(
  df %>% group_by(grupo_real) %>%
    summarise(n = n(), prop = round(n()/nrow(df),3),
              media = round(mean(consultas),2),
              mediana = median(consultas),
              sd = round(sd(consultas),2),
              .groups = "drop")
)

# ------------------------------------------------------------
# 2) stepFlexmix: k = 1..6, NB, multiplos reinicios
# ------------------------------------------------------------
set.seed(42)
fx_steps <- suppressWarnings(
  stepFlexmix(
    consultas ~ 1,
    data    = df,
    k       = 1:5,
    nrep    = 5,
    model   = FLXMRnegbin(),
    control = list(iter.max = 300, minprior = 0.02)
  )
)
cat("\n--- Sumario stepFlexmix ---\n")
print(fx_steps)

# ------------------------------------------------------------
# 3) Tabela de criterios BIC / AIC / ICL para todos os k
# ------------------------------------------------------------
crit_tab <- do.call(rbind, lapply(1:6, function(g) {
  mod <- tryCatch(getModel(fx_steps, which = as.character(g)),
                  error = function(e) NULL)
  if (is.null(mod)) {
    return(data.frame(g = g, k_efetivo = NA, logLik = NA,
                      BIC = NA, AIC = NA, ICL = NA))
  }
  data.frame(
    g         = g,
    k_efetivo = mod@k,
    logLik    = round(as.numeric(logLik(mod)), 1),
    BIC       = round(BIC(mod), 1),
    AIC       = round(AIC(mod), 1),
    ICL       = round(ICL(mod), 1)
  )
}))
cat("\n--- Criterios de selecao por k ---\n")
print(crit_tab, row.names = FALSE)

# ------------------------------------------------------------
# 4) getModel(): recupera modelos por criterio
# ------------------------------------------------------------
best_bic <- getModel(fx_steps, which = "BIC")
best_aic <- getModel(fx_steps, which = "AIC")
best_icl <- getModel(fx_steps, which = "ICL")

cat(sprintf("\nk escolhido  ->  BIC: %d  |  AIC: %d  |  ICL: %d\n",
            best_bic@k, best_aic@k, best_icl@k))

# Para o restante das analises usamos o modelo do BIC
best <- best_bic
cat(sprintf("\n>>> Seguindo com o modelo BIC (k = %d)\n", best@k))

# ------------------------------------------------------------
# 5) parameters(): mu e theta por componente
# ------------------------------------------------------------
pars <- parameters(best)
cat("\n--- parameters(best) ---\n"); print(round(pars, 3))

mu_est    <- exp(pars[1, ])
theta_est <- pars[2, ]
var_teo   <- mu_est + mu_est^2 / theta_est

tab_comp <- data.frame(
  componente   = paste0("C", seq_along(mu_est)),
  peso         = round(prior(best), 3),
  mu           = round(mu_est, 2),
  theta        = round(theta_est, 3),
  var_teorica  = round(var_teo, 1),
  razao_var_mu = round(var_teo / mu_est, 2)
)
cat("\n--- Parametros por componente (escala original) ---\n")
print(tab_comp, row.names = FALSE)

cat("\nMedias-alvo (ordenadas):    ", paste(sort(mu_grupos), collapse = ", "), "\n")
cat("Medias estimadas (ordenadas):", paste(round(sort(mu_est),2), collapse = ", "), "\n")
cat("Proporcoes-alvo (ordenadas):    ", paste(sort(pi_grupos), collapse = ", "), "\n")
cat("Pesos estimados (ordenados):    ", paste(round(sort(prior(best)),3), collapse = ", "), "\n")

# ------------------------------------------------------------
# 6) posterior() e diagnostico de incerteza
# ------------------------------------------------------------
post     <- posterior(best)
max_post <- apply(post, 1, max)
ent      <- -rowSums(post * log(post + 1e-15))
ent_rel  <- mean(ent) / log(best@k)   # entropia relativa (0 = certeza, 1 = uniforme)

cat(sprintf("\nCerteza media: %.3f  |  mediana: %.3f  |  %%>0.9: %.1f%%\n",
            mean(max_post), median(max_post), 100 * mean(max_post > 0.9)))
cat(sprintf("Entropia media: %.4f  (max teorico log(k) = %.4f)  |  ent_rel = %.1f%%\n",
            mean(ent), log(best@k), 100 * ent_rel))

# ------------------------------------------------------------
# 7) clusters(): classificacao MAP e matriz de confusao
# ------------------------------------------------------------
df$cluster_est <- clusters(best)

cat("\n--- Tamanho de cada cluster estimado ---\n")
print(table(df$cluster_est))

cat("\n--- Matriz de confusao: grupo real x cluster estimado ---\n")
conf <- table(real = df$grupo_real, estimado = df$cluster_est)
print(conf)

# Pareamento otimo (Hungarian simples por permutacao) - ate k=4 e barato
if (best@k <= 6) {
  perms  <- gtools::permutations(best@k, best@k)
  acertos <- apply(perms, 1, function(p) {
    sum(diag(conf[, p, drop = FALSE]))
  })
  melhor_perm <- perms[which.max(acertos), ]
  pureza <- max(acertos) / sum(conf)
  cat(sprintf("\nPureza apos pareamento otimo: %.1f%%\n", 100 * pureza))
}

# ------------------------------------------------------------
# 8) refit() + summary(): SE e p-valores por componente
# ------------------------------------------------------------
cat("\n--- refit() + summary() (pode demorar) ---\n")
ref <- tryCatch(refit(best), error = function(e) {
  cat("refit falhou:", conditionMessage(e), "\n"); NULL
})
if (!is.null(ref)) print(summary(ref))

# ------------------------------------------------------------
# 9) Linha-resumo no formato da tabela comparativa
# ------------------------------------------------------------
resumo_linha <- data.frame(
  variavel  = "consultas",
  g_estrela = best@k,
  BIC       = round(BIC(best), 1),
  pi_min    = round(min(prior(best)), 3),
  pi_max    = round(max(prior(best)), 3),
  mu_min    = round(min(mu_est), 2),
  mu_max    = round(max(mu_est), 2),
  post_med  = round(median(max_post), 3),
  pct_gt09  = sprintf("%.1f%%", 100 * mean(max_post > 0.9)),
  ent_rel   = sprintf("%.1f%%", 100 * ent_rel)
)
cat("\n--- Linha-resumo (formato da tabela comparativa) ---\n")
print(resumo_linha, row.names = FALSE)

# ------------------------------------------------------------
# 10) Graficos
# ------------------------------------------------------------
p1 <- ggplot(df, aes(x = consultas, fill = grupo_real)) +
  geom_histogram(bins = 60, position = "stack", color = "white", linewidth = 0.1) +
  labs(title = "Consultas por grupo VERDADEIRO",
       x = "Consultas", y = "Frequencia", fill = "Grupo") +
  theme_minimal()

p2 <- ggplot(df, aes(x = consultas, fill = factor(cluster_est))) +
  geom_histogram(bins = 60, position = "stack", color = "white", linewidth = 0.1) +
  labs(title = paste0("Consultas por cluster ESTIMADO (k=", best@k, ")"),
       x = "Consultas", y = "Frequencia", fill = "Cluster") +
  theme_minimal()

p3 <- ggplot(data.frame(max_post = max_post), aes(x = max_post)) +
  geom_histogram(bins = 40, fill = "steelblue", color = "white") +
  geom_vline(xintercept = 0.9, linetype = "dashed", color = "red") +
  labs(title = "Certeza de classificacao (max posterior por observacao)",
       x = "Maior probabilidade posterior", y = "Frequencia") +
  theme_minimal()

print(p1); print(p2); print(p3)
