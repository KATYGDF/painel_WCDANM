df <- read.csv("C:/Katy/Doutorado/painel_WCDANM/dados/dados_simulados.csv")
vars <- c("consultas","ps","exames","terapias","internacoes")
cat("--- Médias por grupo ---\n")
for(v in vars){
  cat(sprintf("%-14s", v))
  m <- round(tapply(df[[v]], df$grupo_latente, mean), 2)
  cat(paste(names(m),"=",m, collapse="  "), "\n")
}
cat("\n--- Proporção de zeros (%) ---\n")
for(v in vars) cat(sprintf("%-14s %.1f\n", v, mean(df[[v]]==0)*100))
cat("\n--- n por grupo ---\n")
print(table(df$grupo_latente))
cat("\n--- Razão var/média ---\n")
for(v in vars) cat(sprintf("%-14s %.1f\n", v, var(df[[v]])/mean(df[[v]])))
cat("\n--- Gasto médio por grupo ---\n")
m <- round(tapply(df$gasto_total, df$grupo_latente, mean))
cat(paste(names(m),"=",m, collapse="  "),"\n")
cat(sprintf("Media global: %.0f  Mediana: %.0f  P95: %.0f\n",
  mean(df$gasto_total), median(df$gasto_total), quantile(df$gasto_total,.95)))
