############## GRÁFICOS COMPARATIVOS PARA O C-ÍNDEX DAS SIMULAÇÕES

# C-índex de Random Machines com Seleção de subconjunto de variáveis
rm_com_selecao <- c(0.736, 0.69, 0.673, 0.647, 0.733, 0.645, 0.645, 0.766, 0.734,
                    0.714, 0.693, 0.67, 0.752)

# C-índex de Random Machines sem Seleção de subconjunto de variáveis
rm_sem_selecao <- c(0.765, 0.774, 0.766, 0.779, 0.765, 0.798, 0.804, 0.792, 0.752,
                    0.757, 0.808, 0.797, 0.77)

# C-índex de Random Survival Forest
resultados_rsf <- readRDS("outputs/resultado_simulacoes_rsf.rds")
rsf <- resultados_rsf$lista_cindex

# C-índex de SVM linear
resultados_lin <- readRDS("outputs/resultado_simulacoes_svm_lin.rds")
svm_linear <- resultados_lin

# C-índex de SVM Add
resultados_add <- readRDS("outputs/resultado_simulacoes_svm_add.rds")
svm_add <- resultados_add

# C-índex de SVM Rbf
resultados_rbf <- readRDS("outputs/resultado_simulacoes_svm_rbf.rds")
svm_rbf <- resultados_rbf


# Vetor com todos os C-índex
cindex <- c(rm_com_selecao, rm_sem_selecao, rsf, svm_linear, svm_add, svm_rbf)

# Vetor indicando os modelos
modelo <- c(rep("RM com seleção", 13),
            rep("RM", 13),
            rep("RSF", 13),
            rep("SVM linear", 13),
            rep("SVM aditivo", 13),
            rep("SVM gaussiano", 13)
            )

# Vetor indicando o número da simualação
simulacao <- 1:13

df <- data.frame(simulacao, modelo, cindex)

# Boxplot para C-índex das simulações por modelo
library(ggplot2)

bp <- ggplot(df, aes(x = factor(modelo), y = cindex)) +
  geom_boxplot() +
   labs(
    title = "",
    x = "Modelo",
    y = "C-índex"
  ) +
  theme_minimal(base_size = 16)

ggsave(
  filename = "outputs/boxplot_simulacoes.png",
  plot = bp,
  width = 8, height = 6, dpi = 300
)

# Gráfico de linhas para comparar os C-índex
linhas <- ggplot(df, aes(x = simulacao, y = cindex, color = modelo, group = modelo)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_x_continuous(breaks = df$simulacao) +
  labs(
    title = " ",
    x = "Número da simulação",
    y = "C-índex",
    color = "Modelo"
  ) +
  theme_minimal(base_size = 14)

ggsave(
  filename = "outputs/linhas_simulacoes.png",
  plot = linhas,
  width = 8, height = 6, dpi = 300
)
