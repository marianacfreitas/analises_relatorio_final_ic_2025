######################### ANÁLISES DE SHAP

# Define o mirror do CRAN para evitar erro ao instalar pacotes no cluster
options(repos = c(CRAN = "https://cran.r-project.org"))


# Lista de pacotes necessários
pacotes <- c("readr", "survival", "survminer", 
             "dplyr", "survivalsvm", "randomForestSRC", "iml", "ggplot2",
             "ggsci", "patchwork")


# Função para instalar e carregar pacotes com dependências
instalar_e_carregar <- function(pkg) {
  if (!suppressWarnings(require(pkg, character.only = TRUE))) {
    suppressMessages(install.packages(pkg, dependencies = TRUE))
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# Instalação e carregamento silencioso
invisible(sapply(pacotes, instalar_e_carregar))


cat("Pacotes carregados \n")

### Funções usadas em Random Machines
source("funcoes_survival_random_machines.R")

cat("Funções random machines carregadas \n")

### Carregando os dados
dados <- read_csv("data/dados_tratados.csv") |> 
  select(-c(`...1`)) |> 
  mutate(across(-c(time_years, death), as.factor))

cat("Dados tratados carregados \n")

### Definindo dados de treino e teste
cv_data <- split_data(dados, time_var = "time_years", status_var = "death", seed=2402)
treino <- cv_data$train_sample
teste <- cv_data$test_sample

cat("Dados divididos em treino e teste \n")

cat("Começando análise SHAP \n")

options(future.globals.maxSize = 60 * 1024^3)  # 60 GB

teste2 <- teste |>
  select(-c(time_years, death))

modelo_survival_rm <- readRDS("outputs/modelo_final_random_machines.rds")

cat("Modelo final random machines com seleção carregado \n")

# Criando objeto predictor
predictor_rm <- Predictor$new(
  model = modelo_survival_rm,     
  data = teste2,
  y = teste$time_years,
  predict.function = predict_from_model
)

cat("Rodando SHAP para modelo final random machines com seleção \n")

#### Gráfico de importância gregado para todo o conjunto de teste (lento!)
imp <- FeatureImp$new(predictor_rm, loss = "mae")

cat("Importância SHAP finalizado \n")

saveRDS(imp, "outputs/objeto_importancia_shap.rds")

cat("Objeto importância SHAP salvo \n")

# Visualização
grafico_importancia <- plot(imp) +
  labs(
    x = "Importância da variável (erro: MAE)",
    y = NULL
  ) +
  scale_color_jco() +
  theme_minimal(base_size = 14)

ggsave("outputs/grafico_importancia_shap.png", grafico_importancia, width = 8, height = 6, dpi = 300)

cat("Gráfico importância SHAP salvo \n")

### Gráfico de dependência SHAP para top 3 variáveis

top_variaveis <- imp$results |>
  arrange(desc(importance)) |>
  slice(1:3) |>
  pull(feature)

for (var in top_variaveis) {
  efeito <- FeatureEffect$new(predictor_rm, feature = var, method = "pdp")
  saveRDS(efeito, paste0("outputs/objeto_dependencia_shap_", var, ".rds"))
  cat("Objeto de dependência salvo para", var, "\n")
  
  grafico_dependencia <- plot(efeito) +
    labs(
      x = var,
      y = "Tempo até o óbito predito"
    ) +
    scale_color_jco() +
    theme_minimal(base_size = 14)
  
  ggsave(
    filename = paste0("outputs/grafico_dependencia_shap_", var, ".png"),
    plot = grafico_dependencia,
    width = 8, height = 6, dpi = 300
  )
  cat("Gráfico de dependência salvo para", var, "\n")
}

### Gráfico de resumo SHAP

efeitos_todos <- FeatureEffects$new(predictor_rm, method = "pdp+ice")
saveRDS(efeitos_todos, "outputs/objeto_resumo_shap.rds")
cat("Objeto de resumo salvo \n")

# Define as variáveis relevantes manualmente (exemplo com 8)
variaveis_relevantes <- imp$results |>
  arrange(desc(importance)) |>
  slice(1:8) |>
  pull(feature)

# Filtra efeitos
efeitos_filtrados <- efeitos_todos$effects[names(efeitos_todos$effects) %in% variaveis_relevantes]

# Gera os plots individuais
graficos_individuais <- mapply(
  FUN = function(efeito, nome_var) {
    plot(efeito) +
      labs(
        title = nome_var,
        y = "Tempo de sobrevida predito",
        x = NULL
      ) +
      scale_color_jco() +
      theme_minimal(base_size = 12)
  },
  efeito = efeitos_filtrados,
  nome_var = names(efeitos_filtrados),
  SIMPLIFY = FALSE
)

# Divide os gráficos em dois grupos de 4
graficos_1 <- graficos_individuais[1:4]
graficos_2 <- graficos_individuais[5:8]

# Junta todos os gráficos
painel_1 <- wrap_plots(graficos_1, ncol = 2) +
  plot_annotation(title = "Resumo SHAP - Parte 1")

painel_2 <- wrap_plots(graficos_2, ncol = 2) +
  plot_annotation(title = "Resumo SHAP - Parte 2 ")


ggsave("outputs/grafico_resumo_shap_1.png", painel_1, width = 10, height = 7, dpi = 300)
ggsave("outputs/grafico_resumo_shap_2.png", painel_2, width = 10, height = 7, dpi = 300)

cat("Gráfico de resumo SHAP salvo como 'grafico_resumo_shap.png'\n")

rm(list = ls())
cat("O ambiente foi limpo \n")


