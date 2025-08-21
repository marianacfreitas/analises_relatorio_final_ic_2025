
######################  FUNÇÕES DE SURVIVAL RANDOM MACHINES 

source("funcoes_survival_random_machines.R")

cat("Funções random machines carregadas \n")

##################### DEFININDO OS PARÂMETROS PARA SIMULAÇÃO 
## Parâmetros definidos com base no ajuste do modelo exponencial aos dados da FOSP

dados <- read.csv("data/dados_tratados.csv") |>
  select(-c(X))

cat("Dados tratados carregados \n")

dados$time_years <- ifelse(dados$time_years == 0, 0.01, dados$time_years)

modelo_exp <- survreg(Surv(time_years, death) ~ ., data = dados, dist = "exponential")

betas <- coefficients(modelo_exp)

cat("Modelo exponencial foi ajustado aos dados oncológicos \n")

##################### GERANDO DADOS DE SIMULAÇÂO

cat("Executando simulação sem seleção \n")

resultados_sim_sem_selecao <- simulacoes_rm(
  n_sim = 13,        # número de simulações
  n = 2000,           # tamanho de cada conjunto simulado
  p = 23,            # número de variáveis preditoras
  beta = betas, # vetor de parâmetros
  boots_number = 100, # número de modelos bootstrap no ensemble
  n_cores = 6,
  selecionar_variaveis = F, # selecionar subconjunto de variáveis para cada amostra bootstrap
  seed = 2402        # semente base para reprodutibilidade
)

cat("A simulação sem seleção foi finalizada \n")


saveRDS(resultados_sim_sem_selecao, "outputs/resultados_simulacao_sem_selecao.rds")

cat("Os resultados da simulação sem seleção foram salvos \n")

cat("Executando simulação com seleção \n")

resultados_sim_com_selecao <- simulacoes_rm(
  n_sim = 13,        # número de simulações
  n = 2000,           # tamanho de cada conjunto simulado
  p = length(betas),            # número de variáveis preditoras
  beta = betas, # vetor de parâmetros
  boots_number = 100, # número de modelos bootstrap no ensemble
  n_cores = 6,
  selecionar_variaveis = T, # selecionar subconjunto de variáveis para cada amostra bootstrap
  seed = 2402        # semente base para reprodutibilidade
)

cat("A simulação com seleção foi finalizada \n")

saveRDS(resultados_sim_com_selecao, "outputs/resultados_simulacao_com_selecao.rds")

cat("Os resultados da simulação com seleção foram salvos \n")

rm(list = ls())

cat("O ambiente foi limpo \n")
