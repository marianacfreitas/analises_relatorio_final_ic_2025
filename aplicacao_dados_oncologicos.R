############### IMPORTANDO FUNÇÕES SURVIVAL RANDOM MACHINES
source("funcoes_survival_random_machines.R")

cat("Arquivo de funções random machines carregado \n")

################### IMPORTANDO DADOS TRATADOS

dados <- read.csv("data/dados_tratados.csv") |>
   select(-c(X))

cat("dados_tratados carregado \n")

# Criar fórmula de sobrevivência
 formula_dados <-  as.formula(
   paste0("Surv(time_years, death)", " ~ ",
          paste(setdiff(names(dados), c("time_years", "death")), collapse = " + "))
 )
 
 
 # Definindo dados de treino e teste da amostra dos dados
 cv_data <- split_data(dados, time_var = "time_years", status_var = "death", 
                       seed=2402)
 
 treino <- cv_data$train_sample
 teste <- cv_data$test_sample
 
 cat("Dados divididos em treino e teste \n")
 
 cat("O melhor parâmetro foi", cost_escolhido, "\n")
 
 cat("O modelo random machines com seleção está sendo treinado \n")
   
# Treinar modelo com a função survival_random_machines
rm_model <- survival_random_machines(formula = formula_dados, 
                                        train = treino,
                                        test = teste,
                                        n_cores = 6,
                                        boots_number = 100, 
                                        cost = 1,
                                        selecionar_variaveis = TRUE,
                                        seed.bootstrap = 2402,
                                        save = F)
 
 
cat("O modelo ensemble foi gerado \n")
cat("O modelo final está sendo gerado \n")

modelo_survival_rm <- final_surv_rm_model(mod = rm_model, newdata = teste, save = T, save_path = "outputs/modelo_final_random_machines.rds")

cat("O modelo final random machines com seleção foi gerado \n")

# Removendo os modelos com seleção do environment
rm(rm_model, modelo_survival_rm)

# Treinar modelo sem seleção com a função survival_random_machines
rm_model <- survival_random_machines(formula = formula_dados, 
                                     train = treino,
                                     test = teste,
                                     n_cores = 6,
                                     boots_number = 100, 
                                     cost = cost_escolhido,
                                     selecionar_variaveis = FALSE,
                                     seed.bootstrap = 2402,
                                     save = F)


cat("O modelo ensemble sem selação foi gerado \n")
cat("O modelo final sem selação está sendo gerado \n")

modelo_survival_rm <- final_surv_rm_model(mod = rm_model, newdata = teste, save = T, save_path = "outputs/modelo_final_sem_selecao_random_machines.rds")

cat("O modelo final random machines sem seleção foi gerado \n")

# Limpando o environment
rm(list = ls())
cat("O ambiente foi limpo \n")
