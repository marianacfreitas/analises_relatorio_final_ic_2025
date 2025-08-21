# Lista de pacotes necessários
pacotes <- c("readr", "survival", "survminer", 
             "dplyr", "survivalsvm", "randomForestSRC", "iml", "ggplot2",
             "foreach", "doParallel", "fs")


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

################## MODELOS PARA COMPARAÇÃO

### Funções usadas em Random Machines
source("funcoes_survival_random_machines.R")

cat("Funções random machines carregadas \n")

### Carregando os dados
dados <- read_csv("data/dados_tratados.csv") |> 
  select(-c(`...1`)) |> 
  mutate(across(-c(time_years, death), as.factor))

cat("Dados tratados carregados \n")

### Definindo dados de treino e teste dos dados
cv_data <- split_data(dados, time_var = "time_years", status_var = "death", seed=2402)
treino <- cv_data$train_sample
teste <- cv_data$test_sample

cat("Dados divididos em treino e teste \n")

###################### Aplicando o modelo de regressão de cox ###############
modelo_cox <- coxph(Surv(time_years, death) ~ ., data = treino, x = TRUE, y = TRUE)
surv_cox <- survfit(modelo_cox, newdata = teste)

cat("Modelo de regressão de cox aplicado aos dados oncológicos \n")

# Teste para verificar se os riscos são proporcionais
test_ph <- survival::cox.zph(modelo_cox)
test_ph 

saveRDS(test_ph, "outputs/testes_riscos_proporcionais_dados_oncologicos.rds")
cat("Teste de riscos proporcionais feito \n")

######################## Aplicando survival random forest #######################
modelo_rsf <- randomForestSRC::rfsrc(Surv(time_years, death) ~ ., data = as.data.frame(treino), mtry = floor(sqrt(ncol(teste)-2)),
                    splitrule = "logrank", ntree = 100)

cat("Random survival forest foi aplicado \n")

### C-índex do survival random forest
rsf_pred_aux <- predict(modelo_rsf, newdata = teste)
rsf_pred <- rsf_pred_aux$predicted # predição do tempo médio para cada indivíduo

cindex_rsf <- concordance(Surv(teste$time_years, teste$death) ~ rsf_pred)$concordance
cindex_rsf

cat("C-índex do RSF:", cindex_rsf, "\n")

##################### Aplicando SVM com diferentes Kernels ##########################

formula_dados <- as.formula(
  paste0("Surv(time_years, death) ~ ",
         paste(setdiff(names(treino), c("time_years", "death")), collapse = " + "))
)

cat("Executando SVM Rbf \n")

svm_rbf <- survivalsvm(formula_dados,
                                data = treino,
                                type = "regression",
                                gamma.mu = 1,
                                kernel = "rbf_kernel")
 
cat("SVM Rbf executado \n")
 
pred_svm_rbf <- predict(svm_rbf, newdata = teste)$predicted
 
cindex_svm_rbf <- concordance(Surv(teste$time_years, teste$death) ~ as.vector(pred_svm_rbf))$concordance
cindex_svm_rbf
 
cat("C-index para o SVM Rbf:", cindex_svm_rbf, "\n")
 
cat("Executando o SVM aditivo \n")
 
svm_add <- survivalsvm(formula_dados,
                                data = treino,
                                type = "regression",
                                gamma.mu = 1,
                                kernel = "add_kernel")
 
cat("SVM aditivo executado \n")
 
pred_svm_add <- predict(svm_add, newdata = teste)$predicted
 
cindex_svm_add <- concordance(Surv(teste$time_years, teste$death) ~ as.vector(pred_svm_add))$concordance
cindex_svm_add
cat("C-index para SVM aditivo:", cindex_svm_add, "\n")
 
cat("Executando SVM linear \n")
 
svm_lin <- survivalsvm(formula_dados,
                                data = treino,
                                type = "regression",
                                gamma.mu = 1,
                                kernel = "lin_kernel")
 
cat("SVM linear executado \n")
 
pred_svm_lin <- predict(svm_lin, newdata = teste)$predicted
 
cindex_svm_lin <- concordance(Surv(teste$time_years, teste$death) ~ as.vector(pred_svm_lin))$concordance
cindex_svm_lin
cat("C-índex para SVM linear:", cindex_svm_lin, "\n")

c_indexes_dados_oncologicos <- c(cindex_rsf, cindex_svm_add, cindex_svm_lin, 
                                 cindex_svm_rbf)

names(c_indexes_dados_oncologicos) <- c("Random Survival Forest", "SVM Aditivo",
                                        "SVM Linear", "SVM Gaussiano")

saveRDS(c_indexes_dados_oncologicos, "outputs/cindex_modelos_dados_oncologicos.rds")

cat("Os C-índexes dos modelos para comparação foram salvos \n")


######################### COMPARAÇÃO PARA DADOS DE SIMULAÇÃO

cat("Começando a parte de dados de simulação \n")

dados$time_years <- ifelse(dados$time_years == 0, 0.01, dados$time_years)

modelo_exp <- survreg(Surv(time_years, death) ~ ., data = dados, dist = "exponential")

betas <- coefficients(modelo_exp)

cat("Modelo exponencial foi ajustado aos dados oncológicos \n")

# Hiperparâmetros
set.seed(2402)
numero_simulacoes <- 13
n <- 2000         # tamanho de cada conjunto simulado
p <- length(betas)            # número de variáveis preditoras
beta <- betas # vetor de parâmetros
n_cores <- 6
seed <- 2402        # semente base para reprodutibilidade

######## Random Survival Forest

cat("Simulação para Random Forest \n")

resultados <- numeric(numero_simulacoes) # vetor que armazenará os resultados

for (i in 1:numero_simulacoes) { # iteração para cada simulação
  cat("RSF: Simulação", i, "de", numero_simulacoes, "\n")
  
  # Geração dos dados
  dados_sim <- gerar_dados_simulacao(n = n, p = p, beta = beta, seed = seed + i)
  
  # Divisão treino/teste
  split <- split_data(dados_sim, seed = seed + i)
  treino <- split$train_sample
  teste <- split$test_sample
  
  
  # Treinamento do modelo
  modelo_rsf <- randomForestSRC::rfsrc(Surv(time, status) ~ ., data = as.data.frame(treino), mtry = floor(sqrt(ncol(treino)-2)),
                                       splitrule = "logrank", ntree = 100)
  
  ### C-índex do survival random forest
  rsf_pred_aux <- predict(modelo_rsf, newdata = teste)
  rsf_pred <- rsf_pred_aux$predicted # predição do tempo médio para cada indivíduo
  
  cindex_rsf <- concordance(Surv(teste$time, teste$status) ~ rsf_pred)$concordance
  cindex_rsf
  
  
  # Guardar C-index no teste
  resultados[i] <- cindex_rsf
  
  if(i==numero_simulacoes){
    modelo_exemplo <- modelo_rsf
  }
  
  rm(modelo_rsf)
}

resumo_resultados <- list(
  ultimo_modelo = modelo_exemplo,
  lista_cindex = resultados
)

saveRDS(resumo_resultados, "outputs/resultado_simulacoes_rsf.rds")

############# SVM RBF

cat("Simulação para SVM RBF \n")

# número de núcleos disponíveis
numCores <- detectCores() - 1
cl <- makeCluster(numCores)
registerDoParallel(cl)

# exportar funções e objetos criados no script
clusterExport(cl, c("gerar_dados_simulacao", "split_data", "n", "p", "beta", "seed", "numero_simulacoes"))

# rodar as simulações em paralelo
resultados_svm_rbf <- foreach(i = 1:numero_simulacoes,
                      .combine = c, # junta os resultados em um vetor
                      .packages = c("survival", "survivalsvm", "randomForestSRC")) %dopar% {
                        
                        msg <- paste("SVM Rbf: Simulação", i, "de", numero_simulacoes, "\n")
                        
                        # geração dos dados
                        dados_sim <- gerar_dados_simulacao(n = n, p = p, beta = beta, seed = seed + i)
                        
                        # divisão treino/teste
                        split <- split_data(dados_sim, seed = seed + i)
                        treino <- split$train_sample
                        teste  <- split$test_sample
                        
                        formula_surv <-  as.formula(
                          paste0("Surv(time, status)", " ~ ",
                                 paste(setdiff(names(dados_sim), c("time", "status")), collapse = " + "))
                        )
                        
                        # treinamento SVM
                        modelo_svm <- tryCatch({
                          survivalsvm(formula_surv,
                                      data = treino,
                                      type = "regression",
                                      gamma.mu = 1,
                                      kernel = "rbf_kernel")  }, error = function(e) NULL)
                        
                        if (is.null(modelo_svm)) {
                          return(NA_real_) # devolve NA se não treinou
                        }
                        
                        # predição
                        svm_pred_aux <- predict(modelo_svm, newdata = teste)
                        svm_pred <- svm_pred_aux$predicted
                        
                        
                        if (all(is.na(svm_pred))) {
                          return(NA_real_)
                        }
                        
                        # cálculo do C-index
                        cindex_svm <- concordance(Surv(teste$time, teste$status) ~ as.vector(svm_pred))$concordance
                        
                        # retorno da iteração
                        return(cindex_svm)
                      }

stopCluster(cl)

saveRDS(resultados_svm_rbf, "outputs/resultado_simulacoes_svm_rbf.rds")

############# SVM Aditivo

cat("Simulação para SVM Add \n")

# número de núcleos disponíveis
numCores <- detectCores() - 1
cl <- makeCluster(numCores)
registerDoParallel(cl)

# exportar funções e objetos criados no script
clusterExport(cl, c("gerar_dados_simulacao", "split_data", "n", "p", "beta", "seed", "numero_simulacoes"))

# rodar as simulações em paralelo
resultados_svm_add <- foreach(i = 1:numero_simulacoes,
                              .combine = c, # junta os resultados em um vetor
                              .packages = c("survival", "survivalsvm", "randomForestSRC")) %dopar% {
                                
                                msg <- paste("SVM Add: Simulação", i, "de", numero_simulacoes, "\n")
                                
                                # geração dos dados
                                dados_sim <- gerar_dados_simulacao(n = n, p = p, beta = beta, seed = seed + i)
                                
                                # divisão treino/teste
                                split <- split_data(dados_sim, seed = seed + i)
                                treino <- split$train_sample
                                teste  <- split$test_sample
                                
                                formula_surv <-  as.formula(
                                  paste0("Surv(time, status)", " ~ ",
                                         paste(setdiff(names(dados_sim), c("time", "status")), collapse = " + "))
                                )
                                
                                # treinamento SVM
                                modelo_svm <- tryCatch({
                                  survivalsvm(formula_surv,
                                              data = treino,
                                              type = "regression",
                                              gamma.mu = 1,
                                              kernel = "add_kernel")  }, error = function(e) NULL)
                                
                                if (is.null(modelo_svm)) {
                                  return(NA_real_) # devolve NA se não treinou
                                }
                                
                                # predição
                                svm_pred_aux <- predict(modelo_svm, newdata = teste)
                                svm_pred <- svm_pred_aux$predicted
                                
                                
                                if (all(is.na(svm_pred))) {
                                  return(NA_real_)
                                }
                                
                                # cálculo do C-index
                                cindex_svm <- concordance(Surv(teste$time, teste$status) ~ as.vector(svm_pred))$concordance
                                
                                # retorno da iteração
                                return(cindex_svm)
                              }

stopCluster(cl)

saveRDS(resultados_svm_add, "outputs/resultado_simulacoes_svm_add.rds")

############# SVM Linear

cat("Simulação para SVM Linear \n")

# número de núcleos disponíveis
numCores <- detectCores() - 1
cl <- makeCluster(numCores)
registerDoParallel(cl)

# exportar funções e objetos criados no script
clusterExport(cl, c("gerar_dados_simulacao", "split_data", "n", "p", "beta", "seed", "numero_simulacoes"))

# rodar as simulações em paralelo
resultados_svm_lin <- foreach(i = 1:numero_simulacoes,
                              .combine = c, # junta os resultados em um vetor
                              .packages = c("survival", "survivalsvm", "randomForestSRC")) %dopar% {
                                
                                msg <- paste("SVM Linear: Simulação", i, "de", numero_simulacoes, "\n")
                                
                                # geração dos dados
                                dados_sim <- gerar_dados_simulacao(n = n, p = p, beta = beta, seed = seed + i)
                                
                                # divisão treino/teste
                                split <- split_data(dados_sim, seed = seed + i)
                                treino <- split$train_sample
                                teste  <- split$test_sample
                                
                                formula_surv <-  as.formula(
                                  paste0("Surv(time, status)", " ~ ",
                                         paste(setdiff(names(dados_sim), c("time", "status")), collapse = " + "))
                                )
                                
                                # treinamento SVM
                                modelo_svm <- tryCatch({
                                  survivalsvm(formula_surv,
                                              data = treino,
                                              type = "regression",
                                              gamma.mu = 1,
                                              kernel = "lin_kernel")  }, error = function(e) NULL)
                                
                                if (is.null(modelo_svm)) {
                                  return(NA_real_) # devolve NA se não treinou
                                }
                                
                                # predição
                                svm_pred_aux <- predict(modelo_svm, newdata = teste)
                                svm_pred <- svm_pred_aux$predicted
                                
                                
                                if (all(is.na(svm_pred))) {
                                  return(NA_real_)
                                }
                                
                                # cálculo do C-index
                                cindex_svm <- concordance(Surv(teste$time, teste$status) ~ as.vector(svm_pred))$concordance
                                
                                # retorno da iteração
                                return(cindex_svm)
                              }

stopCluster(cl)

saveRDS(resultados_svm_lin, "outputs/resultado_simulacoes_svm_lin.rds")

cat("Os resultados das simulações foram salvos \n")
rm(list = ls())
cat("O ambiente foi limpo \n")




