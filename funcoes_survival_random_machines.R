###################### PACOTES UTILIZADOS

# Define mirror do CRAN (necessário fora do RStudio)
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Lista de pacotes necessários
pacotes <- c("readr", "janitor", "survival", "ggsurvfit", "survminer", 
             "dplyr", "survivalsvm", "caret")

# Função para instalar e carregar pacotes com dependências
instalar_e_carregar <- function(pkg) {
  if (!suppressWarnings(require(pkg, character.only = TRUE))) {
    suppressMessages(install.packages(pkg, dependencies = TRUE))
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# Instalação e carregamento silencioso
invisible(sapply(pacotes, instalar_e_carregar))




################# FUNÇÕES RANDOM MACHINE PARA ANÁLISE DE SOBREVIVÊNCIA

split_data <- function(data, training_ratio = 0.7,
                       test_ratio = 0.3,
                       seed = NULL,
                       time_var = "time",
                       status_var = "status") {
  
  # Verificar se as variáveis de sobrevivência existem
  if(!(time_var %in% names(data)) || !(status_var %in% names(data))) {
    stop("As variáveis de tempo e status devem estar presentes nos dados")
  }
  
  # Verificar se a soma das proporções é 1
  if(sum(training_ratio, test_ratio) != 1) {
    stop("A soma das proporções deve ser igual a 1.")
  }
  
  # Garantir balanceamento de eventos nos conjuntos
  set.seed(seed)
  event_indices <- which(data[[status_var]] == 1)
  nonevent_indices <- which(data[[status_var]] == 0)
  
  # Amostrar eventos e não-eventos separadamente
  train_event <- sample(event_indices, round(length(event_indices) * training_ratio))
  train_nonevent <- sample(nonevent_indices, round(length(nonevent_indices) * training_ratio))
  training_index <- c(train_event, train_nonevent)
  
  # Restante para validação e teste
  remaining_data <- data[-training_index, ]
  remaining_event <- which(remaining_data[[status_var]] == 1)
  remaining_nonevent <- which(remaining_data[[status_var]] == 0)
  
  # Criar os conjuntos
  training_sample <- data[training_index, ]
  test_sample <- remaining_data
  
  # Retornar lista com os conjuntos
  return(list(train_sample = training_sample,
              test_sample = test_sample))
}


survival_random_machines <- function(formula,
                                     train,
                                     test,
                                     boots_number = 100,
                                     cost = 1,
                                     selecionar_variaveis = TRUE,
                                     seed.bootstrap = NULL,
                                     n_cores = NULL,
                                     save = TRUE,
                                     save_path = "modelos_ensemble.rds") {
  
  # Lista de pacotes necessários
  pacotes <- c("survival", "survivalsvm", "foreach", "doParallel", "fs")
  instalar_e_carregar <- function(pkg) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
      suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    } else {
      suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    }
  }
  invisible(sapply(pacotes, instalar_e_carregar))
  
  # Detectar núcleos disponíveis
  if (is.null(n_cores)) n_cores <- max(1, parallel::detectCores() - 2)
  
  # Criar diretório
  dir_create("modelos_amostras_bootstrap")
  
  # Preparar paralelismo
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  
  # Variáveis
  response_var <- all.vars(formula[[2]])
  predictors <- all.vars(formula[[3]])
  p_total <- length(predictors)
  mtry <- floor(sqrt(p_total))
  kernels <- c("lin_kernel", "rbf_kernel", "add_kernel")
  
  print("Treinando modelos iniciais")
  
  # Modelos iniciais para pesos
  modelos_iniciais <- foreach(i = seq_along(kernels), .packages = "survivalsvm") %dopar% {
    survivalsvm(formula, data = train, type = "regression", gamma.mu = cost, kernel = kernels[i])
  }
  names(modelos_iniciais) <- kernels
  
  # Previsões iniciais
  predicted <- foreach(i = seq_along(modelos_iniciais), .packages = "survivalsvm", .combine = 'c') %dopar% {
    list(as.vector(predict(modelos_iniciais[[i]], newdata = test)$predicted))
  }
  names(predicted) <- kernels
  
  print("Calculando C-index dos modelos iniciais")
  
  c_indexes <- foreach(i = seq_along(predicted), .packages = "survival", .combine = 'c') %dopar% {
    concordance(Surv(test[[response_var[1]]], test[[response_var[2]]]) ~ predicted[[i]])$concordance
  }
  names(c_indexes) <- kernels
  
  prob_weights <- c_indexes / sum(c_indexes)
  
  print("Gerando amostras bootstrap")
  
  set.seed(seed.bootstrap)
  boots_index_row_new <- replicate(boots_number, sample(1:nrow(train), nrow(train), replace = TRUE), simplify = FALSE)
  boots_sample <- lapply(boots_index_row_new, function(x) train[x, ])
  out_of_bag <- lapply(boots_index_row_new, function(x) train[-unique(x), ])
  random_kernel <- sample(kernels, boots_number, replace = TRUE, prob = prob_weights)
  
  cat("Treinando nas amostras bootstrap \n")
  
  foreach(b = 1:boots_number, .packages = c("survivalsvm", "survival")) %dopar% {
    kernel_sel <- random_kernel[b]
    data_boot <- boots_sample[[b]]
    
    var_sel <- if(selecionar_variaveis){
      sample(predictors, mtry)
    } else {
      predictors
    }
    
    formula_boot <- as.formula(
      paste0("Surv(", paste(response_var, collapse = ","), ") ~ ",
             paste(var_sel, collapse = " + "))
    )
    
    path_b <- file.path("modelos_amostras_bootstrap", paste0("modelo_bootstrap_", b, ".rds"))
    
      tryCatch({
        modelo_b <- survivalsvm(formula_boot, data = data_boot,
                                type = "regression", gamma.mu = cost, kernel = kernel_sel)
        saveRDS(modelo_b, path_b)
        rm(modelo_b)
        gc()
      }, error = function(e) {
        message(sprintf("Erro no bootstrap %d: %s", b, e$message))
      })
  }
  
  stopCluster(cl)
  
  message("Lendo modelos salvos...")
  model_paths <- dir_ls("modelos_amostras_bootstrap", glob = "*.rds")
  modelos_final <- lapply(model_paths, readRDS)
  
  model_result <- list(
    models = modelos_final,
    boots_sample = boots_sample,
    out_of_bag = out_of_bag,
    kernels = kernels,
    accuracy = c_indexes,
    lambda_values = prob_weights,
    formula = formula
  )
  
  unlink("modelos_amostras_bootstrap", recursive = TRUE)
  
  if (save) {
    saveRDS(model_result, file = file.path(save_path))
  }
  
  class(model_result) <- "rm_surv_model"
  return(model_result)
}



# Calcula os pesos atribuídos a cada modelo do ensemble, retorna o modelo final: ensembles + pesos
final_surv_rm_model <- function(mod, newdata, save = T, save_path = "modelo_random_machines.rds") { 
  # mod é o modelo ensemble (lista de modelos gerados em survival_random_machines) que será usado para fazer as predições
  # newdata é o conjuto de dados teste
  # save: booleano, indica se queremos que o modelo seja salvo
  # save_path: caminho em que o modelo será salvo
  
  # Extrai os nomes das variáveis de tempo e status da fórmula
  time_var <- all.vars(mod$formula)[1]
  status_var <- all.vars(mod$formula)[2]
  
  # Previsão de cada modelo do ensemble para os dados de teste
  predicted <- lapply(mod$models, function(x) {
    as.vector(predict(x, newdata = newdata)$predicted)
  })
  
  # Previsão das amostras fora da amostra (OOB) ---
  # Adaptação: lida com casos de OOB com menos de 2 observações
  predict_oob <- vector("list", length(mod$models)) # lista para guardar previsões OOB
  c_indexes <- numeric(length(mod$models))          # vetor para guardar C-index por modelo
  valid_model_idx <- c()                            # índices dos modelos válidos
  
  for (i in seq_along(mod$models)) {
    oob_data <- mod$out_of_bag[[i]]
    
    # Verifica se o conjunto OOB tem pelo menos 2 observações
    if (nrow(oob_data) < 2) {
      predict_oob[[i]] <- NA
      c_indexes[i] <- NA
      next
    }
    
    # Remove colunas de tempo e status para prever
    X_oob <- as.data.frame(oob_data[, !(names(oob_data) %in% c(time_var, status_var)) ])
    
    # Tenta realizar a predição
    preds <- tryCatch({
      as.vector(predict(mod$models[[i]], newdata = X_oob)$predicted)
    }, error = function(e) {
      NA
    })
    
    # Verifica se a predição foi bem-sucedida
    if (all(is.na(preds))) {
      predict_oob[[i]] <- NA
      c_indexes[i] <- NA
    } else {
      predict_oob[[i]] <- preds
      
      # Cria objeto de sobrevivência
      surv_obj <- survival::Surv(time = oob_data[[time_var]], event = oob_data[[status_var]])
      
      # Calcula concordância (C-index)
      c_indexes[i] <- as.numeric(survival::concordance(surv_obj ~ preds)$concordance)
      
      # Registra índice do modelo válido
      valid_model_idx <- c(valid_model_idx, i)
    }
  }
  
  # Calcula pesos via concordância (C-index), considerando apenas modelos válidos
  valid_c_indexes <- c_indexes[valid_model_idx]
  weights <- valid_c_indexes / sum(valid_c_indexes)
  
  # Combina as predições dos modelos por média ponderada e calcula c-index
  predicted_matrix <- do.call(cbind, predicted[valid_model_idx])
  
  final_prediction <- as.numeric(as.matrix(predicted_matrix) %*% as.numeric(weights))
  
  test_surv <- survival::Surv(time = newdata[[time_var]], 
                              event = newdata[[status_var]])
  cindex_test <- survival::concordance(test_surv ~ final_prediction)$concordance
  
  # Criar objeto do modelo completo para salvar ---
  model_to_save <- list(
    ensemble_models = mod$models,
    bootstrap_samples = mod$boots_sample,
    oob_samples = mod$out_of_bag,
    prediction_rule = list(
      weights = weights,
      aggregation_method = "weighted_average",
      valid_model_idx = valid_model_idx
    ),
    model_formula = mod$formula,
    cindex_oob = mean(valid_c_indexes, na.rm = TRUE),
    cindex_test = cindex_test,
    last_prediction = final_prediction,
    timestamp = Sys.time()
  )
  
  # Salvar modelo completo ---
  if(save){
    saveRDS(model_to_save, file = save_path)
    message("Modelo salvo em: ", save_path)
  }
  
  # 8. Mostrar resumo do modelo ---
  cat("\n=== Modelo Survival Random Machines ===\n")
  cat("Número de modelos no ensemble:", length(mod$models), "\n")
  cat("Modelos válidos usados:", length(valid_model_idx), "\n")
  cat("C-index OOB (média):", round(mean(valid_c_indexes, na.rm = TRUE), 3), "\n")
  cat("C-index no conjunto de teste:", round(cindex_test, 3), "\n")
  cat("Pesos dos modelos válidos:", round(weights, 3), "\n")
  if(save){
  cat("Arquivo do modelo:", save_path, "\n")
  }
  
  # 9. Retornar resultados ---
  return(model_to_save)
}

predict_from_model <- function(saved_model, newdata) {
  # Extrair variáveis de sobrevivência da fórmula
  time_var <- all.vars(saved_model$model_formula[[2]])[1]
  status_var <- all.vars(saved_model$model_formula[[2]])[2]
  
  # Verificar se novas dados têm colunas necessárias
  required_vars <- setdiff(all.vars(saved_model$model_formula[[3]]), 
                           c(time_var, status_var))
  if (!all(required_vars %in% colnames(newdata))) {
    stop("Dados novos não contém todas as variáveis necessárias")
  }
  
  # Verifica se existe o vetor de índices dos modelos válidos
  if (!is.null(saved_model$prediction_rule$valid_model_idx)) {
    valid_idx <- saved_model$prediction_rule$valid_model_idx
  } else {
    valid_idx <- seq_along(saved_model$ensemble_models)
  }
  
  # Extrai índices dos modelos válidos e seus pesos
  weights <- saved_model$prediction_rule$weights
  
  # Subconjunto dos modelos válidos
  valid_models <- saved_model$ensemble_models[valid_idx]
  
  # Fazer previsões com cada modelo válido
  predictions <- lapply(valid_models, function(model) {
    as.vector(predict(model, newdata = newdata)$predicted)
  })
  
  # Monta matriz com as predições válidas
  prediction_matrix <- do.call(cbind, predictions)
  
  # Usa só os pesos correspondentes aos modelos válidos
  weights_valid <- weights
  
  # Calcula predição ponderada
  weighted_prediction <- as.numeric(prediction_matrix %*% weights_valid)
  
  return(weighted_prediction)
}


#################### FUNÇÕES PARA SIMULAÇÃO


gerar_dados_simulacao <- function(n = 300, p=10, prop_censura = 0.3, beta, seed = NULL) {
  # n é o número de observações do conjunto de dados simulados
  # p é o número de variáveis preditoras
  # prop_censura porporção de censura esperada
  # beta são os parâmetros para os dados de simulação
  
  if (!is.null(seed)) {
    set.seed(seed) }
  
  if (length(beta) != p) {
    stop("O vetor beta deve ter comprimento igual ao número de preditores (p).")
  }
  
  # Gerar X_mat: 23 dummies (0/1), cada coluna com proporção aleatória de 1s entre 0.1 e 0.9
  X_mat <- matrix(0, nrow = n, ncol = p)
  for (j in 1:p) {
    prob_1 <- runif(1, 0.1, 0.9)  # proporção de 1s na variável j
    X_mat[, j] <- rbinom(n, size = 1, prob = prob_1)
  }
  colnames(X_mat) <- paste0("X", 1:p)
   
  # taxa da exponencial
  lambda <- exp(X_mat %*% beta)           
  
  # gera os tempos "reais" até o evento de uma exponencial
  tempo_evento <- rexp(n, rate = lambda)
  
  # parâmetro para os tempos censurados
  taxa_censura <- (prop_censura/(1-prop_censura))*lambda
  # gera os tempos de censura também usando a exponencial
  tempo_censura <- rexp(n, rate = taxa_censura)
  
  # define o tempo observado como mínimo entre tempo até o evento e censura (censura à direita)
  tempo_obs <- pmin(tempo_evento, tempo_censura)
  # coluna para definir se o tempo corresponde à falha ou censura (delta)
  status <- as.integer(tempo_evento <= tempo_censura)
  
  # junção de todas as informações em um data frame
  dados <- cbind.data.frame(time = tempo_obs, status = status, X_mat)
  return(dados)
}

simulacoes_rm <- function(n_sim = 10, 
                              n = 300, 
                              p = 10, 
                              beta,
                              seed = 2402,
                              boots_number = 30, 
                              selecionar_variaveis = T,
                              n_cores = NULL) {
  # n_sim: número de simulações
  # n: número de observações por simulação
  # p: número de variáveis preditoras
  # selecionar_variaveis: se TRUE, um conjunto de variáveis será selecionado para cada amostra bootstrap para treinar o modelo 
  # boots_number: número de amostras bootstraps a serem usadas no ensemble
  # n_cores: número de núcleos para o processamento em paralelo
  # save_dir: diretório onde os modelos ensemble de cada simulação serão salvos
  
  resultados <- numeric(n_sim) # vetor que armazenará os resultados
  
  for (i in 1:n_sim) { # iteração para cada simulação
    cat("Simulação", i, "de", n_sim, "\n")
    
    # Geração dos dados
    dados_sim <- gerar_dados_simulacao(n = n, p = p, beta = beta, seed = seed + i)
    
    # Divisão treino/teste
    split <- split_data(dados_sim, seed = seed + i)
    treino <- split$train_sample
    teste <- split$test_sample
    
    # Fórmula
    formula_surv <-  as.formula(
      paste0("Surv(time, status)", " ~ ",
             paste(setdiff(names(dados_sim), c("time", "status")), collapse = " + "))
    )
    
    # Treinamento do ensemble
    modelo_ensemble <- survival_random_machines(
      
      formula = formula_surv,
      train = treino,
      test = teste,
      n_cores = n_cores,
      boots_number = boots_number,
      selecionar_variaveis = selecionar_variaveis,
      seed.bootstrap = seed + i,
      cost = 1,
      save = F
    )
    
    # Aplicação da média ponderada
    modelo_final <- final_surv_rm_model(mod = modelo_ensemble, 
                                        newdata = teste, 
                                        save = F)
    
    # Guardar C-index no teste
    resultados[i] <- modelo_final$cindex_test
    
    if(i==n_sim){
      modelo_exemplo <- modelo_final
    }
    
    rm(modelo_ensemble, modelo_final)
  }
  
  resumo_resultados <- list(
    ultimo_rm = modelo_exemplo,
    lista_cindex = resultados
  )
  
  return(resumo_resultados)
}

