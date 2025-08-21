############## CARREGANDO PACOTES

# Lista completa de pacotes
pacotes <- c("readr", "janitor", "survival", "survminer", 
             "dplyr")

# Define mirror do CRAN (necessário fora do RStudio)
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Função para instalar e carregar pacotes com dependências
instalar_e_carregar <- function(pkg) {
  if (!suppressWarnings(require(pkg, character.only = TRUE))) {
    suppressMessages(install.packages(pkg, dependencies = TRUE))
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# Instalação e carregamento silencioso
invisible(sapply(pacotes, instalar_e_carregar))


############## IMPORTANDO DADOS

dados_aux <- read.csv("data/dataset_clean.csv") |>
  filter(ANODIAG <= 2004)

########### DADOS COM CATEGORIZACAO KAPLAN-MEIER

dados <- dados_aux |>
  clean_names() |>
  mutate(
    death = status_cancer_specific,
    sexo = factor(case_when(
      sexo == 1 ~ "masc",
      sexo == 2 ~ "fem"
    ), levels = c("masc", "fem")),
    cateatend = factor(case_when(
      cateatend == 1 | cateatend == 3 ~ "convenio_ou_particular",
      cateatend == 2 ~ "sus",
      cateatend == 9 ~ "sem_informacao"
    ), levels = c("convenio_ou_particular", "sus", "sem_informacao")),
    diagprev = factor(case_when(
      diagprev == 1 ~ "sem_diag_e_sem_trat",
      diagprev == 2 ~ "com_diag_e_sem_trat"
    ), levels = c("sem_diag_e_sem_trat", "com_diag_e_sem_trat")),
    ecgrup = factor(ecgrup),
    cirurgia = factor(case_when(
      cirurgia == 0 ~ "nao",
      cirurgia == 1 ~ "sim"
    ), levels = c("nao", "sim")),
    hormonio = factor(case_when(
      hormonio == 0 ~ "nao",
      hormonio == 1 ~ "sim"
    ), levels = c("nao", "sim")),
    quimio = factor(case_when(
      quimio == 0 ~ "nao",
      quimio == 1 ~ "sim"
    ), levels = c("nao", "sim")),
    radio = factor(case_when(
      radio == 0 ~ "nao",
      radio == 1 ~ "sim"
    ), levels = c("nao", "sim")),
    outros = factor(case_when(
      outros == 0 ~ "nao",
      outros == 1 ~ "sim"
    ), levels = c("nao", "sim")),
    recnenhum = factor(case_when(
      recnenhum == 0 ~ "nao",
      recnenhum == 1 ~ "sim"
    ), levels = c("nao", "sim")),
    escolari_2 = factor(case_when(
      escolari_2 == 1 ~ "analfabeto",
      escolari_2 == 2 ~ "ens_fund_incompleto",
      escolari_2 == 3 ~ "ens_fund_completo",
      escolari_2 == 4 ~ "ens_medio",
      escolari_2 == 5 ~ "ens_superior"
    ), levels = c("analfabeto", "ens_fund_incompleto", "ens_fund_completo", "ens_medio", "ens_superior")),
    idade_cat = factor(case_when(
      idade <= 49 ~ "0_a_49_anos",
      idade >= 50 & idade <= 74 ~ "50_a_74_anos",
      idade >= 75 ~ "75_anos_mais"
    ), levels = c("0_a_49_anos", "50_a_74_anos", "75_anos_mais")),
    tratcons_cat = factor(case_when(
      tratcons <= 60 ~ "ate_60_dias",
      tratcons > 60 ~ "mais_de_60_dias"
    ), levels = c("ate_60_dias", "mais_de_60_dias"))
  ) |>
  mutate(
    anodiag_cat = factor(
      anodiag,
      levels = c("2000", "2001", "2002", "2003", "2004", "2005")
    ),
    .after = "anodiag"
  ) |>
  mutate( 
    diagtrat_cat = factor(
      ifelse(diagtrat <= 81, "ate_81_dias", "mais_de_81_dias"),
      levels = c("ate_81_dias", "mais_de_81_dias")
    ),
    .after = "diagtrat"
  ) |>
  select(!c(anodiag, diagtrat)) |>
  mutate(
    cateatend = relevel(cateatend, "convenio_ou_particular"),
    cirurgia = relevel(cirurgia, "sim"),
    diagprev = relevel(diagprev, "com_diag_e_sem_trat"),
    diagtrat_cat = relevel(diagtrat_cat, "ate_81_dias"),
    ecgrup = relevel(ecgrup, "I"),
    escolari_2 = relevel(escolari_2, "ens_superior"),
    hormonio = relevel(hormonio, "sim"),
    idade_cat = relevel(idade_cat, "0_a_49_anos"),
    outros = relevel(outros, "sim"),
    quimio = relevel(quimio, "nao"),
    radio = relevel(radio, "nao"),
    recnenhum = relevel(recnenhum, "sim"),
    sexo = relevel(sexo, "fem"),
    tratcons_cat = relevel(tratcons_cat, "mais_de_60_dias")
  ) |>
  select(c(time_years, death, anodiag_cat, cateatend, cirurgia,
           diagprev, diagtrat_cat, ecgrup, escolari_2, hormonio, idade_cat, outros,
           quimio, radio, recnenhum, sexo, tratcons_cat))

write.csv(dados, "data/dados_tratados.csv")

cat("Script finalizado \n")
