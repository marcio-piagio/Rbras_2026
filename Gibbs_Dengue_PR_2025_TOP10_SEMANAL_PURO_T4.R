# ============================================================================
# DENGUE / PARANA / 2025 - MULTINOMIAL E GIBBS EM ESCALA SEMANAL PURA
# REFORMULACAO DE 17/11/2026: SEM PROBABILIDADE DIARIA NA VEROSSIMILHANCA
# ============================================================================
# Indices: r = Regional de Saude, t = SEMANA DE INICIO DOS SINTOMAS,
#          j = SEMANA DE ACOMPANHAMENTO (j=1,2,3,4).
# Y_{r,t,j}: numero de primeiras notificacoes, na semana de acompanhamento j,
#            dos casos cujos sintomas comecaram na coorte semanal t.
# Cada individuo entra em uma unica coorte de inicio t (sem sobreposicao).
# O banco possui datas: usa-se a diferenca de datas APENAS para classificar
# cada registro em uma de quatro semanas. O GIBBS VE SOMENTE Y1,...,Y4.
#
# (Y1,Y2,Y3,Y4,A,B)|N,p,theta ~ Multinomial(N; pi1,...,pi4,piA,piB)
# pi1=p; pi_j=theta*p*(1-p)^(j-1), j=2,3,4;
# piA=(1-p)*(1-theta); piB=theta*(1-p)^4.
# p,theta ~ Beta(1,1); lambda_semana ~ Gamma(2,0.02/7) [forma-taxa];
# N|lambda_semana ~ Poisson(lambda_semana).
#
# AS ESTIMATIVAS DO RELATORIO COM ATRASOS DIARIOS NAO VALEM PARA ESTE MODELO.
# Reestimar e substituir todos os resultados, intervalos e graficos do poster.
# O parametro p varia entre coortes t e regionais r (nao entre j da mesma coorte).
# N/U sao dependentes de pressupostos e de distribuicoes a priori.
# ============================================================================

config <- list(
  arquivo_dbc = "E:/05 - artigo-tese/01-data_base/DENGBR25.dbc",
  arquivo_dbc_2026 = "E:/05 - artigo-tese/01-data_base/DENGBR26.dbc",
  usar_dados_2026 = FALSE,
  corte_2026 = as.Date("2026-03-03"),
  pasta_saida = "E:/banner_rbras/resultados_gibbs_dengue_PR_2025_TOP10_SEMANAL_PURO_T4",
  data_inicio = as.Date("2025-01-01"),
  data_fim_descritiva = as.Date("2025-12-31"),
  data_fim_ajuste = as.Date("2025-05-31"),
  # Para coortes de todo o ano, use data_fim_ajuste=as.Date("2025-12-31").
  # As semanas finais so serao estimadas se houver 28 dias completos de
  # observacao, inclusive quando necessario em arquivo adicional de 2026.
  data_corte = as.Date("2025-12-31"),
  seguimento_max_semanas = 4L,
  somente_coortes_completas = TRUE,  # 7 dias de inicio de sintomas por coorte
  codigos_descartados = c("5", "13"),
  somente_individuais = TRUE,
  excluir_duplicidade_marcada = TRUE,
  n_regionais_ajuste = 10L,  # selecao FIXA pelo total notificado em 2025
  min_notificacoes_ajuste = 30L,
  alerta_min_notificacoes_tardias = 5L,
  iter = 12000L,
  burn = 3000L,
  thin = 10L,
  cadeias = 3L,
  semente = 20261117L,
  # Priori DIRETA da probabilidade de notificacao POR SEMANA.
  ap = 1, bp = 1,       # p_semanal ~ Beta(1,1)
  at = 1, bt = 1,       # theta ~ Beta(1,1)
  # lambda = media de CASOS LATENTES POR COORTE SEMANAL.
  # Ao reescalar a antiga lambda/dia ~ Gamma(2,0.02):
  # lambda_semana=7*lambda_dia ~ Gamma(2,0.02/7); E=700/semana.
  # Essa prior pode ser forte em coortes pequenas: RODAR SENSIBILIDADE.
  al = 2, bl = 0.02/7, # lambda/semana ~ Gamma(forma=2, taxa=0.02/7)
  rodar_sensibilidade = FALSE,
  regional_sensibilidade = 19L,
  semana_sensibilidade = 1L,
  taxas_gamma_sensibilidade = c(0.002, 0.02, 0.2)/7,
  salvar_graficos = TRUE,
  salvar_mapas = TRUE,
  semanas_para_mapas = c(1L, 11L, 21L),
  anos_geobr = 2025L
)

stopifnot(config$seguimento_max_semanas == 4L,
          config$n_regionais_ajuste == 10L,
          config$cadeias >= 2L,
          config$iter > config$burn,
          config$thin >= 1L,
          (config$iter - config$burn) >= config$thin,
          config$data_fim_ajuste <= config$data_fim_descritiva)

pacotes <- c("read.dbc", "geobr", "sf", "dplyr", "ggplot2")
faltam_pacotes <- pacotes[!vapply(pacotes, requireNamespace,
                                FUN.VALUE = logical(1), quietly = TRUE)]
if (length(faltam_pacotes)) {
  stop("Instale os pacotes: install.packages(c(",
       paste(shQuote(faltam_pacotes), collapse = ", "), "))")
}
if (config$usar_dados_2026) {
  if (!file.exists(config$arquivo_dbc_2026))
    stop("Arquivo de 2026 nao encontrado: ", config$arquivo_dbc_2026)
  config$data_corte <- config$corte_2026
  warning("Ao combinar DBC de 2025 e 2026, verifique duplicidades de registros")
}
if (!file.exists(config$arquivo_dbc)) {
  stop("Arquivo de 2025 nao encontrado: ", config$arquivo_dbc)
}
dir.create(config$pasta_saida, recursive = TRUE, showWarnings = FALSE)

nome_regionais <- c(
  "Paranaguá", "Metropolitana (Curitiba)", "Ponta Grossa", "Irati",
  "Guarapuava", "União da Vitória", "Pato Branco", "Francisco Beltrão",
  "Foz do Iguaçu", "Cascavel", "Campo Mourão", "Umuarama", "Cianorte",
  "Paranavaí", "Maringá", "Apucarana", "Londrina", "Cornélio Procópio",
  "Jacarezinho", "Toledo", "Telêmaco Borba", "Ivaiporã"
)
regionais <- data.frame(regional = seq_along(nome_regionais),
                        regional_nome = nome_regionais,
                        stringsAsFactors = FALSE)
regionais$regional_identificacao <- sprintf("%02dª RS - %s",
                                             regionais$regional,
                                             regionais$regional_nome)

salvar_csv <- function(x, nome) {
  utils::write.csv2(x, file.path(config$pasta_saida, nome),
                    row.names = FALSE, fileEncoding = "UTF-8", na = "")
}
parse_data <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXt")) return(as.Date(x))
  z <- trimws(as.character(x))
  z[z %in% c("", "NA", "NULL", "00000000", "0")] <- NA_character_
  d <- rep(as.Date(NA), length(z))
  idx <- !is.na(z) & grepl("^[0-9]{8}$", z)
  d[idx] <- as.Date(z[idx], format = "%Y%m%d")
  idx <- !is.na(z) & grepl("^[0-9]{2}/[0-9]{2}/[0-9]{4}$", z)
  d[idx] <- as.Date(z[idx], format = "%d/%m/%Y")
  idx <- !is.na(z) & grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}", z)
  d[idx] <- as.Date(substr(z[idx], 1L, 10L), format = "%Y-%m-%d")
  d
}
limpa_num <- function(x) gsub("[^0-9]", "", trimws(as.character(x)))
normalizar_nome <- function(z) {
  z <- iconv(as.character(z), to = "ASCII//TRANSLIT")
  z <- tolower(trimws(z))
  z <- gsub("[^a-z0-9 ]", " ", z)
  gsub(" +", " ", z)
}
normalizar_muni <- function(z) {
  z <- limpa_num(z)
  z[is.na(z) | nchar(z) < 6L] <- NA_character_
  substr(z, 1L, 6L)
}
identificar_rs <- function(nome) {
  z <- normalizar_nome(nome)
  if (is.na(z)) return(NA_integer_)
  nms <- normalizar_nome(nome_regionais)
  nms[2L] <- "metropolitana"
  for (i in seq_along(nms)) {
    if (grepl(paste0("(^| )", nms[i], "( |$)"), z)) return(as.integer(i))
  }
  if (grepl("(^| )curitiba( |$)", z)) return(2L)
  prefixo <- regmatches(z, regexpr("^[0-9]{1,2}( |$)", z))
  if (length(prefixo) && nzchar(prefixo)) {
    val <- suppressWarnings(as.integer(trimws(prefixo)))
    if (is.finite(val) && val >= 1L && val <= 22L) return(val)
  }
  NA_integer_
}

# ============================================================================
# 1. MUNICIPIOS DO PARANA -> 22 REGIONAIS DE SAUDE
# ============================================================================
cat("Carregando malha das regioes de saude...\n")
args_geobr <- list(year = config$anos_geobr, code_state = 41,
                   geometry_level = "municipality", simplified = FALSE)
if ("showProgress" %in% names(formals(geobr::read_health_region))) {
  args_geobr$showProgress <- FALSE
}
malha <- do.call(geobr::read_health_region, args_geobr)
if (!all(c("code_muni", "name_health_region") %in% names(malha))) {
  stop("Malha geobr nao trouxe code_muni e name_health_region")
}
malha$municipio_6 <- normalizar_muni(malha$code_muni)
malha$regional_num <- vapply(malha$name_health_region,
                             identificar_rs, integer(1))
municipio_para_rs <- unique(sf::st_drop_geometry(malha)[,
                           c("municipio_6", "regional_num")])
esperados_rs <- c(7,29,12,9,20,9,15,27,9,25,25,21,11,28,30,17,21,21,22,18,7,16)
if (nrow(municipio_para_rs) != 399L ||
    anyNA(municipio_para_rs) ||
    anyDuplicated(municipio_para_rs$municipio_6) ||
    length(unique(municipio_para_rs$regional_num)) != 22L ||
    !identical(as.integer(tabulate(municipio_para_rs$regional_num, nbins = 22L)),
               as.integer(esperados_rs))) {
  stop("Correspondencia municipio/regional invalida: audite a malha geobr")
}
mapa_regional <- malha |>
  dplyr::group_by(regional_num) |>
  dplyr::summarise(do_union = TRUE, .groups = "drop")
mapa_regional <- sf::st_make_valid(mapa_regional)
rm(malha)

# ============================================================================
# 2. LEITURA E FILTROS DOS MICRODADOS SINAN
# ============================================================================
entradas <- c(config$arquivo_dbc,
              if (config$usar_dados_2026) config$arquivo_dbc_2026)
base <- dplyr::bind_rows(lapply(entradas, function(arquivo) {
  cat("Lendo ", arquivo, "...\n", sep = "")
  as.data.frame(read.dbc::read.dbc(arquivo), stringsAsFactors = FALSE)
}))
n_bruto <- nrow(base)
necessarias <- c("DT_SIN_PRI", "DT_NOTIFIC", "ID_MN_RESI")
if (length(setdiff(necessarias, names(base)))) {
  stop("Faltam campos: ",
       paste(setdiff(necessarias, names(base)), collapse = ", "))
}

# A Regional e do MUNICIPIO DE RESIDENCIA, nao de notificacao.
res_pr <- !is.na(normalizar_muni(base$ID_MN_RESI)) &
  grepl("^41", normalizar_muni(base$ID_MN_RESI))
base <- base[res_pr, , drop = FALSE]
n_res_pr <- nrow(base)
if (!n_res_pr) stop("Nenhum registro residente no Parana")

if (config$somente_individuais) {
  if (!"TP_NOT" %in% names(base)) stop("Campo TP_NOT ausente")
  base <- base[limpa_num(base$TP_NOT) == "2", , drop = FALSE]
}
if (config$excluir_duplicidade_marcada && "NDUPLIC_N" %in% names(base)) {
  z <- limpa_num(base$NDUPLIC_N)
  base <- base[is.na(z) | z != "2", , drop = FALSE]
}
if ("CLASSI_FIN" %in% names(base)) {
  classe <- trimws(as.character(base$CLASSI_FIN))
  manter <- is.na(classe) | classe == "" |
    !(classe %in% config$codigos_descartados)
  base <- base[manter, , drop = FALSE]
} else {
  warning("CLASSI_FIN ausente: filtro de descartados nao foi aplicado")
}
n_pos_filtros <- nrow(base)
base$municipio_6 <- normalizar_muni(base$ID_MN_RESI)
base$regional <- municipio_para_rs$regional_num[
  match(base$municipio_6, municipio_para_rs$municipio_6)]
base$dt_sint <- parse_data(base$DT_SIN_PRI)
base$dt_not <- parse_data(base$DT_NOTIFIC)
base <- base[!is.na(base$dt_sint) &
             base$dt_sint >= config$data_inicio &
             base$dt_sint <= config$data_fim_descritiva, , drop = FALSE]
n_sint_2025 <- nrow(base)
n_sem_regional <- sum(is.na(base$regional))
base <- base[!is.na(base$regional), , drop = FALSE]
# Datas usadas apenas para alocar o registro na SEMANA de acompanhamento.
# Nao usamos atrasos em dias na verossimilhanca nem no amostrador Gibbs.
base$semana_atraso <- as.integer(base$dt_not - base$dt_sint) %/% 7L + 1L

# ============================================================================
# 3. SEMANAS CIVIS ANCORADAS EM 1o DE JANEIRO DE 2025
# ============================================================================
inicio_ano <- config$data_inicio
fim_ano <- config$data_fim_descritiva
n_semanas <- ceiling((as.integer(fim_ano - inicio_ano) + 1L)/7)
semanas <- data.frame(semana = seq_len(n_semanas))
semanas$inicio <- inicio_ano + 7L * (semanas$semana - 1L)
semanas$fim <- pmin(fim_ano, semanas$inicio + 6L)
semanas$n_dias <- as.integer(semanas$fim - semanas$inicio) + 1L
semanas$completa <- semanas$n_dias == 7L

# INFERENCIA POR COORTE: uma unica semana de sintomas por ajuste.
# O acompanhamento de cada data individual dura ate 27 dias subsequentes.
# A exigencia abaixo impede estimar com acompanhamento incompleto.
coortes <- semanas[
  semanas$inicio <= config$data_fim_ajuste &
  semanas$fim <= config$data_fim_ajuste &
  (!config$somente_coortes_completas | semanas$completa) &
  (semanas$fim + (7L*config$seguimento_max_semanas - 1L) <=
     config$data_corte),
  , drop = FALSE
]
if (!nrow(coortes)) stop("Nenhuma coorte semanal com seguimento completo")
coortes$data_limite_notificacao <- coortes$fim +
  7L*config$seguimento_max_semanas - 1L

sem_elegiveis_data <- semanas$semana[
  semanas$inicio <= config$data_fim_ajuste &
  semanas$fim <= config$data_fim_ajuste &
  (!config$somente_coortes_completas | semanas$completa)]
cat("DESCRITIVA: ", nrow(semanas), " semanas.\n", sep = "")
cat("INFERENCIA: ", nrow(coortes), " coortes semanais individuais; semanas ",
    paste(range(coortes$semana), collapse = " a "), ".\n", sep = "")
if (length(setdiff(sem_elegiveis_data, coortes$semana))) {
  warning("Semanas finais excluidas por seguimento incompleto: ",
          paste(setdiff(sem_elegiveis_data, coortes$semana), collapse = ", "))
}

base$semana_sint <- as.integer(base$dt_sint - inicio_ano) %/% 7L + 1L
base$semana_not <- ifelse(is.na(base$dt_not), NA_integer_,
                          as.integer(base$dt_not - inicio_ano) %/% 7L + 1L)

# DESCRITIVA DO ANO INTEIRO, INDEPENDENTE DOS AJUSTES.
contagens <- expand.grid(regional = regionais$regional,
                         semana = semanas$semana, KEEP.OUT.ATTRS = FALSE)
contagens <- merge(contagens, semanas, by = "semana", all.x = TRUE, sort = FALSE)
contagens$inicio_sintomas <- 0L
contagens$notificacoes <- 0L
for (r in regionais$regional) {
  x <- base[base$regional == r, , drop = FALSE]
  freq_sint <- tabulate(x$semana_sint, nbins = n_semanas)
  x <- x[!is.na(x$semana_not) & !is.na(x$semana_atraso) &
         x$semana_atraso >= 1L & x$dt_not <= fim_ano &
         x$semana_not %in% semanas$semana, , drop = FALSE]
  freq_not <- tabulate(x$semana_not, nbins = n_semanas)
  idx <- which(contagens$regional == r)
  contagens$inicio_sintomas[idx] <- freq_sint[contagens$semana[idx]]
  contagens$notificacoes[idx] <- freq_not[contagens$semana[idx]]
}
contagens <- contagens[order(contagens$regional, contagens$semana), ]
contagens$regional_nome <- nome_regionais[contagens$regional]
contagens$regional_identificacao <- regionais$regional_identificacao[
  contagens$regional]
ranking_rs <- stats::aggregate(notificacoes ~ regional, data = contagens, FUN = sum)
ranking_rs <- merge(regionais, ranking_rs, by = "regional", sort = FALSE)
ranking_rs <- ranking_rs[order(-ranking_rs$notificacoes, ranking_rs$regional), ]
ranking_rs$posicao <- seq_len(nrow(ranking_rs))
ranking_rs$selecionada <- ranking_rs$posicao <= config$n_regionais_ajuste
regionais_modelo <- ranking_rs[ranking_rs$selecionada, , drop = FALSE]
ids_modelo <- regionais_modelo$regional
contagens_top10 <- contagens[contagens$regional %in% ids_modelo, , drop = FALSE]

cat("\nTOP 10 regionais segundo notificacoes observadas em 2025:\n")
print(regionais_modelo[,c("posicao", "regional_identificacao", "notificacoes")],
      row.names = FALSE)

anual_observado <- stats::aggregate(cbind(inicio_sintomas, notificacoes) ~
                  regional + regional_nome, data = contagens, FUN = sum)
salvar_csv(semanas, "01_semanas_2025.csv")
salvar_csv(coortes, "02_COORTES_SEMANAIS_2025.csv")
salvar_csv(regionais, "03_nomes_regionais.csv")
salvar_csv(ranking_rs, "03A_RANKING_22_REGIONAIS_TOP10.csv")
salvar_csv(contagens, "04_CONTAGENS_SEMANAIS_22_REGIONAIS_DESCRITIVAS.csv")
salvar_csv(contagens_top10, "04B_CONTAGENS_SEMANAIS_TOP10_DESCRITIVAS.csv")
salvar_csv(anual_observado, "09_contagens_anuais_22_RS_descritivo.csv")

# ============================================================================
# 4. VETORES OBSERVADOS POR COORTE DE SINTOMAS
# ============================================================================
# Um registro entra UMA vez: na semana correspondente a DT_SIN_PRI.
# Categoria j = semana de acompanhamento ate a primeira notificacao.
base_inferencia <- base[
  base$regional %in% ids_modelo &
  base$semana_sint %in% coortes$semana &
  !is.na(base$semana_atraso) & base$semana_atraso >= 1L &
  base$semana_atraso <= config$seguimento_max_semanas &
  !is.na(base$dt_not) & base$dt_not <= config$data_corte,
  , drop = FALSE
]
stopifnot(all(base_inferencia$semana_sint %in% coortes$semana))
chaves <- paste(base_inferencia$regional, base_inferencia$semana_sint,
                sep = "_")
semanas_not_por_coorte <- split(as.integer(base_inferencia$semana_atraso), chaves)

vetores_coortes <- expand.grid(
  regional = ids_modelo, semana = coortes$semana,
  j = seq_len(config$seguimento_max_semanas), KEEP.OUT.ATTRS = FALSE
)
vetores_coortes$y <- 0L
for (r in ids_modelo) {
  for (t in coortes$semana) {
    semanas_not <- semanas_not_por_coorte[[paste(r, t, sep = "_")]]
    freq <- if (length(semanas_not))
      tabulate(semanas_not, nbins = 4L) else integer(4)
    idx <- which(vetores_coortes$regional == r &
                 vetores_coortes$semana == t)
    vetores_coortes$y[idx] <- freq[vetores_coortes$j[idx]]
  }
}
vetores_coortes <- vetores_coortes[
  order(vetores_coortes$regional, vetores_coortes$semana, vetores_coortes$j), ]
vetores_coortes$regional_identificacao <-
  regionais$regional_identificacao[vetores_coortes$regional]
vetores_coortes$semana_atraso <- paste0("Semana de acompanhamento ",vetores_coortes$j)
salvar_csv(vetores_coortes, "05_y_j_por_coorte_semanal.csv")

y_largo <- reshape(vetores_coortes[,c("regional", "semana", "j", "y")],
                   idvar = c("regional", "semana"),
                   timevar = "j", direction = "wide")
names(y_largo) <- sub("^y\\.", "y", names(y_largo))
y_largo <- merge(y_largo, coortes[,c("semana", "inicio", "fim", "n_dias")],
                 by = "semana", sort = FALSE)
y_largo$C_observado <- rowSums(y_largo[,paste0("y",1:4)])
y_largo$regional_identificacao <-
  regionais$regional_identificacao[y_largo$regional]
y_largo <- y_largo[order(y_largo$regional, y_largo$semana),]
salvar_csv(y_largo, "06_CONTAGENS_COORTES_SEMANAIS_Y1_Y4.csv")

qualidade <- data.frame(
  indicador = c("registros_brutos", "residentes_PR", "apos_filtros_TP_NOT_DUP_CLASSI",
                "sintomas_2025", "sem_regional_residencia",
                "notificacoes_distintas_utilizadas_ajuste_semana",
                "semanas_descritivas_2025", "coortes_semanais_ajustaveis",
                "coortes_semanais_completas", "combinacoes_regional_semana_previstas",
                "primeira_coorte", "ultima_coorte", "corte_notificacoes"),
  valor = as.character(c(n_bruto, n_res_pr, n_pos_filtros, n_sint_2025,
                         n_sem_regional, nrow(base_inferencia), nrow(semanas),
                         nrow(coortes), sum(coortes$completa),
                         nrow(coortes)*length(ids_modelo),
                         as.character(min(coortes$inicio)),
                         as.character(max(coortes$fim)),
                         as.character(config$data_corte)))
)
salvar_csv(qualidade, "00_auditoria_filtros_COORTES_SEMANAIS.csv")

# ============================================================================
# 5. GIBBS EXATO EM ESCALA SEMANAL POR REGIONAL E COORTE
# ============================================================================
# Para cada (r,t), observado SOMENTE y=(Y1,Y2,Y3,Y4).
# C=Y1+Y2+Y3+Y4; W=Y2+2*Y3+3*Y4; S=Y2+Y3+Y4.
# N|lambda ~ Poisson(lambda), lambda em casos esperados por SEMANA.
# p~Beta(ap,bp); theta~Beta(at,bt); lambda~Gamma(al,bl) [forma-taxa].
# A|- ~ Poi(lambda*(1-p)*(1-theta))
# B|- ~ Poi(lambda*theta*(1-p)^4)
# N=C+A+B
# p|- ~ Beta(ap+C, bp+W+A+4*B)
# theta|- ~ Beta(at+S+B, bt+A)
# lambda|- ~ Gamma(al+N, bl+1)
# Sem h diario, sem F diario, sem fator 7 nas condicionais.
# IMPORTANTE: Yj descreve semana de primeira notificacao APOS inicio
# individual dos sintomas; nao implica semana civil de calendario.
# ============================================================================

resumir <- function(x, prefixo) {
  if (!length(x) || all(!is.finite(x))) {
    ans <- rep(NA_real_, 4)
  } else {
    qs <- unname(stats::quantile(x, c(0.025, 0.5, 0.975), na.rm = TRUE))
    ans <- c(mean(x, na.rm = TRUE), qs)
  }
  names(ans) <- paste0(prefixo, c("_media", "_li95", "_mediana", "_ls95"))
  ans
}
rhat_basico <- function(lista) {
  m <- length(lista)
  n <- min(vapply(lista, length, integer(1)))
  if (m < 2L || n < 4L) return(NA_real_)
  x <- sapply(lista, function(z) z[seq_len(n)])
  W <- mean(apply(x, 2L, stats::var))
  B <- n * stats::var(colMeans(x))
  if (!is.finite(W) || W <= 0) return(NA_real_)
  sqrt(((n-1)/n*W+B/n)/W)
}

probs_multinomiais <- function(p, theta) {
  q <- 1-p
  cbind(pi1 = p,
        pi2 = theta*p*q,
        pi3 = theta*p*q^2,
        pi4 = theta*p*q^3,
        piA = q*(1-theta),
        piB = theta*q^4)
}
# Valide a soma = 1 com parametros genericos.
teste_p <- c(.02, .35, .7, .98)
teste_t <- c(.1, .5, .8, .99)
stopifnot(max(abs(rowSums(probs_multinomiais(teste_p, teste_t))-1)) < 1e-12)

# Cada observacao disponivel possui quatro semanas completas.
gibbs_coorte <- function(y, cfg, semente_base) {
  if (length(y) != 4L || anyNA(y) || any(!is.finite(y)) ||
      any(y < 0L) || any(y != floor(y)))
    stop("y deve ser vetor de quatro contagens semanais inteiras >= 0")
  C <- sum(y)
  W <- sum((0:3)*y)      # Semanas completas de espera (escala semanal)
  S <- sum(y[2:4])       # Notificacoes da 2a a 4a semanas
  if (C < cfg$min_notificacoes_ajuste) return(NULL)
  nsave <- sum((seq_len(cfg$iter) > cfg$burn) &
                 ((seq_len(cfg$iter) - cfg$burn) %% cfg$thin == 0L))
  vars <- c("N", "U", "p", "theta", "lambda", "A", "B", "cobertura")
  cadeias <- vector("list", cfg$cadeias)
  for (ch in seq_len(cfg$cadeias)) {
    set.seed(as.integer(semente_base + ch))
    p <- 0.7
    theta <- 0.5
    lambda <- max(C/0.7, 1)
    saida <- matrix(NA_real_, nsave, length(vars),
                    dimnames = list(NULL, vars))
    pos <- 0L
    for (it in seq_len(cfg$iter)) {
      A <- stats::rpois(1L, lambda*(1-p)*(1-theta))
      B <- stats::rpois(1L, lambda*theta*(1-p)^4)
      N <- C+A+B
      p <- stats::rbeta(1L, cfg$ap+C, cfg$bp+W+A+4*B)
      theta <- stats::rbeta(1L, cfg$at+S+B, cfg$bt+A)
      lambda <- stats::rgamma(1L, shape=cfg$al+N, rate=cfg$bl+1)
      if (it > cfg$burn && ((it-cfg$burn) %% cfg$thin == 0L)) {
        pos <- pos + 1L
        cobertura <- 1 - ((1-p)*(1-theta) + theta*(1-p)^4)
        saida[pos,] <- c(N, A+B, p, theta, lambda, A, B, cobertura)
      }
    }
    cadeias[[ch]] <- saida
  }
  amostras <- as.data.frame(do.call(rbind, cadeias))
  prs <- probs_multinomiais(amostras$p, amostras$theta)
  amostras <- cbind(amostras, as.data.frame(prs))
  stopifnot(all(abs(rowSums(prs)-1)<1e-8))
  rhat <- c(
    rhat_p = rhat_basico(lapply(cadeias, function(m) m[,"p"])),
    rhat_theta = rhat_basico(lapply(cadeias, function(m) m[,"theta"])),
    rhat_lambda = rhat_basico(lapply(cadeias, function(m) m[,"lambda"])),
    rhat_N = rhat_basico(lapply(cadeias, function(m) m[,"N"])),
    rhat_U = rhat_basico(lapply(cadeias, function(m) m[,"U"]))
  )
  # ESS da cadeia MCMC requer pacote 'coda'. NA se indisponivel.
  ess <- c(ess_p=NA_real_, ess_theta=NA_real_, ess_N=NA_real_, ess_U=NA_real_)
  if (requireNamespace("coda", quietly=TRUE)) {
    for (nm in c("p", "theta", "N", "U")) {
      ess[paste0("ess_",nm)] <- tryCatch({
        z <- coda::mcmc.list(lapply(cadeias, function(a)
          coda::mcmc(a[,nm])))
        unname(coda::effectiveSize(z))
      }, error=function(e) NA_real_)
    }
  }
  list(amostras=amostras, rhat=rhat, ess=ess, obs=C, tardios=S,
       esperas_semanais=W)
}

# ============================================================================
# 6. 10 REGIONAIS x COORTES SEMANAIS: MODELOS INDEPENDENTES
# ============================================================================
resumos <- vector("list", length(ids_modelo)*nrow(coortes))
posteriores <- list()
cont <- 0L
n_samp <- config$cadeias * sum((seq_len(config$iter) > config$burn) &
                 ((seq_len(config$iter)-config$burn) %% config$thin == 0L))

for (r in ids_modelo) {
  for (i in seq_len(nrow(coortes))) {
    t <- coortes$semana[i]
    cont <- cont+1L
    ini <- coortes$inicio[i]
    fim <- coortes$fim[i]
    n_dias <- coortes$n_dias[i] # apenas metadado da coorte, NAO entra no Gibbs
    chave <- paste(r,t,sep="_")
    semanas_not <- semanas_not_por_coorte[[chave]]
    if (is.null(semanas_not)) semanas_not <- integer(0)
    y <- tabulate(semanas_not, nbins=4L)
    C <- sum(y)
    S <- sum(y[2:4])
    status <- if (C >= config$min_notificacoes_ajuste) "ajustado" else
      "sem_informacao_suficiente"
    cat(sprintf("%3d/%3d - RS %02d, coorte semanal %02d (%s a %s), C=%d\n",
                cont,length(resumos),r,t,ini,fim,C))
    linha <- data.frame(
      regional = r, regional_nome = nome_regionais[r],
      regional_identificacao = regionais$regional_identificacao[r],
      semana = t, data_inicio = as.character(ini),
      data_fim = as.character(fim), n_dias = n_dias,
      C_observado = C, notificacoes_apos_7dias = S,
      y1 = y[1L], y2 = y[2L], y3 = y[3L], y4 = y[4L],
      modelo_escala = "semanal_puro",
      status = status, stringsAsFactors = FALSE)
    # Criar todas as colunas em branco para preservar estrutura do CSV.
    param <- c("p", "theta", "lambda", "N", "U", "taxa_sub",
               "cobertura", paste0("pi",1:4), "piA", "piB")
    for (nm in param) {
      for (suf in c("media", "li95", "mediana", "ls95")) {
        linha[[paste0(nm,"_",suf)]] <- NA_real_
      }
    }
    for (nm in c("rhat_p","rhat_theta","rhat_lambda","rhat_N","rhat_U",
                 "ess_p", "ess_theta", "ess_N", "ess_U")) {
      linha[[nm]] <- NA_real_
    }
    linha$alerta <- ""

    if (status == "ajustado") {
      fit <- gibbs_coorte(y, config,
                          semente_base = config$semente+10000L*r+10L*t)
      a <- fit$amostras
      a$taxa_sub <- ifelse(a$N > 0, 100*a$U/a$N, NA_real_)
      vals <- unlist(lapply(param, function(nm) resumir(a[[nm]], nm)))
      for (nm in names(vals)) linha[[nm]] <- unname(vals[[nm]])
      for (nm in names(fit$rhat)) linha[[nm]] <- unname(fit$rhat[[nm]])
      for (nm in names(fit$ess)) linha[[nm]] <- unname(fit$ess[[nm]])
      avisos <- character(0)
      if (S < config$alerta_min_notificacoes_tardias) {
        avisos <- c(avisos, "poucas_notificacoes_tardias_theta_fragil")
      }
      if (any(!is.finite(fit$rhat)) || any(fit$rhat > 1.05,na.rm=TRUE)) {
        avisos <- c(avisos, "verificar_convergencia")
      }
      if (any(is.finite(fit$ess) & fit$ess < 400)) {
        avisos <- c(avisos, "ESS_baixo")
      }
      if (length(avisos)) linha$alerta <- paste(avisos, collapse=";")
      posteriores[[chave]] <- a[, c("N","U","p","theta","cobertura",
                                     "taxa_sub","pi1","pi2","pi3","pi4","piA","piB")]
    }
    resumos[[cont]] <- linha
  }
}
res_semanais <- do.call(rbind, resumos)
res_semanais <- res_semanais[order(res_semanais$regional,res_semanais$semana),]
salvar_csv(res_semanais,"07_posteriores_TOP10_COORTES_SEMANAIS_T4.csv")

# Formato longo: pi_{t,r,j} realmente distinto por SEMANA e REGIONAL.
pi_longo <- do.call(rbind,lapply(c("pi1","pi2","pi3","pi4","piA","piB"),
    function(nm) data.frame(
      regional = res_semanais$regional,
      regional_identificacao = res_semanais$regional_identificacao,
      semana = res_semanais$semana,
      data_inicio = res_semanais$data_inicio,
      categoria = nm,
      media = res_semanais[[paste0(nm,"_media")]],
      li95 = res_semanais[[paste0(nm,"_li95")]],
      mediana = res_semanais[[paste0(nm,"_mediana")]],
      ls95 = res_semanais[[paste0(nm,"_ls95")]]
    )))
salvar_csv(pi_longo,"07A_PROBABILIDADES_MULTINOMIAIS_POR_SEMANA.csv")

# ============================================================================
# 7. AGREGACAO ENTRE REGIONAIS NA MESMA SEMANA
#    Diferentemente do script de janelas moveis, os N/U de SEMANAS
#    DISTINTAS sao aditivos, pois cada pessoa pertence a apenas uma coorte.
#    Os agregados consideram somente combinacoes ajustadas.
# ============================================================================
res_agregado <- vector("list",nrow(coortes))
for (i in seq_len(nrow(coortes))) {
  t <- coortes$semana[i]
  tab <- res_semanais[res_semanais$semana==t,,drop=FALSE]
  ks <- paste(ids_modelo,t,sep="_")
  ks <- ks[ks %in% names(posteriores)]
  C_total <- sum(tab$C_observado)
  C_modeladas <- sum(tab$C_observado[tab$status=="ajustado"])
  linha <- data.frame(
    semana = t, data_inicio = as.character(coortes$inicio[i]),
    data_fim = as.character(coortes$fim[i]),
    regionais_ajustadas = length(ks), total_regionais=length(ids_modelo),
    C_observado=C_total, C_modeladas=C_modeladas,
    C_sem_ajuste=C_total-C_modeladas,
    status=if(length(ks)==length(ids_modelo)) "10_regionais_ajustadas"
      else "parcial_TOP10", stringsAsFactors=FALSE)
  for (nm in c("N","U","taxa_sub","p_ponderado_N","cobertura")) {
    for (suf in c("media","li95","mediana","ls95")) {
      linha[[paste0(nm,"_",suf)]] <- NA_real_
    }
  }
  if(length(ks)) {
    Ns <- Us <- ps <- numeric(n_samp)
    for(k in ks) {
      a <- posteriores[[k]]
      stopifnot(nrow(a)==n_samp)
      Ns <- Ns+a$N
      Us <- Us+a$U
      ps <- ps+a$N*a$p
    }
    medidas <- c(resumir(Ns,"N"),resumir(Us,"U"),
                 resumir(ifelse(Ns>0,100*Us/Ns,NA_real_),"taxa_sub"),
                 resumir(ifelse(Ns>0,ps/Ns,NA_real_),"p_ponderado_N"),
                 resumir(ifelse(Ns>0,(Ns-Us)/Ns,NA_real_),"cobertura"))
    for(nm in names(medidas)) linha[[nm]] <- unname(medidas[[nm]])
  }
  res_agregado[[i]] <- linha
}
res_agregado <- do.call(rbind,res_agregado)
salvar_csv(res_agregado,"08_resumo_TOP10_COORTES_SEMANAIS_T4.csv")

# Agregado de TODO o periodo de coortes, sem duplicar uma semana.
# Apenas os ajustes efetivamente realizados entram nesse agregado.
ks_todos <- names(posteriores)
if(length(ks_todos)) {
  Nperiodo <- Uperiodo <- numeric(n_samp)
  for(k in ks_todos) {
    Nperiodo <- Nperiodo+posteriores[[k]]$N
    Uperiodo <- Uperiodo+posteriores[[k]]$U
  }
  lin_total <- data.frame(
    coortes_incluidas=nrow(coortes),
    combinacoes_possiveis=nrow(coortes)*length(ids_modelo),
    combinacoes_ajustadas=length(ks_todos),
    C_observado_modelado=sum(res_semanais$C_observado[
      res_semanais$status=="ajustado"]),
    inicio=as.character(min(coortes$inicio)),
    fim=as.character(max(coortes$fim)),
    status=if(length(ks_todos)==nrow(coortes)*length(ids_modelo))
      "periodo_TOP10_completo" else "periodo_TOP10_parcial",
    as.list(c(resumir(Nperiodo,"N"), resumir(Uperiodo,"U"),
              resumir(ifelse(Nperiodo>0,100*Uperiodo/Nperiodo,NA_real_),
                      "taxa_sub"))), check.names=FALSE)
  salvar_csv(lin_total,"09A_TOTAL_PERIODO_COORTES_NAO_SOBREPOSTAS.csv")
}

# ============================================================================
# 8. GRAFICOS: DESCRITIVA ANUAL E INFERENCIA SEMANAL
# ============================================================================
if(config$salvar_graficos) {
  gg <- ggplot2::ggplot
  tema <- ggplot2::theme_minimal(base_size=11)
  dados <- res_semanais
  dados$regional_identificacao <- factor(dados$regional_identificacao,
    levels=regionais_modelo$regional_identificacao)
  dados$segmento <- ave(as.integer(is.na(dados$p_mediana)),dados$regional,
                       FUN=cumsum)
  dados$grupo <- interaction(dados$regional_identificacao,dados$segmento,
                             drop=TRUE)
  salvar_plot <- function(g,nome,w=13,h=10) {
    ggplot2::ggsave(file.path(config$pasta_saida,nome),g,
                    width=w,height=h,device="pdf")
  }
  eixo <- ggplot2::scale_x_continuous(breaks=sort(unique(c(
    min(coortes$semana),seq(1,max(coortes$semana),by=4),max(coortes$semana)))))
  gp <- gg(dados,ggplot2::aes(x=semana,y=p_mediana,group=grupo)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin=p_li95,ymax=p_ls95),
                         fill="#457B9D",alpha=.19,na.rm=TRUE) +
    ggplot2::geom_line(color="#255072",linewidth=.65,na.rm=TRUE) +
    ggplot2::geom_point(color="#255072",size=.8,na.rm=TRUE) +
    ggplot2::facet_wrap(~regional_identificacao,ncol=2,drop=FALSE) +
    eixo+ggplot2::coord_cartesian(ylim=c(0,1)) +
    ggplot2::labs(title="Probabilidade de notificacao por coorte semanal",
      subtitle="Uma coorte por semana de inicio dos sintomas; ICr95%",
      x="Semana de inicio dos sintomas (1=01 a 07/01/2025)",y="p semanal")+tema
  salvar_plot(gp,"10A_p_COORTES_SEMANAIS_TOP10_paineis.pdf",13,15)

  gs <- gg(dados,ggplot2::aes(x=semana,y=taxa_sub_mediana,group=grupo)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin=taxa_sub_li95,
                                     ymax=taxa_sub_ls95),
                         fill="#B15767",alpha=.18,na.rm=TRUE) +
    ggplot2::geom_line(color="#AE4854",linewidth=.65,na.rm=TRUE) +
    ggplot2::facet_wrap(~regional_identificacao,ncol=2,drop=FALSE) +
    eixo+ggplot2::coord_cartesian(ylim=c(0,100)) +
    ggplot2::labs(title="Nao captacao estimada por coorte semanal",
      subtitle="Quatro semanas de seguimento por coorte; ICr95%",
      x="Semana de inicio dos sintomas",y="U/N (%)")+tema
  salvar_plot(gs,"10B_subnotificacao_COORTES_SEMANAIS_TOP10.pdf",13,15)

  gc <- gg(dados,ggplot2::aes(x=semana,y=regional_identificacao,
                              fill=taxa_sub_mediana))+
    ggplot2::geom_tile(color="white",linewidth=.15)+
    ggplot2::scale_fill_gradient(low="#FDEBE8",high="#AE4854",
                                na.value="grey85")+
    eixo+ggplot2::labs(title="Nao captacao estimada: semana x Regional",
       subtitle="Cinza = semana nao ajustada",x="Coorte semanal",
       y="Regional de Saude",fill="U/N (%)")+tema
  salvar_plot(gc,"10C_heatmap_SUB_COORTES_SEMANAIS.pdf",14,6)

  gh <- gg(dados,ggplot2::aes(x=semana,y=regional_identificacao,
                              fill=p_mediana))+
    ggplot2::geom_tile(color="white",linewidth=.15)+
    ggplot2::scale_fill_gradient(low="#E4EFF7",high="#255072",
                                limits=c(0,1),na.value="grey85")+
    eixo+ggplot2::labs(title="Probabilidade p: semana x Regional",
       subtitle="Coortes semanais independentes; cinza = sem ajuste",
       x="Coorte semanal",y="Regional de Saude",fill="p")+tema
  salvar_plot(gh,"10D_heatmap_P_COORTES_SEMANAIS.pdf",14,6)

  gcomp <- gg(dados,ggplot2::aes(x=semana,y=p_mediana,
          color=regional_identificacao,group=grupo))+
    ggplot2::geom_line(linewidth=.8,na.rm=TRUE)+eixo+
    ggplot2::coord_cartesian(ylim=c(0,1))+
    ggplot2::labs(title="Comparacao de p entre as dez Regionais",
      x="Coorte semanal",y="Mediana posterior de p",color="Regional")+tema
  salvar_plot(gcomp,"10E_p_COORTES_SEMANAIS_COMPARATIVO.pdf",14,7)

  for(v in c("N","U","theta")) {
    nm <- paste0(v,"_mediana")
    baixo <- paste0(v,"_li95")
    alto <- paste0(v,"_ls95")
    g <- gg(dados,ggplot2::aes(x=semana,y=.data[[nm]],group=grupo))+
      ggplot2::geom_ribbon(ggplot2::aes(ymin=.data[[baixo]],
                                      ymax=.data[[alto]]),
                           fill="#457B9D",alpha=.17,na.rm=TRUE)+
      ggplot2::geom_line(color="#255072",linewidth=.6,na.rm=TRUE)+
      ggplot2::facet_wrap(~regional_identificacao,ncol=2,drop=FALSE)+
      eixo+ggplot2::labs(title=paste("Estimativas de",v,"por coorte semanal"),
          x="Coorte semanal",y=v)+tema
    salvar_plot(g,paste0("10F_",v,"_COORTES_SEMANAIS.pdf"),13,15)
  }

  jag <- gg(res_agregado,ggplot2::aes(x=semana,y=taxa_sub_mediana))+
    ggplot2::geom_ribbon(ggplot2::aes(ymin=taxa_sub_li95,
                                      ymax=taxa_sub_ls95),
                          fill="#B15767",alpha=.2,na.rm=TRUE)+
    ggplot2::geom_line(color="#AE4854",linewidth=.85,na.rm=TRUE)+
    ggplot2::geom_point(color="#AE4854",size=1.1,na.rm=TRUE)+
    eixo+ggplot2::coord_cartesian(ylim=c(0,100))+
    ggplot2::labs(title="Nao captacao estimada nas dez Regionais selecionadas",
       subtitle="Agregacao das coortes semanais (sem duplicacao de semanas)",
       x="Coorte semanal",y="Nao captacao (%)")+tema
  salvar_plot(jag,"11A_SUB_AGREGADA_TOP10_COORTES_SEMANAIS.pdf",12,6)

  semanal_pr <- aggregate(cbind(inicio_sintomas,notificacoes)~semana+inicio,
                         data=contagens,FUN=sum)
  gpr <- gg(semanal_pr,ggplot2::aes(x=inicio))+
    ggplot2::geom_line(ggplot2::aes(y=notificacoes,color="Notificacoes"),
                       linewidth=.8)+
    ggplot2::geom_line(ggplot2::aes(y=inicio_sintomas,
                                   color="Inicio dos sintomas"),linewidth=.8)+
    ggplot2::labs(title="Parana: dados observados de dengue em 2025",
       subtitle="Descritiva anual das 22 regionais; nao sao estimativas de N/U",
       x="Semana de calendario",y="Registros",color="Serie")+tema
  salvar_plot(gpr,"11B_DESCRITIVA_22_REGIONAIS_ANO_2025.pdf",12,6)

  gpainel <- gg(contagens_top10,ggplot2::aes(x=inicio,y=notificacoes))+
    ggplot2::geom_line(color="#255072",linewidth=.6)+
    ggplot2::facet_wrap(~regional_identificacao,ncol=2,scales="free_y")+
    ggplot2::labs(title="Contagens semanais observadas nas dez Regionais",
      subtitle="Descritiva do ano inteiro; nao confundir com N latente",
      x="Semana de calendario",y="Notificacoes")+tema
  salvar_plot(gpainel,"11C_DESCRITIVA_TOP10_ANO_2025.pdf",13,15)

  gord <- gg(anual_observado,ggplot2::aes(x=reorder(regional_nome,notificacoes),
                                            y=notificacoes))+
    ggplot2::geom_col(fill="#255072")+
    ggplot2::coord_flip()+
    ggplot2::labs(title="Notificacoes observadas no Parana em 2025",
      subtitle="22 Regionais, ordenadas por volume absoluto",
      x="Regional de Saude",y="Notificacoes em 2025")+tema
  salvar_plot(gord,"11D_RANKING_22_REGIONAIS_2025.pdf",11,9)
}

# ============================================================================
# 9. MAPAS DA ESTIMATIVA POR SEMANA (COORTES DISTINTAS)
# ============================================================================
if(config$salvar_mapas) {
  rotulos <- mapa_regional[
    mapa_regional$regional_num %in% ids_modelo, , drop=FALSE]
  rotulos <- sf::st_transform(rotulos,31982)
  rotulos <- sf::st_point_on_surface(rotulos)
  rotulos <- sf::st_transform(rotulos,4674)
  rotulos$codigo_rs <- sprintf("%02d",rotulos$regional_num)
  mapa_resultado <- function(dados,coluna,titulo,arquivo) {
    dados$regional_num <- dados$regional
    tmp <- dplyr::left_join(mapa_regional,dados,by="regional_num")
    g <- ggplot2::ggplot(tmp)+
      ggplot2::geom_sf(ggplot2::aes(fill=.data[[coluna]]),
                       color="white",linewidth=.25)+
      ggplot2::geom_sf_text(data=rotulos,
                            ggplot2::aes(label=codigo_rs),
                            size=3,fontface="bold",check_overlap=TRUE)+
      ggplot2::scale_fill_gradient(low="#EDF2F6",high="#255072",
                                   na.value="grey85")+
      ggplot2::coord_sf(datum=NA)+
      ggplot2::labs(title=titulo,fill=coluna,
        caption="Numeros identificam as dez RS; cinza = nao selecionada/sem ajuste")+
      ggplot2::theme_void(base_size=11)
    ggplot2::ggsave(file.path(config$pasta_saida,arquivo),
                    g,width=11,height=8,device="pdf")
  }
  map_obs <- anual_observado[
    anual_observado$regional %in% ids_modelo,,drop=FALSE]
  mapa_resultado(map_obs,"notificacoes",
    "Notificacoes observadas em 2025: dez maiores RS",
    "13_MAPA_DESCRITIVO_NOTIFICACOES_2025.pdf")
  semanas_mapas <- unique(c(config$semanas_para_mapas,
                    min(coortes$semana),max(coortes$semana)))
  semanas_mapas <- semanas_mapas[semanas_mapas %in% coortes$semana]
  for(t in semanas_mapas) {
    a <- res_semanais[res_semanais$semana==t,,drop=FALSE]
    mapa_resultado(a,"p_mediana",
      paste0("p: coorte de inicio dos sintomas na semana ",t),
      sprintf("14_MAPA_P_COORTE_SEMANA_%02d.pdf",t))
    mapa_resultado(a,"taxa_sub_mediana",
      paste0("Nao captacao: coorte de sintomas na semana ",t),
      sprintf("15_MAPA_SUB_COORTE_SEMANA_%02d.pdf",t))
  }
}

# ============================================================================
# 10. SENSIBILIDADE OPCIONAL A PRIORI DE LAMBDA (UMA COORTE)
# ============================================================================
if(config$rodar_sensibilidade) {
  r <- config$regional_sensibilidade
  t <- config$semana_sensibilidade
  if(!r %in% ids_modelo || !t %in% coortes$semana)
    stop("Coorte/regional escolhida para sensibilidade nao elegivel")
  semanas_not <- semanas_not_por_coorte[[paste(r,t,sep="_")]]
  y <- tabulate(if(is.null(semanas_not)) integer(0) else semanas_not,
                nbins=4L)
  if(sum(y)<config$min_notificacoes_ajuste) {
    warning("Coorte escolhida nao possui notificacoes suficientes")
  } else {
    sens <- lapply(seq_along(config$taxas_gamma_sensibilidade),function(i) {
      cfg <- config
      cfg$bl <- config$taxas_gamma_sensibilidade[i]
      a <- gibbs_coorte(y,cfg,
                  config$semente+900000L+100L*i)$amostras
      data.frame(regional=r,semana=t,
           taxa_priori_gamma=cfg$bl,
           media_priori_lambda=cfg$al/cfg$bl,
           as.list(c(resumir(a$p,"p"),resumir(a$N,"N"),
                     resumir(a$U,"U"),resumir(a$theta,"theta"))),
           check.names=FALSE)
    })
    salvar_csv(do.call(rbind,sens),"12_SENSIBILIDADE_COORTE_SEMANAL.csv")
  }
}

# ============================================================================
# 11. FINALIZACAO
# ============================================================================
cat("\n===================== CONCLUIDO =====================\n")
cat("Saida: ",normalizePath(config$pasta_saida),"\n",sep="")
cat("Coortes de sintomas: ",nrow(coortes),"\n",sep="")
cat("Ajustes previstos: ",nrow(res_semanais),"\n",sep="")
cat("Ajustes realizados: ",sum(res_semanais$status=="ajustado"),"\n",sep="")
cat("Sem ajuste: ",sum(res_semanais$status!="ajustado"),"\n",sep="")
cat("Agora: uma distribuicao multinomial por (regional, semana de sintomas).\n")
cat("As categorias pi1..pi4 representam SEMANAS DE ACOMPANHAMENTO, nao semanas civis.\n")
cat("Nao reutilize estimativas do modelo DIARIO ou janelas moveis (ex.: 12,6%%): REESTIMAR.\n")
cat("Totais N/U podem ser agregados apenas entre coortes nao sobrepostas.\n")
print(utils::head(res_semanais[,c("regional_identificacao","semana",
          "C_observado","p_mediana","N_mediana","U_mediana",
          "taxa_sub_mediana","status")],12),row.names=FALSE)

