# ============================================================================
# 70a RBras - Dengue PR 2025 - GRAFICOS AZUIS E VERIFICACAO DO BANNER
# Modelo semanal puro: 10 Regionais x 21 coortes de inicio dos sintomas.
#
# PADRAO: LER os 4 CSVs da mesma reproducao MCMC, INCLUIDOS neste pacote.
# OPCIONAL: recalcular_gibbs=TRUE para nova amostragem em R baseada nas 210
# contagens Y1..Y4. A nova execucao e salva EM SUBPASTA, preservando os CSVs
# e os graficos originais. NAO revalida os microdados DENGBR25.dbc.
# ============================================================================
config <- list(
  pasta_arquivos = NA_character_,       # NA: pasta do script se usar source()
  recalcular_gibbs = FALSE,           # TRUE: executar novamente 210 Gibbs
  iter = 100000L, burn = 20000L, thin = 50L, cadeias = 3L,
  semente = 20261117L,
  a_p = 1, b_p = 1, a_theta = 1, b_theta = 1,
  a_lambda = 2, b_lambda = 0.02 / 7,
  azul = "#073A70",
  conferir_total = 78716L
)

if (!requireNamespace("ggplot2", quietly = TRUE))
  stop("Instale o pacote ggplot2: install.packages('ggplot2')")

pasta_do_script <- function() {
  ff <- lapply(sys.frames(), function(z) z$ofile)
  ff <- Filter(function(z) !is.null(z) && length(z) == 1L && nzchar(z), ff)
  if (length(ff)) return(dirname(normalizePath(tail(ff, 1L)[[1L]], mustWork=TRUE)))
  getwd()
}
pasta <- if (is.na(config$pasta_arquivos)) pasta_do_script() else config$pasta_arquivos
leia <- function(nome, em=pasta) {
  arq <- file.path(em,nome)
  if (!file.exists(arq)) stop("Arquivo ausente: ",arq,
    "\nExtraia o ZIP inteiro e execute o script com source('caminho/GERAR_GRAFICOS_E_CONFERIR_BANNER.R').")
  utils::read.csv(arq, check.names=FALSE, fileEncoding="UTF-8-BOM")
}
salve <- function(x,nome,em) utils::write.csv(x,file.path(em,nome),
                                         row.names=FALSE,fileEncoding="UTF-8",na="")
resuma <- function(z,prefixo) {
  q <- unname(stats::quantile(z,c(.025,.5,.975),na.rm=TRUE))
  setNames(q,paste0(prefixo,c("_li95","_mediana","_ls95")))
}

# ---------------------------------------------------------------------------
# MODO OPCIONAL - recalcular TODAS as posteriores a partir de Y1..Y4.
# Os indices r,t sao independentes no modelo, nao ha suavizacao temporal.
# Os intervalos dos agregados SAO calculados por sorteio, nao por somar quantis.
# ---------------------------------------------------------------------------
if (isTRUE(config$recalcular_gibbs)) {
  d <- leia("CONTAGENS_Y1_Y4_210_COORTES_FONTE_RELATORIO.csv")
  d <- d[order(d$rs,d$semana),,drop=FALSE]
  stopifnot(nrow(d)==210L,!anyDuplicated(d[,c("rs","semana")]),
    all(d$C==rowSums(d[,c("y1","y2","y3","y4")])),
    length(unique(d$rs))==10L,length(unique(d$semana))==21L)
  if (sum(d$C)!=config$conferir_total) stop("Total C nao confere")
  C <- as.integer(d$C)
  W <- as.integer(d$y2+2*d$y3+3*d$y4)
  S <- as.integer(d$y2+d$y3+d$y4)
  m <- nrow(d)
  stopifnot(config$iter>config$burn,config$cadeias>=2L,
            (config$iter-config$burn) %% config$thin == 0L)
  por_cadeia <- (config$iter-config$burn) %/% config$thin
  total <- por_cadeia * config$cadeias
  mats <- lapply(c("N","U","p","theta"),function(nm)
    matrix(NA_real_,nrow=total,ncol=m))
  names(mats) <- c("N","U","p","theta")
  for (ch in seq_len(config$cadeias)) {
    set.seed(config$semente+ch)
    p <- rep(.7,m); theta <- rep(.5,m); lambda <- pmax(C/.7,1)
    ix <- 0L
    cat("Executando cadeia ",ch,"/",config$cadeias,"...\n",sep="")
    for (it in seq_len(config$iter)) {
      A <- stats::rpois(m,lambda*(1-p)*(1-theta))
      B <- stats::rpois(m,lambda*theta*(1-p)^4)
      N <- C+A+B
      p <- stats::rbeta(m,config$a_p+C,config$b_p+W+A+4*B)
      theta <- stats::rbeta(m,config$a_theta+S+B,config$b_theta+A)
      lambda <- stats::rgamma(m,shape=config$a_lambda+N,
                                 rate=config$b_lambda+1)
      if (it>config$burn && ((it-config$burn) %% config$thin == 0L)) {
        ix <- ix+1L
        row <- (ch-1L)*por_cadeia+ix
        mats$N[row,] <- N
        mats$U[row,] <- A+B
        mats$p[row,] <- p
        mats$theta[row,] <- theta
      }
    }
  }
  stopifnot(all(sweep(mats$N-mats$U,2L,C,"==")))
  res <- d[,c("rs","nome","semana","C")]
  per_cell <- function(mat,nm) {
    q <- apply(mat,2L,stats::quantile,c(.025,.5,.975),na.rm=TRUE)
    res[[paste0(nm,"_li95")]] <<- unname(q[1L,])
    res[[paste0(nm,"_mediana")]] <<- unname(q[2L,])
    res[[paste0(nm,"_ls95")]] <<- unname(q[3L,])
  }
  per_cell(mats$p,"p"); per_cell(mats$theta,"theta")
  per_cell(mats$N,"N"); per_cell(mats$U,"U")
  per_cell(100*mats$U/mats$N,"taxa_sub")
  # Rhat classico por celula e parametro. ESS pode ser calculado com coda;
  # Rhat sozinho nao garante a qualidade da inferencia.
  rhat <- function(z) {
    partes <- lapply(seq_len(config$cadeias),function(j)
      z[((j-1L)*por_cadeia+1L):(j*por_cadeia),,drop=FALSE])
    medias <- vapply(partes,colMeans,numeric(m))
    Wchain <- rowMeans(vapply(partes,function(v)
      apply(v,2L,stats::var),numeric(m)))
    Bchain <- por_cadeia*apply(medias,1L,stats::var)
    ifelse(Wchain>0,
      sqrt(((por_cadeia-1)/por_cadeia*Wchain+Bchain/por_cadeia)/Wchain),NA_real_)
  }
  for(nm in names(mats)) res[[paste0("rhat_",nm)]] <- rhat(mats[[nm]])
  ordem <- unique(c(17L,14L,15L,19L,10L,16L,12L,18L,2L,13L))
  rr <- vector("list",length(ordem))
  for (i in seq_along(ordem)) {
    id <- ordem[i]; idx <- which(d$rs==id)
    NN <- rowSums(mats$N[,idx,drop=FALSE]); UU <- rowSums(mats$U[,idx,drop=FALSE])
    pp <- rowSums(mats$N[,idx,drop=FALSE]*mats$p[,idx,drop=FALSE])/NN
    tt <- rowSums(mats$N[,idx,drop=FALSE]*mats$theta[,idx,drop=FALSE])/NN
    rr[[i]] <- data.frame(rs=id,nome=d$nome[idx[1]],C=sum(C[idx]),
        as.list(c(resuma(pp,"p"),resuma(tt,"theta"),resuma(NN,"N"),
                  resuma(UU,"U"),resuma(100*UU/NN,"taxa_sub"))),
        check.names=FALSE)
  }
  reg <- do.call(rbind,rr)
  ss <- vector("list",21L)
  for (wk in 1:21) {
    idx <- which(d$semana==wk)
    NN <- rowSums(mats$N[,idx,drop=FALSE]); UU <- rowSums(mats$U[,idx,drop=FALSE])
    ss[[wk]] <- data.frame(semana=wk,C=sum(C[idx]),n_regionais=length(idx),
         as.list(c(resuma(NN,"N"),resuma(UU,"U"),
                   resuma(100*UU/NN,"taxa_sub"))),check.names=FALSE)
  }
  serie <- do.call(rbind,ss)
  NN <- rowSums(mats$N); UU <- rowSums(mats$U)
  glob <- data.frame(C=sum(C),coortes=21L,regionais=10L,ajustes=210L,
     as.list(c(resuma(NN,"N"),resuma(UU,"U"),
        resuma(100*UU/NN,"taxa_sub"))),check.names=FALSE)
  pasta_nova <- file.path(pasta,"REEXECUCAO_R")
  dir.create(pasta_nova,showWarnings=FALSE,recursive=TRUE)
  salve(res,"COORTES_210_MESMA_EXECUCAO.csv",pasta_nova)
  salve(reg,"REGIONAIS_10_MESMA_EXECUCAO.csv",pasta_nova)
  salve(serie,"SERIE_21_SEMANAS_POSTERIOR.csv",pasta_nova)
  salve(glob,"TOTAL_MESMA_EXECUCAO.csv",pasta_nova)
  saveRDS(list(N=mats$N,U=mats$U,p=mats$p,theta=mats$theta,
    contagens=d,config=config),file.path(pasta_nova,"AMOSTRAS_AUDITORIA.rds"),compress="gzip")
  pasta_saida <- pasta_nova
  cat("Nova execucao: nao misture seus numeros com os do banner anterior.\n")
} else {
  pasta_saida <- pasta
  res <- leia("COORTES_210_MESMA_EXECUCAO.csv")
  reg <- leia("REGIONAIS_10_MESMA_EXECUCAO.csv")
  serie <- leia("SERIE_21_SEMANAS_POSTERIOR.csv")
  glob <- leia("TOTAL_MESMA_EXECUCAO.csv")
}

# ---------------------------------------------------------------------------
# AUDITORIA DA ESTRUTURA. Todos os CSVs devem ser da MESMA execucao.
# ---------------------------------------------------------------------------
stopifnot(nrow(res)==210L,nrow(reg)==10L,nrow(serie)==21L,nrow(glob)==1L,
          !anyDuplicated(res[,c("rs","semana")]),
          length(unique(res$rs))==10L,length(unique(res$semana))==21L,
          sum(res$C)==config$conferir_total,
          sum(reg$C)==config$conferir_total,
          all(is.finite(res$p_mediana)), all(res$p_mediana>=0 & res$p_mediana<=1),
          all(res$taxa_sub_mediana>=0 & res$taxa_sub_mediana<=100),
          all(serie$taxa_sub_li95<=serie$taxa_sub_mediana),
          all(serie$taxa_sub_mediana<=serie$taxa_sub_ls95))
if (glob$C[1]!=config$conferir_total) stop("C global nao confere")
cat(sprintf("Conferencia: %d coortes regionais, C=%d, N=%.0f, U=%.0f, nao captacao=%.1f%%\n",
            nrow(res),glob$C[1],glob$N_mediana[1],glob$U_mediana[1],
            glob$taxa_sub_mediana[1]))
if ("rhat_p" %in% names(res)) {
  rh <- max(unlist(res[grep("^rhat_",names(res))]),na.rm=TRUE)
  cat("Maior Rhat classico por coorte: ",round(rh,4),"\n",sep="")
  if (rh>1.05) warning("Rhat >1.05: avalie convergencia antes de apresentar")
}

# ---------------------------------------------------------------------------
# FIGURA 1: p por RS e semana. Sem numeros nas 210 celulas (legibilidade A0).
# FIGURA 2: 100U/N agregado POR SORTEIO e coorte, com ICr95% valido.
# ---------------------------------------------------------------------------
reg$regional_id <- sprintf("%02dª RS - %s",as.integer(reg$rs),reg$nome)
res$regional_id <- sprintf("%02dª RS - %s",as.integer(res$rs),res$nome)
res$regional_id <- factor(res$regional_id,levels=rev(reg$regional_id))
if (anyNA(res$regional_id)) stop("Regionais nao correspondem entre CSVs")

azuis <- c("#EDF3FA","#B4CADF","#5A86AD",config$azul)
base_graf <- ggplot2::theme_minimal(base_size=12) +
  ggplot2::theme(plot.title=ggplot2::element_text(face="bold",color=config$azul,size=15),
   plot.subtitle=ggplot2::element_text(color="#52606D",size=10),
   axis.text.y=ggplot2::element_text(size=9),
   axis.text.x=ggplot2::element_text(size=9),
   legend.position="right",plot.margin=ggplot2::margin(6,8,6,5))

mapa <- ggplot2::ggplot(res,ggplot2::aes(x=semana,y=regional_id,fill=p_mediana)) +
 ggplot2::geom_tile(color="white",linewidth=.40)+
 ggplot2::scale_fill_gradientn(colors=azuis,limits=c(0,1),name="p")+
 ggplot2::scale_x_continuous(breaks=seq_len(21),expand=ggplot2::expansion(add=.5))+
 ggplot2::labs(title="Probabilidade de primeira notificacao por Regional e semana",
     subtitle="Medianas posteriores de p; coortes independentes, sem suavizacao temporal",
     x="Semana de inicio dos sintomas (coorte)",y=NULL)+
 base_graf + ggplot2::theme(panel.grid=ggplot2::element_blank())

serie <- serie[order(serie$semana),,drop=FALSE]
curva <- ggplot2::ggplot(serie,ggplot2::aes(x=semana,y=taxa_sub_mediana)) +
 ggplot2::geom_ribbon(ggplot2::aes(ymin=taxa_sub_li95,ymax=taxa_sub_ls95),
                      fill="#8FB1CF",alpha=.40) +
 ggplot2::geom_line(color=config$azul,linewidth=1.05) +
 ggplot2::geom_point(color=config$azul,size=1.8) +
 ggplot2::scale_x_continuous(breaks=seq_len(21))+
 ggplot2::scale_y_continuous(limits=c(0,max(65,ceiling(max(serie$taxa_sub_ls95)/10)*10)),
                            expand=ggplot2::expansion(mult=c(0,.02)))+
 ggplot2::labs(title="Nao captacao semanal nas dez Regionais de Saude",
     subtitle="Mediana e ICr95% de 100 U/N, agregados por sorteio posterior",
     x="Semana de inicio dos sintomas (coorte)",y="Nao captacao estimada (%)")+
 base_graf + ggplot2::theme(panel.grid.minor=ggplot2::element_blank())

# Mapas de calor de U/N ficam como material complementar, nao no poster.
mapa_sub <- ggplot2::ggplot(res,ggplot2::aes(x=semana,y=regional_id,fill=taxa_sub_mediana)) +
 ggplot2::geom_tile(color="white",linewidth=.4) +
 ggplot2::scale_fill_gradientn(colors=azuis,name="U/N (%)") +
 ggplot2::scale_x_continuous(breaks=seq_len(21),expand=ggplot2::expansion(add=.5))+
 ggplot2::labs(title="Nao captacao por Regional e coorte semanal",
             subtitle="Mapa complementar das medianas posteriores",
             x="Semana de inicio dos sintomas (coorte)",y=NULL)+
 base_graf + ggplot2::theme(panel.grid=ggplot2::element_blank())

salvar_fig <- function(fig,arquivo,w=15,h=6.05) {
  caminho <- file.path(pasta_saida,arquivo)
  ggplot2::ggsave(caminho,plot=fig,width=w,height=h,units="in",device="pdf",bg="white")
  ggplot2::ggsave(sub("\\.pdf$",".png",caminho),plot=fig,width=w,height=h,
                  units="in",dpi=180,bg="white")
  cat("Gravado: ",caminho,"\n",sep="")
}
salvar_fig(mapa,"10D_heatmap_P_COORTES_SEMANAIS.pdf")
salvar_fig(curva,"11A_CURVA_SEMANAL_SUB_TOP10_ICR95.pdf")
salvar_fig(mapa_sub,"10C_heatmap_SUB_COORTES_SEMANAIS.pdf")
cat("\nFinalizado. Se recalcular_gibbs=TRUE, ATUALIZE as tabelas LaTeX\n",
    "com os CSVs da subpasta REEXECUCAO_R antes de submeter o banner.\n")
