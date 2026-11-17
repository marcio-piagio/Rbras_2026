# 70ª RBras — Banner de dengue no Paraná (modelo Gibbs semanal)

## Arquivos principais

- `poster.tex`: código LaTeX do banner, no estilo original.
- `rbras-poster.cls`: classe original preservada, com cartões de conclusão 30% mais baixos.
- `poster.pdf`: banner A0 compilado com `pdflatex`.
- `GERAR_GRAFICOS_E_CONFERIR_BANNER.R`: script R para regenerar os gráficos a partir dos CSVs, com opção de reexecutar o Gibbs.
- `10D_heatmap_P_COORTES_SEMANAIS.pdf`: mapa de calor azul de p; foram retirados os 210 valores internos para melhorar a legibilidade à distância.
- `11A_CURVA_SEMANAL_SUB_TOP10_ICR95.pdf`: curva de não captação agregada por semana, com faixa de credibilidade de 95% calculada **sorteio a sorteio**, não pela soma dos intervalos de células.
- `10C_heatmap_SUB_COORTES_SEMANAIS.pdf`: mapa de calor azul de não captação, apenas complementar.
- `COORTES_210_MESMA_EXECUCAO.csv`, `REGIONAIS_10_MESMA_EXECUCAO.csv`, `SERIE_21_SEMANAS_POSTERIOR.csv`, `TOTAL_MESMA_EXECUCAO.csv`: resultados da **mesma reprodução independente MCMC**.
- `CONTAGENS_Y1_Y4_210_COORTES_FONTE_RELATORIO.csv`: contagens obtidas do relatório para a reprodução (não são uma nova auditoria do arquivo `.dbc`).
- `LaTeX/` e logos na raiz: extraídos do próprio PDF do banner enviado pelo usuário, para manter a apresentação.

## Como gerar os gráficos no RStudio

Descompacte o ZIP por inteiro em uma única pasta. Em seguida, execute:

```r
source("E:/banner_rbras/BANNER_RBRAS_CORRECOES_FINAIS_2025/GERAR_GRAFICOS_E_CONFERIR_BANNER.R", encoding="UTF-8")
```

Substitua o caminho pela pasta real. Com `config$recalcular_gibbs = FALSE` (padrão), o script usa os CSVs incluídos e **não executa o Gibbs**. Os gráficos são gravados ao lado do `.R`.

O comando LaTeX é `pdflatex poster.tex` (executar duas vezes) na pasta do banner. É necessário manter os arquivos de imagem e a classe `.cls` na estrutura do ZIP.

## Recalcular posterior no R (opcional)

Altere `recalcular_gibbs = TRUE`. O R reestimará as 210 combinações com 3 cadeias, 100.000 iterações por cadeia, burn-in de 20.000 e thin=50; usará o vetor de quatro contagens por coorte. Gravará **novos** CSVs e figuras em `REEXECUCAO_R` para não sobrescrever os números já incluídos no banner. Se optar por essa nova execução, **também atualize as duas tabelas e a mensagem principal em `poster.tex`**: não misture resultados de execuções diferentes.

A execução de Gibbs pode levar algum tempo no computador local. Pacote R necessário: `ggplot2`; o pacote `coda` é opcional caso deseje calcular ESS separadamente.

## Resultados desta reprodução independente

- C observado: **78.716**.
- N latente: **89.600** (ICr95%: **87.830–92.097**).
- U não captado: **10.884** (ICr95%: **9.114–13.381**).
- 100U/N: **12,1%** (ICr95%: **10,4%–14,5%**).
- 21 coortes semanais distintas, dez RS e 210 ajustes. Maior Rhat clássico por célula e parâmetro: aproximadamente **1,012**.

*Os valores representam estimativas condicionadas ao modelo e às distribuições a priori. A reprodução foi feita a partir das contagens semanais extraídas do relatório, **não** diretamente do arquivo SINAN `DENGBR25.dbc`. A análise de sensibilidade às prioris (especialmente a de lambda) e uma auditoria independente do `.dbc` continuam recomendadas antes da submissão definitiva.*

O gráfico de linha conecta estimativas de coortes modeladas **independentemente**: não é um modelo temporal suavizado.
