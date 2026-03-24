# ——————————————————————————————————————————————————————————————————————————————
# MODELO A (CASCATA): SELEÇÃO HEURÍSTICA E ALIMENTAÇÃO VIA EXCEL ####
# ——————————————————————————————————————————————————————————————————————————————

pacman::p_load(tidyverse, ggplot2, readxl)

# 1. MOTOR CENTRAL DO ENSEMBLE (NÚCLEO DO MODELO A) ----------------------------
core_ensemble <- function(
    df_h, df_n, var_alvo="Alvo", lim_sde=4, alpha=0.8, teto=1, k_min=2, 
    max_na=2
) {
  y_real <- df_h[[var_alvo]]
  f_raw <- setdiff(names(df_h), c("Paredao", "Candidato", var_alvo))
  
  # Filtro rigoroso de tolerância a NAs no histórico
  qtd_na <- sapply(f_raw, function(x) sum(is.na(df_h[[x]])))
  fontes <- f_raw[qtd_na <= max_na]
  
  n <- length(y_real)
  
  df_h_c <- df_h
  df_n_c <- df_n
  f_cor <- list()
  
  m_b <- sapply(fontes, function(f) mean(abs(df_h[[f]] - y_real), na.rm=TRUE))
  
  for (f in fontes) {
    erros <- y_real - df_h[[f]]
    me <- mean(erros, na.rm=TRUE)
    sde <- sd(erros, na.rm=TRUE)
    
    if ((sde <= lim_sde) && (abs(me) >= 0.5)) {
      df_h_c[[f]] <- df_h[[f]] + me
      df_n_c[[f]] <- df_n[[f]] + me
      f_cor[[f]] <- me
    } else {
      f_cor[[f]] <- 0
    }
  }
  
  m_c <- sapply(fontes, function(f) {
    mean(abs(df_h_c[[f]] - y_real), na.rm=TRUE)
  })
  
  rank <- tibble(Fonte=fontes, MAE_B=m_b, Shift=unlist(f_cor), MAE_C=m_c) %>% 
    arrange(MAE_C)
  
  f_rq <- rank$Fonte
  res_m <- tibble()
  
  # TRATAMENTO DE NAs (IMPUTAÇÃO NEUTRA POR LINHA) -----------------------------
  imputar_na <- function(df) {
    for(i in 1:nrow(df)) {
      lin <- as.numeric(df[i, fontes])
      if(any(is.na(lin))) df[i, fontes][is.na(lin)] <- mean(lin, na.rm=TRUE)
    }
    return(df)
  }
  
  df_h_c <- imputar_na(df_h_c)
  df_n_c <- imputar_na(df_n_c)
  # ----------------------------------------------------------------------------
  
  for (k in 2:length(f_rq)) {
    t_k <- f_rq[1:k]
    pesos <- (1 / m_c[t_k])^alpha
    pesos <- pesos / sum(pesos)
    
    if (max(pesos) > teto) {
      id_m <- which.max(pesos)
      pesos[id_m] <- teto
      id_o <- setdiff(1:k, id_m)
      s_out <- sum(pesos[id_o])
      if (s_out > 0) pesos[id_o] <- pesos[id_o] * ((1 - teto)/s_out)
    }
    
    y_p <- as.numeric(as.matrix(df_h_c[t_k]) %*% pesos)
    mae_e <- mean(abs(y_real - y_p))
    
    res_m <- bind_rows(
      res_m, 
      tibble(
        k=k, 
        MAE=mae_e, 
        IC=n*log(mae_e)+k*log(n),
        Fontes=paste(t_k, collapse=" + "), 
        Pesos=list(pesos)
      )
    )
  }
  
  m_opt <- res_m %>% 
    filter(k >= k_min) %>% 
    filter(IC <= min(IC) + 2) %>% 
    slice(1)
  
  f_o <- unlist(strsplit(m_opt$Fontes, " \\+ "))
  p_o <- m_opt$Pesos[[1]]
  
  pred_h <- as.numeric(as.matrix(df_h_c[f_o]) %*% p_o)
  
  bktest <- df_h_c %>% 
    select(any_of(c("Paredao", "Candidato"))) %>% 
    mutate(Oficial=y_real, Pred=pred_h, Residuo=y_real-pred_h)
  
  pred_f <- as.numeric(as.matrix(df_n_c[f_o]) %*% p_o)
  
  list(Rank=rank, Grade=res_m, Mod=m_opt, F_opt=f_o, P_opt=p_o, 
       Pred=pred_f, Backtest=bktest)
}

# 2. FUNÇÃO CASCATA (WRAPPER HIERÁRQUICO COM DIAGNÓSTICO) ----------------------
imprimir_diag_a <- function(res, nivel) {
  cat(sprintf("\n# %s: DIAGNÓSTICO DETALHADO (MOD. A) --------------\n", nivel))
  cat(">>> ETAPA 1: RANKING E CALIBRAÇÃO\n")
  
  df_rank <- res$Rank %>% 
    mutate(across(is.numeric, ~round(.x, 3)))
  
  print(as.data.frame(df_rank))
  
  cat("\n>>> ETAPA 2/3: SELEÇÃO ÓTIMA (IC-MAE)\n")
  cat(sprintf("Modelo (k=%d) | IC=%.3f\n", res$Mod$k, res$Mod$IC))
  
  for(i in seq_along(res$F_opt)) {
    cat(sprintf(" - %s: %.1f%%\n", res$F_opt[i], res$P_opt[i]*100))
  }
  
  cat("\n>>> ETAPA 4: BACKTESTING IN-SAMPLE LOCAL\n")
  
  df_bk <- res$Backtest %>% 
    mutate(across(is.numeric, ~round(.x, 2)))
  
  print(as.data.frame(df_bk))
  cat(sprintf("RMSE Local: %.3f p.p.\n", sqrt(mean(res$Backtest$Residuo^2))))
}

ensemble_cascata_mod_a <- function(df_h, df_n, tol_na=6) {
  f <- setdiff(names(df_h), c("Paredao", "Candidato", "Posicao", "Oficial"))
  
  cands <- df_n %>% 
    arrange(Posicao) %>% 
    pull(Candidato)
  
  h1 <- df_h %>% 
    filter(Posicao==1) %>% 
    select(Paredao, Candidato, Alvo=Oficial, all_of(f))
  
  n1 <- df_n %>% 
    filter(Posicao==1) %>% 
    select(all_of(f))
  
  r1 <- core_ensemble(h1, n1, k_min = 2, max_na = tol_na)
  p1 <- r1$Pred
  
  r1$Backtest <- r1$Backtest %>% 
    mutate(Paredao = h1$Paredao)
  
  imprimir_diag_a(r1, "NÍVEL 1 (ALVO PRINCIPAL)")
  cat(sprintf("\n[!] PROJEÇÃO NÍVEL 1 (%s): %.2f%%\n", cands[1], p1))
  
  calc_s <- function(x) { x[1] / (x[1] + x[2]) * 100 }
  
  h2 <- df_h %>% 
    filter(Posicao %in% 2:3) %>% 
    arrange(Paredao, Posicao) %>% 
    group_by(Paredao) %>% 
    summarise(
      Candidato=Candidato[1], 
      Alvo=calc_s(Oficial),
      across(all_of(f), calc_s), 
      .groups="drop"
    )
  
  n2 <- df_n %>% 
    filter(Posicao %in% 2:3) %>% 
    arrange(Posicao) %>% 
    summarise(across(all_of(f), calc_s))
  
  r2 <- core_ensemble(h2, n2, k_min = 3, max_na = tol_na)
  s2 <- r2$Pred / 100 
  
  imprimir_diag_a(r2, "NÍVEL 2 (SHARE RELATIVO)")
  
  bk2_adj <- r2$Backtest %>% 
    mutate(Paredao=h2$Paredao, Cand_2=h2$Candidato) %>% 
    left_join(
      r1$Backtest %>% 
        select(Paredao, P_L1=Pred), 
      by="Paredao"
    ) %>% 
    left_join(
      df_h %>% 
        filter(Posicao==2) %>% 
        select(Paredao, Of_Abs=Oficial), 
      by="Paredao"
    ) %>% 
    mutate(P_Abs = (100 - P_L1) * (Pred / 100), Res_Abs = Of_Abs - P_Abs)
  
  cat("\n# NÍVEL 2: DIAGNÓSTICO GLOBAL RECOMPOSTO (MOD. A) ---------\n")
  
  df_pr_a <- bk2_adj %>% 
    select(Paredao, Cand_2, Of_Abs, P_Abs, Res_Abs) %>% 
    rename(Oficial=Of_Abs, Predito=P_Abs, Residuo=Res_Abs) %>% 
    mutate(across(is.numeric, ~round(.x, 2)))
  
  print(as.data.frame(df_pr_a))
  
  cat(sprintf("RMSE Global (2º Colocado): %.3f p.p.\n", 
              sqrt(mean(bk2_adj$Res_Abs^2))))
  
  cat(sprintf("\n[!] SHARE %s SOBRE %s: %.1f%%\n", toupper(cands[2]), 
              toupper(cands[3]), s2 * 100))
  
  p2 <- (100 - p1) * s2
  p3 <- (100 - p1) * (1 - s2)
  t_fin <- tibble(Posicao=1:3, Candidato=cands, Projecao=round(c(p1,p2,p3), 2))
  
  invisible(list(N1=r1, N2=r2, Final=t_fin, Backtest_N1=r1$Backtest))
}

# 3. VISUALIZAÇÃO DINÂMICA (GRÁFICO) -------------------------------------------
plotar_individual <- function(res, num_p, subtitulo = "") {
  df <- res$Final %>% 
    mutate(
      Candidato = factor(Candidato, levels = Candidato),
      Label = paste0(scales::number(Projecao, accuracy=0.01, dec=","), "%")
    )
  
  ggplot(df, aes(x=Candidato, y=Projecao, fill=Candidato)) +
    geom_bar(stat="identity", width=0.7) +
    geom_text(aes(label=Label), size=10, vjust=-0.5) +
    labs(
      title=sprintf("BBB 26 - %dº Paredão", num_p),
      subtitle=subtitulo,
      caption=sprintf("Última atualização: %s", format(Sys.time(), "%A, %Hh%M")), 
      x="", 
      y=""
    ) +
    theme_minimal() + 
    theme(
      legend.position="none", 
      plot.title=element_text(size=20, face="bold"),
      axis.text.x=element_text(size=20, face="bold"),
      axis.text.y=element_blank(), 
      panel.grid=element_blank()
    ) +
    scale_y_continuous(limits=c(0, max(df$Projecao)+15)) +
    scale_fill_brewer(palette="Pastel1")
}

# 4. LEITURA DO EXCEL E EXECUÇÃO ===============================================

# 1. Carrega os dados do Excel
df_completo <- read_excel("dados_paredoes.xlsx") %>% 
  mutate(Media = round((Insta + X + YT) / 3, 2))

# 2. Divisão automática: Histórico (tem dado Oficial) vs Novo (sem Oficial)
hist_cascata <- df_completo %>% 
  filter(!is.na(Oficial))

novo_cascata <- df_completo %>% 
  filter(is.na(Oficial))

# 3. Inferência e Plotagem (Ajuste tol_na conforme o rigor desejado)
res_mod_A <- ensemble_cascata_mod_a(hist_cascata, novo_cascata)

plotar_individual(res_mod_A, num_p = 10, subtitulo = "Estimativa (Modelo A)")

# Olhando as melhores e piores fontes
view(res_mod_A$N1$Rank)








