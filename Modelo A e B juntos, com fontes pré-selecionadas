# ——————————————————————————————————————————————————————————————————————————————
# SISTEMA UNIFICADO DE PREVISÃO BBB (ENSEMBLE DE ENSEMBLES) 
# ——————————————————————————————————————————————————————————————————————————————

# Carregando
pacman::p_load(tidyverse, glmnet, ggplot2)


# 1. FUNÇÕES BASE E TRANSFORMADORES ============================================

logit <- function(p) {
  p <- pmax(pmin(p, 0.999), 0.001)
  log(p / (1 - p))
}

inv_logit <- function(x) exp(x) / (1 + exp(x))

# 2. MOTOR DO MODELO A (SELEÇÃO HEURÍSTICA VIA IC-MAE) =========================

core_ensemble <- function(
    df_h, df_n, var_alvo="Alvo", lim_sde = 4, alpha = 0.9, teto = 1, k_min = 2
) {
  y_real <- df_h[[var_alvo]]
  fontes <- setdiff(names(df_h), c("Paredao", "Candidato", var_alvo))
  n <- length(y_real)
  
  df_h_c <- df_h; df_n_c <- df_n; f_cor <- list()
  
  m_b <- sapply(fontes, function(f) mean(abs(df_h[[f]] - y_real), na.rm=TRUE))
  
  for (f in fontes) {
    erros <- y_real - df_h[[f]]
    me <- mean(erros, na.rm=TRUE); sde <- sd(erros, na.rm=TRUE)
    if ((sde <= lim_sde) && (abs(me) >= 0.5)) {
      df_h_c[[f]] <- df_h[[f]] + me; df_n_c[[f]] <- df_n[[f]] + me
      f_cor[[f]] <- me
    } else f_cor[[f]] <- 0
  }
  
  m_c <- sapply(fontes, function(f) mean(abs(df_h_c[[f]] - y_real), na.rm=TRUE))
  
  rank <- tibble(Fonte=fontes, MAE_B=m_b, Shift=unlist(f_cor), MAE_C=m_c) %>% 
    arrange(MAE_C)
  
  f_rq <- rank$Fonte; res_m <- tibble()
  
  for (k in 2:length(f_rq)) {
    t_k <- f_rq[1:k]
    pesos <- (1 / m_c[t_k])^alpha; pesos <- pesos / sum(pesos)
    if (max(pesos) > teto) {
      id_m <- which.max(pesos); pesos[id_m] <- teto
      id_o <- setdiff(1:k, id_m); s_out <- sum(pesos[id_o])
      if (s_out > 0) pesos[id_o] <- pesos[id_o] * ((1 - teto)/s_out)
    }
    y_p <- as.numeric(as.matrix(df_h_c[t_k]) %*% pesos)
    mae_e <- mean(abs(y_real - y_p))
    res_m <- bind_rows(res_m, tibble(k=k, MAE=mae_e, IC=n*log(mae_e)+k*log(n),
                                     Fontes=paste(t_k, collapse=" + "), 
                                     Pesos=list(pesos)))
  }
  
  m_opt <- res_m %>% filter(k >= k_min) %>% filter(IC <= min(IC) + 2) %>% 
    slice(1)
  
  f_opt <- unlist(strsplit(m_opt$Fontes, " \\+ ")); p_opt <- m_opt$Pesos[[1]]
  
  pred_h <- as.numeric(as.matrix(df_h_c[f_opt]) %*% p_opt)
  bktest <- df_h_c %>% select(any_of(c("Paredao", "Candidato"))) %>%
    mutate(Oficial=y_real, Pred=pred_h, Residuo=y_real-pred_h)
  pred_f <- as.numeric(as.matrix(df_n_c[f_opt]) %*% p_opt)
  
  list(Rank=rank, Grade=res_m, Mod=m_opt, F_opt=f_opt, P_opt=p_opt, 
       Pred=pred_f, Backtest=bktest)
}

imprimir_diag_a <- function(res, nivel) {
  cat(sprintf("\n# %s: DIAGNÓSTICO DETALHADO (MOD. A) --------------\n", nivel))
  cat(">>> ETAPA 1: RANKING E CALIBRAÇÃO\n")
  print(as.data.frame(res$Rank %>% 
                        mutate(across(where(is.numeric), ~ round(.x, 3)))))
  cat("\n>>> ETAPA 2/3: SELEÇÃO ÓTIMA (IC-MAE)\n")
  cat(sprintf("Modelo (k=%d) | IC=%.3f\n", res$Mod$k, res$Mod$IC))
  for(i in seq_along(res$F_opt)) 
    cat(sprintf(" - %s: %.1f%%\n", res$F_opt[i], res$P_opt[i]*100))
  cat("\n>>> ETAPA 4: BACKTESTING IN-SAMPLE LOCAL\n")
  print(as.data.frame(res$Backtest %>% mutate(across(is.numeric,~round(.x,2)))))
  cat(sprintf("RMSE Local: %.3f p.p.\n", sqrt(mean(res$Backtest$Residuo^2))))
}

ensemble_cascata_mod_a <- function(df_h, df_n) {
  fontes <- setdiff(names(df_h), c("Paredao", "Candidato", "Posicao", "Oficial"))
  df_n <- df_n %>% arrange(Posicao); cands <- df_n$Candidato
  
  # NÍVEL 1
  h1 <- df_h %>% filter(Posicao==1) %>% select(Paredao,Candidato,Alvo=Oficial,
                                               all_of(fontes))
  n1 <- df_n %>% filter(Posicao==1) %>% select(all_of(fontes))
  r1 <- core_ensemble(h1, n1, k_min = 2)
  p1 <- r1$Pred
  r1$Backtest <- r1$Backtest %>% mutate(Paredao = h1$Paredao)
  
  imprimir_diag_a(r1, "NÍVEL 1 (ALVO PRINCIPAL)")
  cat(sprintf("\n[!] PROJEÇÃO NÍVEL 1 (%s): %.2f%%\n", cands[1], p1))
  
  # NÍVEL 2
  calc_s <- function(x) { x[1] / (x[1] + x[2]) * 100 }
  h2 <- df_h %>% filter(Posicao %in% 2:3) %>% arrange(Paredao, Posicao) %>% 
    group_by(Paredao) %>% summarise(Candidato=Candidato[1], 
                                    Alvo=calc_s(Oficial),
                                    across(all_of(fontes), calc_s), .groups="drop")
  n2 <- df_n %>% filter(Posicao %in% 2:3) %>% arrange(Posicao) %>% 
    summarise(across(all_of(fontes), calc_s))
  r2 <- core_ensemble(h2, n2, k_min = 3)
  s2 <- r2$Pred / 100 
  
  imprimir_diag_a(r2, "NÍVEL 2 (SHARE RELATIVO)")
  
  # Diagnóstico Global Recomposto
  bk2_adj <- r2$Backtest %>% mutate(Paredao = h2$Paredao, Cand_2 = h2$Candidato) %>%
    left_join(r1$Backtest %>% select(Paredao, P_L1=Pred), by="Paredao") %>%
    left_join(df_h %>% filter(Posicao==2) %>% select(Paredao, Of_Abs=Oficial), 
              by="Paredao") %>%
    mutate(P_Abs = (100 - P_L1) * (Pred / 100), Res_Abs = Of_Abs - P_Abs)
  
  cat("\n# NÍVEL 2: DIAGNÓSTICO GLOBAL RECOMPOSTO (MOD. A) ---------\n")
  df_pr_a <- bk2_adj %>% select(Paredao, Cand_2, Of_Abs, P_Abs, Res_Abs) %>%
    rename(Oficial=Of_Abs, Predito=P_Abs, Residuo=Res_Abs) %>%
    mutate(across(is.numeric, ~round(.x, 2)))
  print(as.data.frame(df_pr_a))
  cat(sprintf("RMSE Global (2º Colocado): %.3f p.p.\n", 
              sqrt(mean(bk2_adj$Res_Abs^2))))
  
  cat(sprintf("\n[!] SHARE %s SOBRE %s: %.1f%%\n", toupper(cands[2]), 
              toupper(cands[3]), s2 * 100))
  
  p2 <- (100 - p1) * s2; p3 <- (100 - p1) * (1 - s2)
  t_fin <- tibble(Posicao=1:3, Candidato=cands, Projecao=round(c(p1,p2,p3), 2))
  
  invisible(list(N1=r1, N2=r2, Final=t_fin, Backtest_N1=r1$Backtest))
}

# 3. MOTOR DO MODELO B (RIDGE ESTRUTURAL CASCATA) ==============================

ensemble_cascata_mod_b <- function(df_h, df_n) {
  df_n <- df_n %>% arrange(Posicao); cands <- df_n$Candidato
  vars <- c("L_I", "L_X", "L_Y", "Pol")
  
  # NÍVEL 1
  prep_n1 <- function(df) {
    df %>% filter(Posicao==1) %>% mutate(
      I_p=Insta/100, X_p=X/100, Y_p=YT/100, Pol=abs(I_p-X_p)*Y_p,
      L_I=logit(I_p), L_X=logit(X_p), L_Y=logit(Y_p))
  }
  h1 <- prep_n1(df_h) %>% mutate(L_A=logit(Oficial/100)); n1 <- prep_n1(df_n)
  
  set.seed(42)
  cv1 <- cv.glmnet(x=as.matrix(h1[vars]), y=h1$L_A, alpha=0, 
                   nfolds=nrow(h1), standardize=TRUE)
  p1 <- inv_logit(as.numeric(predict(cv1, newx=as.matrix(n1[vars]), 
                                     s="lambda.min"))) * 100
  h1$Predito <- inv_logit(as.numeric(predict(cv1, newx=as.matrix(h1[vars]), 
                                             s="lambda.min"))) * 100
  h1$Residuo <- h1$Oficial - h1$Predito
  
  cat("\n# NÍVEL 1: DIAGNÓSTICO DETALHADO (MOD. B) --------------------\n")
  print(as.data.frame(h1 %>% select(Paredao, Candidato, Oficial, Predito, 
                                    Residuo) %>% 
                        mutate(across(is.numeric, ~round(.x, 2)))))
  cat(sprintf("RMSE: %.3f p.p.\n", sqrt(mean(h1$Residuo^2))))
  
  # NÍVEL 2
  prep_n2 <- function(df) {
    if (!"Paredao" %in% names(df)) df$Paredao <- "Novo"
    df %>% filter(Posicao %in% 2:3) %>% arrange(Paredao, Posicao) %>% 
      group_by(Paredao) %>% summarise(
        Cand_2=Candidato[1],
        A_S=if("Oficial" %in% names(.)) Oficial[1]/sum(Oficial) else NA,
        I_S=Insta[1]/sum(Insta), X_S=X[1]/sum(X), Y_S=YT[1]/sum(YT), 
        .groups="drop"
      ) %>% mutate(Pol=abs(I_S-X_S)*Y_S, L_I=logit(I_S), L_X=logit(X_S), 
                   L_Y=logit(Y_S))
  }
  h2 <- prep_n2(df_h) %>% mutate(L_A=logit(A_S)); n2 <- prep_n2(df_n)
  
  set.seed(42)
  cv2 <- cv.glmnet(x=as.matrix(h2[vars]), y=h2$L_A, alpha=0, 
                   nfolds=nrow(h2), standardize=TRUE)
  
  s2 <- inv_logit(as.numeric(predict(cv2, newx=as.matrix(n2[vars]), 
                                     s="lambda.min")))
  
  h2$Predito <- inv_logit(as.numeric(predict(cv2, newx=as.matrix(h2[vars]), 
                                             s="lambda.min")))
  h2$Residuo <- (h2$A_S - h2$Predito) * 100
  h2$Predito <- h2$Predito * 100
  h2$Oficial <- h2$A_S * 100
  
  cat("\n# NÍVEL 2: DIAGNÓSTICO LOCAL SHARE (MOD. B) ----------------\n")
  df_pb <- h2 %>% select(Paredao, Cand_2, Oficial, Predito, Residuo) %>% 
    mutate(across(is.numeric, ~round(.x, 2)))
  print(as.data.frame(df_pb))
  cat(sprintf("RMSE Local: %.3f p.p.\n", sqrt(mean(h2$Residuo^2))))
  
  # Diagnóstico Global Recomposto
  bk2_adj <- h2 %>% select(Paredao, Cand_2, Predito) %>%
    left_join(h1 %>% select(Paredao, P_L1=Predito), by="Paredao") %>%
    left_join(df_h %>% filter(Posicao==2) %>% select(Paredao, Of_Abs=Oficial), 
              by="Paredao") %>%
    mutate(P_Abs = (100 - P_L1) * (Predito / 100), Res_Abs = Of_Abs - P_Abs)
  
  cat("\n# NÍVEL 2: DIAGNÓSTICO GLOBAL RECOMPOSTO (MOD. B) ---------\n")
  df_pr_b <- bk2_adj %>% select(Paredao, Cand_2, Of_Abs, P_Abs, Res_Abs) %>%
    rename(Oficial=Of_Abs, Predito=P_Abs, Residuo=Res_Abs) %>%
    mutate(across(is.numeric, ~round(.x, 2)))
  print(as.data.frame(df_pr_b))
  cat(sprintf("RMSE Global (2º Colocado): %.3f p.p.\n", 
              sqrt(mean(bk2_adj$Res_Abs^2))))
  
  p2 <- (100 - p1) * s2; p3 <- (100 - p1) * (1 - s2)
  t_fin <- tibble(Posicao=1:3, Candidato=cands, Projecao=round(c(p1,p2,p3), 2))
  
  invisible(list(N1=cv1, N2=cv2, Final=t_fin, Backtest_N1=h1))
}
# 4. CONSOLIDAÇÃO E VISUALIZAÇÃO (O CONSENSO) ==================================

#' Gera a média estrutural e plota o gráfico final unificado
# 4. CONSOLIDAÇÃO E GRÁFICO (CONSENSO) -----------------------------------------
gerar_consenso <- function(res_a, res_b, num_p) {
  df <- res_a$Final %>% rename(P_A=Projecao) %>% 
    left_join(res_b$Final %>% rename(P_B=Projecao), by="Candidato") %>%
    mutate(P_Fin=(P_A+P_B)/2, Candidato=factor(Candidato, levels=Candidato),
           Label=paste0(scales::number(P_Fin, accuracy=0.01, dec=","), "%"))
  
  cat("\n# PROJEÇÃO CONSOLIDADA ------------------------------------------\n")
  print(as.data.frame(df %>% select(Candidato, P_A, P_B, P_Fin)))
  
  ggplot(df, aes(x=Candidato, y=P_Fin, fill=Candidato)) +
    geom_bar(stat="identity", width=0.7) +
    geom_text(aes(label=Label), size=9, vjust=-0.5) +
    labs(title=sprintf("BBB 26 - %dº Paredão", num_p),
         subtitle="Projeção Consolidada de Eliminação",
         caption=sprintf("Estimativa baseada em dois modelos independentes.\nAtualização: %s", 
                         format(Sys.time(), "%A, %Hh%M")), x="", y="") +
    theme_minimal() + theme(legend.position="none", 
                            plot.title=element_text(size=20, face="bold"),
                            plot.subtitle=element_text(size=16),
                            axis.text.x=element_text(size=20),
                            axis.text.y=element_blank(), panel.grid=element_blank()) +
    scale_y_continuous(limits=c(0, max(df$P_Fin)+15)) +
    scale_fill_brewer(palette="Pastel1")
}

# 5. PLOTAR INDIVIDUAL ========================================================
plotar_individual <- function(res, titulo, num_p, subtitulo = "") {
  df <- res$Final %>%
    mutate(Candidato = factor(Candidato, levels = Candidato),
           Label = paste0(scales::number(Projecao, accuracy = 0.01, dec = ","), "%"))
  
  ggplot(df, aes(x=Candidato, y=Projecao, fill=Candidato)) +
    geom_bar(stat="identity", width=0.7) +
    geom_text(aes(label=Label), size=9, vjust=-0.5) +
    labs(title=sprintf("BBB 26 - %dº Paredão", num_p),
         subtitle=subtitulo, x="", y="",
         caption=sprintf("Última atualização: %s", 
                         format(Sys.time(), "%A, %Hh%M")), x="", y="") +
    theme_minimal() + 
    theme(legend.position="none", 
          plot.title=element_text(size=22, face="bold"),
          plot.subtitle=element_text(size=16),
          axis.text.x=element_text(size=20),
          axis.text.y=element_blank(), panel.grid=element_blank()) +
    scale_y_continuous(limits=c(0, max(df$Projecao)+15)) +
    scale_fill_brewer(palette="Pastel1")
}

# 6. ALIMENTAÇÃO DE DADOS (BASE ÚNICA) ========================================

hist_cascata <- bind_rows(
  tibble(Paredao="P2", Candidato="Matheus", Posicao=1, Oficial=79.48,
         Insta=77.81, X=80.74, YT=66.42),
  tibble(Paredao="P2", Candidato="Leandro", Posicao=2, Oficial=15.55,
         Insta=12.63, X=14.12, YT=20.54),
  tibble(Paredao="P2", Candidato="Brigido", Posicao=3, Oficial=4.97,
         Insta=9.56,  X=5.14,  YT=13.04),
  
  tibble(Paredao="P3", Candidato="Brigido", Posicao=1, Oficial=77.88,
         Insta=73.20, X=74.55, YT=58.40),
  tibble(Paredao="P3", Candidato="Leandro", Posicao=2, Oficial=12.04,
         Insta=16.40, X=21.19, YT=23.60),
  tibble(Paredao="P3", Candidato="A_Paula", Posicao=3, Oficial=10.08,
         Insta=10.40, X=4.27,  YT=18.00),
  
  tibble(Paredao="P4", Candidato="Sarah",   Posicao=1, Oficial=69.13,
         Insta=69.33, X=66.04, YT=50.24),
  tibble(Paredao="P4", Candidato="Babu",    Posicao=2, Oficial=28.49,
         Insta=23.07, X=25.44, YT=42.92),
  tibble(Paredao="P4", Candidato="Sol",     Posicao=3, Oficial=2.38,
         Insta=7.60,  X=8.52,  YT=6.84),
  
  tibble(Paredao="P5", Candidato="Marcelo", Posicao=1, Oficial=68.56,
         Insta=65.25, X=78.25, YT=66.00),
  tibble(Paredao="P5", Candidato="Samira",  Posicao=2, Oficial=16.25,
         Insta=17.08, X=10.27, YT=20.09),
  tibble(Paredao="P5", Candidato="Solange", Posicao=3, Oficial=15.39,
         Insta=17.67, X=11.47, YT=13.95),
  
  tibble(Paredao="P6", Candidato="Maxiane", Posicao=1, Oficial=63.21,
         Insta=65.07, X=68.20, YT=52.17),
  tibble(Paredao="P6", Candidato="Milena",  Posicao=2, Oficial=36.11,
         Insta=33.27, X=26.93, YT=42.55),
  tibble(Paredao="P6", Candidato="Chaiany", Posicao=3, Oficial=0.68,
         Insta=1.67,  X=4.86,  YT=5.33),
  
  tibble(Paredao="P7", Candidato="Breno",   Posicao=1, Oficial=54.66,
         Insta=53.43, X=48.99, YT=43.55),
  tibble(Paredao="P7", Candidato="Cowboy",  Posicao=2, Oficial=43.12,
         Insta=41.29, X=42.72, YT=53.14),
  tibble(Paredao="P7", Candidato="Jordana", Posicao=3, Oficial=2.22,
         Insta=5.29,  X=8.29,  YT=3.25),
  
  tibble(Paredao="P8", Candidato="Babu",    Posicao=1, Oficial=68.62,
         Insta=70.33, X=77.77, YT=57.63),
  tibble(Paredao="P8", Candidato="Milena",  Posicao=2, Oficial=30.91,
         Insta=27.93, X=13.49, YT=38.24),
  tibble(Paredao="P8", Candidato="Chaiany", Posicao=3, Oficial=0.47,
         Insta=1.73,  X=8.75,  YT=4.58),
  
  tibble(Paredao="P9", Candidato="Breno",   Posicao=1, Oficial=58.96,
         Insta=55.38, X=63.33, YT=48.16),
  tibble(Paredao="P9", Candidato="A_Paula", Posicao=2, Oficial=25.17,
         Insta=16.13, X=13.84, YT=31.32),
  tibble(Paredao="P9", Candidato="Leandro", Posicao=3, Oficial=15.87,
         Insta=28.50, X=22.82, YT=20.18)
) %>%
  mutate(Media = round((Insta + X + YT) / 3, 2))

novo_cascata <- tibble(
  Candidato = c("Jonas", "Juliano", "Gabriela"),
  Posicao   = c(1, 2, 3), 
  Insta     = c(49.94, 35.0, 15.06),
  X         = c(55.70, 38.73,  5.56),
  YT        = c(38.08, 50.36, 12.12)
) %>%
  mutate(Media = round((Insta + X + YT) / 3, 2))


# 7. EXECUÇÃO UNIFICADA ======================================================

# Executa as inferências em paralelo usando a mesma base
res_mod_A <- ensemble_cascata_mod_a(hist_cascata, novo_cascata)
res_mod_B <- ensemble_cascata_mod_b(hist_cascata, novo_cascata)

# Consolida, imprime a tabela do consenso e gera o gráfico final
gerar_consenso(res_mod_A, res_mod_B, num_p = 10)

# 1. Gráfico do Modelo A
plotar_individual(res_mod_A, num_p = 10, subtitulo = "Estimativa dos votos para eliminação")

# 2. Gráfico do Modelo B
plotar_individual(res_mod_B, num_p = 10, subtitulo = "Estimativa dos votos para eliminação")

