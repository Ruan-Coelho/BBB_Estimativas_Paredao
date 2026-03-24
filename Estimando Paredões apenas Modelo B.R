# ——————————————————————————————————————————————————————————————————————————————
# MODELO B (CASCATA): RIDGE ESTRUTURAL COM LOGIT CONDICIONAL ####
# ——————————————————————————————————————————————————————————————————————————————
pacman::p_load(tidyverse, glmnet, ggplot2)

# BLOCO 1: DEFINIÇÃO DE FUNÇÕES ------------------------------------------------

# 1. FUNÇÕES DE TRANSFORMAÇÃO DE BASE 
logit <- function(p) {
  p <- pmax(pmin(p, 0.999), 0.001)
  log(p / (1 - p))
}
inv_logit <- function(x) exp(x) / (1 + exp(x))

# 2. ORQUESTRADOR DO MODELO EM CASCATA 
ridge_cascata_bbb <- function(df_hist, df_novo) {
  
  if (!"Posicao" %in% names(df_novo)) {
    df_novo <- df_novo %>%
      mutate(Media = (Insta + X + YT) / 3,
             Posicao = rank(desc(Media), ties.method = "first")) %>%
      select(-Media)
  }
  
  cat("\n# ——————————————————————————————————————————————————————————————————\n")
  cat("# DIAGNÓSTICO DO MODELO B EM CASCATA (RIDGE ESTRUTURAL) \n")
  cat("# ——————————————————————————————————————————————————————————————————\n\n")
  
  vars <- c("L_I", "L_X", "L_Y", "Pol")
  
  # NÍVEL 1: O ALVO PRINCIPAL (QUEM SAI) 
  prep_n1 <- function(df) {
    df %>% filter(Posicao == 1) %>% mutate(
      I_p = Insta/100, X_p = X/100, Y_p = YT/100,
      Pol = abs(I_p - X_p) * Y_p,
      L_I = logit(I_p), L_X = logit(X_p), L_Y = logit(Y_p)
    )
  }
  
  h1 <- prep_n1(df_hist) %>% mutate(L_Alvo = logit(Oficial/100))
  n1 <- prep_n1(df_novo)
  
  X_tr1 <- as.matrix(h1[vars]); y_tr1 <- h1$L_Alvo
  X_te1 <- as.matrix(n1[vars])
  
  set.seed(42)
  cv1 <- cv.glmnet(x=X_tr1, y=y_tr1, alpha=0, nfolds=nrow(X_tr1), 
                   standardize=TRUE)
  coef1 <- coef(cv1, s="lambda.min")
  
  pred_h1 <- as.numeric(predict(cv1, newx=X_tr1, s="lambda.min"))
  pred_h1 <- inv_logit(pred_h1) * 100
  pred_n1 <- as.numeric(predict(cv1, newx=X_te1, s="lambda.min"))
  pred_n1 <- inv_logit(pred_n1) * 100
  
  h1$Predito <- pred_h1
  h1$Residuo <- h1$Oficial - h1$Predito
  
  cat(">>> NÍVEL 1: COEFICIENTES RIDGE (1º COLOCADO)\n")
  cat(sprintf(" - Intercepto : %+.4f\n", coef1[1,1]))
  cat(sprintf(" - Logit(Inst): %+.4f\n", coef1[2,1]))
  cat(sprintf(" - Logit(X)   : %+.4f\n", coef1[3,1]))
  cat(sprintf(" - Logit(YT)  : %+.4f\n", coef1[4,1]))
  cat(sprintf(" - Polarização: %+.4f\n\n", coef1[5,1]))
  
  cat(">>> NÍVEL 1: BACKTESTING (ERRO HISTÓRICO)\n")
  tab_h1 <- h1 %>% select(Paredao, Candidato, Oficial, Predito, Residuo) %>%
    mutate(across(where(is.numeric), ~round(.x, 2)))
  print(as.data.frame(tab_h1))
  cat(sprintf("RMSE Nível 1: %.3f p.p.\n\n", sqrt(mean(h1$Residuo^2))))
  
  # NÍVEL 2: A DISPUTA RESIDUAL (SHARE 2º VS 3º) 
  prep_n2 <- function(df) {
    if (!"Paredao" %in% names(df)) df$Paredao <- "Novo"
    df %>% filter(Posicao %in% c(2,3)) %>%
      arrange(Paredao, Posicao) %>% group_by(Paredao) %>%
      summarise(
        Cand_2 = Candidato[1],
        Alvo_S = if("Oficial" %in% names(.)) Oficial[1]/sum(Oficial) else NA,
        I_S = Insta[1]/sum(Insta), X_S = X[1]/sum(X), Y_S = YT[1]/sum(YT),
        .groups = "drop"
      ) %>% mutate(
        Pol = abs(I_S - X_S) * Y_S,
        L_I = logit(I_S), L_X = logit(X_S), L_Y = logit(Y_S)
      )
  }
  
  h2 <- prep_n2(df_hist) %>% mutate(L_Alvo = logit(Alvo_S))
  n2 <- prep_n2(df_novo)
  
  X_tr2 <- as.matrix(h2[vars]); y_tr2 <- h2$L_Alvo
  X_te2 <- as.matrix(n2[vars])
  
  set.seed(42)
  cv2 <- cv.glmnet(x=X_tr2, y=y_tr2, alpha=0, nfolds=nrow(X_tr2), 
                   standardize=TRUE)
  coef2 <- coef(cv2, s="lambda.min")
  
  pred_h2 <- as.numeric(predict(cv2, newx=X_tr2, s="lambda.min"))
  pred_h2 <- inv_logit(pred_h2) * 100
  pred_n2 <- as.numeric(predict(cv2, newx=X_te2, s="lambda.min"))
  pred_n2 <- inv_logit(pred_n2) * 100
  
  h2$Oficial_S <- h2$Alvo_S * 100
  h2$Predito_S <- pred_h2
  h2$Residuo <- h2$Oficial_S - h2$Predito_S
  
  cat(">>> NÍVEL 2: COEFICIENTES RIDGE (SHARE DO 2º COLOCADO)\n")
  cat(sprintf(" - Intercepto : %+.4f\n", coef2[1,1]))
  cat(sprintf(" - Logit(Inst): %+.4f\n", coef2[2,1]))
  cat(sprintf(" - Logit(X)   : %+.4f\n", coef2[3,1]))
  cat(sprintf(" - Logit(YT)  : %+.4f\n", coef2[4,1]))
  cat(sprintf(" - Polarização: %+.4f\n\n", coef2[5,1]))
  
  cat(">>> NÍVEL 2: BACKTESTING (ERRO HISTÓRICO DO SHARE)\n")
  tab_h2 <- h2 %>% select(Paredao, Cand_2, Oficial_S, Predito_S, Residuo) %>%
    mutate(across(where(is.numeric), ~round(.x, 2)))
  print(as.data.frame(tab_h2))
  cat(sprintf("RMSE Nível 2: %.3f p.p.\n\n", sqrt(mean(h2$Residuo^2))))
  
  # RECOMPOSIÇÃO DO SIMPLEX (PROJEÇÃO FINAL) 
  cands <- df_novo %>% arrange(Posicao) %>% pull(Candidato)
  
  p_1 <- pred_n1
  p_2 <- (100 - p_1) * (pred_n2 / 100)
  p_3 <- (100 - p_1) * (1 - (pred_n2 / 100))
  
  t_final <- tibble(
    Posicao = c(1, 2, 3), Candidato = cands,
    Projecao = round(c(p_1, p_2, p_3), 2)
  )
  
  cat("# ——————————————————————————————————————————————————————————————————\n")
  cat("# PROJEÇÃO OFICIAL (SIMPLEX 100% RECOMPOSTO) \n")
  cat("# ——————————————————————————————————————————————————————————————————\n")
  print(as.data.frame(t_final), row.names = FALSE)
  cat("\n")
  
  t_coef <- tibble(
    Variavel = rownames(coef1),
    Nivel_1_Alvo = as.numeric(coef1),
    Nivel_2_Share = as.numeric(coef2)
  )
  
  invisible(list(N1 = cv1, N2 = cv2, Coeficientes = t_coef, 
                 Backtest_N1 = h1, Final = t_final))
}

# 3. RASTREADOR DE ESTABILIDADE ESTRUTURAL 
auditar_coeficientes <- function(res_modelo) {
  cat("\n# AUDITORIA DE ESTABILIDADE (RIDGE PENALIZADO) ———————————————————\n")
  print(as.data.frame(res_modelo$Coeficientes %>% 
                        mutate(across(where(is.numeric), ~round(.x, 4)))))
  cat("\n[!] Variações bruscas de sinal indicam quebra do padrão histórico.\n")
}

# 4. VISUALIZAÇÃO DINÂMICA (PROJEÇÃO DO PAREDÃO) 
plotar_projecao_bbb <- function(res_modelo, num_paredao) {
  df_plot <- res_modelo$Final %>%
    mutate(Candidato = factor(Candidato, levels = Candidato))
  
  erro_max <- max(abs(res_modelo$Backtest_N1$Residuo), na.rm = TRUE)
  data_atual <- format(Sys.time(), "%A, %Hh%M")
  teto_y <- max(df_plot$Projecao) + 15 
  
  p <- ggplot(df_plot, aes(x = Candidato, y = Projecao, fill = Candidato)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    geom_text(
      aes(label = paste0(scales::number(Projecao, decimal.mark = ",", 
                                        accuracy = 0.01), "%")),
      size = 8, vjust = -0.5, fontface = "bold", color = "gray20"
    ) +
    labs(
      title = sprintf("BBB 26 - %dº Paredão", num_paredao),
      subtitle = "Estimativa Estrutural de Votos para Eliminação",
      caption = sprintf(
        "Erro retroativo máximo histórico do eliminado: %.2f p.p.\nAtualização: %s", 
        erro_max, data_atual
      ),
      x = "", y = ""
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.caption = element_text(size = 10, face = "italic", hjust = 1, 
                                  color = "gray40"),
      plot.subtitle = element_text(size = 16, margin = margin(b = 20), 
                                   color = "gray30"),
      plot.title = element_text(size = 22, face = "bold"),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 20, face = "bold"),
      panel.grid = element_blank()
    ) +
    scale_y_continuous(limits = c(0, teto_y)) +
    scale_fill_brewer(palette = "Pastel1")
  
  print(p)
}


# BLOCO 2: ALIMENTAÇÃO DE DADOS -----------------------------------------------
historico_empilhado <- bind_rows(
  tibble(Paredao="P2", Candidato="Matheus",   
         Posicao=1, Oficial=79.48, Insta=77.81, X=80.74, YT=66.42),
  tibble(Paredao="P2", Candidato="Leandro",   
         Posicao=2, Oficial=15.55, Insta=12.63, X=14.12, YT=20.54),
  tibble(Paredao="P2", Candidato="Brigido",   
         Posicao=3, Oficial=4.97,  Insta=9.56,  X=5.14,  YT=13.04),
  
  tibble(Paredao="P3", Candidato="Brigido",   
         Posicao=1, Oficial=77.88, Insta=73.20, X=74.55, YT=58.40),
  tibble(Paredao="P3", Candidato="Leandro",   
         Posicao=2, Oficial=12.04, Insta=16.40, X=21.19, YT=23.60),
  tibble(Paredao="P3", Candidato="Ana Paula", 
         Posicao=3, Oficial=10.08, Insta=10.40, X=4.27,  YT=18.00),
  
  tibble(Paredao="P4", Candidato="Sarah",     
         Posicao=1, Oficial=69.13, Insta=69.33, X=66.04, YT=50.24),
  tibble(Paredao="P4", Candidato="Babu",      
         Posicao=2, Oficial=28.49, Insta=23.07, X=25.44, YT=42.92),
  tibble(Paredao="P4", Candidato="Sol",       
         Posicao=3, Oficial=2.38,  Insta=7.60,  X=8.52,  YT=6.84),
  
  tibble(Paredao="P5", Candidato="Marcelo",   
         Posicao=1, Oficial=68.56, Insta=65.25, X=78.25, YT=66.00),
  tibble(Paredao="P5", Candidato="Samira",    
         Posicao=2, Oficial=16.25, Insta=17.08, X=10.27, YT=20.09),
  tibble(Paredao="P5", Candidato="Solange",   
         Posicao=3, Oficial=15.39, Insta=17.67, X=11.47, YT=13.95),
  
  tibble(Paredao="P6", Candidato="Maxiane",   
         Posicao=1, Oficial=63.21, Insta=65.07, X=68.20, YT=52.17),
  tibble(Paredao="P6", Candidato="Milena",    
         Posicao=2, Oficial=36.11, Insta=33.27, X=26.93, YT=42.55),
  tibble(Paredao="P6", Candidato="Chaiany",   
         Posicao=3, Oficial=0.68,  Insta=1.67,  X=4.86,  YT=5.33),
  
  tibble(Paredao="P7", Candidato="Breno",     
         Posicao=1, Oficial=54.66, Insta=53.43, X=48.99, YT=43.55),
  tibble(Paredao="P7", Candidato="Cowboy",    
         Posicao=2, Oficial=43.12, Insta=41.29, X=42.72, YT=53.14),
  tibble(Paredao="P7", Candidato="Jordana",   
         Posicao=3, Oficial=2.22,  Insta=5.29,  X=8.29,  YT=3.25),
  
  tibble(Paredao="P8", Candidato="Babu",      
         Posicao=1, Oficial=68.62, Insta=70.33, X=77.77, YT=57.63),
  tibble(Paredao="P8", Candidato="Milena",    
         Posicao=2, Oficial=30.91, Insta=27.93, X=13.49, YT=38.24),
  tibble(Paredao="P8", Candidato="Chaiany",   
         Posicao=3, Oficial=0.47,  Insta=1.73,  X=8.75,  YT=4.58),
  
  tibble(Paredao="P9", Candidato="Breno",     
         Posicao=1, Oficial=58.96, Insta=55.38, X=63.33, YT=48.16),
  tibble(Paredao="P9", Candidato="Ana Paula", 
         Posicao=2, Oficial=25.17, Insta=16.13, X=13.84, YT=31.32),
  tibble(Paredao="P9", Candidato="Leandro",   
         Posicao=3, Oficial=15.87, Insta=28.50, X=22.82, YT=20.18)
)

novo_paredao_empilhado <- tibble(
  Candidato = c("Jonas", "Juliano", "Gabriela"),
  Insta     = c(50.56, 32.06, 17.38),
  X         = c(55.92, 38.19, 5.09),
  YT        = c(35.87, 51.05, 13.52)
)



# BLOCO 3: EXECUÇÃO DO FLUXO -------------------------------------------------

# 1. Roda a matemática
resultado_modelo_B <- ridge_cascata_bbb(historico_empilhado, 
                                        novo_paredao_empilhado)

# 2. Imprime a tabela de auditoria
auditar_coeficientes(resultado_modelo_B)

# 3. Desenha o gráfico final no RStudio (Ajuste o número do paredão aqui)
plotar_projecao_bbb(resultado_modelo_B, num_paredao = 10)



