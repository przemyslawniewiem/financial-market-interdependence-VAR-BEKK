library(quantmod)
library(vars)
library(urca)
library(tseries)
library(lmtest)
library(sandwich)
library(ggplot2)
library(gridExtra)
library(moments)
library(tidyr)
library(dplyr)
library(ggcorrplot)

# Wczytanie danych
spx <- read.csv("^spx_d.csv", stringsAsFactors = FALSE)
ndq <- read.csv("^ndq_d.csv", stringsAsFactors = FALSE)
gold <- read.csv("xauusd_d.csv", stringsAsFactors = FALSE)
bond <- read.csv("10yusy_b_d.csv", stringsAsFactors = FALSE)
btc <- read.csv("btc_v_d.csv", stringsAsFactors = FALSE)

# Przygotowanie danych
spx$Data <- as.Date(spx$Data)
ndq$Data <- as.Date(ndq$Data)
gold$Data <- as.Date(gold$Data)
bond$Data <- as.Date(bond$Data)
btc$Data <- as.Date(btc$Data)

# Merge wszystkich szeregów
data <- merge(spx[,c("Data","Zamkniecie")], ndq[,c("Data","Zamkniecie")], by.x="Data", by.y="Data", all=FALSE)
names(data) <- c("Date","SPX","NDQ")
data <- merge(data, gold[,c("Data","Zamkniecie")], by.x="Date", by.y="Data", all=FALSE)
names(data)[4] <- "GOLD"
data <- merge(data, bond[,c("Data","Zamkniecie")], by.x="Date", by.y="Data", all=FALSE)
names(data)[5] <- "BOND"
data <- merge(data, btc[,c("Data","Zamkniecie")], by.x="Date", by.y="Data", all=FALSE)
names(data)[6] <- "BTC"

# Logarytmiczne stopy zwrotu
returns <- data.frame(Date = data$Date[-1])
returns$SPX <- diff(log(data$SPX))
returns$NDQ <- diff(log(data$NDQ))
returns$GOLD <- diff(log(data$GOLD))
returns$BOND <- diff(log(data$BOND))
returns$BTC <- diff(log(data$BTC))

# Usunięcie braków
returns <- na.omit(returns)

# ===============================================
# NOWE: STATYSTYKI OPISOWE
# ===============================================
cat("\n=== STATYSTYKI OPISOWE ===\n")
returns_ts <- returns[,-1]

stats_table <- data.frame(
  Mean = colMeans(returns_ts) * 252,  # Zannualizowana średnia
  SD = apply(returns_ts, 2, sd) * sqrt(252),  # Zannualizowane odch. std.
  Skewness = apply(returns_ts, 2, skewness),
  Kurtosis = apply(returns_ts, 2, kurtosis),
  Min = apply(returns_ts, 2, min),
  Max = apply(returns_ts, 2, max)
)
print(round(stats_table, 4))

# ===============================================
# NOWE: WYKRES SKUMULOWANYCH ZWROTÓW
# ===============================================
cat("\n=== WYKRES SKUMULOWANYCH ZWROTÓW ===\n")
cum_returns <- apply(returns_ts, 2, cumsum)

# Wykres za pomocą base R (bardziej stabilny)
plot(returns$Date, cum_returns[,1], type="l", col=1, lwd=2,
     ylim=range(cum_returns), 
     xlab="Data", ylab="Skumulowany zwrot",
     main="Skumulowane log-stopy zwrotu (2018-2025)")
abline(h=0, lty=2, col="black")
for(i in 2:ncol(cum_returns)) {
  lines(returns$Date, cum_returns[,i], col=i, lwd=2)
}
legend("topleft", legend=colnames(cum_returns), col=1:ncol(cum_returns), 
       lwd=2, bty="n", cex=0.8)

# ===============================================
# NOWE: MACIERZ KORELACJI (HEATMAP)
# ===============================================
cat("\n=== MACIERZ KORELACJI ===\n")
corr_matrix <- cor(returns_ts)
print(round(corr_matrix, 4))

# Heatmapa za pomocą base R
library(RColorBrewer)
col_palette <- colorRampPalette(c("#E46726", "white", "#6D9EC1"))(100)
heatmap(corr_matrix, 
        Rowv=NA, Colv=NA, 
        col=col_palette,
        scale="none",
        main="Macierz korelacji stóp zwrotu",
        margins=c(8,8),
        cexRow=1, cexCol=1)

cat("\n=== ANALIZA CAŁEGO OKRESU ===\n")

# Stacjonarność - testy ADF
cat("\nTesty stacjonarności (ADF):\n")
adf_spx <- adf.test(returns$SPX)
adf_ndq <- adf.test(returns$NDQ)
adf_gold <- adf.test(returns$GOLD)
adf_bond <- adf.test(returns$BOND)
adf_btc <- adf.test(returns$BTC)

cat("SPX p-value:", adf_spx$p.value, "\n")
cat("NDQ p-value:", adf_ndq$p.value, "\n")
cat("GOLD p-value:", adf_gold$p.value, "\n")
cat("BOND p-value:", adf_bond$p.value, "\n")
cat("BTC p-value:", adf_btc$p.value, "\n")

# Dobór opóźnienia - ZMIENIONO NA HQ
lag_selection <- VARselect(returns_ts, lag.max=20, type="const")
cat("\nKryteria doboru opóźnienia:\n")
print(lag_selection$selection)
optimal_lag <- lag_selection$selection["HQ(n)"]  # ZMIANA: HQ zamiast AIC

# Estymacja VAR
var_model <- VAR(returns_ts, p=optimal_lag, type="const")
cat("\nModel VAR z", optimal_lag, "opóźnieniami\n")

# Diagnostyka reszt
cat("\nDiagnostyka reszt:\n")
serial_test <- serial.test(var_model, lags.pt=16, type="PT.asymptotic")
cat("Test autokorelacji (Portmanteau) p-value:", serial_test$serial$p.value, "\n")

arch_test <- arch.test(var_model, lags.multi=5)
cat("Test ARCH p-value:", arch_test$arch.mul$p.value, "\n")

normality_test <- normality.test(var_model, multivariate.only=TRUE)
cat("Test normalności (JB) p-value:", normality_test$jb.mul$JB$p.value, "\n")

# Test stabilności VAR
cat("\nTest stabilności VAR:\n")
var_roots <- roots(var_model)
cat("Pierwiastki charakterystyczne:\n")
print(round(var_roots, 4))

if(all(var_roots < 1)) {
  cat("Model jest stabilny\n")
} else {
  cat("Model może być niestabilny\n")
}

# ===============================================
# ROBUST GRANGER TESTS (Newey-West HAC)
# ===============================================
cat("\n=== TESTY PRZYCZYNOWOŚCI GRANGERA (ROBUST - Newey-West HAC) ===\n")

run_robust_granger_matrix <- function(data, p_lag) {
  tickers <- colnames(data)
  n <- length(tickers)
  
  # Macierz wyników
  p_matrix <- matrix(NA, nrow = n, ncol = n)
  rownames(p_matrix) <- tickers
  colnames(p_matrix) <- tickers
  
  core_data <- as.matrix(data)
  n_rows <- nrow(core_data)
  
  for (i in 1:n) { # Target (Skutek)
    for (j in 1:n) { # Source (Przyczyna)
      if (i != j) {
        y_vec <- core_data[(p_lag + 1):n_rows, i]
        
        target_lags <- embed(core_data[, i], p_lag + 1)[, -1, drop = FALSE]
        source_lags <- embed(core_data[, j], p_lag + 1)[, -1, drop = FALSE]
        
        # Model Pełny (Unrestricted)
        df_unres <- data.frame(Y = y_vec, target_lags, source_lags)
        model_unres <- lm(Y ~ ., data = df_unres)
        
        # Model Ograniczony (Restricted)
        df_res <- data.frame(Y = y_vec, target_lags)
        model_res <- lm(Y ~ ., data = df_res)
        
        # Test Walda z Newey-West
        wald_result <- waldtest(model_res, model_unres, vcov = NeweyWest(model_unres, lag = p_lag, prewhite = FALSE))
        
        p_val <- wald_result$`Pr(>F)`[2]
        p_matrix[i, j] <- round(p_val, 4)
      }
    }
  }
  
  return(p_matrix)
}

granger_matrix <- run_robust_granger_matrix(returns_ts, optimal_lag)

cat("\nMacierz przyczynowości (wiersz powoduje kolumnę) - Robust HAC:\n")
print(granger_matrix)

# ===============================================
# ANALIZA PODOKRESÓW
# ===============================================
analyze_subperiod <- function(returns, start_date, end_date, period_name){
  
  cat("\n\n=== ANALIZA PODOKRESU:", period_name, "===\n")
  
  sub_returns <- returns[returns$Date >= start_date & returns$Date <= end_date, ]
  sub_returns_ts <- sub_returns[,-1]
  
  cat("Liczba obserwacji:", nrow(sub_returns), "\n")
  
  # Stacjonarność
  cat("\nTesty stacjonarności (ADF):\n")
  adf_results <- sapply(sub_returns_ts, function(x) adf.test(x)$p.value)
  print(round(adf_results, 4))
  
  # Dobór opóźnienia - ZMIANA: HQ zamiast AIC
  lag_sel <- VARselect(sub_returns_ts, lag.max=20, type="const")
  cat("\nKryteria doboru opóźnienia:\n")
  print(lag_sel$selection)
  opt_lag <- as.numeric(lag_sel$selection["HQ(n)"])  # ZMIANA: HQ
  
  # VAR
  var_sub <- VAR(sub_returns_ts, p=opt_lag, type="const")
  cat("\nModel VAR z", opt_lag, "opóźnieniami\n")
  
  # Diagnostyka
  cat("\nDiagnostyka reszt:\n")
  serial <- serial.test(var_sub, lags.pt=16, type="PT.asymptotic")
  cat("Portmanteau p-value:", serial$serial$p.value, "\n")
  
  arch <- arch.test(var_sub, lags.multi=5)
  cat("ARCH p-value:", arch$arch.mul$p.value, "\n")
  
  norm <- normality.test(var_sub, multivariate.only=TRUE)
  cat("Normalność p-value:", norm$jb.mul$JB$p.value, "\n")
  
  # Stabilność
  cat("\nStabilność:\n")
  roots_sub <- roots(var_sub)
  print(round(roots_sub, 4))
  if(all(roots_sub < 1)) {
    cat("Model stabilny\n")
  } else {
    cat("Model niestabilny\n")
  }
  
  # ZMIANA: Robust Granger z Newey-West HAC
  cat("\nMacierz przyczynowości (Robust HAC):\n")
  granger_sub <- run_robust_granger_matrix(sub_returns_ts, opt_lag)
  print(granger_sub)
  
  return(list(var=var_sub, granger=granger_sub))
}

# Podokresy
period1 <- analyze_subperiod(returns, "2018-03-01", "2020-01-31", "Marzec 2018 - Styczeń 2020")
period2 <- analyze_subperiod(returns, "2020-03-01", "2022-01-31", "Marzec 2020 - Styczeń 2022")
period3 <- analyze_subperiod(returns, "2022-03-01", "2024-01-31", "Marzec 2022 - Styczeń 2024")

cat("\n\n=== PODSUMOWANIE WYNIKÓW ===\n")

cat("\n--- CAŁY OKRES (2018-2025) ---\n")
print(round(granger_matrix, 4))

cat("\n--- OKRES 1: Marzec 2018 - Styczeń 2020 ---\n")
print(round(period1$granger, 4))

cat("\n--- OKRES 2: Marzec 2020 - Styczeń 2022 ---\n")
print(round(period2$granger, 4))

cat("\n--- OKRES 3: Marzec 2022 - Styczeń 2024 ---\n")
print(round(period3$granger, 4))

# WYKRES 
plot_heatmap <- function(mat, title) {
  mat_plot <- mat
  
  mat_plot[is.na(mat_plot)] <- -1
  
  mat_plot_rot <- mat_plot[nrow(mat_plot):1, ]
  
  breaks <- c(-1, 0, 0.01, 0.05, 0.1, 1)
  colors <- c("white", "darkred", "orange", "yellow", "lightgray")
  
  image(1:ncol(mat_plot), 1:nrow(mat_plot), t(mat_plot_rot),
        col = colors,
        breaks = breaks,
        xlab = "Do", 
        ylab = "Od", 
        main = title,
        axes = FALSE)
  
  # Osie
  axis(1, at = 1:ncol(mat_plot), labels = colnames(mat_plot), las = 2)
  axis(2, at = 1:nrow(mat_plot), labels = rev(rownames(mat_plot)), las = 1)
  
  # Siatka
  abline(h = 0.5:(nrow(mat_plot) + 0.5), col = "gray30", lty = 1, lwd = 0.8)
  abline(v = 0.5:(ncol(mat_plot) + 0.5), col = "gray30", lty = 1, lwd = 0.8)
  
  # Dodaj p-values jako tekst
  for(i in 1:nrow(mat)) {
    for(j in 1:ncol(mat)) {
      if(!is.na(mat[i, j])) {
        p_val <- sprintf("%.3f", mat[i, j])
        # Dobór koloru tekstu w zależności od tła
        text_col <- if(mat[i, j] < 0.01) "white" else "black"
        text(j, nrow(mat) - i + 1, p_val, 
             col = text_col, cex = 0.75, font = 2)
      }
    }
  }
  
  # Legenda 
  par(xpd = TRUE)
  legend(x = ncol(mat) + 0.8, y = nrow(mat), 
         legend = c("p < 0.01", "0.01 ≤ p < 0.05", "0.05 ≤ p < 0.10", "p ≥ 0.10", "N/A"),
         fill = c("darkred", "orange", "yellow", "lightgray", "white"),
         bty = "n", cex = 0.65, xjust = 0, yjust = 1)
  par(xpd = FALSE)
}

# Wykresy  
par(mfrow=c(2,2), mar=c(5,5,3,7), oma=c(0,0,0,0))
plot_heatmap(granger_matrix, "Cały okres (2018-2025)")
plot_heatmap(period1$granger, "Okres 1 (2018-2020)")
plot_heatmap(period2$granger, "Okres 2 (2020-2022)")
plot_heatmap(period3$granger, "Okres 3 (2022-2024)")
par(mfrow=c(1,1), mar=c(5,4,4,2))

# ===============================================================================
# CZĘŚĆ B: ANALIZA SPILLOVER Z WYKORZYSTANIEM BEKK(1,1)
# ===============================================================================

library(BEKKs)
library(xts)

cat("\n\n")
cat("================================================================================\n")
cat("                    CZĘŚĆ B: SPILLOVER BEKK(1,1)                                \n")
cat("================================================================================\n")

# ===============================================================================
# 1. FUNKCJA: EKSTRAKCJA RESZT Z VAR
# ===============================================================================

reszty_z_var <- function(returns, start_date, end_date, period_name, lag_list) {
  cat("\n>>> Ekstrakcja reszt VAR dla okresu:", period_name, "\n")
  
  sub <- returns[returns$Date >= start_date & returns$Date <= end_date, ]
  sub_ts <- sub[,-1]
  
  # Dobór opóźnienia (HQ)
  lag_sel <- VARselect(sub_ts, lag.max=20, type="const")
  p_lag <- as.numeric(lag_sel$selection["HQ(n)"])
  
  if (is.null(p_lag) || is.na(p_lag) || p_lag < 1) p_lag <- 1
  
  cat("  - Liczba obserwacji:", nrow(sub_ts), "\n")
  cat("  - Opóźnienie VAR (HQ):", p_lag, "\n")
  
  # Estymacja VAR
  var_model <- VAR(sub_ts, p = p_lag, type = "const")
  
  # Ekstrakcja reszt
  eps_mat <- residuals(var_model)
  eps_mat <- na.omit(eps_mat)
  
  # Dopasowanie indeksu dat
  idx <- sub$Date[(nrow(sub) - nrow(eps_mat) + 1):nrow(sub)]
  eps_xts <- xts(eps_mat, order.by = idx)
  colnames(eps_xts) <- colnames(sub_ts)
  
  cat("  - Liczba reszt:", nrow(eps_xts), "\n")
  
  list(
    eps = eps_xts,
    var_model = var_model,
    p_lag = p_lag
  )
}

# ===============================================================================
# 2. FUNKCJA: ESTYMACJA BEKK(1,1)
# ===============================================================================

dopasuj_bekk <- function(eps_xts, period_name,
                         max_iter = 150,
                         crit = 1e-9) {
  
  cat("\n>>> Estymacja BEKK(1,1) dla okresu:", period_name, "\n")
  
  # Specyfikacja modelu BEKK(1,1) - PEŁNY BEKK (nie diagonal!)
  spec <- bekk_spec(
    model = list(type = "bekk", asymmetric = FALSE),
    init_values = "simple"  # Stabilniejsza inicjalizacja
  )
  
  cat("  - Rozpoczynam optymalizację...\n")
  
  # Estymacja z większą liczbą iteracji i ściślejszym kryterium
  fit <- tryCatch({
    bekk_fit(spec, eps_xts, 
             QML_t_ratios = TRUE, 
             max_iter = max_iter, 
             crit = crit)
  }, error = function(e) {
    cat("  [BŁĄD]", e$message, "\n")
    cat("  >>> Próbuję z luźniejszym kryterium zbieżności...\n")
    
    # Próba z łagodniejszym kryterium
    tryCatch({
      bekk_fit(spec, eps_xts, 
               QML_t_ratios = TRUE, 
               max_iter = max_iter * 2, 
               crit = 1e-6)
    }, error = function(e2) {
      cat("  [BŁĄD KRYTYCZNY] Nie udało się dopasować BEKK.\n")
      return(NULL)
    })
  })
  
  if (is.null(fit)) {
    stop("BEKK nie zbiegł dla okresu ", period_name)
  }
  
  cat("  - BEKK_valid:", fit$BEKK_valid, "\n")
  
  if (!is.null(fit$Portmanteau.test)) {
    cat("  - Portmanteau p-value:", round(fit$Portmanteau.test$p.value, 4), "\n")
  }
  
  return(fit)
}

# ===============================================================================
# 3. FUNKCJA: TESTY SPILLOVER (A i G) + ŁĄCZNY TEST CHI-KWADRAT
# ===============================================================================

test_spillover <- function(fit_bekk, nazwy, period_name) {
  
  cat("\n>>> Testy spillover dla okresu:", period_name, "\n")
  
  A  <- fit_bekk$A
  G  <- fit_bekk$G
  tA <- fit_bekk$A_t
  tG <- fit_bekk$G_t
  
  if (is.null(tA) || is.null(tG)) {
    stop("Brak A_t lub G_t w obiekcie BEKK. Sprawdź QML_t_ratios=TRUE.")
  }
  
  k <- length(nazwy)
  
  wyniki <- data.frame(
    Okres = character(),
    Do = character(),
    Z  = character(),
    A  = numeric(),
    Wartosc_p_A = numeric(),
    G  = numeric(),
    Wartosc_p_G = numeric(),
    Stat_chi2 = numeric(),
    Wartosc_p_AG = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:k) {        # i = Do (cel)
    for (j in 1:k) {      # j = Z  (źródło)
      if (i == j) next
      
      # UWAGA: kierunek Z -> Do w BEKKs to element [j,i]
      tA_ji <- tA[j, i]
      tG_ji <- tG[j, i]
      
      # P-values z rozkładu normalnego (dwustronny test)
      pA <- 2 * (1 - pnorm(abs(tA_ji)))
      pG <- 2 * (1 - pnorm(abs(tG_ji)))
      
      # Test łączny chi-kwadrat (df=2)
      chi2 <- tA_ji^2 + tG_ji^2
      pAG  <- 1 - pchisq(chi2, df = 2)
      
      wyniki[nrow(wyniki) + 1, ] <- list(
        Okres = period_name,
        Do = nazwy[i],
        Z  = nazwy[j],
        A  = A[j, i],
        Wartosc_p_A = pA,
        G  = G[j, i],
        Wartosc_p_G = pG,
        Stat_chi2 = chi2,
        Wartosc_p_AG = pAG
      )
    }
  }
  
  # Sortowanie po p-value łącznym
  wyniki <- wyniki[order(wyniki$Wartosc_p_AG), ]
  
  # Podsumowanie
  istotne <- sum(wyniki$Wartosc_p_AG < 0.05, na.rm = TRUE)
  cat("  - Liczba istotnych spillover (p < 0.05):", istotne, "\n")
  cat("  - Minimalna wartość p (test łączny):", 
      round(min(wyniki$Wartosc_p_AG, na.rm = TRUE), 6), "\n")
  
  return(wyniki)
}

# ===============================================================================
# 4. FUNKCJA: HEATMAPA SPILLOVER
# ===============================================================================

plot_spillover_heatmap <- function(tab_spill, nazwy, period_name) {
  
  # Macierz p-values (test łączny A+G)
  M <- matrix(NA_real_, nrow = length(nazwy), ncol = length(nazwy),
              dimnames = list(nazwy, nazwy))
  
  for (r in 1:nrow(tab_spill)) {
    M[tab_spill$Do[r], tab_spill$Z[r]] <- tab_spill$Wartosc_p_AG[r]
  }
  
  # Konwersja do formatu długiego
  df <- as.data.frame(as.table(M))
  names(df) <- c("Do", "Z", "p")
  df <- df[!is.na(df$p), ]
  
  # Transformacja -log10(p) dla lepszej wizualizacji
  df$score <- -log10(df$p)
  
  # Gwiazdki istotności
  df$star <- ""
  df$star[df$p < 0.001] <- "***"
  df$star[df$p >= 0.001 & df$p < 0.01] <- "**"
  df$star[df$p >= 0.01 & df$p < 0.05] <- "*"
  
  df$Do <- factor(df$Do, levels = nazwy)
  df$Z  <- factor(df$Z,  levels = nazwy)
  
  # Wykres (base R)
  mat_plot <- M
  mat_plot[is.na(mat_plot)] <- 1  # NA jako nieistotne
  
  # Odwróć macierz dla wyświetlenia
  mat_plot_rot <- mat_plot[nrow(mat_plot):1, ]
  
  # Paleta kolorów (czerwony = istotne, biały = nieistotne)
  breaks <- c(0, 0.001, 0.01, 0.05, 0.1, 1)
  colors <- c("darkred", "red", "orange", "yellow", "white")
  
  image(1:ncol(mat_plot), 1:nrow(mat_plot), t(mat_plot_rot),
        col = colors,
        breaks = breaks,
        xlab = "Źródło (Z)", 
        ylab = "Cel (Do)", 
        main = paste("Spillover BEKK(1,1) -", period_name),
        axes = FALSE)
  
  axis(1, at = 1:ncol(mat_plot), labels = colnames(mat_plot), las = 2)
  axis(2, at = 1:nrow(mat_plot), labels = rev(rownames(mat_plot)), las = 1)
  
  abline(h = 0.5:(nrow(mat_plot) + 0.5), col = "gray30", lty = 1, lwd = 0.5)
  abline(v = 0.5:(ncol(mat_plot) + 0.5), col = "gray30", lty = 1, lwd = 0.5)
  
  # Dodaj gwiazdki istotności
  for(i in 1:nrow(M)) {
    for(j in 1:ncol(M)) {
      if(!is.na(M[i, j])) {
        star_text <- ""
        if(M[i, j] < 0.001) star_text <- "***"
        else if(M[i, j] < 0.01) star_text <- "**"
        else if(M[i, j] < 0.05) star_text <- "*"
        
        if(star_text != "") {
          text_col <- if(M[i, j] < 0.01) "white" else "black"
          text(j, nrow(M) - i + 1, star_text, 
               col = text_col, cex = 1.2, font = 2)
        }
      }
    }
  }
  
  # Legenda
  par(xpd = TRUE)
  legend(x = ncol(mat_plot) + 0.8, y = nrow(mat_plot), 
         legend = c("p < 0.001 (***)", "0.001 ≤ p < 0.01 (**)", 
                    "0.01 ≤ p < 0.05 (*)", "p ≥ 0.05"),
         fill = c("darkred", "red", "orange", "white"),
         bty = "n", cex = 0.65, xjust = 0, yjust = 1)
  par(xpd = FALSE)
}

# ===============================================================================
# 5. FUNKCJA: WYKRES TOP WARUNKOWYCH KORELACJI
# ===============================================================================

plot_top_condcorr <- function(fit_bekk, tab_spill, top_n = 3, period_name) {
  
  st <- fit_bekk$sigma_t
  if (!inherits(st, "xts")) {
    cat("  [INFO] Brak sigma_t - pomijam wykres korelacji warunkowych\n")
    return(NULL)
  }
  
  # Wybierz top pary
  top_pairs <- tab_spill[1:min(top_n, nrow(tab_spill)), ]
  
  if(nrow(top_pairs) == 0) return(NULL)
  
  # Pomocnicza funkcja do znajdowania kolumny korelacji
  pick_corr_col <- function(st_cols, a, b) {
    ix <- which(grepl("Conditional correlation", st_cols) &
                  grepl(a, st_cols, fixed = TRUE) &
                  grepl(b, st_cols, fixed = TRUE))
    if (length(ix) == 0) return(NA_character_)
    st_cols[ix[1]]
  }
  
  cols <- colnames(st)
  
  # Przygotuj wykres dla każdej pary
  par(mfrow=c(min(top_n, nrow(top_pairs)), 1), mar=c(3,4,2,1))
  
  for (i in 1:nrow(top_pairs)) {
    a <- top_pairs$Do[i]
    b <- top_pairs$Z[i]
    pair_name <- paste(b, "→", a)
    
    colname <- pick_corr_col(cols, a, b)
    if (is.na(colname)) colname <- pick_corr_col(cols, b, a)
    
    if (!is.na(colname)) {
      plot(index(st), as.numeric(st[, colname]),
           type = "l", col = "steelblue", lwd = 1.5,
           xlab = "", ylab = "Korelacja",
           main = paste(pair_name, "| p =", 
                        round(top_pairs$Wartosc_p_AG[i], 4)))
      abline(h = 0, lty = 2, col = "gray50")
      grid()
    }
  }
  
  par(mfrow=c(1,1), mar=c(5,4,4,2))
}

# ===============================================================================
# 6. DEFINICJA OKRESÓW (IDENTYCZNE JAK W CZĘŚCI A)
# ===============================================================================

periods <- list(
  Full = c("2018-01-01", "2025-12-31"),
  P1 = c("2018-03-01", "2020-01-31"),
  P2 = c("2020-03-01", "2022-01-31"),
  P3 = c("2022-03-01", "2024-01-31")
)

# ===============================================================================
# 7. GŁÓWNA PĘTLA ANALIZY
# ===============================================================================

bekk_wyniki <- list()
wszystkie_spillover <- list()

for (n in names(periods)) {
  
  cat("\n")
  cat("================================================================================\n")
  cat("OKRES:", n, "(", periods[[n]][1], "do", periods[[n]][2], ")\n")
  cat("================================================================================\n")
  
  # Krok 1: Ekstrakcja reszt z VAR
  tmp <- tryCatch({
    reszty_z_var(returns, periods[[n]][1], periods[[n]][2], n, NULL)
  }, error = function(e) {
    cat("[BŁĄD] Nie udało się wyekstrahować reszt dla okresu", n, "\n")
    return(NULL)
  })
  
  if (is.null(tmp)) next
  
  # Krok 2: Estymacja BEKK(1,1)
  fit <- tryCatch({
    dopasuj_bekk(tmp$eps, n, max_iter = 150, crit = 1e-9)
  }, error = function(e) {
    cat("[BŁĄD] Nie udało się dopasować BEKK dla okresu", n, "\n")
    cat("Szczegóły:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(fit)) next
  
  # Krok 3: Testy spillover
  nazwy <- colnames(tmp$eps)
  tab_spill <- test_spillover(fit, nazwy, n)
  
  # Zaokrąglenie dla wydruku
  tab_print <- tab_spill
  num_cols <- c("A", "Wartosc_p_A", "G", "Wartosc_p_G", "Stat_chi2", "Wartosc_p_AG")
  tab_print[, num_cols] <- round(tab_print[, num_cols], 6)
  
  cat("\n>>> TOP 10 NAJSILNIEJSZYCH SPILLOVER (test łączny A+G):\n")
  print(head(tab_print, 10))
  
  # Zapisz wyniki
  bekk_wyniki[[n]] <- list(
    var_lag = tmp$p_lag,
    var_model = tmp$var_model,
    eps = tmp$eps,
    bekk = fit,
    spillover = tab_spill
  )
  
  wszystkie_spillover[[n]] <- tab_spill
}

# ===============================================================================
# 8. PODSUMOWANIE ZBIORCZE
# ===============================================================================

cat("\n\n")
cat("================================================================================\n")
cat("                     PODSUMOWANIE ZBIORCZE - CZĘŚĆ B                            \n")
cat("================================================================================\n")

# Tabela: ile istotnych spillover w każdym okresie
podsumowanie <- data.frame(
  Okres = character(),
  N_reszt = integer(),
  BEKK_valid = logical(),
  Min_p_AG = numeric(),
  Liczba_istotnych = integer(),
  stringsAsFactors = FALSE
)

for (ok in names(bekk_wyniki)) {
  sp <- bekk_wyniki[[ok]]$spillover
  sp <- sp[!is.na(sp$Wartosc_p_AG), ]
  
  podsumowanie[nrow(podsumowanie) + 1, ] <- list(
    Okres = ok,
    N_reszt = nrow(bekk_wyniki[[ok]]$eps),
    BEKK_valid = bekk_wyniki[[ok]]$bekk$BEKK_valid,
    Min_p_AG = min(sp$Wartosc_p_AG, na.rm = TRUE),
    Liczba_istotnych = sum(sp$Wartosc_p_AG < 0.05, na.rm = TRUE)
  )
}

cat("\n>>> PODSUMOWANIE OKRESÓW:\n")
print(podsumowanie)

# Zbiorcza tabela istotnych spillover
istotne_spill <- do.call(rbind, lapply(names(bekk_wyniki), function(n) {
  x <- bekk_wyniki[[n]]$spillover
  x[!is.na(x$Wartosc_p_AG) & x$Wartosc_p_AG < 0.05, ]
}))

cat("\n>>> WSZYSTKIE ISTOTNE SPILLOVER (p < 0.05):\n")
if (is.null(istotne_spill) || nrow(istotne_spill) == 0) {
  cat("Brak istotnych spillover dla progu 0.05.\n")
} else {
  istotne_spill <- istotne_spill[order(istotne_spill$Okres, istotne_spill$Wartosc_p_AG), ]
  istotne_print <- istotne_spill
  num_cols <- sapply(istotne_print, is.numeric)
  istotne_print[, num_cols] <- round(istotne_print[, num_cols], 6)
  print(istotne_print)
}

# ===============================================================================
# 9. WIZUALIZACJE
# ===============================================================================

cat("\n>>> Generowanie wizualizacji...\n")

# Heatmapy spillover
par(mfrow=c(2,2), mar=c(5,5,3,7), oma=c(0,0,2,0))

for (ok in names(bekk_wyniki)) {
  sp <- bekk_wyniki[[ok]]$spillover
  nazwy <- colnames(bekk_wyniki[[ok]]$eps)
  plot_spillover_heatmap(sp, nazwy, ok)
}

mtext("Heatmapy Spillover BEKK(1,1) - Test łączny (A+G)", 
      outer = TRUE, cex = 1.3, font = 2)

par(mfrow=c(1,1), mar=c(5,4,4,2), oma=c(0,0,0,0))

# Wykresy warunkowych korelacji dla każdego okresu
for (ok in names(bekk_wyniki)) {
  if (sum(bekk_wyniki[[ok]]$spillover$Wartosc_p_AG < 0.05, na.rm = TRUE) > 0) {
    cat("\n>>> Wykresy korelacji warunkowych dla okresu:", ok, "\n")
    plot_top_condcorr(bekk_wyniki[[ok]]$bekk, 
                      bekk_wyniki[[ok]]$spillover, 
                      top_n = 3, 
                      period_name = ok)
  }
}




