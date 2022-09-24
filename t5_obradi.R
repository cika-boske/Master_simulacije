moci_ampute_t5 <- rep(0, length(proporcije))
moci_impute_mean_t5 <- rep(0, length(proporcije))
moci_impute_svd_t5 <- rep(0, length(proporcije))
moci_impute_3nn_t5 <- rep(0, length(proporcije))
moci_impute_6nn_t5 <- rep(0, length(proporcije))
moci_impute_10nn_t5 <- rep(0, length(proporcije))
moci_impute_pmm_t5 <- rep(0, length(proporcije))
moci_impute_cart_t5 <- rep(0, length(proporcije))



for (i in 1:length(proporcije)) 
{
  T1_mean <- rep(0, N)
  T1_svd <- rep(0, N)
  T1_3nn <- rep(0, N)
  T1_6nn <- rep(0, N)
  T1_10nn <- rep(0, N)
  T1_pmm <- rep(0, N)
  T1_cart <- rep(0, N)
  T2 <- rep(0, N)
  for (j in 1:N) 
  {
    uzorak <- rmvt(n,
                   mu = rep(0, 2),
                   sigma = kovmat,
                   v = 5)
    uzorak_izbaceno <- delete_MCAR(uzorak, p = proporcije[i])
    uzorak_amputirano <- na.omit(uzorak_izbaceno)
    
    means <- c(mean(uzorak_izbaceno[, 1], na.rm = TRUE), mean(uzorak_izbaceno[, 2], na.rm = TRUE))
    sds <- c(sd(uzorak_izbaceno[, 1], na.rm = TRUE), sd(uzorak_izbaceno[, 2], na.rm = TRUE))
    uzorak_skalirano <- uzorak_izbaceno
    uzorak_skalirano[, 1] <- scale(uzorak_izbaceno[, 1])
    uzorak_skalirano[, 2] <- scale(uzorak_skalirano[, 2])
    
    # statistika za amputaciju
    T2[j] <- BHEP(uzorak_amputirano, a = 1)
    
    # statistika za mean
    uzorak_ubaceno_mean <- impute_mean(uzorak_skalirano, type = "columnwise")
    uzorak_ubaceno_mean[, 1] <- sds[1] * uzorak_ubaceno_mean[, 1] + means[1]
    uzorak_ubaceno_mean[, 2] <- sds[2] * uzorak_ubaceno_mean[, 2] + means[2]
    uzorak_ubaceno_mean <- scale(uzorak_ubaceno_mean)
    T1_mean[j] <- BHEP(uzorak_ubaceno_mean, a = 1)
    
    # statistika za SVD
    uzorak_ubaceno_svd <- eimpute(data.matrix(uzorak_skalirano), r = 1)$x.imp
    uzorak_ubaceno_svd[, 1] <- sds[1] * uzorak_ubaceno_svd[, 1] + means[1]
    uzorak_ubaceno_svd[, 2] <- sds[2] * uzorak_ubaceno_svd[, 2] + means[2]
    uzorak_ubaceno_svd <- scale(uzorak_ubaceno_svd)
    T1_svd[j] <- BHEP(uzorak_ubaceno_svd, a = 1)
    
    # statistika za 3nn
    uzorak_ubaceno_3nn <- knn.impute(uzorak_skalirano, k = 3)
    uzorak_ubaceno_3nn[, 1] <- sds[1] * uzorak_ubaceno_3nn[, 1] + means[1]
    uzorak_ubaceno_3nn[, 2] <- sds[2] * uzorak_ubaceno_3nn[, 2] + means[2]
    uzorak_ubaceno_3nn <- scale(uzorak_ubaceno_3nn)
    T1_3nn[j] <- BHEP(uzorak_ubaceno_3nn, a = 1)
    
    # statistika za 6nn
    uzorak_ubaceno_6nn <- knn.impute(uzorak_skalirano, k = 6)
    uzorak_ubaceno_6nn[, 1] <- sds[1] * uzorak_ubaceno_6nn[, 1] + means[1]
    uzorak_ubaceno_6nn[, 2] <- sds[2] * uzorak_ubaceno_6nn[, 2] + means[2]
    uzorak_ubaceno_6nn <- scale(uzorak_ubaceno_6nn)
    T1_6nn[j] <- BHEP(uzorak_ubaceno_6nn, a = 1)
    
    # statistika za 10nn
    uzorak_ubaceno_10nn <- knn.impute(uzorak_skalirano, k = 10)
    uzorak_ubaceno_10nn[, 1] <- sds[1] * uzorak_ubaceno_10nn[, 1] + means[1]
    uzorak_ubaceno_10nn[, 2] <- sds[2] * uzorak_ubaceno_10nn[, 2] + means[2]
    uzorak_ubaceno_10nn <- scale(uzorak_ubaceno_10nn)
    T1_10nn[j] <- BHEP(uzorak_ubaceno_10nn, a = 1)
    
    uzorak_ubaceno_cart <- data.matrix(complete(mice(uzorak_skalirano, m = 1, method = "cart")))
    uzorak_ubaceno_cart[, 1] <- sds[1] * uzorak_ubaceno_cart[, 1] + means[1]
    uzorak_ubaceno_cart[, 2] <- sds[2] * uzorak_ubaceno_cart[, 2] + means[2]
    uzorak_ubaceno_cart <- scale(uzorak_ubaceno_cart)
    T1_cart[j] <- BHEP(uzorak_ubaceno_cart, a = 1)
    
    uzorak_ubaceno_pmm <- data.matrix(complete(mice(uzorak_skalirano, m = 1)))
    uzorak_ubaceno_pmm[, 1] <- sds[1] * uzorak_ubaceno_pmm[, 1] + means[1]
    uzorak_ubaceno_pmm[, 2] <- sds[2] * uzorak_ubaceno_pmm[, 2] + means[2]
    uzorak_ubaceno_pmm <- scale(uzorak_ubaceno_pmm)
    T1_pmm[j] <- BHEP(uzorak_ubaceno_pmm, a = 1)
  }
  n_new <- nrow(uzorak_amputirano)
  moci_ampute_t5[i] <- moc_testa(nulte1[[n_new]], T2, 0.05)
  
  moci_impute_mean_t5[i] <- moc_testa(T, T1_mean, alpha_calib_mean[i])
  moci_impute_svd_t5[i] <- moc_testa(T, T1_svd, alpha_calib_svd[i])
  moci_impute_3nn_t5[i] <- moc_testa(T, T1_3nn, alpha_calib_3nn[i])
  moci_impute_6nn_t5[i] <- moc_testa(T, T1_6nn, alpha_calib_6nn[i])
  moci_impute_10nn_t5[i] <- moc_testa(T, T1_10nn, alpha_calib_10nn[i])
  moci_impute_cart_t5[i] <- moc_testa(T, T1_cart, alpha_calib_cart[i])
  moci_impute_pmm_t5[i] <- moc_testa(T, T1_pmm, alpha_calib_pmm[i])
  
}

df_amp <- matrix(rep(0, 3*length(proporcije)), ncol = 3)
df_amp <- data.frame(df_amp)
colnames(df_amp) <- c("Пропорције", "Моћи", "Метод")
df_amp$"Пропорције" <- proporcije
df_amp$"Моћи" <- moci_ampute_t5
df_amp$"Метод" <- "CC"

df_mean <- matrix(rep(0, 3*length(proporcije)), ncol = 3)
df_mean <- data.frame(df_mean)
colnames(df_mean) <- c("Пропорције", "Моћи", "Метод")
df_mean$"Пропорције" <- proporcije
df_mean$"Моћи" <- moci_impute_mean_t5
df_mean$"Метод" <- "Mean"

df_svd <- matrix(rep(0, 3*length(proporcije)), ncol = 3)
df_svd <- data.frame(df_svd)
colnames(df_svd) <- c("Пропорције", "Моћи", "Метод")
df_svd$"Пропорције" <- proporcije
df_svd$"Моћи" <- moci_impute_svd_t5
df_svd$"Метод" <- "SVD"

df_3nn <- matrix(rep(0, 3*length(proporcije)), ncol = 3)
df_3nn <- data.frame(df_3nn)
colnames(df_3nn) <- c("Пропорције", "Моћи", "Метод")
df_3nn$"Пропорције" <- proporcije
df_3nn$"Моћи" <- moci_impute_3nn_t5
df_3nn$"Метод" <- "3NN"

df_6nn <- matrix(rep(0, 3*length(proporcije)), ncol = 3)
df_6nn <- data.frame(df_6nn)
colnames(df_6nn) <- c("Пропорције", "Моћи", "Метод")
df_6nn$"Пропорције" <- proporcije
df_6nn$"Моћи" <- moci_impute_6nn_t5
df_6nn$"Метод" <- "6NN"

df_10nn <- matrix(rep(0, 3*length(proporcije)), ncol = 3)
df_10nn <- data.frame(df_10nn)
colnames(df_10nn) <- c("Пропорције", "Моћи", "Метод")
df_10nn$"Пропорције" <- proporcije
df_10nn$"Моћи" <- moci_impute_10nn_t5
df_10nn$"Метод" <- "10NN"


df_pmm <- matrix(rep(0, 3*length(proporcije)), ncol = 3)
df_pmm <- data.frame(df_pmm)
colnames(df_pmm) <- c("Пропорције", "Моћи", "Метод")
df_pmm$"Пропорције" <- proporcije
df_pmm$"Моћи" <- moci_impute_pmm_t5
df_pmm$"Метод" <- "PMM"

df_cart <- matrix(rep(0, 3*length(proporcije)), ncol = 3)
df_cart <- data.frame(df_cart)
colnames(df_cart) <- c("Пропорције", "Моћи", "Метод")
df_cart$"Пропорције" <- proporcije
df_cart$"Моћи" <- moci_impute_cart_t5
df_cart$"Метод" <- "CART"

library(tikzDevice)

df <- rbind(df_amp, df_mean, df_svd, df_3nn, df_6nn, df_10nn, df_pmm, df_cart)

df_melt <- melt(df, id.vars = c("Метод"), variable.name = "Пропорције", value.name = "Моћи")



options(tz="UTC")
options(tikzDefaultEngine = "xetex")

tikz(file = "NOVO_plot_moci_t5_norm1.tex", width = 5, height = 3, Encoding(UTF-8))

ggplot(data = df, aes(x = df[, 1], y = df[, 2], group = df[, 3], colour = df[, 3]))+
  geom_line()+
  xlab("Пропорције")+
  ylab("Моћи")+
  labs(col="Метод")+
  ggtitle("Стандардна t5 алтернатива")

dev.off()






