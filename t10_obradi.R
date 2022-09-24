moci_ampute_t10 <- rep(0, length(proporcije))
moci_impute_mean_t10 <- rep(0, length(proporcije))
moci_impute_svd_t10 <- rep(0, length(proporcije))
moci_impute_3nn_t10 <- rep(0, length(proporcije))
moci_impute_6nn_t10 <- rep(0, length(proporcije))
moci_impute_10nn_t10 <- rep(0, length(proporcije))



for (i in 1:length(proporcije)) 
{
  T1_mean <- rep(0, N)
  T1_svd <- rep(0, N)
  T1_3nn <- rep(0, N)
  T1_6nn <- rep(0, N)
  T1_10nn <- rep(0, N)
  T2 <- rep(0, N)
  for (j in 1:N) 
  {
    uzorak <- rmvt(n,
                   mu = rep(0, 2),
                   sigma = kovmat,
                   v = 10)
    uzorak_izbaceno <- delete_MCAR(uzorak, p = proporcije[i])
    uzorak_amputirano <- na.omit(uzorak_izbaceno)
    T2[j] <- BHEP(uzorak_amputirano, a = 1)
    
    uzorak_ubaceno_mean <- impute_mean(uzorak_izbaceno, type = "columnwise")
    T1_mean[j] <- BHEP(uzorak_ubaceno_mean, a = 1)
    
    uzorak_ubaceno_svd <- eimpute(data.matrix(uzorak_izbaceno), r = 1)$x.imp
    T1_svd[j] <- BHEP(uzorak_ubaceno_svd, a = 1)
    
    uzorak_ubaceno_3nn <- knn.impute(uzorak_izbaceno, k = 3)
    T1_3nn[j] <- BHEP(uzorak_ubaceno_3nn, a = 1)
    
    uzorak_ubaceno_6nn <- knn.impute(uzorak_izbaceno, k = 6)
    T1_6nn[j] <- BHEP(uzorak_ubaceno_6nn, a = 1)
    
    uzorak_ubaceno_10nn <- knn.impute(uzorak_izbaceno, k = 10)
    T1_10nn[j] <- BHEP(uzorak_ubaceno_10nn, a = 1)
  }
  n_new <- nrow(uzorak_amputirano)
  moci_ampute_t10[i] <- moc_testa(nulte1[[n_new]], T2, 0.05)
  
  moci_impute_mean_t10[i] <- moc_testa(T, T1_mean, alpha_calib_mean[i])
  moci_impute_svd_t10[i] <- moc_testa(T, T1_svd, alpha_calib_svd[i])
  moci_impute_3nn_t10[i] <- moc_testa(T, T1_3nn, alpha_calib_3nn[i])
  moci_impute_6nn_t10[i] <- moc_testa(T, T1_6nn, alpha_calib_6nn[i])
  moci_impute_10nn_t10[i] <- moc_testa(T, T1_10nn, alpha_calib_10nn[i])
  
}

df_amp <- matrix(rep(0, 3*length(proporcije)), ncol = 3)
df_amp <- data.frame(df_amp)
colnames(df_amp) <- c("Пропорције", "Моћи", "Метод")
df_amp$"Пропорције" <- proporcije
df_amp$"Моћи" <- moci_ampute_t10
df_amp$"Метод" <- "CC"

df_mean <- matrix(rep(0, 3*length(proporcije)), ncol = 3)
df_mean <- data.frame(df_mean)
colnames(df_mean) <- c("Пропорције", "Моћи", "Метод")
df_mean$"Пропорције" <- proporcije
df_mean$"Моћи" <- moci_impute_mean_t10
df_mean$"Метод" <- "Mean"

df_svd <- matrix(rep(0, 3*length(proporcije)), ncol = 3)
df_svd <- data.frame(df_svd)
colnames(df_svd) <- c("Пропорције", "Моћи", "Метод")
df_svd$"Пропорције" <- proporcije
df_svd$"Моћи" <- moci_impute_svd_t10
df_svd$"Метод" <- "SVD"

df_3nn <- matrix(rep(0, 3*length(proporcije)), ncol = 3)
df_3nn <- data.frame(df_3nn)
colnames(df_3nn) <- c("Пропорције", "Моћи", "Метод")
df_3nn$"Пропорције" <- proporcije
df_3nn$"Моћи" <- moci_impute_3nn_t10
df_3nn$"Метод" <- "3NN"

df_6nn <- matrix(rep(0, 3*length(proporcije)), ncol = 3)
df_6nn <- data.frame(df_6nn)
colnames(df_6nn) <- c("Пропорције", "Моћи", "Метод")
df_6nn$"Пропорције" <- proporcije
df_6nn$"Моћи" <- moci_impute_6nn_t10
df_6nn$"Метод" <- "6NN"

df_10nn <- matrix(rep(0, 3*length(proporcije)), ncol = 3)
df_10nn <- data.frame(df_10nn)
colnames(df_10nn) <- c("Пропорције", "Моћи", "Метод")
df_10nn$"Пропорције" <- proporcije
df_10nn$"Моћи" <- moci_impute_10nn_t10
df_10nn$"Метод" <- "10NN"

library(tikzDevice)

df <- rbind(df_amp, df_mean, df_svd, df_3nn, df_6nn, df_10nn)

df_melt <- melt(df, id.vars = c("Метод"), variable.name = "Пропорције", value.name = "Моћи")



options(tz="UTC")
options(tikzDefaultEngine = "xetex")

tikz(file = "plot_moci_t10_stdn.tex", width = 5, height = 3, Encoding(UTF-8))

ggplot(data = df, aes(x = df[, 1], y = df[, 2], group = df[, 3], colour = df[, 3]))+
  geom_line()+
  xlab("Proporcija")+
  ylab("Moc")+
  labs(col="Метод")+
  ggtitle("Standardna t10 alternativa")

dev.off()
