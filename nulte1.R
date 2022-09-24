########################################################################
# R nema ugradjenu raspodelu BHEP test statistike, tako da cemo za svaki obim
# uzorka simulirati tu raspodelu da bismo mogli sa njome da poredimo

# Svi moguci obimi uzorka
obimi <- 1:n

# Lista u koju cemo upisivati test statistike
nulte1 <- list()

# Upisujemo ih u petlji
for (obim in obimi)
{
  T <- rep(0, N)
  for (i in 1:N)
  {
    uzorak <- rmvnorm(obim,
                      mu = rep(0, 2),
                      sigma = kovmat)
    T[i] <- BHEP(uzorak, a = 1)
  }
  nulte1[[obim]] <- T
}

# Sminka, zbog manjeg broja vrsta od kolona
nulte1[[1]] <- NA
nulte1[[2]] <- NA


##############################################################################
# Cesto cemo koristiti T za pun obim n, pa cemo ga posebno izdvojiti
T <- nulte1[[n]]

# Trebace nam i funkcija koja racuna moc testa
# Levi argument: ono sto je tacno
# Desni argument: Ono sto nas zanima
moc_testa <- function(t1, t2, alpha)
{
  # Verovatnoca da je t2 u kriticnoj oblasti koja je formirana na osnovu t1
  1 - ecdf(t2)(quantile(t1, probs = 1 - alpha))
}
