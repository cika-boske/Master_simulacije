# KALIBRACIJA PMM

# Posmatracemo nedostajanja izmedju 3% i 15% sa razmakom od po pola procenta
proporcije <- seq(0.03, 0.15, 0.005)

# Ideja je da smanjujemo nivo znacajnosti sve dok ne dobijemo dobar size testa
# tj. da namerno testiramo s pogresnim nivoom znacajnosti da bismo korigovali
# Trazimo takvo alfa za svaku proporciju
alpha_calib_pmm <- rep(0, length(proporcije))

# Koristicemo klasican grid-search: isprobavacemo razne alfe
alpha_probno <- seq(0.0001, 0.05, 0.0001)
# Dodacemo jos jedan mali broj na pocetak
privr <- c(0.000001)
alpha_probno <- c(privr, alpha_probno)

# U ovaj niz belezicemo vracenu moc za svaku proporciju
vracene_moci_pmm <- rep(0, length(alpha_calib_pmm))

# Ulazimo u petlju za svaku proporciju, brisemo, imputiramo, racunamo
for (i in 1:length(proporcije))
{
  # u petlji nam je fiksno proporcije[i]
  T_ubaceno <- rep(0, N)
  for (j in 1:N)
  {
    # Sada generisemo uzorak svaki put, brisemo pa imputiramo
    uzorak <- rmvnorm(n,
                      mu = rep(0, 2),
                      sigma = kovmat)
    # vrsimo izbacivanje
    uzorak_izbaceno <- delete_MCAR(uzorak, p = proporcije[i])
    # skaliramo
    mins <- c(min(uzorak_izbaceno[, 1], na.rm = TRUE), min(uzorak_izbaceno[, 2], na.rm = TRUE))
    maxs <- c(max(uzorak_izbaceno[, 1], na.rm = TRUE), max(uzorak_izbaceno[, 2], na.rm = TRUE))
    uzorak_skalirano <- uzorak_izbaceno
    uzorak_skalirano[, 1] <- (uzorak_izbaceno[, 1] - mins[1]) / (maxs[1] - mins[1])
    uzorak_skalirano[, 2] <- (uzorak_izbaceno[, 2] - mins[2]) / (maxs[2] - mins[2])
    # vrsimo imputaciju
    uzorak_ubaceno <- data.matrix(complete(mice(uzorak_skalirano, m = 1)))
    # reskaliramo
    uzorak_ubaceno[, 1] <- (maxs[1] - mins[1]) * uzorak_ubaceno[, 1] + mins[1]
    uzorak_ubaceno[, 2] < (maxs[2] - mins[2]) * uzorak_ubaceno[, 2] + mins[2]
    uzorak_ubaceno[, 1] <- (uzorak_ubaceno[, 1] - min(uzorak_ubaceno[, 1])) / (max(uzorak_ubaceno[, 1]) - min(uzorak_ubaceno[, 1]))
    uzorak_ubaceno[, 2] <- (uzorak_ubaceno[, 2] - min(uzorak_ubaceno[, 2])) / (max(uzorak_ubaceno[, 2]) - min(uzorak_ubaceno[, 2]))
    # Na osnovu ovako napravljenog imputiranog uzorka,
    # racunamo BHEP statistiku
    T_ubaceno[j] <- BHEP(uzorak_ubaceno, a = 1)
  }
  # Sada smo izasli iz unutrasnje petlje i imamo generisanu raspodelu test
  # statistike nakon imputacije, za datu proporciju
  # Sada treba naci sa kojim alfa ona pravi gresku prve vrste od 0.05
  # u naredni vektor upisujemo povratne moci
  povratne_moci <- rep(0, length(alpha_probno))
  # Sada racunamo povratnu moc za svako probno alfa
  for (l in 1:length(povratne_moci)) {
    povratne_moci[l] <- moc_testa(T, T_ubaceno, alpha_probno[l])
  }
  # Sada imamo povratne moci za sve probne alfa, treba naci najmanju
  # Pravimo apsolutne razlike sa 0.05
  tmp <- rep(0.05, length(povratne_moci))
  abs_razlike <- abs(tmp - povratne_moci)
  indeksi_min <- which.min(abs_razlike)
  # Konacno, u kalibrisano alfa upisujemo kandidata na tom indeksu
  alpha_calib_pmm[i] <- alpha_probno[indeksi_min]
  # Zbog kontrole da je zaista vraceno 0.05
  vracene_moci_pmm[i] <- moc_testa(T, T_ubaceno, alpha_calib_pmm[i])
}