# Broj simulacija za Monte karlo
N <- 10000

# Obim uzorka na osnovu kojeg racunamo test statistiku
n <- 50


# kovarijaciona matrica podataka
kovmat <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)

sigma1 <- diag(2)
sigma2 <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
sigme <- list(sigma1, sigma2)
