qcondiglcop <- function(tau, u, cpar) {
    1 - cnstr_Psi((1 - u)/qgamma(tau, cpar-1), cpar)
}

pcondiglcop <- function(v, u, cpar) {
    pgamma((1 - u) / cnstr_Psiinv(1 - v, cpar), cpar - 1)
}

qcondiglcop21 <- qcondiglcop
pcondiglcop21 <- pcondiglcop

pcondiglcop12 <- function(u, v, cpar) {
    pkinv <- cnstr_Psiinv(1 - v, cpar)
    1 - pgamma((1 - u) / pkinv, cpar) / pgamma(1 / pkinv, cpar)
}

qcondiglcop12 <- function(tau, v, cpar) {
    pkinv <- cnstr_Psiinv(1 - v, cpar)
    1 - pkinv * qgamma((1 - tau) * pgamma(1/pkinv, cpar), cpar)
}

diglcop <- function(u, v, cpar) {
    pkinv <- cnstr_Psiinv(1 - v, cpar)
    (1-u)^(cpar-1) / pkinv^k * exp(-(1-u)/pkinv) / (gamma(k) * pgamma(1/pkinv, cpar))
}
