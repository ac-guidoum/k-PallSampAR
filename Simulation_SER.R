##################################################
rho_0     <- 0.3
##################################################
# Example 1 — n_i equals 

# Table 1 (n_i = 50) — k = 4

n_vec     <- rep(50, 4)
deltat_k4 <- c(0, 0.2, 0.6, 1, 2, 3, 8, 9)            # k = 4
sre_table <- PTE_KSAMPLE_parallel(n_vec = n_vec, rho_1 = rho_0,
                                   deltat_grid = deltat_k4,
                                   M = 10000, ncores = 12, p = 0.5, seed = 123)
show4(sre_table) 
plot_sre(sre_table, n_vec, rho_0)

# Table 1 (n_i = 50) — k = 9
n_vec     <- rep(50, 9)
deltat_k9 <- c(0, 0.2, 0.6, 1, 5, 6, 7, 12, 13, 14)   # k = 9
sre_table <- PTE_KSAMPLE_parallel(n_vec = n_vec, rho_1 = rho_0,
                                   deltat_grid = deltat_k9,
                                   M = 10000, ncores = 12, p = 0.5, seed = 123)
show4(sre_table)
plot_sre(sre_table, n_vec, rho_0)

# Table 2 (n_i = 200) — k = 4
n_vec     <- rep(200, 4)
deltat_k4 <- c(0, 0.2, 0.6, 1, 2, 3, 6, 7, 8, 14, 16)            # k = 4
sre_table <- PTE_KSAMPLE_parallel(n_vec = n_vec, rho_1 = rho_0,
                                   deltat_grid = deltat_k4,
                                   M = 10000, ncores = 12, p = 0.5, seed = 123)
show4(sre_table)
plot_sre(sre_table, n_vec, rho_0)

# Table 2 (n_i = 200) — k = 9
n_vec     <- rep(200, 9)
deltat_k9 <- c(0, 0.2, 0.6, 1, 5, 6, 7, 12, 13, 14, 20, 22, 24)   # k = 9
sre_table <- PTE_KSAMPLE_parallel(n_vec = n_vec, rho_1 = rho_0,
                                   deltat_grid = deltat_k9,
                                   M = 10000, ncores = 12, p = 0.5, seed = 123)
show4(sre_table)
plot_sre(sre_table, n_vec, rho_0)

##################################################
# Example 2 — n_i DIFFERENT 

# Table 3 k = 4

n_vec     <- c(30, 50, 100, 300)
deltat_k4 <- c(0, 0.2, 0.6, 1, 2, 3, 8, 9, 14, 15)            # k = 4
sre_table <- PTE_KSAMPLE_parallel(n_vec = n_vec,rho_1 = rho_0,
                                  deltat_grid = deltat_k4,M = 10000, 
                                  ncores = 12, p = 0.5, seed = 123)
show4(sre_table)
plot_sre(sre_table, n_vec, rho_0)

# Table 3 k = 9

n_vec     <- c(30, 50, 80, 100, 120, 180, 220, 300, 330)
deltat_k9 <- c(0, 0.2, 0.6, 1, 2, 3, 8, 9, 14, 15, 24, 25, 26)       # k = 9
sre_table <- PTE_KSAMPLE_parallel(n_vec = n_vec,rho_1 = rho_0,
                                  deltat_grid = deltat_k9,M = 10000, 
                                  ncores = 12, p = 0.5, seed = 123)
show4(sre_table)
plot_sre(sre_table, n_vec, rho_0)
