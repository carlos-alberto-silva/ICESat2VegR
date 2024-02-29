calculatePAI <- function(vegetationPhotons, groundPhotons, rho_v_rho_g_factor = 1) {
  nv <- sum(vegetationPhotons)
  ng <- sum(groundPhotons)
  cover <- 1 / (1 + rho_v_rho_g_factor * (ng / nv))

  P0 <- 1 - 1 * cover
  PAI <- -log(P0) * 2
  PAVD <- Pai / length(vegetationPhotons)

  Rv <- rev(cumsum(vegetationPhotons))
  Pz <- 1 - (Rv / nv) * cover
  PAVDz <- -log(Pz[-length(Pz)] / (Pz[-1])) * 2

  list(
    PAI = PAI,
    PAVD = PAVD,
    PAVDz = PAVDz
  )

}

exampleVegPhotons <- c(5, 3, 1, 2, 4)
exampleGroundPhotons <- c(4)
rho_v_rho_g_factor <- 1.187
(lai <- calculatePAI(exampleVegPhotons, exampleGroundPhotons, rho_v_rho_g_factor))
