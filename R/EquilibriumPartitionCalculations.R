# Code developed to calculate partition coefficients of selected PCBs
# from Aroclor 1254 in well experiments and create figures

# No packages need

# Read data ---------------------------------------------------------------
# Ulrich, N., Endo, S., Brown, T.N., Watanabe, N., Bronner, G., Abraham, M.H.,
# Goss, K.-U., UFZ-LSER database v 3.2.1 [Internet], Leipzig, Germany,
# Helmholtz Centre for Environmental Research-UFZ. 2017 [accessed on 16.03.2023].
# Available from http://www.ufz.de/lserd
pcb <- read.csv("Data/PCB_Abraham_descriptors.csv")
# Read Henry's law data
henry <- read.csv("Data/PCB_Henry_constant.csv")

# Calculate partition coefficients ----------------------------------------
calc_logK <- function(E, S, A, B, V, L,
                      e = 0, s = 0, a = 0, b = 0,
                      v = 0, l = 0, c = 0) {
  e * E + s * S + a * A + b * B + v * V + l * L + c
}


# Geisler, A.; Endo, S.; Goss, K.-U. Environ. Sci. Technol. 2012, 46, 9519-9524.
pcb$logK_lipid_water <- calc_logK(E = pcb$E, S = pcb$S, A = pcb$A,
                                  B = pcb$B, V = pcb$V, L = pcb$L,
                                  e =  0.70, s = -1.08, a = -1.72,
                                  b = -4.14, v =  4.11, l =  0.00,
                                  c = -0.07)

# Endo, S.; Bauerfeind, J.; Goss, K.-U. Environ. Sci. Technol. 2012, 46 (22), 12697-12703.
pcb$logK_protein_water <- calc_logK(E = pcb$E, S = pcb$S, A = pcb$A,
                                    B = pcb$B, V = pcb$V, L = pcb$L,
                                    e =  0.51, s = -0.51, a =  0.26,
                                    b = -2.98, v =  3.01, l =  0.00,
                                    c = -0.65)

# Temperatures
T_ref <- 298.15   # K (25 C)
T     <- 310.15   # K (37 C)

# Gas constant
R <- 8.314        # J/(mol K)

# Convert log10(Kair/water) to Kair/water
henry$Kair_water_ref <- 10^henry$logKair.water

# Temperature correction
henry$Kair_water_37C <- henry$Kair_water_ref * exp(
  -(henry$dUaw / R) * (1 / T - 1 / T_ref))

# Convert corrected K back to log10
henry$logKair_water_37C <- log10(henry$Kair_water_37C)

# Create final dataset ----------------------------------------------------
# Make sure the congener column has the same name in both data frames
names(pcb)[names(pcb) == "Congener"] <- "congener"

# Merge partition coefficients
pcb_partition <- merge(pcb[, c("congener", "logK_lipid_water",
                               "logK_protein_water")],
                       henry[, c("congener", "logKair_water_37C")],
                       by = "congener", all = TRUE)

# Save dataset
write.csv(pcb_partition, "Data/PCB_partition_coefficients.csv",
          row.names = FALSE)



