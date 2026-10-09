# Code developed to calculate fractions of selected PCBs and
# their metabolites (i.e., OH-PCBs and OH-PCB sulfates) in
# well experiments and create figures (Figure 7 A and B)

# Install packages and load libraries -------------------------------------
# Packages
{
  install.packages("gridExtra")
  install.packages("ggplot2")
  install.packages("extrafont")
}

# Libraries
{
  library(reshape2)
  library(ggplot2)
  library(extrafont) 
}

# Physical-chemical Properties --------------------------------------------
# Create matrix to storage physical-chemical properties
pc <- matrix(NA, nrow = 12, ncol = 5)
# Add column names
colnames(pc) <- c('congener', 'logKlipid/water', 'logKprotein/water',
                  'logKair/water', 'dUaw')
# Add row names
pc[, 1] <- c("PCB3", "4-OH-PCB3", "4-OH-PCB3 sulfate", "PCB11",
             "4-OH-PCB11", "4-OH-PCB11 sulfate", "PCB25",
             "4-OH-PCB25", "4-OH-PCB25 sulfate", "PCB52",
             "4-OH-PCB52", "4-OH-PCB52 sulfate")

# Add physical-chemical properties
# SMILES was used to determine OH-PCBs and OH-PCB sulfates
# Reference:
# Ulrich, N., Endo, S., Brown, T.N., Watanabe, N., Bronner, G., Abraham, M.H.,
# Goss, K.-U., UFZ-LSER database v 3.2.1 [Internet], Leipzig, Germany,
# Helmholtz Centre for Environmental Research-UFZ. 2017 [accessed on 16.03.2023].
# Available from http://www.ufz.de/lserd

# logKlipid/water
pc[, 2] <- c(5.05, 3.19, -0.33, 5.60, 4.85, 0.35, 6.06, 5.39, 0.88,
             6.50, 5.22, 1.44)
# logKprotein/water
pc[, 3] <- c(3.40, 3.02, 0.70, 3.84, 3.53, 1.21, 4.22, 4.36, 1.65,
             4.58, 4.34, 2.09)
# logKair/water
pc[, 4] <- c(-1.95, -6.13, -13.42, -1.93, -4.57, -13.48, -1.89,
             -4.64, -13.55, -1.88, -5.59, -13.39)
# dUaw (no data available for OH-PCBs and OH-PCB sulfates)
pc[, 5] <- c(48734.74, 0.00, 0.00, 51662.48, 0.00, 0.00, 53590.22,
             0.00, 0.00, 55517.96, 0.00, 0.00)
# Transform pc matrix into data.frame
pc <- as.data.frame(pc)
# Transform values into numeric values
pc$`logKlipid/water` <- as.numeric(pc$`logKlipid/water`)
pc$`logKprotein/water` <- as.numeric(pc$`logKprotein/water`)
pc$`logKair/water` <- as.numeric(pc$`logKair/water`)
pc$dUaw <- as.numeric(pc$dUaw)

# Calculations ------------------------------------------------------------
# Name parameters
{
  congener <- pc$congener
  logKa.w <- pc$`logKair/water`
  logKlip.w <- pc$`logKlipid/water`
  logKpro.w <- pc$`logKprotein/water`
  dUaw <- pc$dUaw
}

# Function to calculate fractions
fraction = function(logKa.w, dUaw, logKlip.w, logKpro.w)
{
  # Multiwell volumes
  # https://www.corning.com/catalog/cls/documents/drawings/MD_Microplate_Dimension_Sheets_Multiple_Well.pdf
  Vt <- 3.47/1000 # Total volume L
  Vm <- 0.5/1000 # Medium volume L
  Va <- Vt-Vm # Air volume L
  
  # Concentrations in the wells
  Cell <- 4*10^4 # # cell/well
  C.cell <- Cell/Vm # Concentration of cell per well, cell/Vm (L)
  C.lipid.cell <- 9.57*10^-5/10^9 # Lipid content kg/cell
  C.lipid.L <- C.lipid.cell*C.cell # kg/L
  dlipids <- 0.905 # kg/L Reference: https://pubmed.ncbi.nlm.nih.gov/8148928/
  C.lipid <- C.lipid.L/dlipids # Lipid concentration Llipid/lwater
  C.prot.cell.0 <- 1.18*10^-4/10^9 # Protein content kg/#cell
  C.prot.L <- C.prot.cell.0*C.cell # kg/L
  dprot <- 1.43 # kg/L Reference: https://pubmed.ncbi.nlm.nih.gov/10930825/
  C.prot <- C.prot.L/dprot # Lipid concentration Lprot/lwater
  C.water.cell <- 2.84*10^-6 # uL/cell
  V.water.cell <- C.water.cell*Cell/10^6 # L water inside cell
  
  # Temperature correction for Kaw for just PCBs
  R <- 8.3144 # J/(mol K)
  Tst <- 25 # Standard air temperature C
  Tst.1 <- 273.15 + Tst # K
  T2 <- 37 # Incubation air temperature C
  T2 <- 273.15 + T2
  Ka.w.t <- 10^(logKa.w)*exp(-dUaw/R*(1/T2-1/Tst.1)) # Ka.w corrected by incubation air temperature
  logKa.w.t <- log10(Ka.w.t)
  
  # Fraction calculations
  den <- 1 + 10^(logKlip.w)*C.lipid + 10^(logKpro.w)*C.prot +
    10^(logKa.w.t)*Va/(Vm-V.water.cell)
  # Dissolved fraction
  f.dis <- 1/den
  # Dissolved fraction in the medium
  f.dis.m <- f.dis*(Vm-V.water.cell)/Vm
  # Dissolved fraction in the cell
  f.dis.c <- f.dis*V.water.cell/Vm
  # Lipid fraction
  f.lip <- 10^(logKlip.w)*C.lipid/den
  # Protein fraction
  f.prot <- 10^(logKpro.w)*C.prot/den
  # Air (head space) fraction
  f.air <- 10^(logKa.w.t)*Va/Vm/den
  # Sum of cell fraction (dissolved + lipid + protein)
  f.cell <- f.dis.c + f.lip + f.prot
  # Fraction for dissolved in cell
  f.dis.c.c <- f.dis.c/f.cell
  # Fraction of lipids in the cell
  f.lip.c <- f.lip/f.cell
  # Fraction of protein in the cell
  f.prot.c <- f.prot/f.cell
  # Combine all in frac
  frac <- c(f.dis.m, f.dis.c, f.lip, f.prot, f.air, f.cell, f.dis.c.c,
            f.lip.c, f.prot.c)
}

num.congener <- length(congener)
result <- NULL
for (i in 1:num.congener)
{
  result <- rbind(result, fraction(logKa.w[i], dUaw[i],
                                   logKlip.w[i], logKpro.w[i]))
}

final.result <- data.frame(congener, result)
# Name the columns
names(final.result) <- c("congener", "fract.diss.m", "fract.diss.cell",
                         "fract.lip", "fract.prot", "fract.air",
                         "fract.cell", "fract.diss.cell.cell",
                         "fract.lip.cell", "fract.prot.cell")
# Export results
write.csv(final.result, file = "Output/Data/csv/Fractions.csv")

# Plots -------------------------------------------------------------------
# Figure 7 (A). Just fractions in cell, medium and air
# Create data.frame with needed fractions
p.1 <- final.result[,!names(final.result) %in% c("fract.diss.cell",
                              "fract.lip", "fract.prot",
                              "fract.diss.cell.cell",
                              "fract.lip.cell", "fract.prot.cell")]
# Transform data.frame p.1 to 3 column data.frame
p.1.plot <- melt(p.1, id.var = c("congener"),
                variable.name = "phase", value.name = "fraction")
# Name the compounds
p.1.plot$congener <- factor(p.1.plot$congener,
                           levels = c('PCB3', '4-OH-PCB3', '4-OH-PCB3 sulfate',
                                      'PCB11', '4-OH-PCB11', '4-OH-PCB11 sulfate',
                                      'PCB25', '4-OH-PCB25', '4-OH-PCB25 sulfate',
                                      'PCB52', '4-OH-PCB52', '4-OH-PCB52 sulfate'))
# Organize fraction to be displayed in plot
p.1.plot$phase <- factor(p.1.plot$phase,
                        levels = c('fract.air', 'fract.diss.m',
                                   'fract.cell'))
# Plot (Figure 7 (A))
ggplot(p.1.plot, aes(x = congener, y = fraction, fill = phase)) + 
  geom_bar(stat = "identity", col = "white", width = 0.9) +
  scale_fill_manual(labels = c("air" , "medium", "cell"),
                    values = c("#fee8c8", "#fdbb84", "#e34a33")) +
  theme_classic() +
  theme(aspect.ratio = 10/10,
        text = element_text(size = 14, family = "Helvetica",
                            face = "bold", color = "black")) +
  xlab(expression(bold(""))) +
  ylab(expression("Fraction in well")) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

# Figure 7 (B). Just fractions inside cell, i.e., protein, lipids and liquid/cytosol
# Create data.frame with needed fractions
p.2 <- final.result[,!names(final.result) %in% c("fract.diss.m", "fract.diss.cell",
                                                 "fract.lip", "fract.prot",
                                                 "fract.air", "fract.cell")]
# Transform data.frame p.1 to 3 column data.frame
p.2.plot <- melt(p.2, id.var = c("congener"),
                 variable.name = "phase", value.name = "fraction")
# Name the compounds
p.2.plot$congener <- factor(p.2.plot$congener,
                            levels = c('PCB3', '4-OH-PCB3', '4-OH-PCB3 sulfate',
                                       'PCB11', '4-OH-PCB11', '4-OH-PCB11 sulfate',
                                       'PCB25', '4-OH-PCB25', '4-OH-PCB25 sulfate',
                                       'PCB52', '4-OH-PCB52', '4-OH-PCB52 sulfate'))
# Organize fraction to be displayed in plot
p.2.plot$phase <- factor(p.2.plot$phase,
                         levels = c('fract.lip.cell', 'fract.diss.cell.cell',
                                    'fract.prot.cell'))
# Plot (Figure 7 (B))
ggplot(p.2.plot, aes(x = congener, y = fraction, fill = phase)) + 
  geom_bar(stat = "identity", col = "white", width = 0.9) +
  scale_fill_manual(labels = c("lipid" , "liquid/cytosol", "protein"),
                    values = c("#fee0d2", "#fc9272", "#de2d26")) +
  theme_classic() +
  theme(aspect.ratio = 10/10,
        text = element_text(size = 14, family = "Helvetica",
                            face = "bold", color = "black")) +
  xlab(expression(bold(""))) +
  ylab(expression("Fraction in cell")) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
