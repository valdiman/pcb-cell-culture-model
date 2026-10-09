# Code developed to calculate fractions of PCBs in Aroclor 1254
# in well experiments and create figures

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
# Read PCB partition coefficients
pc <- read.csv("Output/Data/PCB_partition_coefficients.csv",
               check.names = FALSE)
# Extract PCB number and sort numerically
pc$PCB_number <- as.numeric(sub("PCB", "", pc$congener))
pc <- pc[order(pc$PCB_number), ]

# Check imported data
head(pc)
str(pc)

# Name parameters ---------------------------------------------------------
congener <- pc$congener
logKa.w  <- pc$logKair_water_37C
logKlip.w <- pc$logK_lipid_water
logKpro.w <- pc$logK_protein_water

fraction <- function(logKa.w, logKlip.w, logKpro.w)
{
  # 35-mm dish experimental conditions
  Vt <- 8.2 / 1000   # Total volume, L
  Vm <- 2.0 / 1000   # Medium volume, L
  Va <- Vt - Vm       # Air/headspace volume, L
  
  # Cells
  Cell <- 1.0e6       # cells/dish
  C.cell <- Cell / Vm # cells/L medium
  
  # Cellular lipid
  C.lipid.cell <- 9.57 * 10^-5 / 10^9  # kg/cell
  C.lipid.L <- C.lipid.cell * C.cell    # kg/L
  dlipids <- 0.905                       # kg/L
  C.lipid <- C.lipid.L / dlipids        # Llipid/Lwater
  
  # Cellular protein
  C.prot.cell.0 <- 1.18 * 10^-4 / 10^9  # kg/cell
  C.prot.L <- C.prot.cell.0 * C.cell     # kg/L
  dprot <- 1.43                           # kg/L
  C.prot <- C.prot.L / dprot              # Lprotein/Lwater
  
  # Intracellular water
  C.water.cell <- 2.84 * 10^-6           # uL/cell
  V.water.cell <- C.water.cell * Cell / 10^6  # L
  
  # Overall denominator
  den <- 1 +
    10^logKlip.w * C.lipid +
    10^logKpro.w * C.prot +
    10^logKa.w * Va / (Vm - V.water.cell)
  
  # Fractions in dish
  f.dis <- 1 / den
  
  f.dis.m <- f.dis * (Vm - V.water.cell) / Vm
  f.dis.c <- f.dis * V.water.cell / Vm
  f.lip <- 10^logKlip.w * C.lipid / den
  f.prot <- 10^logKpro.w * C.prot / den
  f.air <- 10^logKa.w * Va / (Vm - V.water.cell) / den
  
  # Total fraction associated with cells
  f.cell <- f.dis.c + f.lip + f.prot
  
  # Fractions within cells
  f.dis.c.c <- f.dis.c / f.cell
  f.lip.c <- f.lip / f.cell
  f.prot.c <- f.prot / f.cell
  
  return(c(
    f.dis.m,
    f.dis.c,
    f.lip,
    f.prot,
    f.air,
    f.cell,
    f.dis.c.c,
    f.lip.c,
    f.prot.c
  ))
}

# Calculate fractions for all PCBs ----------------------------------------
num.congener <- length(congener)
result <- NULL
for (i in seq_len(num.congener)) {
  
  result <- rbind(
    result,
    fraction(
      logKa.w  = logKa.w[i],
      logKlip.w = logKlip.w[i],
      logKpro.w = logKpro.w[i]
    )
  )
}

final.result <- data.frame(
  congener = congener,
  PCB_number = pc$PCB_number,
  result
)

names(final.result) <- c(
  "congener",
  "PCB_number",
  "fract.diss.m",
  "fract.diss.cell",
  "fract.lip",
  "fract.prot",
  "fract.air",
  "fract.cell",
  "fract.diss.cell.cell",
  "fract.lip.cell",
  "fract.prot.cell"
)

head(final.result)
tail(final.result)

# Check 1: mass-balance check for the whole dish:
final.result$sum.dish <-
  final.result$fract.diss.m +
  final.result$fract.cell +
  final.result$fract.air

range(final.result$sum.dish)

# Check 2: within-cell mass balance:
final.result$sum.cell <-
  final.result$fract.diss.cell.cell +
  final.result$fract.lip.cell +
  final.result$fract.prot.cell

range(final.result$sum.cell)

# Export results
write.csv(final.result, file = "Output/Data/PCB_fractions.csv", row.names = FALSE)

# Aroclor 1254 model ------------------------------------------------------

# Read individual PCB composition
A1254 <- read.csv(
  "Data/Aroclor1254_Frame.csv",
  stringsAsFactors = FALSE
)


# Clean Aroclor composition -----------------------------------------------

# Replace blank values and NM with zero
A1254$A1254_fractionV1[
  trimws(A1254$A1254_fractionV1) %in% c("", "NM")
] <- "0"

A1254$A1254_fractionV2[
  trimws(A1254$A1254_fractionV2) %in% c("", "NM")
] <- "0"

# Convert to numeric
A1254$A1254_fractionV1 <- as.numeric(A1254$A1254_fractionV1)
A1254$A1254_fractionV2 <- as.numeric(A1254$A1254_fractionV2)

# Convert percent composition to fraction
A1254$A1254_fractionV1 <- A1254$A1254_fractionV1 / 100
A1254$A1254_fractionV2 <- A1254$A1254_fractionV2 / 100

# Normalize original Aroclor compositions to sum to 1
A1254$A1254_fractionV1 <-
  A1254$A1254_fractionV1 /
  sum(A1254$A1254_fractionV1)

A1254$A1254_fractionV2 <-
  A1254$A1254_fractionV2 /
  sum(A1254$A1254_fractionV2)

# Check
sum(A1254$A1254_fractionV1)
sum(A1254$A1254_fractionV2)


# Match Aroclor composition with PCB fraction model -----------------------

idx <- match(
  A1254$congener,
  final.result$congener
)

# Check that all congeners were matched
A1254$congener[is.na(idx)]

# Combine composition with PCB-specific fraction results
A1254.model <- cbind(
  A1254,
  final.result[
    idx,
    setdiff(
      names(final.result),
      c("congener", "PCB_number")
    ),
    drop = FALSE
  ]
)


# Analytical reporting exclusions ----------------------------------------

# PCB128 is used as an analytical surrogate.
# PCB166 coelutes with PCB128.
# Neither is included in the reported PCB composition.

exclude_PCBs <- c("PCB128", "PCB166")

# Set excluded congeners to zero BEFORE reportable normalization
A1254.model$A1254_fractionV1_reported <-
  A1254.model$A1254_fractionV1

A1254.model$A1254_fractionV2_reported <-
  A1254.model$A1254_fractionV2

A1254.model$A1254_fractionV1_reported[
  A1254.model$congener %in% exclude_PCBs
] <- 0

A1254.model$A1254_fractionV2_reported[
  A1254.model$congener %in% exclude_PCBs
] <- 0

# Normalize analytically reportable composition to sum to 1
A1254.model$A1254_fractionV1_reported <-
  A1254.model$A1254_fractionV1_reported /
  sum(A1254.model$A1254_fractionV1_reported)

A1254.model$A1254_fractionV2_reported <-
  A1254.model$A1254_fractionV2_reported /
  sum(A1254.model$A1254_fractionV2_reported)

# Check
sum(A1254.model$A1254_fractionV1_reported)
sum(A1254.model$A1254_fractionV2_reported)

# Calculate Aroclor 1254 distribution ------------------------------------

# Version 1
A1254.model$V1_diss_m <-
  A1254.model$A1254_fractionV1_reported *
  A1254.model$fract.diss.m

A1254.model$V1_diss_cell <-
  A1254.model$A1254_fractionV1_reported *
  A1254.model$fract.diss.cell

A1254.model$V1_lipid <-
  A1254.model$A1254_fractionV1_reported *
  A1254.model$fract.lip

A1254.model$V1_protein <-
  A1254.model$A1254_fractionV1_reported *
  A1254.model$fract.prot

A1254.model$V1_air <-
  A1254.model$A1254_fractionV1_reported *
  A1254.model$fract.air

A1254.model$V1_cell <-
  A1254.model$A1254_fractionV1_reported *
  A1254.model$fract.cell


# Version 2
A1254.model$V2_diss_m <-
  A1254.model$A1254_fractionV2_reported *
  A1254.model$fract.diss.m

A1254.model$V2_diss_cell <-
  A1254.model$A1254_fractionV2_reported *
  A1254.model$fract.diss.cell

A1254.model$V2_lipid <-
  A1254.model$A1254_fractionV2_reported *
  A1254.model$fract.lip

A1254.model$V2_protein <-
  A1254.model$A1254_fractionV2_reported *
  A1254.model$fract.prot

A1254.model$V2_air <-
  A1254.model$A1254_fractionV2_reported *
  A1254.model$fract.air

A1254.model$V2_cell <-
  A1254.model$A1254_fractionV2_reported *
  A1254.model$fract.cell

# Analytical PCB grouping -------------------------------------------------

A1254.groups <- read.csv(
  "Data/Aroclor1254_composition.csv",
  stringsAsFactors = FALSE
)

# Function to combine coeluting congeners ---------------------------------

combine_PCB_group <- function(group_name, data) {
  
  # Example:
  # "PCB18/30" -> c("PCB18", "PCB30")
  pcb_numbers <- strsplit(
    sub("^PCB", "", group_name),
    "/"
  )[[1]]
  
  pcb_names <- paste0("PCB", pcb_numbers)
  
  # Select individual congeners in analytical group
  x <- data[data$congener %in% pcb_names, ]
  
  data.frame(
    congener = group_name,
    
    # Reportable Aroclor composition
    A1254_fractionV1 =
      sum(x$A1254_fractionV1_reported),
    
    A1254_fractionV2 =
      sum(x$A1254_fractionV2_reported),
    
    # Version 1
    V1_diss_m    = sum(x$V1_diss_m),
    V1_diss_cell = sum(x$V1_diss_cell),
    V1_lipid     = sum(x$V1_lipid),
    V1_protein   = sum(x$V1_protein),
    V1_air       = sum(x$V1_air),
    V1_cell      = sum(x$V1_cell),
    
    # Version 2
    V2_diss_m    = sum(x$V2_diss_m),
    V2_diss_cell = sum(x$V2_diss_cell),
    V2_lipid     = sum(x$V2_lipid),
    V2_protein   = sum(x$V2_protein),
    V2_air       = sum(x$V2_air),
    V2_cell      = sum(x$V2_cell)
  )
}

# Combine congeners according to analytical method ------------------------

A1254.final <- do.call(
  rbind,
  lapply(
    A1254.groups$congener,
    combine_PCB_group,
    data = A1254.model
  )
)

# Mass-balance checks -----------------------------------------------------

A1254.final$V1_sum <-
  A1254.final$V1_diss_m +
  A1254.final$V1_cell +
  A1254.final$V1_air

A1254.final$V2_sum <-
  A1254.final$V2_diss_m +
  A1254.final$V2_cell +
  A1254.final$V2_air


# Per analytical group
range(
  A1254.final$V1_sum -
    A1254.final$A1254_fractionV1
)

range(
  A1254.final$V2_sum -
    A1254.final$A1254_fractionV2
)


# Entire Aroclor composition
sum(A1254.final$A1254_fractionV1)
sum(A1254.final$A1254_fractionV2)

# Entire modeled Aroclor distribution
sum(A1254.final$V1_sum)
sum(A1254.final$V2_sum)

# Congener profile within cells -------------------------------------------

A1254.final$V1_cell_profile <-
  A1254.final$V1_cell /
  sum(A1254.final$V1_cell)

A1254.final$V2_cell_profile <-
  A1254.final$V2_cell /
  sum(A1254.final$V2_cell)

# Check
sum(A1254.final$V1_cell_profile)
sum(A1254.final$V2_cell_profile)

exclude_PCBs <- c("PCB128", "PCB166")

A1254.model$A1254_fractionV1_reported <-
  A1254.model$A1254_fractionV1

A1254.model$A1254_fractionV2_reported <-
  A1254.model$A1254_fractionV2

# Exclude PCB128 and PCB166
A1254.model$A1254_fractionV1_reported[
  A1254.model$congener %in% exclude_PCBs
] <- 0

A1254.model$A1254_fractionV2_reported[
  A1254.model$congener %in% exclude_PCBs
] <- 0

# Renormalize reportable composition
A1254.model$A1254_fractionV1_reported <-
  A1254.model$A1254_fractionV1_reported /
  sum(A1254.model$A1254_fractionV1_reported)

A1254.model$A1254_fractionV2_reported <-
  A1254.model$A1254_fractionV2_reported /
  sum(A1254.model$A1254_fractionV2_reported)

sum(A1254.model$A1254_fractionV1_reported)
sum(A1254.model$A1254_fractionV2_reported)

sum(A1254.final$A1254_fractionV1)
sum(A1254.final$A1254_fractionV2)

sum(A1254.final$V1_sum)
sum(A1254.final$V2_sum)

A1254.model[
  A1254.model$congener %in% c("PCB128", "PCB166"),
  c(
    "congener",
    "A1254_fractionV1",
    "A1254_fractionV1_reported",
    "A1254_fractionV2",
    "A1254_fractionV2_reported"
  )
]

sum(
  A1254.model$A1254_fractionV1_reported[
    !A1254.model$congener %in% c("PCB128", "PCB166")
  ]
)

sum(
  A1254.model$A1254_fractionV2_reported[
    !A1254.model$congener %in% c("PCB128", "PCB166")
  ]
)

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
