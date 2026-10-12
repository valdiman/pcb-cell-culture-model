# Code developed to calculate fractions of PCBs in Aroclor 1254
# in 35-mm disc experiments (24 hours) and create figures

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
  den <- 1 + 10^logKlip.w * C.lipid + 10^logKpro.w * C.prot +
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
  
  return(c(f.dis.m, f.dis.c, f.lip, f.prot, f.air, f.cell, f.dis.c.c,
           f.lip.c, f.prot.c))
}

# Calculate fractions for all PCBs ----------------------------------------
num.congener <- length(congener)
result <- NULL
for (i in seq_len(num.congener)) {
  result <- rbind(result, fraction(logKa.w  = logKa.w[i],
                                   logKlip.w = logKlip.w[i],
                                   logKpro.w = logKpro.w[i]))
}

final.result <- data.frame(congener = congener, result)

names(final.result) <- c("congener", "fract.diss.m", "fract.diss.cell",
                         "fract.lip", "fract.prot", "fract.air", "fract.cell",
                         "fract.diss.cell.cell", "fract.lip.cell",
                         "fract.prot.cell")

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
write.csv(final.result, file = "Output/Data/PCB_fractions.csv",
          row.names = FALSE)

# Aroclor 1254 model ------------------------------------------------------
# Read individual PCB composition from Frame
# This is needed to get the fraction of each congener, so then
# the coeluting congeners can be calculated.
A1254 <- read.csv("Data/Aroclor1254_Frame.csv", stringsAsFactors = FALSE)

# Clean Aroclor composition -----------------------------------------------
# Replace blank values and NM with zero
A1254$A1254_fractionV1[
  trimws(A1254$A1254_fractionV1) %in% c("", "NM")] <- "0"

A1254$A1254_fractionV2[
  trimws(A1254$A1254_fractionV2) %in% c("", "NM")] <- "0"

# Convert to numeric
A1254$A1254_fractionV1 <- as.numeric(A1254$A1254_fractionV1)
A1254$A1254_fractionV2 <- as.numeric(A1254$A1254_fractionV2)

# Convert percent composition to fraction
A1254$A1254_fractionV1 <- A1254$A1254_fractionV1 / 100
A1254$A1254_fractionV2 <- A1254$A1254_fractionV2 / 100

# Normalize original Aroclor compositions to sum to 1
A1254$A1254_fractionV1 <- A1254$A1254_fractionV1 / sum(A1254$A1254_fractionV1)

A1254$A1254_fractionV2 <- A1254$A1254_fractionV2 / sum(A1254$A1254_fractionV2)

# Check
sum(A1254$A1254_fractionV1)
sum(A1254$A1254_fractionV2)

# Match Aroclor composition with PCB fraction model -----------------------
idx <- match(A1254$congener, final.result$congener)

# Check that all congeners were matched
A1254$congener[is.na(idx)]

# Combine composition with PCB-specific fraction results
A1254.model <- cbind(A1254, final.result[idx, setdiff(names(final.result),
                                                      c("congener")),
                                         drop = FALSE])

# Analytical reporting exclusions ----------------------------------------
# PCB204 is used as an internal standard.
A1254_Le <- read.csv("Data/Aroclor1254_Lehmler.csv", stringsAsFactors = FALSE)
A1254.groups <- A1254_Le$congener

combine_PCB_group <- function(group_name, data, weight_col) {
  # Example:
  # "PCB4/10" -> c("PCB4", "PCB10")
  pcb_numbers <- strsplit(sub("^PCB", "", group_name), "/")[[1]]
  
  pcb_names <- paste0("PCB", pcb_numbers)
  
  # Select congeners belonging to analytical group
  x <- data[data$congener %in% pcb_names, ]
  
  # Check that all expected congeners were found
  if (nrow(x) != length(pcb_names)) {
    warning(paste("Not all congeners found for group:", group_name))
  }
  
  # Aroclor abundance
  w <- x[[weight_col]]
  w_total <- sum(w)
  
  # Individual congener
  
  if (length(pcb_names) == 1) {
    
    return(
      data.frame(
        congener = group_name,
        A1254_fraction = w_total,
        
        fract.diss.m = x$fract.diss.m,
        fract.diss.cell = x$fract.diss.cell,
        fract.lip = x$fract.lip,
        fract.prot = x$fract.prot,
        fract.air = x$fract.air,
        fract.cell = x$fract.cell,
        fract.diss.cell.cell = x$fract.diss.cell.cell,
        fract.lip.cell = x$fract.lip.cell,
        fract.prot.cell = x$fract.prot.cell))
  }
  
  # Coeluting congeners with nonzero Aroclor abundance
  
  if (w_total > 0) {
    
    return(
      data.frame(
        congener = group_name,
        A1254_fraction = w_total,
        
        fract.diss.m =
          sum(w * x$fract.diss.m) / w_total,
        
        fract.diss.cell =
          sum(w * x$fract.diss.cell) / w_total,
        
        fract.lip =
          sum(w * x$fract.lip) / w_total,
        
        fract.prot =
          sum(w * x$fract.prot) / w_total,
        
        fract.air =
          sum(w * x$fract.air) / w_total,
        
        fract.cell =
          sum(w * x$fract.cell) / w_total,
        
        fract.diss.cell.cell =
          sum(w * x$fract.diss.cell.cell) / w_total,
        
        fract.lip.cell =
          sum(w * x$fract.lip.cell) / w_total,
        
        fract.prot.cell =
          sum(w * x$fract.prot.cell) / w_total))
  }
  
  # Coeluting group with zero abundance
  # Weighted averaging is impossible because all weights are zero.
  # Use unweighted mean of the individual PCB partition fractions.
  
  data.frame(
    congener = group_name,
    A1254_fraction = 0,
    
    fract.diss.m = mean(x$fract.diss.m),
    fract.diss.cell = mean(x$fract.diss.cell),
    fract.lip = mean(x$fract.lip),
    fract.prot = mean(x$fract.prot),
    fract.air = mean(x$fract.air),
    fract.cell = mean(x$fract.cell),
    fract.diss.cell.cell = mean(x$fract.diss.cell.cell),
    fract.lip.cell = mean(x$fract.lip.cell),
    fract.prot.cell = mean(x$fract.prot.cell))
}

A1254.liquid <- do.call(rbind, lapply(A1254.groups,
                                      combine_PCB_group,
                                      data = A1254.model,
                                      weight_col = "A1254_fractionV1"))

A1254.gas <- do.call(rbind, lapply(A1254.groups,
                                   combine_PCB_group,
                                   data = A1254.model,
                                   weight_col = "A1254_fractionV2"))

# Quality-control checks --------------------------------------------------
# 1. Check analytical group names and order
identical(A1254.liquid$congener, A1254_Le$congener)
identical(A1254.gas$congener, A1254_Le$congener)

# 2. Check total Aroclor composition
sum(A1254.liquid$A1254_fraction)
sum(A1254.gas$A1254_fraction)

# 3. Whole-dish mass balance
# medium + cell + air = 1
A1254.liquid$sum.dish <- A1254.liquid$fract.diss.m + A1254.liquid$fract.cell +
  A1254.liquid$fract.air

A1254.gas$sum.dish <- A1254.gas$fract.diss.m + A1254.gas$fract.cell +
  A1254.gas$fract.air

range(A1254.liquid$sum.dish)
range(A1254.gas$sum.dish)

# 4. Within-cell mass balance
# dissolved intracellular + lipid + protein = 1

A1254.liquid$sum.cell <- A1254.liquid$fract.diss.cell.cell +
  A1254.liquid$fract.lip.cell + A1254.liquid$fract.prot.cell

A1254.gas$sum.cell <- A1254.gas$fract.diss.cell.cell +
  A1254.gas$fract.lip.cell + A1254.gas$fract.prot.cell

range(A1254.liquid$sum.cell)
range(A1254.gas$sum.cell)

# 5. Check for missing values
anyNA(A1254.liquid)
anyNA(A1254.gas)

# 6. Check that all fractions are between 0 and 1

fraction_cols <- c(
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

range(A1254.liquid[, fraction_cols])
range(A1254.gas[, fraction_cols])

# Export results
write.csv(A1254.liquid, file = "Output/Data/PCB_fractions_A1254_liquid.csv",
          row.names = FALSE)
write.csv(A1254.gas, file = "Output/Data/PCB_fractions_A1254_gas.csv",
          row.names = FALSE)

# Plots -------------------------------------------------------------------
# Select PCB congeners to plot
selected_PCBs <- c(
  "PCB3",
  "PCB11",
  "PCB25",
  "PCB52",
  "PCB101",
  "PCB153")

# Check that all selected congeners exist
missing_PCBs <- selected_PCBs[!selected_PCBs %in% final.result$congener]

if (length(missing_PCBs) > 0) {
  stop(paste("Congeners not found:", paste(missing_PCBs, collapse = ", ")))
}

# Figure A: fractions in medium, cells and air ----------------------------
# Select congeners and fractions
p.1 <- final.result[
  final.result$congener %in% selected_PCBs,
  c("congener", "fract.diss.m", "fract.cell", "fract.air")]

# Preserve selected PCB order
p.1$congener <- factor(p.1$congener, levels = selected_PCBs)

# Convert to long format
p.1.plot <- melt(p.1, id.var = "congener", variable.name = "phase",
                 value.name = "fraction")

# Set stacking order
p.1.plot$phase <- factor(p.1.plot$phase,
                         levels = c("fract.air", "fract.diss.m",
                                    "fract.cell"))

# Plot
plot.dish <- ggplot(p.1.plot, aes(x = congener, y = fraction, fill = phase)) +
  geom_bar(stat = "identity", col = "white", width = 0.9) +
  scale_fill_manual(labels = c("air", "medium", "cell"),
                    values = c("#fee8c8", "#fdbb84", "#e34a33")) +
  theme_classic() +
  theme(aspect.ratio = 1, text = element_text(size = 14, family = "Helvetica",
                                              face = "bold", color = "black"),
        axis.text.x = element_text(angle = 60, hjust = 1)) +
  xlab("") +
  ylab("Fraction in dish")

plot.dish

# Figure B: fractions within cells ----------------------------------------
p.2 <- final.result[
  final.result$congener %in% selected_PCBs,
  c("congener", "fract.diss.cell.cell", "fract.lip.cell", "fract.prot.cell")]

# Preserve selected PCB order
p.2$congener <- factor(p.2$congener, levels = selected_PCBs)

# Convert to long format
p.2.plot <- melt(p.2, id.var = "congener", variable.name = "phase",
                 value.name = "fraction")

# Set stacking order
p.2.plot$phase <- factor(p.2.plot$phase,
                         levels = c("fract.lip.cell",
                                    "fract.diss.cell.cell",
                                    "fract.prot.cell"))

# Plot
plot.cell <- ggplot(p.2.plot, aes(x = congener, y = fraction, fill = phase)) +
  geom_bar(stat = "identity", col = "white", width = 0.9) +
  scale_fill_manual(labels = c("lipid", "liquid/cytosol", "protein"),
                    values = c("#fee0d2", "#fc9272", "#de2d26")) +
  theme_classic() +
  theme(aspect.ratio = 1, text = element_text(size = 14, family = "Helvetica",
                                              face = "bold", color = "black"),
        axis.text.x = element_text(angle = 60, hjust = 1)) +
  xlab("") +
  ylab("Fraction in cell")

plot.cell








