##############################
# Preparing data from DPMUnc #
##############################

# Date: 2022/03/07
# Guillermo Reales

# This script will prepare the data for DPMUnc clustering.

# In this case, we excluded FinnGen and include PanUKBB (R2) as the major panel set. Then we'll select significant traits at FDR 1%.
# Next we'll keep IMDs and BMK, removing anything obviously redundant, giving preference to European datasets and our own/independent IMD.
# Then, we'll output IMD projections for clustering, keeping the list of BMK for later.


# Load packages
library(data.table)
library(magrittr)



# Import projection table after QC, processing and FDR.
pt.sig <- fread("../../data/PT_sig_202303_UKBB.tsv")

# Not interested in blood cells for clustering
pt.cl <- pt.sig[ Trait_class != "BC"]

# We're interested on PC2, 5, and 8, so we can keep just them
pt.cl <- pt.cl[ PC %in% paste0("PC", c(2,5,8))]
sigany <- unique(pt.cl[FDR.PC < 0.01, Trait]) # Get traits that are significant for any, and keep only those
pt.cl  <- pt.cl[Trait %in% sigany]

pt.st  <- reshape(pt.cl[,.(Trait, PC, stars)], idvar="Trait", timevar = "PC", direction = "wide") # Get significant info, if problematic

# Get metadata to help decide
m <- fread("../../data/Metadata_20230307-v1.tsv")

pt.st  <- merge(pt.st, m[, .(Trait, Trait_long, Trait_class, N0, N1, N, Year, Population)], by="Trait")

# Save to do an easier selection of diseases to exclude manually
fwrite(pt.st, "../../data/tmp_redundancies_202303_UKBBsig.tsv", sep="\t")

# Retrieve non-redundant list
nr <- fread("../../data/tmp_processed_redundancies_202303_UKBBsig.tsv")  %>% .[is.na(Redundant), Trait]  

pt.cl  <- pt.cl[ Trait %in% nr ] # Keep non-redundant only

# Now we have the basic pt.sig, with updated delta and var.delta.
# We can now select the PCs we're interested in

# 1. PC5 & PC8
# We want all significant IMD at FDR.PC 1% (stars == "●") for at least one of the selected PCs.
sig <- pt.cl[ PC %in% c("PC5", "PC8") & FDR.PC < 0.01 & Trait_class == "IMD", unique(Trait)]
pt58  <- pt.cl[ PC %in% c("PC5", "PC8") & Trait %in% sig]

length(unique(pt58$Trait))
# 22 IMD datasets after filtering
length(unique(pt58$Label))
# 22 unique labels, too!

# Get data in shape
delta58 <- reshape(pt58[  PC %in% c("PC5", "PC8") , .(PC, Delta, Trait)], idvar="Trait", timevar = "PC", direction = "wide")
delta_var58 <-  reshape(pt58[  PC %in% c("PC5", "PC8") , .(PC, Var.Delta, Trait)], idvar="Trait", timevar = "PC", direction = "wide")

# Save data
fwrite(delta58,     "../data/UKBBsig_PC58_IMD_Delta.tsv", sep="\t")
fwrite(delta_var58, "../data/UKBBsig_PC58_IMD_Var.tsv", sep="\t")

# 2. PC2

pt2 <- pt.cl[ PC == "PC2" & FDR.PC < 0.01 & Trait_class == "IMD"]

# Get data in shape
delta2 <- reshape(pt2[, .(PC, Delta, Trait)], idvar="Trait", timevar = "PC", direction = "wide")
delta_var2 <-  reshape(pt2[, .(PC, Var.Delta, Trait)], idvar="Trait", timevar = "PC", direction = "wide")

# Save data
fwrite(delta2,     "../data/UKBBsig_PC2_IMD_Delta.tsv", sep="\t")
fwrite(delta_var2, "../data/UKBBsig_PC2_IMD_Var.tsv", sep="\t")
