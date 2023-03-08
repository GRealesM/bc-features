#####################################################################
###    Blood cell features - Main analyses and figure creation    ###
#####################################################################

# Guillermo Reales
# Date of last update: 07/03/2023

# This script contains the code for the main analyses and figure and table production for the paper.

#############################
###     Load packages     ###
#############################

library(data.table)
library(magrittr)
library(pheatmap)
library(reshape2)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(gridExtra)
library(stringr)

#############################
###     Load datasets     ###
#############################


load("../data/cell-basis-sparse-3.0.RData")
dictionary <- read.table("../data/Astle_dict.txt", header = T, sep = "\t")

p.table <- fread("../data/Projection_cell_basis_v3_20220804-v1.tsv")
qc.table <- fread("../data/QC_cell_basis_v3_20220804-v1.tsv")
metadata <- fread("../data/Metadata_20230307-v1.tsv")
qc.table <- merge.data.table(qc.table, metadata, all.x = TRUE)
qc.table[, c("Collection", "Chip", "File_ID"):=NULL]


if(is.character(p.table$Var.Delta)){
  p.table$Var.Delta <- as.numeric(p.table$Var.Delta)
}

# Load the SNPs required to create the QC plots
snp.qc <- fread("../data/BCB3_SNPfiltered_20220809.tsv")


# Before we start filtering, save (supplementary) tables with info about *all* datasets

# fwrite(qc.table, "../tables/Table_Sx_Full_QC.tsv", sep = "\t")
# fwrite(p.table,  "../tables/Table_Sx_Full_projection.tsv", sep = "\t")

##################################
### Define helper functions    ###
##################################


make.2labels <- function(x){
  x[,Label:=Trait_long] %>% 
    .[, Label:=gsub(" \\(UKBB\\)", "", Label)] %>% 
    .[, Label:=gsub(" \\(FinnGen\\)", "", Label)] %>% 
    .[, Label:=gsub(" \\(FG\\)", "", Label)] %>% 
    .[, Label:=gsub("Crohn disease \\(strict definition, all UC cases excluded\\)", "Crohn's, strict", Label)] %>% 
    .[, Label:=gsub("Eosinophilic Granulomatosis with Polyangiitis", "EGPA", Label)] %>% 
    .[, Label:=paste0(str_trunc(Label, width = 50), " / ", First_Author)] %>% 
    .[, Label:=gsub(" levels", "", Label)] %>% 
    ## Here we fix some biomarker labels for clarity in the figures
    .[, Label:=gsub("Alanine aminotransferase", "Alanine aminotransferase (ALT)", Label)] %>% 
    .[, Label:=gsub("Aspartate aminotransferase", "Aspartate aminotransferase (AST)", Label)] %>% 
    .[, Label:=gsub("Interleukin-27", "Interleukin-27 (IL27)", Label)] %>% 
    .[, Label:=gsub("Gamma glutamyltransferase", "Gamma glutamyltransferase (GGT)", Label)] %>% 
    .[, Label:=gsub("Eosinophil cationic protein", "Eosinophil cationic protein (ECP)", Label)] %>%
    .[, Label:=gsub("(LOX)", "(LOX1)", fixed = T, Label)] %>% 
    
    .[, Label:=gsub( ".+\\((.+)\\)", "\\1", Label, perl = TRUE)] %>% 
    .[, Label:=gsub("Macrophage inflammatory protein-1α (CCL3, MIP1a...", "MIP1a", fixed = T, Label)] %>%  
    .[, Label:=gsub("C-X3-C motif ligand 1 ...", "CX3CL1", fixed = T, Label)] %>%
    .[, Label:=gsub("Macrophage inflammatory protein-1β (CCL4, MIP1b...", "MIP1b", fixed = T, Label)] %>%  
    .[, Label:=gsub("Matrix metalloproteinase-10 (MMP-10, Stromelysi...", "MMP-10", fixed = T, Label)] %>% 
    .[, Label:=gsub("Platelet endothelial cell adhesion molecule (PE...", "PECAM-1", fixed = T, Label)] %>% 
    .[, Label:=gsub("S100 calcium-binding protein A12 (S100A12, ENRA...", "PECAM-1", fixed = T, Label)] %>%   
    .[, Label:=gsub("Platelet-derived growth factor subunit B (PDGFb...", "PDGFb", fixed = T, Label)] %>%  
    .[, Label:=gsub("Tumor necrosis factor superfamily member 14 (TN...", "TNFSF14", fixed = T, Label)] %>%  
    .[, Label:=gsub("Proteinase-activated receptor 1", "PAR1", Label)] %>%  
    .[, Label:=gsub("Proheparin-binding EGF-like growth factor (HB-E...", "HB-EGF", fixed = T, Label)] %>% 
    .[, Label:=gsub("Dickkopf-related protein 1", "DKK1", Label)] %>%
    .[, Label:=gsub("C-reactive protein", "CRP", Label)] %>% 
    .[, Label:=gsub("HGF, SF", "HGF", Label)] %>%
    .[, Label:=gsub("ST2, IL1RL1", "ST2", Label)] %>% 
    .[, Label:=gsub("IL8, CXCL8", "IL8", Label)]
  # Now create a second label including population
    fl <- tstrsplit(x$Label, split = " / ", fixed = T)[[1]]
    x[, Label2 := paste0(fl, " (", Population, ") / ", First_Author)] 
  
  # Note: this function assumes we're using PanUKBB, so it doesn't replace PanUKBB from labels
}



##################################
### Filtering datasets and QC  ###
##################################

# Remove unnecessary columns from qc.table
qc.table[, c("URL", "Retrieval_date", "Public", "Notes"):=NULL]


# Note HBC_Chen_6 and EOSC_Chen_3 have some issues that we couldn't resolve so far, so we'll exclude them
##### & Trait!="EOSC_Chen_32888493_3"
qc.table <- qc.table[ Trait!="HBC_Chen_32888493_6" ]
p.table <- p.table[ Trait!="HBC_Chen_32888493_6" ]

# Some East Asian Chen datasets showed some weird projection values as well. We ran a PCA on the basis SNP BETAs for all Chen datasets (see Checking_Chen_raw_data.R)
# This PCA confirmed that, for some reason, East Asian Chen had some issue, so we'll remove them
qc.table <- qc.table[!grepl("Chen_32888493_3", Trait)]
p.table <- p.table[!grepl("Chen_32888493_3", Trait)]


# Checking missing projections
summary(qc.table)
qc.table[is.na(qc.table$overall_p),] 
# This time we have only one file MDD_Wray_2, Major depression disorder, that failed to project by low SNP match.
# MDD_Wray is a Meta-analysis, containing data from multiple studies.
# This particular file contained 23andMe data, which due to privacy reasons, they only made public 10k SNPs, which explains the low match.
# Fortunately, we have another file of the same meta-analysis, including many more SNPs.

qc.table <- qc.table[!is.na(qc.table$overall_p),] # Remove datasets with missing overall p-values

# Check SNP match
c(lessthan95 = nrow(qc.table[qc.table$nSNP < nrow(SNP.manifest)*.95,]), lessthan80 = nrow(qc.table[qc.table$nSNP < nrow(SNP.manifest)*.8,]), lessthan50 = nrow(qc.table[qc.table$nSNP < nrow(SNP.manifest)*.5,]))
c(lessthan95 = nrow(qc.table[qc.table$nSNP < nrow(SNP.manifest)*.95,])/nrow(qc.table), lessthan80 = nrow(qc.table[qc.table$nSNP < nrow(SNP.manifest)*.8,])/nrow(qc.table), lessthan50 = nrow(qc.table[qc.table$nSNP < nrow(SNP.manifest)*.5,])/nrow(qc.table))
# Only 1.96% of datasets (207) have <80% SNP match


# Dataset filtering

# This is a crucial point.
# As we have many datasets, many of them redundant or not so interesting for our purposes, and the set of datasets included in the adjustment will influence the resulting set of overall significant datasets, so we must choose carefully which datasets we remove *before* performing the FDR adjustment on the next step.
# 
# At this step, we'll:
# 
# 1.  Extract basis projections from proj.table, for them to be used later.
# 2.  Focus on files that came from a Genome-wide array (ie. not Immunochip), and have at least **80%** of matching SNPs, so we can have confidence in our observations.
# 3.  Remove all FinnGen, Neale UKBB and PanUKBBR1 (as we have more recent PanUKBBR2) datasets, as Astle (basis datasets).
# 4.  Create an additional category for blood cells.

# 1. Extract basis projections for later
pbasis <- p.table[grepl("_Astle_", Trait)][,Trait:=gsub("([A-Z]+)_.+", "\\1", Trait, perl = TRUE)]
pbasis <- merge(pbasis, unique(qc.table[, c("Trait_ID_2.0", "Trait_long")]), by.x = "Trait", by.y = "Trait_ID_2.0")
pbasis[,Label:=Trait_long][, stars:=""][, c("Var.Delta", "z", "P") := NULL]

# 2. Remove datasets with <80% SNP match
qc.filt <- qc.table[nSNP >= max(nSNP) * 0.8]

# 3. Remove Astle (since we used them to build the basis), UKBB (Neale) and FinnGen
qc.filt <- qc.filt[!Reference %in% c("PanUKBBR1", "UKBB", "FinnGenR5", "FinnGenR7") & First_Author != "Astle"]

# 4. Define an extra class (Blood cells) to distinguish them from the rest of biomarkers
btraits <- c(dictionary$Trait, "RBC")
# We'll consider here, in order: Chen, Kanai, Ferreira (19...), Ferreira (20...), Vuckovic (32...).
# Consider all traits that match a basis trait, so we can catch those with mixed traits (blood cells and cytokines)
qc.filt[ Reference %in% c("32888493", "29403010", "19853236", "20045101", "32888494"), Trait_class:="BC" ][ Trait_ID_2.0 %in% c(btraits), Trait_class:="BC"]

# We need to coerce PanUKBBR2 blood cells to BC class. I identified 31 of the 35 "c[0-9]+" traits to be blood-cell related, basically all those below c30500
cnum <- gsub("c", "", qc.filt[grepl("c[0-9]+", Trait_ID_2.0, perl = TRUE) & Reference == "PanUKBBR2", Trait_ID_2.0]) %>% as.numeric %>% .[. < 30500] %>% paste0("c", .)
qc.filt[ Trait_ID_2.0 %in% cnum & Reference == "PanUKBBR2", Trait_class := "BC"]

# Check traits by class
table(qc.filt$Trait_class)
# Create filtered projection table too
pt.filt <- merge(p.table, qc.filt[,c("First_Author","Trait", "Trait_ID_2.0", "Trait_long","Trait_class", "Population")], by = "Trait")



##################################
####   Apply FDR threshold    ####
##################################

# Apply 1% FDR correction to overall p for all remaining datasets
qc.filt[, FDR.overall := p.adjust(overall_p, method = "BH"), by="Trait_class"]
qc.sig <- qc.filt[FDR.overall < 0.01,]
table(qc.sig$Trait_class)

# Apply 5% FDR correction by trait and PC to projections, then filter by overall significant traits
pt.filt[, FDR.PC:=p.adjust(P, method = "BH"), by = c("PC", "Trait_class")][, stars:=ifelse(!is.na(FDR.PC) & FDR.PC<0.05,"○","")][ FDR.PC < 0.01 , stars:="●"]
# Filter PT80 by overall significant traits
pt.sig <- pt.filt[Trait %in% qc.sig$Trait ]

# Make labels for projection table
pt.sig <- make.2labels(pt.sig)

# Save a modified version of qc.filt as Table SXX, which will provide information on the datasets used in the study.


# Save our datasets so far
# These datasets will be used in the InteractiveBCB app
# fwrite(pbasis, "../data/pbasis_202303_UKBB.tsv", sep = "\t")
# fwrite(qc.sig, "../data/QC_sig_202303_UKBB.tsv", sep = "\t")
# fwrite(pt.sig, "../data/PT_sig_202303_UKBB.tsv", sep = "\t")


#-------------------------------------------------------------------------

##################################
####        QC plots          ####
##################################


# Keep only significant datasets
snp.qc <- snp.qc[filename %in% qc.sig$Trait]
snp.qc <- merge(snp.qc, unique(pt.sig[, .(Trait, Label)]), by.x="filename", by.y="Trait")

seplots <- lapply(unique(snp.qc$SNPID), function(snp){
  snp.subset <- snp.qc[SNPID == snp,]
  snp.subset <- na.omit(snp.subset, cols= c("SE","N0N1Calc"))
  snp.subset <- snp.subset[SE != 0,]
  lm.snp.subset <- lm(log(SE) ~ log(N0N1Calc), data = snp.subset)
  p  <- ggplot(snp.subset, aes(x = log(N0N1Calc), y = log(SE), colour = BETA)) +
    geom_point()+
    labs(title = snp)+
    xlab("log(Nc)")+
    theme_classic()
})

plot_grid(plotlist=seplots, ncol=2)

# ggsave("../figures/Figure_S1_QCplot1.png", bg="white", width = 7, height = 7)




##################################
####  Chen heatmap (Figure 1) ####
##################################

# Prepare data for heatmap
pt.chen <- pt.sig[ First_Author == "Chen"][!grepl("Chen_32888493_6", Trait)] # Remove second meta-analysis

PCorder <- paste0("PC", 1:14)
hmcol <- rev(colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#F7F7F7", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))(100))
Mchen.sig <- acast(pt.chen[,c("PC", "Label2", "Delta")], Label2 ~ PC) # PC, Trait, and Delta columns only
Mchen.sig.stars <- acast(pt.chen[,c("PC","Label2","stars")], Label2 ~ PC)
Mchen.sig <- Mchen.sig[,PCorder]
Mchen.sig.stars <- Mchen.sig.stars[,PCorder]
range <- max(abs(Mchen.sig))

# Colour-code by population
annchenviz <- data.frame(pop = gsub(pattern = " \\(.*", replacement = "", pt.chen[PC == "PC8"]$Population), row.names = pt.chen[PC == "PC8"]$Label2) # Note: I chose PC8 arbitrarily, just to avoid row name repetition if including many PCs
anc_colours <- c("#F7EF81", "#63a375", "#586994", "#ffab44", "#e86668", "#C5D5EA")
annchencols <- list(pop = anc_colours[c(3, 5,  2,  4)]) # keep consistency with doughnut plot. adapted to remove East Asians
names(annchencols$pop) <- unique(annchenviz$pop)

# Create heatmap
chmp <- pheatmap(Mchen.sig,  breaks = seq(-range, range, length.out = 100), cluster_cols = FALSE, display_numbers = Mchen.sig.stars, fontsize_row = 8.4, fontsize_number = 11, color = hmcol, annotation_row = annchenviz, annotation_colors = annchencols, annotation_names_row = FALSE, annotation_legend = TRUE)

# Saving function
save_pheatmap_svg <- function(x, filename, width=16, height=16) { 
  svg(filename, width = width, height = height, bg = "white") # Warning: Size in inches, can't change it to pixels
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

# Save raw heatmap, which will be manually modified to make Figure 1
# save_pheatmap_svg(chmp, "../Figures/Figure_1_Chen_heatmap_raw.svg", width = 9.5, height = 12.5)


########################################
####  Trait distribution (Figure 2) ####
########################################



dsp <- c("#001C7F", "#B1400D", "#12711C", "#8C0800", "#591E71", "#592F0D", "#A23582", "#3C3C3C", "#B8850A", "#006374") # Dark sea palette
#msp <- c("#4878D0", "#EE854A", "#6ACC64", "#D65F5F", "#956CB4", "#8C613C", "#DC7EC0", "#797979", "#D5BB67", "#82C6E2") # Muted sea palette
dark.traits <- c(BC = dsp[5], BMK = dsp[3], IMD = dsp[4], INF = dsp[9], CAN = dsp[7], OTH =dsp[1], PSD=dsp[2])
#muted.traits <- c(BC = msp[5], BMK = msp[3], IMD = msp[4], INF = msp[9], CAN = msp[7], OTH =msp[1], PSD=msp[2])

## Text-friendly palette attempt

#tfp <- c(BC = "#AA0000", BMK = "#008040", IMD = "#0F4880", INF = "#AF851A", CAN = "#E73C4E", OTH ="#1C2833", PSD="#708090")
tfp <- c(BC = "#CF000F", BMK = "#2E8856", IMD = "#1460AA", INF = "#B8860B", CAN = "#E65722", OTH ="#1C2833", PSD="#708090")


tcsig <- ggplot(data = qc.filt, aes(x=Trait_class, fill = Trait_class)) +
    geom_bar()+
    scale_y_continuous(breaks = seq(0, 5000, 200))+
    scale_fill_manual(values = tfp)+
    scale_x_discrete(labels=c("BC" = "Blood Cells", "BMK" = "Biomarkers", "CAN" = "Cancer", "IMD" = "IMD", "INF" = "Infectious", "OTH"= "Other", "PSD" = "Psych"))+
    theme_classic(base_size = 11)+
    ylab("Datasets") +
    theme(axis.ticks = element_blank(), axis.line = element_blank(), panel.grid.major.y = element_line(colour = "lightgray", size =  0.3),axis.text.x = element_text(angle = 270, hjust = 0, vjust = 0.5),legend.position = "none",axis.title.x = element_blank())



# A little modification to include OR values in the RHS plot (see how values are computed below)

qcsigOR <- merge(qc.filt[, .(count.filt=.N), by=Trait_class], qc.sig[, .(count.sig=.N), by=Trait_class], all.x = TRUE)
qcsigOR[,prop:=count.sig/count.filt][, `OR vs Other`:= (count.sig/(count.filt - count.sig)) / qcsigOR[Trait_class == "OTH", count.sig/(count.filt - count.sig)]]
qcsigOR[, ORvalue:=format(`OR vs Other`, digits=2)]
qcsigOR <- qcsigOR[Trait_class != "PSD"] 
# qcsigOR[ Trait_class %in% c("OTH", "BC"), ORvalue:=""] # Adapt for plotting


tcsigOR <- ggplot(data = qcsigOR, aes(x=Trait_class, y=count.sig, fill = Trait_class)) +
  geom_bar(stat = "identity")+
  scale_y_continuous(breaks = seq(0, 500, 50), limits = c(0, 180))+
  scale_fill_manual(values =tfp)+
  scale_x_discrete(labels=c("BC" = "Blood Cells", "BMK" = "Biomarkers", "CAN" = "Cancer", "IMD" = "IMD", "INF" = "Infectious", "OTH"= "Other", "PSD" = "Psych"))+
  geom_text(aes(label = ORvalue), vjust = -0.8, hjust = 0.58, size=4)+
  theme_classic(base_size = 11)+
  theme(axis.ticks = element_blank(), axis.line = element_blank(), panel.grid.major.y = element_line(colour = "lightgray", size =  0.3),axis.text.x = element_text(angle=270, hjust = 0, vjust = 0.5),legend.position = "none",axis.title.x = element_blank(), axis.title.y = element_blank())


sigct <- pt.sig[ FDR.PC < 0.01 , .(count = .N), by=c("PC", "Trait_class")][, PC:=factor(PC, levels = paste0("PC", 1:14))][ Trait_class !="BC"]

sig.per.pc <- ggplot(sigct, aes(fill=Trait_class, y=count, x=PC))+
              geom_bar(position="stack", stat="identity")+
              scale_y_continuous(breaks = seq(0, 100, 25), limits = c(0, 105))+
              scale_fill_manual(values =tfp)+
              theme_classic(base_size = 11)+
              theme(axis.ticks = element_blank(), axis.line = element_blank(), panel.grid.major.y = element_line(colour = "lightgray", size =  0.3),axis.text.x = element_text(angle=270, hjust = 0, vjust = 0.5),legend.position = "none",axis.title.x = element_blank(), axis.title.y = element_blank())

top <- plot_grid(tcsig, # Before FDR, should be the same traits for all bases. For comparison purposes
                 tcsigOR, nrow = 1, labels = "AUTO")
bottom <- plot_grid(sig.per.pc, nrow = 1, labels = "C")
plot_grid(top, bottom, nrow = 2)

# ggsave("../figures/Figure_2_traitdistro.png", bg="white", width = 8, height = 8)


#########################################
####  Delta plots PC5+PC8 (Figure 3) ####
#########################################


# Define a more customisable function for individual plotting
make.custom.forest.plots <- function(ps, feature, palette = tfp, fsize = NULL, threshold_basis=0, threshold_gen = 0, remove_sum = FALSE, remove_perc = FALSE, remove_side = "none"){
  ds.col <- data.table(Trait_class = names(tfp), colours=tfp)
  # Need to update names in tfp for scale_colour_manual
  tfp.manual <- tfp
  names(tfp.manual) <- tfp.manual
  
  pdelta <- merge(ps, ds.col, by = "Trait_class", all.x=TRUE)
  
  ft <- paste0("PC", feature)
  dt <- pdelta[PC == ft][order(Delta, decreasing = TRUE)][!(Trait_class == "BC" & abs(Delta) < threshold_basis)][ abs(Delta) > threshold_gen]
  
  if(remove_sum) dt <- dt[!grepl("Sum ", Label, ignore.case = TRUE)]
  if(remove_perc) dt <- dt[!grepl("percentage of", Label, ignore.case = TRUE)]
  if(remove_side == "positive") dt <- dt[Delta <= 0]
  if(remove_side == "negative") dt <- dt[Delta >= 0]
  
  fplot <- ggplot(dt, aes(x = reorder(Label, -Delta), y = Delta, ymin=Delta-ci, ymax=Delta+ci, colour = colours))+
    geom_pointrange()+
    #scale_colour_manual(values = c("red" = "red", "#26547c" = "#26547c", "#049F76" = "#049F76", "#E09D00" = "#E09D00", "#F3B8A5" = "#F3B8A5", "black" = "black"))+
    scale_colour_manual(values = tfp.manual)+
    geom_hline(yintercept = 0, col="red", lty=2)+
    coord_flip()+
    #xlab("Traits")+
    ylab("Delta")+
    ggtitle(ft)+
    theme_minimal()+
    theme(axis.text.y = element_text(colour = dt$colours, size = fsize), legend.position = "none", axis.title.y = element_blank(), panel.grid.minor.x = element_blank())
  fplot
}

tfp <- c(BC = "#CF000F", BMK = "#2E8856", IMD = "#1460AA", INF = "#B8860B", CAN = "#E65722", OTH ="#1C2833", PSD="#708090")

mm <- fread("../data/tmp_processed_redundancies_202303_UKBBsig.tsv") %>% .[is.na(Redundant), Trait]

pchenb <- pt.sig[ grepl("32888493_5", Trait)][ FDR.PC < 0.01] # Trans-ethnic Chen only
pchenb[, Label:=paste0(Trait_long, " / Chen")]

ptr <- pt.sig[Trait %in% mm]

# Prepare dataset
pdelta <- ptr[stars == "●"][ Trait_class != "BC"] # Remove Blood cells to avoid cluttering
pdelta <- rbind(pdelta, pchenb, fill = TRUE)
pdelta[Var.Delta != 0, ci:=sqrt(Var.Delta) * 1.96][is.na(ci),ci:=0][is.na(Label), Label:=Trait]


cp5 <- make.custom.forest.plots(pdelta, feature = 5, palette = tfp, fsize = 9, threshold_basis=0, remove_sum = TRUE, remove_perc = TRUE, remove_side = "none")
cp8 <- make.custom.forest.plots(pdelta, feature = 8, palette = tfp, fsize = 9,  threshold_basis=0, remove_sum = TRUE, remove_perc = TRUE, remove_side = "none")

ggsave("../figures/Figure_S2_PC5_deltaplot.png", cp5, bg="white", width = 8, height = 8)
ggsave("../figures/Figure_S3_PC8_deltaplot.png", cp8, bg="white", width = 8, height = 10)

# Main figure with selected traits
excl <- "HDL|LDL|Glucose|Testosterone|Alkaline phosphatase|ALT |AST |MMP-10|Creatinine|Apolipoprotein|bilirubin|GGT|Phosphate|Cathepsin|Calcium|Cystatin|eGFR|SHBG|Total protein|Non-albumin|Triglycerides|Mean platelet|Mean corpuscular|Distribution Width|Type 2 Diabetes|High Cholesterol|Hypertension|Hemoglobin concentration|Hematocrit|HbA1c|Cholesterol|Albumin|Red blood cell|Urate|Leukocyte"

pdr <- pdelta[!grepl(excl, Label, ignore.case = TRUE)]

cp5r <- make.custom.forest.plots(pdr, feature = 5, palette = tfp, fsize = 9, threshold_basis=0, remove_sum = TRUE, remove_perc = TRUE, remove_side = "none")
cp8r <- make.custom.forest.plots(pdr, feature = 8, palette = tfp, fsize = 9, threshold_basis=0, remove_sum = TRUE, remove_perc = TRUE, remove_side = "none")

c58r <- plot_grid(cp5r, cp8r, nrow = 1, labels = "AUTO")
# ggsave("../figures/Figure_3_deltaplots58.png", c58r, bg="white", width = 11, height = 7.5)

#########################################
####   Delta plots PC2 (Figure 5)    ####
#########################################

cp2 <- make.custom.forest.plots(pdelta, feature = 2 , palette = tfp, threshold_basis=0, remove_sum = TRUE, remove_perc = TRUE, remove_side = "none")

# ggsave("../figures/Figure_S4_PC2_deltaplot.png", cp2, bg="white", width = 8, height = 7.5)

cp2r <- make.custom.forest.plots(pdr, feature = 2, palette = tfp, fsize = 9, threshold_basis=0, remove_sum = TRUE, remove_perc = TRUE, remove_side = "none")

# ggsave("../figures/Figure_5_deltaplot2.png", cp2r, bg="white", width = 6, height =5)
