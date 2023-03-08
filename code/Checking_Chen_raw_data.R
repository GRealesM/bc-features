###############################################################
#####     Checking Chen Raw data across populations    ########
###############################################################

# Date: 2022-11-02
# Author: Guillermo Reales

# Description: Some Chen projections (especially East Asians) seem off. We'll take a look at the raw data by taking the BETA estimates
# for all datasets and basis SNPs and applying a PCA on them

library(data.table)
setDTthreads(10)
library(magrittr)
library(ggplot2)
library(ggrepel)
library(IMDtools)
library(cowplot)

load("../basis_building/cell-basis-sparse-3.0.RData")

man <- SNP.manifest

d <- data.table()
ptofiles <- "../../cell_basis_v2/reduced_datasets/"

chenfiles <- grep("Chen_", dir(ptofiles), value = TRUE)

ds <- lapply(chenfiles, function(x){
       
       dd <- fread(paste0(ptofiles, x))
       dd[, dataset:=x]
})
ds <- rbindlist(ds)

# Check if all SNPs match the manifest and vice versa (they should!)
all(unique(ds$pid) %in% man$pid)
all(man$pid %in% unique(ds$pid))

# Count SNPs and remove those SNPs that aren't present in all datasets
length(unique(ds$dataset)) # 89 datasets
incomp.SNP  <- ds[, .(count = .N), by=pid][order(count, decreasing = TRUE)][ count < max(count), pid ] # 42/1661 are absent in at least one dataset.

ds <- ds[!pid %in% incomp.SNP] # Remove incomplete SNPs
ds[, .(count = .N), by=pid][order(count, decreasing = TRUE)]

# Create appropriate Labels
ds[, label:=gsub("-ft.tsv", "", dataset)][, label:=gsub("Chen_32888493_", "", label)]
ds[, label:=gsub("1", "AFR", label)][, label:=gsub("2", "EUR", label)][, label:=gsub("3", "EAS", label)][, label:=gsub("4", "HIS", label)][, label:=gsub("5", "META1", label)][, label:=gsub("6", "META2", label)]

# We saw that HBC_META2 (HBC_Chen_6) that was problematic, let's take a look at it.
ds[ label == "HBC_META2"]
system("zcat ~/rds/rds-cew54-basis/02-Processed/HBC_Chen_32888493_6-hg38.tsv.gz | head | awk '{ print $10}'") # Estimates are weird in origin! 
# These are sdY-transformed values, which in origin some of them had extremely low SE, resulting in absurd BETA estimates. 
# See /home/gr440/rds/rds-cew54-basis/GWAS_tools/special_cases/exploring_HBC6.R for more details.

# We'll remove these datasets for now
ds <- ds[ label != "HBC_META2"]

# The function below will take the dataset. Turn it into a matrix, perform a PCA and plot it

pca.plot  <- function(ds, plot.legend = TRUE){
        # Create matrix 
        matd <- dcast(ds, label ~ pid, value.var = "BETA" )
        mat <- as.matrix(matd[, 2:ncol(matd)])
        rownames(mat) <- matd$label
        dim(mat)

        # Perform PCA
        pca <- stats::prcomp(mat, center = TRUE, scale = FALSE)

        plot.DT <- data.table(label = rownames(pca$x), pca$x)
        plot.DT <- plot.DT[, 1:11] # Take the first 10 PCs

        plot.DT[, ancestry := tstrsplit(label, "_")[2]][, celltype:=tstrsplit(label, "_")[1]]

        # Plot!
        biplot <- ggplot(plot.DT, aes(PC1,PC2,  label= celltype, colour = ancestry))+
          geom_point()+
          geom_label_repel(data = plot.DT[ancestry == "EAS" | PC1 > 6], size=3,force = 20, seed = 1) + 
          theme_minimal()

       if(plot.legend == FALSE){
       biplot <- ggplot(plot.DT, aes(PC1,PC2,  label= celltype, colour = ancestry))+
          geom_point()+
          geom_label_repel(data = plot.DT[ancestry == "EAS" | PC1 > 6], size=3,force = 20, seed = 1) + 
          theme_minimal()+
          theme(legend.position = "none")
       }

        biplot 

}

biplot1 <- pca.plot(ds)
biplot1
# Here we observe that PLAV_META2 is an outlier, too. We'll remove it and try again.

ds2 <- ds[label != "PLAV_META2"]
biplot2 <- pca.plot(ds2, plot.legend = FALSE)
biplot2
# This confirms that EAS projections are very different from the rest, for whatever reason.

#ggsave( "../Figures/Manuscript/Figure_SXX_Chen_rawdata_biplot.png", biplot2, units = "in", height = 8, width = 8, bg="white")


# What if we tried to do the same but using pre-sdY correction data?
# We repeated the reduction procedure but keeping the original values, instead of sdY-transformed ones.

porigfiles <- "~/rds/rds-cew54-basis/03-Bases/cell_basis_v2/reduced_orig_chen/"
chenrawfiles <- grep("Chen_", dir(porigfiles), value = TRUE)

dso <- lapply(chenrawfiles, function(x){
       
       dd <- fread(paste0(porigfiles, x))
       dd[, dataset:=x]
})
dso <- rbindlist(dso)

# Check if all SNPs match the manifest and vice versa (they should!)
all(unique(dso$pid) %in% man$pid)
all(man$pid %in% unique(dso$pid))

# Count SNPs and remove those SNPs that aren't present in all datasets
length(unique(dso$dataset)) # 89 datasets
incomp.SNP  <- dso[, .(count = .N), by=pid][order(count, decreasing = TRUE)][ count < max(count), pid ] # 42/1661 are absent in at least one dataset.

dso <- dso[!pid %in% incomp.SNP] # Remove incomplete SNPs
dso[, .(count = .N), by=pid][order(count, decreasing = TRUE)]

# Create appropriate Labels
dso[, label:=gsub("-ft.tsv", "", dataset)][, label:=gsub("Chen_32888493_", "", label)]
dso[, label:=gsub("1", "AFR", label)][, label:=gsub("2", "EUR", label)][, label:=gsub("3", "EAS", label)][, label:=gsub("4", "HIS", label)][, label:=gsub("5", "META1", label)][, label:=gsub("6", "META2", label)]


biplot3 <- pca.plot(dso)
biplot3

togplot <- plot_grid( biplot2, biplot3, nrow = 1, labels = "AUTO")

ggsave( "../Figures/Manuscript/Figure_SXX_Chen_rawdata_biplot.png", togplot, units = "in", height = 7, width = 11, bg="white")
