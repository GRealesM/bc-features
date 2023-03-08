################################################
### Allocating biomarkers to IMD clusters   ####
################################################

# Background: Biomarkers are invariably in a different scale than IMD, making joint clustering using DPMUnc challenging. 
# Here, we computed ellipses around IMD clusters to allocate biomarkers to their closest IMD cluster. See paper's Methods for more details.

# Set working directory. Adapt as needed. All paths will be relative to this.
setwd("/home/gr440/rds/rds-cew54-basis/03-Bases/cell_basis_v3_varimax/DPMUnc/code/ellipses")

## Load required packages

source("ellipse-tangents.R") # Custom functions to compute ellipses
library(data.table)
library(magrittr)
library(Matrix)
# devtools::install_github("chr1swallace/random-functions")
library(randomFunctions)    # Custom Chris Wallace's functions
library(rmeta)
library(mvmeta)             # For multivariate meta-analysis
library(ggplot2)
library(ggalluvial)
library(pheatmap)


## Create paths

proj.path <-  "../../../Projections/"
DPres.path <- "../../psm_results/"
fig.path <- "../../../Figures/Manuscript/"
res.path <- "../../ellipses_results/" # Where to store results

## Load datasets

pdata <- fread(file.path(proj.path, "PT_sig_202211_UKBB.tsv")) # Filtered Projection table, after selecting files for downstream analyses,
# see https://grealesm.github.io/Bases/cell_basis_v3_varimax/Reports/Cellv3_report_20230118_UKBB.html#Which_datasets_are_significant_overall for details.
# It needs a bit of processing to remove redundant datasets and ensure it matches the DPMUnc files
nr.ds <- fread("../../data/tmp_processed_redundancies_202301_UKBBsig.tsv") %>% .[ is.na(Redundant), Trait] # Redundant were marked with 1, and non-redundant were NA. 
# 87 datasets

pdata <- pdata[ Trait %in% nr.ds & Trait_class %in% c("IMD", "BMK") & PC %in% paste0("PC", c(5, 8))]  
# Keep non-redundant traits only.
# Keep IMD and BMK only
# Keep only PC5/8 projections
sig <- unique(pdata[ FDR.PC < 0.01, Trait]) # PC-significant for at least one PC
pdata  <- pdata[Trait %in% sig]

# DPMUnc IMD results
pc58 <- readRDS(file.path(DPres.path, "UKBBsig_PC58_IMD_psm_data.rds"))


# Check we're considering the same stuff
IMDtraits  <- unique(pdata[ Trait_class == "IMD", Label])
all(IMDtraits == names(pc58$calls$cl)) # TRUE


################################################################################
#                    USING ELLIPSES TO ASSIGN CLUSTERS                         #
#                                                                              #
################################################################################

## FUNCTIONS

## Projection
## Some Helper functions for projecting
calc_var.proj <- function (beta, seb, pids, pcs_wanted) {
    if (length(beta) != length(seb) || length(beta) != length(pids) ||
        !length(beta))
        stop("arguments must be equal length vectors > 0")
    if (!all(pids %in% SNP.manifest$pid))
        stop("all pids must be members of sparse basis (SNP.manifest$pid)")
    if (length(pids) < 0.95 * nrow(rot.pca))
        warning("more than 5% sparse basis snps missing")
    b <- beta * shrinkage[pids] - beta.centers[pids]
    proj <- b %*% rot.pca[pids, ]
    v <- seb * shrinkage[pids] * rot.pca[pids, ]
    var.proj <- t(v) %*% LD[pids, pids] %*% v
    var.proj[pcs_wanted, pcs_wanted]
}

# Wrapper function for projecting
get_var.proj  <- function(f, ...) {
  sm <- fread(file.path(mpath,f))
    # Some checks
  sm <- unique(sm)
  sm[sm == ""] <- NA # Some missing data might pass as empty string. This will fix that
  sm <- na.omit(sm, cols = c("pid", "BETA", "SE", "P"))
  dups <- sm$pid[duplicated(sm$pid)]

  if(length(dups) > 0){
    dupmessage= "This file has duplicated pids. I removed them prior to projection. You might want to check it."
    message(dupmessage)
    sm <- sm[!pid %in% dups] # Remove all duplicated instances, to be safe
  }
  calc_var.proj(beta=sm$BETA, seb=sm$SE, pids=sm$pid, ...)
}

## Plotting function
plot_heatmap  <- function(mat,INPUT,...) {
   ## pheatmap likes data.frames
  mat2=dcast(mat[!is.na(inferred_group)], cl + inferred_group + Label ~ PC, value.var="Delta") %>% as.data.frame()
  rownames(mat2)=mat2$Label
  mat2$cl  %<>%  as.factor()
  mat2$inferred_group  %<>%  as.factor()
  
  myColor <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#F7F7F7", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))(100)
  mat2.hm   <- mat2[, grep("PC", colnames(mat2))]
  range     <- max(abs(mat2.hm))
  myBreaks  <- seq(-range, range, length.out = 100)

  mat2=mat2[order(mat2$inferred_group),]
  pheatmap(mat2[,INPUT$pcs_wanted],
           annotation_row=mat2[,c("cl","inferred_group")],
           breaks=myBreaks, color=myColor, cluster_rows=FALSE, ...)
 }

# Calculate meta v
calc_meta_v <- function(ss, pcs_wanted, Vars) {
  tmp=dt2mat(ss[PC %in% pcs_wanted], Label ~ PC, value.var="Delta")
  d=mvmeta(tmp, lapply(Vars[rownames(tmp)], function(v) v[pcs_wanted, pcs_wanted]))
  conf_slopes=find_tangent_slopes(mu=d$coefficients, V=d$vcov)
  ## d[c("coefficients","vcov")]
  tmp=dcast(ss, Label ~ PC, value.var="Delta")
  sc=sqrt(sum(d$coefficients^2))
  for(j in seq_along(pcs_wanted)) {
    tmp[[ paste0(pcs_wanted[j],".mu") ]] = d$coefficients[[j]]/sc
    tmp[[ paste0(pcs_wanted[j],".se") ]] = sqrt(diag(d$vcov))[j]/sc
  }
  tmp$lower_slope <- conf_slopes[1]
  tmp$upper_slope <- conf_slopes[2]
  attr(tmp, "meta") <- d[c("coefficients","vcov")]
  tmp
}

## Find matches based on overlap
choose_vector <- function(label, ss, AVG, pcs_wanted=INPUT$pcs_wanted,  Vars, do_plot=FALSE) {

  sstmp=ss[[label]][ match(pcs_wanted, PC) ]
 # if(!any(sstmp$fdr < fdr_threshold)) # No need for new fdr
 #   return(rep(FALSE,length(AVG)))
  Vtest=as.matrix(Vars[[label]][pcs_wanted, pcs_wanted])
  mutest=sstmp$Delta
  Etest=ellipse( Vtest, centre=mutest)
  ## disease cluster vectors
  overlap=sapply(seq_along(AVG), function(j) {
    E=cbind(Etest, yl=AVG[[j]]$lower_slope[1] * Etest[,1], yu=AVG[[j]]$upper_slope[1] * Etest[,1])
    ## check for overlap
    !(all(E[,2] < E[,"yl"] & E[,2] < E[,"yu"]) || all(E[,2] > E[,"yl"] & E[,2] > E[,"yu"]))
  })
   angle_ok=sapply(seq_along(AVG), function(j) {
    cf=attr(AVG[[j]],"meta")$coefficients[,pcs_wanted]
    angle_disease=atan2(cf[2],cf[1])
    angle_test=atan2(mutest[2],mutest[1])
    min(2*pi - abs(angle_disease - angle_test), 
        abs(angle_disease - angle_test)) < 2 * pi/3
  })
  if(do_plot) {
    Elist=lapply(AVG, function(A)
      ellipse(attr(A,"meta")$vcov, centre=attr(A,"meta")$coefficients))
    limits=c(min(c(Etest,unlist(Elist))), max(c(Etest,unlist(Elist))))*1.1
    plot(Etest, type="l", asp=1, xlim=limits, ylim=limits)
    points(0,0,pch="o")
    for(j in seq_along(Elist)) {
      lines(Elist[[j]], col=j+1)
      abline(a=0, b=AVG[[j]]$lower_slope[1], col=j+1)
      abline(a=0, b=AVG[[j]]$upper_slope[1], col=j+1)
     }
  }
  overlap & angle_ok
  }

# Function to do everything in one step
# I introduced some changes and comments for clarity.
# Now I included the VARS in the Input object.
assign_biomarkers <- function(INPUT){

  # Run some checks
  stopifnot(is.list(INPUT))
  stopifnot(all(names(INPUT) == c("data", "pcs_wanted", "calls", "VARS"))) 
  # INPUT must contain
  # data: a data.table containing the projections (Delta, Var.Delta)
  # pcs_wanted: a vector of the PCs that we want to include in the run (ie. c("PC5","PC8") or c("PC5", "PC8", "PC2", "PC5"))
  # calls: a named integer with the DPMUnc-allocated clusters.
  # VARS: Projection variance for each of the traits of interest for each of the PC of interest. See calc_var.proj for more details.

  ## Calculate meta v
  mat=INPUT$data[, .(Label, PC, Delta, SE=sqrt(Var.Delta))] # We don't need filters now because the datasets were filtered above
  mat$cl=INPUT$calls[ mat$Label ]
  if(length(INPUT$pcs_wanted)==2) # Special case for PC5/PC8
    mat[cl==3,cl:=1]
  ss=split(mat,mat[,c("cl")])
  if(length(INPUT$pcs_wanted)==2) {  # Special case for PC5/PC8
    AVG=lapply(ss, calc_meta_v, INPUT$pcs_wanted, INPUT$VARS)
  } else {                           # Special case for PCgrac
    AVG1=lapply(ss, calc_meta_v, INPUT$pcs_wanted[1:2], INPUT$VARS) 
    AVG2=lapply(ss, calc_meta_v, INPUT$pcs_wanted[3:4],INPUT$VARS)
  }

  ## Find matches based on overlap
  ss=split(mat, mat$Label) # Split dataset by Label
  choices=if(length(INPUT$pcs_wanted)==2) { # Special case for PC5/PC8
            lapply(names(ss), choose_vector, ss, AVG, INPUT$pcs_wanted, INPUT$VARS)
          } else {
            lapply(names(ss), function(nm) { # Special case for PCgrac (PC2/11?)
              c(FALSE, choose_vector(nm, ss, AVG2[-1], INPUT$pcs_wanted[3:4], INPUT$VARS)) #, do_plot=TRUE)
            })
          }
  exact1=sapply(choices,sum)==1
  final_choice=unlist(ifelse(exact1==TRUE, sapply(choices,which), NA))
  names(final_choice)=names(ss)
  mat[,inferred_group:=final_choice[Label]]
  mat[!is.na(cl), inferred_group:=cl] #keep original group if disease.

  # Plot results
  plot_heatmap(mat, INPUT)
  invisible(mat)

}


## PREPARING INPUTS

## First, we need full var.proj of the selected datasets
load("~/rds/rds-cew54-basis/03-Bases/cell_basis_v3_varimax/basis_building/cell-basis-sparse-3.0.RData")
mpath <- "../../../../cell_basis_v2/reduced_datasets/"

## Prepare inputs. Keep just BMK traits that are significant for at least one of the PCs
data58 <- pdata[PC %in% c("PC5","PC8")]
sigbmk <- unique(data58[ Trait_class =="BMK" & FDR.PC < 0.01, Trait]) # Setting FDR threshold to 0.01 for at least one of the involved PCs
data58  <- data58[Trait_class == "IMD" | Trait %in% sigbmk ]
labels_wanted_58  <- unique(data58$Label)

datagr <- pdata[PC %in% c("PC2","PC5","PC8", "PC11")]
sigbmk <- unique(datagr[ Trait_class =="BMK" & FDR.PC < 0.01, Trait]) # Setting FDR threshold to 0.01
datagr  <- datagr[Trait_class == "IMD" | Trait %in% sigbmk ]
labels_wanted_grac  <- unique(datagr$Label)


table(unique(pdata[, .(Trait, Trait_class)])$Trait_class)
table(unique(data58[, .(Trait, Trait_class)])$Trait_class)
table(unique(datagr[, .(Trait, Trait_class)])$Trait_class)
# Same IMD for all datasets (ok!), 80 BMK in full data (ok!), 58 in pcgrac (ok), 51 in pc58 (ok!)

## Get projection variances
lookup=with(unique(pdata[,.(Label,Trait)]), structure(Label,names=Trait))
traits_wanted=unique(pdata$Trait)
files_wanted=paste0(traits_wanted, "-ft.tsv")

VARS=lapply(files_wanted, get_var.proj, pcs_wanted=c("PC5","PC8","PC2","PC11"))
names(VARS)=lookup[ traits_wanted ]

VARS58 = VARS[ labels_wanted_58 ]
VARS58 = lapply(VARS58, function(v) v[ c("PC5","PC8"), c("PC5","PC8") ])

VARSgr = VARS[ labels_wanted_grac ]
VARSgr = lapply(VARSgr, function(v) v[ c("PC5","PC8","PC2","PC11"), c("PC5","PC8","PC2","PC11") ])


# Then we can move on and prepare input lists for assign_biomarkers()
INPUT58=list(data = data58, pcs_wanted=c("PC5","PC8"), calls=pc58$calls$cl, VARS=VARS58)
INPUTgr=list(data = datagr, pcs_wanted=c("PC5","PC8","PC2","PC11"), calls=pcgr$calls$cl, VARS=VARSgr) # Important that the PCs are in this order

# Run assignments using assign biomarkers
result58=assign_biomarkers(INPUT58)
resultgr=assign_biomarkers(INPUTgr)


plot_heatmap(result58,INPUT58,file=file.path(fig.path,"Figure_SXX_Heatmap_BMK_assignment_PC58.png"),cluster_cols=FALSE, height=8, width=7)
plot_heatmap(resultgr,INPUTgr,file=file.path(fig.path,"Figure_SXX_Heatmap_BMK_assignment_PCgrac.png"), cluster_cols=FALSE, height=10, width=7.5)

saveRDS(result58, file.path(res.path, "result58.RDS"))
saveRDS(resultgr, file.path(res.path, "resultgr.RDS"))

