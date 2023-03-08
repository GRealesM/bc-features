#############################
# Creating PSM and heatmaps #
#############################

# Date: 2022/08/25
# Guillermo Reales

# After checking the trace plots, we'll import the results, remove the first half of the run, get all seeds together and generate the PSM and heatmaps.
# This is a more interactive script, intended to make publication-ready figures.

# Again, we'll use functions written by Kath Nicholls and Chris Wallace, but I'll adapt them to work with my file structure system.
# As a reminder, all DPMUnc results are in '../results/', with an directory naming system in the form of {exp}_{seed} (eg. PC58_ALL_1).
# The idea here is to incorporate all seeds in every experiment and generate joint traceplots. See below for function adaptations.

setwd("/home/gr440/rds/rds-cew54-basis/03-Bases/cell_basis_v3_varimax/DPMUnc/code")

# Load packages
library(data.table)
setDTthreads(10)
library(clue)
library(magrittr)
library(cluster)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(mcclust)
library(mclust)
library(pheatmap)
library(R.cache)

# Adapted helper functions from Kath's repo (https://github.com/nichollskc/clusteringPublicGeneExpr/blob/26e1fad109009563c8a0dc1254a54f796176f6c5/scripts/psm_plots.R)
# Also from scripts/utils.R (https://github.com/nichollskc/clusteringPublicGeneExpr/blob/26e1fad109009563c8a0dc1254a54f796176f6c5/scripts/utils.R)


palette <- c("#77AADD", "#000000", "#9E0142", "#D53E4F", "#F46D43", "#FEE08B", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2", "#FDAE61", "#E6F598", "#771155", "#AA4488", "#CC99BB", "#114477", "#774411", "#EEEEEE", "#117777", "#117744", "#44AA77", "#88CCAA", "#777711", "#44AAAA", "#AAAA44", "#77CCCC", "#DDDD77", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")
 
# NOTE: next function was in scripts/utils.R. Modified from original to remove burnin
calc_psms <- function(datasets, burnin) {  
    allocs=lapply(paste0(datasets, "/clusterAllocations.csv"), fread) ## read the allocations
    # This line is essential for some reason
    allocs %<>% lapply(., function(x) as.matrix(x[1:nrow(x),]))
    message("Removing burnin.")
    allocs %<>% lapply(., function(x) x[-c(1:burnin),]) # Remove burnin at this stage

    bigalloc = do.call(rbind, allocs)    
    bigpsm=calc_psm(bigalloc,burn=0) ## make a psm, don't discard any burn in because already discarded
    psms = lapply(allocs, function(x) calc_psm(x, burn=0))
    return(list(bigpsm=bigpsm, psms=psms))
}

raw_calc_mclust <- function(obsData) {
    BIC <- mclustBIC(obsData)
    Mclust(obsData, x=BIC)
}

raw_calc_psm=function(x,burn=0.5) {
  n=nrow(x)
  if(burn>0)
    x=x[ (burn*n) : n, , drop=FALSE]
  if(any(is.na(x)))
    x=x[ apply(!is.na(x),1,all), ]
  unq=unique(as.vector(x))
  ## print(unq)
  m=matrix(0,ncol(x),ncol(x))
  for(k in unq) {
    xk=matrix(as.numeric(x==k),nrow(x),ncol(x))
    ## m=m + t(xk) %*% xk
    m=m + crossprod(xk)
  }
  psm=m/nrow(x)
  psm
}

calc_psm <- addMemoization(raw_calc_psm)

calc_mclust <- addMemoization(raw_calc_mclust)

get_ann_colors=function(calls, mclust_calls, obsData, verbose=TRUE) {
  # from spectral, plus some extras
  #palette= c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788") %>%
  # matrix(.,7,3,byrow=TRUE) %>% as.vector()

  counts = data.frame(table(calls))
  cluster_labels = paste0(LETTERS[1:length(counts$Freq)], " (", counts$Freq, ")")

  mclust_calls = adjust_labels_B_to_match_A(calls, mclust_calls)
  mclust_counts = data.frame(table(mclust_calls))
  mclust_cluster_labels = paste0(letters[1:length(mclust_counts$Freq)], " (", mclust_counts$Freq, ")")

  ann=data.frame(row.names=rownames(obsData),
                 cluster=factor(calls, labels=cluster_labels),
                 mclust=factor(mclust_calls, labels=mclust_cluster_labels))
  ncalls=length(unique(calls))
  ncalls_mclust=length(unique(mclust_calls))
  ann_colors=list(cluster = structure(palette[1:ncalls], names=cluster_labels),
                  mclust = structure(palette[1:ncalls_mclust], names=mclust_cluster_labels))
  list(ann=ann,colors=ann_colors)
}

adjust_labels_B_to_match_A <- function(calls_A, calls_B) {
    K_A = length(unique(calls_A))
    K_B = length(unique(calls_B))

    if (K_A < max(calls_A)) {
        print("WARNING: assumptions about cluster labels violated")
    }
    if (K_B < max(calls_B)) {
        print("WARNING: assumptions about cluster labels violated")
    }

    jaccard_mat = matrix(0,
                         nrow=max(K_A, K_B),
                         ncol=max(K_A, K_B))
    for (i in 1:K_A) {
        for (j in 1:K_B) {
            in_A = calls_A == i
                in_B = calls_B == j
                jacc = sum(in_A & in_B) / sum(in_A | in_B)
                jaccard_mat[i, j] = jacc
        }
    }

    new_labels_for_B = c(solve_LSAP(t(jaccard_mat), maximum=TRUE))[1:K_B]

    return(plyr::mapvalues(calls_B, from=1:K_B, to=new_labels_for_B))
}

# Modified to adapt to the input (../results/) and output (../plots/) directory structure, require experiment name and burnin

# This function is modified to remove the hclust annotation and make fonts a bit larger.
psm_plots <- function(exp, burnin, focus_dataset=NULL, update_label=NULL) {

    message("Loading data...")
    datasets = dir("../results/")[grepl(exp, dir("../results/"))] # Little adaptation to original function. From experiment name we retrieve all directories for different seeds.
    datasets = paste0("../results/", datasets)

    obsData = read.table(paste0("../data/", exp, "_Delta.tsv"),
                         header=1, row.names=1, quote="", sep="\t")
    obsVars = read.table(paste0("../data/", exp, "_Var.tsv"),
                         header=1, row.names=1, quote="", sep="\t")
    
    # Keep Trait names for later match
    Trait_ids = rownames(obsData)

    # If you want to update the labels, you can do it at this point. Simply create a data.table with two columns: one "Trait" matching the "Trait"
    # column in obsData and another "Label" with the desired labels.
   if(!is.null(update_label)){
       stopifnot(is.data.table(update_label))
       stopifnot(names(update_label) == c("Trait", "Label"))
       stopifnot(all(rownames(obsData) == rownames(obsVars)))
       md <- update_label[match(rownames(obsData), Trait)]
       rownames(obsData) <- md$Label
       rownames(obsVars) <- md$Label
   }

    message("Preparing PSMs.")
    result = calc_psms(datasets, burnin)
    bigpsm = result$bigpsm
    psms = result$psms

    print(isSymmetric(bigpsm))
    print(max(bigpsm))
    print(min(bigpsm))
    print(diag(bigpsm))
    rownames(bigpsm) = rownames(obsData)
    colnames(bigpsm) = rownames(obsData)

    calls = minbinder(bigpsm, method = "comp") ## calls
    # calls=maxpear(bigpsm, method="comp") ## calls

    # MClust solution
    print("Calculating mclust solution.")
    mclust_solution <- calc_mclust(obsData)
    print(mclust_solution)

    # if(nrow(obsData) >= 20) {
    #     maxK = 20 } else { maxK = nrow(obsData)
    #     }

    # gap_stat <- clusGap(obsData,
    #                     FUN = kmeans,
    #                     nstart = 25,
    #                     K.max = maxK,
    #                     B = 50)
    # print(gap_stat)
    # best_K = which.max(gap_stat$Tab[, 3])
    # kmeans_solution = kmeans(obsData, centers=best_K, nstart=25)

    annotations = get_ann_colors(calls$cl, mclust_solution$classification, obsData)

    # Remove mclust
    annotations = list(ann = annotations$ann[1], colors = list(cluster = annotations$colors$cluster))

    # Same hclust calculation that maxpear does prior to choosing optimal cut point
    hclust.comp <- hclust(as.dist(1 - bigpsm), method = "complete")
    message("Creating PSM heatmap.")
    psm_heatmap = pheatmap(bigpsm,
                           show_rownames = TRUE,
                           show_colnames = FALSE,
                           cluster_rows = hclust.comp,
                           cluster_cols = hclust.comp,
                           annotation_names_row = FALSE,
                           treeheight_col=0,
                           fontsize_row=8,
                           annotation_row = annotations$ann,
                           annotation_col = annotations$ann,
                           color=colorRampPalette((RColorBrewer::brewer.pal(n = 7,
                                                                            name = "Blues")))(100),
                           annotation_colors = annotations$colors,
                           width=14,
                           height=7,
                           filename=paste0("../plots/", exp, "_psm_heatmap.png"))

    generate_balanced_colours <- function(obj) {
        paletteLength = 100
        customColours <- colorRampPalette(c("firebrick3", "white", "navy"))(paletteLength)
        customBreaks <- c(seq(min(obj), 0, length.out=ceiling(paletteLength/2) + 1),
                          seq(max(obj)/paletteLength, max(obj), length.out=floor(paletteLength/2)))
        return(list("colours"=customColours, "breaks"=customBreaks))
    }

    customColours = generate_balanced_colours(obsData)
    message("Creating obsData (raw projections) heatmap.")
    obs_heatmap = pheatmap(obsData,
                           clustering_method="complete",
                           cluster_rows = hclust.comp,
                           annotation_colors = annotations$colors,
                           annotation_row = annotations$ann,
                           cluster_col = FALSE,
                           color = customColours$colours,
                           fontsize_col = 8,
                           fontsize_row = 8,
                           width=6,
                           height = 7,
                           breaks = customColours$breaks,
                           filename=paste0("../plots/", exp, "_obs_heatmap.png"))

    if(grepl("novar", exp)) {
        min_val = min(obsVars)
        max_val = max(obsVars)
        customColours = generate_balanced_colours(c(min_val / 10, max_val * 10))
    } else {
        customColours = generate_balanced_colours(obsVars)
    }
    message("Creating obsVars (raw projection variances) heatmap.")
    obs_vars_heatmap = pheatmap(obsVars,
                                clustering_method="complete",
                                cluster_rows = hclust.comp,
                                cluster_cols = FALSE,
                                annotation_colors = annotations$colors,
                                annotation_row = annotations$ann,
                                color = customColours$colours,
                                fontsize_col = 8,
                                fontsize_row = 8,
                                width=6,
                                height = 7,
                                breaks = customColours$breaks,
                                filename=paste0("../plots/", exp, "_obs_vars_heatmap.png"))

    print("Entries from different PSMs")
    print(sapply(psms, function(x) x[4, 5]))
    message("Creating PSM heatmaps per seed.")
    heatmaps = lapply(psms, function(x) pheatmap(x,
                                                 legend=FALSE,
                                                 color=colorRampPalette((RColorBrewer::brewer.pal(n = 7,
                                                                                                  name = "Blues")))(100),
                                                 cluster_rows = hclust.comp,
                                                 cluster_cols = hclust.comp,
                                                 treeheight_col=0,
                                                 border_color = NA,
                                                 treeheight_row=0))
    for (i in seq_along(datasets)) {
        ggsave(paste0("../plots/", exp, "_psm_heatmap_seed_", i, ".png"), heatmaps[[i]])
    }
    heatmap_grid = grid.arrange(grobs=lapply(heatmaps, function(x) x[[4]]))
    message("Creating PSM heatmaps per seed, but on a grid.")
    ggsave(paste0("../plots/", exp, "_psm_heatmap_grid.png"), heatmap_grid)

    customColours = generate_balanced_colours(bigpsm - psms[[1]])
    message("Creating difference PSMs per seed (Overall - seed PSM).")
    diff_heatmaps = lapply(psms, function(x) pheatmap(bigpsm - x,
                                                 legend=FALSE,
                                                 show_colnames = FALSE,
                                                 show_rownames = FALSE,
                                                 color=customColours$colours,
                                                 breaks=customColours$breaks,
                                                 cluster_rows = hclust.comp,
                                                 cluster_cols = hclust.comp,
                                                 treeheight_col=0,
                                                 border_color = NA,
                                                 treeheight_row=0)[[4]])
    colours = customColours$colours[c(100, 76, 51, 26, 1)]
    breaks = sapply(customColours$breaks[c(101, 76, 51, 26, 1)], function(x) signif(x, 3))
    dummy_plot = ggplot(data.frame(x=as.factor(breaks))) +
        scale_color_manual(values=colours, breaks=breaks) +
        labs(color="Overall PSM - Seed PSM") +
        geom_point(aes(x, x, color=x))
    dummy_grobs <- ggplot_gtable(ggplot_build(dummy_plot))
    leg_index <- which(sapply(dummy_grobs$grobs, function(x) x$name) == "guide-box")
    legend = dummy_grobs$grobs[[leg_index]]

    print(class(legend))
    diff_heatmaps[[length(diff_heatmaps) + 1]] = legend
    heatmap_grid = grid.arrange(grobs=diff_heatmaps)
    ggsave(paste0("../plots/", exp, "_psm_heatmap_diff_grid.png"), heatmap_grid)

    individual_calls = do.call(cbind, lapply(psms, function(x) adjust_labels_B_to_match_A(calls$cl, maxpear(x, method="comp")$cl)))
    all_calls = data.frame(cbind(calls$cl, individual_calls))
    colnames(all_calls) = c("Overall", paste("Seed", seq_along(datasets)))

    print(all_calls)

    map_to_call_counts <- function(calls) {
        call_labels = paste0(LETTERS[1:max(calls)], " (", table(calls), ")")
        return(plyr::mapvalues(calls, from=1:max(calls), to=call_labels))
    }
    cell_labels = apply(all_calls, MARGIN=2, map_to_call_counts)

    message("Creating calls heatmap.")
    calls_heatmap = pheatmap(all_calls,
                             display_numbers = cell_labels,
                             fontsize_number = 4,
                             number_color = "#CCCCCC",
                             show_rownames = TRUE,
                             show_colnames = TRUE,
                             color = palette[1:max(all_calls)],
                             breaks = seq(0.5, max(all_calls) + 0.5, by=1),
                             cluster_col = FALSE,
                             cluster_row = hclust.comp,
                             fontsize_row = 8,
                             width=9,
                             height=7,
                             filename=paste0("../plots/", exp, "_calls_heatmap.png"))
    generate_balanced_colours <- function(obj) {
        paletteLength = 100
        customColours <- colorRampPalette(c("firebrick3", "white", "navy"))(paletteLength)
        customBreaks <- c(seq(min(obj), 0, length.out=ceiling(paletteLength/2) + 1),
                          seq(max(obj)/paletteLength, max(obj), length.out=floor(paletteLength/2)))
        return(list("colours"=customColours, "breaks"=customBreaks))
    }
    message("Done, all good!")
    if(!is.null(update_label)){
        return(list(name=exp, bigpsm=bigpsm, calls=calls, hclust.comp=hclust.comp, mclust_calls=mclust_solution$classification))
    } else {
        return(list(name=exp, bigpsm=bigpsm, calls=calls, hclust.comp=hclust.comp, mclust_calls=mclust_solution$classification, Trait_ID = Trait_ids)) # In case we change labelling, we can use Trait_ids to match datasets
    }
    
}


## Define experiments and generate PSM plots

#exp = c("PC58_ALL", "PC58_IMD", "PC58_BMK", "PCgrac_ALL", "PCgrac_IMD", "PCgrac_BMK")
#exp = c("UKBB_PC58_ALL", "UKBB_PC58_IMD", "UKBB_PC58_BMK", "UKBB_PCgrac_ALL", "UKBB_PCgrac_IMD", "UKBB_PCgrac_BMK")
#exp = c("UKBB_PCgrac_ALL", "UKBB_PCgrac_IMD", "UKBB_PCgrac_BMK")
#exp = c("UKBB_PC58_ALL", "UKBB_PC58_IMD", "UKBB_PC58_BMK")
#exp="UKBB_PC258_IMD"

#exp="UKBB_PC2_IMD"

exp=c("UKBBsig_PC2_IMD", "UKBBsig_PC58_IMD")


labeldict  <- fread("../../Projections/PT_sig_202211_UKBB.tsv") %>% .[,.(Trait, Label)] %>% unique

lapply(exp, function(x){
    result <- psm_plots(x, burnin = 250000, update_label = labeldict)
    saveRDS(result, paste0("../psm_results/", x, "_psm_data.rds"))
})



