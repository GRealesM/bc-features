## Running DPMUnc in Blood cell basis 3


Guillermo Reales

2023/01/18

----

After visualising the significant datasets across different components, we want to cluster diseases using their projection values across multiple components.

For that we'll use shiny new DPMUnc method, with defaults.

We'll use a set of overall significant projections (FDR < 1%), including PanUKBB (not Neale or FinnGen). From this, we'll focus on IMD, and we'll manually remove redundancies prior to clustering. 

We'll cluster the following sets of PCs:

-   **PC58sig**: PC5 & PC8 (PC58), since those PCs accumulate the largest numbers of significant IMD and biomarkers datasets.
-   **PC2sig**: PC2, as it seem to capture signals from neutrophils.

We followed these steps for dataset selection and prior processing:


1.  Select traits that were significant overall at FDR 1%.
2.  Remove all datasets from the blood cell category.
3.  From those, select traits that are FDR.PC \< 1% for either PC of interest (eg. PC5 and PC8, or PC2). See [code](https://github.com/GRealesM/Bases/blob/master/cell_basis_v3_varimax/DPMUnc/code/01-Prepare_data.R) for details.
4.  Remove obvious redundant datasets (manually).
5.  Cluster using IMD traits only.
6.  Apply ellipse assignment (see https://github.com/GRealesM/Bases/blob/master/cell_basis_v3_varimax/DPMUnc/ellipses/biomarkers-to-nearest-cluster-guille.R) on biomarkers for PC5+8.

Then, we ran DPMUnc using delta and var.delta for selected PCs and final selected datasets.

We ran DPMUnc with 5,000,000 iterations, sampled 1/10, with 5 seeds (1-5), and default parameters: kappa0 = 0.01, alpha0 = 2, beta0 = 0.1. Code [here](https://github.com/GRealesM/Bases/blob/master/cell_basis_v3_varimax/DPMUnc/code/02-Running_DPMUnc.R).


Then we created **traceplots** for diagnostics, and **PSM plots** to visualise the clustering.

Please note, if you're reusing the code, that we used seeds 1-5, and plotting scripts might not function properly if using other seeds.


## What's in the box

* `code/` contains all the necessary code to run DPMUnc.
  - `01a-Prepare_data_Aug2022.R` will prepare the projection data from the projection table of significant overall traits and the quantOR transformed. We'll select the PCs we're interested in and keep datasets that are significant for at least one of those PCs from the projection table. We'll also select specific trait classes (eg. IMD and biomarkers) for clustering. This will use FinnGen (not UKBB) data.
  - `01b-Prepare_data_Nov2022.R` will prepare data similarly to `01a` but using PanUKBB projections (not FinnGen) instead.
  - `02-Running_DPMUnc.R` will run DPMUnc with the defaults above. Adapted to be run in parallel using array jobs (one job per seed). We used `slurm_02_RunDPMUNC_*` to that end.
  - `03-Create_traceplots.R` will use the run results in `../results/` to plot traceplots for each seed, and plot them combined. Results from this step in `../plots/`.
  - `04-Create_PSM_heatmaps.R` will remove burn-in iterations (default: first half) and create multiple PSM plots, including general PSM plots and difference PSM plots per seed, using adapted Kath's functions.
  - `04-Create_custom_PSM_heatmaps.R` will remove burn-in iterations (default: first half) and create multiple PSM plots, including general PSM plots and difference PSM plots per seed, using adapted Kath's functions.

* `data/` contains the inputs for DPMUnc, in the format `{exp}_{subset}_Delta.tsv` for projections and `{exp}_{subset}_Var.tsv` for variance, where `exp` is the name we gave the set of PCs (eg. PC58), `subset` refers to the selected datasets (ie. All, IMD only, or BMK only).
  - `tmp_redundancies.tsv` file containing all datasets in both PC58 and PCgrac, to be inspected and selected redundant datasets to remove.
  - `tmp_processed_redundancies.tsv` contains a selection of datasets to remove (Remove = 1).
* `plots/` contains all plots generated from DPMUnc run:
  - `{exp}_{subset}_calls_heatmap.png`: Calls heatmap.
  - `{exp}_{subset}_obs_heatmap.png`: Raw projection heatmap.
  - `{exp}_{subset}_obs_vars_heatmap.png`: Raw variance heatmap.
  - `{exp}_{subset}_psm_heatmap_diff_grid.png`: Difference heatmap (average - individual seed) for each seed, in a grid.
  - `{exp}_{subset}_psm_heatmap_grid.png`: PSM heatmap for each seed, in a grid.
  - `{exp}_{subset}_psm_heatmap_seed_{seed}.png`: PSM heatmap for each seed.
  - `{exp}_{subset}_psm_heatmap.png`: **Main PSM heatmap**.
  - `{exp}_{subset}_trace_medians.png`: Traceplots using medians for each selected block.
  - `{exp}_{subset}_trace_quantiles.png`: Traceplots using quantiles for each selected block.
  - `DPMUnc_Report_202208.Rmd`: A report in Rmarkdown format, with all plots above. See it online [here](https://grealesm.github.io/Bases/cell_basis_v3_varimax/DPMUnc/plots/DPMUnc_Report_202208.html).
* `psm_results/` contains results from the PSM plotting stage, for reproducibility.
* `results/` contains the results from DPMUnc, in the format `{exp}_{subset}_{seed}`. Only available at the HPC, for space reasons.




