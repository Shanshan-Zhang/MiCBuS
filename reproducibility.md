# Reproducibility Information for MiCBuS

This document captures the full software environment used to render the
MiCBuS vignette (`Intro_to_MiCBuS.Rmd`). Including this information ensures that
analyses performed using MiCBuS can be reproduced accurately by other users.

---

## R Environment

```
R version 4.3.1 (2023-06-16 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19045)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.utf8
[2] LC_CTYPE=English_United States.utf8
[3] LC_MONETARY=English_United States.utf8
[4] LC_NUMERIC=C
[5] LC_TIME=English_United States.utf8

time zone: America/New_York
tzcode source: internal
```

---

## Attached Packages During Vignette Execution

```
Biobase_2.60.0
BiocGenerics_0.46.0
knitr_1.46
ggVennDiagram_1.5.2
ggplot2_3.5.1
dplyr_1.1.4
MiCBuS_0.1.0
```

---

## Loaded (Namespace Only) Packages

These packages were loaded automatically by dependencies during vignette execution:

```
RColorBrewer_1.1-3 rstudioapi_0.16.0
jsonlite_1.8.8 magrittr_2.0.3
spatstat.utils_3.2-0 farver_2.1.1
rmarkdown_2.26 zlibbioc_1.46.0
vctrs_0.6.5.9000 ROCR_1.0-11
spatstat.explore_3.2-7 RCurl_1.98-1.14
S4Arrays_1.0.6 htmltools_0.5.8.1
sass_0.4.9 sctransform_0.4.1
parallelly_1.37.1 KernSmooth_2.23-22
bslib_0.7.0 htmlwidgets_1.6.4
ica_1.0-3 plyr_1.8.9
plotly_4.10.4 zoo_1.8-12
cachem_1.0.8 alabama_2023.1.0
igraph_2.0.3 mime_0.12
lifecycle_1.0.4 pkgconfig_2.0.3
Matrix_1.6-4 R6_2.5.1
fastmap_1.1.1 GenomeInfoDbData_1.2.10
MatrixGenerics_1.12.3 fitdistrplus_1.1-11
future_1.33.2 shiny_1.8.1.1
digest_0.6.35 numDeriv_2016.8-1.1
colorspace_2.1-0 S4Vectors_0.38.2
patchwork_1.2.0 DESeq2_1.40.2
Seurat_5.0.3 tensor_1.5
RSpectra_0.16-1 irlba_2.3.5.1
GenomicRanges_1.52.1 labeling_0.4.3
progressr_0.14.0 fansi_1.0.6
spatstat.sparse_3.0-3 httr_1.4.7
polyclip_1.10-6 abind_1.4-5
compiler_4.3.1 withr_3.0.0
BiocParallel_1.34.2 fastDummies_1.7.3
highr_0.10 MASS_7.3-60
DelayedArray_0.26.7 proxyC_0.4.1
tools_4.3.1 lmtest_0.9-40
httpuv_1.6.15 future.apply_1.11.2
goftest_1.2-3 glue_1.6.2
nlme_3.1-162 promises_1.3.0
grid_4.3.1 Rtsne_0.17
cluster_2.1.4 reshape2_1.4.4
generics_0.1.3 gtable_0.3.5
spatstat.data_3.0-4 tidyr_1.3.1
data.table_1.14.8 XVector_0.40.0
sp_2.1-3 utf8_1.2.4
spatstat.geom_3.2-9 RcppAnnoy_0.0.22
ggrepel_0.9.4 RANN_2.6.1
pillar_1.9.0 stringr_1.5.1
spam_2.10-0 RcppHNSW_0.6.0
limma_3.56.2 later_1.3.2
splines_4.3.1 lattice_0.21-8
survival_3.5-5 deldir_2.0-4
tidyselect_1.2.1 locfit_1.5-9.9
miniUI_0.1.1.1 pbapply_1.7-2
gridExtra_2.3 IRanges_2.34.1
SummarizedExperiment_1.30.2 scattermore_1.2
stats4_4.3.1 xfun_0.43
SimBu_1.2.0 matrixStats_1.0.0
stringi_1.7.12 lazyeval_0.2.2
yaml_2.3.8 evaluate_0.23
codetools_0.2-19 SECRET_0.0.0.9000
tibble_3.2.1 cli_3.6.1
uwot_0.2.2 xtable_1.8-4
reticulate_1.36.1 munsell_0.5.1
jquerylib_0.1.4 GenomeInfoDb_1.36.4
Rcpp_1.0.11 globals_0.16.3
spatstat.random_3.2-3 png_0.1-8
parallel_4.3.1 dotCall64_1.1-1
bitops_1.0-7 sparseMatrixStats_1.12.2
listenv_0.9.1 viridisLite_0.4.2
scales_1.3.0 ggridges_0.5.6
crayon_1.5.2 SeuratObject_5.0.1
leiden_0.4.3.1 purrr_1.0.2
rlang_1.1.1 cowplot_1.1.3
```

---

## Key Dependencies in the MiCBuS Workflow

- **Seurat 5.0.3**
- **SECRET 0.0.0.9000**, **alabama 2023.1.0**
- **SimBu 1.2.0**
- **DESeq2 1.40.2**
- **Biobase 2.60.0**, **SummarizedExperiment 1.30.2**, **Matrix 1.6-4**

---

## Recommendations for Reproducible Use

Use R version 4.3.1, match package versions where possible, save your own
sessionInfo(), and cite this reproducibility file when publishing.

---

## Contact

Author: Shanshan Zhang  
GitHub: https://github.com/Shanshan-Zhang/MiCBuS
