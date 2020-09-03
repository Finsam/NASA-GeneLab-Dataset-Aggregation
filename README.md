# Dataset Correction for NASA GeneLab Datasets

Scripts utilized in correction of combined RNA-seq datasets from the NASA GeneLab repository.

**Files and Descriptions:**

- [ComBat_save_custom.R](https://github.com/Finsam/NASA-GeneLab-Dataset-Aggregation/blob/master/ComBat_save_custom.R)
  - Modified version of ComBat from the [Surrogate Variable Analysis](https://bioconductor.org/packages/release/bioc/html/sva.html) package. Allows for extraction of correction factors from ComBat correction (with the `ComBat_save` method), as well as application of custom correction factors to new datasets with correspoding groups (with the `ComBat_custom` method).
  
- [ComBat_seq_save_custom.R](https://github.com/Finsam/NASA-GeneLab-Dataset-Aggregation/blob/master/ComBat_seq_save_custom.R)
  - Modified version of ComBat_seq from the [Surrogate Variable Analysis](https://bioconductor.org/packages/release/bioc/html/sva.html) package. Allows for extraction of correction factors from ComBat_seq correction (with the `ComBat_seq_save` method), as well as application of custom correction factors to new datasets with correspoding groups (with the `ComBat_seq_custom` method).

- [apply_corrections.R](https://github.com/Finsam/NASA-GeneLab-Dataset-Aggregation/blob/master/apply_corrections.R)
  - Script utilized to apply various correction methods to raw RNA transcript counts files and exports corrected counts.

- [batchQC_liver.R](https://github.com/Finsam/NASA-GeneLab-Dataset-Aggregation/blob/master/batchQC_liver.R)
  - Script utilized to apply batchQC diagnostics to corrected and uncorrected datasets.

- [dge_filter.R](https://github.com/Finsam/NASA-GeneLab-Dataset-Aggregation/blob/master/dge_filter.R)
  - Script utilized to filter differential gene expression (DGE) files for corrected and uncorrected datasets. The DGE files were filtered to include flight vs. ground control sample comparisons.

- [mouse_liver_pca.R](https://github.com/Finsam/NASA-GeneLab-Dataset-Aggregation/blob/master/mouse_liver_pca.R)
 - Script utilized to conduct and plot PCA analyses of corrected and uncorrected datasets.
