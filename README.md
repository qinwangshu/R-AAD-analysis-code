[README.md](https://github.com/user-attachments/files/34181500/README.md)
# R code and data for the figures

This directory contains the R scripts and CSV data for Figures 2c, 2d, 3a–c, 4a–b and Supplementary Figures 3a–c. Each figure folder includes a README with running instructions and output filenames.

## Code and input data

| Figure | Folder | R script | Input CSV |
|---|---|---|---|
| 2c | `Fig 2c-R` | `Fig2c_plot.R` | `Fig2c_all_pathways.csv` |
| 2d | `Fig 2d-R` | `Fig2d_77_heatmap.R` | `Fig2d_77_data.csv` |
| 3a | `Fig 3a-R` | `plot_map.R` | `plot_data.csv` |
| 3b | `Fig 3b-R` | `plot.R` | `plot_data.csv` |
| 3c | `Fig 3c-R` | `plot.R` | `plot_data.csv` |
| 4a | `Fig 4a-R` | `Fig4a_plot.R` | `Fig4a_data.csv` |
| 4b | `Fig 4b-R` | `Fig4b_plot.R` | `Fig4b_data.csv` |
| S3a | `Fig s3a-R` | `FigS3a_plot.R` | `FigS3a_data.csv` |
| S3b | `Fig s3b-R` | `FigS3b_plot.R` | `FigS3b_data.csv` |
| S3c | `Fig s3c-R` | `FigS3c_plot.R` | `FigS3c_data.csv` |

## Requirements

R 4.6.1. Figures 2c, 2d, 3a–c, S3a and S3b use packages included with R. Figures 4a and 4b require `brglm2` and `pROC`; Figure S3c requires `Mfuzz` and its dependencies.

Install the additional packages in R:

```r
install.packages(c("brglm2", "pROC", "BiocManager"),
                 repos = "https://cloud.r-project.org")
BiocManager::install("Mfuzz", ask = FALSE, update = FALSE)
```

Figure 2c requires Cairo support for SVG/PNG export. Cairo is also needed for optional PNG export in Figures 2d and S3a. Check support with `capabilities("cairo")`.

## Run a figure

Keep each script and its input CSV in the same folder. Open a terminal in that figure folder and run its script. For example, for Figure 2c:

```text
Rscript --vanilla Fig2c_plot.R
```

Alternatively, set the R/RStudio working directory to the figure folder and run:

```r
source("Fig2c_plot.R", encoding = "UTF-8")
```

For other figures, substitute the script name listed above.

## Outputs

Outputs are saved beside the input CSV. Each script generates an SVG; some also export numerical results, PDF or PNG files. Figure 3a generates `Fig3A_R.svg`; Figures 3b and 3c each generate `figure_R.svg` in their respective folders. Existing outputs with the same names are overwritten. See each folder's README for the complete output list and optional export settings.
