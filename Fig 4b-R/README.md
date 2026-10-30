# Figure 4b

## Code and data

| File | Description |
|---|---|
| `Fig4b_plot.R` | R code for Figure 4b. |
| `Fig4b_data.csv` | Input data. Group labels and four metabolite abundances for 313 participants. |
| `Fig4b_editable.svg` | Reference figure. |

## Requirements

R 4.6.1, brglm2 1.1.1 and pROC 1.19.1. Install the required packages in R:

```r
install.packages(c("brglm2", "pROC"), repos = "https://cloud.r-project.org")
```

## Run

Keep the R script and input CSV together. Open a terminal in this folder and run:

```text
Rscript --vanilla Fig4b_plot.R
```

Alternatively, set the R/RStudio working directory to this folder and run:

```r
source("Fig4b_plot.R", encoding = "UTF-8")
```

## Outputs

- `Fig4b_editable.svg`
- `Fig4b_summary.csv`
- `Fig4b_predictions.csv`
- `Fig4b_ROC_coordinates.csv`

Outputs are saved beside the input CSV. Existing outputs with the same names are overwritten.
