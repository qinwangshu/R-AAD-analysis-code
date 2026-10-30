# Figure 4a

## Code and data

| File | Description |
|---|---|
| `Fig4a_plot.R` | R code for Figure 4a. |
| `Fig4a_data.csv` | Input data. Group labels and four metabolite abundances for 152 participants. |
| `Fig4a_editable.svg` | Reference figure. |

## Requirements

R 4.6.1, brglm2 1.1.1 and pROC 1.19.1. Install the required packages in R:

```r
install.packages(c("brglm2", "pROC"), repos = "https://cloud.r-project.org")
```

## Run

Keep the R script and input CSV together. Open a terminal in this folder and run:

```text
Rscript --vanilla Fig4a_plot.R
```

Alternatively, set the R/RStudio working directory to this folder and run:

```r
source("Fig4a_plot.R", encoding = "UTF-8")
```

## Outputs

- `Fig4a_editable.svg`
- `Fig4a_summary.csv`
- `Fig4a_predictions.csv`
- `Fig4a_ROC_coordinates.csv`

Outputs are saved beside the input CSV. Existing outputs with the same names are overwritten.
