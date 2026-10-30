# Figure 3c

## Code and data

| File | Description |
|---|---|
| `plot.R` | R code for Figure 3c. |
| `plot_data.csv` | Input data. Group means, fold changes and q values for seven LPE features across five comparisons. |
| `original_LPE_rows.csv` | Source-row index for the LPE features; not required for plotting. |
| `figure.svg` | Reference figure. |

## Requirements

R 4.6.1; no additional R packages are required.

## Run

Keep the R script and input CSV together. Open a terminal in this folder and run:

```text
Rscript --vanilla plot.R
```

Alternatively, set the R/RStudio working directory to this folder and run:

```r
source("plot.R", encoding = "UTF-8")
```

## Outputs

- `figure_R.svg`

Outputs are saved beside the input CSV. Existing outputs with the same names are overwritten.
