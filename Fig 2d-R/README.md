# Figure 2d

## Code and data

| File | Description |
|---|---|
| `Fig2d_77_heatmap.R` | R code for Figure 2d. |
| `Fig2d_77_data.csv` | Input data. Feature annotations and abundances for 77 features across 152 samples. |
| `Fig2d_77_editable.svg` | Reference figure. |

## Requirements

R 4.6.1; no additional R packages are required.

## Run

Keep the R script and input CSV together. Open a terminal in this folder and run:

```text
Rscript --vanilla Fig2d_77_heatmap.R
```

Alternatively, set the R/RStudio working directory to this folder and run:

```r
source("Fig2d_77_heatmap.R", encoding = "UTF-8")
```

## Outputs

- `Fig2d_77_editable.svg`

Outputs are saved beside the input CSV. Existing outputs with the same names are overwritten.

To additionally export `Fig2d_77.pdf` or `Fig2d_77_600dpi.png`, set `EXPORT_PDF` or `EXPORT_PNG` to `TRUE` near the top of the script. PNG export requires Cairo support.
