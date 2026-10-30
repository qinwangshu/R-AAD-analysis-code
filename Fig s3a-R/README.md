# Figure S3a

## Code and data

| File | Description |
|---|---|
| `FigS3a_plot.R` | R code for Figure S3a. |
| `FigS3a_data.csv` | Input data. Feature membership in the three sets used in the Venn diagram. |
| `FigS3a_editable.svg` | Reference figure. |

## Requirements

R 4.6.1; no additional R packages are required.

## Run

Keep the R script and input CSV together. Open a terminal in this folder and run:

```text
Rscript --vanilla FigS3a_plot.R
```

Alternatively, set the R/RStudio working directory to this folder and run:

```r
source("FigS3a_plot.R", encoding = "UTF-8")
```

## Outputs

- `FigS3a_editable.svg`

Outputs are saved beside the input CSV. Existing outputs with the same names are overwritten.

To additionally export `FigS3a.pdf` or `FigS3a_600dpi.png`, set `EXPORT_PDF` or `EXPORT_PNG` to `TRUE` near the top of the script. PNG export requires Cairo support.
