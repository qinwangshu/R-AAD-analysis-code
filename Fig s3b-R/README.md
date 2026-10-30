# Figure S3b

## Code and data

| File | Description |
|---|---|
| `FigS3b_plot.R` | R code for Figure S3b. |
| `FigS3b_data.csv` | Input data. Feature annotations and abundances for 81 features across 152 samples. |
| `FigS3b_editable.svg` | Reference figure. |

## Requirements

R 4.6.1; no additional R packages are required.

## Run

Keep the R script and input CSV together. Open a terminal in this folder and run:

```text
Rscript --vanilla FigS3b_plot.R
```

Alternatively, set the R/RStudio working directory to this folder and run:

```r
source("FigS3b_plot.R", encoding = "UTF-8")
```

## Outputs

- `FigS3b_editable.svg`

Outputs are saved beside the input CSV. Existing outputs with the same names are overwritten.
