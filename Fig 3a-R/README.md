# Figure 3a

## Code and data

| File | Description |
|---|---|
| `plot_map.R` | R code for Figure 3a. |
| `plot_data.csv` | Input data. Counts and plotting coordinates for 21 lipid classes. |
| `Fig3A_reference_format_editable.svg` | Reference figure. |

## Requirements

R 4.6.1; no additional R packages are required.

## Run

Keep the R script and input CSV together. Open a terminal in this folder and run:

```text
Rscript --vanilla plot_map.R
```

Alternatively, set the R/RStudio working directory to this folder and run:

```r
source("plot_map.R", encoding = "UTF-8")
```

## Outputs

- `Fig3A_R.svg`

Outputs are saved beside the input CSV. Existing outputs with the same names are overwritten.
