# Figure 3b

## Code and data

| File | Description |
|---|---|
| `plot.R` | R code for Figure 3b. |
| `plot_data.csv` | Input data. One row per lipid species (131 species across six classes), with log₂ fold changes and adjusted q values for AAD vs HP and AAD vs AMI. |
| `Fig3b_editable.svg` | Editable reference figure. |

Each bar shows the proportion of measured species in a lipid class with a negative log₂ fold change and an adjusted q value < 0.05 in both comparisons. The q values come from the original Benjamini–Hochberg correction across 2,330 tests (466 species across five pairwise comparisons).

## Requirements

Tested with R 4.6.1. No additional R packages are required.

## Run

Keep the R script and input CSV together. Open a terminal in this folder and run:

```sh
Rscript --vanilla plot.R
```

Alternatively, set the R or RStudio working directory to this folder and run:

```r
source("plot.R", encoding = "UTF-8")
plot_fig3b()
```

## Outputs

- `Fig3b_editable.svg`

The SVG is written to the working directory. Running the script again overwrites the existing file of the same name.
