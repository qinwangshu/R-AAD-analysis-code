# Figure 3c

## Code and data

| File | Description |
|---|---|
| `plot.R` | R code for the boxplots and comparison labels. |
| `plot_data.csv` | Participant-level relative abundances for seven LPE species. Columns N150–N156 identify the species. The table includes 152 participants across CTRL, HP, AMI and AAD. |
| `comparison_results.csv` | Original summary statistics and adjusted q values for seven LPE species across five pairwise comparisons. |
| `Fig3c_boxplots_editable.svg` | Editable reference figure. |

Each point represents one participant. Boxes span the interquartile range, center lines indicate medians, and whiskers extend to the most extreme observed values within 1.5 times the interquartile range of the box edges. The y-axis scale varies by species.

Brackets compare AAD with each of CTRL, HP and AMI. An asterisk indicates an adjusted q value < 0.05, and `ns` indicates an adjusted q value ≥ 0.05. The script uses the original Benjamini–Hochberg adjustment across 2,330 tests (466 species across five pairwise comparisons).

## Requirements

Tested with R 4.6.1. No additional R packages are required.

## Run

Keep the R script and the two CSV files together. Open a terminal in this folder and run:

```sh
Rscript --vanilla plot.R
```

Alternatively, set the R or RStudio working directory to this folder and run:

```r
source("plot.R", encoding = "UTF-8")
```

## Outputs

- `Fig3c_boxplots_editable.svg`

The SVG is written to the working directory. Running the script again overwrites the existing file of the same name.
