# Figure 3b

## 代码与数据

| 文件 | 说明 |
|---|---|
| `plot.R` | Figure 3b 作图代码。 |
| `plot_data.csv` | 作图数据。六类脂质共 131 个分子，每个分子一行，包含 AAD vs HP 和 AAD vs AMI 的 log₂FC 与校正后 q 值。 |
| `Fig3b_editable.svg` | 可编辑的参考图。 |

横条表示每类实测分子中，在两次比较中均满足 log₂FC < 0 且 q < 0.05 的分子比例。q 值沿用原分析对 2,330 项检验（466 个脂质分子 × 五组两两比较）的 Benjamini–Hochberg 校正结果。

## 运行环境

已用 R 4.6.1 验证，无需额外 R 包。

## 运行方法

将 R 脚本与输入 CSV 放在同一文件夹，在该目录打开终端并运行：

```sh
Rscript --vanilla plot.R
```

也可将 R 或 RStudio 的工作目录设为该文件夹并运行：

```r
source("plot.R", encoding = "UTF-8")
plot_fig3b()
```

## 输出

- `Fig3b_editable.svg`

SVG 保存到工作目录。重新运行会覆盖同名文件。
