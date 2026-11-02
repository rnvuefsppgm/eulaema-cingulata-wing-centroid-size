# Wing centroid size of *Eulaema cingulata* in Atlantic Forest remnants: analysis scripts

Scripts that reproduce the analyses, supplementary tables and figures of the manuscript
**"Wing centroid size variation among *Eulaema cingulata* (Fabricius) populations from Atlantic Forest remnants of different sizes"**
(submitted to the *Journal of Insect Conservation*). They are deposited here at the request of the reviewers; the manuscript
cites this repository in its Data analysis, Data availability and Code availability sections.

## Scripts (run in this order)

| Script | Task |
|---|---|
| `scripts/01_build_individual_dataset.py` | Merges the source spreadsheets (landmark coordinates, centroid size, capture dates) into one table per individual. |
| `scripts/02_analyse_centroid_size_and_wing_shape.py` | Cleans the data set and runs all main analyses: centroid size tests (Kruskal–Wallis, pairwise post hoc tests), regressions and mixed models with assumption checks, Procrustes superimposition, shape PCA, MANOVA/CVA, Mahalanobis distances and Mantel tests. |
| `scripts/03_test_pairwise_shape_differences.py` | Pairwise permutation tests of wing-shape differences between remnants (Holm correction). |
| `scripts/04_plot_figures.py` | Draws the main and supplementary figures (PNG 600 dpi and editable SVG). |
| `scripts/05_crosscheck_models_in_R.R` | Independent re-fit of the main models in R (lme4, vegan, MASS) to verify the Python results. |

## Usage

```bash
pip install -r requirements.txt
python scripts/01_build_individual_dataset.py
python scripts/02_analyse_centroid_size_and_wing_shape.py
python scripts/03_test_pairwise_shape_differences.py
python scripts/04_plot_figures.py
Rscript scripts/05_crosscheck_models_in_R.R .
```

Each script accepts an optional project root as its first argument (default: the repository folder).
Inputs are read from `data/raw/` and all results are written to `outputs/` (tables as CSV, `results.json`, figures in `outputs/figures/`).

## Input data

The source spreadsheets are not distributed in this repository. The individual-level data used in the analyses
(field dictionary, individual records, raw and Procrustes-aligned landmark coordinates, and geographic distances
between remnants) are published with the article as Supplementary Data S1a–e. The original spreadsheets and wing
images are available from the corresponding author on reasonable request. Expected files in `data/raw/`:

- `DADOS MORFOMETRIA (1).xlsx` and `DADOS MORFOMETRIA_TODOS_MOD.xlsx` (landmark coordinates)
- `DADOS MORFOMETRIA_TODOS.xlsx` (centroid size)
- `Data_coleta_El_cingulata_ morfometria.xlsx` (capture dates)
- `wing_photo.png` (wing image used to draw the landmark figure)

## Main outputs

| File in `outputs/` | Content |
|---|---|
| `individuals_clean.csv`, `fragments.csv` | Analysed individuals and remnant attributes |
| `cs_descriptive_by_fragment.csv` | Descriptive statistics of centroid size by remnant |
| `cs_pairwise_posthoc_tests.csv` | Pairwise centroid size comparisons |
| `shape_pairwise_hotelling_mahalanobis.csv` | Pairwise wing-shape tests and Mahalanobis distances |
| `capture_records.csv` | Capture records of the analysed males |
| `procrustes_coords.csv`, `shape_PCs.csv`, `cva_scores.csv` | Shape variables |
| `results.json` | All test statistics reported in the manuscript |
| `R_crosscheck_output.txt` | R verification of the main models |

## Software

Python 3.14.6 (numpy 2.5.1, scipy 1.18.0, pandas 3.0.5, statsmodels 0.15.0, matplotlib 3.11.2, openpyxl, Pillow);
R 4.6.0 (lme4 2.0.6, vegan 2.7.5, MASS). Random seeds are fixed in the scripts, so results are reproducible.

## License

Code released under the MIT License (see `LICENSE`).
