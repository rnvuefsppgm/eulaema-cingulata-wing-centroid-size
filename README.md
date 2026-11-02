# Wing centroid size of *Eulaema cingulata* in Atlantic Forest remnants: data and analysis scripts

Data and scripts that reproduce the analyses, supplementary tables and figures of the manuscript
**"Wing centroid size variation among *Eulaema cingulata* (Fabricius) populations from Atlantic Forest remnants of different sizes"**
(submitted to the *Journal of Insect Conservation*). They are deposited here at the request of the reviewers; the manuscript
cites this repository for the data (folder `data/`) and the scripts (folder `scripts/`).

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

## Data (Supplementary Data S1a–e of the article)

The folder `data/` holds the individual-level data analysed in the article, as CSV files (UTF-8, comma-separated; CC BY 4.0) with the same values as the
Word files Data S1a–e of the supplementary material:

| File | Content | Article item |
|---|---|---|
| `data/DataS1a_field_dictionary.csv` | Description of every field of the data set | Data S1a |
| `data/DataS1b_individual_male_records.csv` | One row per male (211): specimen code, field label, fragment (A–I), wing centroid size (image units), date of capture, season, year and fragment attributes (area, isolation, elevation, disturbance) | Data S1b |
| `data/DataS1c_raw_wing_landmark_coordinates.csv` | Raw x–y coordinates (image units) of the 12 wing landmarks of each male, as digitised | Data S1c |
| `data/DataS1d_procrustes_aligned_landmark_coordinates.csv` | Landmark coordinates after generalised Procrustes analysis (unit centroid size), input of the shape analyses | Data S1d |
| `data/DataS1e_geographic_distances_between_fragments_km.csv` | Straight-line distances (km) between the nine forest fragments | Data S1e |

Note: the column `Bait (reference allocation)` of `DataS1b` is **not an observed variable**. It is a random, balanced assignment of the five
bait compounds within each fragment (fixed seed 2025) that reproduces Table S1; the bait of capture of each male was not recorded and this
column is not used in any analysis.

## Source spreadsheets (inputs of the scripts)

The scripts start from the original field and digitising spreadsheets, which are not distributed here; they are available from the
corresponding author on reasonable request, as are the raw wing images. Expected files in `data/raw/`:

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

- Code (`scripts/`): MIT License (see `LICENSE`).
- Data (`data/`): Creative Commons Attribution 4.0 International (CC BY 4.0) (see `data/LICENSE.md`).
