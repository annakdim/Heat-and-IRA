[README.md](https://github.com/user-attachments/files/31694802/README.md)
# Heat-and-IRA

Code to reproduce the analyses in:

> Dimitrova, A. *et al.* Extreme heat exposure and intrapartum-related asphyxia in The Gambia and Burkina Faso: a case-crossover study. *Nature Communications* (submitted).

This repository contains the R scripts used to clean the clinical cohort data, link births to ERA5 temperature data, and run the time-stratified case-crossover analysis of heat exposure and intrapartum-related asphyxia (IRA) described in the manuscript.

The individual-level clinical data used in the manuscript are **not included** in this repository because they contain sensitive participant information and are subject to ethical and data-sharing restrictions (see the manuscript's Data Availability statement). To allow the code to be installed and run end-to-end, a small **simulated** demo dataset and demo script are provided in [`demo/`](demo/).

## Repository structure

```
1_Clean the data/                    Clean and merge the three clinical cohorts
2_Generate temperature exposure-ERA5/  Link births to ERA5 temperature data and derive exposure percentiles
3_Analysis case-crossover/           Time-stratified case-crossover analysis (main, healthy-case,
                                      full-cohort, and negative-control models)
demo/                                Simulated dataset + self-contained demo of the analysis pipeline
```

## 1. System requirements

**Software dependencies**

The scripts are written in R and require the following packages (all available from CRAN unless noted):

| Package | Used for |
|---|---|
| `dplyr`, `tidyr`, `tidyverse`, `purrr`, `tibble`, `stringr`, `readr`, `readxl`, `haven`, `zoo`, `lubridate` | Data wrangling |
| `sf`, `sp`, `raster`, `rgdal`, `rworldmap` | Spatial data handling and linkage to ERA5 grids |
| `dlnm`, `survival` | Exposure-lag structures and conditional logistic regression (`clogit`) |
| `ggplot2`, `viridis` | Figures |
| `here` | Reproducible file paths |

> **Note:** `rgdal` and `rworldmap` are retired/archived on CRAN as of late 2023. On a recent R installation, install them from the CRAN archive (e.g. via `remotes::install_version()`) or, if preferred, adapt the affected lines in `1_Clean the data/` and `2_Generate temperature exposure-ERA5/` to use `sf`/`terra` equivalents.

**Operating systems tested:** *Windows 11*

**R version tested:** *R version 4.4.2*

**Package versions tested:** 
*  dlnm: 2.4.10
*  dplyr: 1.1.4
  ggplot2: 4.0.0
  haven: 2.5.4
  here: 1.0.2
  lubridate: 1.9.4
  purrr: 1.2.1
  raster: 3.6.30
  readr: 2.1.5
  readxl: 1.4.3
  rworldmap: 1.3.8
  sf: 1.0.21
  sp: 2.1.4
  stringr: 1.5.1
  survival: 3.7.0
  tibble: 3.2.1
  tidyr: 1.3.2
  tidyverse: 2.0.0
  viridis: 0.6.5
  zoo: 1.8.12

**Non-standard hardware:** None required. All analyses run on a standard desktop or laptop computer (no GPU or high-memory cluster needed).

## 2. Installation guide

1. Install [R](https://cran.r-project.org/) (version 4.0 or later recommended) and, optionally, [RStudio](https://posit.co/download/rstudio-desktop/).
2. Install the required packages:

   ```r
   install.packages(c("dplyr", "tidyr", "tidyverse", "purrr", "tibble", "stringr",
                       "readr", "readxl", "haven", "zoo", "lubridate", "sf", "sp",
                       "raster", "rworldmap", "dlnm", "survival", "ggplot2",
                       "viridis", "here"))

   # rgdal is archived on CRAN; install an archived version if needed:
   # remotes::install_version("rgdal", version = "1.6-7")
   ```

3. Clone or download this repository:

   ```bash
   git clone https://github.com/annakdim/Heat-and-IRA.git
   cd Heat-and-IRA
   ```

**Typical install time:** approximately 15-30 minutes on a normal desktop computer, mainly due to compiling the spatial packages (`sf`, `raster`, `rgdal`).

## 3. Demo

Because the clinical data cannot be shared, `demo/` provides a small **simulated** dataset with the same variable names and structure expected by the analysis scripts, so that the pipeline can be installed and verified independently of the real data.

**Instructions to run the demo:**

```r
setwd("Heat-and-IRA")     # repository root
source("demo/simulate_demo_data.R")   # creates demo/demo_data.RData (synthetic data)
source("demo/run_demo.R")             # fits the case-crossover models and prints results
```

**Expected output:** a table of odds ratios (OR), 95% confidence intervals, and p-values for a small set of simulated heat-exposure definitions, printed to the console and saved to `demo/output/demo_results.csv`. Because the input data are randomly simulated, exact values will differ each run and are for illustration only — they are **not** a reproduction of the manuscript's results.

**Expected run time for demo:** less than 1 minute on a normal desktop computer.

## 4. Instructions for use

To run the analysis on your own data:

1. Prepare your dataset with the same structure as `demo/demo_data.RData` (columns `sample`, `country`, `strata_id`, `case`, `stillbirth`, `apgar_score`, `birth_weight`, and temperature-percentile columns named `<metric>_pct_lag<L>` for `metric` in `tmax`, `tmean`, `tmin`, `dtr` and lag `L` in 0-4).
2. Run the scripts in order:
   - `1_Clean the data/` — clean and merge the raw cohort files (adapt the input file paths to your own data).
   - `2_Generate temperature exposure-ERA5/` — link births to ERA5 temperature grids and compute exposure percentiles (requires ERA5 reanalysis data, freely available from the [Copernicus Climate Data Store](https://cds.climate.copernicus.eu/)).
   - `3_Analysis case-crossover/` — fit the time-stratified case-crossover models (main analysis, healthy-case comparison, full-cohort comparison, and negative-control test).
3. Outputs (model estimates and figures) are written to the working directory / subfolders referenced at the top of each script.

## License

This code is released under the [MIT License](LICENSE), a permissive license approved by the Open Source Initiative.

## Citation

If you use this code, please cite the manuscript above.
