# songbird-for-R — Shiny App Workflow

Reimplementation of the Songbird algorithm by Morton et al. (2019) in R with additional functionalities and a **point-and-click Shiny interface**.

Original publication:  
https://doi.org/10.1038/s41467-019-10656-5

Python reference repository:  
https://github.com/biocore/songbird

BMDD preprint (zero-imputation via proportion modeling):  
https://doi.org/10.1101/2025.05.08.652808  
BMDD repo: https://github.com/zhouhj1994/BMDD

Selbal (balance selection) repo:  
https://github.com/malucalle/selbal

Songbird-for-R repo:  
https://github.com/RobertasTupikas/songbird-for-R

---

## Installation

```r
install.packages(c("shiny","shinyjs","DT","plotly","remotes"))

# Core packages from GitHub
remotes::install_github("RobertasTupikas/songbird-for-R")  # songbirdR
remotes::install_github("malucalle/selbal")                # selbal
remotes::install_github("zhouhj1994/BMDD")                 # BMDD (optional)
```

**Conda / Python environment (recommended):** create an environment that includes at least:
- `tensorflow==2.14.0`
- `tensorflow-probability==0.23.0`
- `numpy==1.26.4`

You will load this environment from inside the app (via **reticulate**).

---

## Launch

```r
# If your file defines shinyApp(ui, server):
source("shinyapp.R")

# Or, if the app is a directory:
shiny::runApp("path/to/appdir")
```

The app expects these R packages to be available in the session:
```r
library(shiny)
library(shinyjs)
library(DT)
library(plotly)
```

---

## Data format

- **OTU table (CSV):** **rows = samples**, **columns = OTUs/features**; first column = sample IDs.  
- **Metadata (CSV):** **rows = the same samples**, columns = covariates; first column = sample IDs.  
- Sample IDs must match between files.  
- Categorical variables should be **character/factor**. Numeric codes like `0/1/2` are treated as continuous unless converted.

---

## Workflow

Below is a step-by-step guide to using the **Shiny interface** that comes with `songbird-for-R`. It assumes you have installed the packages and launched the app as described above.

> **Data orientation reminder:** both **OTU table** and **metadata** must have **samples as rows** (first column = sample IDs). Sample IDs should match between files.

### 1) Load (optional) Conda environment
- Open the **Model Settings** tab (left panel).
- Enter your Conda environment name (e.g., `songbird`) and click **Load Conda Environment**.
- You should see a **Loading…** message followed by either a **success** or **error** message.
- This step ensures the correct Python (TensorFlow) environment is used through **reticulate** (needed for GPU and for specific TF versions).

### 2) Upload data
- **OTU Table (CSV):** rows = samples; columns = OTUs/features; **first column = sample IDs**.
- **Metadata (CSV):** rows = the same samples; columns = covariates; **first column = sample IDs**.
- After upload, the app will expose your **metadata columns** for model selection.

### 3) Choose the model variables
- In **Model Formula Columns**, select one or more metadata variables you want the regression to learn from.
- For **categorical** variables, the UI will show a **Reference level** dropdown for each selected variable. Pick your baseline category there.

### 4) Configure training parameters
The model parameters you can set in the interface include:
- `differential_prior`
- `learning_rate`
- `clipnorm`
- `batch_size`
- `seed`
- `num_test_samples`
- `epochs`
- `n_bootstrap` (number of bootstrap resamples)
- `n_cores` (CPU parallelization for bootstraps)
- `use_gpu` (if your Conda env has TF-GPU)

> **Tip (ALR denominator drop):** Internally, the first OTU in the OTU table is used as the ALR denominator and thus **dropped** from the ranking. If that OTU might be important, consider **running the regression twice**, changing which OTU is first.

### 5) Run the regression
- Click **Run Regression**.
- You will see a **running** message; on completion, a **done** message appears.
- Results are exposed in the top-right **browse table**:
  - **beta_mean**, **beta_sd**, **beta_pval** (selectable via a dropdown)
  - You can **search** and **sort** features.
  - Use **Download CSV** to export the currently viewed matrix or **Download whole RDS** to save all results.

### 6) Plot Settings → Rank plot (interactive)
- Switch to the **Plot Settings** tab (left panel).
- Choose the **Coefficient to Plot** (e.g., a level of a categorical variable).
- Optionally **Highlight OTUs** (comma-separated list).
- Set **Y-axis Min/Max** and **P-value Threshold**.
- Click **Generate Plot**.

**What you’ll see:**
- A **rank plot** of OTUs (x = rank by `beta_mean`, y = `beta_mean`).
- Points with **p ≤ threshold** are drawn as a “significant” trace; others are faded.
- Your **highlighted OTUs** appear with a distinct diamond-outline marker.
- **Hover** shows OTU, rank, beta mean, and p-value.
- **Click** any point to **select that OTU in the browse table** on the right.

### 7) Build numerator/denominator sets from the table
- In the top-right table, **select one or more rows** (OTUs).
- Use **Add to Numerator** or **Add to Denominator** to append the selected OTUs to the inputs below.
- You can also **type OTU IDs manually** (comma-separated). The buttons append the **currently selected** table features.

### 8) Log‑ratio visualization & statistics
- Confirm the **Numerator OTUs** and **Denominator OTUs** in their text boxes.
- Pick a **Metadata Variable** for analysis.
- Choose **Analysis Type**:
  - **categorical** → the app renders grouped plots and annotates a suitable test (e.g., t/Wilcoxon/ANOVA/Kruskal–Wallis depending on distribution and group number).
  - **continuous** → the app renders a regression/correlation view with a fitted trend line.
- Click **Visualize Log‑Ratio** to compute and plot.
- Use **Download Log‑Ratios** to export the per‑sample ratios as CSV for use elsewhere.

> If you want to incorporate **zero‑imputation** using **BMDD** in your workflow, run BMDD on your counts first and then feed the imputed abundances into your log‑ratio analyses.

### 9) SELBAL balances (feature selection for ratios)
- At the bottom panel, choose the **Metadata Variable**.
- Toggle **Categorical Variable** if it’s a category.
- Set **Number of folds** and **Number of iterations**.
- Click **Get Balances**.
- You’ll see a progress message; when finished, the app prints **Numerator OTUs** and **Denominator OTUs** selected by SELBAL.
- Copy these into the **log‑ratio inputs** above to visualize and compute statistics.

**Performance note:** SELBAL is compute‑heavier than the regression; in practice it can be ~**6× slower** when the regression runs on **4 cores**.

---

## References
- Morton JT, et al. (2019). Establishing microbial composition measurement standards with reference frames. *Nature Communications*, 10, 2719. https://doi.org/10.1038/s41467-019-10656-5
- Calle ML, et al. (2019). Predictive microbiome analysis with balances. *mSystems*, 4(1), e00016-19. https://doi.org/10.1128/mSystems.00016-19

