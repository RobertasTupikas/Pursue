# Pursue

### Microbiome **P**rofiling **U**sing ALR **R**egression with resample-supported **S**upported **U**ncertainty **E**stimation

Re-implementation of the **Songbird** algorithm (Morton *et al.*, 2019) in **R**, with additional functionality and a point-and-click **Shiny** interface.

- **Original paper:** [https://doi.org/10.1038/s41467-019-10656-5](https://doi.org/10.1038/s41467-019-10656-5)  
- **Python reference repo:** [biocore/songbird](https://github.com/biocore/songbird)  
- **Pursue repo:** [RobertasTupikas/Pursue](https://github.com/RobertasTupikas/Pursue)

---

## Running Pursue: server-side or local

You can run Pursue via the **Nature Research Center (NRC)** server or on your **local machine**.

---

### **Option 1 — NRC Server**

Connect through your browser at  
**[http://192.168.5.86:3838/pursue-simple/](http://192.168.5.86:3838/pursue-simple/)**  

> Currently accessible only on the NRC Wi-Fi network.  
> Skip to the **User Interface Tutorial** below if using this route.

---

### **Option 2 — Local installation**

```r
# Install Pursue core packages
remotes::install_github("RobertasTupikas/songbird-for-R")  # songbirdR
```

**Conda / Python environment requirements:**

| Package | Version |
|----------|----------|
| tensorflow | 2.14.0 |
| tensorflow-probability | 0.23.0 |
| numpy | 1.26.4 |

Load this environment from inside the UI after starting Pursue.

---

### **Launch the Shiny interface**

```r
# Run from file
source("shinyapp.R")

# Or from directory
shiny::runApp("path/to/appdir")
```

**Required R packages:**

```r
library(shiny)
library(shinyjs)
library(DT)
library(plotly)
```

---

## User Interface Tutorial

### **1 — Load Conda environment**
- Open the **Model Settings** tab (left panel).  
- Enter your Conda environment name and click **Load Conda Environment**.  
- On the NRC server, the environment is always named **pursue-tf**.  
- Locally, specify the name of your custom environment.  
- A status message will indicate success or error.

---

### **2 — Upload data**
- **OTU Table (CSV):** rows = samples; columns = OTUs/features; first column = sample IDs.  
- **Metadata (CSV):** rows = same samples; columns = covariates; first column = sample IDs.  
- Sample IDs **must match** between files.  
- Categorical variables must be **character/factor**. Numeric codes (e.g., 0/1/2) are treated as continuous unless converted.

---

### **3 — Select model variables**
- In **Model Formula Columns**, choose one or more metadata variables.  
- For categorical variables, select their **Reference level** in the generated dropdowns.

---

### **4 — Configure training parameters**
- Scroll to the **Model Settings** section to view and edit parameter values.  
- Default settings are suitable for a first run; see the **Parameters Tutorial** below for detailed explanations.

---

### **5 — Run regression**
- Click **Run Regression**.  
- A progress message will appear (“Running…” then “Done.”).  
- Results become accessible in the **Results Table** (top-right panel).

| Output matrix | Description |
|----------------|-------------|
| `beta_mean` | Mean regression coefficients |
| `beta_sd` | Coefficient standard deviations |
| `beta_pval` | Empirical p-values |
| `beta_fdr` | FDR-adjusted q-values |

Use **Download CSV** to export the active table or **Download RDS** to save all outputs.

---

### **6 — Plot results**
- Switch to the **Plot Settings** tab.  
- Choose a **Coefficient** to plot.  
- Optionally specify **Highlighted OTUs** (comma-separated).  
- Adjust **Y-axis limits** and **q-value threshold**.  
- Click **Generate Plot**.

**Displayed output:**
- Rank plot of OTUs ordered by coefficient magnitude.  
- Points with *q ≤ threshold* appear as “significant.”  
- Highlighted OTUs show diamond-outline markers.  
- Hover text includes OTU, rank, β-mean, and q-value.  
- Clicking any point highlights the same OTU in the browse table.

---

## Parameters Tutorial

### **Quick start**
For a high-quality first run:
```r
epochs     = 10000
n_boot     = 200
n_perms    = 2000
```
Leave other parameters at default values.

---

### **Differential_prior**

Defines the standard deviation of the Normal prior on regression coefficients (β). It controls how freely differential effects can vary by determining the strength of regularization (shrinkage) toward zero.  

**Practical:**  
A smaller `differential_prior` (e.g., 0.1–0.5) heavily penalizes large coefficients, resulting in smoother but more conservative outputs; useful for small or noisy datasets where overfitting is a risk.  
Larger values (e.g., 5–10) loosen that constraint, allowing stronger log-fold changes to emerge—useful when you expect strong biological signals or have many samples.  
The Songbird reference implementation uses a comparable scale, typically with σ = 1 as a balanced default.  

**Guideline:**  
- Start with 1.  
- If effects look muted (few significant taxa), raise to 3–5.  
- If results fluctuate across runs, lower to 0.5–1 for stability.

---

### **Clipnorm**

Specifies the maximum global L2 norm for gradient clipping during training. Prevents the optimizer from taking excessively large steps by scaling gradients that exceed this value.  

**Practical:**  
Smaller values (e.g., 1–5) enforce stricter stability but may slow convergence; higher values (10–50) speed learning but risk numerical instability. Start at 10 (Songbird default) and adjust if loss oscillates or diverges.  

**Guideline:**  
- Begin with 10.  
- If loss spikes or training diverges, lower to 5.  
- If learning is extremely slow or gradients vanish, raise slightly (15–20).  

---

### **Random-effect λ**

Regularization coefficient for random-effect terms, controlling their variance penalty. Applies L2 shrinkage to random offsets such as batches or subjects.  

**Practical:**  
Higher values (e.g., 1e-2) force random effects closer to zero, simplifying the model and reducing overfitting. Lower values (1e-4–1e-6) allow more flexibility across clusters. Adjust depending on dataset size and number of random groups.  

**Guideline:**  
- Start at 1e-3.  
- If model underfits (too conservative), lower to 1e-4.  
- If random-effect variance seems excessive, raise to 1e-2.  

---

### **Learning rate**
 
Step size for the Adam optimizer, determining how far weights move in response to each gradient update.  

**Practical:**  
Small learning rates (1e-4–1e-3) ensure stable convergence; large ones (≥1e-2) may cause divergence. Begin with 1e-3 and decrease by half if loss oscillates, or increase slightly if learning is too slow.  

**Guideline:**  
- Start at 1e-3.  
- If training loss fluctuates or diverges, lower to 5e-4.  
- If validation loss plateaus too early, increase to 2e-3.  

---

### **Batch size**

Number of samples processed per training step before gradient averaging. Controls the stochasticity of gradient descent.  

**Practical:**  
Smaller batches (5–10) introduce noise that acts as regularization, helping generalization; larger batches (20+) yield smoother but slower updates. Typical microbiome datasets benefit from batch sizes around 5–15.  

**Guideline:**  
- Default to 5–10 for most datasets.  
- Use smaller batches if overfitting.  
- Use larger batches if training is unstable or memory allows.  

---

### **Number of test samples**

Number of subsampled test replicates evaluated per training iteration to assess generalization.  

**Practical:**  
Higher values (10–20) yield more stable validation metrics but increase runtime. For small datasets (<100 samples), 5 is sufficient; increase if validation loss is noisy.  

**Guideline:**  
- Use 5 for quick tests.  
- Use 10–20 for final results or noisy validation curves.  
- Increase gradually until validation loss stabilizes.  

---

### **Seed**

Random initialization seed ensuring deterministic model behavior across runs.  

---

### **Parallel mode**

Defines the parallel processing backend for resampling operations (`"no"`, `"multicore"`, or `"snow"`).  

**Guideline:**  
- Use `"multicore"` for local Linux/Mac machines.  
- Use `"snow"` for Windows or SLURM clusters.  
- Use `"no"` when debugging serially.  

---

### **Number of cores**

Number of CPU cores allocated for parallel execution during resampling. Speeds up runtime.

**Practical:**  
Match to the number of available physical cores minus one to avoid saturation. Setting too high may slow overall performance due to thread contention.  

**Guideline:**  
- For the NRC server, cap at 16. Right now, this is the soft cap because of optimization issues.
- On your local machine or institutional server, set to as many as there are minus one. 

---

### **Epochs**

Maximum number of epochs (complete training passes) per resample iteration.  

**Practical:**  
Acts as an upper limit for training; used in conjunction with early stopping. If early stopping is disabled, this determines total training length. Raise when convergence is slow; lower for quick diagnostics.  

**Guideline:**  
- Start with 1000 for test runs.
- When early stopping is turned on, set this to a large value (20000). It wont matter much, since early stopping will automatically stop at a number of epochs where conveargence is achieved, however you should leave a reasonable ceiling so it does not underfit.

---

### **Number of bootstraps**

Number of bootstrap resamples used to estimate uncertainty (confidence intervals) around coefficient estimates.  

**Practical:**  
Increasing this number improves the precision of the estimated variability but significantly increases computation time. Too few (e.g., <50) can make intervals noisy; too many (e.g., >1000) may waste resources with little added benefit.  

**Guideline:**  
- Use 100–200 for standard analyses.  
- Use ≥500 for publication-quality uncertainty intervals.  
- For testing or debugging, use 10–20 to save time.  

---

### **Number of permutations**

Number of permutation replicates for constructing empirical null distributions used in p-value estimation.  

**Practical:**  
Low numbers give coarse p-value resolution (minimum p ≈ 1/(n_perms + 1)), while high numbers improve accuracy but slow computation.  

**Guideline:**  
- Minimum 200 for stable significance testing.  
- 1000–4000 for final or published analyses.  
- Fewer than 100 only for exploratory testing.  

---

### **Early stopping**

Toggle controlling whether the training halts automatically when validation loss ceases to improve.  

**Practical:**  
Helps prevent overfitting and shortens training time. If disabled, model trains until maximum `epochs`.  

**Guideline:**  
- Keep `TRUE` for most use cases. 
- Set to `FALSE` only for diagnostics or when manually monitoring convergence.

---

### **Patience**

Number of consecutive evaluation intervals (`eval_every`) without improvement before early stopping triggers.  

**Practical:**  
Lower values make early stopping respond quickly to stagnation; higher values give the model more time to recover from temporary plateaus.  

**Guideline:**  
- Start with 20.  
- Increase to 50 for noisy validation loss.  
- Decrease to 10 for faster stopping.  

---

### **Minimum epochs**

Minimum number of epochs required before early stopping can activate. Ensures sufficient training time before evaluation begins.  

**Practical:**  
Prevents premature stopping during early fluctuations in validation loss.  

**Guideline:**  
- Keep around 200 by default. This will practically never change.
- Lower (100) for quick experiments.  
- Raise (300–500) for complex models or high noise. 

---

### **Evaluate every N steps**

Defines how frequently (in training steps) validation loss is checked for improvement.  

**Practical:**  
Smaller values give more granular feedback but slightly increase overhead; larger values reduce overhead but might miss early signs of convergence.  

**Guideline:**  
- Default 10–20 for most datasets.  
- Increase if training logs are too verbose.  
- Decrease if validation loss is erratic.  

---

### **Tolerance**

Minimum relative improvement in validation loss to be considered progress during early stopping.  

**Practical:**  
Lower `tol` values make convergence stricter and prolong training; higher values stop earlier.  

**Guideline:**  
- Default 1e-4 for general use.  
- Use 1e-5 for precise convergence on stable datasets.  
- Use 1e-3 for faster but less precise training.  

---

### **Learning rate plateau**

Switch to enable automatic learning-rate decay when validation loss stagnates.  

**Practical:**  
When enabled, helps optimization escape shallow minima by reducing step size. Keeps learning steady and prevents oscillations.  

**Guideline:**  
- Keep `TRUE` for stability.  
- Disable if you prefer to control learning rate manually.  

---

### **Learning rate decay factor**

Multiplicative factor applied to the learning rate when plateau detection triggers.  

**Practical:**  
Controls how aggressively the learning rate drops. Too small (<0.3) may stall training; too large (>0.9) may not reduce learning rate enough to stabilize.  

**Guideline:**  
- Default 0.5.  
- Use 0.3 for strong decay on unstable models.  
- Use 0.8–0.9 for gentler schedules.  

---

### **Learning rate minimum**

Minimum allowable learning rate after decays. Prevents the optimizer from becoming effectively frozen.  

**Practical:**  
If training halts too early or becomes static, the minimum learning rate may be too low.  

**Guideline:**  
- Default 1e-6.  
- Increase to 1e-5 if training stops before full convergence.  
- Decrease to 1e-7 only if extremely slow, precise convergence is required.  

---


### **Performance notes**
A regression on an **OTU table (32 × 204)** with  
`epochs = 10000` (typically 4000–8000 with early stopping),  
`n_boot = 200`, `n_perms = 2000`, and `n_cores = 16`  
runs in **≈ 20–50 minutes**.

> Parallel scaling is currently sub-linear: using 16 vs 32 cores yields minor runtime differences.  
> Please limit NRC server runs to **≤ 16 cores** to conserve resources.

---
