# DP-FDA: Privacy and Smoothing in Functional Mean Estimation

Simulation code and results accompanying [**On the Order of Privacy and Smoothing in Functional Mean Estimation**](https://arxiv.org/abs/2611.09799), by Haotian Chen, Yu-Wei Chen, Qifan Song, and Jordan Awan.

The project compares two procedures for functional mean estimation under a common Gaussian differential privacy (GDP) budget:

- **Fit-Then-Protect (FTP):** fit each observed trajectory, average the fitted functions, and add Gaussian privacy noise.
- **Protect-Then-Fit (PTF):** clip and privatize the coordinate means on the observation grid, then fit a function to the privatized vector.

The experiments study fitting dimension, fitting radius, and observation noise in finite-dimensional polynomial spaces and a Gaussian-kernel reproducing kernel Hilbert space (RKHS). Accuracy is measured relative to the latent sample mean: integrated mean squared error (IMSE) for polynomial spaces and mean squared error on the observation grid (MSE) for the RKHS.

## Repository Structure

```text
.
├── Polynomial.R                # Generate polynomial simulation results
├── RKHS.R                      # Generate RKHS simulation results
├── Poly_Visualization.R        # Recreate the three polynomial figures
├── RKHS_Visualization.R        # Recreate the three RKHS figures and inspect a summary table
├── data/
│   ├── Poly/                   # Saved polynomial simulation summaries (.rds)
│   └── RKHS/                   # Saved RKHS simulation summaries (.rds)
├── figures/
│   ├── Poly/                   # Polynomial figures (.pdf)
│   └── RKHS/                   # RKHS figures (.pdf)
└── README.md
```

## Setup

Use **R 4.1 or newer**. Clone or download the repository:

```bash
git clone https://github.com/Hao-tian-Chen/DP-FDA.git
cd DP-FDA
```

Install the required R packages:

```r
install.packages(c("ggplot2", "patchwork", "dplyr", "tidyr", "foreach", "doParallel"))
```

Set your R working directory to the repository root, which contains the four `.R` scripts:

```r
setwd("/path/to/DP-FDA") # Replace with your local path
```

## Usage

### Reproduce the figures

The saved simulation results are included, so the figures can be regenerated without rerunning the simulations. Run in an interactive R/RStudio session:

```r
source("Poly_Visualization.R")
source("RKHS_Visualization.R")
```

The scripts read from `data/` and save PDFs to `figures/`, replacing existing files. 

`RKHS_Visualization.R` also constructs `table_wide` and `table_display` for the noiseless RKHS results, including errors, ratios, and PTF constraint-activation probabilities. It opens these objects using `View()`. 

### Run the simulations

Run the desired script in a fresh R session:

```r
source("Polynomial.R") # Polynomial experiments
# OR
source("RKHS.R")       # RKHS experiments
```

Each script runs three experiments: noiseless fitting-dimension comparisons, fitting-radius comparisons, and observation-noise comparisons. Results overwrite the corresponding files in `data/`.

The experiment settings are specified at the bottom of each script, using 1,000 Monte Carlo repetitions and seed 42. Computation is parallelized with `foreach` and `doParallel`. For a smaller trial run, reduce `n_vec` and `reps`, and set `n_cores` explicitly. Full simulations can take substantial time, and results may vary with the worker count and software environment.

## Results

Polynomial experiments report integrated mean squared error (IMSE); RKHS experiments report MSE on the observation grid. Both measure error relative to the latent sample mean. Error ratios compare FTP to PTF, with values below one favoring FTP.

The `.rds` files contain simulation summaries and can be inspected with `readRDS()`. The main comparison uses `IMSE_A3_FTP_constrained`, `IMSE_A2_PTF_constrained`, and `Ratio_A3_A2`. RKHS columns retain the `IMSE_` prefix but measure grid MSE. 

The dashed `1/m` and `1/T` lines are motivated by the noiseless theoretical comparisons. In the noisy plots they are reference levels, not theoretical guarantees for those noisy settings.

See the paper for the theoretical results and detailed simulation settings.
