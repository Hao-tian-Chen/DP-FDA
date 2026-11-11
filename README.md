# DP-FDA: Privacy and Smoothing in Functional Mean Estimation

Simulation code and results accompanying [**On the Order of Privacy and Smoothing in Functional Mean Estimation**](https://arxiv.org/abs/2611.09799), by Haotian Chen, Yu-Wei Chen, Qifan Song, and Jordan Awan.

The project compares two procedures for functional mean estimation under a common Gaussian differential privacy (GDP) budget:

- **Fit-Then-Protect (FTP):** fit each observed trajectory, average the fitted functions, and add Gaussian privacy noise.
- **Protect-Then-Fit (PTF):** clip and privatize the coordinate means on the observation grid, then fit a function to the privatized vector.

The experiments study fitting dimension, fitting radius, and observation noise in finite-dimensional polynomial spaces and a Gaussian-kernel reproducing kernel Hilbert space (RKHS). Accuracy is measured relative to the latent sample mean: integrated mean squared error (IMSE) for polynomial spaces and mean squared error on the observation grid (MSE) for the RKHS.

## Repository contents

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

The saved data and figures are included. To inspect the results, open the PDFs in `figures/`. To regenerate the plots, use the saved data; rerunning the simulations is not required.

## Setup

Use R 4.1 or newer (the visualization scripts use the native `|>` pipe).

Install the required packages once:

```r
install.packages(c("ggplot2", "patchwork", "dplyr", "tidyr", "foreach", "doParallel"))
```

Clone or download this repository. Set the R working directory to its root, where the four `.R` scripts and the `data/` and `figures/` folders are located:

```r
setwd("/path/to/DP-FDA") # Replace with your local repository path
```

All input and output paths are relative to this directory. Running from a different directory can cause a file-not-found error.

## Recreate the figures from the included results

Run in an interactive R/RStudio session:

```r
source("Poly_Visualization.R")
source("RKHS_Visualization.R")
```

These scripts read the `.rds` files in `data/` and overwrite the corresponding PDFs in `figures/`.

`RKHS_Visualization.R` also constructs `table_wide` and `table_display` for the noiseless RKHS results, including errors, ratios, and PTF constraint-activation probabilities. It opens these objects using `View()`. 

| Experiment | Polynomial input (`data/Poly/`) | Polynomial output (`figures/Poly/`) |
| --- | --- | --- |
| Noiseless fitting-dimension comparison | `Poly_sigma0_0_B1p1.rds` | `IMSE_Ratio_mu_e-1-1.pdf` |
| Fitting-radius comparison | `Poly_Bprime.rds` | `Poly_Bprime_choices.pdf` |
| Observation-noise comparison | `Poly_noisy_sqrtlog.rds` | `IMSE_Ratio_noisy.pdf` |

| Experiment | RKHS input (`data/RKHS/`) | RKHS output (`figures/RKHS/`) |
| --- | --- | --- |
| Noiseless fitting-dimension comparison | `RKHS_sigma0_0_B1p1.rds` | `MSE_Ratio_mu_e-1-1.pdf` |
| Fitting-radius comparison | `RKHS_Bprime.rds` | `RKHS_Bprime_choices.pdf` |
| Observation-noise comparison | `RKHS_noisy.rds` | `MSE_Ratio_noisy.pdf` |

## Simulation settings

The experiment calls at the bottom of each simulation script use sample sizes `2^(2:12)`, latent norm bound `B = 1`, 1,000 Monte Carlo repetitions per setting, and seed 42.

| Experiment | Polynomial spaces | Gaussian-kernel RKHS |
| --- | --- | --- |
| Fitting dimension, without observation noise | `m = 2, 4, 6, 8`; `T = 20`; `mu = 0.1, 1`; `B' = 1.1` | `T = 2, 4, 6, 8`; `mu = 0.1, 1`; `B' = 1.1` |
| Fitting radius | `m = 4`; `T = 20`; `mu = 1`; noise SD = 1 | `T = 4`; `mu = 1`; noise SD = 1 |
| Observation noise | `m = 2, 4, 6, 8`; `T = 20`; `mu = 1`; `B' = sqrt(log(n))` | `T = 2, 4, 6, 8`; `mu = 1`; `B' = log(n)` |

The radius experiment compares `B`, `sqrt(log(n))`, `log(n)`, `sqrt(n)`, and `n^(2/3)`. The noisy experiment compares observation-noise standard deviations 0.1, 1, and 10. Here `log` is the natural logarithm. The polynomial and RKHS noisy experiments intentionally use different fitting-radius sequences.

In the RKHS experiments, latent functions are generated using 20 kernel centers and a Gaussian kernel with bandwidth 0.2. The numerical fitting routines use a regularization floor of `1e-8`.

## Rerun the simulations

From the repository root, run these commands in separate R sessions:

```r
source("Polynomial.R")
```

```r
source("RKHS.R")
```

Each script runs all three experiment blocks at its end and saves new summaries to the corresponding `data/` folder, replacing the included `.rds` files. The full experiments involve many repeated constrained fits and can take substantial time. To try a smaller run, edit the experiment calls first: reduce `n_vec`, set a smaller `reps`, and set `n_cores` explicitly.

The scripts parallelize settings using `foreach` and `doParallel`. With `n_cores = NULL`, they select the worker count automatically. They assign random-number streams to workers; reruns with different worker counts, software versions, or numerical environments need not reproduce the saved results exactly. Use the included summaries to recreate the displayed figures directly.

After generating new results, rerun the visualization scripts to update the PDFs. Keep the expected experiment settings and filenames when using the supplied plotting scripts without modification.

## Read and interpret the saved results

The `.rds` files contain data frames of Monte Carlo summaries, rather than all individual simulated trajectories. For example:

```r
results <- readRDS("data/RKHS/RKHS_sigma0_0_B1p1.rds")
names(results)
head(results)
```

Latent functions are held fixed across repetitions within a setting; observation and privacy noise are regenerated. Reported error ratios are ratios of Monte Carlo mean squared errors, not averages of repetition-level ratios. Ratios below one favor FTP; ratios above one favor PTF.

The main constrained comparison uses the following columns:

| Column | Meaning |
| --- | --- |
| `IMSE_A3_FTP_constrained` | Error of constrained FTP |
| `IMSE_A2_PTF_constrained` | Error of constrained PTF |
| `Ratio_A3_A2` | FTP error divided by PTF error |
| `A3_prop_constrained` | FTP constraint-activation summary |
| `A2_prop_constrained` | PTF constraint-activation summary |

For RKHS results, the columns retain the `IMSE_` prefix for historical naming consistency, but the reported error is **grid MSE**. Additional columns store auxiliary procedures, such as norm-clamped FTP and unconstrained PTF, and diagnostic summaries.

The dashed `1/m` and `1/T` lines are motivated by the noiseless theoretical comparisons. In the noisy plots they are reference levels, not theoretical guarantees for those noisy settings.
