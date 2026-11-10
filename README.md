# Reproducibility Artifact: On the Uniform Size Validity of Welch's t-Test

This repository contains the complete replication code and certified arithmetic verification scripts for the paper:

**Title:** On the Uniform Size Validity of Welch's t-Test and the Consistency of Its Studentizing Variance Ratio  
**Author:** Mohsen Salehi

## Repository Structure

- `verify_certificate.py` — Python script performing the certified computation of both asymptotic quantities:
  - `L`, the exact limiting rejection probability of Welch's test along the diagonal sequence (Theorem 1);
  - `M`, the limit of the expected Gaussian-calibrated rejection probability `E[g(R)]` (Proposition 2);
  - certified interval bounds for both quantities, and the certified margin `M_lo - L_hi > 0.005`, establishing failure of the uniform calibration condition (Proposition 2; Appendix A).
  The script also verifies the boundary exclusion for `k ≤ 100` (underlying `q_k`) and certifies the lower bound `B_100 > 0.05` (Theorem 1).
- `simulate_welch.R` — R script executing the Monte Carlo simulations (200,000 replications, seed 110) across sample sizes {50, 200, 1000, 5000} to assess finite-sample rejection rates.
- `requirements.txt` — Python package dependencies.

## Software Versions & Environment

- Python: 3.10+ with `mpmath` (for high-precision and interval arithmetic) and standard-library modules `fractions` and `math`.
- R: 4.3.0+ (base packages only).

## Instructions for Replication

### 1. Certified verification of L and M (Python)

Run:

```bash
python verify_certificate.py
```

Expected output (recorded in Appendix A of the manuscript):

```text
L_lo = 0.05191344502566980459926
L_hi = 0.05191344502566980459926 + 6.9e-38
M_lo = 0.05731969721554166423991
M_hi = 0.05731969721554166423991 + 1e-44
M_lo - L_hi = 0.005406252189871859639...
certified: M_lo > L_hi : True
q_26 (sanity check) = 0.075518697500228881836
```

The certified margin `M_lo - L_hi > 0.005` establishes the strict inequality `M > L` and hence the failure of the uniform calibration condition (Proposition 2; Appendix A).

### 2. Monte Carlo simulation (R)

Run:

```R
source("simulate_welch.R")
```

This reproduces Table 1 of the manuscript.
