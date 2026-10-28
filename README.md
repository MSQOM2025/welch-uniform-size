# Reproducibility Artifact: On the Uniform Size Validity of Welch's $t$-Test

This repository contains the complete replication code and certified arithmetic verification scripts for the paper:
> **Title:** On the Uniform Size Validity of Welch's $t$-Test and the Consistency of Its Studentizing Variance Ratio  
> **Author:** Mohsen Salehi  

## Repository Structure

- `verify_certificate.py`: Python script utilizing exact rational arithmetic (`fractions.Fraction`) to bound the normal quantile, verify the boundary conditions for $k \le 100$, and certify the rational lower bound $B_{100} > 0.05$.
- `simulate_welch.R`: R script executing the Monte Carlo simulations ($200,000$ replications, seed $110$) across various sample sizes to verify the finite-sample empirical size inflation.
- `requirements.txt`: Python package dependencies.

---

## Software Versions & Environment

The code has been tested under the following environment:
- **Python:** Version 3.10+ (utilizing only standard libraries like `fractions` and `math`).
- **R:** Version 4.3.0+ (utilizing base packages).

---

## Instructions for Replication

### 1. Certified Rational Arithmetic Verification (Python)
To run the exact rational verification of the lower bound $B_{100}$ and cross-multiply the inequality, execute:
```bash
python verify_certificate.py
