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
## Computational verification notes

The certified computations in `verify_certificate.py` use high-precision interval arithmetic via the `mpmath` library, with a working precision of 60 decimal digits. All intermediate quantities are computed with conservative directed-rounding margins, so that the reported intervals provably contain the exact values.

### What is certified

- The normal quantile `z_0.975` is bounded by a rational interval obtained from Machin's formula for π and a Maclaurin expansion of Φ with a rigorous alternating-series remainder.
- For each `k ≤ 100`, the indicator in `q_k` is verified to be constant on the certified bracket, so `q_k` is determined exactly.
- The Poisson weights `π_k` are computed with directed rounding, and the series are truncated at `k = 200` with a geometric tail bound (< 10⁻¹²⁴ beyond `k = 200`, and < 6.9 × 10⁻³⁸ beyond `k = 100` for `L`).

### Reproducibility

The script is fully deterministic: it uses no random number generation, no time-dependent inputs, and fixed precision settings. Running it repeatedly produces identical output.

### Implementation note: `mpmath` vs. `fractions.Fraction`

The current implementation uses `mpmath` (high-precision floating-point arithmetic with interval margins) rather than Python's `fractions.Fraction` for exact rational arithmetic. For the precision required here (final accuracy ~10⁻²⁵), the 60-digit `mpmath` precision is more than sufficient, and the directed-rounding margins ensure the reported intervals are rigorous.

If exact rational arithmetic is preferred, the script can be adapted to use `fractions.Fraction` throughout; the mathematical structure of the computation is unchanged.

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
