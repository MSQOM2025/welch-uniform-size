from decimal import Decimal, localcontext
from fractions import Fraction as F
from math import comb, factorial
import sys


def print_step(msg):
    print(f"\n[+] {msg}")


# =============================================================================
# 1. Alternating Series Framework with Parity-Aware Bounds
# =============================================================================
def alternating_series_bounds(term_func, N):
    """Rigorous lower and upper bounds for an alternating series sum_{n=0}^N (-1)^n a_n

    where a_n > 0 is strictly decreasing.
    """
    s = sum(((-1) ** n) * term_func(n) for n in range(N + 1))
    a_next = term_func(N + 1)
    if N % 2 == 0:
        # N is even: the last added term was positive.
        # s is an UPPER bound, and s - a_next is a LOWER bound.
        return s - a_next, s
    else:
        # N is odd: the last added term was negative.
        # s is a LOWER bound, and s + a_next is an UPPER bound.
        return s, s + a_next


# =============================================================================
# 2. Certified Bounds for pi and sqrt(2*pi) via Machin's Formula
# =============================================================================
print_step("Step 1: Certifying pi and sqrt(2*pi) via Machin's identity...")


def atan_bounds(q, N=36):
    return alternating_series_bounds(
        lambda n: F(1, (2 * n + 1) * (q ** (2 * n + 1))), N
    )


# pi = 16 * arctan(1/5) - 4 * arctan(1/239)
a5_lo, a5_hi = atan_bounds(5, 36)
a239_lo, a239_hi = atan_bounds(239, 36)
pi_lo = 16 * a5_lo - 4 * a239_hi
pi_hi = 16 * a5_hi - 4 * a239_lo
assert pi_lo < pi_hi, "Assertion failed: pi bounds inconsistent"

# Rational bounds for sqrt(2*pi)
sqrt2pi_lo = F(250662827463100050, 10**17)
sqrt2pi_hi = F(250662827463100051, 10**17)
assert (
    sqrt2pi_lo**2 < 2 * pi_lo
), "Assertion failed: sqrt2pi_lo is not a valid lower bound"
assert (
    sqrt2pi_hi**2 > 2 * pi_hi
), "Assertion failed: sqrt2pi_hi is not a valid upper bound"
print("    pi in [pi_lo, pi_hi] and sqrt(2*pi) bracket certified.")


# =============================================================================
# 3. Certified Rational Bracket for z_{0.975}
# =============================================================================
print_step("Step 2: Certifying z_{0.975} quantile bracket...")


def phi_integral_bounds(x, N=36):
    # Taylor series for int_0^x exp(-t^2/2) dt
    return alternating_series_bounds(
        lambda n: (x ** (2 * n + 1)) / (F(2**n * factorial(n) * (2 * n + 1))),
        N,
    )


def Phi_bounds(x, N=36):
    I_lo, I_hi = phi_integral_bounds(x, N)
    phi_lo = F(1, 2) + I_lo / sqrt2pi_hi
    phi_hi = F(1, 2) + I_hi / sqrt2pi_lo
    return phi_lo, phi_hi


z_lo = F(1959963984540054, 10**15)
z_hi = F(1959963984540055, 10**15)

phi_lo_at_zlo, phi_hi_at_zlo = Phi_bounds(z_lo, 36)
phi_lo_at_zhi, phi_hi_at_zhi = Phi_bounds(z_hi, 36)

# Prove 0.975 = 39/40 is strictly bracketed: Phi(z_lo) < 0.975 < Phi(z_hi)
assert phi_hi_at_zlo < F(
    39, 40
), "Assertion failed: z_lo is not strictly below z_{0.975}"
assert phi_lo_at_zhi > F(
    39, 40
), "Assertion failed: z_hi is not strictly above z_{0.975}"
print(f"    z_{{0.975}} in ({float(z_lo)}, {float(z_hi)}) certified.")


# =============================================================================
# 4. Strict Boundary Avoidance and Computation of q_k (k = 1 to 100)
# =============================================================================
print_step(
    "Step 3: Checking boundary avoidance and computing q_k (k = 1 to 100)..."
)

q_values = {}
for k in range(1, 101):
    num_rejected = 0
    kz_lo_sq = k * (z_lo**2)
    kz_hi_sq = k * (z_hi**2)

    for j in range(k + 1):
        d_sq = (2 * j - k) ** 2

        # Check: d^2 is strictly outside the bracket [k*z_lo^2, k*z_hi^2]
        assert (
            d_sq < kz_lo_sq or d_sq > kz_hi_sq
        ), f"Boundary atom ambiguity detected at k={k}, j={j}"

        if d_sq > kz_hi_sq:
            num_rejected += comb(k, j)

    q_values[k] = F(num_rejected, 2**k)

# Audit specific check for k = 26
assert q_values[26] == F(
    2533987, 33554432
), f"Mismatch at k=26: {q_values[26]}"
print("    q_26 = 2533987 / 33554432 verified.")
print(
    "    All (2j - k)^2 != k * z_{0.975}^2 certified strictly for 1 <= k <= 100."
)


# =============================================================================
# 5. Finite-Truncation Summation and B_100 Lower Bound Certification
# =============================================================================
print_step(
    "Step 4: Computing exact intermediate sum and certifying B_100 > 0.05..."
)

sum_terms = F(1, 20)  # k=0 term: 0.05 = 1/20
for k in range(1, 101):
    sum_terms += F(20**k, factorial(k)) * q_values[k]

with localcontext() as ctx:
    ctx.prec = 25
    dec_sum = Decimal(sum_terms.numerator) / Decimal(sum_terms.denominator)
    print(f"    Exact Intermediate Sum: {dec_sum}")
    # Verify the corrected transcription value
    expected_sum = Decimal("25186596.700274514")
    assert (
        abs(dec_sum - expected_sum) < Decimal("1e-6")
    ), "Intermediate sum deviates from audit"

# Rigorous upper bound for exp(20) via Taylor series up to N=200 + geometric tail
N_exp = 200
S_exp = sum(F(20**n, factorial(n)) for n in range(N_exp + 1))
tail_geom = F(20**201, factorial(201)) * F(202, 182)  # 1 / (1 - 20/202) = 202/182
exp20_upper = S_exp + tail_geom

# Certified Lower Bound B_{100}
B100_lower = sum_terms / exp20_upper

with localcontext() as ctx:
    ctx.prec = 32
    b100_dec = Decimal(B100_lower.numerator) / Decimal(B100_lower.denominator)
    print(f"    Certified B_100 Lower Bound: {b100_dec}")
    excess = b100_dec - Decimal("0.05")
    print(f"    Excess above alpha = 0.05:  {excess}")

# Final mathematical assertion of Theorem 1
assert (
    B100_lower > F(1, 20)
), "CRITICAL: B_100 lower bound is NOT strictly greater than 0.05!"

print_step("SUCCESS: All mathematical assertions certified without errors.")
sys.exit(0)
