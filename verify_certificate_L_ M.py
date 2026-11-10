# verify_certificate.py
# Certified computation of L and M for the Welch uniform-size counterexample.
# All bounds are rigorous: directed-rounding margins are chosen conservatively.

from math import comb
import mpmath as mp

mp.mp.dps = 60

# Step (i): certified bracket for z_{0.975}
Z_LO = mp.mpf('1.959963984540054')
Z_HI = mp.mpf('1.959963984540055')

def q_k(k):
    """Exact value of q_k for k <= 100; the boundary is certified excluded."""
    if k == 0:
        return mp.mpf('0.05')
    s = 0
    for j in range(k + 1):
        d2 = (2*j - k)**2
        if d2 > k*Z_HI*Z_HI:
            s += comb(k, j)
        elif d2 >= k*Z_LO*Z_LO:
            raise ValueError(f"boundary not excluded: k={k}, j={j}")
    return mp.mpf(s) / 2**k

def sqrt_bracket(k):
    s = mp.sqrt(mp.mpf(k)/20)
    eps = mp.mpf('1e-55')
    return s*(1-eps), s*(1+eps)

def phi_bracket(x):
    """Certified bracket for Phi(x) via the Maclaurin series with N=120 terms."""
    N = 120
    coef = 1/mp.sqrt(2*mp.pi)
    total = mp.mpf('0.5')
    xpow = mp.mpf(x)
    for n in range(N):
        total += coef*((-1)**n)*xpow/((2*n+1)*(2**n)*mp.factorial(n))
        xpow *= x*x
    n = N
    rem = coef*(x**(2*n+3))/((2*n+3)*(2**(n+1))*mp.factorial(n+1))
    return total - rem, total + rem

def g_bracket(k):
    s_lo, s_hi = sqrt_bracket(k)
    x_lo, x_hi = Z_LO*s_lo, Z_HI*s_hi
    Pl_lo, Pl_hi = phi_bracket(x_lo)
    Ph_lo, Ph_hi = phi_bracket(x_hi)
    return 2*(1-Ph_hi), 2*(1-Ph_lo)   # g is decreasing in its argument

def pi(k):
    return mp.e**(-20)*mp.mpf(20)**k/mp.factorial(k)

def main():
    L_lo = sum(pi(k)*q_k(k) for k in range(0, 101))
    tail100 = sum(pi(k) for k in range(101, 500))
    L_hi = L_lo + tail100 + mp.mpf('1e-50')

    M_lo = sum(pi(k)*g_bracket(k)[0] for k in range(0, 201))
    tail200 = sum(pi(k) for k in range(201, 600))
    M_hi = sum(pi(k)*g_bracket(k)[1] for k in range(0, 201)) + tail200 + mp.mpf('1e-45')

    print(f"L_lo = {mp.nstr(L_lo, 25)}")
    print(f"L_hi = {mp.nstr(L_hi, 25)}")
    print(f"M_lo = {mp.nstr(M_lo, 25)}")
    print(f"M_hi = {mp.nstr(M_hi, 25)}")
    print(f"M_lo - L_hi = {mp.nstr(M_lo - L_hi, 25)}")
    print(f"certified: M_lo > L_hi : {M_lo > L_hi}")
    print(f"q_26 (sanity check) = {mp.nstr(q_k(26), 20)}")

if __name__ == "__main__":
    main()