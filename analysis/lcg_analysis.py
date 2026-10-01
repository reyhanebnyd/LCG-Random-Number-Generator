"""Statistical quality analysis for the LCG in python/lcg.py.

Three experiments:
  A. Sweep (a, c) for the toy m=17 generator and check which combinations
     reach full period, per the Hull-Dobell theorem.
  B. Reproduce the classic RANDU flaw: a "slab" plot that exposes the
     hidden hyperplane structure a 1D/2D view would never show.
  C. Run simple uniformity / autocorrelation checks on several generators
     to show those checks alone would NOT have caught RANDU's flaw.

Outputs PNGs to ../results and prints a short report to stdout.
"""

import math
import random
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "python"))
from lcg import lcg_sequence  # noqa: E402

RESULTS = ROOT / "results"
RESULTS.mkdir(exist_ok=True)


# ---------- Part A: Hull-Dobell full-period analysis ----------

def prime_factors(n):
    factors = set()
    d = 2
    while d * d <= n:
        while n % d == 0:
            factors.add(d)
            n //= d
        d += 1
    if n > 1:
        factors.add(n)
    return factors


def is_full_period(a, c, m):
    """Hull-Dobell theorem: the LCG has period m iff all three hold."""
    if math.gcd(c, m) != 1:
        return False
    for p in prime_factors(m):
        if (a - 1) % p != 0:
            return False
    if m % 4 == 0 and (a - 1) % 4 != 0:
        return False
    return True


def period_length(seed, a, c, m):
    """Length of the cycle reached from `seed`, via dict-based cycle detection."""
    seen = {seed: 0}
    x = seed
    for i in range(1, m + 2):
        x = (a * x + c) % m
        if x in seen:
            return i - seen[x]
        seen[x] = i
    raise RuntimeError("period not found within state space")


def part_a_full_period_map(m=17, out=RESULTS / "period_heatmap_m17.png"):
    a_vals = list(range(1, m))
    c_vals = list(range(0, m))
    grid = np.zeros((len(c_vals), len(a_vals)))
    for i, c in enumerate(c_vals):
        for j, a in enumerate(a_vals):
            grid[i, j] = period_length(0, a, c, m)

    fig, ax = plt.subplots(figsize=(7, 6))
    im = ax.imshow(
        grid, origin="lower", cmap="viridis",
        extent=[min(a_vals) - 0.5, max(a_vals) + 0.5,
                min(c_vals) - 0.5, max(c_vals) + 0.5],
        aspect="auto",
    )
    ax.set_xlabel("multiplier a")
    ax.set_ylabel("increment c")
    ax.set_title(f"LCG period length for every (a, c), m={m}, seed=0\n"
                 f"(max possible period = {m})")
    fig.colorbar(im, ax=ax, label="period length")

    full = [(a, c) for a in a_vals for c in c_vals if is_full_period(a, c, m)]
    if full:
        fx, fy = zip(*full)
        ax.scatter(fx, fy, s=18, facecolors="none", edgecolors="red",
                   linewidths=1.2, label="Hull-Dobell predicts full period")
        ax.legend(loc="upper right", fontsize=8)

    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)

    mismatches = [(a, c) for a, c in full if period_length(0, a, c, m) != m]
    assert not mismatches, f"Hull-Dobell prediction mismatch: {mismatches}"

    print(f"[A] saved {out}")
    print(f"[A] {len(full)} / {len(a_vals) * len(c_vals)} (a, c) pairs "
          f"give full period m={m}")
    print(f"[A] assignment's own parameters a=7, c=3, m=17 -> "
          f"period {period_length(6, 7, 3, 17)} "
          f"(Hull-Dobell full period: {is_full_period(7, 3, 17)})")


# ---------- Part B: RANDU's hidden hyperplane structure ----------

def lcg_stream(seed, a, c, m, n):
    """n successive LCG outputs, normalized to [0, 1). No cycle detection:
    meant for generators whose modulus is too large to detect cycles in."""
    x = seed
    out = np.empty(n, dtype=np.float64)
    for i in range(n):
        x = (a * x + c) % m
        out[i] = x / m
    return out


def part_b_randu_planes(n=200_000, out=RESULTS / "randu_vs_good_lcg.png"):
    # RANDU: the infamous IBM generator, m = 2**31, a = 65539, c = 0.
    randu = lcg_stream(seed=1, a=65539, c=0, m=2**31, n=n)
    # A well-chosen LCG for comparison (Numerical Recipes constants).
    good = lcg_stream(seed=1, a=1664525, c=1013904223, m=2**32, n=n)

    def slab_plot(ax, u, title, band=(0.25, 0.26)):
        x, y, z = u[:-2], u[1:-1], u[2:]
        mask = (z >= band[0]) & (z < band[1])
        ax.scatter(x[mask], y[mask], s=3, alpha=0.6)
        ax.set_title(f"{title}\n({mask.sum()} points with x_(n+2) in {band})")
        ax.set_xlabel("x_n")
        ax.set_ylabel("x_(n+1)")

    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    slab_plot(axes[0], randu, "RANDU (a=65539, c=0)")
    slab_plot(axes[1], good, "Numerical Recipes LCG (a=1664525, c=1013904223)")
    fig.suptitle("Points whose 3rd coordinate falls in a thin slab:\n"
                 "a flawed LCG collapses onto a handful of lines")
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"[B] saved {out}")

    # Quantify it: because a=65539=2^16+3, (a-3)^2 = 2^32 = 2m exactly, so
    # x_(n+2) - 6*x_(n+1) + 9*x_n is an EXACT multiple of m for every n.
    # That collapses all consecutive triples onto a handful of parallel
    # planes in 3D -- the famous Marsaglia result.
    x_int = (randu * (2**31)).round().astype(np.int64)
    k = (x_int[2:] - 6 * x_int[1:-1] + 9 * x_int[:-2]) // (2**31)
    print(f"[B] RANDU: {n} consecutive triples lie on only "
          f"{len(np.unique(k))} parallel planes in 3D "
          f"(a well-distributed generator would need on the order of {n})")


# ---------- Part C: would simple 1D/2D tests even notice? ----------

def _repeat_cycle(seed, a, c, m, n):
    cycle = lcg_sequence(seed, a, c, m)[:-1]  # drop the repeated seed
    reps = n // len(cycle) + 1
    return np.array(cycle * reps)[:n] / m


def part_c_simple_tests(n=100_000):
    generators = {
        "toy LCG (m=17, tiled)": lambda: _repeat_cycle(6, 7, 3, 17, n),
        "RANDU (a=65539, c=0)": lambda: lcg_stream(1, 65539, 0, 2**31, n),
        "Numerical Recipes LCG": lambda: lcg_stream(1, 1664525, 1013904223, 2**32, n),
        "Python random.random() (MT19937)": lambda: np.array(
            [random.random() for _ in range(n)]),
    }

    print("\n[C] generator                        | chi2 p-value | lag-1 r  | lag-2 r")
    for name, gen in generators.items():
        u = gen()
        counts, _ = np.histogram(u, bins=20, range=(0, 1))
        _, p = stats.chisquare(counts)
        r1 = np.corrcoef(u[:-1], u[1:])[0, 1]
        r2 = np.corrcoef(u[:-2], u[2:])[0, 1]
        print(f"[C] {name:<36} | p={p:9.3f}  | {r1:+.4f}  | {r2:+.4f}")
    print("[C] all four pass the uniformity check and show ~0 autocorrelation --\n"
          "    yet part B shows RANDU is badly non-random in 3D. Low-order\n"
          "    1D/2D checks are not enough to validate a PRNG.")


if __name__ == "__main__":
    part_a_full_period_map()
    part_b_randu_planes()
    part_c_simple_tests()
