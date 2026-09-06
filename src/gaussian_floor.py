"""Exact Schur-complement audit of the continuous auxiliary-noise proof."""

import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
from math import prod
from pathlib import Path
import random


def det(matrix):
    a = [[F(x) for x in row] for row in matrix]
    result = F(1)
    for j in range(len(a)):
        pivot = next((i for i in range(j, len(a)) if a[i][j]), None)
        if pivot is None:
            return F(0)
        if pivot != j:
            a[j], a[pivot] = a[pivot], a[j]
            result = -result
        value = a[j][j]
        result *= value
        for i in range(j + 1, len(a)):
            factor = a[i][j] / value
            for k in range(j + 1, len(a)):
                a[i][k] -= factor * a[j][k]
    return result


def inverse(matrix):
    n = len(matrix)
    a = [[F(x) for x in row] + [F(i == j) for j in range(n)] for i, row in enumerate(matrix)]
    for j in range(n):
        pivot = next((i for i in range(j, n) if a[i][j]), None)
        if pivot is None:
            raise ValueError("singular inverse")
        a[j], a[pivot] = a[pivot], a[j]
        value = a[j][j]
        a[j] = [x / value for x in a[j]]
        for i in range(n):
            if i != j:
                value = a[i][j]
                a[i] = [x - value * y for x, y in zip(a[i], a[j])]
    return [row[n:] for row in a]


def audit(sigma, dims):
    a, b, _, _ = dims
    n = len(sigma)
    ib = inverse([row[a:a + b] for row in sigma[a:a + b]])
    k = [[sigma[i][j] - sum(sigma[i][a + u] * ib[u][v] * sigma[a + v][j]
                              for u in range(b) for v in range(b))
          for j in range(a)] for i in range(a)]
    joint = [row[:] + row[:a] for row in sigma]
    joint += [sigma[i][:] + [sigma[i][j] + k[i][j] for j in range(a)] for i in range(a)]
    if det(sigma) <= 0 or det(k) <= 0 or det(joint) <= 0:
        raise AssertionError("nondegenerate Gaussian gate failed")
    labels = [1 << block for block, size in enumerate(dims + (a,)) for _ in range(size)]

    @lru_cache(None)
    def minor(mask):
        indices = [i for i, label in enumerate(labels) if mask & label]
        return det([[joint[i][j] for j in indices] for i in indices])

    def conditional(x, z):
        return minor(x | z) / minor(z)

    def information(x, y, z=0):
        return minor(x | z) * minor(y | z) / (minor(z) * minor(x | y | z))

    # These rational quantities are exp(2I), not floating logarithms.
    for z in (4, 8):
        assert information(1, 2, z) >= information(16, 2, z)
        assert conditional(16, 2 | z) <= conditional(16, 2)
    lhs = conditional(16, 4) * conditional(16, 8) * information(4, 8)
    rhs = minor(16) * conditional(16, 12) * information(4, 8, 16)
    assert lhs == rhs
    assert information(4, 8, 16) >= 1
    assert conditional(16, 2) == 2**a * det(k)
    assert minor(16) >= minor(1)
    assert conditional(16, 12) >= det(k)
    ratio = prod(minor(j) for j in (3, 5, 9, 6, 10)) / prod(minor(j) for j in (1, 2, 7, 11, 12))
    auxiliary = minor(16) * conditional(16, 12) / (conditional(16, 2)**2 * information(1, 2))
    assert ratio >= auxiliary >= F(1, 4**a)
    singular = [row[:] for row in joint]
    for i in range(a):
        for j in range(a):
            singular[n + i][n + j] -= k[i][j]
    assert det(singular) == 0  # Omitting independent noise is not an eligible Gaussian.
    return {"dims": list(dims), "ratio": str(ratio), "auxiliary_lower": str(auxiliary),
            "floor": str(F(1, 4**a)), "strict_conditional_entropy_decrease":
            conditional(16, 6) < conditional(16, 2)}


def nonunit_fixture():
    vectors = [(F(4, 300000), F(0)), (F(4, 300000), F(0)),
               (F(0), F(100)), (F(1, 1000), F(1, 2000)), (F(1000), F(1000))]
    s = [[sum(v[i] * v[j] for v in vectors) for j in range(2)] for i in range(2)]
    inv = inverse(s)
    epsilon = F(1, 100**8)
    return [[F(i == j) + sum(vectors[i][r] * inv[r][t] * vectors[j][t]
                              for r in range(2) for t in range(2)) / epsilon
             for j in range(4)] for i in range(4)]


def build_report():
    rng = random.Random(2026090602)
    records = []
    for dims in ((1, 1, 1, 1), (2, 1, 2, 1), (2, 2, 2, 2), (3, 2, 1, 2)):
        n = sum(dims)
        for case in range(10):
            g = [[F(rng.randrange(-3, 4)) * F(2)**rng.randrange(-5, 6)
                  for _ in range(n)] for _ in range(n)]
            noise = [F(2)**rng.randrange(-12, 13) for _ in range(n)]
            sigma = [[sum(g[i][r] * g[j][r] for r in range(n))
                      + (noise[i] if i == j else 0) for j in range(n)] for i in range(n)]
            records.append(audit(sigma, dims))
    fixture = audit(nonunit_fixture(), (1, 1, 1, 1))
    assert F(fixture["ratio"]) < 1
    assert any(row["strict_conditional_entropy_decrease"] for row in records)
    return {"schema": "gaussian-auxiliary-exact-v1", "status": "PASS", "seed": 2026090602,
            "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "random_cases": len(records), "records": records, "nonunit_fixture": fixture,
            "negative_controls": ["independent-noise omission gives singular joint covariance",
                                  "reversed conditioning inequality fails on strict cases",
                                  "unit determinant ratio fails on the exact five-vector fixture"],
            "scope": "Exact finite formula and sign regression; not the universal proof or a sharpness claim."}


def main():
    if not __debug__:
        raise SystemExit("EXACT_AUDIT_REQUIRES_ASSERTIONS")
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--output", type=Path)
    mode.add_argument("--check", type=Path)
    args = parser.parse_args()
    report = build_report()
    data = (json.dumps(report, indent=2, sort_keys=True) + "\n").encode("ascii")
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_bytes(data)
    elif args.check.read_bytes() != data:
        raise SystemExit("GAUSSIAN_AUXILIARY_REPLAY_MISMATCH")
    print("GAUSSIAN_AUXILIARY_EXACT_PASS", report["random_cases"], "plus nonunit fixture")


if __name__ == "__main__":
    main()
