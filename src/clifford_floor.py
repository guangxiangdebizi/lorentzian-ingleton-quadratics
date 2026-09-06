"""Exact, standard-library Clifford checks; not a substitute for the proof."""

import argparse
from fractions import Fraction as F
import hashlib
import itertools
import json
import math
from pathlib import Path
import random


QUERIES = (3, 5, 9, 6, 10, 1, 2, 7, 11, 12)
SEED = 2026090603


def eye(n):
    return [[F(i == j) for j in range(n)] for i in range(n)]


def zero(n):
    return [[F(0) for _ in range(n)] for _ in range(n)]


def linear_sum(matrices, weights):
    n = len(matrices[0])
    return [[sum(w * a[i][j] for w, a in zip(weights, matrices))
             for j in range(n)] for i in range(n)]


def multiply(a, b):
    n = len(a)
    return [[sum(a[i][k] * b[k][j] for k in range(n))
             for j in range(n)] for i in range(n)]


def determinant(matrix):
    a = [[F(x) for x in row] for row in matrix]
    answer = F(1)
    for k in range(len(a)):
        p = next((i for i in range(k, len(a)) if a[i][k]), None)
        if p is None:
            return F(0)
        if p != k:
            a[k], a[p] = a[p], a[k]
            answer = -answer
        pivot = a[k][k]
        answer *= pivot
        for i in range(k + 1, len(a)):
            factor = a[i][k] / pivot
            for j in range(k + 1, len(a)):
                a[i][j] -= factor * a[k][j]
            a[i][k] = F(0)
    return answer


def rank(matrix):
    a = [[F(x) for x in row] for row in matrix]
    row = 0
    for col in range(len(a[0])):
        p = next((i for i in range(row, len(a)) if a[i][col]), None)
        if p is None:
            continue
        a[row], a[p] = a[p], a[row]
        pivot = a[row][col]
        a[row] = [v / pivot for v in a[row]]
        for i in range(row + 1, len(a)):
            factor = a[i][col]
            a[i] = [x - factor * y for x, y in zip(a[i], a[row])]
        row += 1
        if row == len(a):
            break
    return row


def clifford(r, unsigned=False):
    n = 1 << r
    matrices = []
    for j in range(r):
        a = zero(n)
        for subset in range(n):
            parity = (subset & ((1 << j) - 1)).bit_count() % 2
            a[subset ^ (1 << j)][subset] = F(1 if unsigned else (-1) ** parity)
        matrices.append(a)
    return matrices


def check_clifford(matrices):
    n = len(matrices[0])
    for a in matrices:
        if a != [list(row) for row in zip(*a)]:
            raise AssertionError("Clifford generator is not symmetric")
        if sum(a[i][i] for i in range(n)) != 0:
            raise AssertionError("Clifford generator is not traceless")
        if multiply(a, a) != eye(n):
            raise AssertionError("Clifford generator does not square to I")
    for j, a in enumerate(matrices):
        for b in matrices[j + 1:]:
            if linear_sum([multiply(a, b), multiply(b, a)], [1, 1]) != zero(n):
                raise AssertionError("Clifford anticommutator is nonzero")


def polynomial_multiply(a, b):
    out = {}
    for ka, va in a.items():
        for kb, vb in b.items():
            key = tuple(x + y for x, y in zip(ka, kb))
            out[key] = out.get(key, 0) + va * vb
    return {k: v for k, v in out.items() if v}


def symbolic_determinant(r):
    """Subset expansion of det(t I + sum x_j Gamma_j), exactly in Z[t,x]."""
    n = 1 << r
    matrices = clifford(r)
    degree_zero = (0,) * (r + 1)
    states = {0: {degree_zero: 1}}
    for row in range(n):
        following = {}
        for mask, polynomial in states.items():
            entries = [(row, 0, 1)]
            entries.extend((row ^ (1 << j), j + 1,
                            int(matrices[j][row][row ^ (1 << j)])) for j in range(r))
            for col, variable, coefficient in entries:
                if mask & (1 << col):
                    continue
                sign = (-1) ** (mask >> (col + 1)).bit_count()
                target = following.setdefault(mask | (1 << col), {})
                for monomial, value in polynomial.items():
                    powers = list(monomial)
                    powers[variable] += 1
                    key = tuple(powers)
                    target[key] = target.get(key, 0) + sign * coefficient * value
        states = {mask: {k: v for k, v in poly.items() if v}
                  for mask, poly in following.items()}
    actual = states[(1 << n) - 1]
    quadratic = {}
    for variable in range(r + 1):
        powers = [0] * (r + 1)
        powers[variable] = 2
        quadratic[tuple(powers)] = 1 if variable == 0 else -1
    expected = {degree_zero: 1}
    for _ in range(n // 2):
        expected = polynomial_multiply(expected, quadratic)
    if actual != expected:
        raise AssertionError("Exact multivariate Clifford determinant identity failed")
    return {"r": r, "matrix_size": n, "monomial_count": len(actual),
            "identity": "det(t I + sum x_j Gamma_j) = (t^2 - sum x_j^2)^(2^(r-1))"}


def examples():
    yield "rank_zero", [F(1), F(2), F(0)], [[], [], []]
    for r in range(1, 5):
        weights, vectors = [], []
        for j in range(r):
            v = [F(i == j) for i in range(r)]
            vectors.extend([v, [-x for x in v]])
            weights.extend([F(1), F(1)])
        vectors.extend([[F(0)] * r, [F(0)] * r])
        weights.extend([F(0), F(2)])
        yield "axis_boundary_r%d" % r, weights, vectors
    rng = random.Random(SEED)
    for case in range(12):
        r = 1 + case % 3
        weights, vectors = [], []
        for j in range(r):
            v = [F(i == j) for i in range(r)]
            vectors.extend([v, [-x for x in v]])
            weights.extend([F(1), F(1)])
        for _ in range(4):
            v = [F(rng.randrange(-5, 6), rng.randrange(1, 6)) for _ in range(r)]
            a = sum(abs(x) for x in v)
            vectors.extend([v, [-x for x in v]])
            weights.extend([a, a])
        vectors.extend([[F(0)] * r, [F(0)] * r])
        weights.extend([F(0), F(2)])
        yield "rational_r%d_case%d" % (r, case), weights, vectors


def encode(value):
    return {"numerator": str(value.numerator), "denominator": str(value.denominator)}


def quotient(values):
    return math.prod(values[:5]) / math.prod(values[5:])


def check_example(name, weights, vectors):
    input_rank = len(vectors[0])
    r = max(1, input_rank)
    if input_rank == 0:
        vectors = [[F(0)] for _ in vectors]
    matrices = clifford(r)
    n, size = len(weights), 1 << r
    gram = [[sum(x * y for x, y in zip(v, w)) for w in vectors] for v in vectors]
    h = [[weights[i] * weights[j] - gram[i][j] for j in range(n)] for i in range(n)]
    if not all(x >= 0 for row in h for x in row):
        raise AssertionError("Hessian coefficient nonnegativity failed")
    s = [sum(row) for row in h]
    total = sum(s)
    a_sum = sum(weights)
    if total <= 0 or total != a_sum ** 2:
        raise AssertionError("T positivity or canonical square root failed")
    if [x / a_sum for x in s] != weights:
        raise AssertionError("Canonical a_i reconstruction failed")
    k = [[s[i] * s[j] / total - h[i][j] for j in range(n)] for i in range(n)]
    if k != gram or any(sum(row) for row in k) or rank(k) != input_rank:
        raise AssertionError("Canonical Gram matrix reconstruction/rank/null vector failed")
    atoms = [linear_sum([eye(size)] + matrices, [a] + v) for a, v in zip(weights, vectors)]
    if linear_sum(atoms, [1] * n) != linear_sum([eye(size)], [a_sum]):
        raise AssertionError("Common strictly positive sum failed")
    for a, v, p in zip(weights, vectors, atoms):
        norm2 = sum(x * x for x in v)
        if a < 0 or a * a < norm2:
            raise AssertionError("PSD atom eigenvalue certificate failed")
        identity = linear_sum([multiply(p, p), p, eye(size)], [1, -2 * a, a * a - norm2])
        if identity != zero(size):
            raise AssertionError("PSD atom quadratic minimal-polynomial identity failed")
    noncommuting = any(multiply(a, b) != multiply(b, a)
                       for i, a in enumerate(atoms) for b in atoms[i + 1:])
    if input_rank >= 2 and not noncommuting:
        raise AssertionError("Intended noncommuting example became commuting")
    def value(z):
        return sum(z[i] * h[i][j] * z[j] for i in range(n) for j in range(n))
    points = [[F(0)] * n] + [[F((i * 3 + shift) % 7 - 3) for i in range(n)] for shift in range(4)]
    for z in points:
        if determinant(linear_sum(atoms, z)) != value(z) ** (size // 2):
            raise AssertionError("Power identity failed at an unrestricted real point")
    records = []
    masks_by_mode = {
        "disjoint": [0 if i % 5 == 4 else 1 << (i % 5) for i in range(n)],
        "overlap": [i % 16 for i in range(n)],
    }
    for mode, masks in masks_by_mode.items():
        for epsilon in (F(1, 10 ** 6), F(1, 7), F(10 ** 6)):
            q_values, f_values = [], []
            for query in QUERIES:
                z = [epsilon + int(bool(mask & query)) for mask in masks]
                qv = value(z)
                fv = determinant(linear_sum(atoms, z))
                if qv <= 0 or fv != qv ** (size // 2):
                    raise AssertionError("Common-regularizer determinant power failed")
                q_values.append(qv)
                f_values.append(fv)
            rq, rf = quotient(q_values), quotient(f_values)
            if rf != rq ** (size // 2):
                raise AssertionError("Five-by-five scalar cancellation failed")
            base = 4 if mode == "disjoint" else 128
            if rf < F(1, base ** size) or rq < F(1, base ** 2):
                raise AssertionError("PSD or transferred quadratic floor failed")
            records.append({"mode": mode, "epsilon": encode(epsilon), "quadratic_ratio": encode(rq)})
    return {"name": name, "variables": n, "rank_K": input_rank, "matrix_size": size,
            "noncommuting_atoms": noncommuting, "zero_atoms": weights.count(F(0)),
            "boundary_atoms": sum(a > 0 and a * a == sum(x * x for x in v)
                                  for a, v in zip(weights, vectors)),
            "unrestricted_real_power_checks": len(points), "regularized_cases": records}


def negative_controls():
    broken = clifford(2, unsigned=True)
    rejected = False
    try:
        check_clifford(broken)
    except AssertionError:
        rejected = True
    if not rejected:
        raise AssertionError("Unsigned-generator negative control escaped")
    bad_determinant = determinant(linear_sum([eye(4)] + broken, [3, 1, 1]))
    if bad_determinant == (3 ** 2 - 1 ** 2 - 1 ** 2) ** 2:
        raise AssertionError("Broken power-identity negative control escaped")
    two_positive_h = [[F(1), F(0)], [F(0), F(1)]]
    s = [sum(row) for row in two_positive_h]
    bad_k = [[s[i] * s[j] / sum(s) - two_positive_h[i][j] for j in range(2)] for i in range(2)]
    if not any(bad_k[i][i] < 0 for i in range(2)):
        raise AssertionError("Two-positive-directions negative control escaped")
    h = [[F(-1), F(2)], [F(2), F(1)]]
    s = [sum(row) for row in h]
    total = sum(s)
    if s[0] ** 2 / total >= s[0] ** 2 / total - h[0][0]:
        raise AssertionError("Negative diagonal atom obstruction escaped")
    return {"unsigned_generator_rejected": True,
            "unsigned_generator_determinant": encode(bad_determinant),
            "correct_identity_determinant": 49,
            "two_positive_hessian_rejected": True,
            "negative_diagonal_atom_obstruction_detected": True}


def noncentered_canonical_checks():
    inputs = [
        ([2, 3, 0], [[1], [0], [0]]),
        ([2, 2, 2, 0], [[1, 0], [0, 1], [1, 1], [0, 0]]),
        ([4, 4, 4, 4, 0], [[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 2, 3], [0, 0, 0]]),
    ]
    records = []
    for weights, vectors in inputs:
        n = len(weights)
        h = [[F(weights[i] * weights[j] - sum(x * y for x, y in zip(vectors[i], vectors[j])))
              for j in range(n)] for i in range(n)]
        s = [sum(row) for row in h]
        total = sum(s)
        k = [[s[i] * s[j] / total - h[i][j] for j in range(n)] for i in range(n)]
        if any(sum(row) for row in k):
            raise AssertionError("Noncentered input failed canonical centering")
        minor_count = 0
        for count in range(1, n + 1):
            for indices in itertools.combinations(range(n), count):
                minor = [[k[i][j] for j in indices] for i in indices]
                if determinant(minor) < 0:
                    raise AssertionError("Canonical K is not PSD on a noncentered input")
                minor_count += 1
        for shift in range(5):
            z = [F((i * 3 + shift) % 7 - 3) for i in range(n)]
            original = sum(z[i] * h[i][j] * z[j] for i in range(n) for j in range(n))
            canonical = sum(x * y for x, y in zip(s, z)) ** 2 / total
            canonical -= sum(z[i] * k[i][j] * z[j] for i in range(n) for j in range(n))
            if original != canonical:
                raise AssertionError("Noncentered Lorentz representation identity failed")
        records.append({"variables": n, "rank_K": rank(k), "T": encode(total),
                        "all_principal_minors_checked": minor_count,
                        "unrestricted_real_reconstruction_checks": 5})
    return records


def build_report():
    for r in range(1, 5):
        check_clifford(clifford(r))
    symbolic = [symbolic_determinant(r) for r in range(1, 4)]
    records = [check_example(*example) for example in examples()]
    return {"schema": "clifford-quadratic-audit-v1", "status": "PASS", "seed": SEED,
            "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "scope": "Exact finite algebra checks; the universal proof is in paper/main.tex.",
            "quadratic_examples": len(records), "regularized_cases": sum(len(x["regularized_cases"]) for x in records),
            "symbolic_identities": symbolic, "negative_controls": negative_controls(),
            "noncentered_canonical_checks": noncentered_canonical_checks(), "records": records}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check", type=Path)
    args = parser.parse_args()
    report = build_report()
    encoded = (json.dumps(report, indent=2, sort_keys=True) + "\n").encode("utf-8")
    if args.check and args.check.read_bytes() != encoded:
        raise SystemExit("Certificate mismatch")
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_bytes(encoded)
    print(json.dumps({k: report[k] for k in ("status", "quadratic_examples", "regularized_cases")}))


if __name__ == "__main__":
    main()
