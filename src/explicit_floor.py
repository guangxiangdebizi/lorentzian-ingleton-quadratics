"""Exact regression for the spherical-ultrametric quadratic Ingleton proof.

This checks finite instances, not a replacement for the proof in
paper/main.tex. All construction, polarization, bottleneck and ratio checks use
fractions.Fraction. No floating eigensolver or internet dependency is used.
"""

import argparse
from decimal import Decimal, localcontext
from fractions import Fraction as F
import json
from pathlib import Path
import random


LEFT = (3, 5, 9, 6, 10)
RIGHT = (1, 2, 7, 11, 12)


def zeros(n):
    return [[F(0) for _ in range(n)] for _ in range(n)]


def product(values):
    result = F(1)
    for value in values:
        result *= value
    return result


def decimal(value):
    with localcontext() as context:
        context.prec = 20
        return str(Decimal(value.numerator) / Decimal(value.denominator))


def sphere_point(parameters):
    norm = sum((x * x for x in parameters), F(0))
    denominator = 1 + norm
    return [(1 - norm) / denominator] + [2 * x / denominator for x in parameters]


def sphere_matrix(points, weights):
    n = len(points)
    result = zeros(n)
    for i in range(n):
        assert sum(x * x for x in points[i]) == 1
        for j in range(i):
            inner = sum(a * b for a, b in zip(points[i], points[j]))
            result[i][j] = result[j][i] = weights[i] * weights[j] * (1 - inner)
    return result


def aggregate_and_polarize(matrix, groups, n):
    coefficients = zeros(n)
    for i in range(len(matrix)):
        for j in range(i):
            a, b = sorted((groups[i], groups[j]))
            coefficients[a][b] += matrix[i][j]
    result = zeros(2 * n)
    for i in range(n):
        result[2 * i][2 * i + 1] = result[2 * i + 1][2 * i] = coefficients[i][i]
        for j in range(i + 1, n):
            for a in (2 * i, 2 * i + 1):
                for b in (2 * j, 2 * j + 1):
                    result[a][b] = result[b][a] = coefficients[i][j] / 4

    # S_i=e_(i,0)+e_(i,1), D_i=e_(i,0)-e_(i,1): exact congruence blocks.
    for i in range(n):
        for j in range(n):
            symmetric = sum(result[2 * i + a][2 * j + b] for a in range(2) for b in range(2))
            expected = 2 * coefficients[i][i] if i == j else coefficients[min(i, j)][max(i, j)]
            assert symmetric == expected
            antisymmetric = sum((-1) ** (a + b) * result[2 * i + a][2 * j + b] for a in range(2) for b in range(2))
            assert antisymmetric == (-2 * coefficients[i][i] if i == j else 0)
            mixed = sum((-1) ** b * result[2 * i + a][2 * j + b] for a in range(2) for b in range(2))
            assert mixed == 0
    return result


def evaluate(matrix, values):
    return sum((matrix[i][j] * values[i] * values[j] for i in range(len(matrix)) for j in range(i)), F(0))


def check(matrix, masks, delta):
    active = [i for i, row in enumerate(matrix) if any(row)]
    matrix = [[matrix[i][j] for j in active] for i in active]
    masks = [masks[i] for i in active]
    n = len(matrix)
    assert n >= 2 and 0 < delta < 1
    s = [sum(row) for row in matrix]
    total = sum(s)
    distance2 = [[2 * total * matrix[i][j] / (s[i] * s[j]) for j in range(n)] for i in range(n)]
    bottleneck2 = [row[:] for row in distance2]
    for k in range(n):
        for i in range(n):
            for j in range(n):
                bottleneck2[i][j] = min(bottleneck2[i][j], max(bottleneck2[i][k], bottleneck2[k][j]))
    template = [[s[i] * s[j] * bottleneck2[i][j] / (2 * total) for j in range(n)] for i in range(n)]
    for i in range(n):
        for j in range(n):
            assert (template[i][j] == 0) == (matrix[i][j] == 0)
            assert template[i][j] <= matrix[i][j] <= (n - 1) ** 2 * template[i][j]
            for k in range(n):
                assert bottleneck2[i][j] <= max(bottleneck2[i][k], bottleneck2[k][j])

    m = {}
    q = {}
    bound = n * (n - 1) // 2 * (n - 1) ** 2
    for mask in set(LEFT + RIGHT):
        values = [F(1) if member & mask else delta for member in masks]
        m[mask] = max(template[i][j] * values[i] * values[j] for i in range(n) for j in range(i))
        q[mask] = evaluate(matrix, values)
        assert 0 < m[mask] <= q[mask] <= bound * m[mask]
    max_ratio = product(m[j] for j in LEFT) / product(m[j] for j in RIGHT)
    ratio = product(q[j] for j in LEFT) / product(q[j] for j in RIGHT)
    assert max_ratio >= 1
    assert ratio >= F(1, bound ** 5)
    return ratio, max_ratio, n


def run(seed, cases):
    rng = random.Random(seed)
    records = []
    minimum = F(100)
    maximum_n = 0
    for case in range(cases):
        n = 16 if case < 8 else rng.randrange(2, 9)
        count = n + rng.randrange(2, n + 2)
        dimension = rng.randrange(1, 5)
        points = []
        for i in range(count):
            if i and rng.randrange(5) == 0:
                points.append(points[rng.randrange(i)])
            else:
                points.append(sphere_point([F(rng.randrange(-9, 10), rng.randrange(1, 6)) for _ in range(dimension)]))
        if all(point == points[0] for point in points):
            points[-1] = [-x for x in points[0]]
        weights = [F(2) ** rng.randrange(-20, 21) for _ in range(count)]
        raw = sphere_matrix(points, weights)
        groups = [i % n for i in range(count)]
        rng.shuffle(groups)
        polarized = aggregate_and_polarize(raw, groups, n)
        atom_masks = list(range(16)) if n == 16 else [rng.randrange(16) for _ in range(n)]
        clone_masks = [mask for mask in atom_masks for _ in range(2)]
        epsilon = F(2) ** rng.randrange(-40, 41)
        delta = epsilon / (1 + epsilon)
        for mask in set(LEFT + RIGHT):
            atom_values = [F(1) if member & mask else delta for member in atom_masks]
            assert evaluate(raw, [atom_values[group] for group in groups]) == evaluate(polarized, [x for x in atom_values for _ in range(2)])
        ratio, max_ratio, active = check(polarized, clone_masks, delta)
        minimum = min(minimum, ratio)
        maximum_n = max(maximum_n, active)
        if case < 8:
            records.append({"case": case, "active_clones": active, "epsilon": str(epsilon), "ratio": decimal(ratio), "max_ratio": decimal(max_ratio)})

    # Exact PSD witness below one: a known five-vector family at a finite scale.
    scale = F(100)
    vectors = [
        (F(4, 3) / 100000, F(0)),
        (F(4, 3) / 100000, F(0)),
        (F(0), scale),
        (F(1, 1000), F(1, 2000)),
        (F(1000), F(1000)),
    ]
    matrix = zeros(5)
    for i in range(5):
        for j in range(i):
            determinant = vectors[i][0] * vectors[j][1] - vectors[i][1] * vectors[j][0]
            matrix[i][j] = matrix[j][i] = determinant ** 2
    epsilon = scale ** -8
    witness, witness_max, _ = check(matrix, [1, 2, 4, 8, 0], epsilon / (1 + epsilon))
    assert witness < 1
    return {
        "status": "PASS",
        "arithmetic": "exact fractions.Fraction; decimal strings for presentation only",
        "seed": seed,
        "random_cases": cases,
        "largest_active_clone_count": maximum_n,
        "all_16_membership_atoms_cases": min(cases, 8),
        "assertions": [
            "exact diagonal polarization recovery for every queried value",
            "exact symmetric, antisymmetric, mixed polarization congruence blocks",
            "exact zero support preservation",
            "coefficient distortion at most (N-1)^2",
            "bottleneck ultrametric inequalities",
            "max-product Ingleton for every random instance",
            "sum/max bounds and explicit final floor",
            "finite PSD instance has ratio strictly below one",
        ],
        "universal_constant": "476656^(-5)",
        "universal_constant_decimal": decimal(F(1, 476656 ** 5)),
        "minimum_random_ratio": decimal(minimum),
        "large_case_samples": records,
        "finite_psd_witness": {"ratio": decimal(witness), "max_ratio": decimal(witness_max), "epsilon": str(epsilon)},
        "scope_limit": "Finite regression is not a proof or an optimality certificate. See the complete symbolic proof.",
    }


def build_report():
    return run(20260905, 160)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=20260905)
    parser.add_argument("--cases", type=int, default=160)
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parents[1] / "reports" / "explicit_floor.json")
    args = parser.parse_args()
    report = run(args.seed, args.cases)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
