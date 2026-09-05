"""Independent exact checks of the quadratic floor candidate.

No imports from the candidate implementation. The two principal checks are
enumeration of input matrices (not generated Gram instances) and computation
of lattice-quotient dimensions in a truncated Laurent coefficient space.
"""

from fractions import Fraction as F
from itertools import combinations, product
import json
from pathlib import Path
import random


LEFT = (3, 5, 9, 6, 10)
RIGHT = (1, 2, 7, 11, 12)


def psd_exact(matrix):
    a = [list(map(F, row)) for row in matrix]
    while a:
        if any(a[i][i] < 0 for i in range(len(a))):
            return False
        pivots = [i for i in range(len(a)) if a[i][i] > 0]
        if not pivots:
            return not any(any(row) for row in a)
        p = pivots[0]
        rest = [i for i in range(len(a)) if i != p]
        a = [[a[i][j] - a[i][p] * a[p][j] / a[p][p]
              for j in rest] for i in rest]
    return True


def certify_and_template(matrix):
    n0 = len(matrix)
    assert all(len(row) == n0 for row in matrix)
    assert all(matrix[i][i] == 0 for i in range(n0))
    assert all(matrix[i][j] == matrix[j][i] and matrix[i][j] >= 0
               for i in range(n0) for j in range(n0))
    active = [i for i in range(n0) if any(matrix[i])]
    if not active:
        return None
    b = [[F(matrix[i][j]) for j in active] for i in active]
    n = len(b)
    s = [sum(row) for row in b]
    total = sum(s)
    k = [[s[i] * s[j] / total - b[i][j] for j in range(n)]
         for i in range(n)]
    if not psd_exact(k):
        return None
    d2 = [[2 * total * b[i][j] / (s[i] * s[j]) for j in range(n)]
          for i in range(n)]
    for i, j, k0 in product(range(n), repeat=3):
        excess = d2[i][j] - d2[i][k0] - d2[k0][j]
        assert excess <= 0 or excess ** 2 <= 4 * d2[i][k0] * d2[k0][j]

    # Kruskal plus unique tree paths, deliberately not Floyd-Warshall.
    parent = list(range(n))

    def root(x):
        while parent[x] != x:
            x = parent[x]
        return x

    adjacency = [[] for _ in range(n)]
    for weight, i, j in sorted((d2[i][j], i, j) for i, j in combinations(range(n), 2)):
        if root(i) == root(j):
            continue
        parent[root(i)] = root(j)
        adjacency[i].append((j, weight))
        adjacency[j].append((i, weight))
    u2 = [[F(0) for _ in range(n)] for _ in range(n)]
    for source in range(n):
        stack = [(source, -1, F(0))]
        while stack:
            vertex, predecessor, peak = stack.pop()
            u2[source][vertex] = peak
            stack.extend((neighbor, vertex, max(peak, weight))
                         for neighbor, weight in adjacency[vertex]
                         if neighbor != predecessor)
    c = [[s[i] * s[j] * u2[i][j] / (2 * total) for j in range(n)]
         for i in range(n)]
    for i, j in product(range(n), repeat=2):
        assert (b[i][j] == 0) == (c[i][j] == 0)
        assert c[i][j] <= b[i][j] <= (n - 1) ** 2 * c[i][j]
    for i, j, k0 in product(range(n), repeat=3):
        assert u2[i][j] <= max(u2[i][k0], u2[k0][j])
    for i, j, k0, ell in combinations(range(n), 4):
        products = [c[i][j] * c[k0][ell], c[i][k0] * c[j][ell], c[i][ell] * c[j][k0]]
        assert products.count(max(products)) >= 2
    return b, c, active


def evaluate_ratio(matrix, masks, delta, maximum=False):
    n = len(matrix)
    values = {}
    for query in set(LEFT + RIGHT):
        terms = [matrix[i][j] * (1 if masks[i] & query else delta)
                 * (1 if masks[j] & query else delta)
                 for i, j in combinations(range(n), 2)]
        values[query] = max(terms) if maximum else sum(terms)
    result = F(1)
    for query in LEFT:
        result *= values[query]
    for query in RIGHT:
        result /= values[query]
    return result


def audit_integer_inputs():
    edges = list(combinations(range(4), 2))
    accepted = rejected = loop_cases = parallel_cases = 0
    rng = random.Random(731209)
    for weights in product(range(4), repeat=6):
        matrix = [[F(0) for _ in range(4)] for _ in range(4)]
        for (i, j), weight in zip(edges, weights):
            matrix[i][j] = matrix[j][i] = F(weight)
        answer = certify_and_template(matrix)
        if answer is None:
            rejected += 1
            continue
        accepted += 1
        b, c, active = answer
        n = len(active)
        loop_cases += n < 4
        parallel_cases += any(b[i][j] == 0 for i, j in combinations(range(n), 2))
        for delta in (F(1, 257), F(1, 2), F(256, 257)):
            masks = [rng.randrange(16) for _ in active]
            assert evaluate_ratio(c, masks, delta, maximum=True) >= 1
            assert evaluate_ratio(b, masks, delta) >= F(1, (n * (n - 1) // 2 * (n - 1) ** 2) ** 5)
    # Two disjoint positive edges have two positive eigenvalues.
    disconnected = [[0, 1, 0, 0], [1, 0, 0, 0],
                    [0, 0, 0, 1], [0, 0, 1, 0]]
    assert certify_and_template(disconnected) is None
    assert not psd_exact([[0, 1], [1, 0]])
    return {"enumerated_matrices": 4 ** 6, "accepted": accepted,
            "rejected_including_zero": rejected, "accepted_with_loops": loop_cases,
            "accepted_with_parallel_pairs": parallel_cases,
            "query_checks": accepted * 3,
            "negative_controls": "two-positive matrix and zero-pivot indefinite matrix rejected"}


def audit_two_clones():
    b, c, _ = certify_and_template([[0, 7], [7, 0]])
    assert b == c
    count = 0
    for masks in product(range(16), repeat=2):
        for epsilon in (F(1, 2 ** 40), F(1), F(2 ** 40)):
            delta = epsilon / (1 + epsilon)
            assert evaluate_ratio(b, masks, delta) >= 1
            count += 1
    return {"all_memberships_and_regularizers": count, "constant": 1}


def poly_add(a, b, sign=1):
    answer = a.copy()
    for exponent, value in b.items():
        answer[exponent] = answer.get(exponent, F(0)) + sign * value
        if answer[exponent] == 0:
            del answer[exponent]
    return answer


def poly_mul(a, b):
    answer = {}
    for x, vx in a.items():
        for y, vy in b.items():
            answer[x + y] = answer.get(x + y, F(0)) + vx * vy
    return {key: value for key, value in answer.items() if value}


def shift(poly, exponent):
    return {key + exponent: value for key, value in poly.items()}


def determinant(a, b):
    return poly_add(poly_mul(a[0], b[1]), poly_mul(a[1], b[0]), -1)


def sparse_rank(vectors):
    basis = {}
    for original in vectors:
        vector = original.copy()
        while vector:
            pivot = min(vector)
            coefficient = vector[pivot]
            if pivot not in basis:
                basis[pivot] = {key: value / coefficient for key, value in vector.items()}
                break
            for key, value in basis[pivot].items():
                vector[key] = vector.get(key, F(0)) - coefficient * value
                if vector[key] == 0:
                    del vector[key]
    return len(basis)


def truncated_module_rank(columns, lower, upper):
    generators = []
    width = upper - lower
    for column in columns:
        least = min(exponent for poly in column for exponent in poly)
        for power in range(max(0, upper - least)):
            vector = {}
            for coordinate, poly in enumerate(column):
                for exponent, value in poly.items():
                    degree = exponent + power
                    if lower <= degree < upper:
                        vector[coordinate * width + degree - lower] = value
            if vector:
                generators.append(vector)
    return sparse_rank(generators)


def audit_laurent_modules():
    rng = random.Random(916024)
    cases = 24
    quotients = 0
    negative_valuations = 0
    parallel_cases = 0
    for case in range(cases):
        # Root with two binary children; the fifth leaf repeats leaf 0.
        root = rng.randrange(-4, 2)
        child_left = root + rng.randrange(1, 4)
        child_right = root + rng.randrange(1, 4)
        x = [{root: F(1), child_left: F(1)},
             {root: F(1), child_left: F(3)},
             {root: F(4), child_right: F(2)},
             {root: F(4), child_right: F(5)}]
        if case % 2 == 0:
            x.append(x[0].copy())
            parallel_cases += 1
        alpha = [rng.randrange(-3, 4) for _ in x]
        columns = [({a: F(1)}, shift(poly, a)) for a, poly in zip(alpha, x)]
        n = len(columns)
        p = {}
        for i, j in combinations(range(n), 2):
            det = determinant(columns[i], columns[j])
            if not det:
                assert x[i] == x[j]
                continue
            h = min(poly_add(x[i], x[j], -1))
            p[i, j] = min(det)
            assert p[i, j] == alpha[i] + alpha[j] + h
            negative_valuations += p[i, j] < 0
        selected = next(iter(p))
        valdet = p[selected]
        adj_min = min(exponent for i in selected for poly in columns[i] for exponent in poly)
        upper = valdet - adj_min
        regularizer = rng.randrange(1, 5)
        lower = min(exponent for column in columns for poly in column for exponent in poly) - regularizer
        assert upper > lower
        base_rank = truncated_module_rank(columns, lower, upper)
        r0 = max(-value for value in p.values())
        masks = [rng.randrange(16) for _ in columns]
        lengths = {}
        for query in range(16):
            included = [i for i in range(n) if masks[i] & query]
            extra = [tuple(shift(poly, -regularizer) for poly in columns[i]) for i in included]
            measured = truncated_module_rank(columns + extra, lower, upper) - base_rank
            predicted = max(regularizer * ((i in included) + (j in included)) - value
                            for (i, j), value in p.items()) - r0
            assert measured == predicted
            assert 0 <= measured <= 2 * regularizer
            lengths[query] = measured
            quotients += 1
        assert sum(lengths[q] for q in LEFT) >= sum(lengths[q] for q in RIGHT)
        ambient_columns = [tuple(shift(poly, -regularizer) for poly in column) for column in columns]
        assert truncated_module_rank(ambient_columns, lower, upper) - base_rank == 2 * regularizer
    return {"tree_cases": cases, "direct_quotient_dimension_checks": quotients,
            "parallel_class_cases": parallel_cases, "negative_determinant_valuations": negative_valuations,
            "ambient_dimension_checks": cases,
            "method": "exact rational Gaussian rank in one finite Laurent coefficient space"}


def build_report():
    return {"status": "PASS", "arithmetic": "fractions.Fraction only",
              "independence": "no imports or generated instances from the candidate tests",
              "integer_inputs": audit_integer_inputs(),
              "two_clone_boundary": audit_two_clones(),
              "laurent_modules": audit_laurent_modules(),
              "limit": "Finite checks supplement the symbolic audit; they do not prove a universal theorem."}


def main():
    report = build_report()
    output = Path(__file__).resolve().parents[1] / "reports" / "hostile_audit.json"
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
