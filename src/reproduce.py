#!/usr/bin/env python3
"""Exact certificates for the Lorentzian-quadratic Ingleton note.

This file deliberately uses only the Python standard library.  The finite
valuated-matroid enumeration is a convention regression, not the proof of the
general theorem.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Iterable, Mapping, Sequence


ROOT = Path(__file__).resolve().parents[1]
EXPECTED_VALUATED_HASH = (
    "9eb0e7e8d5edf3cde90e21556308bb5475e8039243bb779a2f5233cb1b8a8321"
)
EXPECTED_SEPARATION_HASH = (
    "860c9ae48e6e9b616771ec9d22a926e44f1d3c9822a1b883720c60b1c130b1df"
)
EXPECTED_FINITE_WITNESS = Fraction(
    1151375679760487544161993156626743807199200025995873,
    1309813763635875854527900038158486602837796044630625,
)

Polynomial = dict[tuple[int, ...], Fraction]
Univariate = list[Fraction]

INGLETON_LEFT = (("A", "B"), ("A", "C"), ("A", "D"), ("B", "C"), ("B", "D"))
INGLETON_RIGHT = (
    ("A",),
    ("B",),
    ("A", "B", "C"),
    ("A", "B", "D"),
    ("C", "D"),
)


def _trim(polynomial: Univariate) -> Univariate:
    while len(polynomial) > 1 and polynomial[-1] == 0:
        polynomial.pop()
    return polynomial


def multiply(first: Univariate, second: Univariate) -> Univariate:
    product = [Fraction(0) for _ in range(len(first) + len(second) - 1)]
    for i, left in enumerate(first):
        for j, right in enumerate(second):
            product[i + j] += left * right
    return _trim(product)


def evaluate_univariate(polynomial: Univariate, value: Fraction) -> Fraction:
    result = Fraction(0)
    for coefficient in reversed(polynomial):
        result = result * value + coefficient
    return result


def regularized_evaluation(
    polynomial: Mapping[tuple[int, ...], Fraction], subset: Iterable[int]
) -> Univariate:
    """Return f(1_S + epsilon 1_E) as a polynomial in epsilon."""
    chosen = frozenset(subset)
    degree = sum(next(iter(polynomial)))
    result = [Fraction(0) for _ in range(degree + 1)]
    for exponent, coefficient in polynomial.items():
        inside_degree = sum(exponent[index] for index in chosen)
        outside_degree = degree - inside_degree
        for added_degree in range(inside_degree + 1):
            result[outside_degree + added_degree] += coefficient * math.comb(
                inside_degree, added_degree
            )
    return _trim(result)


def _union_groups(
    groups: Mapping[str, frozenset[int]], names: Sequence[str]
) -> frozenset[int]:
    return frozenset().union(*(groups[name] for name in names))


def ratio_polynomials(
    polynomial: Mapping[tuple[int, ...], Fraction],
    groups: Mapping[str, frozenset[int]],
) -> tuple[Univariate, Univariate]:
    numerator = [Fraction(1)]
    denominator = [Fraction(1)]
    for names in INGLETON_LEFT:
        numerator = multiply(
            numerator,
            regularized_evaluation(polynomial, _union_groups(groups, names)),
        )
    for names in INGLETON_RIGHT:
        denominator = multiply(
            denominator,
            regularized_evaluation(polynomial, _union_groups(groups, names)),
        )
    return numerator, denominator


def ratio_at(
    polynomial: Mapping[tuple[int, ...], Fraction],
    groups: Mapping[str, frozenset[int]],
    epsilon: Fraction,
) -> Fraction:
    numerator, denominator = ratio_polynomials(polynomial, groups)
    return evaluate_univariate(numerator, epsilon) / evaluate_univariate(
        denominator, epsilon
    )


def five_vector_polynomial(x_value: Fraction, scale: int) -> Polynomial:
    """Cauchy--Binet coefficients of the five-vector PSD pencil."""
    x = Fraction(x_value)
    ell = Fraction(scale)
    return {
        (1, 0, 1, 0, 0): x * ell**-3,
        (1, 0, 0, 1, 0): x * Fraction(1, 4) * ell**-8,
        (1, 0, 0, 0, 1): x * ell**-2,
        (0, 1, 1, 0, 0): x * ell**-3,
        (0, 1, 0, 1, 0): x * Fraction(1, 4) * ell**-8,
        (0, 1, 0, 0, 1): x * ell**-2,
        (0, 0, 1, 1, 0): ell**-1,
        (0, 0, 1, 0, 1): ell**5,
        (0, 0, 0, 1, 1): Fraction(1, 4),
    }


def five_vector_limit(x_value: Fraction) -> Fraction:
    x = Fraction(x_value)
    return (x + 1) ** 4 / (2 * x * (2 * x + 1) ** 2)


@dataclass(frozen=True)
class Qsqrt17:
    """An exact a + b*sqrt(17) element."""

    a: Fraction
    b: Fraction = Fraction(0)

    @staticmethod
    def coerce(value: int | Fraction | "Qsqrt17") -> "Qsqrt17":
        if isinstance(value, Qsqrt17):
            return value
        return Qsqrt17(Fraction(value))

    def __add__(self, other: int | Fraction | "Qsqrt17") -> "Qsqrt17":
        rhs = self.coerce(other)
        return Qsqrt17(self.a + rhs.a, self.b + rhs.b)

    __radd__ = __add__

    def __neg__(self) -> "Qsqrt17":
        return Qsqrt17(-self.a, -self.b)

    def __sub__(self, other: int | Fraction | "Qsqrt17") -> "Qsqrt17":
        return self + (-self.coerce(other))

    def __rsub__(self, other: int | Fraction | "Qsqrt17") -> "Qsqrt17":
        return self.coerce(other) - self

    def __mul__(self, other: int | Fraction | "Qsqrt17") -> "Qsqrt17":
        rhs = self.coerce(other)
        return Qsqrt17(
            self.a * rhs.a + 17 * self.b * rhs.b,
            self.a * rhs.b + self.b * rhs.a,
        )

    __rmul__ = __mul__

    def inverse(self) -> "Qsqrt17":
        norm = self.a * self.a - 17 * self.b * self.b
        if norm == 0:
            raise ZeroDivisionError
        return Qsqrt17(self.a / norm, -self.b / norm)

    def __truediv__(self, other: int | Fraction | "Qsqrt17") -> "Qsqrt17":
        return self * self.coerce(other).inverse()

    def __rtruediv__(self, other: int | Fraction | "Qsqrt17") -> "Qsqrt17":
        return self.coerce(other) / self

    def __pow__(self, exponent: int) -> "Qsqrt17":
        if exponent < 0:
            return (self.inverse()) ** (-exponent)
        result = Qsqrt17(Fraction(1))
        base = self
        power = exponent
        while power:
            if power & 1:
                result = result * base
            base = base * base
            power //= 2
        return result

    def encode(self) -> dict[str, str]:
        return {"rational": str(self.a), "sqrt17_coefficient": str(self.b)}


def upper_bound_report() -> dict:
    groups = {name: frozenset({index}) for index, name in enumerate("ABCD")}
    x_rational = Fraction(16, 9)
    samples = []
    for scale in (4, 10, 100, 1_000):
        ratio = ratio_at(
            five_vector_polynomial(x_rational, scale),
            groups,
            Fraction(1, scale**8),
        )
        samples.append(
            {
                "scale": scale,
                "epsilon": str(Fraction(1, scale**8)),
                "ratio": str(ratio),
            }
        )

    if Fraction(samples[0]["ratio"]) != EXPECTED_FINITE_WITNESS:
        raise AssertionError("finite rational witness changed")

    sqrt17 = Qsqrt17(Fraction(0), Fraction(1))
    x_star = (3 + sqrt17) / 4
    derivative_polynomial = 2 * x_star**2 - 3 * x_star - 1
    optimized = (x_star + 1) ** 4 / (2 * x_star * (2 * x_star + 1) ** 2)
    expected_optimized = Qsqrt17(Fraction(-107, 128), Fraction(51, 128))
    if derivative_polynomial != Qsqrt17(Fraction(0)):
        raise AssertionError("critical point identity failed")
    if optimized != expected_optimized:
        raise AssertionError("optimized radical identity failed")

    return {
        "finite_rational_witness": {
            "L": 4,
            "x": "16/9",
            "epsilon": "1/65536",
            "ratio": str(EXPECTED_FINITE_WITNESS),
            "below_nine_tenths": EXPECTED_FINITE_WITNESS < Fraction(9, 10),
        },
        "rational_parameter_limit": {
            "x": "16/9",
            "ratio": str(five_vector_limit(x_rational)),
        },
        "optimized_limit": {
            "x_star": {"rational": "3/4", "sqrt17_coefficient": "1/4"},
            "ratio": optimized.encode(),
            "display": "(-107+51*sqrt(17))/128",
            "decimal_display": "0.8068623977070366",
        },
        "convergence_samples": samples,
    }


def principal_minor_separation_report() -> dict:
    """Exact two-variable obstruction to a common 4x4 principal-minor model."""
    epsilon = Fraction(1, 3)
    groups = {
        "A": frozenset({0}),
        "B": frozenset({1}),
        "C": frozenset(),
        "D": frozenset({0, 1}),
    }

    def value(subset: frozenset[int]) -> Fraction:
        x0 = epsilon + (0 in subset)
        x1 = epsilon + (1 in subset)
        return x0 * x1

    unions = {
        "A": groups["A"],
        "B": groups["B"],
        "AB": groups["A"] | groups["B"],
        "AC": groups["A"] | groups["C"],
        "AD": groups["A"] | groups["D"],
        "BC": groups["B"] | groups["C"],
        "BD": groups["B"] | groups["D"],
        "ABC": groups["A"] | groups["B"] | groups["C"],
        "ABD": groups["A"] | groups["B"] | groups["D"],
        "CD": groups["C"] | groups["D"],
    }
    values = {name: value(subset) for name, subset in unions.items()}
    order = ["AB", "AC", "AD", "BC", "BD", "A", "B", "ABC", "ABD", "CD"]

    a, b, u = values["A"], values["B"], values["AB"]
    p_value, q_value = values["ABC"], values["ABD"]
    x = u * values["AC"] / (a * p_value)
    y = u * values["BC"] / (b * p_value)
    z = u * values["AD"] / (a * q_value)
    w = u * values["BD"] / (b * q_value)
    h = u**3 * values["CD"] / (a * b * p_value * q_value)
    ratio = x * y * z * w / h

    expected_values = [
        Fraction(16, 9),
        Fraction(4, 9),
        Fraction(16, 9),
        Fraction(4, 9),
        Fraction(16, 9),
        Fraction(4, 9),
        Fraction(4, 9),
        Fraction(16, 9),
        Fraction(16, 9),
        Fraction(16, 9),
    ]
    if [values[name] for name in order] != expected_values:
        raise AssertionError("principal-minor separation values changed")
    if (x, y, z, w, h) != (1, 1, 4, 4, 16) or ratio != 1:
        raise AssertionError("principal-minor canonical invariants changed")
    if h - 13 != 3:
        raise AssertionError("principal-minor separator gap changed")

    certificate = {
        "f": "x0*x1",
        "epsilon": "1/3",
        "A": [0],
        "B": [1],
        "C": [],
        "D": [0, 1],
        "V_order": order,
        "V": [str(values[name]) for name in order],
        "canonical": [str(entry) for entry in (x, y, z, w, h)],
        "ratio": str(ratio),
        "separator_bound": "13",
        "target_h": "16",
        "gap": "3",
    }
    stream = json.dumps(
        certificate, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    )
    digest = hashlib.sha256(stream.encode("ascii")).hexdigest()
    if digest != EXPECTED_SEPARATION_HASH:
        raise AssertionError(f"principal-minor certificate changed: {digest}")
    return {
        **certificate,
        "certificate_stream_sha256": digest,
        "scope": (
            "Rules out a simultaneous common-scalar plus diagonal modular-gauge "
            "4x4 PSD principal-minor representation; it does not rule out an "
            "indirect use of principal-minor inequalities."
        ),
    }


VALUATED_PAIRS = tuple(itertools.combinations(range(4), 2))
VALUATED_LEFT = ({0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3})
VALUATED_RIGHT = ({0}, {1}, {0, 1, 2}, {0, 1, 3}, {2, 3})


def satisfies_basis_exchange(finite_pairs: tuple[frozenset[int], ...]) -> bool:
    bases = set(finite_pairs)
    if not bases:
        return False
    for first in bases:
        for second in bases:
            for removed in first - second:
                if not any(
                    frozenset((first - {removed}) | {added}) in bases
                    for added in second - first
                ):
                    return False
    return True


def satisfies_four_point(values: tuple[int | float, ...]) -> bool:
    sums = (
        values[0] + values[5],
        values[1] + values[4],
        values[2] + values[3],
    )
    minimum = min(sums)
    return sum(value == minimum for value in sums) >= 2


def regularized_rank(
    values: tuple[int | float, ...], subset: set[int], lambda_value: int
) -> int:
    return max(
        lambda_value * len(set(pair) & subset) - int(value)
        for pair, value in zip(VALUATED_PAIRS, values)
        if math.isfinite(value)
    )


def finite_ingleton_gap(
    values: tuple[int | float, ...], lambda_value: int
) -> int:
    left = sum(
        regularized_rank(values, subset, lambda_value) for subset in VALUATED_LEFT
    )
    right = sum(
        regularized_rank(values, subset, lambda_value) for subset in VALUATED_RIGHT
    )
    return left - right


def valuated_rank_two_report() -> dict:
    digest = hashlib.sha256(b"rank2-valuated-ingleton-n4-v1\n")
    input_count = 0
    case_count = 0
    negative_count = 0
    minimum_gap: int | None = None
    maximum_gap: int | None = None
    for values in itertools.product((0, 1, 2, 3, math.inf), repeat=6):
        finite_values = [value for value in values if math.isfinite(value)]
        if not finite_values or min(finite_values) != 0:
            continue
        finite_pairs = tuple(
            frozenset(pair)
            for pair, value in zip(VALUATED_PAIRS, values)
            if math.isfinite(value)
        )
        if not satisfies_basis_exchange(finite_pairs):
            continue
        if not satisfies_four_point(values):
            continue
        input_count += 1
        encoded = [
            "inf" if not math.isfinite(value) else int(value) for value in values
        ]
        for lambda_value in (1, 2, 3):
            gap = finite_ingleton_gap(values, lambda_value)
            digest.update(
                (json.dumps([encoded, lambda_value, gap], separators=(",", ":")) + "\n").encode(
                    "ascii"
                )
            )
            case_count += 1
            negative_count += gap < 0
            minimum_gap = gap if minimum_gap is None else min(minimum_gap, gap)
            maximum_gap = gap if maximum_gap is None else max(maximum_gap, gap)

    report = {
        "pair_order": [list(pair) for pair in VALUATED_PAIRS],
        "valuation_alphabet": [0, 1, 2, 3, "inf"],
        "lambda_values": [1, 2, 3],
        "valuated_rank_two_inputs": input_count,
        "input_lambda_cases": case_count,
        "negative_ingleton_gaps": negative_count,
        "minimum_gap": minimum_gap,
        "maximum_gap": maximum_gap,
        "certificate_stream_sha256": digest.hexdigest(),
        "scope": (
            "Finite n=4 convention regression only. The general theorem uses "
            "rank-two realization and the module-length proof."
        ),
    }
    expected = (1_929, 5_787, 0, 0, 9, EXPECTED_VALUATED_HASH)
    actual = (
        input_count,
        case_count,
        negative_count,
        minimum_gap,
        maximum_gap,
        digest.hexdigest(),
    )
    if actual != expected:
        raise AssertionError(f"valuated certificate changed: {actual!r}")
    return report


def build_report() -> dict:
    source_path = Path(__file__).resolve()
    paper_path = ROOT / "paper" / "main.tex"
    bibliography_path = ROOT / "paper" / "references.bib"
    return {
        "schema": "lorentzian-ingleton-quadratics-certificates-v2",
        "provenance": {
            "reproduce_py_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
            "paper_main_tex_sha256": hashlib.sha256(paper_path.read_bytes()).hexdigest(),
            "paper_references_bib_sha256": hashlib.sha256(
                bibliography_path.read_bytes()
            ).hexdigest(),
            "runtime": "Python >=3.11; standard library only",
        },
        "theorem_scope": (
            "The computations certify formulas and finite conventions; they do "
            "not replace the proof of the uniform positive lower bound."
        ),
        "upper_bound": upper_bound_report(),
        "principal_minor_model_separation": principal_minor_separation_report(),
        "finite_valuated_rank_two_audit": valuated_rank_two_report(),
    }


def render_report(report: dict) -> str:
    return json.dumps(report, ensure_ascii=True, indent=2, sort_keys=True) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--output", type=Path, help="write a freshly generated report")
    mode.add_argument("--check", type=Path, help="compare against a committed report")
    args = parser.parse_args()

    report = build_report()
    rendered = render_report(report)
    rendered_bytes = rendered.encode("ascii")
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_bytes(rendered_bytes)
    elif args.check:
        expected = args.check.read_bytes()
        if rendered_bytes != expected:
            raise SystemExit(f"certificate mismatch: {args.check}")
        print(f"CERTIFICATES_OK {args.check}")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
