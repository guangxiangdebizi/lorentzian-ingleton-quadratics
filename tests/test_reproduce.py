from __future__ import annotations

import json
import sys
import unittest
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import reproduce  # noqa: E402


class ReproductionTests(unittest.TestCase):
    def test_finite_upper_bound_witness(self) -> None:
        report = reproduce.upper_bound_report()
        witness = report["finite_rational_witness"]
        self.assertEqual(Fraction(witness["ratio"]), reproduce.EXPECTED_FINITE_WITNESS)
        self.assertTrue(witness["below_nine_tenths"])

    def test_optimized_radical_identity(self) -> None:
        report = reproduce.upper_bound_report()["optimized_limit"]
        self.assertEqual(report["display"], "(-107+51*sqrt(17))/128")
        self.assertEqual(report["decimal_display"], "0.8068623977070366")

    def test_finite_valuated_rank_two_regression(self) -> None:
        report = reproduce.valuated_rank_two_report()
        self.assertEqual(report["valuated_rank_two_inputs"], 1_929)
        self.assertEqual(report["input_lambda_cases"], 5_787)
        self.assertEqual(report["negative_ingleton_gaps"], 0)
        self.assertEqual(
            report["certificate_stream_sha256"], reproduce.EXPECTED_VALUATED_HASH
        )

    def test_principal_minor_model_separation(self) -> None:
        report = reproduce.principal_minor_separation_report()
        self.assertEqual(report["canonical"], ["1", "1", "4", "4", "16"])
        self.assertEqual(report["ratio"], "1")
        self.assertEqual(report["gap"], "3")
        self.assertEqual(
            report["certificate_stream_sha256"], reproduce.EXPECTED_SEPARATION_HASH
        )

    def test_committed_report_is_current(self) -> None:
        committed = json.loads(
            (ROOT / "reports" / "certificates.json").read_text(encoding="utf-8")
        )
        self.assertEqual(committed, reproduce.build_report())


if __name__ == "__main__":
    unittest.main()
