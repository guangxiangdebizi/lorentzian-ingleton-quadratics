from __future__ import annotations

from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import explicit_floor  # noqa: E402
import hostile_audit  # noqa: E402


class ExplicitFloorTests(unittest.TestCase):
    def test_universal_constant_and_provenance_dependencies(self):
        report = json.loads((ROOT / "reports" / "certificates.json").read_text(encoding="utf-8"))
        self.assertEqual(report["schema"], "lorentzian-ingleton-quadratics-certificates-v3")
        self.assertEqual(496 * 31 ** 2, 476656)
        self.assertEqual(report["explicit_quadratic_floor"]["universal_constant"], "476656^(-5)")
        self.assertEqual(report["explicit_quadratic_floor"]["largest_active_clone_count"], 32)
        for key, path in (
            ("explicit_floor_py_sha256", "src/explicit_floor.py"),
            ("hostile_audit_py_sha256", "src/hostile_audit.py"),
            ("test_explicit_floor_py_sha256", "tests/test_explicit_floor.py"),
        ):
            self.assertEqual(report["provenance"][key], hashlib.sha256((ROOT / path).read_bytes()).hexdigest())

    def test_signature_negative_controls(self):
        disconnected = [[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]]
        self.assertIsNone(hostile_audit.certify_and_template(disconnected))
        self.assertIsNone(hostile_audit.certify_and_template([[0, 0], [0, 0]]))
        self.assertFalse(hostile_audit.psd_exact([[0, 1], [1, 0]]))
        self.assertTrue(hostile_audit.psd_exact([[1, 1], [1, 1]]))

    def test_adversarial_certificate_counts(self):
        report = json.loads((ROOT / "reports" / "certificates.json").read_text(encoding="utf-8"))
        independent = report["independent_explicit_floor_audit"]["integer_inputs"]
        self.assertEqual(independent["enumerated_matrices"], 4096)
        self.assertEqual(independent["accepted"], 1317)
        self.assertEqual(independent["rejected_including_zero"], 2779)
        self.assertEqual(independent["accepted_with_loops"], 234)
        self.assertEqual(independent["accepted_with_parallel_pairs"], 531)
        self.assertEqual(independent["query_checks"], 3951)
        self.assertEqual(report["explicit_quadratic_floor"]["random_cases"], 160)

    def test_exact_parallel_support_and_loop_deletion(self):
        matrix = [[0, 0, 7, 0], [0, 0, 7, 0], [7, 7, 0, 0], [0, 0, 0, 0]]
        b, c, active = hostile_audit.certify_and_template(matrix)
        self.assertEqual(active, [0, 1, 2])
        self.assertEqual(b[0][1], 0)
        self.assertEqual(c[0][1], 0)
        self.assertGreater(c[0][2], 0)

    def test_single_square_exact_polarization(self):
        polarized = explicit_floor.aggregate_and_polarize([[F(0), F(3)], [F(3), F(0)]], [0, 0], 1)
        self.assertEqual(explicit_floor.evaluate(polarized, [F(2), F(2)]), 12)
        for epsilon in (F(1, 2 ** 40), F(1), F(2 ** 40)):
            ratio, max_ratio, n = explicit_floor.check(polarized, [15, 15], epsilon / (1 + epsilon))
            self.assertEqual((ratio, max_ratio, n), (1, 1, 2))

    def test_two_clone_boundary_exhaustion(self):
        report = hostile_audit.audit_two_clones()
        self.assertEqual(report["all_memberships_and_regularizers"], 768)
        self.assertEqual(report["constant"], 1)

    def test_laurent_quotients_measured_as_dimensions(self):
        report = hostile_audit.audit_laurent_modules()
        self.assertEqual(report["direct_quotient_dimension_checks"], 384)
        self.assertEqual(report["ambient_dimension_checks"], 24)
        self.assertEqual(report["parallel_class_cases"], 12)
        self.assertEqual(report["negative_determinant_valuations"], 102)


if __name__ == "__main__":
    unittest.main()
