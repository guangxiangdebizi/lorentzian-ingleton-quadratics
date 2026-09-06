from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import clifford_floor
import gaussian_floor


class GaussianCliffordTests(unittest.TestCase):
    def test_current_constants_and_provenance(self):
        report = json.loads((ROOT / "reports/certificates.json").read_text())
        bound = report["quadratic_lower_bound"]
        self.assertEqual(bound["all_overlaps"], "1/16384")
        self.assertEqual(bound["disjoint_queries"], "1/16")
        self.assertFalse(bound["sharpness_claimed"])
        for key, path in (
            ("gaussian_floor_py_sha256", "src/gaussian_floor.py"),
            ("clifford_floor_py_sha256", "src/clifford_floor.py"),
            ("test_gaussian_clifford_py_sha256", "tests/test_gaussian_clifford.py"),
        ):
            self.assertEqual(report["provenance"][key], hashlib.sha256((ROOT / path).read_bytes()).hexdigest())

    def test_exact_nonunit_gaussian_fixture(self):
        result = gaussian_floor.audit(gaussian_floor.nonunit_fixture(), (1, 1, 1, 1))
        self.assertLess(F(result["ratio"]), 1)
        self.assertGreaterEqual(F(result["ratio"]), F(result["auxiliary_lower"]))
        self.assertGreaterEqual(F(result["auxiliary_lower"]), F(1, 4))

    def test_exterior_generators_and_negative_controls(self):
        for rank in range(1, 5):
            clifford_floor.check_clifford(clifford_floor.clifford(rank))
        controls = clifford_floor.negative_controls()
        self.assertTrue(controls["unsigned_generator_rejected"])
        self.assertTrue(controls["two_positive_hessian_rejected"])
        self.assertTrue(controls["negative_diagonal_atom_obstruction_detected"])

    def test_complete_symbolic_identity_not_only_samples(self):
        clifford_floor.symbolic_determinant(2)

    def test_rank_zero_and_full_certificate_counts(self):
        report = json.loads((ROOT / "reports/certificates.json").read_text())
        audit = report["clifford_power_audit"]
        self.assertEqual(audit["quadratic_examples"], 17)
        self.assertEqual(audit["regularized_cases"], 102)
        self.assertEqual(len(audit["symbolic_identities"]), 3)
        self.assertTrue(any(row["rank_K"] == 0 and row["matrix_size"] == 2 for row in audit["records"]))
        self.assertEqual(report["gaussian_auxiliary_audit"]["random_cases"], 40)

    def test_determinant_power_exponent_cancellation(self):
        for rank in range(1, 9):
            size = 2**rank
            self.assertEqual(F(1, 16384)**(size // 2), F(1, 128**size))
            self.assertEqual(F(1, 16)**(size // 2), F(1, 4**size))


if __name__ == "__main__":
    unittest.main()
