"""Regression tests for the single-docking validation status mapping.

The UI previously checked Top-N recovery before ``protocol_validated``, so a
correct top-1 pose was always labelled "Top-10 recovery only" and the UI warned
that top-1 failed even when it had reproduced the crystal pose.
"""

import unittest

from docking.validation import derive_single_docking_status


def _details(**overrides):
    base = {
        "top1_rmsd_pass": True,
        "top1_cluster_pass": True,
        "recovery_rmsd_pass": True,
        "recovery_cluster_pass": True,
        "protocol_validated": True,
        "recovery": {"rmsd_angstrom": 1.4},
    }
    base.update(overrides)
    return base


class TestDeriveSingleDockingStatus(unittest.TestCase):
    def test_protocol_validated_wins_over_recovery(self):
        status, warnings = derive_single_docking_status(_details(), seed_count=3)
        self.assertEqual(status, "validation passed")
        self.assertEqual(warnings, [])

    def test_top10_only_when_top1_rmsd_fails(self):
        details = _details(protocol_validated=False, top1_rmsd_pass=False)
        status, warnings = derive_single_docking_status(details, seed_count=3)
        self.assertEqual(status, "validation passed (Top-10 recovery)")
        self.assertTrue(any("did not rank it first" in warning for warning in warnings))

    def test_top10_only_when_top1_cluster_unsupported(self):
        details = _details(protocol_validated=False, top1_cluster_pass=False)
        status, warnings = derive_single_docking_status(details, seed_count=3)
        self.assertEqual(status, "validation passed (Top-10 recovery)")
        self.assertTrue(any("lacks support" in warning for warning in warnings))

    def test_failed_when_recovery_rmsd_fails(self):
        details = _details(
            protocol_validated=False,
            top1_rmsd_pass=False,
            recovery_rmsd_pass=False,
        )
        status, warnings = derive_single_docking_status(details, seed_count=3)
        self.assertEqual(status, "validation failed")
        self.assertTrue(warnings)

    def test_failed_when_recovery_cluster_unsupported(self):
        details = _details(
            protocol_validated=False,
            top1_rmsd_pass=False,
            recovery_cluster_pass=False,
        )
        status, warnings = derive_single_docking_status(details, seed_count=3)
        self.assertEqual(status, "validation failed")
        self.assertTrue(any("lacks support" in warning for warning in warnings))

    def test_insufficient_seeds_is_exploratory(self):
        status, warnings = derive_single_docking_status(_details(), seed_count=2)
        self.assertEqual(status, "exploratory—insufficient independent seeds")
        self.assertEqual(warnings, [])

    def test_missing_details_is_unavailable(self):
        status, warnings = derive_single_docking_status(None, seed_count=3)
        self.assertEqual(status, "validation unavailable")
        self.assertEqual(warnings, [])


if __name__ == "__main__":
    unittest.main()
