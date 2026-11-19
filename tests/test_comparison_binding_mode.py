"""Regression tests for reference-guided consensus selection and the non-fatal
binding-mode outcome used by the matched comparison."""

import unittest
from pathlib import Path
from unittest.mock import patch

from docking.comparison import verify_comparison_clusters
from docking.models import DockingPose, PoseCluster


def _cluster(cluster_id: int, seeds: list[int], median: float) -> PoseCluster:
    poses = [
        DockingPose(
            seed=seed,
            rank=1,
            score=median,
            rmsd_lb=0.0,
            rmsd_ub=0.0,
            pdbqt_block=f"MODEL 1\nREMARK VINA RESULT: -8.0 0.0 0.0\nATOM      1  C   LIG     1       0.000   0.000   0.000\n",
            cluster_id=cluster_id,
        )
        for seed in seeds
    ]
    return PoseCluster(
        cluster_id=cluster_id,
        seed_count=len(seeds),
        supporting_seeds=seeds,
        top_score=median,
        median_score=median,
        representative_pose=poses[0],
        representative_pose_source={"seed": poses[0].seed, "rank": 1},
    )


class TestBindingModeOutcome(unittest.TestCase):
    def test_incompatible_binding_mode_is_recorded_not_raised(self):
        healthy = _cluster(1, [42, 101, 2024], -7.0)
        mutant = _cluster(1, [42, 101, 2024], -6.9)
        details = {}
        with patch(
            "docking.comparison.evaluate_binding_mode_compatibility",
            return_value=(
                False,
                {"evaluated": True, "ligand_pose_rmsd_angstrom": 3.9, "reason": "different pose"},
            ),
        ):
            h, m, delta, warnings = verify_comparison_clusters(
                [healthy],
                [mutant],
                healthy_receptor_pdb=Path("healthy.pdb"),
                mutant_receptor_pdb=Path("mutant.pdb"),
                details=details,
                strict_binding_mode=False,
            )
        self.assertIs(delta, None)
        self.assertEqual(details["binding_mode_status"], "not_evaluated")
        self.assertTrue(any("ΔScore was not evaluated" in item for item in warnings))
        self.assertEqual(h.cluster_id, 1)
        self.assertEqual(m.cluster_id, 1)

    def test_incompatible_binding_mode_still_raises_in_strict_mode(self):
        healthy = _cluster(1, [42, 101], -7.0)
        mutant = _cluster(1, [42, 101], -6.9)
        with patch(
            "docking.comparison.evaluate_binding_mode_compatibility",
            return_value=(False, {"evaluated": True, "reason": "different pose"}),
        ):
            with self.assertRaises(ValueError):
                verify_comparison_clusters(
                    [healthy],
                    [mutant],
                    healthy_receptor_pdb=Path("healthy.pdb"),
                    mutant_receptor_pdb=Path("mutant.pdb"),
                )

    def test_reference_guided_healthy_cluster_overrides_dominant(self):
        dominant = _cluster(1, [42, 101, 2024], -7.2)
        reference = _cluster(2, [42, 101], -6.8)
        mutant = _cluster(1, [42, 101, 2024], -6.9)
        with patch(
            "docking.comparison.evaluate_binding_mode_compatibility",
            return_value=(True, {"evaluated": True}),
        ):
            h, _, _, _ = verify_comparison_clusters(
                [dominant, reference],
                [mutant],
                healthy_receptor_pdb=Path("healthy.pdb"),
                mutant_receptor_pdb=Path("mutant.pdb"),
                healthy_cluster=reference,
            )
        self.assertEqual(h.cluster_id, 2)


if __name__ == "__main__":
    unittest.main()
