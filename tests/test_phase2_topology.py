"""Phase 2 ligand-topology regression tests.

Scientific references: ``wiki/open-force-field-and-ligand-parameterization.md``,
``wiki/force-fields.md``, and ``wiki/gromacs-protein-ligand-protocol.md``.
"""

from __future__ import annotations

import importlib.util
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest.mock import patch

from md.topology import TopologyBuildError, _enforce_ligand_charge
from acpype_wrapper import (
    EXIT_ACPYPE_NO_AMBER,
    EXIT_CHARGE,
    EXIT_NO_OUTPUT,
    main as acpype_wrapper_main,
)

ROOT = Path(__file__).resolve().parents[1]


def _write_itp(path: Path, charges: list[float]) -> None:
    lines = ["[ moleculetype ]", "LIG 3", "", "[ atoms ]"]
    lines.extend(
        f"{index} c3 1 LIG C{index} 1 {charge:.6f} 12.011"
        for index, charge in enumerate(charges, start=1)
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_gro(path: Path, atom_count: int) -> None:
    rows = ["test ligand", str(atom_count)]
    rows.extend(
        f"    1LIG     C{index:<2d}{index:5d}   0.000   0.000   0.000"
        for index in range(1, atom_count + 1)
    )
    rows.append("   2.00000   2.00000   2.00000")
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")


def _load_converter_with_stubbed_parmed():
    """Import pure converter helpers without requiring AmberTools/ParmEd."""
    module_name = "_phase2_converter_under_test"
    fake_parmed = types.ModuleType("parmed")
    fake_parmed.load_file = lambda *args, **kwargs: None
    with patch.dict(sys.modules, {"parmed": fake_parmed}):
        spec = importlib.util.spec_from_file_location(
            module_name, ROOT / "convert_amber_to_gmx.py"
        )
        module = importlib.util.module_from_spec(spec)
        assert spec and spec.loader
        spec.loader.exec_module(module)
    return module


class TestACPYPEChargeBalancing(unittest.TestCase):
    def test_charge_rounding_residual_is_distributed_over_extreme_atoms(self):
        converter = _load_converter_with_stubbed_parmed()
        balanced, residual = converter.balance_charges([0.8, 0.8, -0.3], 1.0)
        self.assertAlmostEqual(residual, 0.3)
        self.assertAlmostEqual(sum(balanced), 1.0)
        self.assertAlmostEqual(balanced[0], 0.65)
        self.assertAlmostEqual(balanced[1], 0.65)
        self.assertAlmostEqual(balanced[2], -0.3)

    def test_equal_positive_and_negative_extremes_follow_acpype_tie_rule(self):
        converter = _load_converter_with_stubbed_parmed()
        balanced, _ = converter.balance_charges([-0.8, 0.8, 0.2], 0.0)
        self.assertAlmostEqual(sum(balanced), 0.0)
        self.assertAlmostEqual(balanced[0], -0.8)
        self.assertAlmostEqual(balanced[1], 0.6)
        self.assertAlmostEqual(balanced[2], 0.2)


class TestTopologyChargeGate(unittest.TestCase):
    def test_requested_neutral_ligand_charge_passes(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            itp = Path(temp_dir) / "ligand.itp"
            _write_itp(itp, [-0.2, 0.2])
            total, atom_count = _enforce_ligand_charge(itp, 0)
            self.assertAlmostEqual(total, 0.0)
            self.assertEqual(atom_count, 2)

    def test_unbalanced_ligand_charge_fails_before_grompp(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            itp = Path(temp_dir) / "ligand.itp"
            _write_itp(itp, [0.7, 0.3])
            with self.assertRaisesRegex(TopologyBuildError, r"sum\(q\)=\+1\.000000"):
                _enforce_ligand_charge(itp, 0)


class TestACPYPEWrapperFileStateMachine(unittest.TestCase):
    def _run(self, acpype_rc, amber_present, itp_charges, net_charge=0):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        old_cwd = Path.cwd()
        import os

        os.chdir(temp.name)
        self.addCleanup(os.chdir, old_cwd)
        Path("ligand.mdl").write_text("test input\n", encoding="utf-8")
        out = Path("LIG.acpype")
        out.mkdir()
        if amber_present:
            (out / "LIG_AC.prmtop").write_text(
                "%FLAG POINTERS\n%FLAG ATOM_NAME\n%FLAG CHARGE\n", encoding="utf-8"
            )
            (out / "LIG_AC.inpcrd").write_text("title\n1\n", encoding="utf-8")
        if itp_charges is not None:
            _write_itp(out / "LIG_GMX.itp", itp_charges)
            _write_gro(out / "LIG_GMX.gro", len(itp_charges))
        def rebuild(root, requested_charge):
            _write_itp(out / f"{root}_GMX.itp", [float(requested_charge or 0), 0.0])
            _write_gro(out / f"{root}_GMX.gro", 2)
            return 0

        with patch("acpype_wrapper._run_acpype", return_value=acpype_rc), patch(
            "acpype_wrapper._run_converter", side_effect=rebuild
        ) as converter:
            rc = acpype_wrapper_main(
                ["-i", "ligand.mdl", "-b", "LIG", "-n", str(net_charge), "-o", "gmx"]
            )
        return rc, converter

    def test_valid_acpype_files_are_accepted_without_overwrite(self):
        rc, converter = self._run(0, True, [-0.2, 0.2])
        self.assertEqual(rc, 0)
        converter.assert_not_called()

    def test_valid_files_survive_acpype_nonzero_exit(self):
        rc, converter = self._run(134, True, [-0.2, 0.2])
        self.assertEqual(rc, 0)
        converter.assert_not_called()

    def test_amber_files_trigger_fallback_for_bad_gmx_charge(self):
        rc, converter = self._run(134, True, [0.7, 0.302])
        self.assertEqual(rc, 0)
        converter.assert_called_once()

    def test_failed_acpype_without_amber_files_fails(self):
        rc, converter = self._run(1, False, None)
        self.assertEqual(rc, EXIT_ACPYPE_NO_AMBER)
        converter.assert_not_called()

    def test_success_without_any_output_fails(self):
        rc, converter = self._run(0, False, None)
        self.assertEqual(rc, EXIT_NO_OUTPUT)
        converter.assert_not_called()

    def test_bad_charge_without_amber_fallback_fails(self):
        rc, converter = self._run(0, False, [0.7, 0.3])
        self.assertEqual(rc, EXIT_CHARGE)
        converter.assert_not_called()


if __name__ == "__main__":
    unittest.main()
