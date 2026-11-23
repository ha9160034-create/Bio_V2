"""Unit tests for PDF report generation in reporting/pdf_generator.py."""

import base64
import io
import unittest
from PIL import Image

from reporting.pdf_generator import (
    build_analysis_pdf,
    build_docking_pdf,
    build_md_pdf,
    sanitize_text,
)


class MockFigure:
    """Mock matplotlib Figure object providing a savefig method."""

    def __init__(self, color: str = "blue", width: int = 300, height: int = 150) -> None:
        self.color = color
        self.width = width
        self.height = height

    def savefig(self, buf, format: str = "png", dpi: int = 150, **kwargs) -> None:
        img = Image.new("RGB", (self.width, self.height), color=self.color)
        img.save(buf, format=format.upper())



def create_sample_png_base64(width: int = 100, height: int = 80, color: str = "blue") -> str:
    """Generate a valid base64-encoded PNG data URL."""
    img = Image.new("RGB", (width, height), color=color)
    buf = io.BytesIO()
    img.save(buf, format="PNG")
    raw_b64 = base64.b64encode(buf.getvalue()).decode("ascii")
    return f"data:image/png;base64,{raw_b64}"


class TestPDFReportGenerator(unittest.TestCase):
    """Test suite verifying PDF generation across analysis, docking, and MD tabs."""

    @classmethod
    def setUpClass(cls):
        cls.sample_b64_url = create_sample_png_base64(120, 90, color="teal")
        # Also raw base64 without prefix
        raw_only = cls.sample_b64_url.split(",", 1)[1]
        cls.sample_raw_b64 = raw_only

    def test_sanitize_text(self):
        """Verify that common biophysical symbols and emojis are safely sanitized."""
        raw = "RMSD ≤ 2.0 Å, ΔG: -8.5 kcal/mol, ±0.3, ⚯ ligand, Cα atoms, • item, 🟢 Pass"
        clean = sanitize_text(raw)
        self.assertIn("Angstrom", clean)
        self.assertIn("dG", clean)
        self.assertIn("+/-", clean)
        self.assertIn("<=", clean)
        self.assertIn("[PASS]", clean)
        # Verify clean can be encoded to latin-1 without exception
        clean.encode("latin-1")

    # ── Tab 1: Analysis PDF Tests ─────────────────────────────────────────────

    def test_build_analysis_pdf_empty_screenshots(self):
        """Verify build_analysis_pdf returns valid PDF bytes with empty screenshots."""
        data = {
            "wt_id": "1M17",
            "wt_chain": "A",
            "mt_id": "1M17_mutant",
            "mt_chain": "A",
            "coordinate_mode": "Deposited",
            "alignment": {
                "score": 1250.5,
                "identity_pct": 98.4,
                "mutations": 2,
                "gaps": 0,
                "ca_rmsd": 0.412,
                "text": "CLUSTAL O(1.2.4) multiple sequence alignment\nWT  ACDEFGHIKLMNPQRSTVWY\nMT  ACDEFGHIKLMNPQRSTVWY",
            },
            "comparison_rows": [
                {
                    "res_num_h": 45,
                    "res_num_m": 45,
                    "WT_Res": "VAL",
                    "MT_Res": "ALA",
                    "SASA_H": 32.5,
                    "SASA_M": 48.1,
                    "SASA_Delta": 15.6,
                    "Status": "[MUTATION]",
                }
            ],
            "notes": "Mutation is located away from catalytic pocket.",
        }
        pdf_bytes = build_analysis_pdf(data, screenshots=[])
        self.assertIsInstance(pdf_bytes, bytes)
        self.assertTrue(pdf_bytes.startswith(b"%PDF-"))
        self.assertGreater(len(pdf_bytes), 500)

    def test_build_analysis_pdf_with_screenshot(self):
        """Verify build_analysis_pdf returns valid PDF bytes with sample base64 screenshot."""
        data = {
            "wt_id": "1M17",
            "wt_chain": "A",
            "mt_id": "1M17_mut",
            "mt_chain": "A",
            "alignment": {
                "identity_pct": 100.0,
                "ca_rmsd": 0.0,
            },
        }
        screenshots = [self.sample_b64_url, self.sample_raw_b64]
        pdf_bytes = build_analysis_pdf(data, screenshots=screenshots)
        self.assertIsInstance(pdf_bytes, bytes)
        self.assertTrue(pdf_bytes.startswith(b"%PDF-"))
        self.assertGreater(len(pdf_bytes), 1000)

    # ── Tab 2: Docking PDF Tests ──────────────────────────────────────────────

    def test_build_docking_pdf_empty_screenshots(self):
        """Verify build_docking_pdf returns valid PDF bytes with empty screenshots."""
        data = {
            "target_id": "1M17",
            "ligand_name": "Gefitinib",
            "workflow_mode": "Single Target",
            "validation_status": "Protocol validated on WT reference",
            "validation_details": "Redocking RMSD is 0.72 Angstrom <= 2.0 Angstrom threshold.",
            "lowest_vina_score": -9.14,
            "poses": [
                {"rank": 1, "score": -9.14, "rmsd_lb": 0.0, "rmsd_ub": 0.0},
                {"rank": 2, "score": -8.82, "rmsd_lb": 1.2, "rmsd_ub": 1.8},
            ],
            "clusters": [{"id": 1, "size": 8}, {"id": 2, "size": 2}],
            "properties": {
                "mw": 446.9,
                "logp": 3.75,
                "hbd": 1,
                "hba": 7,
                "rotatable_bonds": 6,
                "lipinski_violations": 0,
            },
        }
        pdf_bytes = build_docking_pdf(data, screenshots=[])
        self.assertIsInstance(pdf_bytes, bytes)
        self.assertTrue(pdf_bytes.startswith(b"%PDF-"))
        self.assertGreater(len(pdf_bytes), 500)

    def test_build_docking_pdf_with_screenshot(self):
        """Verify build_docking_pdf returns valid PDF bytes with sample base64 screenshot."""
        data = {
            "target_id": "1M17",
            "ligand_name": "Erlotinib",
            "lowest_vina_score": -8.9,
            "poses": [{"rank": 1, "score": -8.9}],
        }
        screenshots = [self.sample_b64_url]
        pdf_bytes = build_docking_pdf(data, screenshots=screenshots)
        self.assertIsInstance(pdf_bytes, bytes)
        self.assertTrue(pdf_bytes.startswith(b"%PDF-"))
        self.assertGreater(len(pdf_bytes), 1000)

    # ── Tab 3: MD PDF Tests ───────────────────────────────────────────────────

    def test_build_md_pdf_empty_screenshots(self):
        """Verify build_md_pdf returns valid PDF bytes with empty screenshots and plots."""
        # Create sample mock figures
        fig1 = MockFigure(color="blue", width=400, height=200)
        fig2 = MockFigure(color="green", width=400, height=200)

        plots = {"rmsd_trajectory": fig1, "hbond_count": fig2}

        data = {
            "target_pdb": "1M17",
            "ligand_name": "Gefitinib",
            "gromacs_version": "GROMACS 2026.3",
            "hardware": "CUDA Acceleration",
            "unit_cell": "Rhombic Dodecahedron",
            "duration_ns": 10.0,
            "ensemble": "NPT",
            "temperature_k": 300.0,
            "mean_ligand_rmsd": 1.45,
            "mean_hbond_count": 2.8,
            "mmpbsa": {
                "dg_bind": -24.8,
                "van_der_waals": -38.2,
                "electrostatic": -12.4,
                "polar_solvation": 29.5,
                "nonpolar_solvation": -3.7,
            },
        }

        pdf_bytes = build_md_pdf(data, plots=plots, screenshots=[])
        self.assertIsInstance(pdf_bytes, bytes)
        self.assertTrue(pdf_bytes.startswith(b"%PDF-"))
        self.assertGreater(len(pdf_bytes), 1000)

    def test_build_md_pdf_with_screenshot(self):
        """Verify build_md_pdf returns valid PDF bytes with sample base64 screenshot."""
        fig = MockFigure(color="purple", width=300, height=150)
        plots = {"stability": fig}

        data = {
            "target_pdb": "1M17",
            "ligand_name": "Gefitinib",
            "mean_ligand_rmsd": 1.2,
        }
        screenshots = [self.sample_b64_url]

        pdf_bytes = build_md_pdf(data, plots=plots, screenshots=screenshots)
        self.assertIsInstance(pdf_bytes, bytes)
        self.assertTrue(pdf_bytes.startswith(b"%PDF-"))
        self.assertGreater(len(pdf_bytes), 1500)

    # ── Edge Cases & Robustness ───────────────────────────────────────────────

    def test_handle_empty_data_gracefully(self):
        """Verify all three builders handle empty data dicts and None values cleanly."""
        pdf1 = build_analysis_pdf({}, [])
        pdf2 = build_docking_pdf({}, [])
        pdf3 = build_md_pdf({}, {}, [])

        for b in [pdf1, pdf2, pdf3]:
            self.assertIsInstance(b, bytes)
            self.assertTrue(b.startswith(b"%PDF-"))

    def test_special_characters_in_all_builders(self):
        """Verify that strings containing Å, ±, ΔG, ⚯, emojis do not trigger encoding errors."""
        tricky_data = {
            "title": "Simulation: ⚯ Complex Å ± ΔG ≥ 99%",
            "notes": "Backbone Cα RMSD: 1.25 Å; ΔSASA: -12.4 Å²; Status: 🟢 Pass / 🔴 Mutation",
            "validation_details": "Target ⚯ Ligand: ± 0.05 Å error threshold.",
            "mmpbsa": {"dg_bind": -35.2},
        }
        pdf1 = build_analysis_pdf(tricky_data, [])
        pdf2 = build_docking_pdf(tricky_data, [])
        pdf3 = build_md_pdf(tricky_data, {}, [])

        for b in [pdf1, pdf2, pdf3]:
            self.assertIsInstance(b, bytes)
            self.assertTrue(b.startswith(b"%PDF-"))


if __name__ == "__main__":
    unittest.main()
