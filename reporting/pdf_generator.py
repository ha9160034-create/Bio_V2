"""PDF Generator for BioBio Molecular Simulation Platform.

Uses fpdf2 to build professional PDF reports for:
1. Structural Analysis & Variant Comparison (Tab 1)
2. Molecular Docking Simulation (Tab 2)
3. Molecular Dynamics Simulation (Tab 3)

Handles missing/optional fields cleanly without raw 'None' values,
sanitizes Unicode for standard PDF latin-1 font compatibility,
and embeds base64 PNG screenshots and high-DPI matplotlib figures.
"""

from __future__ import annotations

import base64
import io
import unicodedata
from datetime import datetime
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

from fpdf import FPDF
from PIL import Image

# ── Unicode Sanitization Mapping ─────────────────────────────────────────────

UNICODE_REPLACEMENTS: dict[str, str] = {
    # Biophysical units & symbols
    "Å": "Angstrom",
    "Å": "Angstrom",
    "±": "+/-",
    "ΔG": "dG",
    "Δ": "d",
    "≥": ">=",
    "≤": "<=",
    "⚯": "x",
    "•": "-",
    "·": "-",
    "×": "x",
    "÷": "/",
    "μ": "u",
    "µ": "u",
    "α": "alpha",
    "β": "beta",
    "γ": "gamma",
    "δ": "delta",
    "ε": "epsilon",
    "Cα": "C-alpha",
    "Cβ": "C-beta",
    "²": "^2",
    "³": "^3",
    "°": " deg",
    "℃": " deg C",
    "℉": " deg F",
    "—": "--",
    "–": "-",
    "“": '"',
    "”": '"',
    "‘": "'",
    "’": "'",
    "…": "...",
    "≈": "~=",
    "≠": "!=",
    "∞": "inf",
    # Status badges and symbols
    "🟢": "[PASS]",
    "🔴": "[MUTATION]",
    "🟡": "[WARN]",
    "⚪": "[N/A]",
    "🔵": "[INFO]",
    "✅": "[OK]",
    "❌": "[FAIL]",
    "⚠️": "[WARN]",
    "ℹ️": "[INFO]",
    "🏆": "[VALIDATED]",
    "🎯": "[TARGET]",
    "🔬": "",
    "📊": "",
    "💊": "",
    "🧬": "",
    "👥": "",
    "🔢": "",
    "📝": "",
    "💡": "",
    "📸": "",
    "🗑️": "",
    "📁": "",
    "💎": "",
    "⚡": "",
}


def sanitize_text(text: Any) -> str:
    """Sanitize arbitrary text or objects for standard PDF Latin-1 compatibility.
    
    Replaces scientific notation (Å, ±, ΔG, etc.) and emojis, decomposing or
    dropping any remaining non-Latin-1 codepoints to prevent encoding errors.
    """
    if text is None:
        return ""
    s = str(text)
    if not s:
        return ""

    for src, dst in UNICODE_REPLACEMENTS.items():
        if src in s:
            s = s.replace(src, dst)

    clean_chars: list[str] = []
    for ch in s:
        try:
            ch.encode("latin-1")
            clean_chars.append(ch)
        except UnicodeEncodeError:
            # Attempt NFKD decomposition
            decomposed = unicodedata.normalize("NFKD", ch)
            encoded_part = False
            for d in decomposed:
                try:
                    d.encode("latin-1")
                    clean_chars.append(d)
                    encoded_part = True
                except UnicodeEncodeError:
                    pass
            if not encoded_part:
                clean_chars.append("?")

    return "".join(clean_chars)


def decode_image_data(image_input: Any) -> Optional[io.BytesIO]:
    """Decode base64 string, data URL, or bytes into a valid PIL-verified BytesIO stream."""
    if not image_input:
        return None

    try:
        raw_bytes: bytes
        if isinstance(image_input, str):
            clean_str = image_input.strip()
            if not clean_str:
                return None
            if "," in clean_str and "base64" in clean_str.split(",", 1)[0]:
                clean_str = clean_str.split(",", 1)[1]
            raw_bytes = base64.b64decode(clean_str)
        elif isinstance(image_input, (bytes, bytearray)):
            raw_bytes = bytes(image_input)
        elif isinstance(image_input, io.BytesIO):
            image_input.seek(0)
            raw_bytes = image_input.read()
        else:
            return None

        # Verify image using PIL
        buf = io.BytesIO(raw_bytes)
        with Image.open(buf) as img:
            img.verify()
        buf.seek(0)
        return buf
    except Exception:
        return None


# ── PDF Document Base ─────────────────────────────────────────────────────────


class BioBioPDF(FPDF):
    """Customized PDF document styling for BioBio reports."""

    def __init__(self, doc_category: str = "Molecular Report") -> None:
        super().__init__(orientation="P", unit="mm", format="A4")
        self.doc_category = sanitize_text(doc_category)
        self.set_auto_page_break(auto=True, margin=18)
        self.set_margins(left=15, top=18, right=15)
        self.alias_nb_pages()

    def header(self) -> None:
        """Render top page header."""
        self.set_font("Helvetica", "B", 8)
        self.set_text_color(100, 116, 139)  # Slate grey
        self.cell(0, 5, sanitize_text(f"BioBio Platform  |  {self.doc_category}"), align="L")
        self.set_font("Helvetica", "", 8)
        self.cell(0, 5, datetime.now().strftime("%Y-%m-%d %H:%M"), align="R", new_x="LMARGIN", new_y="NEXT")
        self.set_draw_color(226, 232, 240)
        self.set_line_width(0.3)
        self.line(15, self.get_y() + 1, self.w - 15, self.get_y() + 1)
        self.ln(4)

    def footer(self) -> None:
        """Render bottom page footer with page numbers."""
        self.set_y(-14)
        self.set_draw_color(226, 232, 240)
        self.set_line_width(0.3)
        self.line(15, self.get_y(), self.w - 15, self.get_y())
        self.ln(2)
        self.set_font("Helvetica", "", 8)
        self.set_text_color(148, 163, 184)
        self.cell(0, 5, "BioBio Molecular Simulation & Variant Analysis", align="L")
        self.cell(0, 5, f"Page {self.page_no()} of {{nb}}", align="R")

    # ── High-Level Content Layout Helpers ─────────────────────────────────────

    def add_report_banner(self, title: str, subtitle: Optional[str] = None) -> None:
        """Draw top banner with title and optional subtitle."""
        self.set_fill_color(248, 250, 252)
        self.set_draw_color(203, 213, 225)
        self.set_line_width(0.4)
        banner_start = self.get_y()

        self.set_font("Helvetica", "B", 16)
        self.set_text_color(15, 23, 42)  # Dark slate
        sanitized_title = sanitize_text(title)
        self.cell(0, 8, sanitized_title, new_x="LMARGIN", new_y="NEXT")

        if subtitle:
            self.set_font("Helvetica", "", 9)
            self.set_text_color(71, 85, 105)
            self.cell(0, 5, sanitize_text(subtitle), new_x="LMARGIN", new_y="NEXT")

        self.ln(3)

    def add_section_header(self, title: str) -> None:
        """Add styled section header with blue accent bar."""
        # Ensure we don't start header right at the bottom of page
        if self.get_y() > self.page_break_trigger - 25:
            self.add_page()
        else:
            self.ln(3)

        self.set_draw_color(37, 99, 235)  # Royal blue accent
        self.set_line_width(0.8)
        y = self.get_y()
        self.line(self.l_margin, y, self.l_margin + 5, y)

        self.set_font("Helvetica", "B", 11)
        self.set_text_color(30, 41, 59)
        self.cell(0, 6, sanitize_text(f"  {title}"), new_x="LMARGIN", new_y="NEXT")
        self.ln(1)

    def add_key_value_grid(
        self, pairs: Sequence[Tuple[str, Any]], cols: int = 2
    ) -> None:
        """Display clean key-value summary metrics, skipping None/empty gracefully."""
        valid_pairs = [(k, v) for k, v in pairs if v is not None and str(v).strip() != ""]
        if not valid_pairs:
            return

        col_w = self.epw / cols
        row_height = 6

        self.set_font("Helvetica", "B", 8)
        self.set_text_color(100, 116, 139)

        for i in range(0, len(valid_pairs), cols):
            chunk = valid_pairs[i : i + cols]
            # Print keys
            self.set_font("Helvetica", "B", 8)
            self.set_text_color(100, 116, 139)
            for key, _ in chunk:
                self.cell(col_w, 4, sanitize_text(key.upper()))
            self.ln(4)

            # Print values
            self.set_font("Helvetica", "", 9)
            self.set_text_color(15, 23, 42)
            for _, val in chunk:
                val_str = sanitize_text(str(val))
                self.cell(col_w, row_height, val_str)
            self.ln(row_height + 2)

    def add_table(
        self,
        headers: list[str],
        rows: Sequence[Sequence[Any]],
        col_widths: Optional[list[float]] = None,
        max_rows: int = 30,
    ) -> None:
        """Render a clean structured table with alternating row colors."""
        if not rows:
            return

        if not col_widths or len(col_widths) != len(headers):
            col_widths = [self.epw / len(headers)] * len(headers)

        # Table Header
        self.set_font("Helvetica", "B", 8)
        self.set_fill_color(241, 245, 249)
        self.set_text_color(51, 65, 85)
        self.set_draw_color(226, 232, 240)
        self.set_line_width(0.2)

        for i, header in enumerate(headers):
            self.cell(
                col_widths[i],
                6,
                sanitize_text(header),
                border=1,
                align="C",
                fill=True,
            )
        self.ln(6)

        # Table Rows
        self.set_font("Helvetica", "", 8)
        self.set_text_color(15, 23, 42)

        rendered_rows = list(rows)[:max_rows]
        for row_idx, row in enumerate(rendered_rows):
            if self.get_y() > self.page_break_trigger - 10:
                self.add_page()
                # Re-render header on new page
                self.set_font("Helvetica", "B", 8)
                self.set_fill_color(241, 245, 249)
                self.set_text_color(51, 65, 85)
                for i, header in enumerate(headers):
                    self.cell(
                        col_widths[i],
                        6,
                        sanitize_text(header),
                        border=1,
                        align="C",
                        fill=True,
                    )
                self.ln(6)
                self.set_font("Helvetica", "", 8)
                self.set_text_color(15, 23, 42)

            fill = (row_idx % 2 == 1)
            if fill:
                self.set_fill_color(248, 250, 252)

            for col_idx, cell_val in enumerate(row):
                w = col_widths[col_idx]
                val_str = "" if cell_val is None else sanitize_text(str(cell_val))
                self.cell(
                    w,
                    5,
                    val_str,
                    border=1,
                    align="C",
                    fill=fill,
                )
            self.ln(5)

        if len(rows) > max_rows:
            self.set_font("Helvetica", "I", 7)
            self.set_text_color(100, 116, 139)
            self.cell(
                0,
                5,
                sanitize_text(f"(Showing first {max_rows} of {len(rows)} entries)"),
                align="R",
            )
            self.ln(4)

        self.ln(2)

    def add_callout(self, text: str, label: str = "NOTE") -> None:
        """Render a shaded info/warning callout block."""
        if not text or not str(text).strip():
            return
        sanitized = sanitize_text(text)
        self.set_fill_color(241, 245, 249)
        self.set_draw_color(203, 213, 225)
        self.set_line_width(0.3)
        self.set_font("Helvetica", "B", 8)
        self.set_text_color(30, 41, 59)

        callout_w = self.epw
        self.cell(callout_w, 5, sanitize_text(f"[{label}]"), fill=True, border="TLR")
        self.ln(5)
        self.set_font("Helvetica", "", 8)
        self.multi_cell(callout_w, 4.5, sanitized, fill=True, border="BLR")
        self.ln(3)

    def add_screenshots_section(
        self,
        screenshots: Sequence[Any],
        section_title: str = "3D Structure Views & Snapshots",
        default_caption: str = "3D Structure View",
    ) -> None:
        """Decode base64 PNG data URLs and embed images with sequential captions."""
        if not screenshots:
            return

        valid_buffers: list[io.BytesIO] = []
        for item in screenshots:
            buf = decode_image_data(item)
            if buf:
                valid_buffers.append(buf)

        if not valid_buffers:
            return

        self.add_section_header(section_title)

        target_w = min(130.0, self.epw)
        for i, buf in enumerate(valid_buffers, start=1):
            try:
                buf.seek(0)
                with Image.open(buf) as pil_img:
                    orig_w, orig_h = pil_img.size
                ratio = orig_h / max(orig_w, 1)
                rendered_h = target_w * ratio

                # Check if image + caption will fit on current page
                if self.get_y() + rendered_h + 16 > self.page_break_trigger:
                    self.add_page()

                x_pos = self.l_margin + (self.epw - target_w) / 2
                y_pos = self.get_y()

                buf.seek(0)
                self.image(buf, x=x_pos, y=y_pos, w=target_w)
                self.set_y(y_pos + rendered_h + 2)

                self.set_font("Helvetica", "I", 8)
                self.set_text_color(71, 85, 105)
                caption = sanitize_text(f"Figure {i}: {default_caption}")
                self.cell(0, 5, caption, align="C", new_x="LMARGIN", new_y="NEXT")
                self.ln(4)
            except Exception:
                continue


# ── Builders for Application Tabs ─────────────────────────────────────────────


def build_analysis_pdf(data: dict, screenshots: list[str]) -> bytes:
    """Build PDF report for Tab 1: Structural Analysis & Variant Comparison.
    
    Args:
        data: Analysis dictionary containing protein IDs, alignment, SASA and contact metrics.
        screenshots: List of base64 PNG data URLs of 3D viewer angles.

    Returns:
        Raw PDF document bytes.
    """
    data = data or {}
    screenshots = screenshots or []

    pdf = BioBioPDF(doc_category="Structural Analysis & Variant Comparison")
    pdf.add_page()

    title = data.get("title") or "Structural Analysis & Variant Comparison Report"
    subtitle = data.get("subtitle") or "Pairwise Alignment, Residue Contacts & SASA Profile"
    pdf.add_report_banner(title, subtitle)

    # 1. Target & Input Metadata
    meta_pairs: list[Tuple[str, Any]] = []
    if data.get("wt_id") or data.get("healthy_id") or data.get("h_id"):
        meta_pairs.append(("Wild-Type (WT)", data.get("wt_id") or data.get("healthy_id") or data.get("h_id")))
    if data.get("wt_chain") or data.get("h_chain"):
        meta_pairs.append(("WT Chain", data.get("wt_chain") or data.get("h_chain")))
    if data.get("mt_id") or data.get("mutant_id") or data.get("m_id"):
        meta_pairs.append(("Mutant Variant (MT)", data.get("mt_id") or data.get("mutant_id") or data.get("m_id")))
    if data.get("mt_chain") or data.get("m_chain"):
        meta_pairs.append(("MT Chain", data.get("mt_chain") or data.get("m_chain")))
    if data.get("coordinate_mode"):
        meta_pairs.append(("Coordinate Mode", data.get("coordinate_mode")))

    if meta_pairs:
        pdf.add_section_header("Protein Target Information")
        pdf.add_key_value_grid(meta_pairs, cols=2)

    # 2. Sequence & Structural Alignment Metrics
    align_data = data.get("alignment") or {}
    align_pairs: list[Tuple[str, Any]] = []
    if "score" in align_data or "alignment_score" in data:
        align_pairs.append(("Alignment Score", round(float(align_data.get("score", data.get("alignment_score", 0))), 1)))
    if "identity_pct" in align_data or "sequence_identity" in data:
        pct = align_data.get("identity_pct", data.get("sequence_identity"))
        align_pairs.append(("Sequence Identity", f"{float(pct):.1f}%"))
    if "mutations" in align_data or "mutations_count" in data:
        align_pairs.append(("Mutations", align_data.get("mutations", data.get("mutations_count"))))
    if "gaps" in align_data or "indels_count" in data:
        align_pairs.append(("Gaps & Indels", align_data.get("gaps", data.get("indels_count"))))
    if "ca_rmsd" in align_data or "rmsd" in data:
        rmsd_val = align_data.get("ca_rmsd", data.get("rmsd"))
        if rmsd_val is not None:
            align_pairs.append(("Backbone C-alpha RMSD", f"{float(rmsd_val):.3f} Angstrom"))

    if align_pairs:
        pdf.add_section_header("Sequence & Structural Alignment")
        pdf.add_key_value_grid(align_pairs, cols=3)

    # Text alignment snippet if provided
    align_text = align_data.get("text") or data.get("alignment_text")
    if align_text and str(align_text).strip():
        pdf.set_font("Courier", "", 7)
        pdf.set_fill_color(248, 250, 252)
        pdf.set_text_color(30, 41, 59)
        pdf.multi_cell(pdf.epw, 3.5, sanitize_text(str(align_text)[:2000]), fill=True, border=1)
        pdf.ln(2)

    # 3. Structural Comparison Table (WT vs MT)
    comp_rows = data.get("comparison_rows") or data.get("comparison_data") or []
    if hasattr(comp_rows, "to_dict"):
        comp_rows = comp_rows.to_dict(orient="records")

    if comp_rows and isinstance(comp_rows, list):
        pdf.add_section_header("Residue-Level Comparative Profile (WT vs. MT)")
        headers = ["WT #", "MT #", "WT Res", "MT Res", "WT SASA", "MT SASA", "dSASA", "Status"]
        col_w = [18, 18, 20, 20, 25, 25, 25, 29]
        formatted_rows = []
        for r in comp_rows:
            if isinstance(r, dict):
                formatted_rows.append([
                    r.get("res_num_h", r.get("WT Residue #", "-")),
                    r.get("res_num_m", r.get("MT Residue #", "-")),
                    r.get("WT_Res", r.get("WT Residue", "-")),
                    r.get("MT_Res", r.get("MT Residue", "-")),
                    r.get("SASA_H", r.get("WT SASA", "-")),
                    r.get("SASA_M", r.get("MT SASA", "-")),
                    r.get("SASA_Delta", r.get("Delta SASA", "-")),
                    r.get("Status", "-"),
                ])
        if formatted_rows:
            pdf.add_table(headers, formatted_rows, col_widths=col_w, max_rows=25)

    # 4. Single Chain Structural Results (if comparison not present)
    single_res = data.get("results") or data.get("structural_results") or []
    if hasattr(single_res, "to_dict"):
        single_res = single_res.to_dict(orient="records")
    if not comp_rows and single_res and isinstance(single_res, list):
        pdf.add_section_header("Residue Structural Analysis")
        headers = ["Chain", "Res #", "Residue", "Class", "Contact Dist (Angstrom)", "SASA (Angstrom^2)"]
        col_w = [20, 25, 30, 35, 35, 35]
        s_rows = []
        for r in single_res:
            if isinstance(r, dict):
                s_rows.append([
                    r.get("chain", "-"),
                    r.get("res_num", "-"),
                    r.get("res_name", "-"),
                    r.get("class", "-"),
                    r.get("min_dist", "-"),
                    r.get("sasa", "-"),
                ])
        if s_rows:
            pdf.add_table(headers, s_rows, col_widths=col_w, max_rows=25)

    # 5. Summary / Notes callout
    if data.get("summary") or data.get("notes"):
        pdf.add_callout(data.get("summary") or data.get("notes"), label="SUMMARY")

    # 6. Screenshots
    pdf.add_screenshots_section(
        screenshots,
        section_title="3D Structural Views",
        default_caption="3D Structure View",
    )

    return bytes(pdf.output())


def build_docking_pdf(data: dict, screenshots: list[str]) -> bytes:
    """Build PDF report for Tab 2: Molecular Docking Simulation (AutoDock Vina).
    
    Args:
        data: Docking results dictionary (poses, scores, validation status, properties).
        screenshots: List of base64 PNG data URLs of 3D binding poses.

    Returns:
        Raw PDF document bytes.
    """
    data = data or {}
    screenshots = screenshots or []

    pdf = BioBioPDF(doc_category="Molecular Docking Simulation")
    pdf.add_page()

    title = data.get("title") or "Molecular Docking Simulation Report"
    lig_name = data.get("ligand_name") or data.get("ligand_id")
    target_id = data.get("target_id") or data.get("target_pdb")
    subtitle = f"AutoDock Vina Simulation: {lig_name} in {target_id}" if lig_name and target_id else "AutoDock Vina Simulation & Pose Analysis"
    pdf.add_report_banner(title, subtitle)

    # 1. Docking Setup & Target Metadata
    kpi_pairs: list[Tuple[str, Any]] = []
    if target_id:
        kpi_pairs.append(("Target Receptor", target_id))
    if lig_name:
        kpi_pairs.append(("Ligand Compound", lig_name))
    if data.get("workflow_mode"):
        kpi_pairs.append(("Workflow Mode", data.get("workflow_mode")))
    if data.get("validation_status"):
        kpi_pairs.append(("Validation Status", data.get("validation_status")))

    best_score = data.get("lowest_vina_score", data.get("best_score"))
    if best_score is not None:
        kpi_pairs.append(("Lowest Vina Score", f"{float(best_score):.2f} kcal/mol"))

    poses = data.get("poses") or []
    if poses:
        kpi_pairs.append(("Total Poses", len(poses)))

    clusters = data.get("clusters") or []
    if clusters:
        kpi_pairs.append(("Pose Clusters", len(clusters)))

    if kpi_pairs:
        pdf.add_section_header("Docking Simulation Summary")
        pdf.add_key_value_grid(kpi_pairs, cols=2)

    # Validation Details / Rationale Callout
    val_reason = data.get("validation_details") or data.get("validation_reason")
    if val_reason:
        pdf.add_callout(val_reason, label="VALIDATION RATIONALE")

    # 2. Physicochemical / Lipinski Rule of 5 Properties
    props = data.get("properties") or data.get("drug_properties") or {}
    prop_pairs: list[Tuple[str, Any]] = []
    if "mw" in props:
        prop_pairs.append(("Molecular Weight", f"{props['mw']:.1f} g/mol"))
    if "logp" in props:
        prop_pairs.append(("LogP", f"{props['logp']:.2f}"))
    if "hbd" in props:
        prop_pairs.append(("H-Bond Donors", props["hbd"]))
    if "hba" in props:
        prop_pairs.append(("H-Bond Acceptors", props["hba"]))
    if "rotatable_bonds" in props:
        prop_pairs.append(("Rotatable Bonds", props["rotatable_bonds"]))
    if "lipinski_violations" in props:
        prop_pairs.append(("Lipinski Violations", props["lipinski_violations"]))

    if prop_pairs:
        pdf.add_section_header("Ligand Physicochemical & Drug-Likeness Profile")
        pdf.add_key_value_grid(prop_pairs, cols=3)

    # 3. Poses Summary Table
    if poses and isinstance(poses, list):
        pdf.add_section_header("Binding Pose Conformations")
        headers = ["Rank", "Binding Affinity (kcal/mol)", "RMSD l.b. (Angstrom)", "RMSD u.b. (Angstrom)"]
        col_w = [25, 60, 50, 45]
        p_rows = []
        for idx, p in enumerate(poses, start=1):
            if isinstance(p, dict):
                p_rows.append([
                    p.get("rank", idx),
                    f"{float(p.get('score', p.get('affinity', 0))):.2f}",
                    p.get("rmsd_lb", p.get("rmsd_lower", "-")),
                    p.get("rmsd_ub", p.get("rmsd_upper", "-")),
                ])
        if p_rows:
            pdf.add_table(headers, p_rows, col_widths=col_w, max_rows=15)

    # 4. Pharmacological Interactions
    inter = data.get("interactions")
    if inter and isinstance(inter, (dict, list)):
        pdf.add_section_header("Receptor-Ligand Interaction Profile")
        if isinstance(inter, dict):
            i_pairs = [(k, v) for k, v in inter.items() if v]
            pdf.add_key_value_grid(i_pairs, cols=2)
        elif isinstance(inter, list):
            i_rows = [[i.get("type", "-"), i.get("residue", "-"), i.get("distance", "-")] for i in inter if isinstance(i, dict)]
            if i_rows:
                pdf.add_table(["Interaction Type", "Receptor Residue", "Distance (Angstrom)"], i_rows)

    # 5. Screenshots
    pdf.add_screenshots_section(
        screenshots,
        section_title="3D Binding Pose Visualizations",
        default_caption="3D Docked Pose View",
    )

    return bytes(pdf.output())


def build_md_pdf(data: dict, plots: dict, screenshots: list[str]) -> bytes:
    """Build PDF report for Tab 3: Molecular Dynamics Simulation (GROMACS).
    
    Args:
        data: MD simulation metrics, pose metadata, and thermodynamic summaries.
        plots: Dict of matplotlib Figure objects (or bytes/streams) for RMSD, RMSF, etc.
        screenshots: List of base64 PNG data URLs of 3D trajectory snapshots.

    Returns:
        Raw PDF document bytes.
    """
    data = data or {}
    plots = plots or {}
    screenshots = screenshots or []

    pdf = BioBioPDF(doc_category="Molecular Dynamics Simulation")
    pdf.add_page()

    title = data.get("title") or "Molecular Dynamics Simulation Report"
    lig_name = data.get("ligand_name") or data.get("ligand_id")
    target_id = data.get("target_pdb") or data.get("target_id")
    subtitle = f"GROMACS MD Production: {lig_name} in {target_id}" if lig_name and target_id else "GROMACS 2026.3 GPU-Accelerated Production Dynamics"
    pdf.add_report_banner(title, subtitle)

    # 1. Simulation Parameters & Environment
    sim_pairs: list[Tuple[str, Any]] = []
    if target_id:
        sim_pairs.append(("Target Receptor", target_id))
    if lig_name:
        sim_pairs.append(("Ligand Name", lig_name))
    if data.get("gromacs_version") or data.get("engine"):
        sim_pairs.append(("MD Engine", data.get("gromacs_version") or data.get("engine")))
    if data.get("hardware"):
        sim_pairs.append(("Hardware Acceleration", data.get("hardware")))
    if data.get("unit_cell"):
        sim_pairs.append(("Unit Cell Geometry", data.get("unit_cell")))
    if data.get("duration_ns") or data.get("production_duration"):
        dur = data.get("duration_ns") or data.get("production_duration")
        sim_pairs.append(("Production Duration", f"{dur} ns" if not str(dur).endswith("ns") else str(dur)))
    if data.get("ensemble"):
        sim_pairs.append(("Ensemble", data.get("ensemble")))
    if data.get("temperature_k"):
        sim_pairs.append(("Temperature", f"{data['temperature_k']} K"))
    if data.get("pressure_bar"):
        sim_pairs.append(("Pressure", f"{data['pressure_bar']} bar"))

    if sim_pairs:
        pdf.add_section_header("Simulation Parameters & Setup")
        pdf.add_key_value_grid(sim_pairs, cols=3)

    # 2. Key Dynamic & Stability Metrics
    kpi_pairs: list[Tuple[str, Any]] = []
    if "mean_ligand_rmsd" in data or "mean_rmsd" in data:
        rmsd = data.get("mean_ligand_rmsd", data.get("mean_rmsd"))
        kpi_pairs.append(("Mean Ligand RMSD", f"{float(rmsd):.3f} Angstrom"))
    if "max_rmsd" in data:
        kpi_pairs.append(("Max Ligand RMSD", f"{float(data['max_rmsd']):.3f} Angstrom"))
    if "mean_hbond_count" in data or "mean_hbonds" in data:
        hb = data.get("mean_hbond_count", data.get("mean_hbonds"))
        kpi_pairs.append(("Mean H-Bond Count", f"{float(hb):.1f}"))
    if "duration_ps" in data:
        kpi_pairs.append(("Total Time Steps", f"{data['duration_ps']} ps"))

    if kpi_pairs:
        pdf.add_section_header("Trajectory Stability Metrics")
        pdf.add_key_value_grid(kpi_pairs, cols=2)

    # 3. MM-PB(GB)SA Energetics Summary
    mmpbsa = data.get("mmpbsa") or data.get("mmpbsa_res") or {}
    mmpbsa_pairs: list[Tuple[str, Any]] = []
    if "dg_bind" in mmpbsa or "total_energy" in mmpbsa:
        dg = mmpbsa.get("dg_bind", mmpbsa.get("total_energy"))
        mmpbsa_pairs.append(("Binding Free Energy (dG)", f"{float(dg):.2f} kcal/mol"))
    if "van_der_waals" in mmpbsa or "vdw" in mmpbsa:
        vdw = mmpbsa.get("van_der_waals", mmpbsa.get("vdw"))
        mmpbsa_pairs.append(("Van der Waals (vdW)", f"{float(vdw):.2f} kcal/mol"))
    if "electrostatic" in mmpbsa or "elec" in mmpbsa:
        elec = mmpbsa.get("electrostatic", mmpbsa.get("elec"))
        mmpbsa_pairs.append(("Electrostatic (Elec)", f"{float(elec):.2f} kcal/mol"))
    if "polar_solvation" in mmpbsa or "polar" in mmpbsa:
        polar = mmpbsa.get("polar_solvation", mmpbsa.get("polar"))
        mmpbsa_pairs.append(("Polar Solvation", f"{float(polar):.2f} kcal/mol"))
    if "nonpolar_solvation" in mmpbsa or "sasa" in mmpbsa:
        nonpolar = mmpbsa.get("nonpolar_solvation", mmpbsa.get("sasa"))
        mmpbsa_pairs.append(("Non-polar Solvation (SASA)", f"{float(nonpolar):.2f} kcal/mol"))

    if mmpbsa_pairs:
        pdf.add_section_header("MM-PB(GB)SA Binding Free Energy Decomposition")
        pdf.add_key_value_grid(mmpbsa_pairs, cols=2)

    # 4. Embedded Matplotlib Figures from `plots` dict
    if plots and isinstance(plots, dict):
        valid_plots = []
        for name, plot_obj in plots.items():
            if plot_obj is None:
                continue
            # If matplotlib Figure object: convert via fig.savefig(buf, format='png', dpi=150)
            if hasattr(plot_obj, "savefig"):
                buf = io.BytesIO()
                plot_obj.savefig(buf, format="png", dpi=150, bbox_inches="tight")
                buf.seek(0)
                valid_plots.append((name, buf))
            elif isinstance(plot_obj, (bytes, bytearray)):
                valid_plots.append((name, io.BytesIO(bytes(plot_obj))))
            elif isinstance(plot_obj, io.BytesIO):
                plot_obj.seek(0)
                valid_plots.append((name, plot_obj))
            elif isinstance(plot_obj, str):
                buf = decode_image_data(plot_obj)
                if buf:
                    valid_plots.append((name, buf))

        if valid_plots:
            pdf.add_section_header("Trajectory Analysis Charts")
            target_w = min(140.0, pdf.epw)
            for plot_name, buf in valid_plots:
                try:
                    buf.seek(0)
                    with Image.open(buf) as pil_img:
                        pw, ph = pil_img.size
                    ratio = ph / max(pw, 1)
                    rendered_h = target_w * ratio

                    if pdf.get_y() + rendered_h + 16 > pdf.page_break_trigger:
                        pdf.add_page()

                    x_pos = pdf.l_margin + (pdf.epw - target_w) / 2
                    y_pos = pdf.get_y()

                    buf.seek(0)
                    pdf.image(buf, x=x_pos, y=y_pos, w=target_w)
                    pdf.set_y(y_pos + rendered_h + 2)

                    pdf.set_font("Helvetica", "I", 8)
                    pdf.set_text_color(71, 85, 105)
                    label = sanitize_text(str(plot_name).replace("_", " ").title())
                    pdf.cell(0, 5, sanitize_text(f"Chart: {label}"), align="C", new_x="LMARGIN", new_y="NEXT")
                    pdf.ln(3)
                except Exception:
                    continue

    # 5. Screenshots
    pdf.add_screenshots_section(
        screenshots,
        section_title="3D Molecular Conformations",
        default_caption="MD Trajectory Conformation",
    )

    return bytes(pdf.output())
