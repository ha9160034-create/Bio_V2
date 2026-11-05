# Protein-Ligand Interaction and Pharmacological Profiler
"""
وحدة تحليل التفاعلات الكيميائية بين الدواء والبروتين، وكشف الروابط الهيدروجينية،
وتقييم الخصائص الصيدلانية ومعايير ليبينسكي (Lipinski's Rule of 5).
"""

from io import StringIO
import numpy as np
from scipy.spatial.distance import cdist
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors
from rdkit.Chem import FilterCatalog

from md.topology import autodock_type_to_element


# Transparent structural alerts for common colloidal-aggregator chemotypes.
# These are heuristics, not the Shoichet Aggregator Advisor model.
_AGGREGATOR_SMARTS = {
    "rhodanine": "O=C1NC(=S)SC1",
    "quinone": "O=C1C=CC(=O)C=C1",
    "catechol": "Oc1ccccc1O",
    "flavone": "O=c1cc(-c2ccccc2)oc2ccccc12",
}


def _aggregator_alerts(mol) -> list[str]:
    alerts = []
    for name, smarts in _AGGREGATOR_SMARTS.items():
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is not None and mol.HasSubstructMatch(pattern):
            alerts.append(name)
    return alerts

def calculate_drug_properties(smiles: str) -> dict:
    """
    حساب الخصائص الصيدلانية والكيميائية للدواء وفق قواعد ليبينسكي (Lipinski Rule of Five)
    لتقييم مدى ملائمة الجزيء كدواء فموي (Drug-likeness / ADME).
    """
    if not smiles or not smiles.strip():
        return {}
    
    mol = Chem.MolFromSmiles(smiles.strip())
    if mol is None:
        return {"error": "صيغة SMILES غير صالحة كيميائياً"}

    mw = round(Descriptors.MolWt(mol), 2)
    logp = round(Descriptors.MolLogP(mol), 2)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    rotb = Lipinski.NumRotatableBonds(mol)
    tpsa = round(rdMolDescriptors.CalcTPSA(mol), 2)
    formula = rdMolDescriptors.CalcMolFormula(mol)

    # التحقق من قواعد ليبينسكي الخمس (Lipinski Violations)
    violations = 0
    if mw > 500: violations += 1
    if logp > 5.0: violations += 1
    if hbd > 5: violations += 1
    if hba > 10: violations += 1
    
    lipinski_pass = (violations <= 1)

    pains_name = ""
    try:
        params = FilterCatalog.FilterCatalogParams()
        params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
        catalog = FilterCatalog.FilterCatalog(params)
        match = catalog.GetFirstMatch(mol)
        if match is not None:
            pains_name = match.GetDescription()
    except Exception:
        pains_name = ""

    # Wiki-aligned, explicitly heuristic aggregation screen.  It is not the
    # Shoichet Aggregator Advisor and must not be presented as an experimental
    # result.  logP > 3.5 follows wiki/chemical-screening-artifacts.md.
    aggregation_alerts = _aggregator_alerts(mol)
    aggregation_flags = []
    if logp > 3.5:
        aggregation_flags.append(f"logP {logp:.2f} > 3.5")
    if mw >= 300:
        aggregation_flags.append(f"MW {mw:.1f} ≥ 300")
    if tpsa < 75:
        aggregation_flags.append(f"TPSA {tpsa:.2f} < 75")
    aggregation_risk = bool(aggregation_alerts) or (
        logp > 3.5 and mw >= 300 and tpsa < 75
    )
    aggregation_note = ""
    if aggregation_risk:
        parts = []
        if aggregation_alerts:
            parts.append("structural alerts: " + ", ".join(aggregation_alerts))
        if aggregation_flags:
            parts.append("properties: " + ", ".join(aggregation_flags))
        aggregation_note = (
            "Heuristic only (" + "; ".join(parts) + "). "
            "Not Aggregator Advisor; confirm experimentally (e.g. detergent sensitivity)."
        )

    return {
        "formula": formula,
        "molecular_weight": mw,
        "logp": logp,
        "hbd": hbd,
        "hba": hba,
        "rotatable_bonds": rotb,
        "tpsa": tpsa,
        "lipinski_violations": violations,
        "lipinski_pass": lipinski_pass,
        "pains_alert": bool(pains_name),
        "pains_name": pains_name,
        "aggregation_risk": aggregation_risk,
        "aggregation_note": aggregation_note,
        "aggregation_alerts": aggregation_alerts,
        "aggregation_flags": aggregation_flags,
        "aggregation_logp_threshold": 3.5,
        "drug_likeness": "مطابق لقواعد ليبينسكي" if lipinski_pass else f"تنبيه: يوجد {violations} تجاوزات لقواعد ليبينسكي"
    }

def _parse_atoms_from_block(block_text: str, is_ligand: bool = False) -> list[dict]:
    """استخراج إحداثيات وخصائص الذرات من نصوص PDB أو PDBQT.

    PDBQT stores an AutoDock atom *type* in columns 78-79 (for example ``A``,
    ``OA``, ``NA``, ``HD``), not a PDB element.  Reading those columns as an
    element turns aromatic carbon into ``A`` and can drop it from carbon-contact
    counts.  Ligand records therefore go through the documented AutoDock type
    mapping; receptor PDB records keep their real element columns.
    """
    atoms = []
    if not block_text:
        return atoms

    for line in block_text.splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            if len(line) < 54:
                continue
            try:
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                chain_id = line[21].strip() or "A"
                res_num = line[22:26].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])

                element = ""
                if is_ligand:
                    raw_type = line[77:].strip() if len(line) >= 78 else ""
                    if not raw_type:
                        raw_type = line[76:78].strip() if len(line) >= 78 else ""
                    element = autodock_type_to_element(raw_type, atom_name).upper()
                elif len(line) >= 78:
                    element = line[76:78].strip().upper()
                if not element:
                    alpha = "".join(c for c in atom_name if c.isalpha())
                    element = alpha[:1].upper() if alpha else "C"

                atoms.append({
                    "atom_name": atom_name,
                    "res_name": res_name,
                    "chain": chain_id,
                    "res_num": res_num,
                    "coord": np.array([x, y, z], dtype=np.float32),
                    "element": element,
                    "is_ligand": is_ligand
                })
            except (ValueError, IndexError):
                continue

    return atoms

CONTACTS_DISCLAIMER = (
    "Distance-based contact candidates only. They are not hydrogen-bond assignments, "
    "binding free-energy components, or experimental interaction evidence."
)

def analyze_protein_ligand_interactions(
    receptor_pdb_text: str,
    ligand_pdbqt_block: str,
    polar_cutoff: float = 3.5,
    contact_cutoff: float = 4.0
) -> dict:
    """
    استكشاف أولي قائم على المسافات الهندسية للتماسات بين وضعية الدواء والمستقبل البروتيني:
    1. مرشحات التماس القطبي (Polar contact candidates): ذرات O/N/F بمسافة <= 3.5 Å.
    2. مرشحات التماس الكربوني (Carbon contact candidates): ذرات كربون بمسافة <= 4.0 Å.
    3. أحماض الجيب المحيطة المتماسة (Interacting Residues): أحماض بمسافة تماس <= 4.0 Å.
    """
    rec_atoms = _parse_atoms_from_block(receptor_pdb_text, is_ligand=False)
    lig_atoms = _parse_atoms_from_block(ligand_pdbqt_block, is_ligand=True)

    if not rec_atoms or not lig_atoms:
        return {
            "polar_contact_candidates": [],
            "carbon_contact_candidates": [],
            "interacting_residues": [],
            "polar_contact_count": 0,
            "carbon_contact_count": 0,
            "disclaimer": CONTACTS_DISCLAIMER
        }

    rec_coords = np.array([a["coord"] for a in rec_atoms])
    lig_coords = np.array([a["coord"] for a in lig_atoms])

    dist_matrix = cdist(rec_coords, lig_coords)

    polar_contacts = []
    carbon_contacts = []
    interacting_res_map = {}

    polar_elements = {"O", "N", "F"}

    for r_idx, r_atom in enumerate(rec_atoms):
        for l_idx, l_atom in enumerate(lig_atoms):
            dist = float(dist_matrix[r_idx, l_idx])

            # تسجيل الأحماض المتفاعلة بشكل عام ضمن النطاق القريب
            if dist <= contact_cutoff:
                res_key = f"{r_atom['res_name']} {r_atom['res_num']} (Chain {r_atom['chain']})"
                if res_key not in interacting_res_map or dist < interacting_res_map[res_key]["min_distance"]:
                    interacting_res_map[res_key] = {
                        "res_name": r_atom["res_name"],
                        "res_num": r_atom["res_num"],
                        "chain": r_atom["chain"],
                        "min_distance": round(dist, 2)
                    }

            # 1. كشف مرشحات التماس القطبي (Polar Contact Candidates)
            if dist <= polar_cutoff:
                if r_atom["element"] in polar_elements and l_atom["element"] in polar_elements:
                    polar_contacts.append({
                        "residue": f"{r_atom['res_name']}{r_atom['res_num']}",
                        "res_num": r_atom["res_num"],
                        "chain": r_atom["chain"],
                        "rec_atom": r_atom["atom_name"],
                        "lig_atom": l_atom["atom_name"],
                        "distance_angstrom": round(dist, 2)
                    })

            # 2. كشف مرشحات التماس الكربوني (Carbon-Carbon Contacts)
            if dist <= contact_cutoff:
                if r_atom["element"] == "C" and l_atom["element"] == "C":
                    carbon_contacts.append({
                        "residue": f"{r_atom['res_name']}{r_atom['res_num']}",
                        "res_num": r_atom["res_num"],
                        "chain": r_atom["chain"],
                        "distance_angstrom": round(dist, 2)
                    })

    # إزالة التكرارات
    unique_polar = []
    seen_polar = set()
    for pc in polar_contacts:
        key = (pc["residue"], pc["rec_atom"], pc["lig_atom"])
        if key not in seen_polar:
            seen_polar.add(key)
            unique_polar.append(pc)

    unique_interacting_res = sorted(
        list(interacting_res_map.values()),
        key=lambda x: x["min_distance"]
    )

    return {
        "polar_contact_candidates": unique_polar,
        "carbon_contact_candidates": carbon_contacts[:25],
        "interacting_residues": unique_interacting_res,
        "polar_contact_count": len(unique_polar),
        "carbon_contact_count": len(carbon_contacts),
        "disclaimer": CONTACTS_DISCLAIMER
    }

def _contact_key(contact: dict) -> tuple[str, str]:
    """Stable contact identity: chain + residue label (resname+number)."""
    return str(contact.get("chain", "")).strip(), str(contact.get("residue", "")).strip()


def _format_contact_key(key: tuple[str, str]) -> str:
    chain, residue = key
    return f"{chain}:{residue}" if chain else residue


def _parse_mutation_residue_specs(residue_specs: list | None) -> list[tuple[str | None, str]]:
    """Accept both ``A:102`` and legacy bare ``102`` mutation-residue specs."""
    parsed: list[tuple[str | None, str]] = []
    for raw in residue_specs or []:
        spec = str(raw).strip()
        if not spec:
            continue
        if ":" in spec:
            chain, residue = (part.strip() for part in spec.split(":", 1))
            parsed.append((chain or None, residue))
        else:
            parsed.append((None, spec))
    return parsed


def compare_interactions(
    healthy_interactions: dict,
    mutant_interactions: dict,
    mutation_residue_nums: list = None
) -> dict:
    """
    مقارنة التماسات الاستكشافية بين السليم والمصاب لتحديد:
    - مرشحات التماس القطبي المفقودة في المصاب.
    - مرشحات التماس القطبي الجديدة في المصاب.
    - هل يدخل حمض الطفرة في نطاق التماس المباشر مع الدواء؟

    The returned keys are the canonical ``*_contact_*`` names used by the UI;
    deprecated aliases are also kept so older callers do not break.
    """
    h_polar_keys = {
        _contact_key(c) for c in healthy_interactions.get("polar_contact_candidates", [])
    }
    m_polar_keys = {
        _contact_key(c) for c in mutant_interactions.get("polar_contact_candidates", [])
    }

    lost_polar_keys = sorted(h_polar_keys - m_polar_keys)
    gained_polar_keys = sorted(m_polar_keys - h_polar_keys)
    conserved_polar_keys = sorted(h_polar_keys.intersection(m_polar_keys))
    lost_polar_residues = [_format_contact_key(k) for k in lost_polar_keys]
    gained_polar_residues = [_format_contact_key(k) for k in gained_polar_keys]
    conserved_polar_residues = [_format_contact_key(k) for k in conserved_polar_keys]

    # Mutation engagement must compare (chain, resnum), not a CHAIN:RESIDUE
    # string against a bare residue number.  The UI passes pocket specs such as
    # {"A:102"} while residue records store res_num="102", chain="A".
    mutation_targets = _parse_mutation_residue_specs(mutation_residue_nums)

    def _is_engaged(interactions: dict) -> bool:
        for residue in interactions.get("interacting_residues", []):
            res_num = str(residue.get("res_num", "")).strip()
            chain = str(residue.get("chain", "")).strip()
            for wanted_chain, wanted_residue in mutation_targets:
                if res_num == wanted_residue and (wanted_chain is None or chain == wanted_chain):
                    return True
        return False

    mutation_engaged_healthy = _is_engaged(healthy_interactions)
    mutation_engaged_mutant = _is_engaged(mutant_interactions)

    result = {
        "lost_polar_contact_residues": lost_polar_residues,
        "gained_polar_contact_residues": gained_polar_residues,
        "conserved_polar_contact_residues": conserved_polar_residues,
        "healthy_polar_contact_count": healthy_interactions.get("polar_contact_count", 0),
        "mutant_polar_contact_count": mutant_interactions.get("polar_contact_count", 0),
        "mutation_engaged_in_contacts": mutation_engaged_healthy or mutation_engaged_mutant,
        "mutation_residues_checked": [str(n).strip() for n in (mutation_residue_nums or [])],
        "disclaimer": CONTACTS_DISCLAIMER,
    }
    # Backwards-compatible aliases (the UI used these names before the fix).
    result.update({
        "lost_polar_residues": lost_polar_residues,
        "gained_polar_residues": gained_polar_residues,
        "conserved_polar_residues": conserved_polar_residues,
        "healthy_polar_count": result["healthy_polar_contact_count"],
        "mutant_polar_count": result["mutant_polar_contact_count"],
    })
    return result

class ScoreDifferenceInterpretation(str):
    """Exploratory score difference interpretation string with helper accessors."""

    @property
    def summary(self) -> str:
        return str(self)

    @property
    def direction(self) -> str:
        if "less favorable" in str(self):
            return "less favorable"
        elif "more favorable" in str(self):
            return "more favorable"
        return "unchanged"

    @property
    def category(self) -> str:
        if self.direction == "less favorable":
            return "فارق درجة أقل ملاءمة للمصاب (Less Favorable)"
        elif self.direction == "more favorable":
            return "فارق درجة أفضل للمصاب (More Favorable)"
        return "فارق درجة غير متغير (Unchanged)"

    def __getitem__(self, key):
        if isinstance(key, str):
            if key == "summary":
                return str(self)
            elif key == "direction":
                return self.direction
            elif key == "category":
                return self.category
            raise KeyError(key)
        return super().__getitem__(key)


def interpret_score_difference(delta_score: float) -> ScoreDifferenceInterpretation:
    sign = "less favorable" if delta_score > 0 else "more favorable" if delta_score < 0 else "unchanged"
    msg = (
        f"Exploratory protocol-specific score difference: {delta_score:+.2f} kcal/mol "
        f"({sign} for the mutant under this docking setup). "
        "This is not binding affinity, resistance, efficacy, or clinical guidance."
    )
    return ScoreDifferenceInterpretation(msg)

def generate_2d_depiction_svg(smiles: str, width: int = 340, height: int = 220) -> str:
    """
    Generate an SVG string for 2D depiction of a molecule from SMILES string using RDKit.
    Returns empty string if invalid SMILES.
    """
    clean_s = smiles.strip()
    if not clean_s:
        return ""
    try:
        from rdkit import Chem
        from rdkit.Chem.Draw import rdMolDraw2D
        mol = Chem.MolFromSmiles(clean_s)
        if not mol:
            return ""
        drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
        opts = drawer.drawOptions()
        opts.clearBackground = True
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        return drawer.GetDrawingText()
    except Exception:
        return ""

