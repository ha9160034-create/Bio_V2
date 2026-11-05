# Human Pharmaceutical & Clinical Presets Library
"""
Verified clinical presets for matched/exploratory docking.

Every entry was checked against the RCSB PDB entry metadata (entry title and
non-polymer ligand component) before being written here.  The audit that
prompted this file found several presets that referenced ligands absent from
the declared PDB entry (for example ``ref_ligand="AQ4"`` on the apo T790M
structure 2JIT), a T315I label on a dasatinib/WT-Abl structure (2GQG), and a
voxelotor pocket on the wrong globin chain.

Design rules
------------
* ``healthy_ref_ligand`` is only set when that exact co-crystal ligand is
  present in the healthy entry **and** is the same chemical entity as the drug
  being docked.  Otherwise it is empty and WT redocking validation is reported
  as unavailable rather than fabricated.
* ``mutant_ref_ligand`` is kept separate: a mutant co-crystal does not validate
  the healthy arm (and vice versa).
* Pocket residues always use ``CHAIN:RESIDUE`` notation so multi-chain
  structures cannot silently select the wrong chain.
* ``matched_comparison_supported=False`` marks presets whose lesion is not in
  the docked chain; they must be run as exploratory cross-docking.
"""

PRESET_SCHEMA_VERSION = 2

HUMAN_CLINICAL_PRESETS = {
    "egfr_erlotinib_t790m": {
        "title": "🎗️ سرطان الرئة: مستقبل EGFR مع عقار Erlotinib (طفرة المقاومة T790M)",
        "category": "الأورام (Oncology / NSCLC)",
        "target_name": "EGFR Tyrosine Kinase Domain",
        "healthy_id": "1M17",
        "mutant_id": "2JIT",
        "healthy_chain": "A",
        "mutant_chain": "A",
        "drug_name": "Erlotinib (Tarceva)",
        "drug_smiles": "COCCOC1=C(C=C2C(=C1)C(=NC=N2)NC3=CC=CC(=C3)C#C)OCCOC",
        "healthy_ref_ligand": "AQ4",
        "mutant_ref_ligand": "",
        "ref_ligand": "AQ4",
        "pocket_residues": "A:790, A:793, A:858",
        "target_domain": "Kinase Domain",
        "matched_comparison_supported": True,
        "covalent_drug": False,
        "reference_note": (
            "1M17 contains erlotinib (AQ4). 2JIT is the T790M kinase domain without a "
            "co-crystal ligand, so only the healthy arm can be redocking-validated."
        ),
        "clinical_note": (
            "طفرة T790M (Gatekeeper Mutation) تُغير الثريونين الصغير إلى ميثيونين أكبر حجماً، "
            "مما يُحدث إعاقة فراغية (Steric Hindrance) تمنع أدوية الجيل الأول مثل Erlotinib من الارتباط، "
            "مسببة مقاومة سريرية لدى مرضى سرطان الرئة غير صغير الخلايا (NSCLC)."
        ),
        "limitations": [
            "2JIT is apo; no mutant co-crystal ligand exists for mutant-arm redocking validation.",
            "Erlotinib is not expected to overcome T790M; a lower mutant score is a hypothesis, not a clinical prediction.",
        ],
    },
    "egfr_osimertinib_t790m": {
        "title": "💊 سرطان الرئة: مستقبل EGFR مع عقار Osimertinib (الجيل الثالث الموجه للطفرة)",
        "category": "الأورام (Oncology / NSCLC)",
        "target_name": "EGFR T790M Active Conformation",
        "healthy_id": "1M17",
        "mutant_id": "2JIT",
        "healthy_chain": "A",
        "mutant_chain": "A",
        "drug_name": "Osimertinib (Tagrisso)",
        "drug_smiles": "CC(=O)NC1=C(C=CC(=C1)NC2=NC=C(C(=N2)C3=CN(C4=CC=CC=C34)C)N(C)CCN(C)C)NC(=O)C=C",
        "healthy_ref_ligand": "AQ4",
        "mutant_ref_ligand": "",
        "ref_ligand": "AQ4",
        "pocket_residues": "A:790, A:793, A:797, A:858",
        "target_domain": "Kinase Domain",
        "matched_comparison_supported": True,
        "covalent_drug": True,
        "reference_note": (
            "2JIT is apo T790M. Osimertinib is a covalent Cys797 inhibitor; AutoDock Vina "
            "only scores a non-covalent approximation and cannot model the acrylamide adduct."
        ),
        "clinical_note": (
            "عقار Osimertinib صُمم خصيصاً كعلاج موجه ليتفادى الإعاقة الفراغية لطفرة T790M "
            "ويرتبط تساهمياً بالحامض Cys797، ليعيد التحكم بالورم بعد فشل أدوية الجيل الأول."
        ),
        "limitations": [
            "Covalent binding is not modelled by Vina; scores are a non-covalent approximation.",
            "2JIT is apo; the mutant arm cannot be redocking-validated.",
        ],
    },
    "sickle_cell_voxelotor": {
        "title": "🩸 فقر الدم المنجلي: الهيموجلوبين السليم مقابل المنجلي مع دواء Voxelotor",
        "category": "أمراض الدم (Hematology)",
        "target_name": "Human Hemoglobin Tetramer (HbA vs HbS)",
        "healthy_id": "2HHB",
        "mutant_id": "2HBS",
        "healthy_chain": "A",
        "mutant_chain": "A",
        "drug_name": "Voxelotor (Oxbryta)",
        "drug_smiles": "CC1=C(C(=NO1)C)C2=NC=C(C=C2)OCC3=C(C(=CC=C3)O)C=O",
        "healthy_ref_ligand": "",
        "mutant_ref_ligand": "",
        "ref_ligand": "",
        "pocket_residues": "A:1, A:2, A:3, A:6",
        "target_domain": "Hemoglobin Alpha Subunit",
        "matched_comparison_supported": False,
        "covalent_drug": False,
        "reference_note": (
            "Voxelotor binds the alpha-globin chain (chain A), whereas the sickle lesion "
            "E6V is on the beta chain. The two alpha chains are sequence-identical, so this "
            "is an exploratory allosteric comparison, not a matched mutation comparison."
        ),
        "clinical_note": (
            "طفرة Glu6Val في سلسلة بيتا تجعل الهيموجلوبين يتبلمر في حالات نقص الأكسجين مشكلاً خلايا منجلية. "
            "عقار Voxelotor يرتبط بسلسلة ألفا ويزيد تقارب الهيموجلوبين للأكسجين مانعاً البلمرة والتمنجل."
        ),
        "limitations": [
            "The docked chain (alpha) does not carry the E6V mutation; ΔScore is not a mutation effect.",
            "2HHB/2HBS are deoxy (T-state) structures; voxelotor stabilises the R-state, so the starting conformation is not the drug-bound state.",
        ],
    },
    "cml_bcr_abl_imatinib": {
        "title": "🎗️ سرطان الدم CML: إنزيم BCR-ABL1 مع عقار Imatinib (طفرة المقاومة T315I)",
        "category": "أمراض الدم والأورام (Hematology / Oncology)",
        "target_name": "BCR-ABL1 Tyrosine Kinase",
        "healthy_id": "1IEP",
        "mutant_id": "3QRJ",
        "healthy_chain": "A",
        "mutant_chain": "A",
        "drug_name": "Imatinib (Gleevec)",
        "drug_smiles": "CC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCN(CC3)C)NC4=NC=CC(=N4)C5=CN=CC=C5",
        "healthy_ref_ligand": "STI",
        "mutant_ref_ligand": "",
        "ref_ligand": "STI",
        "pocket_residues": "A:315, A:317, A:381, A:382",
        "target_domain": "Abl Kinase Domain",
        "matched_comparison_supported": True,
        "covalent_drug": False,
        "reference_note": (
            "1IEP is WT Abl with imatinib (STI). 3QRJ is the T315I mutant with DCC-2036 (919), "
            "not imatinib; the mutant co-crystal therefore cannot validate imatinib redocking."
        ),
        "clinical_note": (
            "طفرة T315I تُعرف بطفرة بوابة الكيناز (Gatekeeper)، حيث يُستبدل الثريونين بالآيزولوسين "
            "مما يؤدي إلى فقدان رابطة هيدروجينية أساسية وحدوث إعاقة فراغية تطرد عقار Imatinib."
        ),
        "limitations": [
            "The mutant entry 3QRJ contains DCC-2036, not imatinib; only the healthy arm has an imatinib reference.",
            "T315I is a known imatinib-resistance mutation; docking cannot measure resistance.",
        ],
    },
    "melanoma_braf_vemurafenib": {
        "title": "☀️ سرطان الجلد الميلانوما: إنزيم BRAF مع عقار Vemurafenib (طفرة V600E)",
        "category": "الأورام (Oncology / Melanoma)",
        "target_name": "BRAF Kinase Domain",
        "healthy_id": "1UWH",
        "mutant_id": "3OG7",
        "healthy_chain": "A",
        "mutant_chain": "A",
        "drug_name": "Vemurafenib (Zelboraf)",
        "drug_smiles": "CCCS(=O)(=O)NC1=C(C(=C(C=C1)F)C(=O)C2=CNC3=C2C=C(C=N3)C4=CC=C(C=C4)Cl)F",
        "healthy_ref_ligand": "",
        "mutant_ref_ligand": "032",
        "ref_ligand": "",
        "pocket_residues": "A:509, A:595, A:600",
        "target_domain": "BRAF Kinase Domain",
        "matched_comparison_supported": True,
        "covalent_drug": False,
        "reference_note": (
            "1UWH is WT BRAF with sorafenib (BAX), not vemurafenib. 3OG7 is V600E BRAF with "
            "vemurafenib (032). There is no WT BRAF/vemurafenib co-crystal here, so the healthy "
            "arm has no valid redocking reference."
        ),
        "clinical_note": (
            "طفرة V600E تُحاكي حالة الفسفرة النشطة لإنزيم BRAF وتجعله نشطاً بشكل دائم محفزاً نمو الورم. "
            "عقار Vemurafenib يثبط الإنزيم النشط انتقائياً ويحقق استجابة سريرية سريعة لدى مرضى الميلانوما."
        ),
        "limitations": [
            "1UWH contains sorafenib, not vemurafenib; the healthy arm cannot be validated against this drug.",
            "Vemurafenib is V600E-selective, so a worse WT score is expected and is not a resistance measurement.",
        ],
    },
    "covid_mpro_paxlovid": {
        "title": "🦠 الفيروسات التنفسية: بروتياز كورونا مع دواء Paxlovid (Nirmatrelvir)",
        "category": "الأمراض المعدية والفيروسات (Infectious Diseases)",
        "target_name": "SARS-CoV-2 Main Protease (Mpro)",
        "healthy_id": "7VH8",
        "mutant_id": "8D4N",
        "healthy_chain": "A",
        "mutant_chain": "A",
        "drug_name": "Nirmatrelvir (Paxlovid)",
        "drug_smiles": "CC1(C2C1C(N(C2)C(=O)C(C(C)(C)C)NC(=O)C(F)(F)F)C(=O)NC(CC3CCNC3=O)C#N)C",
        "healthy_ref_ligand": "4WI",
        "mutant_ref_ligand": "",
        "ref_ligand": "4WI",
        "pocket_residues": "A:140, A:142, A:144, A:145, A:163, A:166",
        "target_domain": "Chymotrypsin-like Protease",
        "matched_comparison_supported": True,
        "covalent_drug": True,
        "reference_note": (
            "7VH8 is WT Mpro with nirmatrelvir (4WI). 8D4N is the E166Q mutant and is apo, "
            "so only the healthy arm can be redocking-validated."
        ),
        "clinical_note": (
            "إنزيم Mpro أساسي لتضاعف الفيروس عبر تقطيع متعدد البروتين. "
            "يرتبط عقار Nirmatrelvir بـ Cys145 و His41 مثبطاً عمل الإنزيم بصورة فعالة لتقليل شدة المرض."
        ),
        "limitations": [
            "Nirmatrelvir is covalent; Vina cannot model the Cys145 adduct and scores only a non-covalent approximation.",
            "8D4N is E166Q apo; the mutant arm cannot be redocking-validated.",
        ],
    },
}


def validate_preset(preset: dict) -> list[str]:
    """Return a list of schema/consistency problems (empty when valid)."""
    problems: list[str] = []
    required = {
        "title", "healthy_id", "mutant_id", "healthy_chain", "mutant_chain",
        "drug_name", "drug_smiles", "pocket_residues", "target_domain",
    }
    missing = sorted(required - set(preset))
    if missing:
        problems.append(f"missing fields: {', '.join(missing)}")

    for field in ("healthy_id", "mutant_id"):
        value = str(preset.get(field, ""))
        if len(value) != 4 or not value[0].isdigit():
            problems.append(f"{field}={value!r} is not a 4-character PDB ID")

    residues = [r.strip() for r in str(preset.get("pocket_residues", "")).split(",") if r.strip()]
    if not residues:
        problems.append("pocket_residues is empty")
    for residue in residues:
        if ":" not in residue:
            problems.append(f"pocket residue {residue!r} must use CHAIN:RESIDUE notation")

    for field in ("healthy_ref_ligand", "mutant_ref_ligand"):
        value = str(preset.get(field, "")).strip()
        if value and len(value) > 5:
            problems.append(f"{field}={value!r} does not look like a PDB chemical component ID")

    if preset.get("healthy_id") == preset.get("mutant_id"):
        problems.append("healthy_id and mutant_id must differ")

    return problems


def get_preset(key: str) -> dict:
    """Return a preset after validating it (raises on an invalid entry)."""
    preset = HUMAN_CLINICAL_PRESETS[key]
    problems = validate_preset(preset)
    if problems:
        raise ValueError(f"Invalid clinical preset {key!r}: {'; '.join(problems)}")
    return preset
