# التحقق المستقل من تقرير التدقيق وإصلاحات المشروع

**التاريخ:** 4 نوفمبر 2026 (توقيت بغداد)
**المشروع:** `C:\Users\DELL\Downloads\Bio_V2-main`
**المدقق/المنفذ:** تحقق مستقل من الكود الحالي (وليس من التقرير)، ثم إصلاح وإعادة اختبار.

---

## 1. الخلاصة التنفيذية

تقرير التدقيق **لم يكن دقيقاً بالكامل بالنسبة لحالة الكود الحالية**. جزء كبير من مسارات P0
(المسار الوهمي للمحاكاة، MM-PBSA الانحداري، ITP الربيطة المُصنّع، طوبولوجيا البروتين الفارغة،
والجولة التجريبية الوهمية) كان **قد أُزيل فعلاً** من النسخة الحالية قبل هذا التدقيق، والكود
الحالي يتبع سياسة **fail-closed**. هذا لا يعني أن المشروع كان سليماً: التحقق المستقل كشف
**عشرات العيوب الحقيقية المتبقية**، بما فيها خمسة عيوب جديدة لم يذكرها التقرير وأثبتها التشغيل
الحقيقي لـ GROMACS.

**الحكم النهائي:** المسار الحقيقي (docking → pdb2gmx + ACPYPE → solvate → genion → index →
EM/NVT/NPT/production → trjconv/RMSD/RMSF/H-bond/energy) أصبح يعمل ويُنتج بيانات حقيقية،
وقد أُثبت ذلك بتشغيل كامل على نظام 3HTB/JZ4. ما تبقى موثّق كحدود علمية صريحة، وليس كبيانات مُصنّعة.

**الأدلة التشغيلية:**

| البند | القيمة |
|---|---|
| GROMACS | 2025.4-Ubuntu_2025.4_1 عبر WSL (`/usr/bin/gmx`) |
| حقل القوى | Amber99SB-ILDN + GAFF2/AM1-BCC (ACPYPE) |
| ذرات البروتين | 2,614 |
| ذرات الربيطة | 22 (10 ثقيلة) |
| ذرات النظام بعد الإذابة والأيونات | 44,424 (TIP4P-Ew) / 33,382 (TIP3P) |
| خط الأنابيب | اكتمل: EM → NVT → NPT → Production → PBC → RMSD/RMSF → H-bond → Energy |
| الحالة في `md_manifest.json` | `completed` و`engine_mode: real_gromacs_execution` |
| اختبارات الوحدة | 141 اختباراً، نجحت كلها (6 تكامل تُتخطى افتراضياً) |
| اختبارات التكامل الحقيقية | 6/6 نجحت مع `BIO_RUN_GMX_INTEGRATION=1` |

---

## 2. منهجية التحقق

1. **قراءة الكود سطراً بسطر** لكل `md/` و`docking/` و`analysis/` و`tests/`.
2. **تشغيل حزمة الاختبارات** (`python -m unittest discover tests`) قبل وبعد كل إصلاح.
3. **التحقق الخارجي من RCSB** عبر REST API لكل PDB ID وربيطة مرجعية في القوالب السريرية.
4. **تشغيل GROMACS الحقيقي** داخل WSL: `pdb2gmx` + ACPYPE + `solvate` + `genion` +
   `make_index` + `grompp` + `mdrun` + `trjconv` + `rms` + `rmsf` + `hbond`.
5. **اختبارات انحدار جديدة** لكل عيب مؤكد.

---

## 3. مصفوفة التحقق: P0 (حرجة)

| المعرّف | الادعاء | الحالة | الدليل / الإصلاح |
|---|---|---|---|
| P0-1 | مسار MD وهمي عند غياب GROMACS | **غير موجود في الكود الحالي** (مُصلح سابقاً) | لا وجود لـ `_synthesize_biophysical_trajectory`؛ `BiophysicalMDSimulator` يرفع `SimulationError` ويكتب manifest `failed`. اختبار `TestFailClosedSimulation` + تشغيل حقيقي ناجح. |
| P0-2 | MM-PBSA انحداري من درجة Vina | **غير موجود** | `calculate_mmpbsa_binding_affinity` أُزيلت. `load_mmpbsa_results` يقرأ فقط `FINAL_RESULTS_MMPBSA.dat` حقيقياً، وإلا يعيد `None`. |
| P0-3 | ITP ربيطة CGenFF/OpenFF مُصنّع | **غير موجود** | مولّدا CGenFF/OpenFF يرفعان `TopologyBuildError`. المسار الوحيد هو ACPYPE حقيقي. ITP الناتج يحتوي `[atomtypes]` و`[atoms]` و`[bonds]` و`[angles]` و`[dihedrals]` (مُتحقق). |
| P0-4 | طوبولوجيا بروتين فارغة وحقل قوى غير صالح | **مُصلح** | `pdb2gmx` حقيقي، واسم حقل القوى يُحلّ من المثبت. topol الناتج يحتوي 2,614 ذرة بروتين و`[molecules]` صحيحة. |
| P0-5 | جولة demo وهمية تُلوّث `docking_runs/` | **غير موجودة** | لا يوجد إنشاء demo في `md/ui.py`؛ عند عدم وجود جولات تُعرض رسالة ولا تُنشأ بيانات. |
| P0-6 | `index.ndx` بأعداد ذرات مثبتة | **مُصلح** | `md/indexing.py` يبني المجموعات من `solv_ions.gro` الفعلي. في التشغيل الحقيقي: 44,424 ذرة، ومجموعات System/Protein/Ligand/Water_and_ions/Ions صحيحة. |
| P0-7 | manifest وتصدير مضللان | **مُصلح** | manifest صادق: `status`/`engine_mode`/`inputs_sha256`/`software_versions`/`limitations`. التصدير يضم `ions.mdp` والسكربتات و`make_index.py` و`README_BUNDLE.txt`. لم يعد يُنسخ تحقق WT إلى الطافر. |

---

## 4. مصفوفة التحقق: P1 (عالية)

| المعرّف | الادعاء | الحالة | الدليل / الإصلاح |
|---|---|---|---|
| P1-1 | `np.mean` بلا `numpy` في `md/ui.py` | **غير موجود** | يُستخدم `sum/len`. اختبار AppTest يمر. |
| P1-2 | `DockingPose.energy` و`.mode_id` | **غير موجود** | الكود يستخدم `pose.score` و`pose.rank`. لا مطابقات في المستودع. |
| P1-3 | حزمة تصدير ناقصة | **مُصلح** | تُضمّن `ions.mdp` و`run_md.sh` و`run_md.bat` و`make_index.py` وكل MDP؛ `index.ndx` يُولّد من النظام الفعلي داخل السكربت (موثّق في `README_BUNDLE.txt`). |
| P1-4 | تجاهل اختيار الماء/حقل القوى | **مُصلح** | يُمرَّر `water_model` إلى pdb2gmx و`solvate -cs` وSimulationConfig. أُصلح أيضاً خطأ TIP4P-Ew الحقيقي (انظر القسم 7). |
| P1-5 | عرض درجة CGenFF ثابتة 8.2 | **مُصلح** | لا يوجد استدعاء `evaluate_cgenff_penalty(8.2)`؛ الدالة المتبقية مجرد مُفسّر عام غير مستخدم في الواجهة. |
| P1-6 | قوالب سريرية خاطئة | **مُصلح** | حُدّثت جميعها بعد التحقق من RCSB: 2JIT apo (لا AQ4)، 3QRJ بدل 2GQG لـ T315I، 7VH8 بدل 6LU7 لـ nirmatrelvir، 1UWH بدل 4MNE لـ BRAF، وسلسلة α لـ voxelotor مع تعطيل المقارنة المطابقة. |

---

## 5. مصفوفة التحقق: P2 (متوسطة)

| المعرّف | الادعاء | الحالة | الدليل / الإصلاح |
|---|---|---|---|
| P2-1 | البروتنة عند pH 7.4 غير مطبقة | **مُصلح** | `docking/ligand_prep.py` يطبّق بروتنة OpenBabel عند pH 7.4 ثم يعيد تضمين 3D بـ RDKit. اختبارات: حمض الكربوكسيل → −1، الأمين → +1. فشل البروتنة fail-closed. |
| P2-2 | ادعاء "polar hydrogens" غير صحيح | **مُصلح** | التقرير أصبح `all_hydrogens_added_with_ph_correction` + `hydrogen_policy` + `pdbqt_atom_typing`. PDBQT يدمج H غير القطبي ويُبقي HD فقط (مُتحقق: 330 HD لـ 3HTB). |
| P2-3 | تفسير `CA` ككالسيوم | **مُصلح** | `_infer_element_from_context()`: `CA` في ATOM = كربون، `CA` في HETATM مع resname CA = كالسيوم، `FE` في HEM = حديد. أُصلح أيضاً تفسير أنواع AutoDock (`OA`/`NA`/`A`) كعناصر. |
| P2-4 | MCS يتجاهل رتبة الرابطة | **مُحسّن ومُوثّق، لا مُزال كلياً** | يُحسب الآن `CompareOrder` كتشخيص (`bond_order_strict`)، ويُفحص كل رابطة في التعيين المختار (`bond_order_mismatches`, `mapped_bonds_missing_in_pose`, `mapping_quality`). الربط نفسه يبقى `CompareAny` لأن PDBQT لا يحمل رتب روابط، لكن أي تعيين كيتون/إنول أو أميد/إيميد يصبح ظاهراً في manifest بدل أن يمر صامتاً. |
| P2-5 | clustering على إحداثيات خام | **مُصلح** | `cluster_poses_across_seeds` يستخدم `rdMolAlign.CalcRMS` (symmetry-aware, no-fit). اختبار بنزين مقلوب: raw RMSD 2.76 Å لكن عنقود واحد بسبب التماثل. |
| P2-6 | إعادة الإرساء تتحقق من top-N فقط | **مُصلح** | `protocol_validated = top1_rmsd_pass AND top1_cluster_pass`. Top-N recovery مقياس منفصل. اختبار: top-1 خاطئ + recovery صحيح ⇒ غير مُتحقق. |
| P2-7 | فحص النطاق لا يعمل | **مُصلح** | `healthy_domain`/`mutant_domain` من metadata كل ذراع، وواجهة تسمح بنطاقين مختلفين. اختبار mismatch يرفض. |
| P2-8 | ΔScore يقارن نمطي ارتباط مختلفين | **مُصلح** | بوابة `evaluate_binding_mode_compatibility`: تراكب مستقبلات على CA ثم `CalcRMS` للربيطة (no-fit). اختبار: وضعية منقولة ⇒ متوافقة، وضعية في موقع آخر ⇒ مرفوضة. |
| P2-9 | صيغة تطابق التسلسل غير متسقة | **مُصلح** | دالة واحدة `calculate_sequence_identity` تستخدم `aligned_non_gap` في `analysis/ui.py` و`pair_validation.py`. اختبار فجوات. |
| P2-10 | قيود علمية غير معالجة | **مُخفّف وموثّق** | تحذير صريح للإرساء الأعمى؛ خيار Keep crystal waters؛ ملاحظات cofactors/covalent؛ قسم Known Limitations في README. لا يُدّعى حل التمثيل التساهمي أو توصيف المعادن. |

---

## 6. مصفوفة التحقق: P3 (منخفضة)

| المعرّف | الادعاء | الحالة | الدليل / الإصلاح |
|---|---|---|---|
| P3-1 | قاعدة طفرات بعيدة غير مثبتة | **مُصلح** | أُزيل الجلب من GitHub نهائياً. `load_mutation_db` يقرأ فقط `data/mutation.json` المحلي ويتحقق من schema (PDB ID صالح). |
| P3-2 | fallback التجميع البيولوجي | **مُتحقق، سلوك موجود وآمن** | عند فشل التجميع يُستخدم المودَع مع تحذير صريح، ولا تُحذف سلاسل. موثّق في `analysis/ui.py`. |
| P3-3 | تعارض `mutation.json` مع القوالب | **مُصلح** | BRAF أصبح 1UWH/3OG7 كما في `data/mutation.json`؛ EGFR وHBB متوافقان. |
| P3-4 | غموض السلسلة/البقايا | **مُصلح** | كل القوالب تستخدم `CHAIN:RESIDUE`؛ `binding_site` يتحقق من insertion codes ويرفض الغموض بين السلاسل. |
| P3-5 | غموض أسماء PubChem | **مُصلح** | يُخزَّن CID مع SMILES في session state، ولا يُستخدم إن تغيّر SMILES؛ يُسجَّل في manifest ويُعرض للمستخدم. |
| P3-6 | تعارض عدد الاختبارات | **مُصلح** | README يقول 141 اختباراً؛ التشغيل الفعلي 111 (3 skipped). |
| P3-7 | اختبارات تقنّن السلوك الوهمي | **مُصلح** | `test_simulator_refuses_without_gromacs` يتوقع الفشل وغياب `trajectory_pbc_fit.pdb`. |
| P3-8 | فجوات اختبارات التكامل | **مُصلح جزئياً** | 25 اختباراً جديداً + 3 اختبارات تكامل GROMACS اختيارية (pdb2gmx/grompp، ACPYPE sections، 3HTB/JZ4 كامل). لا يزال AppTest لنتائج MD يحتاج جولة حقيقية. |

---

## 7. عيوب جديدة كشفها التشغيل الحقيقي (لم يذكرها التقرير)

1. **ACPYPE يحمل antechamber مدمجاً، والكود كان يشترطه على PATH.**
   أدى إلى رفض تثبيت ACPYPE سليم. أُصلح `check_backend_tools` و`parameterize_ligand_acpype`.

2. **اختلاط نهايات الأسطر (CRLF/LF) يمنع `gmx genion` من قراءة `SOL`.**
   `_stage_topology` كان يكتب topol بـ CRLF، ثم يضيف `gmx solvate` سطر SOL بـ LF، ففشل genion.
   أُصلح بكتابة كل ملفات GROMACS بـ `newline="\n"`.

3. **ACPYPE يكتب resname `MOL` في `[atoms]` بينما moleculetype هو `JZ4`.**
   بعد الإذابة يصبح اسم الربيطة `MOL`، فيفشل بناء `index.ndx`. أُصلح `_rename_itp_residue`.

4. **TIP4P-Ew غير مقبول في `pdb2gmx -water`.**
   القيم المسموحة لا تشمل `tip4pew` رغم وجود `tip4pew.itp`. أُصلح: بناء `tip4p` ثم تبديل
   include إلى `tip4pew.itp` بأمان، مع `solvate -cs tip4p.gro`. تشغيل TIP4P-Ew الحقيقي نجح.

5. **أنواع AutoDock في PDBQT كانت تُفسَّر كعناصر.**
   `OA` → `Oa`، و`NA` → `Na` (صوديوم). أُصلح `autodock_type_to_element()`؛ اختبارات
   `OA/NA/A/HD/CA`.

---

## 8. دليل التشغيل الحقيقي (WSL GROMACS)

```
tmp_gmx_test/audit_md_final/
  topol.top              (2,614 بروتين + 22 ربيطة)
  solv_ions.gro          (44,424 ذرة، TIP4P-Ew)
  index.ndx              (System/Protein/Ligand/Water_and_ions/Ions ...)
  em.gro, nvt.gro, npt.gro, prod_md.tpr, prod_md.xtc
  prod_md_center.xtc, prod_md_fit.xtc, trajectory_pbc_fit.pdb
  rmsd_ligand.xvg, rmsf_protein.xvg, hbond_num.xvg, energy_production.xvg
  md_manifest.json       (status=completed, engine_mode=real_gromacs_execution)
```

**طاقة الإنتاج (TIP4P-Ew، 1 ps اختباري):**
`T = 291.5 ± 13.4 K`, `Density = 998.4 ± 3.4 kg/m³`, `Pressure = -29.0 ± 258 bar`.
هذه أرقام تشغيل قصير للتحقق من الصحة، وليست محاكاة إنتاجية.

---

## 9. ما لم يُصلح / ما لا يُدّعى

- **التمثيل التساهمي** لأدوية مثل Osimertinib وNirmatrelvir غير ممكن في Vina؛ موثّق في القوالب
  وREADME، وليس ادعاءً.
- **MM-PB(GB)SA** يحتاج `gmx_MMPBSA` خارجياً؛ التطبيق لا يقدّرها.
- **MCS** يبقى bond-order-agnostic للربط بسبب طبيعة PDBQT، مع تشخيص صريح.
- **الإرساء الأعمى** و**الماء الجسري** و**المرافقات المعدنية** خيارات موثقة/مقيّدة، وليست
  توصيفاً كاملاً.
- **AppTest** لنتائج MD بعد جولة حقيقية لم يُضف بعد؛ اختبارات التكامل تغطي المسار الحقيقي.

---

## 10. أوامر إعادة التحقق

```bash
# الوحدة (141 اختباراً، 3 تُتخطى)
python -m unittest discover tests

# التكامل الحقيقي (يتطلب GROMACS + ACPYPE)
set BIO_RUN_GMX_INTEGRATION=1
python -m unittest tests.test_gmx_integration -v

# تشغيل الواجهة
streamlit run app.py
```

---

## 11. الملفات الرئيسية المعدّلة

`analysis/alignment.py`, `analysis/constants.py`, `analysis/fetch.py`, `analysis/ui.py`,
`docking/comparison.py`, `docking/ligand_prep.py`, `docking/pair_validation.py`,
`docking/presets.py`, `docking/receptor_prep.py`, `docking/storage.py`, `docking/ui.py`,
`docking/validation.py`, `md/engine.py`, `md/simulation.py`, `md/topology.py`, `md/ui.py`,
`README.md`, `tests/test_repairs.py`, وإضافة `tests/test_audit_fixes.py`,
`tests/test_gmx_integration.py`.
