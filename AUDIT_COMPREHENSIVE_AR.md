# تقرير الفحص الشامل لمشروع Bio_V2-main

**تاريخ الفحص:** 4 نوفمبر 2026
**المشروع:** `C:\Users\DELL\Downloads\Bio_V2-main`
**المرجعيات:** `AGENTS.md`، صفحات `wiki/` العلمية، الكود الحالي، الاختبارات، وتوثيق `gmx_MMPBSA` الرسمي.
**نوع الجولة:** فحص وتشخيص — **لم أجرِ تعديلات جديدة على الكود في هذه الجولة**.

---

## 1. الحكم العام

المشروع **يعمل فعلياً** في المسار الحقيقي: `pdb2gmx + ACPYPE + GROMACS`، ولا يحتوي على مسار المحاكاة الوهمي أو MM-PBSA الانحداري الذي وصفه التقرير القديم. الاختبارات الحالية **111 اختباراً تنجح** (6 تكامل تُتخطى افتراضياً)، واختبارات التكامل الأربعة تنجح مع WSL GROMACS.

لكن الفحص المستقل كشف **مشكلات حقيقية متبقية**، أهمها:

1. **لوحة مقارنة التفاعلات في الإرساء تعرض أصفاراً دائماً** بسبب اختلاف أسماء المفاتيح بين `compare_interactions` و`docking/ui.py`.
2. **قارئ نتائج MM-PBSA لا يقرأ ملف `gmx_MMPBSA` الحقيقي** (يفشل مع صيغة `ΔVDWAALS` و`ΔG binding`)، وأمر MM-PBSA المعروض في الواجهة **غير صالح** وفق توثيق `gmx_MMPBSA`.
3. **RMSD للربيطة يُحسب نسبةً إلى بداية الإنتاج (`prod_md.tpr`) وليس إلى البنية المصغّرة بالطاقة (`em.tpr`)**، خلافاً لبروتوكول Lemkul في `wiki/gromacs-protein-ligand-protocol.md`.
4. **المحاكاة مسار واحد فقط**، بلا replicates أو NVE drift check أو burn-in أو أشرطة خطأ، بينما `wiki/simulation-and-benchmark-best-practices.md` يشترط 3–5 مسارات مستقلة.
5. **MD يرمي كل HETATM** (بما فيها MSE الحقيقي والمرافقات والماء الجسري) قبل `pdb2gmx`، رغم أن `wiki/receptor-preparation.md` يشترط الحفاظ على الماء الجسري والمعادن.
6. **`prepare_receptor` لا يحل altloc ولا يحدد model 0**، ما قد يُدخل ذرتين متزامنتين في PDBQT.
7. **الـ backend الأصلي على Windows معطّل عملياً** لأن اكتشاف حقول القوى والأدوات يستخدم bash، و`shell_command` يرفض win32.
8. **تحليل عناصر PDBQT في `docking/interactions.py`** يقرأ نوع AutoDock كأنه عنصر؛ الكربون العطري `A` لا يُحتسب ككربون، والهالوجينات قد تُقرأ خطأً.
9. **`run_matched_docking_comparison` لا يفرض `pair_validation_meta` داخلياً**؛ الحماية موجودة في الواجهة فقط، ويمكن استدعاء الدالة برمجياً على بروتينين غير متجانسين.
10. **قوالب `docking/presets.py` غير موصولة بالواجهة**؛ صُححت بياناتها لكن لا يوجد مسار يستخدمها.

---

## 2. المنهجية

1. قراءة `AGENTS.md` بالكامل، وربط كل قاعدة بالكود.
2. قراءة صفحات wiki الأساسية:
   - `wiki/gromacs-protein-ligand-protocol.md`
   - `wiki/simulation-and-benchmark-best-practices.md`
   - `wiki/mm-pbsa-and-mm-gbsa.md`
   - `wiki/force-fields.md`
   - `wiki/open-force-field-and-ligand-parameterization.md`
   - `wiki/autodock-vina.md`
   - `wiki/docking-failure-modes.md`
   - `wiki/ligand-preparation-and-protonation.md`
   - `wiki/receptor-preparation.md`
   - `wiki/binding-site-detection.md`
   - `wiki/protein-ligand-interactions.md`
3. قراءة كل وحدات `md/` و`docking/` و`analysis/` و`tests/`.
4. تشغيل `python -m unittest discover tests` → **111 اختباراً، OK (skipped=3)**.
5. تشغيل `BIO_RUN_GMX_INTEGRATION=1 python -m unittest tests.test_gmx_integration` → **4/4 OK**.
6. اختبارات إعادة إنتاج صغيرة للمشكلات المشتبه بها (مفاتيح التفاعلات، قارئ MM-PBSA، MSE، altloc، Windows native).
7. الرجوع إلى توثيق `gmx_MMPBSA` الرسمي للتحقق من صحة الأمر:
   <https://valdes-tresanco-ms.github.io/gmx_MMPBSA/1.7.0/gmx_MMPBSA_command-line/>

---

## 3. ما هو سليم ويستحق الحفاظ عليه

- **Fail-closed MD**: لا مسار وهمي؛ غياب GROMACS أو فشل أي خطوة يكتب `md_manifest.json` بحالة `failed`.
- **طوبولوجيا حقيقية**: `pdb2gmx` + ACPYPE؛ ITP يحتوي `[atomtypes]`, `[atoms]`, `[bonds]`, `[angles]`, `[dihedrals]`.
- **Index من النظام الفعلي**: `md/indexing.py` يبني المجموعات من `solv_ions.gro`؛ لا أعداد ذرات مثبتة.
- **Guardrails**: صندوق dodecahedron، `v-rescale`، `Parrinello-Rahman` في الإنتاج، `tc-grps` مدمج، قيود موضعية على ذرات الربيطة الثقيلة، PBC rewrap.
- **تحقق إرساء صارم نسبياً**: top-1 + دعم بذرتين لحالة `protocol validated`، وتجميع symmetry-aware، وبوابة تطابق نمط الارتباط قبل ΔScore.
- **بروتنة pH 7.4**: OpenBabel pH-aware ثم 3D RDKit؛ فشلها fail-closed.
- **شفافية**: SHA-256 للمدخلات، إصدارات GROMACS/RDKit/OpenBabel، manifest صادق.
- **قوالب سريرية مصححة**: 2JIT apo، 3QRJ لـ T315I، 7VH8 لـ nirmatrelvir، 1UWH لـ BRAF، سلسلة α لـ voxelotor (مع تعطيل matched).
- **اختبارات**: 111 اختباراً، منها اختبارات انحدار وتكامل حقيقي.

---

## 4. المشكلات عالية الأولوية (P1)

### P1-1: لوحة مقارنة التفاعلات تعرض أصفاراً دائماً

**الموقع:**
- `docking/interactions.py:223-231` يُرجع:
  `healthy_polar_contact_count`, `mutant_polar_contact_count`,
  `lost_polar_contact_residues`, `gained_polar_contact_residues`,
  `conserved_polar_contact_residues`.
- `docking/ui.py:1956-1960` يقرأ:
  `healthy_polar_count`, `mutant_polar_count`,
  `lost_polar_residues`, `gained_polar_residues`, `conserved_polar_residues`.

**الأثر:** كل مقاييس التماس القطبي في تبويب المقارنة تظهر `0`، ولا تُعرض أي روابط مفقودة/مكتسبة، حتى لو كانت موجودة فعلاً.

**الدليل التجريبي:**
```text
actual key: ['GLN102']
UI key expected: None
```

**الإصلاح المقترح:** توحيد أسماء المفاتيح في `compare_interactions` أو في الواجهة، وإضافة اختبار يتحقق من القيم غير الصفرية.

---

### P1-2: قارئ MM-PBSA لا يقرأ مخرجات `gmx_MMPBSA` الحقيقية

**الموقع:** `md/analysis.py:326-388`.

**السبب:** القارئ يبحث عن:
- `Delta G binding` (نص ASCII)، بينما مخرجات `gmx_MMPBSA` تستخدم `ΔG binding`.
- مفاتيح `VDWAALS`, `EEL`, `EGB`, `ESURF`، بينما مخرجات `gmx_MMPBSA` تستخدم `ΔVDWAALS`, `ΔEEL`, `ΔEGB`, `ΔESURF`.

**الأثر:** وجود `FINAL_RESULTS_MMPBSA.dat` حقيقي لا يكفي؛ `load_mmpbsa_results` يعيد `None`، وتظل الواجهة تقول "MM-PBSA غير متاح".

**الدليل:** اختبار على عينة من صيغة gmx_MMPBSA الرسمية أعاد `None`.

**الإصلاح المقترح:** دعم البادئة `Δ` و`Delta`، وقراءة `ΔTOTAL` و`ΔG binding`، والتحقق من `Calculations performed using N complex frames`، وتسجيل `-TΔS` إن وُجد.

---

### P1-3: أمر MM-PBSA المعروض في الواجهة غير صالح

**الموقع:** `md/ui.py:559-562`.

الأمر المعروض:
```bash
gmx_MMPBSA -O -f prod_md.xtc -s prod_md.tpr -n index.ndx \
  -ci index.ndx -cg Protein Ligand -cs prod_md.tpr -o FINAL_RESULTS_MMPBSA.dat
```

وفق توثيق `gmx_MMPBSA` الرسمي، الخيارات هي:
`-cs <complex structure/tpr>`, `-ci <index>`, `-cg <receptor> <ligand>`, `-ct <trajectory>`, و`-cp <topology>` عند الحاجة. لا توجد خيارات `-f`, `-s`, `-n` في `gmx_MMPBSA`.

**الأثر:** المستخدم الذي ينسخ الأمر سيحصل على خطأ، ولن تظهر نتائج MM-PBSA.

**الإصلاح المقترح:** توليد أمر مثل:
```bash
gmx_MMPBSA -O -cs prod_md.tpr -ci index.ndx -cg Protein Ligand \
  -ct prod_md_center.xtc -cp topol.top -o FINAL_RESULTS_MMPBSA.dat
```
مع التحقق من أن `prod_md_center.xtc` مُعاد لفّه ومتوافق مع `-cp`.

---

### P1-4: RMSD الربيطة يُحسب نسبةً إلى بداية الإنتاج لا إلى بنية EM

**الموقع:** `md/simulation.py:739-744`:
```python
["rms", "-s", "prod_md.tpr", "-f", "prod_md_fit.xtc", ...]
```

**المرجع:** `wiki/gromacs-protein-ligand-protocol.md:204-211` يوضح أن RMSD الربيطة يُحسب نسبةً إلى البنية المصغّرة بالطاقة (`gmx rms -s em.tpr ...`).

**الأثر:** `rmsd_ligand.xvg` يقيس الانزياح عن بداية الإنتاج (نهاية NPT)، وليس عن الوضعية الأصلية بعد EM؛ قد يبدو الدواء "مستقراً" حتى لو انزاح خلال NPT.

**الإصلاح المقترح:** استخدام `-s em.tpr` أو أول إطار EM كمرجع، وتوثيق المرجع في manifest.

---

### P1-5: مسار إنتاج واحد بلا replicates/NVE/burn-in/أشرطة خطأ

**الموقع:** `md/simulation.py:706-713` يُشغّل `prod_md` مرة واحدة؛ `md/ui.py:517-528` يعرض متوسط RMSD ومتوسط H-bonds كأرقام نهائية.

**المرجع:**
- `wiki/simulation-and-benchmark-best-practices.md:18` يطلب NVE drift check للتحقق من dt/cutoffs.
- السطور 35-37 تطلب 3–5 مسارات مستقلة وبذر velocities مختلفة وburn-in.
- `wiki/molecular-dynamics-simulation.md:89-90` يطلب replicates والتحقق من convergence.

**الأثر:** النتائج تصلح كفحص استقرار نسبي، لكنها **لا تصلح** لادعاء تقارب عينات أو استنتاج طاقة/ديناميكا موثوقة. عرض "متوسط RMSD" بلا خطأ معياري أو burn-in مضلل.

**الإصلاح المقترح:** إضافة `n_replicates` (3–5) وبذور مختلفة، NVE قصير، تقرير block averages/CI، واستبعاد أول 10–20% من المسار.

---

### P1-6: `run_matched_docking_comparison` لا يفرض التحقق من الزوج

**الموقع:** `docking/comparison.py:375` يستقبل `pair_validation_meta`، لكنه **يخزّنه فقط** في manifest (`:445`, `:468`, `:481`) ولا يتحقق من `is_valid`. التحقق الفعلي يتم في `docking/ui.py` قبل الاستدعاء.

**الأثر:** أي استدعاء برمجي مباشر (API/سكربت) يمكنه حساب ΔScore بين بروتينين غير متجانسين، رغم أن `wiki/receptor-preparation.md` و`docking/pair_validation.py` يشترطان ≥90% تطابق.

**الإصلاح المقترح:** داخل الدالة:
```python
if not (pair_validation_meta or {}).get("is_valid"):
    raise ValueError("Matched comparison requires a validated homologous pair.")
```
وكذلك التحقق من تطابق pocket mapping.

---

### P1-7: MD يرمي MSE الحقيقي والمرافقات والماء الجسري

**الموقع:** `md/topology.py:381-428`.

`clean_receptor_for_md` يعالج `record.startswith("ATOM")` فقط؛ في PDB الحقيقي غالباً `MSE` يُسجّل `HETATM`، فيُرمى قبل التحويل إلى MET. كذلك تُرمى كل `HETATM` (HEM, MG, ZN, HOH) بلا خيار للحفاظ عليها.

**الدليل:** إدخال `HETATM MSE` أعطى:
```text
ATOM 1 CA ALA ...
Removed 2 HETATM records ...
```
أي أن MSE حُذف ولم يُحوّل، رغم أن docstring يقول "converts selenomethionine (MSE) to MET".

**المرجع:** `wiki/gromacs-protein-ligand-protocol.md:18` و`wiki/receptor-preparation.md:27-33` يطلبان الحفاظ على الماء الجسري والمعادن/المرافقات الوظيفية.

**الأثر:** بروتينات الهيم/المعادن أو MSE تفقد مرافقاتها أو residue كامل، وقد يفشل `pdb2gmx` أو ينتج بروتيناً ناقصاً.

**الإصلاح المقترح:** معالجة `HETATM` لـ MSE قبل الفلترة، وإضافة خيار `keep_cofactors`/`keep_bridging_waters` مع توثيق القرار في manifest.

---

### P1-8: altloc وموديلات PDB المتعددة غير محسومة في receptor prep

**الموقع:** `docking/receptor_prep.py:97-125`.

- لا يوجد اختيار `model 0`؛ `PDBIO.save` يكتب كل الموديلات، وOpenBabel يقرأ أول موديل فقط في اختبارنا.
- لا يوجد اختيار altloc؛ PDB يحتوي `CA A` و`CA B` يمران معاً إلى PDBQT (تحققنا: ذرتان كربون بموضعين مختلفين).

**المرجع:** `wiki/receptor-preparation.md:30-31` يطلب اختيار altloc الأعلى occupancy أو المتوافق مع الموقع.

**الإصلاح المقترح:** تمرير `model=0`، واختيار altloc حسب occupancy أو `A`، وتسجيل عدد altloc المحذوفة.

---

### P1-9: backend Windows native معطّل عملياً

**الموقع:** `md/engine.py:115-120`, `:338-366`, `:447-482`.

`GromacsBackend.shell_command` يرفع استثناءً على `sys.platform == "win32"`، لكن `list_force_fields` و`check_backend_tools` يعتمدان عليه. النتيجة على native Windows:
```text
force fields: []
resolve error: Amber99SB-ILDN is not installed ... installed: none
```
مع أن `gmx.exe` موجود.

**الأثر:** مشروع يدعم `run_app.bat` وvina.exe لا يستطيع تشغيل MD على GROMACS Windows الأصلي؛ يعمل فقط عبر WSL.

**الإصلاح المقترح:** استخدام `Path`/`os.listdir` لاكتشاف `share/gromacs/top/*.ff` على native Windows، وعدم الاعتماد على bash للأدوات؛ أو توثيق WSL كمتطلب صريح.

---

## 5. مشكلات متوسطة (P2)

### P2-1: عناصر PDBQT في تحليل التفاعلات

**الموقع:** `docking/interactions.py:75-80`.

يُقرأ `line[76:78]` كنوع AutoDock لا كعنصر. النتائج:
- `A` (كربون عطري) يبقى `A`، و`carbon_contacts` يشترط `element == "C"` → الكربونات العطرية لا تُحتسب.
- `Cl`/`Br` قد تُقرأ `C`/`B`؛ المعادن قد تُقرأ خطأً.

**الدليل:** erlotinib PDBQT الحقيقي أعطى `A=14`, `C=8`؛ 14 كربوناً عطرياً خارج عدّ التماس الكربوني.

**الإصلاح المقترح:** استخدام `autodock_type_to_element()` من `md/topology.py` أو دالة مشتركة، وإضافة اختبار لـ `A/OA/NA/HD/Cl`.

---

### P2-2: `mutation_engaged_in_contacts` يفشل مع `CHAIN:RESIDUE`

**الموقع:** `docking/interactions.py:215-220`; `docking/ui.py:1843` يمرر `h_pocket_list` بصيغة `A:102`.

`h_all_res_nums` يحتوي `"102"`، و`mut_nums_clean` يحتوي `"A:102"` → `any(...)` false دائماً.

**الأثر:** المؤشر "هل الطفرة مشاركة في التماس؟" يفشل في الواجهة التي تستخدم CHAIN:RESIDUE.

**الإصلاح المقترح:** تطبيع الطرفين إلى `(chain,resnum)` قبل المقارنة.

---

### P2-3: حجم الصندوق والتوفير يُحسبان بقطر ثابت 5 nm

**الموقع:** `md/simulation.py:637`:
```python
compute_rhombic_dodecahedron_volume(cfg.box_distance_nm, 5.0)
```
`5.0` قطر افتراضي وليس قطر البروتين الفعلي. لذلك `dodecahedron_volume_nm3` و`cubic_volume_nm3` و`savings_percent` في manifest ليست لهذا النظام.

**الإصلاح المقترح:** حساب قطر الغلاف من `complex.gro` (max-min) أو من `editconf` output.

---

### P2-4: استخراج الطاقة لا يغطي Kinetic/Total ولا ينبّه عند الفقد

**الموقع:** `md/simulation.py:899-902`.

يطلب `Kinetic En.` و`Total Energy`، لكن ملف `energy_production.xvg` الحقيقي احتوى فقط:
`Potential`, `Temperature`, `Pressure`, `Density`.
لا يوجد تحذير بأن اختيارين لم يُستخرجا.

**الإصلاح المقترح:** التحقق من أسماء legends بعد `gmx energy`، وإضافة تحذير/إعادة محاولة، وتسجيل `n_frames`.

---

### P2-5: `load_mmpbsa_results` يثبّت dielectric/model/frames

**الموقع:** `md/analysis.py:365-387`.

- `internal_dielectric=2.0` دائماً.
- `solvation_model` = "MM-GBSA (GB-OBC2)" لأي EGB، حتى لو استُخدم igb=5/8 أو PB.
- `frames_evaluated=0` دائماً.
- لا يقرأ `-TΔS` ولا يتحقق من entropy.

**الأثر:** manifest/الواجهة يعرضان طريقة غير مؤكدة.

**الإصلاح المقترح:** قراءة `igb`, `eps_int`, `frames`, `-TΔS` من الملف، وعدم افتراض نموذج GB.

---

### P2-6: SASA/RMSD في التحليل الهيكلي يتجاهلان insertion codes

**الموقع:** `analysis/structure.py:153` (`sasa_map[res.id[1]]`)، `:182-193` (`res_dict` بمفتاح `res.id[1]`).

**الأثر:** بقايا `100A` و`100B` تتصادمان؛ SASA وCA-RMSD قد يعيدان قيماً خاطئة.

**الإصلاح المقترح:** مفتاح `(chain, resseq, icode)`، وإرجاع تسميات كاملة.

---

### P2-7: فشل المحاذاة يمكن أن يُسقط واجهة التحليل

**الموقع:** `analysis/alignment.py:96-108` يرفع `ValueError` عند سلاسل فارغة؛ `analysis/ui.py:278` و`:469` تستدعيانها بلا حماية. `get_alignment` نفسها تعيد `("",0,"","")` عند الاستثناء.

**الأثر:** إذا فشلت المحاذاة (PDB غير قياسي/سلسلة فارغة)، تظهر استثناءات بدل رسالة واضحة.

**الإصلاح المقترح:** إرجاع قيم صفرية مع تحذير عند المدخلات الفارغة، أو التحقق قبل الاستدعاء.

---

### P2-8: شارة "validation passed" للإرساء الفردي

**الموقع:** `docking/ui.py:1029` (`if "passed" in val_stat`)، و`:905` تضع `"validation passed (Top-10 recovery)"` عندما ينجح recovery فقط.

**الأثر:** شارة خضراء "تم التحقق بنجاح" رغم فشل top-1؛ الحالة النصية تقول Top-10 لكن الشارة مضللة.

**الإصلاح المقترح:** حالة/شارة منفصلة لـ "Top-10 recovery only" بلون تحذيري، أو اشتراط top-1 للشارة الخضراء.

---

### P2-9: TIP4P-Ew يُبنى بحيلة استبدال include وغير مُتحقق منه

**الموقع:** `md/topology.py:455-531`, `md/engine.py:428-441`.

`pdb2gmx -water tip4p` ثم استبدال `tip4p.itp` بـ`tip4pew.itp`، و`solvate -cs tip4p.gro`. Run الحقيقي نجح، لكن:
- هل `tip4p.gro` متوافق هندسياً مع `tip4pew.itp`؟
- هل `[moleculetype]`/atom order متطابقان؟

**الأثر:** خطر عدم اتساق نموذج الماء/الهندسة، خصوصاً في الكثافة والطاقة.

**الإصلاح المقترح:** التحقق من تطابق atom names/order وgeometry، أو توثيق أن TIP4P-Ew غير مدعوم رسمياً وإزالته من الخيارات.

---

### P2-10: `pdb2gmx` بلا `-his`/`-ter`/stdin

**الموقع:** `md/topology.py:466-475`.

لبروتينات متعددة السلاسل أو histidines حرجة، قد يطلب `pdb2gmx` اختيار protonation/termini تفاعلياً. لم يظهر ذلك في 3HTB، لكنه خطر reproducibility.

**الإصلاح المقترح:** تمرير `-his` صريح، `-ter`، أو stdin محدد، وتسجيل القيم.

---

### P2-11: `method` ثابت وOpenFF alias مضلل

**الموقع:** `md/topology.py:1049` يكتب `method="Amber99SB-ILDN + GAFF2/AM1-BCC"` دائماً؛ `md/engine.py:387` يحوّل `"openff"` إلى `amber99sb-ildn` بدلاً من الفشل.

**الأثر:** manifest قد يصف مساراً غير الذي نُفّذ؛ طلب OpenFF يُنفّذ كـ Amber بصمت.

**الإصلاح المقترح:** اشتقاق method من `force_field`/`water`/ligand route، ورفض OpenFF غير المدعوم.

---

### P2-12: SASA/contacts تشمل المرافقات والربيطة

**الموقع:** `analysis/structure.py:33-44` يحذف الماء فقط، ويُبقي HETATM/ligands؛ `calculate_all_distances` يستخدم كل الذرات.

**الأثر:** SASA لبقايا البروتين يتأثر بوجود ligand/cofactor؛ `min_dist` قد يكون إلى ligand لا protein. قد يكون مقصوداً (complex SASA) لكنه غير موثق.

**الإصلاح المقترح:** توضيح الوضع أو توفير خيار protein-only.

---

## 6. مشكلات منخفضة (P3)

1. **HETATM MSE docstring/test mismatch**: الاختبار يستخدم `ATOM MSE` فقط، فلا يحمي الحالة الحقيقية `HETATM MSE`.
2. **README architecture** يذكر `test_repairs.py` فقط ولا يذكر `test_md.py`, `test_audit_fixes.py`, `test_gmx_integration.py`.
3. **تعارض 23% vs 29.3%**: `AGENTS.md` وwiki يقولان ~23%، والكود يحسب 29.3% (النسبة الهندسية الصحيحة dodecahedron/cube). التوثيق يحتاج توحيد.
4. **Traceability**: بعض وحدات `docking/` (مثل `receptor_prep.py`, `interactions.py`) لا تحتوي docstring يذكر صفحات wiki المستخدمة، خلافاً لـ AGENTS.md.
5. **AppTest لنتائج MD** بعد جولة حقيقية غير موجود؛ الاختبار الحالي يغطي إقلاع التطبيق فقط.
6. **`calculate_rmsd_from_trajectory` يعيد 0.0** عند عدم تطابق عدد الذرات بين الإطارات (`md/analysis.py:199-214`) بدلاً من رفع خطأ/تخطي.
7. **`calculate_rmsf_from_trajectory`** لا يطبق alignment صريحاً؛ يعتمد على أن المسار مُعاد لفّه، وهو غير مضمون إن استُدعيت الدالة خارج الواجهة.
8. **PAINS/aggregators/PAINS filtering** غير موجود في مسار الإدخال، رغم توصية wiki `ligand-preparation-and-protonation.md:122-132` و`docking-failure-modes.md:41-46`.
9. **Vina MW bias** لا يوجد تحذير صريح في الواجهة؛ wiki `docking-failure-modes.md:22-27` ينبه إلى أن الدرجات تعتمد على عدد الذرات.
10. **جودة البنية** (resolution/Rfree/RSCC) لا تُفحص عند تحميل PDB؛ wiki `receptor-preparation.md:13-19` يطلب معايير جودة.

---

## 7. ملاحظات علمية أوسع (ليست أخطاء برمجية)

- **Vina ليس طاقة ارتباط حرة**: الواجهة تقول ذلك في disclaimer، لكن يجب أن يظل واضحاً أن ΔScore ليس ΔΔG.
- **Covalent drugs**: Osimertinib/Nirmatrelvir لا يمكن تمثيلهما تساهمياً في Vina؛ لا يوجد كشف تلقائي أو حجب.
- **Metalloproteins**: `prepare_receptor` يحتفظ بالمعادن لكن PDBQT/charges عامة؛ wiki `docking-failure-modes.md:35-39` يحذر من فشل COMT/metalloenzymes.
- **Crystal waters**: خيار `keep_crystal_waters` موجود في الإرساء، لكن MD `clean_receptor_for_md` لا يوفر خياراً مماثلاً.
- **Biological assembly fallback**: إذا فشل تحميل assembly، تُستخدم asymmetric unit مع تحذير؛ قد تفقد سلسلة/موقعاً. `wiki/receptor-preparation.md:16`.
- **Sampling**: 1 ns افتراضي قد لا يكفي لتحولات pocket؛ لا يوجد تحليل convergence.

---

## 8. خطة إصلاح مقترحة (مرتبة)

### المرحلة 1 — فوري (قبل أي استخدام بحثي)
1. إصلاح مفاتيح `compare_interactions` ↔ `docking/ui.py`.
2. إصلاح قارئ MM-PBSA لدعم `Δ` و`ΔG binding` و`frames` و`-TΔS`.
3. تصحيح أمر `gmx_MMPBSA` في الواجهة.
4. استخدام `em.tpr` كمرجع RMSD الربيطة.
5. فرض `pair_validation_meta["is_valid"]` داخل `run_matched_docking_comparison`.
6. معالجة `HETATM MSE` قبل الفلترة، وإضافة خيار الحفاظ على cofactors/water.

### المرحلة 2 — أسبوع
1. altloc/model selection في `prepare_receptor`.
2. دعم native Windows لاكتشاف حقول القوى والأدوات.
3. إصلاح PDBQT element mapping في `interactions.py` واستخدامه.
4. إضافة replicates + NVE + burn-in + block averages/CI.
5. إصلاح حجم الصندوق ليُحسب من النظام الفعلي.
6. إصلاح energy selections والتحقق من legends.
7. إصلاح insertion-code keys في SASA/RMSD.

### المرحلة 3 — تحسين
1. ربط القوالب السريرية بالواجهة أو إزالتها.
2. إضافة PAINS/aggregator/structure-quality gates.
3. إضافة H-bond geometry/persistence وinteraction energy decomposition.
4. إضافة AppTest لنتائج MD.
5. توحيد توثيق 23%/29.3% وإضافة traceability docstrings.

---

## 9. حدود هذا الفحص

- لم أُجرِ تشغيلاً كاملاً جديداً لـ MD في هذه الجولة؛ استندت إلى:
  - اختبارات التكامل الحقيقية 3/3،
  - سجل/مخرجات `tmp_gmx_test/audit_md_final` الموجودة في المشروع،
  - اختبارات إعادة إنتاج صغيرة للمشكلات المذكورة.
- لم أتحقق من صحة كل مخرجات `gmx_MMPBSA` على نظام حقيقي (الحزمة غير مثبتة)، لكن تحققت من صيغة المخرجات الرسمية وقارنتها بالقارئ.
- لم أعدّل الكود في هذه الجولة؛ التقرير تشخيصي. يمكن تنفيذ الإصلاحات في جولة تالية.

---

## 10. الملفات المرجعية

- [AGENTS.md](AGENTS.md)
- [wiki/gromacs-protein-ligand-protocol.md](wiki/gromacs-protein-ligand-protocol.md)
- [wiki/simulation-and-benchmark-best-practices.md](wiki/simulation-and-benchmark-best-practices.md)
- [wiki/mm-pbsa-and-mm-gbsa.md](wiki/mm-pbsa-and-mm-gbsa.md)
- [wiki/receptor-preparation.md](wiki/receptor-preparation.md)
- [wiki/ligand-preparation-and-protonation.md](wiki/ligand-preparation-and-protonation.md)
- [wiki/autodock-vina.md](wiki/autodock-vina.md)
- [wiki/docking-failure-modes.md](wiki/docking-failure-modes.md)
- [docking/interactions.py](docking/interactions.py)
- [docking/ui.py](docking/ui.py)
- [docking/comparison.py](docking/comparison.py)
- [docking/receptor_prep.py](docking/receptor_prep.py)
- [md/analysis.py](md/analysis.py)
- [md/simulation.py](md/simulation.py)
- [md/topology.py](md/topology.py)
- [md/engine.py](md/engine.py)
- [md/ui.py](md/ui.py)
- [analysis/structure.py](analysis/structure.py)

---

## 11. تحديث: الإصلاحات المنفذة بعد هذا التقرير

نُفذت الحزمة الأولى من الإصلاحات مع اختبارات انحدار:

### تم إصلاحه
- **مفاتيح مقارنة التفاعلات**: `compare_interactions` يعيد المفاتيح القانونية + aliases،
  و`docking/ui.py` يقرأ `healthy_polar_contact_count` وغيرها. أُضيف تطبيع `CHAIN:RESIDUE`
  لـ `mutation_engaged_in_contacts`.
- **PDBQT element/type**: `docking/interactions.py` يستخدم `autodock_type_to_element`؛
  الكربون العطري `A` أصبح `C`.
- **MM-PBSA**: قارئ جديد يدعم `ΔVDWAALS`/`ΔG binding`/`ΔTOTAL`، ويقرأ `igb` و`intdiel`
  و`frames` و`-TΔS`، ولا يخلط مجمع/مستقبل/ربيطة. أمر `gmx_MMPBSA` في الواجهة صُحح إلى
  `-cs/-ci/-cg/-ct/-cp`.
- **RMSD**: مرجع الربيطة أصبح `em.tpr` (بدلاً من `prod_md.tpr`).
- **pair validation**: `run_matched_docking_comparison` يرفض الآن أي `pair_validation_meta`
  غير صالح أو غائب.
- **MSE**: `HETATM MSE` يُحوّل إلى MET/SD قبل pdb2gmx بدلاً من حذفه.
- **altloc/model**: `prepare_receptor` يحتفظ بأول موديل ويختار altloc الأعلى occupancy.
- **Windows native**: اكتشاف حقول القوى عبر المسار، و`check_backend_tools` عبر PATH،
  وتشغيل ACPYPE مباشرة على Windows.
- **الصندوق**: القطر يُحسب من `complex.gro` الفعلي؛ `solute_bounding_diameter_nm` في manifest.
- **الطاقة**: اختيارات `Kinetic-En.`/`Total-Energy` الصحيحة مع تحذير عند فقد أي حد.
- **burn-in**: متوسطات RMSD/H-bonds تستبعد أول 20% (قابل للضبط).
- **NVE drift check**: تشغيل NVE قصير بعد NPT، مع انحدار خطي و`drift_kj_mol_per_atom_per_ps`
  وتحذير عند تجاوز العتبة.
- **أشرطة الخطأ**: `summarize_timeseries` يحسب SD/SEM وblock SEM، ويعرضها في الواجهة.
- **PAINS**: فحص PAINS عبر RDKit FilterCatalog، مع تحذير في الواجهة.
- **جودة البنية**: `fetch_entry_quality` يجلب resolution/R_free/method من RCSB ويعرض تحذيرات.
- **replicates مستقلة بالكامل**: `n_replicates` (1–5)؛ كل مسار إضافي يعيد NVT/NPT/الإنتاج
  من بنية EM ببذرة سرعات مستقلة، ويُجمع RMSD/H-bonds بين المسارات (mean ± SD).
- **replicates في التصدير**: حزمة `run_md.sh`/`run_md.bat` تحتوي `nvt_repN.mdp` وحلقة NVT/NPT/إنتاج.
- **حارس المرافقات**: حذف HEM/FE/MG/ZN/MN/CU/CA أصبح **fail-closed**؛ يحتاج
  `allow_cofactor_removal=True` صريحاً، وإلا يُرفض MD برسالة واضحة.
- **SASA protein-only**: تُحسب SASA على بقايا البروتين فقط افتراضياً، مع الاحتفاظ
  بالمرافقات/الربيطة في البنية المعروضة و3D viewer.
- **contacts protein-only**: `calculate_all_distances` يستبعد HETATM افتراضياً،
  مع `include_hetero=True` لخريطة تماس المركب.
- **pdb2gmx غير تفاعلي**: `-nointer` و`-chainsep id_or_ter` مضافان صراحةً، ويُتحقق منهما
  في سجل التكامل الحقيقي.
- **اختبار MD كامل**: `test_short_full_md_pipeline` يشغّل pdb2gmx + ACPYPE + GROMACS
  (EM/NVE/NVT/NPT/إنتاج/تحليل) ويتحقق من manifest والطاقة وRMSD/H-bond.
- **insertion codes**: SASA وCA-RMSD يستخدمان `(resnum, icode)`.
- **identity edge case**: المحاذاة الفارغة تعيد metrics غير صالحة بدلاً من رفع استثناء.
- **UI badge**: `Top-10 recovery` لم يعد يُعرض كتحقق أخضر كامل.
- **method/OpenFF**: method مشتق من FF الفعلي، وOpenFF لم يعد يُحوّل بصمت إلى Amber.
- **AppTest لنتائج MD**: `tests/test_md_ui.py` يفتح تبويب النتائج بحالة جلسة محقونة
  ويتحقق من عدم وجود استثناءات ومن ظهور حالة MM-PBSA الصادقة.
- **TIP4P-Ew موثّق ومُتحقق**: `validate_water_model_topology` يفحص `SOL` و4 ذرات
  وتطابق `tip4p.gro`، ويُسجّل النتيجة في `TopologyResult.water_validation` وmanifest.
- **RSCC للربيطة والبروتين**: `parse_ligand_rscc` و`parse_residue_rscc` يقرآن تقرير
  التحقق الرسمي؛ يُعرض RSCC للربيطة وملخص بقايا البروتين (median، <0.8) وقائمة
  البقايا الأقل من 0.8، مع تحذير. القيمة ليست درجة إرساء.
- **Aggregator محسّن**: عتبة logP > 3.5 (وفق wiki) + MW/TPSA + تنبيهات بنيوية
  (rhodanine/quinone/catechol/flavone)، مع توضيح أنها ليست Aggregator Advisor.
- **احتفاظ الأيونات الأحادية**: ZN/MG/CA/K تُضاف إلى complex.gro/topol.top وتمر عبر
  solvate/genion؛ HEM/FE/MN/CU ما زالت fail-closed لعدم توفر moleculetype مناسب.

### لا يزال مفتوحاً
- **توصيف المرافقات متعددة الذرات والماء الجسري**: HEM/المعادن غير المدعومة والماء
  البلوري تحتاج توصيفاً صريحاً؛ الاحتفاظ الحالي يغطي الأيونات الأحادية فقط.
- Aggregator Advisor: التشابه مع قاعدة بيانات المجمعات المعروفة غير متاح؛ الموجود
  تنبيهات بنيوية/فيزيائية شفافة فقط.

**الاختبارات بعد الإصلاح:** 141 اختباراً، OK (skipped=6)، و6/6 تكامل GROMACS.

