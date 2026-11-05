# 🧬 Bio-Impact Analyzer
> **منصة النمذجة الجزيئية المتقدمة ومحاكاة الارتباط الدوائي (Molecular Modeling & Mutational Docking)**

---

## 📌 نظرة عامة (Overview)
منصة **Bio-Impact Analyzer** هي أداة حاسوبية وبيولوجية متكاملة مصممة للصيادلة، أطباء الأورام، والباحثين في علم الأحياء الجزيئي وتصميم الأدوية، تتيح:
1. **التحليل الهيكلي ومقارنة الطفرات**: مقارنة البنية البلورية للبروتين السليم (Wild-Type) مع البروتين المصاب (Mutant)، مع حساب التغير في مساحة السطح المعرضة للمذيب (SASA)، محاذاة التسلسل (Sequence Alignment)، وتصنيف الأثر الكيميائي والفيزيائي للطفرة مع عرض ثلاثي أبعاد تفاعلي بواسطة `py3Dmol`.
2. **محاكاة الإرساء الجزيئي المزدوج (Matched WT-vs-Variant Docking)**: إرساء الجزيئات الدوائية على الجيوب المستهدفة ومقارنة السليم بالطافر مع التحقق الصارم من تطابق التسلسل ($\ge 90\%$) ومحاذاة البقايا المتجانسة وحساب $\Delta \text{Score}$.
3. **الإرساء المتقاطع الاستكشاري (Exploratory Cross-Protein Docking)**: إرساء ومقارنة أي بروتينين مستقلين بشكل متوازٍ ومنفصل دون فرض تجانس جيني أو حساب $\Delta \text{Score}$ مضلل.
4. **الموثوقية والتحقق الإحصائي (Scientific Rigor & Reproducibility)**: دعم البذور العشوائية المتعددة (Multi-seed Docking)، تجميع الوضعيات في عناقيد إحصائية (Pose Clustering)، والتحقق من إعادة الإرساء (Redocking Validation) وتدقيق هوية الربائط عبر قاموس RCSB.

---

## 🏗️ هيكلية المشروع (Project Architecture)

```
BIOBIO/
│
├── app.py                       # نقطة الانطلاق الرئيسية للتطبيق (Streamlit)
│
├── analysis/                    # حزمة التحليل الهيكلي والمقارنة الجزيئية
│   ├── __init__.py              # تصدير دوال الحزمة
│   ├── constants.py             # ثوابت الأحماض الأمينية والخواص الفيزيوكيميائية
│   ├── fetch.py                 # جلب هياكل PDB (Deposited & Assembly) وقاعدة الطفرات
│   ├── structure.py             # قراءة البنية، فحص السلاسل، خوارزميات KDTree، وحسابات SASA
│   ├── alignment.py             # خوارزميات المحاذاة (Biopython) وتحليل أثر الطفرات
│   ├── visualization.py         # العرض ثلاثي الأبعاد (py3Dmol) ورسم Plotly البياني
│   └── ui.py                    # واجهة المستخدم الخاصة بالتحليل الهيكلي وعزل الجلسات
│
├── docking/                     # حزمة الإرساء الجزيئي والتحليل الدوائي
│   ├── __init__.py
│   ├── ui.py                    # واجهة المستخدم للإرساء (Matched & Cross-docking)
│   ├── engine.py                # إدارة تشغيل محرك AutoDock Vina
│   ├── binding_site.py          # حساب صندوق الجيب (Grid Box)
│   ├── receptor_prep.py         # تهيئة وتجهيز البروتين المستقبل وتصحيح أعمدة العناصر
│   ├── ligand_prep.py           # تهيئة الليجاند وتحويل SMILES إلى PDBQT
│   ├── openbabel_io.py          # معالجة آمنة لملفات وتنسيقات OpenBabel
│   ├── pair_validation.py       # التحقق البيولوجي من تطابق الأزواج ومحاذاة بقايا الجيب
│   ├── comparison.py            # المقارنة المزدوجة المتطابقة والإرساء المتقاطع
│   ├── interactions.py          # رصد الروابط الهيدروجينية ومسافات التلامس
│   ├── presets.py               # مكتبة الحالات السريرية البشرية
│   ├── pubchem.py               # جلب الأدوية والمعلومات الجزيئية من PubChem
│   ├── rcsb.py                  # تنزيل ملفات الـ PDB وجلب صيغ الربائط من قاموس RCSB
│   ├── validation.py            # فحص RMSD وتجميع العناقيد والتحقق من هوية الربائط
│   ├── storage.py               # إدارة مسارات التشغيل والمخرجات وتدقيق الـ Manifests
│   └── models.py                # نماذج وهياكل البيانات (Dataclasses)
│
├── data/                        # البيانات والمراجع المحلية
│   └── mutation.json            # قاعدة بيانات ربط البروتينات السليمة بالمصابة تلقائياً
│
├── tools/                       # الأدوات التنفيذية
│   └── vina/
│       ├── vina                 # محرك AutoDock Vina (Linux)
│       └── vina.exe             # محرك AutoDock Vina (Windows)
│
├── tests/                       # حزمة الاختبارات الآلية والتحقق العلمي
│   ├── __init__.py
│   └── test_repairs.py          # اختبارات الإرساء والتحليل، مع ملفات اختبار إضافية لـ MD والتدقيق
│
├── docking_runs/                # مجلد حفظ مخرجات عمليات الإرساء
│   └── .gitkeep
│
├── requirements.txt             # اعتمادات ومكتبات بايثون
└── .gitignore                   # استثناءات Git
```

---

## 🚀 طريقة التشغيل (How to Run)

### 1. تثبيت الحزم والمتطلبات:
على Windows شغّل `setup_app.bat` لإنشاء البيئة المستقلة `.venv` وتثبيت متطلبات المشروع.
ثم شغّل `run_app.bat`؛ يستخدم هذا المشغّل Python الخاص بالمشروع.
البيئة المحلية الحالية تستفيد من حزم Python العامة المطابقة للمتطلبات عبر
`include-system-site-packages`، مع تثبيت NumPy وStreamlit وACPYPE داخل `.venv`.
لم تُغيّر نسخ الحزم العامة. نسخة RDKit المعتمدة في المتطلبات هي `2026.3.6`.
للتثبيت اليدوي من PowerShell:

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe -m streamlit run app.py --server.address 127.0.0.1
```

على Linux:
```bash
pip install -r requirements.txt
```

### 2. تشغيل التطبيق:
```bash
streamlit run app.py
```

---

## 🧪 التحقق والاختبار (Verification)

### 1. اختبار استيراد الوحدات الأساسية:
```bash
python -c "import analysis; import docking; print('All modules loaded successfully!')"
```

### 2. تشغيل حزمة الاختبارات الآلية (149 اختباراً؛ 6 منها تكامل GROMACS اختياري):
```bash
python -m unittest discover tests
```

### 3. محاكاة الديناميكا الجزيئية (GROMACS):
تتطلب المحاكاة تثبيت GROMACS (نظامياً أو داخل WSL على Windows) وACPYPE/AmberTools.
تثبيت ACPYPE في بيئة Windows لا يثبت GROMACS ولا يضمن توفر ثنائيات AmberTools.
يجب أن تتوفر أدوات التوصيف في نفس البيئة التي ينفذ فيها GROMACS.
إذا أعاد WSL خطأ `REGDB_E_CLASSNOTREG` أو `Class not registered`، افتح PowerShell
**كمسؤول** ونفذ تحديث/تثبيت حزمة Microsoft الرسمية، ثم تثبيت Ubuntu:

```powershell
winget install --id Microsoft.WSL --exact --source winget
wsl --install -d Ubuntu
```

قد يلزم إعادة تشغيل Windows وإكمال إنشاء مستخدم Ubuntu. لا يحذف هذا الإجراء توزيعاتك.
بعدها يجب تجهيز GROMACS وACPYPE/AmberTools داخل Ubuntu والتحقق من `gmx --version`
ومن توفر `acpype` و`antechamber` و`parmchk2` و`tleap` أو بدائلها المدمجة في ACPYPE.
التوثيق الرسمي: [WSL](https://learn.microsoft.com/windows/wsl/install)،
[GROMACS](https://manual.gromacs.org/documentation/current/install-guide/index.html)،
[ACPYPE](https://pypi.org/project/acpype/).

لا يوجد مسار بديل يُصنّع مساراً أو طاقات: عند غياب GROMACS تفشل الجولة صراحةً
ويُكتب `md_manifest.json` بحالة `failed`. اختبار التكامل الكامل (pdb2gmx + ACPYPE +
grompp + mdrun قصير) يتم عبر:

```bash
set BIO_RUN_GMX_INTEGRATION=1
python -m unittest tests.test_gmx_integration -v
```

---

## ⚠️ حدود علمية معروفة (Known Limitations)

- **الإرساء** يعطي فرضية ارتباط وترتيباً داخل بروتوكول محدد، وليس طاقة ارتباط حرة مقاسة.
- **الأدوية التساهمية** (Osimertinib، Nirmatrelvir) لا يمكن تمثيل رابطنها التساهمية في Vina؛
  النتائج تقريب غير تساهمي ويجب قراءتها بهذا القيد.
- **الإرساء الأعمى** (whole-protein grid) منخفض الثقة ولا يُعد بديلاً عن جيب متحقق.
- **الماء البلوري** يُزال افتراضياً (يمكن الاحتفاظ به عبر خيار Keep crystal waters)؛
  الجيوب التي تعتمد على ماء جسري تحتاج تفعيل هذا الخيار وتوثيقه في الـ manifest.
- **المرافقات المعدنية**: في الإرساء تُعالَج بشحنات Gasteiger عامة. في MD، حذف
  HEM/FE/MG/ZN/MN/CU/CA **مرفوض افتراضياً** (fail-closed)؛ فعّل `allow_cofactor_removal`
  فقط بعد التأكد من أن المرافق غير مطلوب، ولا يُدّعى توصيف حقل قوى له.
- **احتفاظ الأيونات الأحادية**: ZN/MG/CA/K يمكن الاحتفاظ بها في
  `complex.gro`/`topol.top` عبر `retain_monatomic_ions`؛ HEM/MN/CU/FE والماء البلوري
  ما زالت تحتاج توصيفاً صريحاً.
- **PAINS/aggregation**: يُعرض تنبيه PAINS عبر RDKit، وتحذير تجمّع **heuristic** وفق
  logP > 3.5 + MW/TPSA + تنبيهات بنيوية (rhodanine/quinone/catechol/flavone)؛
  هذا ليس Aggregator Advisor ولا دليلاً تجريبياً على التجمّع.
- **SASA/contacts**: تُحسب SASA وخريطة التماس افتراضياً على بقايا البروتين فقط، مع إبقاء
  المرافقات/الربيطة للعرض؛ استخدم `include_hetero_in_sasa` / `include_hetero` للمركب الكامل.
- **جودة البنية**: تُعرض resolution/R_free/method من RCSB، و**RSCC للربيطة المرجعية** وملخص
  **RSCC لبقايا البروتين** (median، وعدد < 0.8) مع قائمة البقايا الأقل من 0.8، من تقرير
  التحقق الرسمي، مع تحذير عند الانخفاض.
- **MM-PB(GB)SA** لا تُحسب داخل التطبيق؛ تُعرض فقط عند وجود `FINAL_RESULTS_MMPBSA.dat`
  حقيقي من `gmx_MMPBSA`.
- **المحاكاة**: يوجد **فحص NVE drift** بعد NPT، وأشرطة خطأ داخل المسار (SD/SEM وblock SEM)
  مع burn-in 20%. يدعم التطبيق **n_replicates = 1–5**؛ كل مسار إضافي يعيد NVT/NPT/الإنتاج
  من بنية EM ببذرة سرعات مستقلة، مع تجميع RMSD/H-bonds بين المسارات. الأرقام تصلح للمقارنة
  النسبية ضمن نفس البروتوكول.
- **الصندوق**: التوفير الهندسي المحسوب لـ dodecahedron مقابل cube هو ~29.3%؛ قيمة ~23%
  المذكورة في `AGENTS.md` هي قيمة تحفظية للتخطيط، وليست الناتج الحسابي للملف.
