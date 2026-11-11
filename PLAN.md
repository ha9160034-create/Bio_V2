# PLAN.md — خطة تقوية BIOBIO للإصدار 1.0 (نسخة كاملة قابلة للنشر)

**الحالة**: خطة معتمدة للتنفيذ — لم يُعدَّل أي كود بعد.
**النطاق**: المراحل 0 → 5 (نسخة كاملة قابلة للنشر)، مع بوّابة قبول آلية.
**المرجع العلمي**: `wiki/**` — كل قرار فيزيائي مُشتق منها ومقتبس.
**آخر تحديث**: 2026-11-10

---

## كيف تُقرأ هذه الخطة

| الرمز | المعنى |
|:--|:--|
| **P0** | مانع للنشر — لا تُشارك النسخة قبل إصلاحه |
| **P1** | مطلوب للموثوقية العلمية |
| **P2** | تحسين جودة/سرعة، لا يمنع النشر |
| **مُجمَّد** | لا يُلمس إطلاقاً (انظر §8) |
| **قبول** | شرط قابل للقياس يثبت نجاح المهمة |

---

## 0. ملخص تنفيذي

الفيزياء في المشروع **سليمة ومطابقة للويكي** (force field، خطوات الزمن، الضوابط الحرارية والضغطية، الـ cutoffs، إعادة اللف الدوري، معايير التحليل). ما يجعل النسخة غير قابلة للنشر اليوم ليس العلم، بل **قابلية التكرار**:

1. لا يوجد تثبيت إصدارات (RDKit ثلاث نسخ، Vina نسختان، ACPYPE نسختان، GROMACS غير مثبّت).
2. مسار طوبولوجيا الليكند غير حتمي: ParmEd يستبدل مخرجات ACPYPE ويتجاوز معايرة الشحنة (+0.002 e → تحذير PME).
3. الصورة المنشورة `:latest` قد لا تطابق المستودع.
4. `gmx` مربوط مباشرة ببناء AVX2_256 → ينهار على معالجات أقدم.
5. الـ EXE يعتمد على Python + pywebview على جهاز المستخدم.

الخطة تعالج هذه الخمس دون المساس بأي قرار نجح سابقاً (§8).

---

## 1. القرارات المغلقة

### 1.1 القرار الأول — إصدار RDKit: **تثبيت `2026.09.1`**

**السبب المعماري**: نسأل «أي نسخة تعطي نتائج علمية صحيحة؟» والسؤال الصحيح هو «أي نسخة **مُثبَّتة ومُتحقَّق منها بقياس**؟». إصدار RDKit لا يحدّد صحة الفيزياء؛ يحدّد ثلاثة أشياء فقط:

| ما يتأثر بـ RDKit | الملف | الأثر |
|:--|:--|:--|
| الـ conformer الابتدائي (ETKDGv3, seed 42) | `docking/ligand_prep.py:26-27` | يغيّر نقطة انطلاق بحث Vina |
| خريطة الذرات لحساب RMSD (MCS, `CompareAny`) | `docking/validation.py:216` | يغيّر رقم RMSD المُبلَّغ |
| canonical SMILES / الصيغة / الشحنة الرسمية | `docking/ligand_prep.py:86-148` | يغيّر التحقق المدخلي |

**ولا يتأثر بها**: أي شيء في MD (GROMACS لا يرى RDKit)، ولا البروتونيشن (يقوم به OpenBabel عند pH 7.4، `docking/ligand_prep.py:37-63`)، ولا دالة التسجيل (Vina).

**الأدلة المؤيدة لـ 2026.09.1**:

1. **هو الإصدار الوحيد الذي أنتج أدلة متحقَّقة في هذا المشروع**: 7 تشغيلات إرساء + redocking RMSD = 1.598 Å (rank 2). التبديل يُلغي هذه الأدلة بلا مقابل مقاس.
2. **هو خط الإصلاحات الأحدث** (سبتمبر فوق مارس): يحتوي تراكم إصلاحات ستة أشهر على مسارات نستخدمها فعلاً (ETKDG، MCS، sanitization).
3. **`2026.03.x` ليس «أكثر صحة»** — هو فقط أقدم؛ ولا يوجد دليل أن فيه إصلاحاً غائباً عن سبتمبر.
4. **التكلفة**: اختيار 2026.03 يفرض إعادة تحقق كاملة بلا فائدة مقاسة.

**نقطة الأمانة**: لم أتمكن من التحقق من ملاحظات إصدار RDKit أونلاين (لا وصول شبكي في هذه الجلسة). لذلك لا أعتمد على «الحداثة» كحجة، بل على **بوّابة قياس**:

> **قاعدة ملزمة**: يُولَّد «بصمة ذهبية» (golden fingerprint) لأول بناء ناجح، وتُحفظ في `tests/golden/rdkit_fingerprint.json`. أي تغيير في البصمة = فشل البناء، ويجب إعادة التحقق كاملاً أو الرجوع للإصدار السابق.

**بصمة مرجعية فعلية قِستها على `2026.03.5`** (لتثبت أن القياس ممكن ومستقر). الجزيء المرجعي هو **الإرلوتينيب الصحيح** `C#Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCCOC)cc23)c1` (CID 176870، 29 ذرة ثقيلة) — وهو مطابق لـ `prepared_smiles` في manifest الإرساء:

```json
{
  "canonical_smiles": "C#Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCCOC)cc23)c1",
  "formula": "C22H23N3O4",
  "heavy_atoms": 29,
  "formal_charge": 0,
  "hbd_hba": [1, 7],
  "rings": 3,
  "ETKDGv3":      { "embed_code": 0, "mmff_energy": 92.4289, "coord_sha256_3dp": "ebc72427a30d7a13" },
  "srETKDGv3":    { "embed_code": 0, "mmff_energy": 92.4289, "coord_sha256_3dp": "ebc72427a30d7a13" },
  "mcs":          { "num_atoms": 29, "num_bonds": 31, "smarts_sha256": "32a1d956bfce178a" }
}
```

> **تصحيح موثّق (2026-11-10)**: البصمة الأولى في هذه الخطة كانت محسوبة على SMILES ناقص مجموعة إيثوكسي واحدة (`C20H19N3O3`، 26 ذرة)، أي **جزيء مختلف**. الصيغة أعلاه هي الصحيحة والمطابقة لمدخلات الإرساء الفعلية. سبب التصحيح: `docking_runs/.../manifest.json → ligand.prepared_smiles` و`heavy_atom_count = 29`.

**قاعدة الرجوع**: إذا أخفقت البصمة على 2026.09.1، أو أخفق redocking المرجعي، نعود إلى `2026.03.x` **ونعيد التحقق كاملاً**. لا خلط بين إصدارين في أي حال.

**ملاحظة تسمية**: حزمة pip تُسمّى `2026.9.1`، بينما `rdkit.__version__` يُرجع `2026.09.1`. كلاهما نفس الإصدار؛ الوثائق ستستخدم `2026.09.1`.

### 1.2 القرار الثاني — النطاق: **نسخة كاملة قابلة للنشر (المراحل 0–5)**

يشمل: تثبيت الإصدارات، حتمية الطوبولوجيا، نقل وتشغيل، مطابقة علمية، بوّابة قبول، حزمة نشر بدون إنترنت.

### 1.3 المبدأ الحاكم

> **لا شيء يعتمد على `latest`، ولا شيء يُقبل بلا قياس، ولا مسار فاشل يكتب نتيجة صامتاً.**

---

## 2. خط الأساس المقاس (Baseline — مرجع لا جدال فيه)

| القياس | القيمة | الدليل |
|:--|:--|:--|
| GROMACS | `2026.3-conda_forge` | `md_runs/md_20261109_142719_f90275/nvt.log` |
| بناءات SIMD المتاحة | `AVX2_256`, `AVX_256`, `SSE2` | `ls /home/dell/gmx_fast/bin.*` |
| موزّع `gmx` | سكربت يقرأ `/proc/cpuinfo` ويعود AVX2_256→AVX_256→SSE2 | `head /home/dell/gmx_fast/bin/gmx` |
| `amber99sb-ildn` defaults | `1 2 yes 0.5 0.8333` | `amber99sb-ildn.ff/forcefield.itp` |
| الأداء (خط الأساس) | **4.840 ns/day** (8 خيوط) | `artifacts/baseline/em/nvt_perf.log` (5000 خطوة، 10 ps) |
| إعادة إنتاج الإرساء (بذرة 42) | **مطابق بالبايت**: `-6.822 / -6.656 / -6.622`، hash `1585e6e7…fa013` | `artifacts/baseline/docking/` |
| إعادة إنتاج EM | **مطابق**: 1489 خطوة، `Fmax = 9.2863733e+02`، `Ep = -1.6638664e+06` | `artifacts/baseline/em/em.log` |
| شحنة البروتين | −8.000 e | `pdb2gmx.log` |
| شحنة الليكند | **+0.002 e** (خطأ) | `ERLOT_AC.prmtop` + ITP |
| شحنة النظام | +0.002 e → تحذير PME | `grompp_nvt.log` |
| EM | تقارب عند Fmax = 928.6 بعد 1489 خطوة | `em.log` |
| الحالة الأخيرة لـ NVT | توقف عند 18000/50000 (36 ps) | `nvt.log` |
| الصندوق | dodecahedron = 1049.08 nm³ | `editconf.log` |
| solute diameter | 9.405 nm (ذيل C-terminal 977–995) | `editconf.log` |
| بديل triclinic | 706.08 nm³ (−33%) | قياس مباشر |
| بديل cubic | 1483.62 nm³ | قياس مباشر |
| الماء والأيونات | 32283 SOL، 103 NA، 95 CL، 102,128 ذرة | `solvate.log` / `genion.log` |
| Vina في Docker | 1.2.5 (`/app/tools/vina/vina`) | 7 manifests |
| Vina على Windows | 1.2.7 | `tools/vina/vina.exe --version` |
| RDKit في Docker | 2026.09.1 | manifests |
| RDKit في `.venv` | 2026.03.5 | قياس مباشر |
| ACPYPE في الحاوية | 2023.10.27 | `acpype.log` |
| ACPYPE في `.venv` | 2026.10.15 | قياس مباشر |
| OpenBabel | 3.1.0 (الطرفان) | manifests + `.venv` |
| Python | 3.11.18 (Docker) / 3.12.3 (venv) | — |
| أفضل درجة إرساء (1M17/Erlotinib) | −6.822 / −6.520 / −6.565 kcal/mol | manifests |
| Redocking top-1 RMSD | **5.824 Å (فشل > 2.0 Å)** | manifest |
| Redocking rank-2 RMSD | 1.598 Å (نجاح) | manifest |
| حجم EXE | 10,752 بايت | `BIOBIO.exe` |
| تشغيلات MD | 15 مجلداً، 14 بلا نتائج | `md_runs/` |
| تشغيلات الإرساء | 7، كلها 1M17 | `docking_runs/` |

---

## 3. سجل المشاكل الكامل

### 3.1 من تقرير الفحص (مشاكل أثبتها القياس)

| المعرّف | المشكلة | السبب الجذري | الحل | قبول |
|:--|:--|:--|:--|:--|
| **P0-1** | لا تثبيت إصدارات | `FROM ...:latest`، conda/pip بلا pin | `§5` جدول تثبيت + `constraints.txt` | `verify_stack.py` يطابق `expected_versions.json` |
| **P0-2** | Vina مزدوج/ثلاثي | apt `autodock-vina` + binary مرفق 1.2.5 + Windows 1.2.7 | حذف apt، تثبيت 1.2.5 بـ SHA-256 | نفس البصمة من كل مسار |
| **P0-3** | ParmEd يستبدل ACPYPE ويُفسد معايرة الشحنة | الغلاف يكتب فوق النجاح؛ لا معايرة | منطق أربع حالات + `balanceCharges` + حارس | `sum(q)=0`، صفر تحذير PME |
| **P1-4** | `gmx` مربوط بـ AVX2_256 | `ln -sf .../bin.AVX2_256/gmx` | استخدام الموزّع `/opt/conda/bin/gmx` | يعمل على AVX2 وغير AVX2 |
| **P1-5** | الصورة قد لا تطابق المستودع | `:latest` + غلاف مختلف داخل الصورة | وسم ثابت + digest | EXE يرفض digest مخالفاً |
| **P1-6** | OpenBabel مثبّت مرتين (conda + wheel) | `Dockerfile:15` و`:41` | تثبيت واحد + إثبات أن CLI وbindings يعملان | `obabel -V` + `import openbabel` |
| **P1-7** | ترقيع `acpype/utils.py` | `printf >> utils.py` | إلغاء الترقيع؛ `topology.py` يقرأ `$AMBERHOME` | لا تعديل على مكتبة طرف ثالث |
| **P2-8** | EXE يعتمد على Python + pywebview | `Launcher.cs:188` | msedge صريح + بروفايل معزول | نافذة مستقلة بلا Python |
| **P2-9** | بطء I/O عبر `/mnt/c` | bind mount من Windows | مخرجات داخل volume مُدار + نسخ نهائي | زمن I/O أقل (قياس) |

### 3.2 من دليل التحويل — مع تصحيح التشخيص

| # | العطل | الحل المعتمد | هل يعيده إصلاحنا؟ | الحارس |
|:--|:--|:--|:--|:--|
| 1 | MD بطيء (SSE4.1 + `-v` + 1 ns×3) | conda AVX2 + `OMP_NUM_THREADS` + بلا `-v` | **لا** — الموزّع يختار AVX2_256 على نفس المعالج | `gmx --version` + ns/day ≥ الأساس |
| 2 | تعارض conda | بيئتان منفصلتان | **لا** — الحذف يخصّ apt vina فقط | البناء ينجح، الأداتان تعملان |
| 3 | `open Dockerfile` | مسار البناء | **لا** — `build.ps1` بـ `Set-Location $PSScriptRoot` | أمر بناء واحد |
| 4 | PyInstaller `%TEMP%` | `csc.exe` | **لا** — لا نلمسه | لا `*.spec` في المستودع |
| 5 | `Please open Docker Desktop` كذباً | `IsAdmin()` + `runas` | **لا** — نُبقيه | — |
| 6 | `GMXRC.bash: line 85` | `ENTRYPOINT /entrypoint.sh` | **لا** — نُبقيه | الحاوية تقلع بلا stderr |
| 7 | Brave/Edge يدمج النافذة | pywebview | **لا** — مسار msedge صريح + `--user-data-dir` معزول | النافذة مستقلة وBrave مفتوح |
| 8 | ACPYPE غير مثبت + `SyntaxError` | shebang + دالة `bundled_amber_dir` | **لا** — نُلغي الترقيع ونقرأ `$AMBERHOME` | `command -v acpype` صحيح |
| 9 | `free(): invalid pointer` | ParmEd converter | **لا** — لكن بتصحيح: لا كتابة فوق النجاح + معايرة شحنة | `sum(q)=0` |
| 10 | `FROM` parse + حذف `md/indexing.py` | ملفات منفصلة + إبقاء `indexing.py` | **لا** — `.dockerignore` بدل `rm -rf` | `import md.simulation` ينجح |

**تصحيحات على الدليل نفسه (للنسخة 2.0 منه):**

1. **#9 تشخيص ناقص**: «ينهار قبل حفظ `_GMX`» غير صحيح دائماً. في 2026-11-09 طُبعت رسالة `SUCCESS: Generated ..._GMX.itp` ثم استُبدلت المخرجات. الصح: «ينهار أحياناً قبل أو بعد الحفظ ⇒ يجب الحكم على الملفات لا على رمز الخروج».
2. **#9 أثر علمي غير مذكور**: مسار ParmEd يتجاوز `balanceCharges()` → +0.002 e → تحذير PME. يجب توثيقه وحلّه.
3. **#1 استشهاد غير دقيق**: `wiki/poli-2020-mm-pbsa-virtual-screening.md:30` يدعم 0.5–5 ns **كمسار استخراج لقطات MM-PBSA**، لا كدليل استقرار. الصياغة تحتاج تحديد السياق.
4. **#2 `autodock-vina` عبر apt**: يمكن حذفه بأمان؛ وجوده أنشأ نسخة Vina ثالثة غير مثبّتة.
5. **ادعاء «يعمل بنجاح 100%»**: لا يتوافق مع 14 من 15 مجلد MD بلا نتائج، وآخر NVT متوقف عند 36/100 ps. الأدق: «يعمل، مع قيود معروفة».

### 3.3 مشاكل لم تُذكر سابقاً (اكتُشفت أثناء الفحص)

| المعرّف | المشكلة | الأثر | الحل |
|:--|:--|:--|:--|
| **N-1** | `md/simulation.py:81` `prod_duration_ns = 1.0` بينما الواجهة 0.5 | مصدران للحقيقة → نتائج غير متوقعة عند الاستدعاء البرمجي | مصدر واحد: الواجهة تقرأ من dataclass |
| **N-2** | `parmed` مثبّت في `base` و`amber` بلا pin | نسختان قد تختلفان → ITP مختلف | تثبيت واحد موحّد |
| **N-3** | `openbabel` conda + `openbabel-wheel` pip | ABI مزدوج محتمل | واحد فقط (§Phase 1) |
| **N-4** | لا سياسة احتفاظ: 15 مجلد MD + 7 إرساء | تضخّم القرص على جهاز الصديق | `--keep-last N` + زر حذف |
| **N-5** | `docker run` بلا `--rm` | حاويات يتيمة عند انهيار المشغّل | `--rm` + `--label biobio=1` |
| **N-6** | المنفذ 8501 ثابت | فشل إن كان مشغولاً | فحص + اختيار منفذ بديل |
| **N-7** | `normalize_manifest_status` يرفع `ValueError` | تشغيل قديم يكسر الواجهة | قبول القيم القديمة كـ legacy |
| **N-8** | مسارات/أسماء عربية أو بمسافات | قد تكسر `pdb2gmx`/`editconf` | تطهير الأسماء إلى ASCII + تحقق مسبق |
| **N-9** | EXE يكتب بجانب نفسه | فشل في `Program Files` أو OneDrive | `%LOCALAPPDATA%\BIOBIO` افتراضياً |
| **N-10** | لا `md_manifest.json` عند النجاح | لا provenance | كتابته دائماً |
| **N-11** | `RunCmd` بلا مهلة أثناء `pull` | تعليق بلا نهاية عند شبكة بطيئة | مهلة + رسالة تقدم |
| **N-12** | `OMP_NUM_THREADS=$(nproc)` قد يفرط بالاشتراك | بطء على 4c/8t | إمكانية تجاوز عبر متغير بيئة |
| **N-13** | ترويسة ITP تقول «November 1» وmtime 9 نوفمبر | شك في provenance | تسجيل `generated_at` وhash للمدخلات |
| **N-14** | `convert_amber_to_gmx.py` بلا `main()` ولا تحقق مدخلات | يعمل نسبةً لـ CWD | إعادة هيكلة + `argparse` |
| **N-15** | `acpype_wrapper.py` يتجاهل رمز خروج acpype | فشل صامت | تسجيل + قرار صريح |
| **N-16** | `.dockerignore` لا يستثني `build/`, `dist/`, `*.exe`, `*.spec` | صورة متضخمة + تسريب أدوات | allowlist (§Phase 3) |
| **N-17** | لا `HEALTHCHECK` | EXE لا يعرف إن كانت الحاوية مريضة | `HEALTHCHECK` على 8501 |
| **N-18** | `BIOBIO.spec`/`launcher.py`/`native_window.py` باقية | قد يُعاد بناء EXE معطوب | حذفها (§Phase 3) |
| **N-19** | لا تحقق مسبق من WSL2/القرص/الذاكرة | تجربة أولى غامضة للصديق | preflight في الـ EXE |
| **N-20** | `EXPOSE 8501` + `0.0.0.0` | مقبول، لكن يجب توثيقه | توثيق في `SECURITY.md` |

---

## 4. المراحل

### المرحلة 0 — خط الأساس والرجوع (0.5 يوم)

**المهام**
1. `git tag baseline-before-hardening`
2. تسجيل digest الصورة الحالية: `docker image inspect --format '{{index .RepoDigests 0}}' hassancreator/biobio-app:latest`
3. `docker save hassancreator/biobio-app:latest -o artifacts/biobio-baseline.tar`
4. تشغيل مرجعي: إرساء بذرة 42 + EM، وحفظ النواتج في `artifacts/baseline/`
5. إنشاء `baseline.json`: الإصدارات + digest + 4.840 ns/day + شحنة +1M17.

**المخرج**: `baseline.json` + أرشيف رجوع.
**قبول**: `docker load` للأرشيف يعيد نفس digest.

---

### المرحلة 1 — P0: التكرار العلمي (2–3 أيام)

#### 1.1 مصدر واحد للحقيقة
- إنشاء `constraints.txt` (pip) و`environment.lock` (نتيجة `micromamba list --explicit`).
- حذف `requirements.txt` أو جعله مُشتقاً آلياً من `constraints.txt`.
- `Dockerfile` يستخدم `-c constraints.txt` و`gromacs=2026.3=<build>`.

#### 1.2 إزالة ازدواج Vina
- حذف `autodock-vina` من `apt-get` (**لا يمسّ فصل البيئات**).
- `docking/storage.py`: تحقق SHA-256 قبل التشغيل، وتسجيل `vina_sha256` في `software_versions`.
- `tools/vina/vina.exe` (Windows، 1.2.7) يبقى **dev-only** مع تحذير في الواجهة عند اختلاف النسخة.

#### 1.3 تثبيت `verify_stack.py`
ملف `scripts/verify_stack.py` يطبع ويقارن:
```
gromacs, gmx_simd, vina+sha256, rdkit, openbabel, python,
acpype, ambertools, parmed, streamlit
```
يقارن بـ `expected_versions.json` ويرجع 1 عند أي اختلاف.

#### 1.4 تثبيت الصورة المنشورة
- `Launcher.cs:11` → `biobio-app:1.0.0` (لا `latest`).
- بعد السحب: مقارنة `RepoDigests` بـ digest مُضمَّن في الـ EXE.
- توفير `artifacts/biobio-1.0.0.tar` + `.sha256` للنشر بلا إنترنت.

**قبول المرحلة**: `verify_stack.py` ينجح داخل الصورة، وdigest الـ EXE == digest الصورة.

---

### المرحلة 2 — P1: طوبولوجيا حتمية (1–2 يوم)

#### 2.1 منطق الغلاف (أربع حالات)
```python
rc = run_acpype()                      # يُسجَّل كاملاً
have_amber = prmtop.exists() and inpcrd.exists()
have_gmx   = itp.exists() and gro.exists() and itp.stat().st_size > 0

if have_gmx and validate(itp, net_charge):     # ACPYPE نجح (أو انهار بعد الحفظ)
    ok = True
elif have_amber:                                # انهيار بعد ملفات Amber
    convert_amber_to_gmx(resname, net_charge)   # ParmEd + معايرة
    ok = validate(itp, net_charge)
else:
    ok = False
sys.exit(0 if ok else 1)
```
**ثلاثة قواعد**: لا كتابة فوق مخرجات ACPYPE السليمة؛ `free(): invalid pointer` غير قاتل إذا وُجد `prmtop`؛ معايرة شحنة إلزامية في مسار ParmEd.

#### 2.2 معايرة الشحنة في `convert_amber_to_gmx.py`
```python
total = sum(charges)
residual = total - float(net_charge)
lim = max(charges, key=abs)                       # ACPYPE: أقصى |q|
n_lim = sum(1 for q in charges if q == lim)
balanced = [q - residual / n_lim if q == lim else q for q in charges]
assert abs(sum(balanced) - float(net_charge)) < 1e-6
```
**تصحيح مقاس (2026-11-10)**: الكود المقترح أصلاً كان يضع الفارق على ذرة واحدة؛ قراءة `acpype/topol.py:1564` أثبتت أن `balanceCharges` يوزّعه بالتساوي على **كل** ذرة تحمل أكبر `|q|`. طُبّق سلوك ACPYPE حرفياً لضمان تطابق المسارين (ACPYPE الأصلي وParmEd البديل).
**يبقى** حذف `[ defaults ]` (تحقّق: `0.5 / 0.8333` متطابقة مع GAFF2).

#### 2.3 حارس في `md/topology.py`
بعد `_rename_itp_residue` (وليس داخل `_parse_itp_atoms`، حتى لا يخفي خطأ عمود الشحنة):
```python
def _enforce_ligand_charge(itp_path, net_charge):
    q_sum, n_atoms = _ligand_charge_sum(itp_path)   # قراءة مستقلة لعمود charge
    if abs(q_sum - float(net_charge)) > 1e-4:
        raise TopologyBuildError(...)
    return q_sum, n_atoms
```
تُعاد `charge_sum` و`charge_atoms` و`net_charge` في قاموس الليكند لتظهر في السجلات.

**قبول المرحلة**:
- `sum(q) == 0.000000` لـ ERLOT.
- `grompp` بلا تحذير net charge.
- اختبار سلبي: ليكند +1 مع `net_charge=0` ⇒ يفشل بخطأ صريح.

---

### المرحلة 3 — P1: النقل والتشغيل (1–2 يوم)

#### 3.1 `gmx` بلا symlink هش
```dockerfile
ln -sf /opt/conda/bin/gmx /usr/local/bin/gmx   # الموزّع (CPUID)
```
`entrypoint.sh` يسجّل:
```bash
echo "[BIOBIO] GROMACS SIMD: $("$GMX" --version 2>/dev/null | grep -i 'SIMD instructions' | cut -d: -f2)"
```
**قبول**: على جهازك `AVX2_256` (بلا خسارة)، وعلى جهاز بلا AVX2 يعمل بـ SSE2.

#### 3.2 إلغاء ترقيع acpype
- حذف `RUN printf ... >> acpype/utils.py`.
- `md/topology.py:1146` يقرأ `$AMBERHOME` أولاً ثم `bundled_amber_dir()`.
**قبول**: لا تعديل على site-packages؛ `acpype` يعمل.

#### 3.3 `.dockerignore` بدل `rm -rf`
```gitignore
.git
.venv
__pycache__/
*.pyc
build/
dist/
BIOBIO.exe
*.spec
launcher.py
native_window.py
Dockerfile
*.md
md_runs/
docking_runs/
tests/
artifacts/
```
**حُذف نهائياً** أي `rm -rf /app/...` بعد `COPY`.
**قبول**: `md/indexing.py` موجود؛ `python -c "import md.simulation"` ينجح.

#### 3.4 نافذة بلا Python
```csharp
string edge = @"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe";
string profile = Path.Combine(Environment.GetFolderPath(
    Environment.SpecialFolder.LocalApplicationData), "BIOBIO", "Viewer");
Process.Start(edge, $"--app={URL} --user-data-dir=\"{profile}\" " +
                    "--no-first-run --no-default-browser-check --window-size=1440,920");
```
مسار صريح ⇒ لا Brave؛ بروفايل معزول ⇒ لا اندماج مع Edge مفتوح.

#### 3.5 Preflight
`docker info` + `wsl --status` + مساحة ≥ 15 GB + ذاكرة ≥ 8 GB + فحص المنفذ 8501.
إضافة `--rm --label biobio=1` إلى `docker run`، و`HEALTHCHECK` في الصورة.

**قبول المرحلة**: تجربة أولى نظيفة على جهاز بلا Python ولا صورة محلية.

---

### المرحلة 4 — العلم: مطابقة الويكي (2–3 أيام)

#### 4.1 اختيار الصندوق تلقائياً
الويكي يوصي بـ dodecahedron لأنه يوفّر ~23% **للبروتينات الكروية** (`wiki/gromacs-protein-ligand-protocol.md:91`). هنا القياس يخالفه:

| الهندسة | الحجم |
|:--|:--|
| dodecahedron | 1049.08 nm³ |
| **triclinic** | **706.08 nm³ (−33%)** |
| cubic | 1483.62 nm³ |

**القاعدة**: احسب الحجمين من أبعاد الـ solute الفعلية، واختر الأصغر **فقط** إذا كان الفرق > 10%؛ وإلا ابقِ dodecahedron.
**تحذير إضافي**: إن كان `solute_diameter > 7 nm`، اعرض تنبيهاً عن الذيل C-terminal (977–995) مع خيار القص وتوثيق أثره على C-lobe.

#### 4.2 `replicates` الافتراضي = 1
`md/ui.py:381-390`: `value=3` → `value=1`، مع شرح أن 3 مطلوبة إحصائياً عند النشر.
(لا يتعارض مع قرار 0.5 ns.)

#### 4.3 توضيح دلالة 0.5 ns
- **وضع سريع**: 0.5 ns، مع نص: «صالح كمدخل MM-PBSA أو فحص أولي — ليس دليل استقرار».
  (مستند إلى `wiki/poli-2020-mm-pbsa-virtual-screening.md:30` الذي يحدّد 0.5–5 ns لاستخراج 10–100 لقطة MM-PBSA.)
- **وضع نشر**: preset ≥ 10 ns، replicates=3، مع الزمن المتوقع.

#### 4.4 تصحيح تسمية التحقق
المقاس: top-1 = 5.824 Å (فشل حسب `wiki/scoring-functions.md:80`)، rank-2 = 1.598 Å.

| الحالي | الصحيح |
|:--|:--|
| `validation passed (Top-10 recovery)` | `top-1 failed (5.82 Å); top-10 recovery passed (1.60 Å, rank 2)` |

#### 4.5 `md_manifest.json` دائماً
```
gromacs_version, gmx_simd, force_field, water_model,
ligand_charge_method, ligand_net_charge, ligand_charge_sum,
box_type, box_volume_nm3, solute_diameter_nm,
mdp_sha256{em,nvt,npt,md}, acpype_version, vina_version,
production_ns, n_replicates, analysis_reference, generated_at
```

**قبول المرحلة**: أي تشغيل ينتج manifest كاملاً؛ `box_volume ≤ 750 nm³` لهذا النظام.

---

### المرحلة 5 — بوّابة القبول والنشر (1 يوم)

`tests/test_packaging_contract.py` يعمل داخل الصورة ويفشل البناء:

| # | الفحص | الشرط |
|:--|:--|:--|
| 1 | `gmx --version` | `2026.3` + سطر SIMD |
| 2 | `vina --version` + sha256 | `AutoDock Vina 1.2.5` + مطابقة |
| 3 | `import openbabel` + `obabel -V` | كلاهما يعمل (تثبيت واحد) |
| 4 | `antechamber`, `tleap` | يعملان من بيئة `amber` |
| 5 | شحنة ليكند مرجعي | `abs(sum - net) < 1e-6` |
| 6 | `grompp` على نظام صغير | **0 WARNING** |
| 7 | `amber99sb-ildn` defaults | `0.5 / 0.8333` |
| 8 | `md/indexing.py` | `import md.simulation` ينجح |
| 9 | EM 10 خطوات | يكتمل |
| 10 | بصمة RDKit | تطابق `tests/golden/rdkit_fingerprint.json` |

`scripts/build.ps1`: بناء → بوّابة → وسم digest → `docker save` → `.sha256`.

---

## 5. جدول تثبيت الإصدارات

| المكوّن | القيمة المثبّتة | المصدر | ملاحظة |
|:--|:--|:--|:--|
| الصورة الأساس | `mambaorg/micromamba:2.10.0@sha256:e1476182…c1198` | بدل `latest` | **تصحيح مقاس (2026-11-10)**: الصورة العاملة فعلاً هي 2.10.0 على Debian 13؛ 1.5.10 كانت Debian 12 وتغيّر ABI. يُسجَّل الـ digest |
| Python | `3.11.18` | Docker الحالي | — |
| GROMACS | `2026.3` + build string | `nvt.log` | مصدر النتائج الحالية |
| RDKit | `2026.09.1` | manifests | **قرار §1.1** |
| OpenBabel | **wheel فقط** `openbabel-wheel==3.1.1.23` (يُبلّغ 3.1.0) | manifests | **تصحيح مقاس**: حزمة conda `openbabel` تتشارك الملفات نفسها مع الـ wheel؛ إزالة أحدهما بعد تثبيت الآخر تُفسد الاستيراد. الـ wheel وحده يعمل ويطابق خط الأساس |
| ParmEd | إصدار conda الحالي | — | موحّد في البيئتين |
| AmberTools | الإصدار الحالي + build | — | — |
| ACPYPE | `2023.10.27` | `acpype.log` | ما أنتج النتائج |
| Vina | `1.2.5` + SHA-256 | الملف المرفق | — |
| Streamlit | `1.62.0` | `requirements.txt` | — |
| NumPy | `>=1.26,<2.3` | `requirements.txt` | يمنع كسر ABI |

---

## 6. ما يُعدَّل vs ما لا يُعدَّل

### 6.1 ملفات ستُعدَّل

| الملف | التغيير |
|:--|:--|
| `Dockerfile` | pin كل شيء، إزالة apt vina، إزالة ترقيع acpype، موزّع gmx، `HEALTHCHECK` |
| `.dockerignore` | allowlist بدل `rm -rf` |
| `entrypoint.sh` | تسجيل SIMD + إمكانية تجاوز الخيوط |
| `acpype_wrapper.py` | منطق أربع حالات + تسجيل + رموز خروج صريحة |
| `convert_amber_to_gmx.py` | `main()` + `argparse` + معايرة شحنة + تحقق |
| `md/topology.py` | `$AMBERHOME` أولاً + حارس الشحنة |
| `md/simulation.py` | اختيار صندوق تلقائي + manifest دائم + توحيد `prod_duration_ns` |
| `md/ui.py` | `replicates=1` + presets + تحذير الذيل + تسمية التحقق |
| `docking/storage.py` | SHA-256 لـ Vina + إصدارات كاملة |
| `docking/ui.py` | إظهار top-1 vs top-N |
| `Launcher.cs` + `BIOBIO.exe` | وسم ثابت + digest + msedge معزول + preflight + `--rm` |
| `requirements.txt` | يُستبدل بـ `constraints.txt` |
| `app.py` | شريط إصدارات (اختياري) |

### 6.2 ملفات جديدة

```
PLAN.md
constraints.txt
environment.lock
expected_versions.json
scripts/build.ps1
scripts/verify_stack.py
scripts/golden_rdkit.py
scripts/release.ps1
tests/test_packaging_contract.py
tests/golden/rdkit_fingerprint.json
artifacts/baseline.json
CHANGELOG.md
SECURITY.md
```

### 6.3 ملفات **لن تُلمس** (مُجمَّدة)

| العنصر | السبب |
|:--|:--|
| `wiki/**` | مرجع علمي — للقراءة فقط |
| `raw/**` | مصادر غير قابلة للتعديل (قاعدة المشروع) |
| `md/engine.py` | طبقة تنفيذ GROMACS ناجحة |
| `md/indexing.py` | حُذف سابقاً فتعطّل المشروع (#10) |
| `analysis/**` | خارج نطاق التحويل |
| `docking/engine.py`, `ligand_prep.py`, `receptor_prep.py`, `openbabel_io.py`, `validation.py` | منطق علمي متحقَّق |
| `tools/vina/vina`, `tools/vina/vina.exe` | بصمة مرجعية |
| قيم `.mdp` (dt, cutoffs, PME, LINPCS, tcoupl, pcoupl, tau, POSRES, DispCorr) | مطابقة 100% للويكي |
| بنية البيئتين (base / amber) | حلّت تعارض conda |
| `ENTRYPOINT ["/entrypoint.sh"]` | حلّ `GMXRC.bash: line 85` |
| `csc.exe` بدل PyInstaller | حلّ `%TEMP%` |
| `IsAdmin()` + `runas` | حلّ صلاحيات Docker |
| shebang الغلاف `#!/opt/conda/envs/amber/bin/python` | حلّ بايثون الفرعي |
| حذف `[ defaults ]` من ITP الليكند | `0.5/0.8333` متطابقة |
| حذف `-v` من `_mdrun` | قرار السرعة |

---

## 7. سجل المخاطر

| المخاطرة | الاحتمال | الأثر | التخفيف |
|:--|:--|:--|:--|
| تثبيت إصدار يكسر حلاً قديماً | متوسط | مرتفع | مرحلة 0 + بوّابة + أرشيف رجوع |
| `free(): invalid pointer` يتكرر بشكل قاتل | متوسط | متوسط | ParmEd + معايرة + حارس |
| CPU بلا AVX2 عند صديق | منخفض | متوسط | موزّع gmx |
| انقطاع شبكة أثناء `pull` | متوسط | منخفض | مهلة + tar بديل |
| تثبيت RDKit يغيّر conformer | مرتفع | منخفض | بصمة ذهبية + إعادة تحقق |
| صور Docker متضخمة (>8 GB) | متوسط | منخفض | `.dockerignore` + `--no-cache` |
| نتائج قديمة لا تُقارن بالجديدة | مرتفع | متوسط | `schema_version` + `generated_at` + الإصدارات في manifest |

---

## 8. الجدول الزمني

| المرحلة | الزمن | يعتمد على | مخرج قابل للتسليم |
|:--|:--|:--|:--|
| 0 — الأساس | 0.5 يوم | — | `baseline.json` + أرشيف |
| 1 — P0 تثبيت | 2–3 أيام | 0 | `constraints.txt` + `verify_stack` |
| 2 — طوبولوجيا | 1–2 يوم | 1 | غلاف حتمي + حارس شحنة |
| 3 — نقل/تشغيل | 1–2 يوم | 1 | صورة تعمل على أي CPU + EXE بلا Python |
| 4 — علم | 2–3 أيام | 2 | صندوق أمثل + manifest + تسمية صحيحة |
| 5 — بوابة/نشر | 1 يوم | 1–4 | `biobio-1.0.0.tar` + `.sha256` |
| إعادة تحقق علمي | 1 يوم | 5 | redocking + 1 ns MD موثّقان |
| **الإجمالي** | **12–16 يوم عمل** | | |

---

## 9. قائمة تحقق الإصدار (Release Checklist)

1. [ ] كل الاختبارات في §المرحلة 5 خضراء.
2. [ ] `verify_stack.py` يطابق `expected_versions.json`.
3. [ ] بصمة RDKit مطابقة للملف الذهبي.
4. [ ] redocking مرجعي: top-1 أو top-10 كما هو موثّق، بلا تغيير في التسمية.
5. [ ] `grompp` بلا تحذير شحنة.
6. [ ] `docker save` + `.sha256` + digest مضمّن في الـ EXE.
7. [ ] `CHANGELOG.md` يوثّق كل تغيير.
8. [ ] تجربة نظيفة على جهاز بلا Docker image ولا Python.
9. [ ] `README-USER.md` بالعربية: تثبيت Docker Desktop، المساحة، الزمن، حدود النتائج.

---

## 10. خارج النطاق / فجوات معروفة

- **GPU**: لا تسريع CUDA في هذه النسخة؛ يبقى تحسينه بنداً مستقبلاً منفصلاً.
- **HMR / خطوة 4 fs**: يحتاج إعادة معاملة طوبولوجيا وتحققاً مستقلاً — خارج 1.0.
- **CHARMM36/CGenFF**: الويكي يفضّله، لكن المشروع مبني على AMBER/GAFF2؛ التبديل قرار علمي مستقل.
- **حساسية النتائج لتغيير الإصدارات**: أي تغيير لـ GROMACS/RDKit/Vina يستلزم إعادة التحقق كاملاً؛ الخطة تفرض ذلك ولا تستثني منه.
- **حدود علمية ثابتة**: 0.5 ns لا تُثبت استقرار ارتباط (الويكي: أحداث ارتباط تحتاج µs–ms)؛ درجات Vina فرضيات ترتيب لا طاقات ارتباط.

---

## 11. معايير «تم»

| البند | الشرط |
|:--|:--|
| التكرار | بناء من الصفر يعطي نفس الأرقام (بصمات) على نفس الإصدارات |
| الحتمية | نفس المدخلات + نفس البذور ⇒ نفس الشحنة والصندوق والـ manifest |
| السلامة العلمية | صفر تحذير PME؛ شحنة 0.000000؛ معاملات mdp مطابقة للويكي |
| النقل | يعمل على CPU بلا AVX2، وعلى جهاز بلا Python، وبلا إنترنت (بعد التحميل) |
| الشفافية | كل نتيجة تحمل إصداراتها وبصماتها وحدودها |
| الرجوع | أي إخفاق يُعاد إلى `artifacts/biobio-baseline.tar` |

---

## 12. الخطوة التالية

البدء بالمرحلة 0 ثم 1. لن يُعدَّل أي ملف قبل تأكيد صريح على هذه النسخة من `PLAN.md`.
