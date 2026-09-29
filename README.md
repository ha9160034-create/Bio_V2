# 🧬 Bio-Impact Analyzer

A comprehensive computational platform for structural bioinformatics and molecular docking simulations.

## ✨ المميزات الرئيسية (Key Features)

1. **التحليل الهيكلي والمقارنة (Structural Analysis & Comparison)**:
   - مقارنة السلاسل والمحاذاة التتابعية (Sequence Alignment).
   - حساب مساحة السطح المتاح للمذيب (SASA Profiling).
   - عرض ثلاثي الأبعاد تفاعلي للبروتينات باستخدام `py3Dmol`.

2. **محاكاة الارتباط الجزيئي (Molecular Docking Simulation)**:
   - محاكاة الارتباط باستخدام محرك **AutoDock Vina**.
   - تحديد وتعيين جيوب الارتباط (Pocket & Binding Site Mapping).
   - تحليل التفاعلات غير التساهمية وتوليد التقارير.

---

## 🚀 كيفية التثبيت والتشغيل (Installation & Quickstart)

### المتطلبات (Prerequisites)
- Python 3.10+
- Git

### خطوات التشغيل (Run Locally)

1. **استنساخ المستودع (Clone repository):**
   ```bash
   git clone https://github.com/<username>/<repo-name>.git
   cd Bio_V2
   ```

2. **إنشاء وتفعيل البيئة الافتراضية (Create virtual environment):**
   ```bash
   python -m venv venv
   # On Windows:
   venv\Scripts\activate
   # On Linux/macOS:
   source venv/bin/activate
   ```

3. **تثبيت الحزم المطلوبة (Install dependencies):**
   ```bash
   pip install -r requirements.txt
   ```

4. **تشغيل التطبيق (Launch application):**
   ```bash
   streamlit run b.py
   ```

---

## 📁 هيكلية المشروع (Project Structure)

```text
Bio_V2/
├── analysis/         # وحدات التحليل الهيكلي والحسابي (Structure & SASA analysis)
├── docking/          # وحدات محاكاة الارتباط ومحرك Vina وتوليد التفاعلات
├── docking_runs/     # مخرجات عمليات المحاكاة
├── data/             # ملفات البيانات والجينات (Data & mutations)
├── tests/            # الاختبارات الآلية (Automated test suite)
├── tools/            # أدوات خارجية مدمجة (Vina binaries)
├── b.py              # نقطة الدخول الرئيسية للتطبيق (Main Streamlit App)
└── requirements.txt  # الاعتماديات (Dependencies)
```

---

## 🧪 الاختبارات (Running Tests)

لتشغيل الاختبارات والتأكد من سلامة النظام:
```bash
pytest tests/
```
