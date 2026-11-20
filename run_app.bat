@echo off
setlocal
title Bio-Impact Analyzer
echo ==========================================================
echo Starting Bio-Impact Analyzer...
echo ==========================================================
cd /d "%~dp0"

:: 1. Prefer Windows System Python 3.12 (Direct laptop environment)
py -3.12 -c "import streamlit, rdkit, Bio, openbabel" >nul 2>&1
if not errorlevel 1 (
    echo [OK] Using Windows Python 3.12 on this laptop...
    py -3.12 -m streamlit run app.py --server.address 127.0.0.1
    goto :done
)

:: 2. Try default system python
python -c "import streamlit, rdkit, Bio, openbabel" >nul 2>&1
if not errorlevel 1 (
    echo [OK] Using system Python on this laptop...
    python -m streamlit run app.py --server.address 127.0.0.1
    goto :done
)

:: 3. Try virtual environment if present
if exist ".venv\Scripts\python.exe" (
    echo [OK] Using virtual environment .venv...
    set "PATH=%~dp0.venv\Scripts;%PATH%"
    ".venv\Scripts\python.exe" -m streamlit run app.py --server.address 127.0.0.1
    goto :done
)

echo [ERROR] Required libraries not found. Installing into Python...
py -3.12 -m pip install -r requirements.txt
if not errorlevel 1 (
    py -3.12 -m streamlit run app.py --server.address 127.0.0.1
    goto :done
)

pause
:done

