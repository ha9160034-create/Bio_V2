@echo off
setlocal
title Bio-Impact Analyzer
echo ==========================================================
echo Starting Bio-Impact Analyzer...
echo ==========================================================
cd /d "%~dp0"
if not exist ".venv\Scripts\python.exe" (
    echo Project environment is missing. Run setup_app.bat first.
    pause
    exit /b 1
)
set "PATH=%~dp0.venv\Scripts;%PATH%"
".venv\Scripts\python.exe" -m streamlit run app.py --server.address 127.0.0.1
pause
