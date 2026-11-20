@echo off
setlocal
cd /d "%~dp0"
echo ==========================================================
echo Installing Bio_V2 dependencies on this laptop...
echo ==========================================================

:: Install to Windows System Python 3.12
py -3.12 -m pip install --upgrade pip
py -3.12 -m pip install -r requirements.txt
if not errorlevel 1 (
    echo [OK] Setup complete! You can now start the app with run_app.bat.
    pause
    exit /b 0
)

:: Fallback to python
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
if not errorlevel 1 (
    echo [OK] Setup complete! You can now start the app with run_app.bat.
    pause
    exit /b 0
)

echo [ERROR] Installation failed. Please ensure Python 3.12 is installed.
pause
exit /b 1

