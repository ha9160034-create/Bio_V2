@echo off
title Push Bio_V2 to GitHub
echo ==========================================================
echo  Pushing latest update to Bio_V2...
echo ==========================================================
cd /d "%~dp0"
git push origin main
echo.
echo ==========================================================
if errorlevel 1 (
    echo [!] Push failed or cancelled. Please check your credentials.
) else (
    echo [OK] Successfully updated https://github.com/ha9160034-create/Bio_V2 !
)
echo ==========================================================
pause
