@echo off
rem One-time setup: installs the Python packages SG Calibrator needs.
where py >nul 2>nul
if errorlevel 1 (
  echo Python 3 was not found. Install it from https://www.python.org/downloads/ ^(tick "Add python.exe to PATH" and "tcl/tk"^), then run this again.
  pause
  exit /b 1
)
py -3 -m pip install --upgrade -r "%~dp0requirements.txt"
echo.
echo Done. Double-click run_gui.bat to start.
pause
