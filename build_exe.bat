@echo off
rem Optional: build a single SGCalibrator.exe (needs install.bat done first).
cd /d "%~dp0"
py -3 -m pip install pyinstaller
py -3 -m PyInstaller --onefile --windowed --name SGCalibrator sgcal_gui.py
echo The exe is in the dist folder.
pause
