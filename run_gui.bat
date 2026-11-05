@echo off
cd /d "%~dp0"
py -3 sgcal_gui.py
if errorlevel 1 pause
