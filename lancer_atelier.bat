@echo off
rem G2Egrid workbench / atelier G2Egrid : local server and browser.
rem Requires uv / il faut uv : https://docs.astral.sh/uv/
rem The first run installs the environment ; the next ones start in seconds.
cd /d "%~dp0"
where uv >nul 2>nul
if errorlevel 1 (
  echo uv not found / uv est introuvable : https://docs.astral.sh/uv/getting-started/installation/
  pause
  exit /b 1
)
uv run --extra api python -m g2egrid_api %*
if errorlevel 1 pause
