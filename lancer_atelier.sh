#!/bin/sh
# G2Egrid workbench / atelier G2Egrid : local server and browser (macOS, Linux).
# Requires uv / il faut uv : https://docs.astral.sh/uv/
cd "$(dirname "$0")" || exit 1
if ! command -v uv >/dev/null 2>&1; then
  echo "uv not found / uv est introuvable : https://docs.astral.sh/uv/getting-started/installation/"
  exit 1
fi
exec uv run --extra api python -m g2egrid_api "$@"
