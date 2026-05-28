#!/usr/bin/env bash
# Generate the gFlex logo from scratch.
#
# Usage:
#   bash docs/logo/make_logo.sh          # full pipeline
#   bash docs/logo/make_logo.sh --fast   # skip data generation (reuse /tmp cache)
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO"

FAST=0
for arg in "$@"; do [[ "$arg" == "--fast" ]] && FAST=1; done

if [[ $FAST -eq 0 ]]; then
    echo "── Step 1: generate flexure data ──────────────────────────"
    python docs/logo/generate_logo_data.py
else
    echo "── Step 1: skipped (--fast) ───────────────────────────────"
fi

echo "── Step 2: render 3-D scene (Blender) ─────────────────────"
blender --background --python docs/logo/blender_logo.py

echo "── Step 3: composite text ──────────────────────────────────"
python docs/logo/add_logo_text.py

echo "────────────────────────────────────────────────────────────"
echo "Done.  Output: docs/_static/logo.png"
ls -lh docs/_static/logo.png
