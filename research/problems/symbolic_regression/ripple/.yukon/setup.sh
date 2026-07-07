#!/usr/bin/env bash
# Per-leaf setup for the FrontierCS **symbolic_regression** research family. The evaluator.py imports
# numpy, pandas, and sympy — install them into the setup snapshot so run.sh can score. Datasets are
# already bundled in resources/data (no download step). Kept idempotent + quiet; fails closed if pip
# can't provide the deps (run.sh would otherwise crash on import).
set -euo pipefail

# Prefer an already-satisfied environment (snapshot reuse) — only install if an import is missing.
if python3 -c "import numpy, pandas, sympy" >/dev/null 2>&1; then
  echo "setup.sh: numpy/pandas/sympy already present" >&2
  exit 0
fi

echo "setup.sh: installing research eval deps (numpy, pandas, sympy)" >&2
python3 -m pip install --quiet --disable-pip-version-check numpy pandas sympy >&2

# Verify the deps import (fail closed rather than defer a cryptic crash to run.sh).
python3 -c "import numpy, pandas, sympy" || { echo "setup.sh: research eval deps failed to import after install" >&2; exit 1; }
echo "setup.sh: research eval deps ready" >&2
