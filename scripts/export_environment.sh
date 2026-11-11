#!/usr/bin/env bash
# Export the exact conda/pip solution of a running BIOBIO image.
#
# Usage (from the repository root):
#   docker run --rm --entrypoint /bin/bash biobio-app:1.0.0 /app/scripts/export_environment.sh
#
# The script prints a lock file on stdout; redirect it to environment.lock.
set -euo pipefail

echo "# environment.lock - exact conda/pip solution of the verified BIOBIO image"
echo "# Generated inside: $(python -V 2>&1)"
echo "#"
echo "# BASE ENVIRONMENT - micromamba list --explicit"
micromamba list -n base --explicit
echo ""
echo "# AMBER ENVIRONMENT - micromamba list --explicit"
micromamba list -n amber --explicit
echo ""
echo "# PYTHON (pip) - pip freeze"
python -m pip freeze
