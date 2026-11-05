#!/usr/bin/env bash
# Install the MD runtime, without choosing simulation parameters.
# References: wiki/gromacs-protein-ligand-protocol.md, wiki/force-fields.md,
# wiki/open-force-field-and-ligand-parameterization.md.
# Package sources: packages.ubuntu.com and pypi.org/project/acpype/2026.10.15/.
set -euo pipefail
if [ "$(id -u)" -ne 0 ]; then
    echo 'Run this installation in the WSL distribution as root.' >&2
    exit 1
fi

export DEBIAN_FRONTEND=noninteractive
apt-get update
apt-get install -y --no-install-recommends \
    gromacs gromacs-data openbabel python3-venv ca-certificates \
    libgfortran5 libstdc++6 libgomp1 libblas3 liblapack3 libcurl4t64 \
    libxrender1 libxext6 libsm6

MD_VENV=/opt/bio-v2-md
if [ -e "$MD_VENV" ] && [ ! -f "$MD_VENV/pyvenv.cfg" ]; then
    echo "$MD_VENV exists and is not a Python environment; refusing to overwrite it." >&2
    exit 1
fi
if [ ! -x "$MD_VENV/bin/python" ]; then
    python3 -m venv "$MD_VENV"
fi
"$MD_VENV/bin/python" -m pip install --disable-pip-version-check 'acpype==2026.10.15'

md_link_if_safe() {
    local md_source="$1" md_destination="$2"
    if [ -e "$md_destination" ] || [ -L "$md_destination" ]; then
        if [ "$(readlink -f "$md_destination")" != "$(readlink -f "$md_source")" ]; then
            echo "Refusing to replace an existing tool: $md_destination" >&2
            return 1
        fi
    else
        ln -s "$md_source" "$md_destination"
    fi
}
md_link_if_safe "$MD_VENV/bin/acpype" /usr/local/bin/acpype
MD_AMBER_DIR=$("$MD_VENV/bin/python" -c 'from acpype.utils import bundled_amber_dir; print(bundled_amber_dir() or "")')
if [ -z "$MD_AMBER_DIR" ]; then
    echo 'ACPYPE installed without bundled AmberTools; actual MD parameterization is not ready.' >&2
    exit 1
fi
for MD_TOOL in antechamber parmchk2 tleap sqm; do
    if [ ! -x "$MD_AMBER_DIR/bin/$MD_TOOL" ]; then
        echo "Missing AmberTools executable: $MD_TOOL" >&2
        exit 1
    fi
    md_link_if_safe "$MD_AMBER_DIR/bin/$MD_TOOL" "/usr/local/bin/$MD_TOOL"
done

gmx --version
acpype --version
command -v gmx acpype obabel antechamber parmchk2 tleap sqm
echo 'MD_BACKEND_PACKAGES_INSTALLED; real simulation verification is still required.'
