# BIOBIO runtime image — pinned on purpose.
#
# Design references:
#   * PLAN.md §1.1, §1.2, §5 (single source of truth for versions)
#   * wiki/gromacs-protein-ligand-protocol.md (GROMACS/AMBER setup)
#   * wiki/open-force-field-and-ligand-parameterization.md (ACPYPE/GAFF2)
#
# Rules:
#   * No `latest`, no unpinned conda/pip installs.
#   * Do not install autodock-vina from apt: the app ships a SHA-256-pinned
#     binary at /app/tools/vina/vina and docking/storage.py verifies it.
#   * Do not install conda `openbabel`: it collides file-for-file with
#     openbabel-wheel and the verified baseline (3.1.0) is the wheel.

FROM mambaorg/micromamba:2.10.0@sha256:e1476182cd40eb1f7c0cfae4adec9bfd78f3d7f05c61368306abfdc8939c1198

USER root

RUN apt-get update && apt-get install -y --no-install-recommends \
    procps \
    libxrender1 \
    libxext6 \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/*

# Base environment: MD engine + topology tooling.
RUN micromamba install -y -n base -c conda-forge \
    python=3.11.18 \
    gromacs=2026.3=nompi_h26635d9_100 \
    parmed=4.3.1=py311h74232ed_1 \
    numpy=2.2.6 \
    pandas=3.1.1 \
    && micromamba clean --all --yes

# Amber environment: AM1-BCC charges (antechamber/sqm) + GAFF2 topology.
RUN micromamba create -y -n amber -c conda-forge \
    python=3.11.18 \
    ambertools=26.0=cuda_None_nompi_py311h13e84fe_101 \
    acpype=2023.10.27=pyhd8ed1ab_0 \
    parmed=4.3.1=py311h74232ed_1 \
    numpy=2.2.6 \
    && micromamba clean --all --yes

ENV AMBERHOME="/opt/conda/envs/amber"
ENV PATH="/usr/local/bin:/opt/conda/bin.AVX2_256:/opt/conda/bin:/opt/conda/envs/amber/bin:/usr/bin:/bin:${PATH}"

# Python stack: versions are fixed in constraints.txt (single source of truth).
COPY constraints.txt requirements.txt /tmp/lock/
RUN pip install --no-cache-dir --root-user-action=ignore \
    -c /tmp/lock/constraints.txt \
    -r /tmp/lock/requirements.txt

WORKDIR /app
COPY . /app

RUN tr -d '\r' < /app/convert_amber_to_gmx.py > /usr/local/bin/convert_amber_to_gmx.py && \
    tr -d '\r' < /app/acpype_wrapper.py > /usr/local/bin/acpype && \
    tr -d '\r' < /app/acpype_wrapper.py > /opt/conda/envs/amber/bin/acpype && \
    tr -d '\r' < /app/entrypoint.sh > /entrypoint.sh && \
    chmod +x /usr/local/bin/convert_amber_to_gmx.py /usr/local/bin/acpype /opt/conda/envs/amber/bin/acpype /entrypoint.sh && \
    ln -sf /opt/conda/bin.AVX2_256/gmx /usr/local/bin/gmx

# Kept for Phase 3 replacement by .dockerignore (see PLAN.md #10).
RUN rm -rf /app/tests /app/Dockerfile /app/Launcher.cs /app/launcher.py /app/native_window.py /app/convert_amber_to_gmx.py /app/acpype_wrapper.py /app/entrypoint.sh && \
    mkdir -p /app/docking_runs /app/md_runs && \
    if [ -f /app/tools/vina/vina ]; then chmod +x /app/tools/vina/vina; fi

EXPOSE 8501
ENTRYPOINT ["/entrypoint.sh"]
