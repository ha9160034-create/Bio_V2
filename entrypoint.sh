#!/bin/bash
export AMBERHOME="/opt/conda/envs/amber"
export PATH="/usr/local/bin:/opt/conda/bin.AVX2_256:/opt/conda/bin:/opt/conda/envs/amber/bin:/usr/bin:/bin:$PATH"
export OMP_NUM_THREADS=$(nproc)
exec /opt/conda/bin/streamlit run /app/app.py --server.port=8501 --server.address=0.0.0.0 --server.headless=true --client.toolbarMode=minimal
