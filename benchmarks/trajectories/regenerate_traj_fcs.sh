#!/usr/bin/env bash
# Rebuild the eight-site counting benchmark from available trajectory batches.
# If the large jump records are absent, the checked-in estimates under est/ are used.
# Usage: bash benchmarks/trajectories/regenerate_traj_fcs.sh [--validate]
set -e
cd "$(dirname "$0")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}" OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}" MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
if [ -n "${QENV:-}" ] && [ -x "${QENV}/bin/python" ]; then PYTHON_BIN="${QENV}/bin/python"; fi
mkdir -p est
for run in gpu_U0.04_c10 gpu_U0.08_c12 gpu_U0.10_c14; do
  n=$(ls "runs/$run" 2>/dev/null | grep -c '^\(cpu\)\?batch_[0-9]*\.npz$' || true)
  if [ "$n" -gt 0 ]; then "$PYTHON_BIN" traj_fcs_estimator.py "runs/$run" --out "est/est_$run.json"; fi
done
if [ ! -f est/est_gpu_U0.10_cut8.json ] && [ -d runs/gpu_U0.10 ]; then
  "$PYTHON_BIN" traj_fcs_estimator.py runs/gpu_U0.10 --out est/est_gpu_U0.10_cut8.json
fi
if [ "${1:-}" = "--validate" ] || [ ! -f traj_validation.json ]; then "$PYTHON_BIN" traj_validate.py; fi
"$PYTHON_BIN" traj_fcs_final.py
ls -l traj_fcs_results.json traj_validation.json chain8_contour.json fig_fcs.pdf
