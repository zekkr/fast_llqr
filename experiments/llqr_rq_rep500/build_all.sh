#!/usr/bin/env bash
set -euo pipefail

project_root="${SSQR_PROJECT_ROOT:-$(pwd)}"
experiment_dir="${SSQR_EXPERIMENT_DIR:-experiments/llqr_rq_rep500}"
cd "$project_root"

fc="${FC:-gfortran}"
source_dir="reproduction/frozen/fortran"
build_dir="$experiment_dir/build"
mkdir -p "$build_dir"

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1

read -r -a lapack_flags <<< "$(R CMD config LAPACK_LIBS)"
read -r -a blas_flags <<< "$(R CMD config BLAS_LIBS)"

compile_shared() {
  local output="$1"
  shift
  if [[ "$(uname -s)" == Darwin ]]; then
    local sdk
    sdk="$(xcrun --show-sdk-path)"
    "$fc" -shared -fPIC -ffree-line-length-none "$@" \
      -Wl,-syslibroot,"$sdk" "${lapack_flags[@]}" "${blas_flags[@]}" -o "$output"
  else
    "$fc" -shared -fPIC -ffree-line-length-none "$@" \
      "${lapack_flags[@]}" "${blas_flags[@]}" -o "$output"
  fi
}

compile_shared "$build_dir/ssqr_checked.so" \
  -O0 -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow \
  "$source_dir/weighted_qr_core.f90" "$source_dir/kernel_entry.f90"

optimized=(-O3 -march=native -funroll-loops -ffast-math)
compile_shared "$build_dir/ssqr_optimized.so" "${optimized[@]}" \
  "$source_dir/weighted_qr_core.f90" "$source_dir/kernel_entry.f90"
compile_shared "$build_dir/llqr_seq_lean_sortskip.so" "${optimized[@]}" \
  "$source_dir/llqr_seq_lean_sortskip.f90"
compile_shared "$build_dir/tvcqr_seq_lean_nohistory.so" "${optimized[@]}" \
  "$source_dir/tvcqr_seq_lean_nohistory.f90"

"$fc" --version > "$build_dir/compiler_version.txt"
printf '%s\n' "${optimized[*]}" > "$build_dir/compiler_flags.txt"
uname -a > "$build_dir/build_uname.txt"
if command -v lscpu >/dev/null 2>&1; then lscpu > "$build_dir/build_lscpu.txt"; fi
if command -v sha256sum >/dev/null 2>&1; then
  sha256sum "$source_dir"/*.f90 "$experiment_dir"/*.R "$experiment_dir"/*.sh \
    > "$build_dir/source_sha256.txt"
else
  shasum -a 256 "$source_dir"/*.f90 "$experiment_dir"/*.R "$experiment_dir"/*.sh \
    > "$build_dir/source_sha256.txt"
fi
