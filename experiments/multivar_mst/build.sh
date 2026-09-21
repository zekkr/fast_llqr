#!/usr/bin/env bash
set -euo pipefail

project_root="${SSQR_PROJECT_ROOT:-$(pwd)}"
experiment_dir="experiments/multivar_mst"
build_dir="$project_root/$experiment_dir/build"
mkdir -p "$build_dir"
cd "$project_root"

fc="${FC:-gfortran}"
read -r -a lapack_flags <<< "$(R CMD config LAPACK_LIBS)"
read -r -a blas_flags <<< "$(R CMD config BLAS_LIBS)"
common=(-shared -fPIC -ffree-line-length-none)
if [[ "$(uname -s)" == Darwin ]]; then
  sdk="$(xcrun --show-sdk-path)"
  common+=("-Wl,-syslibroot,$sdk")
fi
sources=("$experiment_dir/fortran/weighted_qr_tree_core.f90" "$experiment_dir/fortran/tree_entry.f90")
"$fc" "${common[@]}" -O0 -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow \
  "${sources[@]}" "${lapack_flags[@]}" "${blas_flags[@]}" -o "$build_dir/multivar_mst_checked.so"
"$fc" "${common[@]}" -O3 -march=native -funroll-loops -ffast-math \
  "${sources[@]}" "${lapack_flags[@]}" "${blas_flags[@]}" -o "$build_dir/multivar_mst_optimized.so"
"$fc" --version > "$build_dir/compiler_version.txt"
printf '%s\n' '-O3 -march=native -funroll-loops -ffast-math' > "$build_dir/compiler_flags.txt"
if command -v sha256sum >/dev/null 2>&1; then
  sha256sum "${sources[@]}" "$experiment_dir"/*.R "$experiment_dir"/*.sh > "$build_dir/source_sha256.txt"
else
  shasum -a 256 "${sources[@]}" "$experiment_dir"/*.R "$experiment_dir"/*.sh > "$build_dir/source_sha256.txt"
fi
