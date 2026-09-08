#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC_DIR="${SCRIPT_DIR}/fortran"
BUILD_DIR="${SCRIPT_DIR}/build"
FC="${FC:-gfortran}"

mkdir -p "${BUILD_DIR}"

COMMON_FLAGS=(-O3 -march=native -funroll-loops -ffast-math -fPIC)

build_linux() {
  local source_file="$1"
  local output_file="$2"
  local -a lapack_flags=(-llapack)
  local -a blas_flags=(-lblas)

  if command -v R >/dev/null 2>&1; then
    read -r -a lapack_flags <<< "$(R CMD config LAPACK_LIBS)"
    read -r -a blas_flags <<< "$(R CMD config BLAS_LIBS)"
  fi

  "${FC}" -shared "${COMMON_FLAGS[@]}" -o "${output_file}" "${source_file}" \
    "${lapack_flags[@]}" "${blas_flags[@]}"
}

build_macos() {
  local source_file="$1"
  local output_file="$2"
  local object_file="${output_file%.so}.o"
  local gcc_prefix
  gcc_prefix="$(brew --prefix gcc)"
  "${FC}" -c "${COMMON_FLAGS[@]}" "${source_file}" -o "${object_file}"
  xcrun clang -dynamiclib -undefined dynamic_lookup "${object_file}" \
    -o "${output_file}" -framework Accelerate \
    -L"${gcc_prefix}/lib/gcc/current" -lgfortran -lquadmath
}

build_one() {
  local source_name="$1"
  local output_name="$2"
  if [[ "$(uname -s)" == "Darwin" ]]; then
    build_macos "${SRC_DIR}/${source_name}" "${BUILD_DIR}/${output_name}"
  else
    build_linux "${SRC_DIR}/${source_name}" "${BUILD_DIR}/${output_name}"
  fi
}

build_one llqr_ppro_fast_v41_active_first.f90 llqr_v41.so
build_one llqr_ppro_fast_v42_active_index.f90 llqr_v42.so
build_one llqr_seq_lean_sortskip.f90 llqr_lean_seq.so
build_one tvcqr_ppro_v41_active_first.f90 tvcqr_v41.so
build_one tvcqr_ppro_v42_contract_fast.f90 tvcqr_v42.so
build_one tvcqr_seq_lean_nohistory.f90 tvcqr_lean_seq.so

for binary in llqr_v41.so llqr_v42.so llqr_lean_seq.so tvcqr_v41.so tvcqr_v42.so tvcqr_lean_seq.so; do
  test -s "${BUILD_DIR}/${binary}"
done

printf 'Built six experiment binaries in %s\n' "${BUILD_DIR}"
