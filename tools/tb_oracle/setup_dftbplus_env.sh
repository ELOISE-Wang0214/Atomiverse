#!/usr/bin/env bash
# Source this file to put DFTB+ 24.1 + SK paths on PATH for offline pins.
# Install lives outside the git tree (see software/README.md next to molcrafts).
#
#   source tools/tb_oracle/setup_dftbplus_env.sh
#   dftb+ --version
#
# Override with ATV_DFTBPLUS_ROOT / ATV_DFTB_PARAMS_ROOT if you install elsewhere.

: "${ATV_DFTBPLUS_ROOT:=/nobackup/proj/disk/teoroo/personal/jicli594/work/molcrafts/software/dftbplus-24.1}"
: "${ATV_DFTB_PARAMS_ROOT:=/nobackup/proj/disk/teoroo/personal/jicli594/work/molcrafts/software/dftb-params}"

if [[ ! -x "${ATV_DFTBPLUS_ROOT}/bin/dftb+" ]]; then
  echo "error: dftb+ not found at ${ATV_DFTBPLUS_ROOT}/bin/dftb+" >&2
  echo "  expected prebuilt 24.1 under molcrafts/software/ (see software/README.md)" >&2
  return 1 2>/dev/null || exit 1
fi

export PATH="${ATV_DFTBPLUS_ROOT}/bin:${PATH}"
export ATV_DFTBPLUS_SK_MIO="${ATV_DFTB_PARAMS_ROOT}/mio-1-1"
export ATV_DFTBPLUS_SK_3OB="${ATV_DFTB_PARAMS_ROOT}/3ob-3-1"

if [[ ! -f "${ATV_DFTBPLUS_SK_MIO}/H-H.skf" ]]; then
  echo "warning: mio SK missing at ${ATV_DFTBPLUS_SK_MIO}/H-H.skf" >&2
fi

echo "DFTB+ ready: $(command -v dftb+)"
dftb+ --version 2>&1 | head -3
echo "  ATV_DFTBPLUS_SK_MIO=${ATV_DFTBPLUS_SK_MIO}"
echo "  ATV_DFTBPLUS_SK_3OB=${ATV_DFTBPLUS_SK_3OB}"
