#! /usr/bin/env bash

set -e

if [ -z "$CONDA_PREFIX" ]; then
    echo >&2 "CONDA_PREFIX needs to be set"
    exit 1
fi

# https://docs.anaconda.com/anaconda/install/silent-mode/
# https://docs.anaconda.com/anaconda/user-guide/tasks/shared-pkg-cache/
# https://www.anaconda.com/conda-configuration-engine-power-users/

mkdir -p "$CONDA_PREFIX"
echo >  "$CONDA_PREFIX"/.condarc "# Generated by .ci/setup-conda.sh"

if [ -n "$SHARED_CACHE_DIR" ]; then
    mkdir -p "$SHARED_CACHE_DIR"/conda_pkgs
    echo >> "$CONDA_PREFIX"/.condarc "pkgs_dirs:"
    echo >> "$CONDA_PREFIX"/.condarc "  - $SHARED_CACHE_DIR/conda_pkgs"
fi
if [ -n "$USE_CONDARC" ]; then
    cat $USE_CONDARC >> "$CONDA_PREFIX"/.condarc
fi
if [ ! -x "$CONDA_PREFIX"/bin/conda ]; then
    mkdir -p "$CONDA_PREFIX"/conda-meta
    curl -L "$CONDA_INSTALLER_URL_BASE$CONDA_INSTALLER_FILE" -C - -o "$SHARED_CACHE_DIR"/"$CONDA_INSTALLER_FILE"
    bash "$SHARED_CACHE_DIR"/"$CONDA_INSTALLER_FILE" -b -f -p "$CONDA_PREFIX"
fi