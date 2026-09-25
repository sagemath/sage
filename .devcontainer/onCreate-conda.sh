#!/bin/sh

set -eu

# Interactive shells should enter the development environment rather than
# Conda's base environment.  The environment itself is synchronized by the
# updateContentCommand, which also runs when a prebuild is refreshed.
conda_dir=${CONDA_DIR:-/opt/conda}
"$conda_dir/bin/conda" config --set auto_activate_base false
"$conda_dir/bin/conda" init bash
