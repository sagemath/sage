# Do not keep running on errors
set -e

# Create conda environment
conda install mamba -n base -c conda-forge -y
mamba env create --file src/environment-dev-3.11-linux.yml || mamba env update --file src/environment-dev-3.11-linux.yml
conda init bash

# Build sage
conda run -n sage-dev ./bootstrap
conda run -n sage-dev pip install --no-build-isolation -v -v -e ./src
