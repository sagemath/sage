# Do not keep running on errors
set -e

# Create conda environment
conda config --env --add channels conda-forge
conda config --env --set channel_priority strict
conda update -y --all --override-channels -c conda-forge
conda install mamba=1 -n base -y
mamba env create -y --file environment-3.11-linux.yml || mamba env update -y --file environment-3.11-linux.yml
conda init bash

# Build sage
conda run -n sage-dev pip install --no-build-isolation -v -v -e .
