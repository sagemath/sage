# Do not keep running on errors
set -e

# Create conda environment
./bootstrap-conda
conda install mamba -n base -c conda-forge -y
mamba env create --file src/environment-dev.yml || mamba env update --file src/environment-dev.yml
conda init bash

# Build sage
conda run -n sage-dev ./bootstrap
conda run -n sage-dev ./configure --with-python=/opt/conda/envs/sage-dev/bin/python --prefix=/opt/conda/envs/sage-dev
conda run -n sage-dev pip install --no-build-isolation -v -v -e ./pkgs/sage-conf ./pkgs/sage-setup 
conda run -n sage-dev pip install --no-build-isolation -v -v -e ./src
