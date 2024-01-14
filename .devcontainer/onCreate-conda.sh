# Do not keep running on errors
set -e

# Create conda environment
mkdir /tmp/conda-pkgs && mkdir -p /opt/conda && ln -s /tmp/conda-pkgs /opt/conda/pkgs
conda install mamba -n base -c conda-forge -y
df -h
mamba env create --quiet --file src/environment-dev-3.11-linux.yml || mamba env update --quiet --file src/environment-dev-3.11-linux.yml
conda init bash
df -h

# Build sage
conda run -n sage-dev ./bootstrap
conda run -n sage-dev ./configure --with-python=/opt/conda/envs/sage-dev/bin/python --prefix=/opt/conda/envs/sage-dev
conda run -n sage-dev pip install --no-build-isolation -v -v -e ./pkgs/sage-conf ./pkgs/sage-setup 
conda run -n sage-dev pip install --no-build-isolation -v -v -e ./src
