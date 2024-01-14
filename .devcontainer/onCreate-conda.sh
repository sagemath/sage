# Do not keep running on errors
set -e
# Show commands run
set -x

# Create conda environment
mkdir /tmp/conda-pkgs && mkdir -p /opt/conda && ln -s /tmp/conda-pkgs /opt/conda/pkgs
conda install --quiet mamba -n base -c conda-forge -y
df -h / /tmp
mamba env create --file src/environment-dev-3.11-linux.yml || mamba env update --file src/environment-dev-3.11-linux.yml
conda init bash
df -h / /tmp
mamba clean --all --quiet --yes
df -h / /tmp

# Build sage
conda run -n sage-dev ./bootstrap
conda run -n sage-dev ./configure --with-python=/opt/conda/envs/sage-dev/bin/python --prefix=/opt/conda/envs/sage-dev
conda run -n sage-dev pip install --no-build-isolation -v -v -e ./pkgs/sage-conf ./pkgs/sage-setup 
conda run -n sage-dev pip install --no-build-isolation -v -v -e ./src
