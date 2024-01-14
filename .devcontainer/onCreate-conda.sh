# Do not keep running on errors
set -e
# Show commands run
set -x

# Create conda environment
# A space bottleneck occurs at build time of the container
# when the smallest (2-core, 32 GB) configuration is used.
mkdir /tmp/conda-pkgs && mkdir -p /opt/conda && ln -s /tmp/conda-pkgs /opt/conda/pkgs
conda install --quiet mamba -n base -c conda-forge -y
# The following line shows the available space.
# At the time of writing:
# Filesystem      Size  Used Avail Use% Mounted on
# overlay          32G   23G  7.5G  76% /
# /dev/sda1        44G   52K   42G   1% /tmp
# The 7.5G (also shared with /workspaces) are not enough for our
# conda packages, which is why we redirect them into /tmp
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
