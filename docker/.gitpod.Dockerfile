# Use minimal Ubuntu installation that includes mamba
FROM condaforge/mambaforge

# Some basic system packages
RUN apt update && apt-get install -yq --no-install-recommends sudo gpg curl lsb-release

# Make Docker available, like the default gitpod image does
# from https://github.com/gitpod-io/workspace-images/blob/main/chunks/tool-docker/Dockerfile @ 3f0988f2d06768d22d0aa1454ef0e963b0db65f3
# - removed unnecessary "sudo"
# - replaced use of "install-packages" (https://github.com/gitpod-io/workspace-images/blob/main/base/install-packages)
# https://docs.docker.com/engine/install/ubuntu/
RUN curl -fsSL https://download.docker.com/linux/ubuntu/gpg | gpg --dearmor -o /usr/share/keyrings/docker-archive-keyring.gpg \
    && echo "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/docker-archive-keyring.gpg] https://download.docker.com/linux/ubuntu \
    $(lsb_release -cs) stable" | tee /etc/apt/sources.list.d/docker.list > /dev/null \
    && apt update \
    && apt-get install -yq --no-install-recommends docker-ce docker-ce-cli containerd.io

# Workaround so that vscode internals (such as git) find things installed in the conda env (notably ssh which is required to contribute to trac)
ENV PATH $PATH:/workspace/sagetrac-mirror/venv/bin
# Default to non-admin user
USER gitpod
