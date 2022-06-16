# Use minimal Ubuntu installation that includes mamba
FROM condaforge/mambaforge
# Workaround so that vscode internals (such as git) find things installed in the conda env (notably ssh which is required to contribute to trac)
ENV PATH $PATH:/workspace/sagetrac-mirror/venv/bin
# Default to non-admin user
USER gitpod
