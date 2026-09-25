#!/bin/sh

set -eu

environment_name=sage
conda_dir=${CONDA_DIR:-/opt/conda}
environment_prefix="$conda_dir/envs/$environment_name"

if [ -d "$environment_prefix/conda-meta" ]; then
    "$conda_dir/bin/mamba" install --yes --name "$environment_name" sage
else
    "$conda_dir/bin/mamba" create --yes --name "$environment_name" sage
fi

"$environment_prefix/bin/sage" --version
"$conda_dir/bin/conda" clean --all --yes

activation_command=". $conda_dir/etc/profile.d/conda.sh && conda activate $environment_name"
touch "$HOME/.bashrc"
if ! grep -Fqx "$activation_command" "$HOME/.bashrc"; then
    printf '\n%s\n' "$activation_command" >> "$HOME/.bashrc"
fi
