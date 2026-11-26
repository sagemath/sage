#!/bin/sh

set -eu

environment_name=sage-dev
environment_file=environment-3.12-linux.yml
conda_dir=${CONDA_DIR:-/opt/conda}
conda_command="$conda_dir/bin/conda"
environment_prefix="$conda_dir/envs/$environment_name"
environment_stamp="$environment_prefix/conda-meta/sage-environment.sha256"

expected_python=$(
    sed -n 's/^  - python=\([0-9][0-9]*\.[0-9][0-9]*\).*/\1/p' "$environment_file" |
        head -n 1
)
if [ -z "$expected_python" ]; then
    echo >&2 "cannot determine the Python version from $environment_file"
    exit 1
fi

environment_hash=$(sha256sum "$environment_file" | awk '{print $1}')
environment_changed=yes
current_python=
replace_environment=no

if [ -d "$environment_prefix/conda-meta" ]; then
    if [ -x "$environment_prefix/bin/python" ]; then
        current_python=$(
            "$environment_prefix/bin/python" -c \
                'import sys; print(f"{sys.version_info.major}.{sys.version_info.minor}")' \
                2>/dev/null || true
        )
    fi

    if [ "$current_python" != "$expected_python" ]; then
        echo "Replacing $environment_name: Python $current_python does not match $expected_python"
        replace_environment=yes
    elif [ -r "$environment_stamp" ] &&
         [ "$(sed -n '1p' "$environment_stamp")" = "$environment_hash" ]; then
        environment_changed=no
    fi
fi

if [ "$environment_changed" = yes ]; then
    # Leave the old stamp absent until dependency synchronization, the
    # editable build, and the Sage smoke test have all succeeded.
    rm -f "$environment_stamp"

    # Meson build directories cannot be reused across a Conda lock or Python
    # ABI change.  Remove them before resolving packages so they cannot make
    # the environment repair itself exhaust a Codespaces disk.
    for stale_build_dir in build/conda-cp*; do
        if [ -d "$stale_build_dir" ]; then
            echo "Removing stale build directory $stale_build_dir"
            rm -rf -- "$stale_build_dir"
        fi
    done

    # Free caches before asking Conda to write transaction metadata; a stale
    # prebuild may already have filled the Codespaces disk.
    "$conda_command" clean --all --yes
    if [ "$replace_environment" = yes ]; then
        "$conda_command" env remove --yes --prefix "$environment_prefix"
        "$conda_command" clean --all --yes
    fi

    if [ -d "$environment_prefix/conda-meta" ]; then
        "$conda_command" env update --prune \
            --prefix "$environment_prefix" --file "$environment_file"
    else
        "$conda_command" env create --yes \
            --prefix "$environment_prefix" --file "$environment_file"
    fi

    actual_python=$(
        "$environment_prefix/bin/python" -c \
            'import sys; print(f"{sys.version_info.major}.{sys.version_info.minor}")'
    )
    if [ "$actual_python" != "$expected_python" ]; then
        echo >&2 "expected Python $expected_python in $environment_name, found $actual_python"
        exit 1
    fi

    "$conda_command" clean --all --yes
fi

python_tag=$(
    "$environment_prefix/bin/python" -c \
        'import sys; print(f"cp{sys.version_info.major}{sys.version_info.minor}")'
)
build_dir="build/conda-$python_tag"

# Reinstall the editable metadata and entry point for the current checkout.
"$conda_command" run --no-capture-output --prefix "$environment_prefix" \
    python -m pip install --no-build-isolation --no-deps --verbose --editable . \
    --config-settings="builddir=$build_dir"
"$conda_command" run --no-capture-output --prefix "$environment_prefix" \
    python -m pip check

# Force the editable loader to compile now and verify the installed console
# entry point before Codespaces reports that the environment is ready.
"$conda_command" run --no-capture-output --prefix "$environment_prefix" \
    env MESONPY_EDITABLE_VERBOSE=1 \
    sage -c 'assert ZZ(2) + ZZ(2) == 4'

printf '%s\n' "$environment_hash" > "$environment_stamp"

activation_command="conda activate $environment_name"
touch "$HOME/.bashrc"
if ! grep -Fqx "$activation_command" "$HOME/.bashrc"; then
    printf '\n%s\n' "$activation_command" >> "$HOME/.bashrc"
fi
