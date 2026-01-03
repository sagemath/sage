#! /bin/sh

SYSTEM=$(build/bin/sage-guess-package-system)

if [ "$SYSTEM" = "fedora" ] && [ -f /etc/fedora-release ]; then
    # Need to use --setopt=tsflags="" to avoid errors with gphelp
    dnf5 install -y pari-gp --setopt=tsflags=""

    # Mitigate upstream packaging bug: https://bugzilla.redhat.com/show_bug.cgi?id=2332429
    # by swapping the incorrectly installed OpenCL-ICD-Loader for the expected ocl-icd
    dnf5 -y swap --repo='fedora' OpenCL-ICD-Loader ocl-icd
fi

eval $(build/bin/sage-print-system-package-command $SYSTEM "$@" update)
eval $(build/bin/sage-print-system-package-command $SYSTEM --yes --ignore-missing install $(build/bin/sage-get-system-packages $SYSTEM $(build/bin/sage-package list :standard:)))

# Disable build isolation following the advice of https://mesonbuild.com/meson-python/how-to-guides/editable-installs.html#build-dependencies
# Install build dependencies manually as workaround for https://github.com/astral-sh/uv/issues/1516
uv venv --clear
. $UV_PROJECT_ENVIRONMENT/bin/activate # https://github.com/astral-sh/uv/issues/14022
uv pip install \
    meson-python \
    "cypari2 >=2.2.1" \
    "cysignals >=1.11.2, != 1.12.0" \
    "cython >=3.0, != 3.0.3, < 3.1.0" \
    "gmpy2 >=2.1.5" \
    memory_allocator \
    "numpy >=1.25" \
    jinja2 \
    setuptools
uv sync --frozen --inexact --no-build-isolation -v --config-settings=builddir=build/build-$SYSTEM-$devcontainerId
