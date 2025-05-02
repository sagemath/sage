#! /bin/sh

SYSTEM=$(build/bin/sage-guess-package-system)

if [ "$SYSTEM" = "fedora" ]; then
    # Need to use --setopt=tsflags="" to avoid errors with gphelp
    yum install -y pari-gp --setopt=tsflags=""
fi

eval $(build/bin/sage-print-system-package-command $SYSTEM "$@" update)
eval $(build/bin/sage-print-system-package-command $SYSTEM --yes --ignore-missing install $(build/bin/sage-get-system-packages $SYSTEM $(build/bin/sage-package list :standard:)))

# Disable build isolation following the advice of https://mesonbuild.com/meson-python/how-to-guides/editable-installs.html#build-dependencies
# Install build dependencies manually as workaround for https://github.com/astral-sh/uv/issues/1516
uv pip install \
    meson-python \
    "cypari2 >=2.2.1" \
    "cysignals >=1.11.2, != 1.12.0" \
    "cython >=3.0, != 3.0.3" \
    "gmpy2 ~=2.1.b999" \
    memory_allocator \
    "numpy >=1.25" \
    jinja2
uv sync --frozen --no-build-isolation-package=sagemath
