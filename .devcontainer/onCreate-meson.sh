#! /bin/sh

SYSTEM=$(build/bin/sage-guess-package-system)
eval $(build/bin/sage-print-system-package-command $SYSTEM "$@" update)
eval $(build/bin/sage-print-system-package-command $SYSTEM --yes --ignore-missing install $(build/bin/sage-get-system-packages $SYSTEM $(build/bin/sage-package list :standard:)))

# Disable build isolation following the advice of https://mesonbuild.com/meson-python/how-to-guides/editable-installs.html#build-dependencies
uv sync --frozen
