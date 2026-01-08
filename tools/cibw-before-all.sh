#!/usr/bin/env bash

# This script is run from the root of the Sage source tree before cibuildwheel starts building any wheels.

# Exit on error
set -e

export PATH=$(pwd)/build/bin:$PATH
SYSTEM=$(sage-guess-package-system)

if [ "$SYSTEM" = "homebrew" ]; then
    source .homebrew-build-env
fi

eval $(build/bin/sage-print-system-package-command $SYSTEM "$@" update)
eval $(build/bin/sage-print-system-package-command $SYSTEM --yes --no-install-recommends --spkg install _bootstrap _prereq)
./bootstrap
eval $(build/bin/sage-print-system-package-command $SYSTEM --yes --no-install-recommends --ignore-missing install $(build/bin/sage-get-system-packages $SYSTEM $(build/bin/sage-package list :standard:)))
./configure --enable-build-as-root --without-system-python3
export MAKE="make -j8"
make sagelib-build-deps
