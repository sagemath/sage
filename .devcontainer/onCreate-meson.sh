#! /bin/sh
# Run this script from SAGE_ROOT. Invoke with "--sudo" if sudo is needed.
export PATH=$(pwd)/build/bin:$PATH
SYSTEM=$(sage-guess-package-system)
eval $(sage-print-system-package-command $SYSTEM "$@" update)
eval $(sage-print-system-package-command $SYSTEM --yes "$@" --spkg install python3)
eval $(sage-print-system-package-command $SYSTEM --yes --ignore-missing install $(sage-get-system-packages $SYSTEM $(sage-package list :standard:)))
