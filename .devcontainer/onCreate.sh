#! /bin/sh
# Run this script from SAGE_ROOT. Invoke with "--sudo" if sudo is needed.
export PATH=$(pwd)/build/bin:$PATH
SYSTEM=$(sage-guess-package-system)
eval $(sage-print-system-package-command $SYSTEM "$@" update)
eval $(sage-print-system-package-command $SYSTEM --yes "$@" --spkg install _bootstrap _prereq python3 git $EXTRA_SAGE_PACKAGES)
if [ -n "$EXTRA_SYSTEM_PACKAGES" ]; then
    eval $(sage-print-system-package-command $SYSTEM --yes "$@" install $EXTRA_SYSTEM_PACKAGES)
fi
