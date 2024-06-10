#! /bin/sh
# The portability-.../devcontainer.json configurations run this script after
# the container is started.
#
# The script assumes that it is run from SAGE_ROOT.
#
# If "config.log" or "logs" are symlinks (for example, created by 'tox -e local-...',
# or after https://github.com/sagemath/sage/issues/33262), they might point outside of
# the dev container, so remove them. Likewise for upstream.
for f in config.log logs upstream; do
    if [ -L $f ]; then
        rm -f $f
    fi
done
# If possible (ensured after https://github.com/sagemath/sage/issues/33262), keep the
# logs in the container.
if [ ! -d logs ]; then
    ln -s /sage/logs logs
fi
# Bootstrap, configure, and build the Sage distribution, reusing the Sage
# installation from the prebuilt image.
set -e
set -x
make configure
if [ -x /sage/config.status ]; then
    eval ./configure $(/sage/config.status --config) --enable-build-as-root --prefix=/sage/local --with-sage-venv
else
    ./configure --enable-build-as-root --prefix=/sage/local --with-sage-venv
fi
make build V=0
