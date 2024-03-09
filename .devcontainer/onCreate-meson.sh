#! /bin/sh
SYSTEM=$(build/bin/sage-guess-package-system)
eval $(build/bin/sage-print-system-package-command $SYSTEM "$@" update)
eval $(build/bin/sage-print-system-package-command $SYSTEM --yes --ignore-missing install $(build/bin/sage-get-system-packages $SYSTEM $(build/bin/sage-package list :standard:)))

# Create a virtual environment (some systems do not allow to install packages globally)
# Virtual envionment has to be created outside of the source directory for meson to work
python3 -m venv ../venv --system-site-packages
. ../venv/bin/activate
pip3 install --upgrade pip

# Needed for lrcalc
eval $(build/bin/sage-print-system-package-command $SYSTEM --yes --ignore-missing install wget)
wget math.rutgers.edu/~asbuch/lrcalc/lrcalc-2.1.tar.gz \
    && tar zxvf lrcalc-2.1.tar.gz \
    && cd lrcalc-2.1 \
    && ./configure \
    && make \
    && make install \
    && cd .. \
    && rm -fr lrcalc-2.1 lrcalc-2.1.tar.gz

 # Needed for fpylll
FPLLL_VERSION=5.4.5
wget https://github.com/fplll/fplll/releases/download/${FPLLL_VERSION}/fplll-${FPLLL_VERSION}.tar.gz \
    && tar -xf fplll-${FPLLL_VERSION}.tar.gz \
    && cd fplll-${FPLLL_VERSION} \
    && ./configure \
    && make \
    && make install \
    && cd .. \
    && rm -fr fplll-${FPLLL_VERSION} fplll-${FPLLL_VERSION}.tar.gz

pip3 install --no-build-isolation --config-settings=builddir=builddir --editable . -v
