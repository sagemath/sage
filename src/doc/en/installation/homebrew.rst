.. _sec-installation-homebrew:

Installing Sage with Homebrew
=============================

On Linux and macOS, the community-maintained `dimpase/tap
<https://github.com/dimpase/homebrew-tap#sagemath>`_ provides two ways to use
Sage. Install ``dimpase/tap/sagemath`` to get the Sage command and Jupyter
kernel. To build a Sage checkout yourself, install ``dimpase/tap/sagemath-deps``
for its native libraries, Python packages, datasets, and build tools.

The tap publishes prebuilt packages (bottles) for x86-64 and ARM64 Linux and
for macOS Tahoe on Apple Silicon. Run ``brew info dimpase/tap/sagemath``
after setting up the taps to see the packaged Sage version.

Set up Homebrew
---------------

Install `Homebrew <https://brew.sh/>`_ and follow the installer's instructions
to add ``brew`` to your shell. On Linux, Homebrew (formerly Linuxbrew) uses
``/home/linuxbrew/.linuxbrew`` as its default prefix for bottles. See
`Homebrew on Linux <https://docs.brew.sh/Homebrew-on-Linux>`_ for prerequisites.

Both formulae use packages from three community taps. To add the taps and
allow Homebrew to load their formulae, run::

    $ brew tap dimpase/tap
    $ brew tap dimpase/gap
    $ brew tap macaulay2/tap
    $ brew trust --tap dimpase/tap dimpase/gap macaulay2/tap

Install and run Sage
--------------------

Homebrew installs the dependencies automatically when you install Sage::

    $ brew install dimpase/tap/sagemath
    $ sage
    sage: factor(2026)
    2 * 1013

The formula also installs the SageMath Jupyter kernel. If you use JupyterLab
from another environment, add Homebrew's kernels to its search path before
launching it::

    $ export JUPYTER_PATH="$(brew --prefix)/share/jupyter${JUPYTER_PATH:+:$JUPYTER_PATH}"
    $ jupyter lab

Select the SageMath kernel when creating a notebook.

.. _sec-homebrew-source:

Build a Sage checkout with Homebrew dependencies
------------------------------------------------

For development, start in the root of your Sage source checkout. The dependency
formula supplies ``uv`` and the dependencies for the Sage version packaged by
the tap. Other revisions of Sage may require different dependency versions.

Create a virtual environment using Homebrew's Python 3.14. The
``--system-site-packages`` option makes the bottled Python packages available
inside it::

    $ brew install dimpase/tap/sagemath-deps
    $ uv venv --python "$(brew --prefix python@3.14)/bin/python3.14" --system-site-packages .venv
    $ source .venv/bin/activate
    $ source "$(brew --prefix dimpase/tap/sagemath-deps)/share/sagemath-deps/python3.14-env.sh"

The environment script sets the search paths for the native libraries and
build tools. Install the checkout in editable mode using those dependencies,
then check that Sage imports successfully::

    $ uv pip install --no-build-isolation --no-deps --editable . \
        -Csetup-args=-Dbuild-docs=false -Csetup-args=-Dc_std=gnu17 \
        -Csetup-args=--wrap-mode=nofallback
    $ python -c 'from sage.all import factor; print(factor(2026))'
    2 * 1013

When returning to the checkout in a new shell, repeat the two ``source``
commands before running ``sage`` or rebuilding.
