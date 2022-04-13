.. _sec-installation-conda:

Install from conda-forge
========================

SageMath can be installed via Conda from the
`conda-forge <https://conda-forge.org>`_ conda channel.

This works on Linux and macOS on ``x86_64`` processors,
and on Linux on ``aarch64`` processors (using Miniforge).

This requires a working Conda installation: either Miniforge, Miniconda
or Anaconda. If you don't have one yet, we recommend installing
`Miniforge <https://github.com/conda-forge/miniforge#miniforge3>`_.

Miniforge uses conda-forge as the default channel. If you are
using Miniconda or Anaconda, set it up to use conda-forge:

* Add the conda-forge channel: ``conda config --add channels conda-forge``
* Change channel priority to strict: ``conda config --set channel_priority strict``

Optionally, use `mamba <https://github.com/mamba-org/mamba>`_
which uses a faster dependency solver than `conda`.

.. code-block:: shell

    conda install mamba


Installing all of SageMath from conda (not for development)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Create a new conda environment containing SageMath, either with ``mamba`` or ``conda``:

.. code-block:: shell

    mamba create -n sage sage python=X
    conda create -n sage sage python=X

where ``X`` is version of Python, e.g. ``3.9``.

To use Sage from there,

* Enter the new environment: ``conda activate sage``
* Start SageMath: ``sage``


Using conda to provide system packages for the Sage distribution
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If Conda is installed (check by typing ``conda info``), there are two ways to
prepare for installing SageMath from source:

  - If you are using a git checkout::

      $ ./bootstrap

  - Create a new empty environment and activate::

      $ conda create -n sage-build
      $ conda activate sage-build

  - Install standard packages recognized by sage's ``spkg-configure`` mechanism::

      $ conda env update --file environment.yml -n sage-build

  - Or install all standard and optional packages recognized by sage::

      $ conda env update --file environment-optional.yml -n sage-build

  - Then the SageMath distribution will be built using the compilers provided by Conda
    and using many packages installed by Conda::

      $ ./bootstrap
      $ ./configure --prefix=$CONDA_PREFIX
      $ make

Using conda to provide all dependencies for the Sage library (experimental)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can build and install the Sage library from source, using conda to
provide all of its dependencies. This bypasses most of the build
system of the Sage distribution and is the fastest way to set up an
environment for Sage development.

Note that this is still an experimental feature and may not work as
intended.

Here we assume that you are using a git checkout.

  - Optionally, set the build parallelism for the Sage library. Use
    whatever the meaningful value for your machine is - no more than
    the number of cores::

      $ export SAGE_NUM_THREADS=24

  - As a recommended step, install the ``mamba`` package manager.  If
    you skip this step, replace ``mamba`` by ``conda`` in the
    following steps::

      $ conda install mamba

  - Create and activate a new conda environment that provides the
    bootstrapping prerequisites. You can replace 3.9 by another Python
    version::

      $ mamba create -n sage-build python=3.9 \
            gettext autoconf automake libtool pkg-config
      $ conda activate sage-build

  - Run ``bootstrap``; this generates the files ``src/environment*.yml`` used
    in the next step::

      $ ./bootstrap

  - Populate the conda environment with the dependencies of Sage::

      $ mamba env update -n sage-build -f src/environment.yml  # alternatively, use

    Alternatively, you can use ``src/environment-optional.yml``, which will
    install some additional packages.

  - Activate the conda environment again::

      $ conda activate sage-build

  - Run the ``configure`` script::

      $ ./configure --with-python=$CONDA_PREFIX/bin/python             \
                    --prefix=$CONDA_PREFIX                             \
                    $(for pkg in $(./sage -package list :standard:     \
                                     --has-file spkg-configure.m4      \
                                     --has-file distros/conda.txt); do \
                          echo --with-system-$pkg=force;               \
                      done)

  - Install the build prerequisites of the Sage library::

      $ pip install --no-build-isolation -v -v --editable pkgs/sage-conf pkgs/sage-setup

  - Install the Sage library::

      $ pip install --no-build-isolation -v -v --editable src

  - Verify that Sage has been installed::

      $ sage -c 'print(version())'
      SageMath version 9.6.beta5, Release Date: 2022-03-12

Note that ``make`` is not used at all.  All dependencies
(including all Python packages) are provided by conda.

Thus, you will get a working version of Sage much faster.  However,
note that this will invalidate the use of any Sage-the-distribution
commands such as ``sage -i``. Do not use them.

By using ``pip install --editable`` in the above steps, the Sage
library is installed in editable mode.  This means that when you only
edit Python files, there is no need to rebuild the library; it
suffices to restart Sage.

After editing any Cython files, rebuild by repeating the command::

  $ pip install --no-build-isolation -v -v --editable src
