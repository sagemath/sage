.. _sec-installation-conda:

Install from conda-forge
========================

SageMath can be installed on Linux and macOS via Conda from the
`conda-forge <https://conda-forge.org>`_ conda channel.

Both the ``x86_64`` (Intel) architecture and the ``arm64``/``aarch64``
architectures (including Apple Silicon, M1) are supported.

You will need a working Conda installation: either Mambaforge/Miniforge,
Miniconda or Anaconda. If you don't have one yet, we recommend installing
`Mambaforge <https://github.com/conda-forge/miniforge#mambaforge>`_ as
follows. In a terminal,

.. code-block:: shell

   $ curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh
   $ sh Mambaforge-$(uname)-$(uname -m).sh

* Mambaforge and Miniforge use conda-forge as the default channel.

* If you are using Miniconda or Anaconda, set it up to use conda-forge:

  * Add the conda-forge channel: ``conda config --add channels conda-forge``

  * Change channel priority to strict: ``conda config --set channel_priority strict``

Optionally, use `mamba <https://github.com/mamba-org/mamba>`_,
which uses a faster dependency solver than ``conda``.
If you installed Mambaforge, it is already provided. Otherwise, use

.. code-block:: shell

   $ conda install mamba


.. _sec-installation-conda-binary:

Installing all of SageMath from conda (not for development)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Create a new conda environment containing SageMath, either with ``mamba`` or ``conda``:

.. code-block:: shell

    $ mamba create -n sage sage python=X        # either
    $ conda create -n sage sage python=X        # or

where ``X`` is version of Python, e.g. ``3.9``.

To use Sage from there,

* Enter the new environment: ``conda activate sage``
* Start SageMath: ``sage``

If there are any installation failures, please report them to
the conda-forge maintainers by opening a `GitHub Issue for
conda-forge/sage-feedstock <https://github.com/conda-forge/sage-feedstock/issues>`_.


.. _sec-installation-conda-source:

Using conda to provide system packages for the Sage distribution
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If Conda is installed (check by typing ``conda info``), one can install SageMath
from source as follows:

  - If you are using a git checkout::

      $ ./bootstrap-conda

  - Create a new conda environment including all standard packages
    recognized by sage, and activate it::

      $ conda env create --file environment.yml --name sage-build
      $ conda activate sage-build

    Alternatively, use ``environment-optional.yml`` in place of
    ``environment.yml`` to create an environment with all standard and optional
    packages recognized by sage.

  - Then the SageMath distribution will be built using the compilers provided by Conda
    and using many packages installed by Conda::

      $ ./bootstrap
      $ ./configure --with-python=$CONDA_PREFIX/bin/python \
                    --prefix=$CONDA_PREFIX
      $ make


.. _sec-installation-conda-develop:

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

  - As a recommended step, install the ``mamba`` package manager. If
    you skip this step, replace ``mamba`` by ``conda`` in the
    following steps::

      $ conda install mamba

  - Generate the conda environment files ``src/environment*.yml`` used
    in the next step::

      $ ./bootstrap-conda

  - Create and activate a new conda environment with the dependencies of Sage
    and a few additional developer tools::

      $ mamba env create --file src/environment-dev.yml --name sage-dev
      $ conda activate sage-dev

    Alternatively, you can use ``src/environment.yml`` or
    ``src/environment-optional.yml``, which will only install standard
    (and optional) packages without any additional developer tools.

    By default, the most recent version of Python supported by Sage is
    installed. You can use the additional option ``python=3.9`` in the above
    ``env create`` command to select another Python version (here 3.9). 

  - Run the ``configure`` script::

      $ ./bootstrap
      $ ./configure --with-python=$CONDA_PREFIX/bin/python             \
                    --prefix=$CONDA_PREFIX                             \
                    $(for pkg in $(./sage -package list :standard:     \
                                     --exclude rpy2                    \
                                     --has-file spkg-configure.m4      \
                                     --has-file distros/conda.txt); do \
                          echo --with-system-$pkg=force;               \
                      done)

  - Install the build prerequisites and the Sage library::

      $ pip install --no-build-isolation -v -v --editable ./pkgs/sage-conf ./pkgs/sage-setup
      $ pip install --no-build-isolation -v -v --editable ./src

  - Verify that Sage has been installed::

      $ sage -c 'print(version())'
      SageMath version 9.6.beta5, Release Date: 2022-03-12

Note that ``make`` is not used at all. All dependencies
(including all Python packages) are provided by conda.

Thus, you will get a working version of Sage much faster.  However,
note that this will invalidate the use of any Sage-the-distribution
commands such as ``sage -i``. Do not use them.

By using ``pip install --editable`` in the above steps, the Sage
library is installed in editable mode.  This means that when you only
edit Python files, there is no need to rebuild the library; it
suffices to restart Sage.

After editing any Cython files, rebuild the Sage library using::

  $ pip install --no-build-isolation -v -v --editable src

In order to update the conda environment later, you can run::

  $ mamba env update --file src/environment-dev.yml --name sage-dev
