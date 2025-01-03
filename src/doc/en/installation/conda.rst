.. _sec-installation-conda:

Install from conda-forge
========================

SageMath can be installed on Linux and macOS via Conda from the
`conda-forge <https://conda-forge.org>`_ conda channel.

Both the ``x86_64`` (Intel) architecture and the ``arm64``/``aarch64``
architectures (including Apple Silicon, M1, M2, M3, M4) are supported.

You will need a working Conda installation: either Miniforge (or Mambaforge),
Miniconda or Anaconda. If you don't have one yet, we recommend installing
`Miniforge <https://github.com/conda-forge/miniforge>`_ as
follows. In a terminal,

.. code-block:: shell

   $ curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
   $ bash Miniforge3-$(uname)-$(uname -m).sh

* Miniforge (and Mambaforge) use conda-forge as the default channel.

* If you are using Miniconda or Anaconda, set it up to use conda-forge:

  * Add the conda-forge channel: ``conda config --add channels conda-forge``

  * Change channel priority to strict: ``conda config --set channel_priority strict``

If you installed Miniforge (or Mambaforge), we recommend to use
`mamba <https://mamba.readthedocs.io/en/latest/index.html>`_ in the following,
which uses a faster dependency solver than ``conda``.

.. _sec-installation-conda-binary:

Installing all of SageMath from conda (not for development)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Create a new conda environment containing SageMath, either with ``mamba`` or ``conda``:

.. tab:: mamba

  .. code-block:: shell

      $ mamba create -n sage sage python=X

.. tab:: conda

  .. code-block:: shell

      $ conda create -n sage sage python=X

where ``X`` is version of Python, e.g. ``3.12``.

To use Sage from there,

* Enter the new environment: ``conda activate sage``
* Start SageMath: ``sage``

If there are any installation failures, please report them to
the conda-forge maintainers by opening a `GitHub Issue for
conda-forge/sage-feedstock <https://github.com/conda-forge/sage-feedstock/issues>`_.

.. _sec-installation-conda-develop:

Using conda to provide all dependencies for the Sage library
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can build and install the Sage library from source, using conda to
provide all of its dependencies. This bypasses most of the build
system of the Sage distribution and is the fastest way to set up an
environment for Sage development.

Here we assume that you are using a git checkout.

- Optionally, set the build parallelism for the Sage library. Use
  whatever the meaningful value for your machine is - no more than
  the number of cores::

    $ export SAGE_NUM_THREADS=24

- Create and activate a new conda environment with the dependencies of Sage
  and a few additional developer tools:

  .. tab:: mamba

    .. code-block:: shell

        $ mamba env create --file environment-3.12-linux.yml --name sage-dev
        $ conda activate sage-dev

  .. tab:: conda

    .. code-block:: shell

        $ conda env create --file environment-3.12-linux.yml --name sage-dev
        $ conda activate sage-dev

  Alternatively, you can use ``environment-3.12-linux.yml`` or
  ``environment-optional-3.12-linux.yml``, which will only install standard
  (and optional) packages without any additional developer tools.

  A different Python version can be selected by replacing ``3.12`` with the 
  desired version.

- Bootstrap the source tree and install the build prerequisites and the Sage library::

    $ ./bootstrap
    $ pip install --no-build-isolation --config-settings editable_mode=compat -v -v --editable ./src

  If you encounter any errors, try to install the ``sage-conf`` package first::

    $ pip install --no-build-isolation -v -v --editable ./pkgs/sage-conf_conda

  and then run the last command again.

- Verify that Sage has been installed::

    $ sage -c 'print(version())'
    SageMath version 10.2.beta4, Release Date: 2023-09-24

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

  $ pip install --no-build-isolation --config-settings editable_mode=compat -v -v --editable src

In order to update the conda environment later, you can run::

  $ mamba env update --file environment-3.12-linux.yml --name sage-dev

To build the documentation, use::

  $ pip install --no-build-isolation -v -v --editable ./pkgs/sage-docbuild
  $ sage --docbuild all html

.. NOTE::

   The switch ``--config-settings editable_mode=compat`` restores the
   `legacy setuptools implementation of editable installations
   <https://setuptools.pypa.io/en/latest/userguide/development_mode.html>`_.
   Adventurous developers may omit this switch to try the modern,
   PEP-660 implementation of editable installations, see :issue:`34209`.

.. NOTE::

  You can update the conda lock files by running
  ``.github/workflows/conda-lock-update.py`` or by running
  ``conda-lock --platform linux-64 --filename environment-3.12-linux.yml --lockfile environment-3.12-linux.lock``
  manually.
