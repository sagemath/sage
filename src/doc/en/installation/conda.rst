.. _sec-installation-conda:

Install from conda-forge
========================

SageMath can be installed on Linux and macOS via Conda from the
`conda-forge <https://conda-forge.org>`_ conda channel.

Both the ``x86_64`` (Intel) architecture and the ``arm64``/``aarch64``
architectures (including Apple Silicon, M1, M2, M3, M4) are supported.

You will need a working Conda installation: either Miniforge, Miniconda or
Anaconda. If you don't have one yet, we recommend installing `Miniforge
<https://github.com/conda-forge/miniforge>`_ as follows. In a terminal,

.. code-block:: shell

    $ curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
    $ bash Miniforge3-$(uname)-$(uname -m).sh

* Miniforge uses conda-forge as the default channel. However, if you are using
  Miniconda or Anaconda, set it up to use conda-forge:

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

      $ mamba create -n sage sage

.. tab:: conda

  .. code-block:: shell

      $ conda create -n sage sage

To use Sage from there,

* Enter the new environment: ``conda activate sage``
* Start SageMath: ``sage``

If there are any installation failures, please report them to
the conda-forge maintainers by opening a `GitHub Issue for
conda-forge/sage-feedstock <https://github.com/conda-forge/sage-feedstock/issues>`_.
