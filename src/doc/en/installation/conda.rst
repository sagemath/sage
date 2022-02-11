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

Create a new conda environment containing SageMath, either with ``mamba`` or ``conda``:

.. code-block:: shell

    mamba create -n sage sage python=X
    conda create -n sage sage python=X

where ``X`` is version of Python, e.g. ``3.8``.

To use Sage from there,

* Enter the new environment: ``conda activate sage``
* Start SageMath: ``sage``

Instructions for using Conda for SageMath development are on
`the Conda page of the Sage wiki <https://wiki.sagemath.org/Conda>`__.
