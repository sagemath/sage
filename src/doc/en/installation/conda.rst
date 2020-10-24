.. _sec-installation-conda:

Install from conda-forge
========================

SageMath can be installed from `conda-forge <https://conda-forge.org>`_ on Linux
and macOS running x86-64 that most current desktops and laptops use.

To install SageMath, install `Miniconda <https://conda.io/miniconda.html>`_ and
then type in the following commands in a terminal:

* Add the conda-forge channel: ``conda config --add channels conda-forge``
* Create a new environment containing SageMath: ``conda create -n sage sage python=X``, where
  ``X`` is version of Python, e.g. ``3.7``
* Enter the new environment: ``conda activate sage``
* Start SageMath: ``sage``

.. note:: The dependency resolution process of conda is slow; a faster experimental resolver
   we tried with success is `mamba <https://github.com/QuantStack/mamba>`_. It can be used as follows.

   .. code-block:: shell

      conda install mamba -c conda-forge # installs mamba
      mamba create -n sage sage -c conda-forge # replaces "conda create..."


