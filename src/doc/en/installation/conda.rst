.. _sec-installation-conda:

Install from conda-forge
========================

SageMath can be installed from `conda-forge <https://conda-forge.org>`_ on Linux
and macOS running x86-64 that most current desktops and laptops use.

To install SageMath, install `Miniconda <https://conda.io/miniconda.html>`_ and
then type in the following commands in a terminal:

* Add the conda-forge channel: ``conda config --add channels conda-forge``
* Update all packages: ``conda update --all``
* Create a new environment containing SageMath: ``conda create -n sage sage``
* Enter the new environment: ``source activate sage``
* Start SageMath: ``sage``
