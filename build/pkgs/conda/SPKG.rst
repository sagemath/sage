conda: Installs conda-forge packages in an isolated environment
===============================================================

Description
-----------

conda-forge is a community-led collection of recipes, build infrastructure
and distributions for the conda package manager.

mamba is a fast cross-platform package manager, a compatible
reimplementation of the conda package manager in C++.

micromamba is a tiny version of the mamba package manager. It is a statically linked
C++ executable with a separate command line interface. It does not need a base
environment and does not come with a default version of Python.

Installing this script package installs ``micromamba`` in ``$SAGE_LOCAL/bin``
and creates an isolated conda environment within ``$SAGE_LOCAL/var/lib/sage/conda``.


License
-------

BSD-3-Clause


Upstream Contact
----------------

https://conda-forge.org/

https://github.com/mamba-org/mamba

https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html
