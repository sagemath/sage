rst2ipynb: Convert reStructuredText files to Jupyter notebooks
==============================================================

Description
-----------

The rst2pynb program converts a standalone reStructuredText file to a
Jupyter notebook file.

This is currently achieved by converting to markdown with pandoc and
then to Jupyter notebook using notedown, plus some configuration and
tweaks.

License
-------

BSD 3-Clause License


Upstream Contact
----------------

Authors: Scott Sievert and Nicolas M. Thiéry Home page:
https://github.com/nthiery/rst-to-ipynb

Special Update/Build Instructions
---------------------------------

Fetch tarball from https://pypi.python.org/pypi/rst2ipynb/

As it is written in Haskell, pandoc must be installed from the distro.

The main rationale for having a notedown package in Sage (rather than
just let pip fetch it) is that the version on pipy (1.5.0, 2015-10-07)
is outdated and lacks important features / fixes for us.
