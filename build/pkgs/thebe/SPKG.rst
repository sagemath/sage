thebe: Add live Jupyter interaction to static websites
======================================================

Description
-----------

Jupyter javascript plugin for static sites. Thebe takes the Jupyter
front end, and make it work outside of the notebook context.

This is used by Sage's Sphinx-based documentation build system to
produce html documentation that can be turned live (see
https://trac.sagemath.org/ticket/20690).

License
-------

MIT


Upstream Contact
----------------

- Home page: https://oreillymedia.github.io/thebe/
- Source: https://github.com/oreillymedia/thebe/

Dependencies
------------

None.


Special Update/Build Instructions
---------------------------------

There are no release numbers, hence find the latest commit, download
https://github.com/oreillymedia/thebe/archive/$%7BCOMMIT%7D.zip and
rename it thebe-${COMMIT:0:8}.zip
