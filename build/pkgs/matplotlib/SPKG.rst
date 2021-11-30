matplotlib: Python 2D plotting library
======================================

Description
-----------

From the Matplotlib website: matplotlib is a python 2D plotting library
which produces publication quality figures in a variety of hardcopy
formats and interactive environments across platforms. matplotlib can be
used in python scripts, the python and ipython shell (ala matlab or
mathematica), web application servers, and six graphical user interface
toolkits.

License
-------

The Matplotlib license - see
http://matplotlib.sourceforge.net/users/license.html: Matplotlib only
uses BSD compatible code, and its license is based on the PSF license.
See the Open Source Initiative licenses page for details on individual
licenses. Non-BSD compatible licenses (eg LGPL) are acceptable in
matplotlib Toolkits. For a discussion of the motivations behind the
licencing choice, see Licenses.


Upstream Contact
----------------

https://matplotlib.org

The matplotlib mailing lists: see
http://sourceforge.net/projects/matplotlib

Dependencies
------------

-  python
-  numpy
-  setuptools (>= 0.7)
-  freetype
-  patch (used in spkg-install)
-  dateutil
-  pyparsing
-  tornado
-  kiwisolver


Build Instructions/Changes
--------------------------

-  NOTE: To drastically cut down on spkg size, we delete the internal
   testing images. To do this, we repackage the tarball by removing
   the contents of ``lib/matplotlib/tests/baseline_images/*``, this is
   done by the ``spkg-src`` script.

-  ``setup.py.patch``: disable loading of Tests. Otherwise, ``setup.py``
   raises an error because it can't find the deleted files
   from ``src/lib/matplotlib/tests/baseline_images/*``

-  NOTE: as of matplotlib-1.0.0 and Sage 4.6, Sage does not use
   $HOME/.matplotlib by default. Instead, it sets MPLCONFIGDIR to
   a subdirectory in $DOT_SAGE, see src/bin/sage-env
