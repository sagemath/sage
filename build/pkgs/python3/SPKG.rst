python3: The Python programming language
========================================

Description
-----------

By default, Sage will try to use system's ``python3`` to set up a virtual
environment, a.k.a. `venv <https://docs.python.org/3.10/library/venv.html>`_
rather than building a Python 3 installation from scratch.

Sage will accept versions 3.8.x to 3.10.x.

You can also use ``--with-python=/path/to/python3_binary`` to tell Sage to use
``/path/to/python3_binary`` to set up the venv. Note that setting up the venv requires
a number of Python modules to be available within the Python in question. Currently,
as of Sage 9.7, these modules are as follows: ``sqlite3``, ``ctypes``, ``math``,
``hashlib``, ``crypt``, ``socket``, ``zlib``, ``distutils.core``, ``ssl`` -
they will be checked for by the ``configure`` script.

Use the ``configure`` option ``--without-system-python3`` in case you want Python 3
built from scratch.


Upstream Contact
----------------

https://www.python.org
