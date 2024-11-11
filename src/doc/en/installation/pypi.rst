.. HIGHLIGHT:: shell-session

.. _sec-install-from-pypi:

Install from PyPI
=================

For installing Sage in a Python environment from PyPI, Sage provides the
``pip``-installable package
`sagemath-standard <https://pypi.org/project/sagemath-standard/>`__.

Unless you need to install Sage into a specific existing environment, we
recommend to create and activate a fresh virtual environment, for
example ``~/sage-venv/``:

::

           $ python3 -m venv ~/sage-venv
           $ source ~/sage-venv/bin/activate

As the first installation step, install
`sage_conf <https://pypi.org/project/sage-conf/>`__, which builds
various prerequisite packages in a subdirectory of ``~/.sage/``:

::

           (sage-venv) $ python3 -m pip install -v sage_conf

After a successful installation, a wheelhouse provides various Python
packages. You can list the wheels using the command:

::

           (sage-venv) $ ls $(sage-config SAGE_SPKG_WHEELS)

If this gives an error saying that ``sage-config`` is not found, check
any messages that the ``pip install`` command may have printed. You may
need to adjust your ``PATH``, for example by:

::

           $ export PATH="$(python3 -c 'import sysconfig; print(sysconfig.get_path("scripts", "posix_user"))'):$PATH"

Now install the packages from the wheelhouse and the
`sage_setup <https://pypi.org/project/sage-conf/>`__ package, and
finally install the Sage library:

::

           (sage-venv) $ python3 -m pip install $(sage-config SAGE_SPKG_WHEELS)/*.whl sage_setup
           (sage-venv) $ python3 -m pip install --no-build-isolation -v sagemath-standard

The above instructions install the latest stable release of Sage. To
install the latest development version instead, add the switch ``--pre``
to all invocations of ``python3 -m pip install``.

**NOTE:** PyPI has various other ``pip``-installable packages with the
word “sage” in their names. Some of them are maintained by the Sage
project, some are provided by Sage users for various purposes, and
others are entirely unrelated to Sage. Do not use the packages
``sage`` and ``sagemath``. For a curated list of packages, see the
chapter `Packages and
Features <https://doc.sagemath.org/html/en/reference/spkg/index.html>`__
of the Sage Reference Manual.

