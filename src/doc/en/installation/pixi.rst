.. _sec-installation-pixi:

Developing with Pixi
====================

Pixi is the recommended way to set up a Sage development checkout. It installs
Python, compilers, libraries, and development tools from conda-forge into a
local environment. The repository's ``pixi.lock`` records the package versions
for each supported platform and Python version.

Prerequisites
-------------

Install Git and `Pixi <https://pixi.prefix.dev/latest/installation/>`_.
An existing Conda installation is not required.

The workspace includes environments for Linux (x86-64 and AArch64), macOS
(Intel and Apple Silicon), and Windows (x86-64).
On macOS, version 13.3 or later is required. Install the Xcode Command Line
Tools with ``xcode-select --install``.

Native Windows support is experimental and does not yet provide all Sage
features. For the usual Sage development workflow on Windows, use WSL and
follow the Linux instructions. For a native build, install the
`Visual Studio 2022 Build Tools <https://visualstudio.microsoft.com/downloads/>`_
with the C++ build tools and Windows SDK, then run the commands below from
the x64 Native Tools Command Prompt for Visual Studio 2022.

Build and run Sage
------------------

Clone the repository and install the locked development environment:

.. code-block:: console

    $ git clone --origin upstream https://github.com/sagemath/sage.git
    $ cd sage
    $ pixi install --locked
    $ pixi run build

Pixi selects the appropriate packages for your platform. The ``--locked``
option checks that ``pixi.lock`` matches ``pixi.toml`` and does not change the
lockfile.

The build task installs your checkout in editable mode using the environment's
build dependencies. It does not download a separate Sage binary package.
Start Sage or test a source file with:

.. code-block:: console

    $ pixi run sage
    $ pixi run test src/sage/rings/integer.pyx

Additional arguments are passed to Sage's doctester, for example:

.. code-block:: console

    $ pixi run test -p 4 --long src/sage/graphs

Alternatively, enter the environment with ``pixi shell``. You can then run
``sage``, ``python``, and other installed commands directly. Leave the shell
with ``exit``. Use ``pixi run`` or ``pixi shell`` to activate the compiler
environment as well as selecting the correct Python interpreter.

Python versions and rebuilds
----------------------------

The default environment uses Python 3.13. The named environments ``py313`` and
``py314`` select Python 3.13 and 3.14 respectively. For example:

.. code-block:: console

    $ pixi run -e py314 build
    $ pixi run -e py314 sage
    $ pixi run -e py314 test src/sage/rings/integer.pyx

Each environment keeps its own compiled files in ``sage-build`` under the
environment directory, such as ``.pixi/envs/default/sage-build``. This prevents
build files from being shared between incompatible Python or library versions.
Use this directory instead of ``builddir`` in the
:ref:`Meson commands <section-meson-build-directory>`.

An editable installation uses Python source changes immediately and rebuilds
changed Cython modules when importing them. Rerun ``pixi run build`` after
changing package metadata, adding dependencies, or updating the environment.
Use the same ``-e`` option for the build and subsequent commands.

The build task starts a fresh build when the installed Conda packages change.
To request a clean rebuild manually, use ``pixi run build --clean``. This
clears only the selected environment's Sage build directory.

Updating dependencies
---------------------

After pulling changes, install the updated lockfile and rebuild:

.. code-block:: console

    $ pixi install --locked
    $ pixi run build

When changing Sage's dependencies, edit ``pyproject.toml`` and regenerate the
dependency section of ``pixi.toml``:

.. code-block:: console

    $ pixi run python tools/update-pixi.py
    $ pixi lock

The generator shares the package-name mappings and platform exclusions used by
the Conda environment generator. Keep the generated dependency sections and
``pixi.lock`` in the same commit as the dependency changes. To verify that
the generated sections are current without modifying them, use
``pixi run python tools/update-pixi.py --check``.

Use ``pixi update PACKAGE`` to update a selected dependency within the existing
constraints. Review the lockfile changes and rebuild and test Sage before
committing them. See the `Pixi documentation
<https://pixi.prefix.dev/latest/>`_ for environment and lockfile commands.
