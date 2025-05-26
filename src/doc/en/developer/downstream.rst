=============================
Packaging SageMath Downstream
=============================

This document is intended for downstream maintainers (e.g., Linux distribution
package maintainers) who wish to create redistributable builds of Sage.

.. contents::
   :local:
   :depth: 2

Dependencies
============

SageMath relies on a broad set of Python and system libraries. These must be
provided by the downstream distribution. The definitive list of dependencies is
found in `pyproject.toml <https://github.com/sagemath/sage/blob/develop/pyproject.toml>`_.

These include:
 - `build-system.requires`: Python packages needed for building SageMath,
 - `project.dependencies`: Python packages required at runtime,
 - `project.optional-dependencies`: optional dependencies for additional
    functionality,
 - `external.build-requires` and `external.host-requires`: system dependencies
    needed for building,
 - `external.dependencies`: system libraries required at runtime.

The `external` section follows `PEP 725 <https://peps.python.org/pep-0725/>`_
and specifies dependencies in the form of ̀PURLs.
At the moment, there is no standard interface to translate these PURLs into
system package names. However, the names should be quite self-explanatory.
You may also consult the section :ref:`spkg` for a list of Sage's
dependencies and their corresponding system package names in various
distributions.

Build Procedure
===============

1. **Obtain the Source**:
   Clone the SageMath repository:

   .. code-block:: bash

      git clone https://github.com/sagemath/sage.git

    Alternatively, download the sdist tarball from the
    `SageMath PyPI project<https://pypi.org/project/sagemath/>`_ or from the
    `GitHub releases <https://github.com/sagemath/sage/releases>`_.

1. **Prepare the Build Environment**:
   Ensure a clean and consistent build environment with access to all
   required system libraries and Python packages.

2. **Build**:

    Create a wheel using the `build` module:

   .. code-block:: bash

      python -m build --wheel --no-isolation

    If you are sure that all dependencies are available, you may also add the
    `--skip-dependency-check` option.
    Moreover, if you care about reproducible builds, it is recommended to
    use `-Cbuild-dir=build` to specify a build directory, see this
    `Meson-Python issue <https://github.com/mesonbuild/meson-python/issues/671>`_.

3. **Install**:

    The resulting wheel can be installed using

   .. code-block:: bash

      python -m installer --destdir="<pkgdir>" dist/sagemath-*.whl

   where `<pkgdir>` is the directory where you want to install the package
   (usually a temporary directory for packaging).

4. **Test the Build**:

   Run the Sage tests to ensure functionality:

    .. code-block:: bash

        python -m sage.doctest --all

    However, some tests are known to fail, see :issue:`39872`.


If you maintain a downstream package and encounter build issues or patches
that may benefit others, please consider contributing back by reporting issues
or opening pull requests on the SageMath GitHub repository.

Other considerations:
- **Package naming**: Use `sagemath`, or `python-sagemath` if your distribution
has a convention for Python packages.

Example Downstream Packages
===========================

- `Arch Linux <https://archlinux.org/packages/extra/x86_64/sagemath>`_

