.. _sec-GNU-Linux:

Linux Package Managers
======================

SageMath is available from various distributions and can be installed
by package managers.

As of Sage 10.2, we can recommend the following distributions, which
provide well-maintained and up-to-date SageMath packages:
`Arch Linux <https://archlinux.org/>`_
and `Void Linux <https://voidlinux.org/>`_.
Gentoo users might want to give the
`sage-on-gentoo <https://github.com/cschwan/sage-on-gentoo>`_ overlay
a try.

.. only:: html

   .. |codespace-downstream-archlinux-latest| image:: https://github.com/codespaces/badge.svg
      :target: https://codespaces.new/sagemath/sage?devcontainer_path=.devcontainer%2Fdownstream-archlinux-latest%2Fdevcontainer.json

   .. |package-downstream-archlinux-latest| image:: https://repology.org/badge/version-for-repo/arch/sagemath.svg
      :target: https://repology.org/project/sagemath/versions

   .. |codespace-downstream-voidlinux-latest| image:: https://github.com/codespaces/badge.svg
      :target: https://codespaces.new/sagemath/sage?devcontainer_path=.devcontainer%2Fdownstream-voidlinux-latest%2Fdevcontainer.json

   .. |package-downstream-voidlinux-latest| image:: https://repology.org/badge/version-for-repo/void_x86_64/sagemath.svg
      :target: https://repology.org/project/sagemath/versions

   .. |codespace-downstream-gentoo-overlay| image:: https://github.com/codespaces/badge.svg
      :target: https://codespaces.new/sagemath/sage?devcontainer_path=.devcontainer%2Fdownstream-gentoo-overlay%2Fdevcontainer.json

   .. list-table::
      :widths: 35 30 35
      :header-rows: 0

      * - `downstream-archlinux-latest
          <https://github.com/sagemath/sage/tree/develop/.devcontainer/downstream-archlinux-latest>`_
        - |package-downstream-archlinux-latest|
        - |codespace-downstream-archlinux-latest|

      * - `downstream-voidlinux-latest
          <https://github.com/sagemath/sage/tree/develop/.devcontainer/downstream-voidlinux-latest>`_
        - |package-downstream-voidlinux-latest|
        - |codespace-downstream-voidlinux-latest|

      * - `downstream-gentoo-overlay
          <https://github.com/sagemath/sage/tree/develop/.devcontainer/downstream-gentoo-overlay>`_
        - `cschwan/sage-on-gentoo <https://github.com/cschwan/sage-on-gentoo>`_
        - |codespace-downstream-gentoo-overlay|

See `repology.org: sagemath
<https://repology.org/project/sagemath/versions>`_ for information
about versions of SageMath packages in various other distributions.

**Do not install a version of Sage older than 9.5.**
If you are on an older version of your distribution and a recent
version of SageMath is only available on a newer version of the
distribution, consider upgrading your distribution.

See `the _sagemath dummy package <../reference/spkg/_sagemath.html>`_
for the names of packages that provide a standard installation of
SageMath, including documentation and Jupyter.

The  `GitHub wiki page Distribution <https://github.com/sagemath/sage/wiki/Distribution>`_ collects information
regarding packaging and distribution of SageMath.
