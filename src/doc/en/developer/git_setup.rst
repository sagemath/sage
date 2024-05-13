.. _chapter-git-setup:

==============
Setting Up Git
==============


.. _section-git-install:

Installing Git
--------------

Depending on your platform, use the following to install Git:

Linux
    See :ref:`spkg_git` for the installation command on your
    Linux distribution.

Windows (WSL)
    We strongly recommend to install the package using the Linux
    distribution's package manager.  Native Windows installations of
    Git may also work, but there are possible pitfalls.

macOS
    Install the Xcode Command Line Tools.

Some further resources for installation help are:

* `Section 1.5 of the Git book
  <https://git-scm.com/book/en/v2/Getting-Started-Installing-Git>`_

* The `Git homepage <http://git-scm.com>`_ for the most recent
  information

* `Git install help page on GitHub <https://github.com/git-guides/install-git>`_


.. _section-git-setup-name:

Configuring Git
---------------

Assuming your name ``alice`` and email address ``alice@wonderland.com``,

.. CODE-BLOCK:: shell-session

    [alice@localhost ~]$ git config --global user.name "Alice Adventure"
    [alice@localhost ~]$ git config --global user.email alice@wonderland.com

This will write the settings into your Git configuration file
``~/.gitconfig`` with your name and email:

.. CODE-BLOCK:: text

    [user]
        name = Alice Adventure
        email = alice@wonderland.com

Of course, replace ``Alice Adventure`` and ``alice@wonderland.com`` with your
actual name and email address.

This is the basic Git configuration for now. For further tips on configuring
Git, see :ref:`section-git-configuration`.

