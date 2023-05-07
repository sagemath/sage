.. _chapter-git-setup:

==============
Setting Up Git
==============

To work on the Sage source code, you need a working Git installation,
configured at least to use your name and email address for commits.

For further information about Git, see :ref:`chapter-git-background`, but this
is not required for a beginning Sage developer.


.. _section-git-install:

Installing Git
--------------

First, try Git on the command line by typing ``git``. Most platforms will have it installed by
default if other development tools are installed. If that fails, use the
following to install Git:

Debian / Ubuntu
    Run ``sudo apt-get install git-core``

Fedora
    Run ``sudo yum install git-core``

Windows (Cygwin)
    Install the Cygwin package Git. Do not attempt to use native
    Windows installations of Git.

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

The commit message of any change contains your name and email address to
acknowledge your contribution and to have a point of contact if there are
questions in the future. Filling it in is required if you want to share your
changes. The simplest way to do this is from the command line. Assuming your
name ``alice`` and email address ``alice@wonderland.com``,

.. CODE-BLOCK:: shell-session

    [alice@localhost ~]$ git config --global user.name "Alice Adventure"
    [alice@localhost ~]$ git config --global user.email alice@wonderland.com

This will write the settings into your :ref:`Git configuration file
<section-git-configuration>` with your name and email:

.. CODE-BLOCK:: text

    [user]
        name = Alice Adventure
        email = alice@wonderland.com

Of course, replace ``Alice Adventure`` and ``alice@wonderland.com`` with your
actual name and email address.

