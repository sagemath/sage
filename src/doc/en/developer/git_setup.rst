.. _chapter-git-setup:

==============
Setting Up Git
==============

To work on the Sage source code, you need

* a working ``git`` installation, see :ref:`section-git-install`.

* configure ``git`` to use your name and email address for commits, see
  :ref:`section-git-setup-name`.

The :ref:`chapter-git-background` chapter contains further information
about ``git`` that might be useful to some but are not required.


.. _section-git-install:

Installing Git
--------------

First, try ``git`` on the command line. Most distributions will have
it installed by default if other development tools are installed. If
that fails, use the following to install git:

Debian / Ubuntu
    ``sudo apt-get install git-core``

Fedora
    ``sudo yum install git-core``

Windows (Cygwin)
    Install the Cygwin package ``git``. Do not attempt to use native
    Windows installations of ``git``.

Windows (WSL)
    We strongly recommend to install the package using the Linux
    distribution's package manager.  Native Windows installations of
    ``git`` may also work, but there are possible pitfalls.

macOS
    Install the Xcode Command Line Tools.

Some further resources for installation help are:

* `Section 1.5 of the git book
  <https://git-scm.com/book/en/v2/Getting-Started-Installing-Git>`_

* The `git homepage <http://git-scm.com>`_ for the most recent
  information.

* `Github install help pages <http://help.github.com>`_


.. _section-git-setup-name:

Your Name and Email
-------------------

The commit message of any change contains your name and email address
to acknowledge your contribution and to have a point of contact if
there are questions in the future; filling it in is required if you
want to share your changes. The simplest way to do this is from the
command line:

.. CODE-BLOCK:: shell-session

    [user@localhost ~] git config --global user.name "Your Name"
    [user@localhost ~] git config --global user.email you@yourdomain.example.com

This will write the settings into your :ref:`git configuration file
<section-git-configuration>` with your name and email:

.. CODE-BLOCK:: text

    [user]
        name = Your Name
        email = you@yourdomain.example.com

Of course you'll need to replace ``Your Name`` and ``you@yourdomain.example.com``
with your actual name and email address.
