.. _chapter-git-setup:

==============
Setting Up Git
==============

To work on the Sage source code, you need

* a working git installation, see :ref:`section-git-install`. Sage
  actually comes with git, see below. However, it is recommended that
  you have a system-wide install if only to save you some typing.

* configure git to use your name and email address for commits, see
  :ref:`section-git-setup-name`. The Sage development scripts will
  prompt you if you don't. But, especially if you use git for other
  projects in the future as well, you really should configure git.

The :ref:`chapter-git-background` chapter contains further information
about git that might be useful to some but are not required.


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

Windows
    Download and install `msysGit
    <http://code.google.com/p/msysgit/downloads/list>`_

OS X
    Use the `git OSX installer
    <http://code.google.com/p/git-osx-installer/downloads/list>`_

Finally, Sage includes git. Obviously there is a chicken-and-egg
problem to checkout the Sage source code from its git repository, but
one can always download a Sage source tarball or binary
distribution. You can then run git via the ``sage -git`` command line
switch. So, for example, ``git help`` becomes ``sage -git help`` and
so on. Note that the examples in the developer guide will assume that
you have a system-wide git installation.

Some further resources for installation help are:

* `Chapter 2 of the git book
  <http://book.git-scm.com/2_installing_git.html>`_

* The `git homepage <http://git-scm.com>`_ for the most recent
  information.

* `Github install help pages <http://help.github.com>`_


.. _section-git-setup-name:

Your Name and Email
-------------------

The commit message of any change contains your name and email address
to acknowledge your contribution and to have a point of contact if
there are questions in the future; Filling it in is required if you
want to share your changes. The simplest way to do this is from the
command line::

    [user@localhost ~] git config --global user.name "Your Name"
    [user@localhost ~] git config --global user.email you@yourdomain.example.com

This will write the settings into your :ref:`git configuration file
<section-git-configuration>` with your name and email::

    [user]
        name = Your Name
        email = you@yourdomain.example.com

Of course you'll need to replace ``Your Name`` and ``you@yourdomain.example.com``
with your actual name and email address.
