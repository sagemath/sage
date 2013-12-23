.. _set-up-fork:

==================
Set Up Your Fork
==================

First you follow the instructions for :ref:`forking`.

Overview
========

::

   git clone git@github.com:your-user-name/sagenb.git
   cd sagenb
   git remote add upstream git://github.com/sagemath/sagenb.git

In etail
=========

Clone Your Fork
---------------

#. Clone your fork to the local computer with ``git clone
   git@github.com:your-user-name/sagenb.git``
#. Investigate.  Change directory to your new repo: ``cd sagenb``. Then
   ``git branch -a`` to show you all branches.  You'll get something
   like::

      * master
      remotes/origin/master

   This tells you that you are currently on the ``master`` branch, and
   that you also have a ``remote`` connection to ``origin/master``.
   What remote repository is ``remote/origin``? Try ``git remote -v`` to
   see the URLs for the remote.  They will point to your github fork.

   Now you want to connect to the upstream `Sage Notebook github`_ repository, so
   you can merge in changes from trunk.

.. _linking-to-upstream:

Linking Your Repository to the Upstream Repo
--------------------------------------------

::

   cd sagenb
   git remote add upstream git://github.com/sagemath/sagenb.git

``upstream`` here is just the arbitrary name we're using to refer to the
main `Sage Notebook`_ repository at `Sage Notebook github`_.

Note that we've used ``git://`` for the URL rather than ``git@``.  The
``git://`` URL is read only.  This means we that we can't accidentally
(or deliberately) write to the upstream repo, and we are only going to
use it to merge into our own code.

Just for your own satisfaction, show yourself that you now have a new
'remote', with ``git remote -v show``, giving you something like::

   upstream     git://github.com/sagemath/sagenb.git (fetch)
   upstream     git://github.com/sagemath/sagenb.git (push)
   origin       git@github.com:your-user-name/sagenb.git (fetch)
   origin       git@github.com:your-user-name/sagenb.git (push)

.. include:: links.inc

