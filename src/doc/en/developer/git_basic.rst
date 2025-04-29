.. highlight:: shell-session

.. _chapter-git-basic:

==========
Git Basics
==========

Git is a tool to exchange commits (file changes) and branches (organized of
commits) with other developers.

As a distributed revision control system, Git does not have the notion of a
central server. However, for Sage development, Git communicates with other
developers via `the Sage repository <https://github.com/sagemath/sage>`_ on
GitHub. Hence we assume that throughout this guide.

.. _section-git-ssh:

Git authentication through SSH
==============================

In order to push changes securely to a remote repository, Git uses public-key
cryptography. This section will show you how to set up the necessary
cryptographic keys for the case that you want to use SSH(Secure Shell) protocol
to authenticate your Git to GitHub, instead of HTTPS protocol.


.. _section-github-ssh-key:

Generating your SSH keys
------------------------

Check whether you already have suitable SSH keys by inspecting ``.ssh``
directory in your home directory. If you don't have suitable SSH keys yet, you
can create a key pair with the ``ssh-keygen`` tool.

Follow either `the detailed instructions
<https://git-scm.com/book/en/v2/Git-on-the-Server-Generating-Your-SSH-Public-Key>`_
or the following brief instructions::

    [alice@localhost ~]$ ssh-keygen
    Generating public/private rsa key pair.
    Enter file in which to save the key (/home/alice/.ssh/id_rsa):
    Enter passphrase (empty for no passphrase):
    Enter same passphrase again:
    Your identification has been saved in /home/alice/.ssh/id_rsa.
    Your public key has been saved in /home/alice/.ssh/id_rsa.pub.
    The key fingerprint is:
    ce:32:b3:de:38:56:80:c9:11:f0:b3:88:f2:1c:89:0a alice@localhost
    The key's randomart image is:
    +--[ RSA 2048]----+
    |  ....           |
    |   ..            |
    |   .o+           |
    | o o+o.          |
    |E + .  .S        |
    |+o .   o.        |
    |. o   +.o        |
    |      oB         |
    |     o+..        |
    +-----------------+

This will generate a new random private RSA key
in the ``.ssh`` folder in your home directory. By default, they are

``~/.ssh/id_rsa``
  Your private key. Keep safe. **Never** hand it out to anybody.

``~/.ssh/id_rsa.pub``
  The corresponding public key. This and only this file can be safely
  disclosed to third parties.

The ``ssh-keygen`` tool will let you generate a key with a different
file name, or protect it with a passphrase. Depending on how much you
trust your own computer or system administrator, you can leave the
passphrase empty to be able to login without any human intervention.


Adding your public key for authentication to GitHub
---------------------------------------------------

Follow the procedure `Adding a new SSH key to your GitHub account
<https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account>`_.
Then check that it works by::

    [alice@localhost sage]$ git remote add origin git@github.com:alice/sage.git
    [alice@localhost sage]$ git remote -v
    origin  git@github.com:alice/sage.git (fetch)
    origin  git@github.com:alice/sage.git (push)


.. _section-git-push:

Pushing your changes to a remote
================================

Push your branch to the remote ``origin`` with either ::

    [alice@localhost sage]$ git push --set-upstream origin HEAD:my_branch

or ::

    [alice@localhost sage]$ git push origin HEAD:my_branch

if your branch already has an upstream branch. Here "upstream" means the the
remote ``origin``, which is *upstream* to your local Git repo.

Here, ``HEAD`` means that you are pushing the most recent commit (and, by
extension, all of its parent commits) of the current local branch to the remote
branch.


.. _section-git-checkout:

Checking out a PR
=================

If you want to work with the changes of a PR branch, you must
make a local copy of the branch. In particular, Git has no concept of directly
working with the remote branch, the remotes are only bookmarks for
things that you can get from/to the remote server. Hence, the first
thing you should do is to get everything from the branch
into your local repository. This is achieved by::

    [alice@localhost sage]$ git fetch upstream pull/12345/head
    remote: Enumerating objects: 12, done.
    remote: Counting objects: 100% (12/12), done.
    remote: Compressing objects: 100% (3/3), done.
    remote: Total 12 (delta 9), reused 11 (delta 9), pack-reused 0
    Unpacking objects: 100% (12/12), 2.22 KiB | 206.00 KiB/s, done.
    From https://github.com/sagemath/sage
     * branch                  refs/pull/12345/head -> FETCH_HEAD

The ``pull/12345/head`` branch refers to the branch of the PR #12345 of the
remote ``upstream``. The branch is now temporarily (until you fetch something
else) stored in your local Git database under the alias ``FETCH_HEAD``. In the
second step, we make it available as a new local branch and switch to it. Your
local branch can have a different name, for example::

    [alice@localhost sage]$ git checkout -b my_branch FETCH_HEAD
    Switched to a new branch 'my_branch'

creates a new branch in your local Git repository named ``my_branch``
and modifies your local Sage filesystem tree to the state of the files
in the branch. You can now edit files and commit changes to your
local branch.


.. _section-git-pull:

Getting changes from a remote
=============================

A common task during development is to synchronize your local copy of the
branch with the branch on the GitHub Sage repo. In particular, assume you
downloaded the branch of a PR made by someone else, say Bob, and made some
suggestions for improvements on the PR. Now Bob incorporated your suggestions
into his branch, and you want to get the added changes to complete your review.
Assuming that you originally got your local branch as in
:ref:`section-git-checkout`, you can just issue::

    [bob@localhost sage]$ git pull upstream pull/12345/head
    From https://github.com/sagemath/sage
     * branch                  refs/pull/35608/head -> FETCH_HEAD
    Merge made by the 'ort' strategy.
     src/doc/common/python3.inv          | Bin 98082 -> 131309 bytes
     src/doc/common/update-python-inv.sh |   7 ++++---
     2 files changed, 4 insertions(+), 3 deletions(-)

This command downloads the changes from the branch of the PR and merges
them into your local branch.


.. _section-git-pull-develop:

Updating develop
================

The ``develop`` branch can be updated just like any other branch. However, your
local copy of the develop branch should stay **identical** to the GitHub Sage repo develop
branch.

If you accidentally added commits to your local copy of ``develop``, you must
delete them before updating the branch.

One way to ensure that you are notified of potential problems is to use ``git
pull --ff-only``, which will raise an error if a non-trivial merge would be
required::

    [alice@localhost sage]$ git checkout develop
    [alice@localhost sage]$ git pull --ff-only upstream develop

If this pull fails, then something is wrong with the local copy of the
master branch. To switch to the correct Sage master branch, use::

    [alice@localhost sage]$ git checkout develop
    [alice@localhost sage]$ git reset --hard upstream/develop


.. _section-git-merge:

Merging and rebasing
====================

Sometimes, a new version of Sage is released while you work on a Git branch.

Let us assume you started ``my_branch`` at commit ``B``. After a while, your
branch has advanced to commit ``Z``, but you updated ``develop`` (see
:ref:`section-git-pull-develop`) and now your Git history looks like this:

.. CODE-BLOCK:: text

                     X---Y---Z my_branch
                    /
               A---B---C---D develop

How should you deal with such changes? In principle, there are two ways:

* **Rebase:** The first solution is to **replay** commits ``X,Y,Z`` atop of the
  new ``develop``. This is called **rebase**, and it rewrites your current
  branch:

  .. CODE-BLOCK:: text

      git checkout my_branch
      git rebase -i develop

  In terms of the commit graph, this results in:

  .. CODE-BLOCK:: text

                             X'--Y'--Z' my_branch
                            /
               A---B---C---D develop

  Note that this operation rewrites the history of ``my_branch`` (see
  :ref:`section-git-rewriting-history`). This can lead to problems if somebody
  began to write code atop of your commits ``X,Y,Z``. It is safe otherwise.

  **Alternatively**, you can rebase ``my_branch`` while updating ``develop`` at the
  same time (see :ref:`section-git-pull`):

  .. CODE-BLOCK:: text

    git checkout my_branch
    git pull -r develop

* **Merging** your branch with ``develop`` will create a new commit above the two
  of them:

  .. CODE-BLOCK:: text

      git checkout my_branch
      git merge develop

  The result is the following commit graph:

  .. CODE-BLOCK:: text

                     X---Y---Z---W my_branch
                    /           /
               A---B---C-------D develop

  - **Pros:** you did not rewrite history (see
    :ref:`section-git-rewriting-history`).The additional commit is then easily
    pushed to the git repository and distributed to your collaborators.

  - **Cons:** it introduced an extra merge commit that would
    not be there had you used rebase.

  **Alternatively**, you can merge ``my_branch`` while updating ``develop`` at the
  same time (see :ref:`section-git-pull`):

  .. CODE-BLOCK:: text

    git checkout my_branch
    git pull develop

**In case of doubt** use merge rather than rebase. There is less risk involved,
and rebase in this case is only useful for branches with a very long history.


.. _section-git-mergetool:

Merge tools
===========

Simple conflicts can be easily solved with Git only (see :ref:`section-git-conflict`)

For more complicated ones, a range of specialized programs are
available. Because the conflict marker includes the hash of the most recent
common parent, you can use a three-way diff::

    [alice@laptop]$ git mergetool

    This message is displayed because 'merge.tool' is not configured.
    See 'git mergetool --tool-help' or 'git help config' for more details.
    'git mergetool' will now attempt to use one of the following tools:
    meld opendiff kdiff3 [...] merge araxis bc3 codecompare emerge vimdiff
    Merging:
    fibonacci.py

    Normal merge conflict for 'fibonacci.py':
      {local}: modified file
      {remote}: modified file
    Hit return to start merge resolution tool (meld):

If you don't have a favourite merge tool we suggest you try `meld
<http://meldmerge.org/>`_ (cross-platform). The result looks like the following
screenshot.

.. IMAGE:: static/meld-screenshot.png

The middle file is the most recent common parent; on the right is
Bob's version and on the left is Alice's conflicting version. Clicking
on the arrow moves the marked change to the file in the adjacent
pane.


.. _section-git-conflict:

Conflict resolution
-------------------

Merge conflicts happen if there are overlapping edits, and they are an
unavoidable consequence of distributed development. Fortunately,
resolving them is common and easy with Git. As a hypothetical example,
consider the following code snippet:

.. CODE-BLOCK:: python

    def fibonacci(i):
        """
        Return the `i`-th Fibonacci number
        """
        return fibonacci(i-1) * fibonacci(i-2)

This is clearly wrong. Two developers, namely Alice and Bob, decide to
fix it. Bob corrected the seed values:

.. CODE-BLOCK:: python

    def fibonacci(i):
       """
       Return the `i`-th Fibonacci number
       """
       if i > 1:
           return fibonacci(i-1) * fibonacci(i-2)
       return [0, 1][i]

and turned those changes into a new commit::

    [bob@laptop sage]$ git add fibonacci.py
    [bob@laptop sage]$ git commit -m 'return correct seed values'

He made his changes a PR to the GitHub Sage repo and got it merged to the
``develop`` branch. His ``fibonacci`` function is not yet perfect but is
certainly better than the original.

Meanwhile, Alice changed the multiplication to an addition since that is the
correct recursion formula:

.. CODE-BLOCK:: python

    def fibonacci(i):
        """
        Return the `i`-th Fibonacci number
        """
        return fibonacci(i-1) + fibonacci(i-2)

and merged her branch with the latest ``develop`` branch fetched from the GitHub Sage repo::

    [alice@home sage]$ git add fibonacci.py
    [alice@home sage]$ git commit -m 'corrected recursion formula, must be + instead of *'
    [alice@home sage]$ git fetch upstream develop:develop
    [alice@home sage]$ git merge develop
    ...
    CONFLICT (content): Merge conflict in fibonacci.py
    Automatic merge failed; fix conflicts and then commit the result.

The file now looks like this:

.. skip    # doctester confuses >>> with input marker

.. CODE-BLOCK:: python

    def fibonacci(i):
        """
        Return the `i`-th Fibonacci number
        """
    <<<<<<< HEAD
        return fibonacci(i-1) + fibonacci(i-2)
    =======
        if i > 1:
            return fibonacci(i-1) * fibonacci(i-2)
        return [0, 1][i]
    >>>>>>> 41675dfaedbfb89dcff0a47e520be4aa2b6c5d1b

The conflict is shown between the conflict markers ``<<<<<<<`` and
``>>>>>>>``. The first half (up to the ``=======`` marker) is Alice's
current version, the second half is Bob's version. The 40-digit hex
number after the second conflict marker is the SHA1 hash of the most
recent common parent of both.

It is now Alice's job to resolve the conflict by reconciling their
changes, for example by editing the file. Her result is:

.. CODE-BLOCK:: python

    def fibonacci(i):
        """
        Return the `i`-th Fibonacci number
        """
        if i > 1:
            return fibonacci(i-1) + fibonacci(i-2)
        return [0, 1][i]

And then upload both her original change *and* her merge commit to the GitHub Sage repo::

    [alice@laptop sage]$ git add fibonacci.py
    [alice@laptop sage]$ git commit -m "merged Bob's changes with mine"

The resulting commit graph now has a loop::

    [alice@laptop sage]$ git log --graph --oneline
    *   6316447 merged Bob's changes with mine
    |\
    | * 41675df corrected recursion formula, must be + instead of *
    * | 14ae1d3 return correct seed values
    |/
    * 14afe53 initial commit
    [alice@laptop sage]$ git push origin

This time, there is no merge conflict since Alice's branch already merged the ``develop`` branch.

