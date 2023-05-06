.. highlight:: shell-session

.. _chapter-git-basic:

==========
Git Basics
==========

.. _section-git-ssh:

Git authentication through SSH
==============================

In order to push changes securely to a remote repository, git uses
public-key cryptography. This section will show you how to set up the
necessary cryptographic keys for Secure Shell (SSH).


Generating your SSH Keys
------------------------

Check whether you have already have suitable SSH keys by inspecting ``.ssh``
directory in your home directory. If you don't have suitable SSH keys yet, you
can create a key pair with the ``ssh-keygen`` tool.

Follow either `the detailed instructions
<https://git-scm.com/book/en/v2/Git-on-the-Server-Generating-Your-SSH-Public-Key>`_
or the following brief instructions::

    [alice@localhost ~]$ ssh-keygen
    Generating public/private rsa key pair.
    Enter file in which to save the key (/home/user/.ssh/id_rsa):
    Enter passphrase (empty for no passphrase):
    Enter same passphrase again:
    Your identification has been saved in /home/user/.ssh/id_rsa.
    Your public key has been saved in /home/user/.ssh/id_rsa.pub.
    The key fingerprint is:
    ce:32:b3:de:38:56:80:c9:11:f0:b3:88:f2:1c:89:0a user@localhost
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


.. _section-github-ssh-key:

Linking your Public Key to your GitHub Account
----------------------------------------------

In order to push your code directly to a branch on your remote ``origin``, Your
GitHub account needs to know your public key.

No action needed if you have already contributed to any other project on GitHub
and set up git credentials or SSH keys for this.

New users of GitHub should follow either `Caching your GitHub credentials in
Git
<https://docs.github.com/en/get-started/getting-started-with-git/caching-your-github-credentials-in-git>`_
or generate an SSH keypair, or use an already existing one, and upload the
public key to your GitHub account settings `Connecting to GitHub with SSH
<https://docs.github.com/en/authentication/connecting-to-github-with-ssh>`_


Adding your Public Key for authentication on another server
-----------------------------------------------------------

If you have an account on a lab or department computer that allows you
to log in remotely via SSH, you can now also use your SSH keys to
log in. Just copy the **public** key file (ending in ``.pub``) to
``~/.ssh/authorized_keys`` on the remote computer and make sure that
the file is only read/writable by yourself. Voila, the next time you
ssh into that machine you don't have to provide your password.

.. NOTE::

    In the command above we set up the remote to only track the
    ``master`` branch on the Sage repo (the ``-t master``
    option). This avoids clutter by not automatically downloading all
    branches ever created. But it also means that you will not fetch
    everything that is on the repo by default, and you need to explicitly
    tell git which branch you want to get from the repo. See the
    :ref:`section-git-checkout` section for examples.

Note that write operations (``push``) use the ssh protocol (specified by the ``git@``
part). For this to work, you need to have a GitHub account and to set up your ssh public
key.  Authentication is necessary if you want to upload anything to ensure
that it really is from you.

The above instructions set up the remote to perform read-only operations (``fetch``)
using HTTPS from a mirror of the Sage repo instead. The mirror is faster and
more reliable than our git server. However, this configuration is not recommended if
you use VS Code as an IDE.

If you want to use SSH only for both ``fetch`` and ``push``, use the
following commands instead::

    [alice@localhost sage]$ git remote add trac git@trac.sagemath.org:sage.git -t master
    [alice@localhost sage]$ git remote -v
    origin      https://github.com/sagemath/sage.git (fetch)
    origin      https://github.com/sagemath/sage.git (push)
    trac        git@trac.sagemath.org:sage.git (fetch)
    trac        git@trac.sagemath.org:sage.git (push)


.. _section-git-checkout:

Checking out a branch
---------------------

If you want to work with the changes in a PR, you must
make a local copy. In particular, Git has no concept of directly
working with the remote branch, the remotes are only bookmarks for
things that you can get from/to the remote server. Hence, the first
thing you should do is to get everything from the Sage repo's branch
into your local repository. This is achieved by::

    [user@localhost sage]$ git fetch upstream u/user/description
    remote: Counting objects: 62, done.
    remote: Compressing objects: 100% (48/48), done.
    remote: Total 48 (delta 42), reused 0 (delta 0)
    Unpacking objects: 100% (48/48), done.
    From trac.sagemath.org:sage
    * [new branch]      u/user/description -> FETCH_HEAD

The ``u/user/description`` branch is now temporarily (until you fetch
something else) stored in your local git database under the alias
``FETCH_HEAD``. In the second step, we make it available as a new
local branch and switch to it. Your local branch can have a different
name, for example::

    [user@localhost sage]$ git checkout -b my_branch FETCH_HEAD
    Switched to a new branch 'my_branch'

creates a new branch in your local git repository named ``my_branch``
and modifies your local Sage filesystem tree to the state of the files
in that ticket. You can now edit files and commit changes to your
local branch.


.. _section-git-push:

Pushing your changes to a remote
--------------------------------

Push your branch to the remote with either::

    [alice@localhost sage]$ git push --set-upstream origin HEAD:branch-name

if you started the branch yourself and do not follow any other branch, or use::

    [alice@localhost sage]$ git push origin HEAD:u/user/description

if your branch already has an upstream branch.

Here, ``HEAD`` means that you are pushing the most recent commit (and, by
extension, all of its parent commits) of the current local branch to the remote
branch.


.. _section-git-pull:

Getting changes from a remote
-----------------------------

A common task during development is to synchronize your local copy of
the branch with the branch on the GitHub sage repo. In particular, assume you
downloaded the branch of a PR made by someone else, say Bob, and made some suggestions for
improvements on the PR. Now Bob incorporated
your suggestions into his branch, and you want to get the added
changes to complete your review. Assuming that you originally got
your local branch as in :ref:`section-git-checkout`, you can just
issue::

    [bob@localhost sage]$ git pull upstream pull/12345/head
    From https://github.com/sagemath/sage
     * branch                  refs/pull/35608/head -> FETCH_HEAD
    Merge made by the 'ort' strategy.
     src/doc/common/python3.inv          | Bin 98082 -> 131309 bytes
     src/doc/common/update-python-inv.sh |   7 ++++---
     2 files changed, 4 insertions(+), 3 deletions(-)

This command downloads the changes from the branch of the PR and merges
them into your local branch.


.. _section-git-pull-master:

Updating Master
---------------

The ``master`` branch can be updated just like any other branch. However, your
local copy of the master branch should stay **identical** to the GitHub Sage repo master
branch.

If you accidentally added commits to your local copy of ``master``, you must
delete them before updating the branch.

One way to ensure that you are notified of potential problems is to use ``git
pull --ff-only``, which will raise an error if a non-trivial merge would be
required::

    [user@localhost sage]$ git checkout master
    [user@localhost sage]$ git pull --ff-only upstream master

If this pull fails, then something is wrong with the local copy of the
master branch. To switch to the correct Sage master branch, use::

    [user@localhost sage]$ git checkout master
    [user@localhost sage]$ git reset --hard upstream/master


.. _section-git-merge:

Merging and Rebasing
====================

Sometimes, a new version of Sage is released while you work on a git branch.

Let us assume you started ``my_branch`` at commit ``B``. After a while, your
branch has advanced to commit ``Z``, but you updated ``master`` (see
:ref:`section-git-pull-master`) and now your git history looks like this (see
:ref:`section_walkthrough_logs`):

.. CODE-BLOCK:: text

                     X---Y---Z my_branch
                    /
               A---B---C---D master

How should you deal with such changes? In principle, there are two ways:


* **Rebase:** The first solution is to **replay** commits ``X,Y,Z`` atop of the
  new ``master``. This is called **rebase**, and it rewrites your current
  branch:

  .. CODE-BLOCK:: text

      git checkout my_branch
      git rebase -i master

  In terms of the commit graph, this results in:

  .. CODE-BLOCK:: text

                             X'--Y'--Z' my_branch
                            /
               A---B---C---D master

  Note that this operation rewrites the history of ``my_branch`` (see
  :ref:`section-git-rewriting-history`). This can lead to problems if somebody
  began to write code atop of your commits ``X,Y,Z``. It is safe otherwise.

  **Alternatively**, you can rebase ``my_branch`` while updating master at the
  same time (see :ref:`section-git-pull`):

  .. CODE-BLOCK:: text

    git checkout my_branch
    git pull -r master

* **Merging** your branch with ``master`` will create a new commit above the two
  of them:

  .. CODE-BLOCK:: text

      git checkout my_branch
      git merge master

  The result is the following commit graph:

  .. CODE-BLOCK:: text

                     X---Y---Z---W my_branch
                    /           /
               A---B---C-------D master

  - **Pros:** you did not rewrite history (see
    :ref:`section-git-rewriting-history`).The additional commit is then easily
    pushed to the git repository and distributed to your collaborators.

  - **Cons:** it introduced an extra merge commit that would
    not be there had you used rebase.

  **Alternatively**, you can merge ``my_branch`` while updating master at the
  same time (see :ref:`section-git-pull`):

  .. CODE-BLOCK:: text

    git checkout my_branch
    git pull master

**In case of doubt** use merge rather than rebase. There is less risk involved,
and rebase in this case is only useful for branches with a very long history.


.. _section-git-mergetool:

Merge Tools
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

Conflict Resolution
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

This is clearly wrong; Two developers, namely Alice and Bob, decide to
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

    [alice@laptop sage]$ git add fibonacci.py
    [alice@laptop sage]$ git commit -m 'return correct seed values'

He made her changes a PR to the GitHub sage repo and quickly got merged to the ``develop`` branch. Yes, his `fibonacci` function is not yet perfect but is certainly better than the original. Meanwhile, Alice changed the
multiplication to an addition since that is the correct recursion
formula:

.. CODE-BLOCK:: python

    def fibonacci(i):
        """
        Return the `i`-th Fibonacci number
        """
        return fibonacci(i-1) + fibonacci(i-2)

and merged his branch with the latest ``develop`` branch fetched from the GitHub Sage repo::

    [bob@home sage]$ git add fibonacci.py
    [bob@home sage]$ git commit -m 'corrected recursion formula, must be + instead of *'
    [bob@home sage]$ git merge develop
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
        if i > 1:
            return fibonacci(i-1) * fibonacci(i-2)
        return [0, 1][i]
    =======
        return fibonacci(i-1) + fibonacci(i-2)
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

