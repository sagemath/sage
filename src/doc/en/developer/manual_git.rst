.. highlight:: shell-session

.. _chapter-manual-git:

===================================
Using Git with the Sage Trac Server
===================================

.. WARNING::

    **Sage development moved to GitHub in February 2023.** After the transition,
    some parts of this guide (especially those related with `the Sage Trac
    server <https://trac.sagemath.org>`_) became obsolete and need to be
    updated according to the new workflow on GitHub. See our `transition guide
    from Trac to GitHub
    <https://github.com/sagemath/trac-to-github/blob/master/docs/Migration-Trac-to-Github.md>`_
    for the preliminary version of the workflow.

Now we continue our introduction to git from :ref:`chapter-walkthrough`.
We discuss how to push your local changes to a remote repository
so that your changes can be reviewed for inclusion in Sage.

In the following, we assume that you are in the source directory of Sage (``SAGE_ROOT``),
obtained either from a source tarball or by cloning a Sage git repository
such as https://github.com/sagemath/sage.git, as described in the
`README <https://github.com/sagemath/sage/#readme>`_.
In either case, this source directory is actually the main worktree of
a local git repository.

Using the following command, we can see which remote repository or repositories
are associated with this local repository::

    [user@localhost sage]$ git remote -v
    origin      https://github.com/sagemath/sage.git (fetch)
    origin      https://github.com/sagemath/sage.git (push)

.. _section-git-ssh:

Git authentication through SSH
==============================

In order to push changes securely to a remote repository, git uses
public-key cryptography. This section will show you how to set up the
necessary cryptographic keys for Secure Shell (SSH).


Checking whether you have already have suitable SSH keys
--------------------------------------------------------

Follow the instructions in
https://docs.gitlab.com/ee/user/ssh.html#see-if-you-have-an-existing-ssh-key-pair.


Generating your SSH Keys
------------------------

If you don't have suitable SSH keys yet, you can create a key pair
with the ``ssh-keygen`` tool.

Follow either the detailed instructions at
https://docs.gitlab.com/ee/user/ssh.html#generate-an-ssh-key-pair
or the following brief instructions::

    [user@localhost ~]$ ssh-keygen
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


.. _section-trac-ssh-key:

Linking your Public Key to your Trac Account
--------------------------------------------

In order to push your code directly to a branch on the git repository
trac.sagemath.org, the Sage trac server needs to know your public
key. You can upload it in the preferences, that is

1. Go to https://trac.sagemath.org

2. Log in with your trac username/password

3. Click on "Preferences"

4. Go to the "SSH Keys" tab

5. Paste the content of your public key file
   (e.g. ``~/.ssh/id_rsa.pub``)

6. Click on "Save changes"

Note that this does **not** allow you to ssh into any account on trac,
it is only used to authenticate you to the gitolite installation on
trac. You can test that you are being authenticated correctly by
issuing some basic gitolite commands, for example::

    [user@localhost ~]$ ssh git@trac.sagemath.org info
    hello user, this is git@trac running gitolite3 (unknown) on git 1.7.9.5

     R W      sage
    [user@localhost ~]$ ssh git@trac.sagemath.org help
    hello user, this is gitolite3 (unknown) on git 1.7.9.5

    list of remote commands available:

        desc
        help
        info
        perms
        writable

Adding your Public Key for authentication on another server
-----------------------------------------------------------

If you have an account on a lab or department computer that allows you
to log in remotely via SSH, you can now also use your SSH keys to
log in. Just copy the **public** key file (ending in ``.pub``) to
``~/.ssh/authorized_keys`` on the remote computer and make sure that
the file is only read/writeable by yourself. Voila, the next time you
ssh into that machine you don't have to provide your password.


.. _section-git-trac:

The git repository trac.sagemath.org
====================================

The Sage trac server is another git repository for the Sage source tree, it is
served via the ssh protocol. To add it as a remote repository to your local git
repository, use these commands::

    [user@localhost sage]$ git remote add trac https://github.com/sagemath/sagetrac-mirror.git -t master
    [user@localhost sage]$ git remote set-url --push trac git@trac.sagemath.org:sage.git

.. WARNING::

    **Sage development moved to GitHub in February 2023.** After the move, the
    Sage trac server git@trac.sagemath.org:sage.git no longer operates, but all
    branches are available (in read-only mode) on
    https://github.com/sagemath/sagetrac-mirror.git.

Instead of ``trac`` you can use any other name you want, of course.
To verify that it is set up correctly::

    [user@localhost sage]$ git remote -v
    origin      https://github.com/sagemath/sage.git (fetch)
    origin      https://github.com/sagemath/sage.git (push)
    trac        https://github.com/sagemath/sagetrac-mirror.git (fetch)
    trac        git@trac.sagemath.org:sage.git (push)

It is perfectly fine to have multiple remote repositories for git,
think of them as bookmarks. You can then use ``git pull`` to get
changes and ``git push`` to upload your local changes using::

    [user@localhost sage]$ git <push|pull> trac [ARGS]

.. NOTE::

    In the command above we set up the remote to only track the
    ``master`` branch on the trac server (the ``-t master``
    option). This avoids clutter by not automatically downloading all
    branches ever created. But it also means that you will not fetch
    everything that is on trac by default, and you need to explicitly
    tell git which branch you want to get from trac. See the
    :ref:`section-git-checkout` section for examples.

Note that write operations (``push``) use the ssh protocol (specified by the ``git@``
part). For this to work, you need to have a trac account and to set up your ssh public
key as described in `Trac authentication through ssh
<http://doc.sagemath.org/html/en/developer/trac.html#trac-authentication-through-ssh>`_.
Authentication is necessary if you want to upload anything to ensure
that it really is from you.

The above instructions set up the remote to perform read-only operations (``fetch``)
using HTTPS from a mirror of the trac repository instead. The mirror is faster and
more reliable than our git server. However, this configuration is not recommended if
you use VS Code as an IDE.

If you want to use ssh only for both ``fetch`` and ``push``, use the
following commands instead::

    [user@localhost sage]$ git remote add trac git@trac.sagemath.org:sage.git -t master
    [user@localhost sage]$ git remote -v
    origin      https://github.com/sagemath/sage.git (fetch)
    origin      https://github.com/sagemath/sage.git (push)
    trac        git@trac.sagemath.org:sage.git (fetch)
    trac        git@trac.sagemath.org:sage.git (push)

* The Patch buildbot will automatically test your ticket. See :trac:`wiki/patchbot`
  for more information about its features and limitations. Make sure that you
  look at the log, especially if the patch buildbot did not give you
  the green blob.


.. _section-git-checkout:

Checking Out Tickets
--------------------

Trac tickets that are finished or in the process of being worked on
can have a git branch attached to them. This is the "Branch:" field in
the ticket description. The branch name is generally of the form
``u/user/description``, where ``user`` is the name of the user who
made the branch and ``description`` is some free-form short
description (and can include further slashes).

If you want to work with the changes in that remote branch, you must
make a local copy. In particular, git has no concept of directly
working with the remote branch, the remotes are only bookmarks for
things that you can get from/to the remote server. Hence, the first
thing you should do is to get everything from the trac server's branch
into your local repository. This is achieved by::

    [user@localhost sage]$ git fetch trac u/user/description
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

Pushing Your Changes to a Ticket
--------------------------------

To add your local branch to a trac ticket, you should first decide on
a name on the Sage trac repository.

For read/write permissions on git branches, see
:ref:`section-git_trac-branch-names`

In order to avoid name clashes, you can use
``u/your_username/a_description_of_your_branch`` (the description can contain
slashes, but no spaces). Then:

- **Fill** the ``Branch`` field of the trac ticket with that name.

- **Push** your branch to trac with either::

    [user@localhost sage]$ git push --set-upstream trac HEAD:u/user/description

  if you started the branch yourself and do not follow any other branch,
  or use::

    [user@localhost sage]$ git push trac HEAD:u/user/description

  if your branch already has an upstream branch.

Here, ``HEAD`` means that you are pushing the most recent commit (and, by
extension, all of its parent commits) of the current local branch to the remote
branch.

The ``Branch`` field on the trac ticket can appear in red/green. See
:ref:`section-trac-fields` to learn what it means.


.. _section-git-pull:

Getting Changes
---------------

A common task during development is to synchronize your local copy of
the branch with the branch on trac. In particular, assume you
downloaded somebody else's branch and made some suggestions for
improvements on the trac ticket. Now the original author incorporated
your suggestions into his branch, and you want to get the added
changesets to complete your review. Assuming that you originally got
your local branch as in :ref:`section-git-checkout`, you can just
issue::

    [user@localhost sage]$ git pull trac u/user/description
    From trac.sagemath.org:sage
     * branch            u/user/description -> FETCH_HEAD
    Updating 8237337..07152d8
    Fast-forward
     src/sage/tests/cmdline.py      | 3 ++-
     1 file changed, 2 insertions(+), 1 deletions(-)

where now ``user`` is the other developer's trac username and
``description`` is some description that he chose. This command will
download the changes from the originally-used remote branch and merge
them into your local branch. If you haven't published your local
commits yet then you can also rebase them via::

    [user@localhost sage]$ git pull -r trac u/user/description
    From trac.sagemath.org:sage
     * branch            u/user/description -> FETCH_HEAD
    First, rewinding head to replay your work on top of it...
    Applying: my local commit

See :ref:`section-git-merge` section for an in-depth explanation of
merge vs. rebase.

So far, we assumed that there are no conflicts. It is unavoidable in
distributed development that, sometimes, the same location in a source
source file is changed by more than one person. Reconciling these
conflicting edits is explained in the :ref:`section-git_trac-conflict`
section.


.. _section-git-pull-master:

Updating Master
---------------

The ``master`` branch can be updated just like any other branch. However, your
local copy of the master branch should stay **identical** to the trac master
branch.

If you accidentally added commits to your local copy of ``master``, you must
delete them before updating the branch.

One way to ensure that you are notified of potential problems is to use ``git
pull --ff-only``, which will raise an error if a non-trivial merge would be
required::

    [user@localhost sage]$ git checkout master
    [user@localhost sage]$ git pull --ff-only trac master

If this pull fails, then something is wrong with the local copy of the
master branch. To switch to the correct Sage master branch, use::

    [user@localhost sage]$ git checkout master
    [user@localhost sage]$ git reset --hard trac/master


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

Finally, **do nothing unless necessary:** it is perfectly fine for your branch
to be behind ``master``. You can always merge/rebase if/when your branch's name
appears in red on its trac page (see :ref:`section-trac-fields`), or when you
will really need a feature that is only available in the current master.

.. _section-git-mergetool:

Merge Tools
===========

Simple conflicts can be easily solved with git only (see :ref:`section-git_trac-conflict`)

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
