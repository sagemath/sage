.. _chapter-manual-git:

================
Git the Hard Way
================

For beginners to Sage development with no git experience, we recommend using
the Sage development scripts as explained in :ref:`chapter-walk-through`, which
simplify using git and the trac server. However, you can use git
directly to work on Sage if you want to take off the training
wheels. This chapter will tell you how to do so assuming some
basic familiarity with git.

We assume that you have a copy of the Sage git repository, for example
by running::

    [user@localhost ~]$ git clone git://github.com/sagemath/sage.git
    [user@localhost ~]$ cd sage
    [user@localhost sage]$ make


.. _section-git-branch:

Branching Out
=============

A branch is any set of changes that deviates from the current official
Sage tree. Whenever you start developing some new feature or fix a bug
you should first create a new branch to hold the changes. It is easy
to create a new branch, just check out (switch to) the branch from
where you want to start (that is, ``master``) and use the *git
branch* command::

    [user@localhost sage]$ git checkout master
    [user@localhost sage]$ git branch my_new_branch
    [user@localhost sage]$ git checkout my_new_branch
    [user@localhost sage]$ git branch
      master
    * my_new_branch

The asterisk shows you which branch you are on. Without an argument,
the *git branch* command just displays a list of all local branches
with the current one marked by an asterisk. Also note that *git
branch* creates a new branch, but does not switch to it. To avoid
typing the new branch name twice you can use the shortcut ``git
checkout -b my_new_branch`` to create and switch to the new branch in
one command.


.. _section-git-commit:

Commits (Snapshots)
===================

Once you have your own branch feel free to make any changes as you
like. Whenever you have reached your goal, a milestone towards it, or
just feel like you got some work done you should commit your
changes. That is, snapshot the state of all files in the
repository. First, you need to *stage* the changed files, which tells
git which files you want to be part of the next commit::

    ... edit foobar.txt ...

    [user@localhost sage]$ git status
    # On branch my_branch
    # Untracked files:
    #   (use "git add <file>..." to include in what will be committed)
    #
    #       foobar.txt
    nothing added to commit but untracked files present (use "git add" to track)

    [user@localhost sage]$ git add foobar.txt
    [user@localhost sage]$ git status
    # On branch my_branch
    # Changes to be committed:
    #   (use "git reset HEAD <file>..." to unstage)
    #
    #   new file:   foobar.txt
    #

Once you are satisfied with the list of staged files, you create a new
snapshot with the *commit* command::

    [user@localhost sage]$ git commit
    ... editor opens ...
    [my_branch 31331f7] Added the very important foobar text file
     1 file changed, 1 insertion(+)
      create mode 100644 foobar.txt

This will open an editor for you to write your commit message. The
commit message should generally have a one-line description, followed
by an empty line, followed by further explanatory text::

    Added the very important foobar text file

    This is an example commit message. You see there is a one-line
    summary followed by more detailed description, if necessary.

You can then continue working towards your next milestone, make
another commit, repeat until finished. As long as you do not
*checkout* another branch, all commits that you make will be part of
the branch that you created.



.. _section-git-trac:

The Trac Server
===============

The Sage trac server also holds a copy of the Sage repository, it is
served via ssh. To add it as a remote repository to your local git
repository, use the command::

    [user@localhost sage]$ git remote add trac git@trac.sagemath.org:sage.git -t master
    [user@localhost sage]$ git remote -v
    origin      git://github.com/sagemath/sage.git (fetch)
    origin      git://github.com/sagemath/sage.git (push)
    trac        git@trac.sagemath.org:sage.git (fetch)
    trac        git@trac.sagemath.org:sage.git (push)

Instead of ``trac`` you can use any local name you want, of course. It
is perfectly fine to have multiple remote repositories for git, think
of them as bookmarks. You can then use ``git pull`` to get changes and
``git push`` to upload your local changes using::

    [user@localhost sage]$ git <push|pull> trac [ARGS]

.. note::
   
    In the command above we set up the remote to only track the
    ``master`` branch on the trac server (the ``-t master``
    option). This avoids clutter by not automatically downloading all
    branches ever created. But it also means that you will not fetch
    everything that is on trac by default, and you need to explicitly
    tell git which branch you want to get from trac. See the
    :ref:`section-git-checkout` section for examples.

The way we set up the remote here is via ssh authentication (the
``git@`` part), this requires you to have a trac account and to set up
your ssh public key as described in
:ref:`section-trac-ssh-key`. Authentication is necessary if you want
to upload anything to ensure that it really is from you. However, if
you just want to download branches from the trac server then you can
set up the remote to use the git protocol without authentication::

    [user@localhost sage]$ git remote add trac git://trac.sagemath.org/sage.git -t master

Setting up the remote repository this way allows you to perform all
steps covered this manual (except for :ref:`section-git-push`) without
having a trac account. To switch between the two setups, just remove
the current remote repository with ``git remote remove trac`` and then
run the respective ``git remote add trac ...`` command.
     



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
a name on the Sage trac repository. In order to avoid name clashes,
you have push permissions to branches of the form ``u/user/*`` where
``user`` is your trac username and ``*`` is any valid git branch name.
By default, you do *not* have push permissions
to other user's branches or the Sage master branch. In the following,
we will be using ``u/user/description`` as the branch name, where it
is understood that you replaced

* ``user`` with your trac username, and
* ``description`` with some (short but self-explanatory) description of
  your branch. May contain further slashes, but spaces are not allowed.

Your first step should be to put your chosen name into the "Branch:"
field on the trac ticket. To push your branch to trac you then use
either::

    [user@localhost sage]$ git push --set-upstream trac HEAD:u/user/description

if you started the branch yourself and do not follow any other branch,
or use::

    [user@localhost sage]$ git push trac HEAD:u/user/description

if your branch already has an upstream branch.  The ``HEAD`` means
that you are pushing the most recent commit (and, by extension, all of
its parent commits) of the current local branch to the remote
branch.

The ``Branch`` field on the trac ticket page is color coded:
red means there is an issue,
green means it will merge cleanly into ``master``. If it is red, the
tooltip will tell you what is wrong.  If it is green, then it will
link to a diff of the changes against ``master``.



.. _section-git-pull:

Getting Changes
---------------

A common task during development is to synchronize your local copy of
the branch with the branch on trac. In particular, assume you
downloaded somebody else's branch made some suggestions for
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
conflicting edits is explained in the :ref:`section-git-conflict`
section.


.. _section-git-pull-master:

Updating Master
---------------

The ``master`` branch can be updated just like any other
branch. However, you should be take care to keep your local copy of
the master branch identical to the trac master branch, since this is
the current official Sage version. In particular, if you accidentally
added commits to your local copy of the master then you need to delete
those instead of merging them with the official master branch. One way
to ensure that you are notified of potential problems is to use ``git
pull --ff-only``, which will raise an error if a non-trivial merge
would be required::

    [user@localhost sage]$ git checkout master
    [user@localhost sage]$ git pull --ff-only trac master

If this pull fails, then something is wrong with the local copy of the
master branch. To switch to the correct Sage master branch, use::

    [user@localhost sage]$ git checkout master
    [user@localhost sage]$ git reset --hard trac/master


.. _section-git-merge:

Merging and Rebasing
====================

Invariably, Sage development continues while you are working on your
local branch. For example, let us assume you started ``my_branch`` at
commit ``B``. After a while, your branch has advanced to commit ``Z``
while the Sage master branch has advanced to ``D`` ::

                     X---Y---Z my_branch
                    /
               A---B---C---D master

How should you deal with upstream changes while you are
still developing your code? In principle, there are two ways of
dealing with it:


* The first solution is to change the commits in your local branch to
  start out at the new master. This is called **rebase**, and it
  rewrites your current branch::
   
      git checkout my_branch
      git rebase master

  Here, we assumed that ``master`` is your local and up-to-date copy
  of the master branch. Alternatively, you can pull changes from the
  trac server and rebase the current in one go with the combination
  ``git pull -r master`` command, see :ref:`section-git-pull`. In
  terms of the commit graph, this results in::

                             X'--Y'--Z' my_branch
                            /
               A---B---C---D master

  Since the SHA1 hash includes the hash of the parent, all commits
  change. This means that you should only ever use rebase if nobody
  else has used one of your ``X``, ``Y``, ``Z`` commits to base their
  development on. 


* The other solution is to not change any commits, and instead create
  a new merge commit ``W`` which merges in the changes from the newer
  master. This is called **merge**, and it merges your current branch
  with another branch::

      git checkout my_branch
      git merge master

  Here, we assumed that ``master`` is your local and up-to-date copy
  of the master branch. Alternatively, you can pull changes from the
  trac server and merge them into the current branch with the
  combination ``git pull master`` command, see
  :ref:`section-git-pull`. The result is the following commit graph::

                     X---Y---Z---W my_branch
                    /           /
               A---B---C-------D master

  The downside is that it introduced an extra merge commit that would
  not be there had you used rebase. But that is also the advantage of
  merging: None of the existing commits is changed, only a new commit
  is made. This additional commit is then easily pushed to the git
  repository and distributed to your collaborators.


As a general rule of thumb, use merge if you are in doubt. The
downsides of rebasing can be really severe for other developers, while
the downside of merging is just minor. Finally, and perhaps the most
important advice, do nothing unless necessary. It is perfectly fine
for your branch to be behind the master branch. Just keep developing
your feature. Trac will tell you if it doesn't merge cleanly with the
current master by the color of the "Branch:" field, and the patchbot
(coloured blob on the trac ticket) will test whether your branch still
works on the current master. Unless either a) you really need a
feature that is only available in the current master, or b) there is a
conflict with the current master, there is no need to do anything on
your side.


.. _section-git-conflict:

Conflict Resolution
===================

Merge conflicts happen if there are overlapping edits, and they are an
unavoidable consequence of distributed development. Fortunately,
resolving them is common and easy with git. As a hypothetical example,
consider the following code snippet::

    def fibonacci(i):
        """
        Return the `i`-th Fibonacci number
        """
        return fibonacci(i-1) * fibonacci(i-2)

This is clearly wrong; Two developers, namely Alice and Bob, decide to
fix it. First, in a cabin in the woods far away from any internet
connection, Alice corrects the seed values::

    def fibonacci(i):
       """
       Return the `i`-th Fibonacci number
       """
       if i > 1:
           return fibonacci(i-1) * fibonacci(i-2)
       return [0, 1][i]

and turns those changes into a new commit::

    [alice@laptop]$ git commit -m 'return correct seed values'
    [fibonacci_alice 14ae1d3] return correct seed values
     1 file changed, 3 insertions(+), 1 deletion(-)

However, not having an internet connection, she cannot immediately
send her changes to the trac server. Meanwhile, Bob changes the
multiplication to an addition since that is the correct recursion
formula::

    def fibonacci(i):
        """
        Return the `i`-th Fibonacci number
        """
        return fibonacci(i-1) + fibonacci(i-2)

and immediately uploads his change::

    [bob@home]$ git commit -m 'corrected recursion formula, must be + instead of *'
    [fibonacci_bob 41675df] corrected recursion formula, must be + instead of *
    1 file changed, 1 insertion(+), 1 deletion(-)

    [bob@home]$ git push trac HEAD:u/bob/fibonacci
    Counting objects: 5, done.
    Delta compression using up to 8 threads.
    Compressing objects: 100% (2/2), done.
    Writing objects: 100% (3/3), 320 bytes | 0 bytes/s, done.
    Total 3 (delta 1), reused 0 (delta 0)
    To trac.sagemath.org:sage
       14afe53..41675df  HEAD -> u/bob/fibonacci

Eventually, Alice returns to civilization. In her mailbox, she finds a
trac notification email that Bob has uploaded further changes to their
joint project. Hence, she starts out by getting his changes into her
own local branch::

    [alice@laptop]$ git pull trac u/bob/fibonacci
    From trac.sagemath.org:sage
     * branch            u/bob/fibonacci     -> FETCH_HEAD
    Auto-merging fibonacci.py
    CONFLICT (content): Merge conflict in fibonacci.py
    Automatic merge failed; fix conflicts and then commit the result.

.. skip    # doctester confuses >>> with input marker

The file now looks like this::

    def fibonacci(i):
        """
        Return the `i`-th Fibonacci number
        """
    <<<<<<< HEAD
        if i > 1:
            return fibonacci(i-1) * fibonacci(i-2)
        return i
    =======
        return fibonacci(i-1) + fibonacci(i-2)
    >>>>>>> 41675dfaedbfb89dcff0a47e520be4aa2b6c5d1b

The conflict is shown between the conflict markers ``<<<<<<<`` and
``>>>>>>>``. The first half (up to the ``=======`` marker) is Alice's
current version, the second half is Bob's version. The 40-digit hex
number after the second conflict marker is the SHA1 hash of the most
recent common parent of both.

It is now Alice's job to resolve the conflict by reconciling their
changes, for example by editing the file. Her result is::
    
    def fibonacci(i):
        """
        Return the `i`-th Fibonacci number
        """
        if i > 1:
            return fibonacci(i-1) + fibonacci(i-2)
        return [0, 1][i]
    
And then upload both her original change *and* her merge commit to trac::

    [alice@laptop]$ git commit -m "merged Bob's changes with mine"
    [fibonacci_allice 6316447] merged Bob's changes with mine
    $ git push trac HEAD:u/alice/fibonacci

The resulting commit graph now has a loop::
    
    $ git log --graph --oneline
    *   6316447 merged Bob's changes with mine
    |\  
    | * 41675df corrected recursion formula, must be + instead of *
    * | 14ae1d3 return correct seed values
    |/  
    * 14afe53 initial commit
    
If Bob decides to do further work on the ticket then he will have to
pull from ``u/alice/fibonacci``. However, this time there is no
conflict on his end: git downloads both Alice's conflicting commit and
her resolution.


Merge Tools
-----------

Just editing the file with the conflict markers is often the simplest
solution. However, for more complicated conflicts there is a range of
specialized programs available to help you identify the
conflicts. Because the conflict marker includes the hash of the most
recent common parent, you can use a three-way diff::

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
    
If you don't have a favorite merge tool we suggest you try meld
(cross-platform). The result looks like the following screenshot.

.. image:: static/meld-screenshot.png

The middle file is the most recent common parent; on the right is
Bob's version and on the left is Alice's conflicting version. Clicking
on the arrow moves the marked change to the file in the adjacent
pane. 
