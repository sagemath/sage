.. _chapter-walkthrough:

========================
Sage Development Process
========================

This section is a concise overview of the Sage development process. In
it, we will see how to make changes to the Sage source code and record
them in the git revision control system.

In the following section on :ref:`chapter-git_trac` we will look at
communicating these changes back to the Sage project.  We also have a handy
`one-page "cheat sheet"
<http://github.com/sagemath/git-trac-command/raw/master/doc/git-cheat-sheet.pdf>`_
of commonly used git commands that you can print out and leave on your
desk.  We have some :ref:`recommended references and tutorials
<section-git-tutorials>` as well.

You can alternatively fork and create a pull request at
`github <http://github.com/sagemath/sage>`_ which will automatically fetch
your code and open a ticket on our trac server.


.. _section-walkthrough-setup-git:

Configuring Git
===============

One way or another, ``git`` is what Sage uses for tracking changes.
So first, open a shell (for instance, Terminal on Mac) and check that
``git`` works::

    [user@localhost]$ git
    usage: git [--version] [--help] [-C <path>] [-c name=value]
    ...
    The most commonly used git commands are:
       add        Add file contents to the index
    ...
       tag        Create, list, delete or verify a tag object signed with GPG

    'git help -a' and 'git help -g' lists available subcommands and some
    concept guides. See 'git help <command>' or 'git help <concept>'
    to read about a specific subcommand or concept.


Don't worry about the giant list of subcommands. You really only need
a handful for effective development, and we will walk you through them
in this guide. If you got a "command not found" error, then you don't
have git installed. Now is the time to install it; see
:ref:`chapter-git-setup` for instructions.

Because we also track who does changes in Sage with git, you must tell
git how you want to be known. This only needs to be done once::

    [user@localhost]$ git config --global user.name "Your Name"
    [user@localhost]$ git config --global user.email you@yourdomain.example.com

If you have multiple accounts / computers use the same name on each of
them. This name/email combination ends up in commits, so do it now
before you forget!


.. _section-walkthrough-sage-source:

Obtaining the Sage Source Code
==============================

Obviously one needs the Sage source code to develop.  You can use your
local installation of Sage, or (to start without Sage) download it
from github which is a public read-only mirror (=faster) of our
internal git repository::

    [user@localhost ~]$ git clone git://github.com/sagemath/sage.git
    Cloning into 'sage'...
    [...]
    Checking connectivity... done.

This creates a directory named ``sage`` containing the sources for the
current stable and development releases of Sage. You next need to switch
to the develop branch (latest development release)::

    [user@localhost ~]$ cd sage
    [user@localhost sage]$ git checkout develop

You will then need to `compile Sage
<http://doc.sagemath.org/html/en/installation/source.html>`_ in order to use it. If
you cloned, you will need to remain on the internet for it to download various
packages of Sage::

    [user@localhost sage]$ make

.. NOTE::

    Mac OS X allows changing directories without using exact capitalization.
    Beware of this convenience when compiling for OS X. Ignoring exact
    capitalization when changing into :envvar:`SAGE_ROOT` can lead to build
    errors for dependencies requiring exact capitalization in path names.

For the experts, note that the repository at
`git.sagemath.org <http://git.sagemath.org>`_ is where development
actually takes place.


.. _section-walkthrough-branch:

Branching Out
=============

In order to start modifying Sage, we want to make a *branch* of Sage.
A branch is a copy (except that it doesn't take up twice the space) of
the Sage source code where you can store your modifications to the
Sage source code and which you can upload to trac tickets.

It is easy to create a new branch, just check out (switch to) the branch
from where you want to start (that is, ``master``) and use the ``git
branch`` command::

    [user@localhost sage]$ git checkout master
    [user@localhost sage]$ git branch last_twin_prime
    [user@localhost sage]$ git checkout last_twin_prime

You can list all branches using::

    [user@localhost]$ git branch
      master
    * last_twin_prime

The asterisk shows you which branch you are on. Without an argument,
the ``git branch`` command just displays a list of all local branches
with the current one marked by an asterisk. Also note that ``git
branch`` creates a new branch, but does not switch to it. For this,
you have to use ``git checkout``::

    [user@localhost sage]$ git checkout master
    Switched to branch 'master'
    Your branch is up-to-date with 'github/master'.
    [user@localhost sage]$ git branch
    * master
      last_twin_prime
    [user@localhost sage]$ git checkout last_twin_prime
    Switched to branch 'last_twin_prime'

Note that, unless you explicitly upload ("push") a branch to remote
git repository, the local branch will only be on your computer and not
visible to anyone else.

To avoid typing the new branch name twice you can use the shortcut
``git checkout -b my_new_branch`` to create and switch to the new
branch in one command.



.. _section_walkthrough_logs:

The History
===========

It is always a good idea to check that you are making your edits on
the version that you think you are on. The first one shows you the
topmost commit in detail, including its changes to the sources::

    [user@localhost sage]$ git show

To dig deeper, you can inspect the log::

    [user@localhost sage]$ git log

By default, this lists all commits in reverse chronological order.

- If you find your branch to be in the wrong place, see the
  :ref:`section-git-recovery` section.

- Many programs are available to help you visualize the history tree
  better. ``tig`` is a very nice text-mode such tool.

.. _section-walkthrough-add-edit:

Editing the Source Code
=======================

Once you have your own branch, feel free to make any changes as you
like. :ref:`Subsequent chapters <section-writing-code-for-sage>` of
this developer guide explain how your code should look like to fit
into Sage, and how we ensure high code quality throughout.

*Status* is probably the most important git command. It tells
you which files changed, and how to continue with recording the
changes::

    [user@localhost sage]$ git status
    On branch master
    Changes not staged for commit:
      (use "git add <file>..." to update what will be committed)
      (use "git checkout -- <file>..." to discard changes in working directory)

        modified:   some_file.py
        modified:   src/sage/primes/all.py

    Untracked files:
      (use "git add <file>..." to include in what will be committed)

        src/sage/primes/last_pair.py

    no changes added to commit (use "git add" and/or "git commit -a")

To dig deeper into what was changed in the files you can use::

    [user@localhost sage]$ git diff some_file.py

to show you the differences.



.. _section-walkthrough-make:

Rebuilding Sage
===============

Once you have made any changes you of course want to build Sage and
try out your edits. As long as you only modified the Sage library
(that is, Python and Cython files under ``src/sage/...``) you just
have to run::

    [user@localhost sage]$ ./sage -br

to rebuild the Sage library and then start Sage. This should be quite
fast. If you made changes to
:ref:`third-party packages <chapter-packaging>`, then you have to run ::

    [user@localhost sage]$ make

as if you were `installing Sage from scratch
<http://doc.sagemath.org/html/en/installation/source.html>`_.
However, this time only packages which were changed (or which depend
on a changed package) will be recompiled,
so it shoud be much faster than compiling Sage
the first time. Rarely there are conflicts with other packages,
or with the already-installed older version of the package that you
changed, in that case you do have to recompile everything using::

    [user@localhost sage]$ make distclean && make

Also, don't forget to run the tests (see :ref:`chapter-doctesting`)
and build the documentation (see :ref:`chapter-sage_manuals`).


.. _section-walkthrough-commit:

Commits (Snapshots)
===================

Whenever you have reached your goal, a milestone towards it, or
just feel like you got some work done you should *commit* your
changes. A commit is just a snapshot of the state of all files in
the *repository* (the program you are working on).

Unlike with some other revision control programs, in git you first
need to *stage* the changed files, which tells git which files you
want to be part of the next commit::

    [user@localhost sage]$ git status
    # On branch my_branch
    # Untracked files:
    #   (use "git add <file>..." to include in what will be committed)
    #
    #       src/sage/primes/last_pair.py
    nothing added to commit but untracked files present (use "git add" to track)

    [user@localhost sage]$ git add src/sage/primes/last_pair.py
    [user@localhost sage]$ git status
    # On branch my_branch
    # Changes to be committed:
    #   (use "git reset HEAD <file>..." to unstage)
    #
    #   new file:   src/sage/primes/last_pair.py
    #

Once you are satisfied with the list of staged files, you create a new
snapshot with the ``git commit`` command::

    [user@localhost sage]$ git commit
    ... editor opens ...
    [my_branch 31331f7] Added the very important foobar text file
     1 file changed, 1 insertion(+)
      create mode 100644 foobar.txt

This will open an editor for you to write your commit message. The
commit message should generally have a one-line description, followed
by an empty line, followed by further explanatory text::

    Added the last twin prime

    This is an example commit message. You see there is a one-line
    summary followed by more detailed description, if necessary.

You can then continue working towards your next milestone, make
another commit, repeat until finished. As long as you do not
``git checkout`` another branch, all commits that you make will be part of
the branch that you created.





