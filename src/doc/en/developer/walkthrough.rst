.. highlight:: shell-session

.. _chapter-walkthrough:

========================
Development Walk-through
========================

This section is a concise overview of the Sage development process. We will see
how to make changes to the Sage source code and record them in the Git revision
control system.

In the sections of the following chapter :ref:`section-development-on-github`,
we will look at communicating these changes back to the Sage project. All
changes to Sage source code have to go through `the Sage repository
<https://github.com/sagemath/sage>`_ on GitHub.

For examples, we assume your name Alice. Always replace it with your own name.

.. _section-walkthrough-setup-git:

Checking Git
============

First, open a shell (for instance, Terminal on Mac) and check that Git works::

    [alice@localhost ~]$ git
    usage: git [--version] [--help] [-C <path>] [-c name=value]
    ...
    The most commonly used git commands are:
       add        Add file contents to the index
    ...
       tag        Create, list, delete or verify a tag object signed with GPG

    'git help -a' and 'git help -g' lists available subcommands and some
    concept guides. See 'git help <command>' or 'git help <concept>'
    to read about a specific subcommand or concept.

Don't worry about the giant list of subcommands. You really only need a handful
of them for effective development, and we will walk you through them in this
guide. If you got a "command not found" error, then you don't have Git
installed; now is the time to install it. See
:ref:`section-git-install` for instructions.

Because we also track who does what changes with Git, you must tell
Git how you want to be known. Check if Git knows you::

    [alice@localhost ~]$ git config --global user.name
    Alice Adventure
    [alice@localhost ~]$ git config --global user.email
    alice@wonderland.com

If you have multiple computers, then use the same name on each of them.  This
name/email combination ends up in commits. So if it's not set yet, do it now
before you forget! This only needs to be done once. See
:ref:`section-git-setup-name` for instructions.

.. _section-walkthrough-sage-source:

Obtaining the Sage source code
==============================

Obviously one needs the Sage source code to develop. You can use your
local installation of Sage (if you installed Sage from source), or
(to start from scratch) download it from our Sage repository on GitHub::

    [alice@localhost ~]$ git clone --origin upstream https://github.com/sagemath/sage.git
    Cloning into 'sage'...
    [...]
    Checking connectivity... done.

This creates a directory named ``sage`` containing the sources for the
current stable and development releases of Sage. You next need to switch
to the develop branch (latest development release)::

    [alice@localhost ~]$ cd sage
    [alice@localhost sage]$ git checkout develop

Next, build Sage, following the instruction in the file `README.md
<https://github.com/sagemath/sage/#readme>`_ in ``SAGE_ROOT``. If all
prerequisites to build are in place, the commands ``./configure && make -j4``
will do it.  Additional details can be found in the section on `installation
from source <../installation/source.html>`_ in the Sage installation guide. If
you wish to use conda-forge, see the section on `conda
<../installation/conda.html>`_.

.. NOTE::

    macOS allows changing directories without using exact capitalization.
    Beware of this convenience when compiling for macOS. Ignoring exact
    capitalization when changing into :envvar:`SAGE_ROOT` can lead to build
    errors for dependencies requiring exact capitalization in path names.


.. _section-walkthrough-branch:

Branching out
=============

In order to start modifying Sage, we want to make a new *branch* in the local
Sage repo. A branch is a copy (except that it doesn't take up twice the space)
of the Sage source code where you can store your modifications to the Sage
source code (and which you can push to your fork of the Sage repository on GitHub).

To begin with, type the command ``git branch``. You will see the following::

    [alice@localhost sage]$ git branch
    * develop
      master

The asterisk shows you which branch you are on. Without an argument,
the ``git branch`` command displays a list of all local branches
with the current one marked by an asterisk.

It is easy to create a new branch. First make sure you are on the branch from
which you want to branch out. That is, if you are not currently on the
``develop`` branch, type the command ``git checkout develop``::

    [alice@localhost sage]$ git checkout develop
    Switched to branch 'develop'
    Your branch is up-to-date with 'origin/develop'.

Then use the ``git branch`` command to create a new branch, as follows::

    [alice@localhost sage]$ git branch last_twin_prime

Also note that ``git branch`` creates a new branch, but does not switch
to it. For this, you have to use ``git checkout``::

    [alice@localhost sage]$ git checkout last_twin_prime
    Switched to branch 'last_twin_prime'

Now if you use the command ``git branch``, you will see the following::

    [alice@localhost sage]$ git branch
      develop
    * last_twin_prime
      master

Note that unless you explicitly push a branch to a remote Git repository, the
branch is a local branch that is only on your computer and not visible to
anyone else.

To avoid typing the new branch name twice you can use the shortcut
``git checkout -b last_twin_prime develop`` to create and switch to the new
branch based on ``develop`` in one command.


.. _section_walkthrough_logs:

The history
===========

It is always a good idea to check that you are making your edits on the branch
that you think you are on. The following command shows you the topmost commit
in detail, including its changes to files::

    [alice@localhost sage]$ git show

To dig deeper, you can inspect the log::

    [alice@localhost sage]$ git log

By default, this lists all commits in reverse chronological order.

- If you find your branch to be in the wrong place, see the
  :ref:`section-git-recovery` section.

- Many tools are available to help you visualize the history tree better.
  For instance, ``tig`` is a very nice text-mode tool.

.. _section-walkthrough-add-edit:

Editing the source code
=======================

Once you have your own branch, feel free to make any changes to source files as
you like. The chapter :ref:`section-writing-code-for-sage` explains how your
code should look like to fit into Sage, and how we ensure high code quality
throughout.

The Git command ``git status`` is probably the most important of all Git
commands. It tells you which files changed, and how to continue with recording
the changes::

    [alice@localhost sage]$ git status
    On branch last_twin_prime
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

    [alice@localhost sage]$ git diff some_file.py

to show you the differences.


.. _section-walkthrough-make:

Rebuilding Sage
===============

Once you have made any changes, you of course want to build Sage and try out
your edits. As long as you only modified the Sage library (that is, Python and
Cython files under ``src/sage/...``) you just have to run::

    [alice@localhost sage]$ ./sage -br

to rebuild the Sage library and then start Sage.

.. NOTE::

    All changes to Python files take effect immediately after restarting Sage
    (unless you have used ``./configure --disable-editable`` when you built
    Sage). Hence you can just start Sage instead of ``./sage -br`` if only Python
    files were modified.

If you made changes to :ref:`third-party packages <chapter-packaging>`
installed as part of Sage, then you have to run ::

    [alice@localhost sage]$ make build

as if you were `installing Sage from scratch
<http://doc.sagemath.org/html/en/installation/source.html>`_.  However, this
time only, the packages which were changed (or which depend on a changed package)
will be rebuilt, so it should be much faster than building Sage the first
time.

.. NOTE::

    If you have `pulled a branch from the GitHub Sage repository
    <https://doc.sagemath.org/html/en/developer/workflows.html#checking-out-an-existing-pr>`_,
    it may depend on changes to third-party packages, so ``./sage -br`` may
    fail.  If this happens (and you believe the code in this branch should
    compile), try running ``make build``.

Rarely there are conflicts with other packages,
or with the already-installed older version of the package that you
changed, in that case you do have to recompile everything using::

    [alice@localhost sage]$ make distclean && make build

Also, don't forget to run the tests (see :ref:`chapter-doctesting`)
and build the documentation (see :ref:`chapter-sage_manuals`).

.. NOTE::

    If you switch between branches based on different releases, the timestamps
    of modified files will change. This triggers recythonization and recompilation
    of modified files on subsequent builds, whether or not you have made any
    additional changes to files. To minimize the impact of switching between branches,
    install ccache using the command ::

        [alice@localhost sage]$ ./sage -i ccache

    Recythonization will still occur when rebuilding, but the recompilation stage
    first checks whether previously compiled files are cached for reuse before
    compiling them again. This saves considerable time rebuilding.


.. _section-walkthrough-commit:

Making commits
==============

Whenever you have reached your goal, a milestone towards it, or
just feel like you got some work done you should *commit* your
changes. A commit is just a snapshot of the state of all files in
the repository.

Unlike with some other revision control programs, in Git you first
need to *stage* the changed files, which tells Git which files you
want to be part of the next commit::

    [alice@localhost sage]$ git status
    # On branch my_branch
    # Untracked files:
    #   (use "git add <file>..." to include in what will be committed)
    #
    #       src/sage/primes/last_pair.py
    nothing added to commit but untracked files present (use "git add" to track)

    [alice@localhost sage]$ git add src/sage/primes/last_pair.py
    [alice@localhost sage]$ git status
    # On branch my_branch
    # Changes to be committed:
    #   (use "git reset HEAD <file>..." to unstage)
    #
    #   new file:   src/sage/primes/last_pair.py
    #

Once you are satisfied with the list of staged files, you create a new
snapshot with the ``git commit`` command::

    [alice@localhost sage]$ git commit
    ... editor opens ...
    [my_branch 31331f7] Added the very important foobar text file
     1 file changed, 1 insertion(+)
      create mode 100644 foobar.txt

This will open an editor for you to write your commit message. The
commit message should generally have a one-line description, followed
by an empty line, followed by further explanatory text:

.. CODE-BLOCK:: text

    Added the last twin prime

    This is an example commit message. You see there is a one-line
    summary followed by more detailed description, if necessary.

You can then continue working towards your next milestone, make
another commit, repeat until finished. As long as you do not
``git checkout`` another branch, all commits that you make will be part of
the branch that you created.

