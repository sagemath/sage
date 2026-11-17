.. highlight:: shell-session

.. _chapter-walkthrough:

========================
Development Walk-through
========================

This section is a concise overview of the Sage development process. We will see
how to make changes to the Sage source code and record them in the Git revision
control system.

.. _section-quick-start:

Quick start
===========

If you are in a hurry, you can skip the details and just follow these steps:

1. Install Git (see :ref:`section-git-install`) and
   `Pixi <https://pixi.prefix.dev/latest/installation/>`_.

2. Clone the Sage repository from GitHub

   .. code-block:: console

       $ git clone --origin upstream https://github.com/sagemath/sage.git

3. Change into the directory

   .. code-block:: console

       $ cd sage

4. Install the locked development dependencies

   .. code-block:: console

       $ pixi install --locked

   Pixi selects the environment for your operating system automatically.
   See `the Pixi setup instructions <../installation/pixi.html>`_ for platform
   prerequisites, including the compiler setup for native Windows builds.

5. Build and install Sage

   .. code-block:: console

       $ pixi run build

   Start Sage with ``pixi run sage``. To test a file you changed, use
   ``pixi run test path/to/file.py``.

6. Create a new branch for your changes

   .. code-block:: console

       $ git checkout -b my_branch develop

7. Make your changes, and push them to your fork on GitHub

   .. code-block:: console

       $ git add .
       $ git commit -m "Your commit message here"
       $ git push origin my_branch

8. Create a pull request on GitHub to merge your changes into the Sage repository.

.. _section-walkthrough-setup-git:

Checking Git
============

First, open a shell (for instance, Terminal on Mac) and check that Git works

.. code-block:: console

    $ git --version
    git version xyz

If you got a "command not found" error, then you don't have Git
installed; now is the time to install it. See
:ref:`section-git-install` for instructions.

Because we also track who does what changes with Git, you must tell
Git how you want to be known. Check if Git knows you

.. code-block:: console

    $ git config --global user.name
    Alice Adventure
    $ git config --global user.email
    alice@wonderland.com

If you see your name and email address, then you are all set. This
name/email combination ends up in commits. So if it's not set yet, do it now
before you forget! This only needs to be done once. See
:ref:`section-git-setup-name` for instructions.

.. _section-walkthrough-sage-source:

Obtaining the Sage source code
==============================

Obviously one needs the Sage source code to develop. You can download it
from our Sage repository on GitHub

.. code-block:: console

    $ git clone --origin upstream https://github.com/sagemath/sage.git
    Cloning into 'sage'...
    $ cd sage

This creates a directory named ``sage`` containing the most recent version of
the Sage source code.

Building Sage
=============

Sage is a large project with many dependencies. For development, we recommend
using Pixi to install the dependencies from conda-forge and build your checkout.
You can instead use the approach described in the section
`Sage-the-distribution (classical build) <../installation/source-distro.html>`_
of the Sage installation guide, or use ``uv``, see `System Package Manager <../installation/source.html#using-system-package-manager>`_.

Install Pixi by following its `official instructions
<https://pixi.prefix.dev/latest/installation/>`_. Then, from your Sage checkout,
run:

.. code-block:: console

    $ pixi install --locked
    $ pixi run build
    $ pixi run sage

The first command creates a local environment in ``.pixi/`` using the versions
recorded in ``pixi.lock``. The build task installs Sage in editable mode, so
changes to the source are used by subsequent Sage sessions. Cython changes
are rebuilt when needed.

To work in an activated shell, use ``pixi shell`` and then run ``sage`` or other
development commands directly. Use ``exit`` to leave that shell. Run tests of
an individual file with:

.. code-block:: console

    $ pixi run test src/sage/rings/integer.pyx

The default environment uses Python 3.13. To build and run with Python 3.14,
select the ``py314`` environment:

.. code-block:: console

    $ pixi run -e py314 build
    $ pixi run -e py314 sage

Each environment has its own dependencies and build directory. Rerun the build
task after changing Python package metadata or updating the environment.

See `Developing with Pixi <../installation/pixi.html>`_ for platform setup,
dependency updates, and build-directory details. Other build options are
described in `Building from source <../installation/source.html>`_.

.. _section-walkthrough-branch:

Branching out
=============

In order to start modifying Sage, we want to make a new *branch* in the local
Sage repo. A branch is a copy (except that it doesn't take up twice the space)
of the Sage source code where you can store your modifications to the Sage
source code (and which you can push to your fork of the Sage repository on GitHub).

To begin with, type the command ``git branch``. You will see the following

.. code-block:: console

    $ git branch
    * develop
      master

The asterisk shows you which branch you are on. Without an argument,
the ``git branch`` command displays a list of all local branches
with the current one marked by an asterisk.

It is easy to create a new branch, as follows

.. code-block:: console

    $ git checkout -b last_twin_prime develop

This will create a new branch named ``last_twin_prime`` based on
the ``develop`` branch and switch to it.

Now if you use the command ``git branch``, you will see the following

.. code-block:: console

    $ git branch
      develop
    * last_twin_prime
      master

Note that unless you explicitly push a branch to a remote Git repository, the
branch is a local branch that is only on your computer and not visible to
anyone else.

.. _section-walkthrough-add-edit:

Editing the source code
=======================

Once you have your own branch, feel free to make any changes to source files as
you like. The chapter :ref:`section-writing-code-for-sage` explains how your
code should look like to fit into Sage, and how we ensure high code quality
throughout.

The Git command ``git status`` is probably the most important of all Git
commands. It tells you which files changed, and how to continue with recording
the changes

.. code-block:: console

    $ git status
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

To dig deeper into what was changed in the files you can use

.. code-block:: console

    $ git diff some_file.py

to show you the differences.


.. _section-walkthrough-testing:

Testing changes
===============

Once you have made any changes, you of course want to try out
your edits. All changes to Python and Cython files take effect immediately
after restarting Sage, so there is no need to explicitly rebuild Sage.

The changes can be tested by running Sage and verifying that the modifications
work as expected. For example, if you modified a function, you can call it
directly in Sage to ensure it behaves as intended.

Additionally, you can write or modify doctests in the relevant files to
confirm the correctness of your changes.
To run the doctests for a specific file, use the following command::

    $ pixi run test path/to/your/file.py

This will execute all the doctests in the specified file and report any
failures. Make sure all tests pass before proceeding
(see :ref:`chapter-doctesting` for more details).
Also, don't forget to build the documentation (see :ref:`chapter-sage_manuals`).

.. _section-walkthrough-commit:

Making commits
==============

Whenever you have reached your goal, a milestone towards it, or
just feel like you got some work done you should *commit* your
changes. A commit is just a snapshot of the state of all files in
the repository.

You first need to *stage* the changed files, which tells Git which files you
want to be part of the next commit

.. code-block:: console

    $ git status
    On branch last_twin_prime
    Untracked files:
      (use "git add <file>..." to include in what will be committed)
          src/sage/primes/last_pair.py
    nothing added to commit but untracked files present (use "git add" to track)

    $ git add src/sage/primes/last_pair.py
    $ git status
    On branch last_twin_prime
    Changes to be committed:
      (use "git reset HEAD <file>..." to unstage)
      new file:   src/sage/primes/last_pair.py

Once you are satisfied with the list of staged files, you create a new
snapshot with the ``git commit`` command

.. code-block:: console

    $ git commit
    ... editor opens ...
    [last_twin_prime 31331f7] Added the very important foobar text file
     1 file changed, 1 insertion(+)
      create mode 100644 foobar.txt

This will open an editor for you to write your commit message. The
commit message should generally have a one-line description, followed
by an empty line, followed by further explanatory text:

.. code-block:: text

    Added the last twin prime

    This is an example commit message. You see there is a one-line
    summary followed by more detailed description, if necessary.

You can then continue working towards your next milestone, make
another commit, repeat until finished. As long as you do not
``git checkout`` another branch, all commits that you make will be part of
the branch that you created.

.. tip::

    Consider activating the pre-commit hook so that your code is automatically
    checked before each commit:

    .. code-block:: console

        $ git config core.hooksPath git-hooks

    This only needs to be done once after cloning the repository.

Open pull request
=================

Once you are happy with your changes, you can propose these for review and
integration into the main project.
The first step is to push your branch to your fork of the `the Sage repository
<https://github.com/sagemath/sage>`_ on GitHub. This is done with the command

.. code-block:: console

    $ git push origin last_twin_prime

Now you can go `to GitHub and create a pull request
<https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request-from-a-fork>`_.
See :ref:`chapter-workflows` for more details on the workflow of
creating a pull request and the review process.
