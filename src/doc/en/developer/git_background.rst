.. highlight:: shell-session

.. _chapter-git-background:

=======================
Git Tips and References
=======================

This chapter contains additional material about the Git revision
control system.  See :ref:`chapter-git-setup` for the minimal
steps needed for Sage development.


.. _section-git-configuration:

Configuration tips
==================

Your personal Git configurations are saved in the ``~/.gitconfig``
file in your home directory. Here is an example:

.. code-block:: text

    [user]
        name = Alice Adventure
        email = alice@wonderland.com

    [core]
        editor = emacs

You can edit this file directly or you can use Git to make changes for
you

.. code-block:: console

    $ git config --global user.name "Alice Adventure"
    $ git config --global user.email alice@wonderland.com
    $ git config --global core.editor vim


Aliases
-------

Aliases are personal shortcuts for Git commands. For example, you
might want to be able to shorten ``git checkout`` to ``git co``.  Or
you may want to alias ``git diff --color-words`` (which gives a nicely
formatted output of the diff) to ``git wdiff``. You can do this with

.. code-block:: console

    $ git config --global alias.ci "commit -a"
    $ git config --global alias.co checkout
    $ git config --global alias.st "status -a"
    $ git config --global alias.stat "status -a"
    $ git config --global alias.br branch
    $ git config --global alias.wdiff "diff --color-words"

The above commands will create an ``alias`` section in your ``.gitconfig``
file with contents like this:

.. code-block:: text

    [alias]
        ci = commit -a
        co = checkout
        st = status -a
        stat = status -a
        br = branch
        wdiff = diff --color-words


Editor
------

To set the editor to use for editing commit messages, you can use

.. code-block:: console

    $ git config --global core.editor vim

or set the ``EDITOR`` environment variable.


Merging
-------

To enforce summaries when doing merges (``~/.gitconfig`` file again):

.. code-block:: text

    [merge]
        log = true

Or from the command line

.. code-block:: console

    $ git config --global merge.log true


.. _section-fancy-log:

Fancy log output
----------------

Here is an alias to get a fancy log output. It should go in the
``alias`` section of your ``.gitconfig`` file:

.. code-block:: text

    lg = log --graph --pretty=format:'%Cred%h%Creset -%C(yellow)%d%Creset %s %Cgreen(%cr) %C(bold blue)[%an]%Creset' --abbrev-commit --date=relative

Using this ``lg`` alias gives you the changelog with a colored ASCII graph

.. code-block:: console

    $ git lg
    * 6d8e1ee - (HEAD, origin/my-fancy-feature, my-fancy-feature) NF - a fancy file (45 minutes ago) [Matthew Brett]
    *   d304a73 - (origin/placeholder, placeholder) Merge pull request #48 from hhuuggoo/master (2 weeks ago) [Jonathan Terhorst]
    |\
    | * 4aff2a8 - fixed bug 35, and added a test in test_bugfixes (2 weeks ago) [Hugo]
    |/
    * a7ff2e5 - Added notes on discussion/proposal made during Data Array Summit. (2 weeks ago) [Corran Webster]
    * 68f6752 - Initial implementation of AxisIndexer - uses 'index_by' which needs to be changed to a call on an Axes object - this is all very sketchy right now. (2 weeks ago) [Corr
    *   376adbd - Merge pull request #46 from terhorst/master (2 weeks ago) [Jonathan Terhorst]
    |\
    | * b605216 - updated joshu example to current api (3 weeks ago) [Jonathan Terhorst]
    | * 2e991e8 - add testing for outer ufunc (3 weeks ago) [Jonathan Terhorst]
    | * 7beda5a - prevent axis from throwing an exception if testing equality with non-axis object (3 weeks ago) [Jonathan Terhorst]
    | * 65af65e - convert unit testing code to assertions (3 weeks ago) [Jonathan Terhorst]
    | *   956fbab - Merge remote-tracking branch 'upstream/master' (3 weeks ago) [Jonathan Terhorst]
    | |\
    | |/


.. _section-git-tutorials:

Tutorials and summaries
=======================

There are many, many tutorials and command summaries available online.

Beginner
--------

* `gittutorial <https://git-scm.com/docs/gittutorial>`_ is an introductory tutorial
  from Git project.

* `Git magic
  <http://www-cs-students.stanford.edu/~blynn/gitmagic/index.html>`_
  is an extended introduction with intermediate detail.

* The `Git parable
  <http://tom.preston-werner.com/2009/05/19/the-git-parable.html>`_ is
  an easy read explaining the concepts behind Git.

* Although it also contains more advanced material about branches and
  detached head and the like, the visual summaries of merging and branches
  in `Learn Git Branching <http://pcottle.github.io/learnGitBranching/>`_
  are really quite helpful.


Advanced
--------

* `GitHub help <http://help.github.com>`_ has an excellent series of
  how-to guides.

* The `pro Git book <http://git-scm.com/book>`_ is a good in-depth book on Git.

* `Github Training Kit <http://training.github.com>`_ has an excellent series
  of tutorials as well as videos and screencasts.

* `Git ready <http://www.gitready.com/>`_ is a nice series of
  tutorials.

* A good but technical page on `Git concepts
  <http://www.eecs.harvard.edu/~cduan/technical/git/>`_


Git best practices
==================

There are many ways of working with Git. Here are some posts on the
rules of thumb that other projects have come up with:

* Linus Torvalds on `Git management
  <https://web.archive.org/web/20120511084711/http://kerneltrap.org/Linux/Git_Management>`_.

* Linus Torvalds on `Git workflow
  <http://www.mail-archive.com/dri-devel@lists.sourceforge.net/msg39091.html>`_. Summary:
  use the Git tools to make the history of your edits as clean as
  possible; merge from upstream edits as little as possible in
  branches where you are doing active development.


Manual pages online
===================

You can get these on your own machine with (e.g) ``git help push`` or
(same thing) ``git push --help``, but, for convenience, here are the
online manual pages for some common commands:

* `git add <https://git-scm.com/docs/git-add>`_
* `git branch <https://git-scm.com/docs/git-branch.html>`_
* `git checkout <https://git-scm.com/docs/git-checkout.html>`_
* `git clone <https://git-scm.com/docs/git-clone.html>`_
* `git commit <https://git-scm.com/docs/git-commit.html>`_
* `git config <https://git-scm.com/docs/git-config.html>`_
* `git diff <https://git-scm.com/docs/git-diff.html>`_
* `git log <https://git-scm.com/docs/git-log.html>`_
* `git pull <https://git-scm.com/docs/git-pull.html>`_
* `git push <https://git-scm.com/docs/git-push.html>`_
* `git remote <https://git-scm.com/docs/git-remote.html>`_
* `git status <https://git-scm.com/docs/git-status.html>`_

