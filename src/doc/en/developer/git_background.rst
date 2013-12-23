.. _chapter-git-background:

===================
Tips and References
===================

This chapter contains additional material about the git revision
control system. It is not necessary if you stick with the Sage
development scripts. See :ref:`chapter-git-setup` for the minimal
steps needed for Sage development.






.. _section-git-configuration:

Configuration Tips
==================

Your personal git configurations are saved in the ``~/.gitconfig``
file in your home directory. Here is an example::

    [user]
        name = Your Name
        email = you@yourdomain.example.com

    [core]
        editor = emacs

You can edit this file directly or you can use git to make changes for
you::

    [user@localhost ~] git config --global user.name "Your Name"
    [user@localhost ~] git config --global user.email you@yourdomain.example.com
    [user@localhost ~] git config --global core.editor vim



Aliases
-------

Aliases are personal shortcuts for git commands. For example, you
might want to be able to shorten ``git checkout`` to ``git co``.  Or
you may want to alias ``git diff --color-words`` (which gives a nicely
formatted output of the diff) to ``git wdiff``. You can do this with::

    [user@localhost ~] git config --global alias.ci "commit -a"
    [user@localhost ~] git config --global alias.co checkout
    [user@localhost ~] git config --global alias.st "status -a"
    [user@localhost ~] git config --global alias.stat "status -a"
    [user@localhost ~] git config --global alias.br branch
    [user@localhost ~] git config --global alias.wdiff "diff --color-words"

The above commands will create an ``alias`` section in your ``.gitconfig``
file with contents like this::

    [alias]
        ci = commit -a
        co = checkout
        st = status -a
        stat = status -a
        br = branch
        wdiff = diff --color-words


Editor
------

To set the editor to use for editing commit messages, you can use::

    [user@localhost ~] git config --global core.editor vim

or set the `EDITOR` environment variable.

Merging
-------

To enforce summaries when doing merges (``~/.gitconfig`` file again)::

    [merge]
        log = true

Or from the command line::

    [user@localhost ~] git config --global merge.log true


.. _section-fancy-log:

Fancy Log Output
----------------

Here is an alias to get a fancy log output; it should go in the
``alias`` section of your ``.gitconfig`` file::

    lg = log --graph --pretty=format:'%Cred%h%Creset -%C(yellow)%d%Creset %s %Cgreen(%cr) %C(bold blue)[%an]%Creset' --abbrev-commit --date=relative

Using this ``lg`` alias gives you the changelog with a colored ascii graph::

    [user@localhost ~] git lg
    * 6d8e1ee - (HEAD, origin/my-fancy-feature, my-fancy-feature) NF - a fancy file (45 minutes ago) [Matthew Brett]
    *   d304a73 - (origin/placeholder, placeholder) Merge pull request #48 from hhuuggoo/master (2 weeks ago) [Jonathan Terhorst]
    |\
    | * 4aff2a8 - fixed bug 35, and added a test in test_bugfixes (2 weeks ago) [Hugo]
    |/
    * a7ff2e5 - Added notes on discussion/proposal made during Data Array Summit. (2 weeks ago) [Corran Webster]
    * 68f6752 - Initial implimentation of AxisIndexer - uses 'index_by' which needs to be changed to a call on an Axes object - this is all very sketchy right now. (2 weeks ago) [Corr
    *   376adbd - Merge pull request #46 from terhorst/master (2 weeks ago) [Jonathan Terhorst]
    |\
    | * b605216 - updated joshu example to current api (3 weeks ago) [Jonathan Terhorst]
    | * 2e991e8 - add testing for outer ufunc (3 weeks ago) [Jonathan Terhorst]
    | * 7beda5a - prevent axis from throwing an exception if testing equality with non-axis object (3 weeks ago) [Jonathan Terhorst]
    | * 65af65e - convert unit testing code to assertions (3 weeks ago) [Jonathan Terhorst]
    | *   956fbab - Merge remote-tracking branch 'upstream/master' (3 weeks ago) [Jonathan Terhorst]
    | |\
    | |/



Tutorials and Summaries
=======================

* `Github help <http://help.github.com>`_ has an excellent series of
  how-to guides.

* The `pro git book <http://git-scm.com/book>`_ is a good in-depth book on git.

* `Github Training <http://training.github.com>`_ has an excellent series
  of tutorials as well as videos and screencasts.

* A `git cheat sheet <http://github.com/guides/git-cheat-sheet>`_ is a
  page giving summaries of common commands.

* The `git user manual
  <http://schacon.github.com/git/user-manual.html>`_.

* The `git tutorial <http://schacon.github.com/git/gittutorial.html>`_.

* `Git ready <http://www.gitready.com/>`_ is a nice series of
  tutorials.

* `Git magic
  <http://www-cs-students.stanford.edu/~blynn/gitmagic/index.html>`_
  is an extended introduction with intermediate detail.

* The `git parable
  <http://tom.preston-werner.com/2009/05/19/the-git-parable.html>`_ is
  an easy read explaining the concepts behind git.

* `Git foundation
  <http://matthew-brett.github.com/pydagogue/foundation.html>`_
  expands on the `git parable`_.

* `Fernando Perez' git page
  <http://www.fperez.org/py4science/git.html>`_ contains many links
  and tips.

* A good but technical page on `git concepts
  <http://www.eecs.harvard.edu/~cduan/technical/git/>`_

* `Git svn crash course <http://git-scm.com/course/svn.html>`_: git
  for those of us used to `subversion
  <http://subversion.tigris.org/>`_


Git Best Practices
==================

There are many ways of working with git; here are some posts on the
rules of thumb that other projects have come up with:

* Linus Torvalds on `git management
  <https://web.archive.org/web/20120511084711/http://kerneltrap.org/Linux/Git_Management>`_

* Linus Torvalds on `linux git workflow
  <http://www.mail-archive.com/dri-devel@lists.sourceforge.net/msg39091.html>`_. Summary:
  use the git tools to make the history of your edits as clean as
  possible; merge from upstream edits as little as possible in
  branches where you are doing active development.


Manual Pages Online
===================

You can get these on your own machine with (e.g) ``git help push`` or
(same thing) ``git push --help``, but, for convenience, here are the
online manual pages for some common commands:

* `git add <http://schacon.github.com/git/git-add.html>`_
* `git branch <http://schacon.github.com/git/git-branch.html>`_
* `git checkout <http://schacon.github.com/git/git-checkout.html>`_
* `git clone <http://schacon.github.com/git/git-clone.html>`_
* `git commit <http://schacon.github.com/git/git-commit.html>`_
* `git config <http://schacon.github.com/git/git-config.html>`_
* `git diff <http://schacon.github.com/git/git-diff.html>`_
* `git log <http://schacon.github.com/git/git-log.html>`_
* `git pull <http://schacon.github.com/git/git-pull.html>`_
* `git push <http://schacon.github.com/git/git-push.html>`_
* `git remote <http://schacon.github.com/git/git-remote.html>`_
* `git status <http://schacon.github.com/git/git-status.html>`_



