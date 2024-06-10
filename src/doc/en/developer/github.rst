.. highlight:: shell-session


.. _chapter-github:

=============================
The Sage Repository on GitHub
=============================

The center of Sage development is `the SageMath organization on GitHub
<https://github.com/sagemath>`_, which consists of many repositories related
with Sage. The most important one among them is of course `the Sage repository
<https://github.com/sagemath/sage>`_, which we call "the Sage repo" for short.


.. _section-github-account:

Obtaining a GitHub account
==========================

To share your work on Sage, you need a GitHub account. If you do not have one
yet, choose a username and `create an account <https://github.com/join>`_. In
the following, we assume your username "alice". So you always read your own
username if you see "alice".


.. _section-github-cli:

Using the GitHub CLI
====================

GitHub provides a command-line interface, the GitHub CLI, that can be used
instead of the web interface.  The central component of the GitHub CLI is the
``gh`` command that you can use in your terminal.

Installation
------------

The page :ref:`spkg_github_cli` documents how to install the ``gh`` command for
your platform. Or see `GitHub CLI <https://cli.github.com>`_ from GitHub.

Configuration
-------------

You have to authenticate to your GitHub account to allow ``gh`` command to
interact with GitHub. Typically the authorization proceeds as follows::

    [alice@localhost sage]$ gh auth login
    ? What is your preferred protocol for Git operations? HTTPS
    ? Authenticate Git with your GitHub credentials? Yes
    ? How would you like to authenticate GitHub CLI? Login with a web browser

    ! First copy your one-time code: 3DA8-5ADA
    Press Enter to open github.com in your browser...
    ✓ Authentication complete.
    - gh config set -h github.com git_protocol https
    ✓ Configured git protocol
    ✓ Logged in as sage

where a web browser is used to enter credentials. You can also use an
authentication token instead, in which case you must first generate `a Personal
Access Token here <https://github.com/settings/tokens>`_.

Next set the default repo for the ``gh`` command::

    [alice@localhost sage]$ gh repo set-default sagemath/sage

and check::

    [alice@localhost sage]$ gh repo view
    sagemath/sage
    ...

which will show the default repo along with its readme, which is quite long.

``gh`` extensions
-----------------

``gh`` is extendable; e.g.  a useful extension to ``gh`` allows testing of
Sage's GitHub Actions locally, using Docker. It is called ``act`` and can be
installed by running::

    [alice@localhost sage]$ gh extension install https://github.com/nektos/gh-act

Append ``--force`` flag to the command above to force an upgrade of the extension.
More details on configuring and using ``gh act`` are in :ref:`chapter-portability_testing`.


Linking Git to your GitHub account
==================================

In order for your Git to work with GitHub, your GitHub account needs to be
linked with your Git. No action is needed if you have already contributed to
any other project on GitHub and set up Git credentials or SSH keys for this.

The above dialogue from ``gh auth login`` linked your Git to GitHub by HTTPS
protocol using your GitHub credentials. Alternatively you may like to use
popular `Git Credential Manager
<https://github.com/git-ecosystem/git-credential-manager>`_ which stores your
credentials natively to your platform. For more information, see `Caching your
GitHub credentials in Git
<https://docs.github.com/en/get-started/getting-started-with-git/caching-your-github-credentials-in-git>`_.

If you prefer SSH to HTTPS for authenticating Git to GitHub, then you follow
:ref:`section-git-ssh` to generate an SSH keypair and add the SSH public key to
your GitHub account. A simple way to upload the public key is to choose SSH
protocol for Git operations in the dialogue from ``gh auth login`` command
above. For more details, see `Connecting to GitHub with SSH
<https://docs.github.com/en/authentication/connecting-to-github-with-ssh>`_.

We assume HTTPS protocol in the rest of this guide.


Forking the Sage repository
===========================

The first step is to create `your personal fork
<https://docs.github.com/en/get-started/quickstart/fork-a-repo#forking-a-repository>`_
of the Sage repo on GitHub. After logging in to your GitHub account, visit the
Sage repo https://github.com/sagemath/sage, and simply click "Fork" on the Sage
repo. Then your fork of the Sage repo is created at
https://github.com/alice/sage.

Next if you don't have a local Git repo of Sage, then start afresh `cloning
your fork
<https://docs.github.com/en/get-started/quickstart/fork-a-repo#cloning-your-forked-repository>`_:

.. tab:: By HTTPS protocol

   ::

    [alice@localhost ~]$ git clone https://github.com/alice/sage.git
    Cloning into 'sage'...
    remote: Enumerating objects: 914565, done.
    remote: Counting objects: 100% (2738/2738), done.
    remote: Compressing objects: 100% (855/855), done.
    remote: Total 914565 (delta 1950), reused 2493 (delta 1875), pack-reused 911827
    Receiving objects: 100% (914565/914565), 331.09 MiB | 11.22 MiB/s, done.
    Resolving deltas: 100% (725438/725438), done.
    Updating files: 100% (9936/9936), done.
    [alice@localhost ~]$ cd sage
    [alice@localhost sage]$ git remote -v
    origin  https://github.com/alice/sage.git (fetch)
    origin  https://github.com/alice/sage.git (push)

.. tab:: By SSH protocol

   ::

    [alice@localhost ~]$ git clone git@github.com:alice/sage.git
    Cloning into 'sage'...
    remote: Enumerating objects: 914565, done.
    remote: Counting objects: 100% (2738/2738), done.
    remote: Compressing objects: 100% (855/855), done.
    remote: Total 914565 (delta 1950), reused 2493 (delta 1875), pack-reused 911827
    Receiving objects: 100% (914565/914565), 331.09 MiB | 11.22 MiB/s, done.
    Resolving deltas: 100% (725438/725438), done.
    Updating files: 100% (9936/9936), done.
    [alice@localhost ~]$ cd sage
    [alice@localhost sage]$ git remote -v
    origin  git@github.com:alice/sage.git (fetch)
    origin  git@github.com:alice/sage.git (push)


If you already have a local Git repo and only want to link your fork as ``origin`` remote, then do:

.. tab:: By HTTPS protocol

   ::

    [alice@localhost sage]$ git remote add origin https://github.com/alice/sage.git
    [alice@localhost sage]$ git remote -v
    origin  https://github.com/alice/sage.git (fetch)
    origin  https://github.com/alice/sage.git (push)
    [alice@localhost sage]$ git fetch origin
    remote: Enumerating objects: 1136, done.
    remote: Counting objects: 100% (1084/1084), done.
    remote: Compressing objects: 100% (308/308), done.
    remote: Total 1136 (delta 825), reused 982 (delta 776), pack-reused 52
    Receiving objects: 100% (1136/1136), 2.62 MiB | 5.30 MiB/s, done.
    Resolving deltas: 100% (838/838), completed with 145 local objects.
    From https://github.com/alice/sage
     * [new branch]      develop     -> origin/develop

.. tab:: By SSH protocol

   ::

    [alice@localhost sage]$ git remote add origin git@github.com:alice/sage.git
    [alice@localhost sage]$ git remote -v
    origin  git@github.com:alice/sage.git (fetch)
    origin  git@github.com:alice/sage.git (push)
    [alice@localhost sage]$ git fetch origin
    remote: Enumerating objects: 1136, done.
    remote: Counting objects: 100% (1084/1084), done.
    remote: Compressing objects: 100% (308/308), done.
    remote: Total 1136 (delta 825), reused 982 (delta 776), pack-reused 52
    Receiving objects: 100% (1136/1136), 2.62 MiB | 5.30 MiB/s, done.
    Resolving deltas: 100% (838/838), completed with 145 local objects.
    From git@github.com:alice/sage
     * [new branch]      develop     -> origin/develop

You also add the Sage repo ``sagemath/sage`` as your remote ``upstream``:

.. tab:: By HTTPS protocol

   ::

    [alice@localhost sage]$ git remote add upstream https://github.com/sagemath/sage.git
    [alice@localhost sage]$ git remote -v
    origin  https://github.com/alice/sage.git (fetch)
    origin  https://github.com/alice/sage.git (push)
    upstream    https://github.com/sagemath/sage.git (fetch)
    upstream    https://github.com/sagemath/sage.git (push)

.. tab:: By SSH protocol

   ::

    [alice@localhost sage]$ git remote add upstream git@github.com:sagemath/sage.git
    [alice@localhost sage]$ git remote -v
    origin  git@github.com:alice/sage.git (fetch)
    origin  git@github.com:alice/sage.git (push)
    upstream    git@github.com:sagemath/sage.git (fetch)
    upstream    git@github.com:sagemath/sage.git (push)

To prevent accidental pushes to ``upstream`` (instead of ``origin``), you may want to disable it by running::

    [alice@localhost sage]$ git remote set-url --push upstream DISABLE

Of course, you can give arbitrary names to your Git remotes, but ``origin`` and
``upstream`` are the established defaults, which will make it easier to use tools
such as the GitHub CLI.


.. _section-github-bug-report:

Reporting bugs
==============

If you think you have found a bug in Sage, here is the procedure:

- Search through our Google groups `sage-devel <https://groups.google.com/group/sage-devel>`_, `sage-support <https://groups.google.com/group/sage-support>`_ for postings related to your possible bug (it
  may have been fixed/reported already). You also search `the GitHub issues
  <https://github.com/sagemath/sage/issues>`_ to see if anyone else has already
  opened an issue about your bug.

- If you do not find anything but you are not sure that you have found a bug,
  ask about it on `sage-devel <https://groups.google.com/group/sage-devel>`_.

- If you are sure that you have found a bug, then create on GitHub a new issue about the bug.

  A bug report should contain:

  - An explicit and **reproducible example** illustrating your bug (and/or the
    steps required to reproduce the buggy behavior). It also helps to describe what
    behaviour is expected.

  - The **version of Sage** you run, as well as the version of the optional
    packages that may be involved in the bug.

  - If relevant, describe your **operating system** as accurately as you can and the
    architecture of your CPU (32 bit, 64 bit, ...).

Follow :ref:`section-github-create-issue` for further guide. Thank you in
advance for reporting bugs to improve Sage!


.. _section-github-new-enhancement:

Planning an enhancement
=======================

In addition to bug reports, you should also open an issue if you have some new
code or an idea that makes Sage better. If you have a feature request, start a
discussion on `sage-devel <https://groups.google.com/group/sage-devel>`_ first,
and then if there seems to be a general agreement that you have a good idea,
open an issue describing the idea.

Before opening a new issue, consider the following points:

- Make sure that nobody else has opened an issue (or a PR) about the same
  or closely related issue. Search through the existing issues and PRs with
  some key words.

- It is much better to open several specific issues than one that
  is very broad. Indeed, a single issue which deals with lots of
  different issues can be quite problematic, and should be avoided.

- Be precise: If foo does not work on macOS but is fine on Linux,
  mention that in the title. Use the keyword option so that
  searches will pick up the issue.

- The problem described in the issue must be solvable. For
  example, it would be silly to open an issue whose purpose was
  "Make Sage the best mathematical software in the world". There is
  no metric to measure this properly and it is highly subjective.

- If appropriate, provide URLs to background information or sage-devel
  conversation relevant to the issue you are reporting.


.. _section-github-create-issue:

Opening an issue
================

Whether it's reporting a bug or planning an enhancement, `issue
<https://docs.github.com/en/issues/tracking-your-work-with-issues/about-issues>`_
should be opened on our Sage repo `sagemath/sage
<https://github.com/sagemath/sage/issues>`_ on GitHub.

- Think of an apt title. People scan through the titles of issues to decide
  which ones to look into further. So write a title that concisely describes
  what the issue is about.

- Describe the issue in detail in the issue body. What is the issue? How can we
  solve the issue? Add links to relevant issues/PRs, and other resources.

  You may use GitHub mention ``@USERNAME`` to get attention from the people
  who would be interested in the issue or has expertise in this issue.

- Add appropriate labels to the created issue:

  - **Type** labels with prefix ``t:`` such as ``t: bug``, ``t: enhancement``,
    ``t: feature``, ``t: performance``, ``t: refactoring``,
    ``t: tests``

  - **Component** labels with prefix ``c:`` such as ``c: basic arithmetic``,
    ``c: linear algebra``, ``c: geometry``, etc.

  - **Priority** labels with prefix ``p:`` such as ``p: trivial / 5``,
    ``p: minor / 4``, ``p: major / 3``, ``p: critical / 2``, and ``p: blocker / 1``

  If the issue is not expected to be solved in the near future, you may add
  ``wishlist item`` label.


.. _section-github-create-pr:

Creating a Pull Request
=======================

If you worked on an issue, and prepared a fix for a bug or wrote code for
enhancing Sage, then you create a PR on the Sage repo `sagemath/sage
<https://github.com/sagemath/sage/issues>`_.

In addition to what were said about opening an issue, the following applies:

- The title should concisely describe what the PR does. If the PR solves an
  issue, describe briefly what the PR solves (do not simply put the issue
  number in the title).

- Explain what the PR solves in detail in the body. If the PR solves an issue,
  you may mention the issue here.

- Add type, component, and priority labels. If this PR solves an existing
  issue, please duplicate the labels of the issue to this PR.

- **Dependencies**: Use the phrase ``- Depends on``, followed by the issue or PR
  reference. Repeat this in separate lines if there is more than one
  dependency. This format is understood by various dependency managers.

If you are working on a PR and the PR is not yet quite ready for review, then
`open the PR as draft
<https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-pull-requests#draft-pull-requests>`_.


.. _section-github-pr-status:

The status of a PR
==================

If a PR is in the state of draft, the review process does not start. Otherwise,
review process will start for the PR as soon as a reviewer gets interested with
the PR, and the status of the PR will be indicated by **status** labels with
prefix ``s:``.

- ``s: needs review``: The code is ready to be peer-reviewed. If the code is not
  yours, then you can review it. See :ref:`chapter-review`.

- ``s: needs work``: Something needs to be changed in the code. The reason should
  appear in the comments.

- ``s: needs info``: The author of the PR or someone else should answer to a
  question or provide information to proceed the review process.

- ``s: positive review``: The PR has been reviewed positively, and the release manager
  will merge it to the ``develop`` branch of the Sage repo in due time.

If the PR does not get positive review and it is decided to close the PR, then
the PR will get one of **resolution** labels: ``r: duplicate``, ``r: invalid``,
``r: wontfix``, ``r: worksforme``.


.. _section-github-stopgaps:

The stopgap
===========

When Sage returns wrong results, an issue and a PR should be created:

- A stopgap issue with all available details.
- A stopgap PR (e.g. :issue:`12699`)

The stopgap PR does not fix the problem but adds a warning that will be
printed whenever anyone uses the relevant code, until the problem is
finally fixed.

To produce the warning message, use code like the following:

.. CODE-BLOCK:: python

    from sage.misc.stopgap import stopgap
    stopgap("This code contains bugs and may be mathematically unreliable.",
        ISSUE_NUM)

Replace ``ISSUE_NUM`` by the reference number for the stopgap issue. On the stopgap issue, enter the reference number for the stopgap PR. Stopgap issues and PRs should be marked as critical.

.. NOTE::

    If mathematically valid code causes Sage to raise an error or
    crash, for example, there is no need for a stopgap.  Rather,
    stopgaps are to warn users that they may be using buggy code; if
    Sage crashes, this is not an issue.


Commenting issues and PRs
=========================

Anyone can comment on an issue or a PR. If a PR is linked to an issue,
you may not be sure where the comment should go. Then

- Comments on the reported issue should go on the issue.

- Comments on the submitted code should go on the PR.


Checks on PRs
=============

If you manage to fix a bug or enhance Sage, you are our hero. See
:ref:`chapter-walkthrough` for making changes to the Sage source code and
:ref:`section-github-create-pr` to create a PR for the changes.

For each push to a PR, automated tests for the branch of the PR run on GitHub
Actions.

* A `linting workflow
  <https://github.com/sagemath/sage/blob/develop/.github/workflows/lint.yml>`_
  checks that the code of the current branch adheres to the style guidelines
  using :ref:`section-tools-pycodestyle` (in the ``pycodestyle-minimal``
  configuration) and :ref:`section-tools-relint`.

  In order to see details when it fails, you can click on the check
  and then select the most recent workflow run.

* The `build and test workflow
  <https://github.com/sagemath/sage/blob/develop/.github/workflows/build.yml>`_
  on GitHub Actions builds Sage for the current branch (incrementally
  on top of an installation of the ``develop`` branch) and runs the
  test. Details are again available by clicking on the check.

  The automatic workflow runs on a container based on
  ``ubuntu-focal-standard``.  To request a run of the workflow on a different
  platform, you can issue a `workflow dispatch
  <https://docs.github.com/en/actions/managing-workflow-runs/manually-running-a-workflow#running-a-workflow>`_.
  You can select any of the platforms for which a `prebuilt container image
  <https://github.com/orgs/sagemath/packages?tab=packages&q=with-targets-optional>`_
  exists.

* The `build documentation workflow
  <https://github.com/sagemath/sage/blob/develop/.github/workflows/doc-build.yml>`_
  on GitHub Actions builds the HTML documentation for the current branch.

  A link to the built doc is added in a comment, and so you can easily inspect changes
  to the documentation without the need to locally rebuild the docs yourself.

  If the doc build fails, you can go to Actions tab and examine `documentation
  build workflow
  <https://github.com/sagemath/sage/actions/workflows/doc-build.yml>`_ and
  choose the particular branch to see what went wrong.


Final notes
===========

* Every bug fixed should result in a doctest.

* There are many enhancements possible for Sage and too few developers to
  implement all the good ideas.

* If you are a developer, be nice and try to solve a stale/old issue
  every once in a while.

* Some people regularly do triage. In this context, this means that we
  look at new bugs and classify them according to our perceived
  priority. It is very likely that different people will see
  priorities of bugs very differently from us, so please let us know
  if you see a problem with specific PRs.

