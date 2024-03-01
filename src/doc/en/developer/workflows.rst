.. highlight:: shell-session

.. _chapter-workflows:

=====================
Using Git with GitHub
=====================

We continue our introduction to Sage development from :ref:`chapter-walkthrough`.
We discuss how to push your local changes to your fork of the GitHub Sage repository
so that your changes can be reviewed for inclusion in Sage.

Before proceeding, check that you have ``origin`` and ``upstream`` remotes right::

    [alice@localhost sage]$ git remote -v
    origin  https://github.com/alice/sage.git (fetch)
    origin  https://github.com/alice/sage.git (push)
    upstream    https://github.com/sagemath/sage.git (fetch)
    upstream    https://github.com/sagemath/sage.git (push)


Development workflow at a glance
================================

.. IMAGE:: static/workflow.*
    :align: center

1. Alice creates a :ref:`new local branch <section-walkthrough-branch>` and
   :ref:`commits <section-walkthrough-commit>` changes to the Sage source files.

2. Alice pushes the local branch to the remote ``origin``, her fork of the Sage
   repo on GitHub, and with it :ref:`creates a PR <section-workflows-push>` to
   the Sage repo. When ready, Alice sets the PR to ``needs review`` status.

3. Bob, a developer acting as reviewer, :ref:`examines the PR
   <section-workflows-pr-checkout>`, looks through the changes, leaves comments
   on the PR, and requests fixes (``needs work``).

4. Alice makes more commits on top of her local branch, and pushes the new
   commits to the remote ``origin``. These new commits are reflected in the PR.

5. Bob looks through the changes in the new commits and reviews the changes.

6. After a few of iterations of commenting and fixing, finally the reviewer Bob
   is satisfied, and then he approves the PR and sets it to ``positive review``
   status.


.. _section-workflows-pr-create:

Creating a new PR
=================

Suppose you have written an algorithm for calculating the last twin prime,
committed the code to a local branch based upon ``develop`` branch. Now you
want to add it to Sage. You would first open a PR for that::

    [alice@localhost sage]$ gh pr create
    ? Where should we push the 'last-twin-prime' branch? user/sage

    Creating pull request for user:last-twin-prime into develop in sagemath/sage

    ? Title Last twin prime
    ? Choose a template PULL_REQUEST_TEMPLATE.md
    ? Body <Received>
    ? What's next? Submit as draft
    https://github.com/sagemath/sage/pull/12345

This will create a new PR titled "Last twin prime" in the Sage repo for the
branch pushed to your fork ``alice/sage`` from the local branch on your
desktop. The title is automatically derived from the last commit title. If you
don't like this, then you can use the ``-t`` switch to specify it explicitly.
See the manual page of the command `gh pr create
<https://cli.github.com/manual/gh_pr_create>`_ for details.

If you did not provide enough details about the PR at the prompts, you may want
to edit the PR further via the web interface.


.. _section-workflows-pr-checkout:

Checking out an existing PR
===========================

If you want to base your work on an existing PR or want to review the code of a PR,
then you would run::

    [alice@localhost sage]$ gh pr checkout 12345
    remote: Enumerating objects: 7, done.
    remote: Counting objects: 100% (7/7), done.
    remote: Compressing objects: 100% (7/7), done.
    remote: Total 7 (delta 0), reused 0 (delta 0), pack-reused 0
    Unpacking objects: 100% (7/7), 25.50 KiB | 2.83 MiB/s, done.
    From https://github.com/sagemath/sage
     * [new ref]               refs/pull/12345/head -> last-twin-prime
    Switched to branch 'last-twin-prime'

The command ``gh pr checkout`` downloads the branch of the PR. Just
like the ``create`` command, you can specify the local branch name explicitly using
the ``-b`` switch if you want.


.. _section-workflows-push:

Uploading more changes to GitHub
================================

Once you have created a PR, edit the appropriate files and commit your changes
to your local branch as described in :ref:`section-walkthrough-add-edit` and
:ref:`section-walkthrough-commit`.

If you are ready to share the changes up to now, upload your new commits to
your fork by::

    [alice@localhost sage]$ git push origin
    Enumerating objects: 13, done.
    Counting objects: 100% (13/13), done.
    Delta compression using up to 12 threads
    Compressing objects: 100% (7/7), done.
    Writing objects: 100% (7/7), 1.98 KiB | 1.98 MiB/s, done.
    Total 7 (delta 6), reused 0 (delta 0), pack-reused 0
    remote: Resolving deltas: 100% (6/6), completed with 6 local objects.
    To https://github.com/alice/sage.git
     + 352d842907...56ffdab967 last-twin-prime -> last-twin-prime

Note that you do not push the branch to the remote ``upstream`` the Sage repo.
Instead the new commits pushed to the remote ``origin`` are shown in the PR at
the Sage repo.


.. _section-workflows-finish:

Finishing it up
===============

It is common to go through a few iterations of commits before you
push the branch, and you will probably also have pushed your branch a few
times before your branch is ready for review.

Once you are happy with the changes you pushed, they must be
reviewed by someone else before they can be included in the next
release of Sage. To mark your PR as ready for review, you should
set it to ``needs review`` status.


.. _section-workflows-merge:

Merging the upstream develop branch
===================================

It commonly happens that ``develop`` branch at the remote ``upstream`` was
updated and you need to merge the upstream changes to your local branch. Then
you do::

    [alice@localhost sage]$ git fetch upstream develop:develop

This fast-forwards your local ``develop`` branch to the upstream
``develop`` branch.

Now you go back to your working branch and merge the updated ``develop`` branch::

    [alice@localhost sage]$ git merge develop
    ....

If there was no upstream change conflicting with the changes you made locally,
this merge operation will finish cleanly. Otherwise, you are in *merge
conflict*. This rarely happens since Git is smart in merging changes. However,
once merge conflict occurs, you have to manually resolve the conflicts. The
conflict resolving procedure is explained in :ref:`section-git-conflict`.

