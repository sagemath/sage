.. _chapter-workflows:

=======================
Distributed Development
=======================

Git is a tool to exchange commits (organized into branches) with other
developers. As a distributed revision control system, it does not have
the notion of a central server. The Sage trac server is just one of
many possible remote repositories from your point of view. This lets
you use and experiment with different ways to interact with other
developers. In this chapter, we describe some common ways to develop
for Sage.

For simplicity, let us assume two developers (Alice and Bob) are
collaborating on a ticket. The first step of opening the ticket is
always the same, and could be performed by either Alice or Bob or a
third person.





Simple Workflow
===============

.. image:: static/flowchart.*
    :align: center


1. Alice creates a :ref:`new local branch <section-git-branch>` and
   :ref:`commits <section-git-commit>` changes to the Sage sources.

2. Alice :ref:`uploads her branch <section-git-push>` to the trac
   server and fills in the "Branch:" field with her remote branch name
   ``u/alice/description``.

3. Bob :ref:`downloads Alice's branch <section-git-checkout>`, looks
   through the source, and leaves a comment on the ticket about a
   mistake in Alice's code.

4. Alice fixes the bug on top of her current branch, and uploads the
   updated branch.

5. Bob :ref:`retrieves Alice's updates <section-git-pull>` and reviews
   the changes.

6. Once Bob is satisfied, he sets the ticket to positive review. The
   "Author:" field is set to Alice's full name, and the "Reviewer:"
   field is set to Bob's full name.

Alternatively, Bob might want to make some changes himself. Then,
instead, we would have

3. Bob :ref:`downloads Alice's branch <section-git-checkout>`, makes
   changes, and :ref:`commits <section-git-commit>` them to his local
   branch.

4. Bob :ref:`uploads his branch <section-git-push>` to the trac server
   and fills in the "Branch:" field with his remote branch name
   ``u/bob/description``.

5. Alice :ref:`downloads Bob's branch <section-git-checkout>` and
   reviews his changes.

6. Once Alice is satisfied, she sets the ticket to positive review. If
   both contributions are of comparable size, then the "Author:" and
   "Reviewer:" fields are set to both Alice's and Bob's full name.




Public Repository
=================

In addition to the user branches (``u/<user>/<description>`` on the
Sage trac server with ``<user>`` replaced by your trac user name) that
only you can write to, you can also create a public branch that
everybody with a trac account can write to. These start with
``public/`` plus some description. To avoid branch name collisions it
is a good idea to include your trac user name in the branch name, so
it is recommended that you use ``public/<user>/<description>`` as the
branch name. Now all ticket authors push to the same remote branch.

1. Alice creates a :ref:`new local branch <section-git-branch>` and
   :ref:`commits <section-git-commit>` some changes to the Sage library.

2. Alice :ref:`uploads her branch <section-git-push>` as a public
   branch to the trac server and fills in the "Branch:" field with her
   remote branch name ``public/alice/description``.

3. Bob :ref:`downloads Alice's branch <section-git-checkout>` and
   makes changes to his local copy.

4. Bob :ref:`commits <section-git-commit>` changes to his local branch
   of the Sage sources.

5. Bob uploads his changes to the joint remote repository::

       [bob@localhost sage]$ git push trac local_branch:public/alice/description

6. Alice :ref:`retrieves Bob's updates <section-git-pull>`, makes
   more changes, commits, and pushes them to trac.

7. Charly reviews the final version, and then sets the ticket to
   positive review. The "Author:" field is set to Alice's and Bob's
   full name, and the "Reviewer:" field is set to Charly's full name.




GitHub
======

Yet another possible workflow is to use GitHub (or any other
third-party git repository) to collaboratively edit your new branch,
and only push the result to trac once you and your ticket co-authors
are satisfied.


Fork
----

The first step is to create your own fork of the Sage repository;
simply click "Fork" on the `Sage GitHub repository
<https://github.com/sagemath/sage>`_. Then add it as one of the
remotes to your local Sage repository. In the following, we will use
the label "github" for this remote repository, though you are of
course free to use a different one::

    $ git remote add github git@github.com:github_user_name/sage.git
    $ git remote -v
    github      git@github.com:github_user_name/sage.git (fetch)
    github      git@github.com:github_user_name/sage.git (push)
    trac        git@trac.sagemath.org:sage.git (fetch)
    trac        git@trac.sagemath.org:sage.git (push)
    $ git fetch github
    remote: Counting objects: 107, done.
    remote: Compressing objects: 100% (63/63), done.
    remote: Total 74 (delta 41), reused 40 (delta 10)
    Unpacking objects: 100% (74/74), done.
    From github.com:github_user_name/sage
    * [new branch]      master     -> github/master
    

Develop
-------

You now use the github repository to develop your ticket branch; First
create a new branch::

    $ git checkout -b my_branch --track github/master
    Branch my_branch set up to track remote branch master from github.
    Switched to a new branch 'my_branch'
    $ git push github my_branch
    Total 0 (delta 0), reused 0 (delta 0)
    To git@github.com:github_user_name/sage.git
     * [new branch]      my_branch -> my_branch

Because of the ``--track`` option, the ``git pull`` command will
default to downloading your coauthor's changes from your github
branch. Alternatively, you can create a new branch on your fork's
GitHub webpage.

At this point you can use the GitHub workflow that you prefer. In
particular, your choices are

* Give your coauthors write permissions to your github fork. Every
  author edits/commits to their own local copy and they jointly push
  to your github branch.

* Have every coauthor create their own fork and send you (the lead
  author) pull requests to your GitHub fork.

* Use the GitHub web page editing & commiting feature, that way you
  can make changes without ever using your local machine.


Push to Trac
------------

When you are satisfied with your branch, you push it to the Sage trac
server::

    $ git push trac u/user/description

and then fill in the "Branch" field in the trac ticket description as
explained in :ref:`section-git-push`.

