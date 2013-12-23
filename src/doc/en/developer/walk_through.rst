.. _chapter-walk-through:

========================
Sage Development Process
========================

This section is a concise overview of the Sage development process. In
it, we will see how to make changes to the Sage source code and
communicate these changes back to the Sage project. If you are a
beginner to Sage development, this introductory guide is here to help
you become familiar with the Sage development process.

Sage comes with a set of developer scripts, which help you with common
interactions with the bug tracker (see :ref:`chapter-sage-trac`) and
with handling revisions of your code. The developer scripts use the
git distributed revison control system under the hood which you'll
have to install (see :ref:`chapter-git-setup`), but you do not need to
know anything about it (see :ref:`chapter-manual-git` only if you want
to).

Most of the commands in the following section will not work unless you have an
account on the bug tracker. If you want to contribute to Sage, it is a good
idea to get an account now (see :ref:`section-trac-account`).

We assume here that the ``sage`` executable of your development
installation of Sage is in your ``PATH``. If this is not the case, you
might have to replace ``sage`` by ``./sage`` or ``/path/to/sage`` in
the following. You can also use the developer scripts from the Sage
prompt. All commandline options to ``sage -dev`` are also available as
methods of the ``dev`` object in a Sage session. That is, for example,
to checkout a ticktet you can either run::

    [user@localhost]$ sage -dev checkout --ticket 1729
    On ticket #1729 with associated local branch "ticket/1729".

    #  Use "sage --dev merge" to include another ticket/branch.
    #  Use "sage --dev commit" to save changes into a new commit.

in a terminal or, equivalently, within Sage

.. skip   # don't actually doctest

::

    sage: dev.checkout(1729)
    On ticket #1729 with associated local branch "ticket/1729".
 
    #  Use "dev.merge()" to include another ticket/branch.
    #  Use "dev.commit()" to save changes in a new commit.

Note that the number sign ``#`` (a.k.a. hash or pound sign) is the
comment marker for both the shell and Python. So if you were to input
``#1729``, it will be interpreted as the comment "1729" and not passed
to the development scripts. Always specify the ticket number as a
plain number, without the number sign in front.

.. warning::

    During the transitional period it can happen that you end up
    on a branch where the developer scripts are not available or
    outdated. If this is the case, i.e., if ``sage -dev`` does not
    work properly anymore, run::

        git pull git://trac.sagemath.org/sage.git master
        sage -b

    This will merge the latest version of the developer scripts
    with your current branch. After rebuilding the Sage library,
    the dev scripts will work again.


.. _section-walkthrough-add:

Contributing to the Sage Source Code
====================================

.. _section-walkthrough-add-create:

Create a Ticket
---------------

Suppose you have written an algorithm for calculating the last twin prime, and
want to add it to Sage. You would first open a ticket for that::

    [user@localhost]$ sage -dev create-ticket

This will give you an editor in which you can give a summary and a
description of what you want to do. If you are not sure which values
to put for the other fields, you can leave the defaults or have a look
at :ref:`section-trac-fields`. After you close the editor, a new
ticket will be opened on the trac server. From that point on, everyone
can see what you intend to do which lets us avoid duplicating work. If
you want to cancel the creation of a ticket, then you can just save an
empty file. This will abort the operation.

Alternatively, you can use the `web interface to the Sage trac
development server <http://trac.sagemath.org>`_ to open a new ticket,
just log in and click on "Create Ticket".


.. _section-walkthrough-add-edit:

Editing the Source Code
-----------------------

If you want to work on a ticket which you or somebody else created,
you first need to make a local "branch". The development scripts
maintain a mapping between local branches and trac tickets. Creating a
new local branch for a ticket is easy::

    [user@localhost]$ sage -dev checkout --ticket 1729
    On ticket #1729 with associated local branch "ticket/1729".

    #  Use "sage --dev merge" to include another ticket/branch.
    #  Use "sage --dev commit" to save changes into a new commit.

Essentially, a branch is a copy (except that it doesn't take up twice
the space) of the Sage source code where you can store your
modifications to the Sage source code and which you can upload to trac
tickets. Your new branch is now called ``ticket/<TICKETNUM>``. Unless
you upload ("push") it, see below, it will only be on your local
system and not visible to anyone else.

At this point you can start editing the source code. The subsequent
chapters of this developer guide explain how your code should look
like to fit into Sage, and how we ensure high code quality
throughout. Whenever you have reached one of your goals, you should
make a *commit*. This takes a snapshot of the whole Sage source code
that you have been working on and records the changes into your local
branch::

    [user@localhost]$ sage -dev commit
    Commit your changes to branch "ticket/1729"? [Yes/no] y

    #  Use "sage --dev push" to push your commits to the trac server once you are
    #  done.

You will be asked to write a message describing your changes. It is
common to write a one line summary, then a blank line, and then a 1-2
paragraph explanation of your changes. If your changes are minor, then
just the one-line summary can be enough.

If you are working on a larger project, it can be useful to break up
your work into multiple commits: Each commit is saved, enabling you to
retrieve older versions of files from the repository. So, even if you
accidentally delete something, you can get it back later. Also, if you
find a mistake in one of your earlier commits, then you just correct
it in the Sage source code and then add another commit on top.


.. _section-walkthrough-add-push:

Uploading Changes to Trac
-------------------------

At some point, you may wish to share your changes with the rest of us:
maybe it is ready for review, or maybe you are collaborating with
someone and want to share your changes "up until now". This is simply
done by::

    [user@localhost]$ sage -dev push

On trac, your remote branch will be called
``u/<USERNAME>/ticket/<TICKETNUM>``. This name will automatically be
added to the "Branch:" field on the ticket. Other developers then know
where to find your work in the git repository.

It is common to go through some iterations of ``sage -dev commit``
before you upload, and you will probably also have uploaded a few
times before your changes are ready for review.

If you are happy with the changes you uploaded, you want somebody else
to review them, so they can be included into the next version of
Sage. If your ticket is ready for review, you should set it to
``needs_review`` on the trac server. This can be done though the `web
interface <http://trac.sagemath.org>`_, or, alternatively, using the
development scripts. For the latter, run::

    [user@localhost]$ sage -dev edit-ticket

This will give you an editor in which you can edit the ticket. Change the
status to::

    Status: needs_review

And add yourself as an author for that ticket by inserting the following as the
first line::

    Authors: Your Real Name

If you want to add an additional comment for potential reviewers, run::

    [user@localhost]$ sage -dev comment


.. _section-walkthrough-add-local:

Starting Without a Ticket
-------------------------

You might not want to create a trac ticket for your changes. For
example, if you are only working on your own code or if you are making
experimental changes that you are likely to throw away if they do not
work out. In that case, you can also start a branch that only lives in
your local repository. To do this, you use checkout but specify a
branch name instead of the ticket number. For example, to create a new
branch ``my_branch``, you would run::

    [user@localhost]$ sage -dev checkout --branch my_branch

This is assuming that you do not already have a local branch called
``my_branch``. If that were the case, you would just switch to the
already-existing branch. Once on your branch, you can work with it as
described in :ref:`section-walkthrough-add-edit`.

You can upload your local branch later to an existing ticket. This
works exactly like in the case where you started with a ticket, except
that you have to specify the ticket number. That is::

    [user@localhost]$ sage -dev push --ticket <TICKETNUM>
    
where you have to replace ``<TICKETNUM>`` with the number of the trac
ticket. 


.. _section-walkthrough-merge:

Merging
=======

As soon as you are working on a bigger project that spans multiple
tickets you will want to base your work on branches that have not been
merged into Sage yet. This is natural in collaborative development,
and in fact you are very much encouraged to split your work into
logically different parts. Ideally, each part that is useful on its
own and and can be reviewed independently should be a different
ticket, instead of a huge patch bomb.

For this purpose, you can incorporate branches from other tickets (or
just other local branches) into your current branch. This is called
merging, and all it does is include commits from other branches into
your current branch. In particular, this is done when a new Sage
release is made: the finished tickets are merged with the Sage master
and the result is the next Sage version. Git is smart enough to not
merge commits twice. In particular, it is possible to merge two
branches, one of which had already merged the other branch.

The syntax for merging is easy. If the code that you want to
incorporate is on a trac ticket number ``<TICKETNUM>``, use::

    [user@localhost]$ sage -dev merge --ticket <TICKETNUM>

Optionally, you can add the merged ticket to the trac "Dependency:"
field. Note that the merged commits become part of the current branch,
regardless of whether they are noted on trac. Adding a dependency
implies that the dependency must be reviewed first. After the
dependency is reviewed, the commits that came from the dependency are
no longer listed in the output of ``sage -dev diff``.

.. warning::

    You should avoid merging tickets both ways. Once ticket A merged
    ticket B and ticket B merged ticket A, there is no way to
    distinguish commits that were originally made in ticket A or in
    ticket B. Effectively, merging both ways combines the branches and
    makes individual review impossible.

    In practice, you should only merge when one of the following holds:

    * Either two tickets conflict, then you have to merge one into the
      other in order to resolve the merge conflict.

    * Or you definitely need a feature that has been developed as part
      of another branch.

A special case of merging is merging in the ``master`` branch. This
brings your local branch up to date with the newest Sage version. The
above warning against unnecessary merges still applies, though. Try to
do all of your development with the Sage version that you originally
started with. The only reason for merging in the master branch is if
you need a new feature or if your branch conflicts.


.. _section-walkthrough-review:

Reviewing
=========

Now suppose you want to review the existing work on a ticket, such as the one
you created in the last section.  For definiteness, suppose you want to review
#12270. You would do that as follows::

    [user@localhost]$ sage -dev checkout --ticket 12270

This command will download the branch on Trac in case you do not have any local
work on ticket 12270. (If you do, you may have to merge your changes; see
below). You can now test the ticket; you'll probably want to call ``make`` or
``sage -b`` first to rebuild Sage with the changes. Another important
command is::

    [user@localhost]$ sage -dev diff

which lists all source code changes that are part of the current
branch. That is, it lists the changes from the current master to the
current branch. If the ticket were to be positively reviewed, this is
the code that will be added to Sage. Note that there is no way to
"exclude dependencies", just as there is no guarantee that unreviewed
dependencies will become part of Sage. The best way to exclude
dependencies from the diff output is to review them. Once the
dependency becomes part of the master branch, they are automatically
removed.

Most likely, your will want to add a comment to the ticket as part of
your review::

    [user@localhost]$ sage -dev comment

This will open a text editor in which you can type, and upload the
result to Trac.
    
It is also possible that you make some changes to the code as part of
your review. After you have done that, you can upload your changes
back to trac::

    [user@localhost]$ sage -dev commit
    [user@localhost]$ sage -dev push

This will update the ticket to now point to your branch, including
your changes. Your branch is based on the original author's branch, so
s/he can easily incorporate your changes into his/her own branch (see
below).


.. _section-walkthrough-collaborate:

Collaboration
=============

It is very easy to collaborate by just going through the above steps any number of times::

    # developer 1
    <EDIT EDIT>
    sage -dev commit
    sage -dev push

    # developer 2
    sage -dev pull
    <EDIT EDIT>
    sage -dev commit
    sage -dev push

    # developer 1
    sage -dev pull
    <EDIT EDIT>
    sage -dev commit
    sage -dev push
    (etc)

The obvious problem is when you both work on the same ticket simultaneously::

    # developer 1
    <EDIT EDIT>
    sage -dev commit
    sage -dev push

    # developer 2
    <EDIT EDIT>
    sage -dev commit
    sage -dev push
    Changes not compatible with remote branch
    u/<developer1>/ticket/12270; consider 
    downloading first. Are you sure you want to continue?

Developer 2 should probably select ``No``, and do as suggested::

    sage -dev pull

This will try to merge the changes developer 1 made into the ones that
developer 2 made. The latter should check whether all seems okay, and
if so, upload the changes::

    sage -dev push   # works now

It is possible that the changes cannot be automatically merged. In
that case, developer 2 will have to do some manual fixup after
downloading and before uploading::

    <EDIT EDIT FOR FIXUP>
    sage -dev commit
    sage -dev push


