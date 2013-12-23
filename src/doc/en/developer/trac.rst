.. _chapter-sage-trac:

====================
The Sage Trac Server
====================

All changes to Sage source code have to go through the `Sage trac
development server <http://trac.sagemath.org>`_. The purpose
of the Sage trac server is to

1. Provide a place for discussion on issues and store a permanent
   record.

2. Provide a repository of source code and all proposed changes.

3. Link these two together.

There is also a `wiki <http://trac.sagemath.org/wiki>`_ for more general
organizational web pages, like Sage development workshops.

Thus if you find a bug in Sage, if you have new code to submit, want
to review new code already written but not yet included in Sage, or if
you have corrections for the documentation, you should post on the
trac server. Items on the server are called *tickets*, and anyone may
search or browse the tickets. For a list of recent changes, just visit
the `Sage trac timeline <http://trac.sagemath.org/timeline>`_.


Authentication
==============

There are two avenues to prove to the trac server that you are who you
claim to be. First, to change the ticket web pages you need to log in
to trac using a username/password. Second, there is public key
cryptography used by git when copying new source files to the
repository. This section will show you how to setup both.


.. _section-trac-account:

Obtaining an Account
--------------------

You first need to open an account if you want to *change* anything on
the Sage trac server, even if you just want to comment on a
ticket. Part of the process is to prove that you are a human to keep
spam at a minimum. To get an account read the developer manual (this
document) and then send an email to
``sage-trac-account@googlegroups.com`` that contains all of the
following:

* your full name,
* preferred username,
* contact email,
* and reason for needing a trac account.

Your trac account also grants you access to the sage wiki. Make sure
you understand the review process, and the procedures for opening and
closing tickets before making changes. The remainder of this chapter
contains various guidelines on using the trac server.

Generating and Uploading your SSH Keys
--------------------------------------

The git installation on the development server uses SSH keys to decide if and
where you are allowed to upload code. No SSH key is required to report a bug or
comment on a ticket, but as soon as you want to contribute code yourself you
need to provide trac with the public half of your own personal key. In recent
versions of Sage, you can use Sage to generate an upload an SSH key

.. skip   # do not doctest

::

    sage: dev.upload_ssh_key()
    The trac git server requires your SSH public key to be able to identify you.
    Upload "/home/vbraun/.ssh/id_dsa.pub" to trac? [Yes/no] y
    Trac username: user
    Trac password:
    Your key has been uploaded.

You can also manually generate an SSH key and upload it to trac. This is
described in the following two sections.


Manually Generating your SSH Keys
---------------------------------

If you don't have a private key yet, you can
create it with the ``ssh-keygen`` tool::

    [user@localhost ~]$ ssh-keygen
    Generating public/private rsa key pair.
    Enter file in which to save the key (/home/user/.ssh/id_rsa):
    Enter passphrase (empty for no passphrase):
    Enter same passphrase again:
    Your identification has been saved in /home/user/.ssh/id_rsa.
    Your public key has been saved in /home/user/.ssh/id_rsa.pub.
    The key fingerprint is:
    ce:32:b3:de:38:56:80:c9:11:f0:b3:88:f2:1c:89:0a user@localhost
    The key's randomart image is:
    +--[ RSA 2048]----+
    |  ....           |
    |   ..            |
    |   .o+           |
    | o o+o.          |
    |E + .  .S        |
    |+o .   o.        |
    |. o   +.o        |
    |      oB         |
    |     o+..        |
    +-----------------+

This will generate a new random private RSA key
in the ``.ssh`` folder in your home directory. By default, they are

``~/.ssh/id_rsa``
  Your private key. Keep safe. **Never** hand it out to anybody.

``~/.ssh/id_rsa.pub``
  The corresponding public key. This and only this file can be safely
  disclosed to third parties.

The ``ssh-keygen`` tool will let you generate a key with a different
file name, or protect it with a passphrase. Depending on how much you
trust your own computer or system administrator, you can leave the
passphrase empty to be able to login without any human intervention.

If you have accounts on multiple computers you can use the SSH keys to
log in. Just copy the **public** key file (ending in ``.pub``) to
``~/.ssh/authorized_keys`` on the remote computer and make sure that
the file is only read/writeable by yourself. Voila, the next time you
ssh into that machine you don't have to provide your password.


.. _section-trac-ssh-key:

Manually Linking your Public Key to your Trac Account
-----------------------------------------------------

The Sage trac server needs to know one of your public keys. You can
upload it in the preferences, that is

1. Go to http://trac.sagemath.org

2. Log in with your trac username/password

3. Click on "Preferences"

4. Go to the "SSH Keys" tab

5. Paste the content of your public key file
   (e.g. ``~/.ssh/id_rsa.pub``)

6. Click on "Save changes"

Note that this does **not** allow you to ssh into any account on trac,
it is only used to authenticate you to the gitolite installation on
trac. You can test that you are being authenticated correctly by
issuing some basic gitolite commands, for example::

    [user@localhost ~]$ ssh git@trac.sagemath.org info
    hello user, this is git@trac running gitolite3 (unknown) on git 1.7.9.5

     R W      sage
    [user@localhost ~]$ ssh git@trac.sagemath.org help
    hello user, this is gitolite3 (unknown) on git 1.7.9.5

    list of remote commands available:

        desc
        help
        info
        perms
        writable


Reporting Bugs
==============

If you think you have found a bug in Sage, you should first search
through the following Google groups for postings related to your
possible bug:

* ``sage-devel``: http://groups.google.com/group/sage-devel
* ``sage-support``: http://groups.google.com/group/sage-support

Maybe the problem you have encountered has already been discussed. You
should also search the trac server to see if anyone else has opened a
ticket about your bug.

If you do not find anything, and you are not sure that you have found
a bug, ask about it on ``sage-devel``. You might be asked to open a
new ticket on the trac server. As mentioned above, you need an account
to do this. To report a bug, login and click on the "New ticket"
button. Type a meaningful one-liner in the "Short summary" box, with
more information in the larger box below. You should include at least
one explicit, reproducible example illustrating your bug (and/or the
steps required to reproduce the buggy behavior). You should also
include the version of Sage (and any relevant packages) you are using,
and operating system information, being precise as possible (32-bit,
64-bit, ...).

Between the "Summary" and "Full description" boxes, there is a
place to choose the "Type" of the ticket: "Defect", "Enhancement",
or "Task". Use your best judgment here; a bug should probably be
reported as a "Defect".

Also pick a component for your bug; this is sometimes
straightforward. If your bug deals with Sage's calculus
implementation, choose "calculus". If it is not obvious, do your
best. Choose a milestone; if you are not sure what to choose, just
choose the numbered version of Sage from the menu ("sage-5.10", for
example). Type in some helpful keywords. In the box labeled "Assign
to", type "somebody" if you are not sure what else to do.

Hit the "Preview" button to make sure everything looks okay, and
then hit "Submit ticket".

If you do not have an account on the trac system to report directly,
you are still encouraged to report any possible bug to the
``sage-devel`` mailing list at ``sage-devel@googlegroups.com``.
The list is moderated for new users and requires subscription.
In your bug report to ``sage-devel``, make sure to include the
following information:

- **operating system**: as precise as possible and architecture
  (32-bit, 64-bit, ...)

- affected version: the exact **version number** and the downloaded
  package (source, precompiled, virtual machine image, or an upgrade
  from a previous version (which one?))

- provide a **reproducible example** and/or define the steps to
  reproduce the erroneous behaviour.

Thank you in advance for reporting bugs to improve Sage in the future!


Guidelines for Opening Tickets
==============================

In addition to bug reports, you should also open a ticket if you
have some new code which extends Sage's capabilities. If you have a
feature request, start a discussion on ``sage-devel`` first,
and then if there seems to be general agreement that you have a
good idea, open a ticket describing the idea.

When you consider opening a new ticket, please bear the following
points in mind.

- Before opening a ticket, make sure that nobody else has opened a
  ticket about the same or closely related issue.

- It is much better to open several specific tickets than one that
  is very broad. Indeed, a single ticket which deals with lots of
  different issues can be quite problematic, and should be avoided.

- Be precise: If foo does not work on OS X but is fine on Linux,
  mention that in the title. Use the keyword option so that
  searches will pick up the issue.

- The problem described in the ticket must be solvable. For
  example, it would be silly to open a ticket whose purpose was
  "Make Sage the best mathematical software in the world". There is
  no metric to measure this properly and it is highly subjective.

- If appropriate, provide URLs to background information or email
  threads relevant to the problem you are reporting.


.. _section-trac-fields:

The Ticket Fields
=================

When you open a new ticket or change an existing ticket, you will find
a variety of fields that can be changed. Here is a comprehensive
overview:

* **Reported by:** The trac account name of whoever created the
  ticket. Cannot be changed.

* **Owned by:** Trac account name of owner, by default the person in
  charge of the **Component:**. Generally not used in the Sage trac.

* **Priority:** The priority of the ticket. Keep in mind that the
  "blocker" label should be used very sparingly.

* **Milestone:** Milestones are usually goals to be met while working
  toward a release. In Sage’s trac, we use milestones instead of
  releases. Each ticket must have a milestone assigned. If you are
  unsure, assign it to the current milestone.

* **Component:** A list of components of Sage, pick one that most
  closely matches the ticket.

* **Keywords:** List of keywords. Fill in any keywords that you think
  will make your ticket easier to find. Tickets that have been worked
  on at Sage days ``NN`` (some number) ofter have ``sdNN`` as keyword.

* **Cc:** List of trac user names to Cc (send emails for changes on
  the ticket). Note that users that enter a comment are automatically
  substcribed to email updates and don't need to be listed under Cc.

* **Merged in:** The Sage release where the ticket was merged in. Only
  changed by the release manager.

* **Authors:** Real name of the ticket author (or list of authors).

* **Reviewers:** Real name of the ticket reviewer (or list of
  reviewers).

* **Report Upstream:** If the ticket is a bug in an upstream component
  of Sage, this field is used to summarize the communication with the
  upstream developers.

* **Work issues:** Issues that need to be resolved before the ticket
  can leave the "needs work" status.

* **Branch:** See :ref:`section-git-branch`

* **Dependencies:** Does the ticket depend on another ticket?
  Sometimes, a ticket requires that another ticket be applied
  first. If this is the case, put the dependencies as a
  comma-separated list (``#1234, #5678``) into the "Dependencies:"
  field.

* **Stopgaps:** See :ref:`section-trac-stopgaps`.



.. _section-trac-stopgaps:

Stopgaps
========

If a component of Sage produces a mathematical error, you should open
two tickets: a main ticket with all available details, and also a
"stopgap" ticket. This second ticket should have a patch which will be
merged into Sage if no one fixes the main issue; this patch should print a
warning when anyone uses the relevant code. To produce the warning
message, use code like the following::

    from sage.misc.stopgap import stopgap
    stopgap("This code contains bugs and may be mathematically unreliable.",
        TICKET_NUM)

Replace ``TICKET_NUM`` by the ticket number for the main ticket.  See
:trac:`12699`, for example.  On the main trac ticket, you should also
enter the ticket number for the stopgap ticket in the "Stopgaps"
field. Stopgap tickets should be marked as blockers.

.. note::

    If mathematically valid code causes Sage to raise an error or
    crash, for example, there is no need for a stopgap.  Rather,
    stopgaps are to warn users that they may be using buggy code; if
    Sage crashes, this is not an issue.


Working on Tickets
==================

If you manage to fix a bug or enhance Sage you are our hero. See
:ref:`chapter-walk-through` for making changes to the Sage source
code, uploading them to the Sage trac server, and finally putting your
new branch on the trac ticket. The following are some other relevant
issues:

* The Patch buildbot wil automatically test your ticket. See `the
  patchbot wiki <http://wiki.sagemath.org/buildbot>`_ for more
  information about its features and limitations. Make sure that you
  look at the log, especially if the patch buildbot did not give you
  the green blob.

* Every bug fixed should result in a doctest.

* This is not an issue with defects, but there are many enhancements
  possible for Sage and too few developers to implement all the good
  ideas. The trac server is useful for keeping ideas in a central
  place because in the Google groups they tend to get lost once they
  drop off the first page.

* If you are a developer, be nice and try to solve a stale/old ticket
  every once in a while.

* Some people regularly do triage. In this context, this means that we
  look at new bugs and classify them according to our perceived
  priority. It is very likely that different people will see
  priorities of bugs very differently from us, so please let us know
  if you see a problem with specific tickets.


.. _section-review-patches:

Reviewing Patches
=================

All code that goes into Sage is peer-reviewed, to ensure that the
conventions discussed in this manual are followed, to make sure that
there are sufficient examples and doctests in the documentation, and
to try to make sure that the code does, mathematically, what it is
supposed to.

If someone (other than you) has posted a git branch for a ticket on
the trac server, you can review it! Look at the branch diff (by
clicking on the ) to see if it makes sense.  Download it (see
:ref:`section-walkthrough-review`) and build Sage with the new
code. Now ask yourself questions such as the following:

- Does the new source code make sense?

- When you run it in Sage, does it fix the problem reported on the
  ticket?

- Does it introduce any new problems?

- Is it documented sufficiently, including both explanation and
  doctests? All code in Sage must have doctests, so if the ticket
  author changes code which did not have a doctest before, the new
  version must include one. In particular, all new code must be 100%
  doctested. Use the command ``sage -coverage <files>`` to see the
  coverage percentage of ``<files>``.

- In particular, is there a doctest illustrating that the bug has been
  fixed? If a function used to give the wrong answer and this ticket
  fixes that, then it should include a doctest illustrating its new
  success.  The surrounding docstring shoud contain the ticket number,
  for example ``See :trac:`12345```.

- If the ticket claims to speed up some computation, does the ticket
  contain code examples to illustrate the claim? The ticket should
  explain the speed efficiency before applying the patch. It should
  also explain the speed efficiency gained after applying the patch.

- Does the reference manual build without errors? You can test the
  reference manual using the command ``sage -docbuild reference html``
  to build the HTML version. The PDF version of the reference manual
  must also build without errors. Use the command ``sage -docbuild
  reference pdf`` to test it out. The latter command requires that you
  have LaTeX installed on your system.

- Do all doctests pass without errors? It is difficult to predict
  which components of Sage will be affected by a given patch and you
  should run tests on the whole library---including those flagged as
  ``#long``---before giving a positive review. You can test the Sage
  library with ``make ptestlong``. See :ref:`chapter-doctesting` for
  more information.

- Do the code and documentation follow conventions documented in the
  following sections?

  - :ref:`chapter-code-basics`
  - :ref:`chapter-python`
  - :ref:`chapter-cython`

If the answers to these and other such reasonable questions are yes,
then you might want to give the patch a positive review. On the main
ticket page, write a comment in the box explaining your review. If you
don't feel experienced enough for this, make a comment explaining what
you checked, and end by asking if someone more experienced will take a
look.  If you think there are issues with the patch, explain them in
the comment box and change the status to "needs work". Browse the
tickets on the trac server to see how things are done.

.. note::

    "The perfect is the enemy of the good"

    The point of the review is to ensure that the Sage code guidelines
    are followed and that the the implementation is mathematically
    correct. Please refrain from aditional feature requests or
    open-ended discussion about alternative implementations. If you
    want the patch written differently, your suggestion should be a
    clear and actionable request.



Closing Tickets
===============

Only the Sage release manager will close tickets. Most likely, this is
not you nor will your trac account have the necessary permissions. If
you feel strongly that a ticket should be closed or deleted, then
change the status of the ticket to *needs review* and change the
milestone to *sage-duplictate/invalid/wontfix*. You should also
comment on the ticket, explaining why it should be closed. If another
developer agrees, he sets the ticket to *positive review*.

A related issue is re-opening tickets. You should refrain from
re-opening a ticket that is already closed. Instead, open a new ticket
and provide a link in the description to the old ticket.


Reasons to Invalidate Tickets
=============================

**One Issue Per Ticket**: A ticket must cover only one issue
and should not be a laundry list of unrelated issues. If a ticket
covers more than one issue, we cannot close it and while some of
the patches have been applied to a given release, the ticket would
remain in limbo.

**No Patch Bombs**: Code that goes into Sage is peer-reviewed. If
you show up with an 80,000 lines of code bundle that completely
rips out a subsystem and replaces it with something else, you can
imagine that the review process will be a little tedious. These
huge patch bombs are problematic for several reasons and we prefer
small, gradual changes that are easy to review and apply. This is
not always possible (e.g. coercion rewrite), but it is still highly
recommended that you avoid this style of development unless there
is no way around it.

**Sage Specific**: Sage's philosophy is that we ship everything
(or close to it) in one source tarball to make debugging possible.
You can imagine the combinatorial explosion we would have to deal
with if you replaced only ten components of Sage with external
packages. Once you start replacing some of the more essential
components of Sage that are commonly packaged (e.g. Pari, GAP,
lisp, gmp), it is no longer a problem that belongs in our tracker.
If your distribution's Pari package is buggy for example, file a
bug report with them. We are usually willing and able to solve
the problem, but there are no guarantees that we will help you
out. Looking at the open number of tickets that are Sage specific,
you hopefully will understand why.

**No Support Discussions**: The trac installation is not meant to
be a system to track down problems when using Sage. Tickets should
be clearly a bug and not "I tried to do X and I couldn't get it to
work. How do I do this?" That is usually not a bug in Sage and it
is likely that ``sage-support`` can answer that question for you. If
it turns out that you did hit a bug, somebody will open a concise
and to-the-point ticket.

**Solution Must Be Achievable**: Tickets must be achievable. Many
times, tickets that fall into this category usually ran afoul to
some of the other rules listed above. An example would be to
"Make Sage the best CAS in the world". There is no metric to
measure this properly and it is highly subjective.


