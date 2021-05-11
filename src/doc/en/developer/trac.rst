.. highlight:: shell-session

.. _chapter-sage-trac:

====================
The Sage Trac Server
====================

All changes to Sage source code have to go through the `Sage Trac
development server <https://trac.sagemath.org>`_. The purpose
of the Sage trac server is to

1. Provide a place for discussion on issues and store a permanent
   record.

2. Provide a repository of source code and all proposed changes.

3. Link these two together.

There is also a `wiki <https://trac.sagemath.org/wiki>`_ for more general
organizational web pages, like Sage development workshops.

Thus if you find a bug in Sage, if you have new code to submit, want
to review new code already written but not yet included in Sage, or if
you have corrections for the documentation, you should post on the
trac server. Items on the server are called *tickets*, and anyone may
search or browse the tickets. For a list of recent changes, just visit
the `Sage trac timeline <https://trac.sagemath.org/timeline>`_.

.. _section-trac-account:

Obtaining an Account
====================

**New:** Previously, it was necessary to manually request a Trac account in
order to post anything to Sage's Trac.  Now, if you have a GitHub account, you
may log in using it to create and comment on tickets, and edit wiki pages on
Sage's Trac.

A manual account request is currently only necessary if you prefer not to
use GitHub or if you want to log into the old `Sage Wiki
<https://wiki.sagemath.org>`_.  This may change as well in the future.

To obtain a non-GitHub account, send an email to
``sage-trac-account@googlegroups.com`` containing:

* your full name,
* preferred username,
* contact email,
* and reason for needing a trac account

Your trac account also grants you access to the `sage wiki
<https://wiki.sagemath.org>`_. Make sure you understand the review process, and
the procedures for opening and closing tickets before making changes. The
remainder of this chapter contains various guidelines on using the trac server.

Trac authentication through SSH
===============================

There are two avenues to prove to the trac server that you are who you
claim to be. First, to change the ticket web pages you need to log in
to trac using a username/password. Second, there is public key
cryptography used by git when copying new source files to the
repository. This section will show you how to set up both.

Generating and Uploading your SSH Keys
--------------------------------------

The git installation on the development server uses SSH keys to decide if and
where you are allowed to upload code. No SSH key is required to report a bug or
comment on a ticket, but as soon as you want to contribute code yourself you
need to provide trac with the public half of your own personal key.
Details are described in the following two sections.


Generating your SSH Keys
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

Linking your Public Key to your Trac Account
-----------------------------------------------------

The Sage trac server needs to know one of your public keys. You can
upload it in the preferences, that is

1. Go to https://trac.sagemath.org

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

.. _trac-bug-report:

Reporting Bugs
==============

If you think you have found a bug in Sage, here is the procedure:

- Search through our Google groups for postings related to your possible bug (it
  may have been fixed/reported already):

  * ``sage-devel``: `<https://groups.google.com/group/sage-devel>`_
  * ``sage-support``: `<https://groups.google.com/group/sage-support>`_

  Similarly, you can search :ref:`chapter-sage-trac` to see if anyone else has
  opened a ticket about your bug.

- If you do not find anything, and you are not sure that you have found a bug,
  ask about it on `sage-devel <https://groups.google.com/group/sage-devel>`_. A
  bug report should contain:

  - An explicit and **reproducible example** illustrating your bug (and/or the
    steps required to reproduce the buggy behavior).

  - The **version** of Sage you run, as well as the version of the optional
    packages that may be involved in the bug.

  - Describe your **operating system** as accurately as you can and your
    architecture (32-bit, 64-bit, ...).

- You might be asked to open a new ticket. In this case, follow the
  :ref:`section-trac-new-ticket`.

Thank you in advance for reporting bugs to improve Sage in the future!

.. _section-trac-new-ticket:

Guidelines for Opening Tickets
==============================

In addition to bug reports (see :ref:`trac-bug-report`), you should also open a
ticket if you have some new code that makes Sage a better tool. If you have a
feature request, start a discussion on ``sage-devel`` first, and then if there
seems to be general agreement that you have a good idea, open a ticket
describing the idea.

- Do you already have a **trac account**? If not, :ref:`click here
  <section-trac-account>`.

**Before** opening a new ticket, consider the following points:

- Make sure that nobody else has opened a ticket about the same or closely
  related issue.

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

- For bug reports: the ticket's description should contain the information
  described at :ref:`trac-bug-report`.

- If appropriate, provide URLs to background information or sage-devel
  conversation relevant to the problem you are reporting.

**When creating** the ticket, you may find useful to read
:ref:`section-trac-fields`.

Unless you know what you are doing, leave the milestone field to its default
value.

.. _section-trac-fields:

The Ticket Fields
=================

When you open a new ticket or change an existing ticket, you will find a variety
of fields that can be changed. Here is a comprehensive overview (for the
'status' field, see :ref:`section-trac-ticket-status`):

* **Reported by:** The trac account name of whoever created the
  ticket. Cannot be changed.

* **Owned by:** Trac account name of owner, by default the person in charge of
  the Component (see below). Generally not used in the Sage trac.

* **Type:** One of ``enhancement`` (e.g. a new feature), ``defect`` (e.g. a bug
  fix), or ``task`` (rarely used).

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

* **Authors:** Real name of the ticket author(s).

* **Reviewers:** Real name of the ticket reviewer(s).

* **Report Upstream:** If the ticket is a bug in an upstream component
  of Sage, this field is used to summarize the communication with the
  upstream developers.

* **Work issues:** Issues that need to be resolved before the ticket
  can leave the "needs work" status.

* **Branch:** The Git branch containing the ticket's code (see
  :ref:`section-walkthrough-branch`). It is displayed in green color,
  unless there is a conflict between the branch and the latest beta
  release (red color). In this case, the branch should be merged or
  rebased.

* **Dependencies:** Does the ticket depend on another ticket?
  Sometimes, a ticket requires that another ticket be applied
  first. If this is the case, put the dependencies as a
  comma-separated list (``#1234, #5678``) into the "Dependencies:"
  field.

* **Stopgaps:** See :ref:`section-trac-stopgaps`.

.. _section-trac-ticket-status:

The status of a ticket
======================

The status of a ticket appears right next to its number, at the top-left corner
of its page. It indicates who has to work on it.

- **new** -- the ticket has only been created (or the author forgot to change
  the status to something else).

  If you want to work on it yourself it is better to leave a comment to say
  so. It could avoid having two persons doing the same job.

- **needs_review** -- the code is ready to be peer-reviewed. If the code is not
  yours, then you can review it. See :ref:`chapter-review`.

- **needs_work** -- something needs to be changed in the code. The reason should
  appear in the comments.

- **needs_info** -- somebody has to answer a question before anything else can
  happen. It should be clear from the comments.

- **positive_review** -- the ticket has been reviewed, and the release manager
  will close it.

The status of a ticket can be changed using a form at the bottom of the ticket's
page. Leave a comment explaining your reasons whenever you change it.

.. _section-trac-stopgaps:

Stopgaps
========

When Sage returns wrong results, two tickets should be opened:

- A main ticket with all available details.
- A "stopgap" ticket (e.g. :trac:`12699`)

This second ticket does not fix the problem but adds a warning that will be
printed whenever anyone uses the relevant code. This, until the problem is
finally fixed.

To produce the warning message, use code like the following:

.. CODE-BLOCK:: python

    from sage.misc.stopgap import stopgap
    stopgap("This code contains bugs and may be mathematically unreliable.",
        TICKET_NUM)

Replace ``TICKET_NUM`` by the ticket number for the main ticket. On the main
trac ticket, enter the ticket number for the stopgap ticket in the "Stopgaps"
field (see :ref:`section-trac-fields`). Stopgap tickets should be marked as
blockers.

.. NOTE::

    If mathematically valid code causes Sage to raise an error or
    crash, for example, there is no need for a stopgap.  Rather,
    stopgaps are to warn users that they may be using buggy code; if
    Sage crashes, this is not an issue.


Working on Tickets
==================

If you manage to fix a bug or enhance Sage you are our hero. See
:ref:`chapter-walkthrough` for making changes to the Sage source
code, uploading them to the Sage trac server, and finally putting your
new branch on the trac ticket. The following are some other relevant
issues:

* The Patch buildbot will automatically test your ticket. See `the
  patchbot wiki <https://wiki.sagemath.org/buildbot>`_ for more
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

Reviewing and closing Tickets
=============================

Tickets can be closed when they have positive review or for other reasons. To
learn how to review, please see :ref:`chapter-review`.

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


