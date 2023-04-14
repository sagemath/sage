.. nodoctest

.. _chapter-review:

=========================
The reviewer's check list
=========================

.. WARNING::

    **Sage development moved to GitHub in February 2023.** After the transition,
    some parts of this guide (especially those related with `the Sage Trac
    server <https://trac.sagemath.org>`_) became obsolete and need to be
    updated according to the new workflow on GitHub. See our `transition guide
    from Trac to GitHub
    <https://github.com/sagemath/trac-to-github/blob/master/docs/Migration-Trac-to-Github.md>`_
    for the preliminary version of the workflow.

All code that goes into Sage is peer-reviewed. Two reasons for this are:

- Because a developer cannot think of everything at once;
- Because a fresh pair of eyes may spot a mathematical error,
  a corner-case in the code, insufficient documentation, a missing
  consistency check, etc.

Anybody (e.g. you) can do this job for somebody else's ticket. This document
lists things that the reviewer must check before deciding that a ticket is
ready for inclusion into Sage.

- Do you know what the **trac server** is? If not, :ref:`click here
  <chapter-sage-trac>`.

- Do you have a **trac account**? If not, :ref:`click here
  <section-trac-account>`.

You can now begin the review by reading the diff code.

**Read the diff:** the diff (i.e. the ticket's content) can be obtained by
clicking on the (green) branch's name that appears on the trac ticket. If that
name appears in red (see :ref:`section-trac-fields`) you can say so in a comment
and set the ticket to ``needs_work`` (see :ref:`section-trac-ticket-status`).

**Build the code:** while you read the code, you can :ref:`rebuild Sage with the
new code <section-walkthrough-make>`. If you do not know how to **download the
code**, :ref:`click here <section-git_trac-review>` (with git trac) or
:ref:`here <section-git-checkout>` (git only).


The following should generally be checked while reading and testing the code:

- **The purpose**: Does the code address the ticket's stated aim? Can it
  introduce any new problems? Does testing the new or fixed functionality
  with a variety of input, not just the examples in the documentation,
  give expected and robust output (and no unexpected errors or crashes)?

- **User documentation**: Is the use of the new code clear to a user? Are all
  mathematical notions involved standard, or is there explanation (or a link
  to one) provided? Can he/she find the new code easily if he/she needs it?

- **Code documentation**: Is the code sufficiently commented so that a developer
  does not have to wonder what exactly it does?

- **Conventions**: Does the code respect :ref:`Sage's conventions
  <chapter-code-basics>`? :ref:`Python's convention <chapter-python>`?
  :ref:`Cython's convention <chapter-cython>`?

- **Doctest coverage**: Do all functions contain doctests? Use ``sage -coverage
  <files>`` to check it. Are all aspects of the new/modified methods and classes
  tested (see :ref:`section-doctest-writing`)?

- **Bugfixes**: If the ticket contains a bugfix, does it add a doctest
  illustrating that the bug has been fixed? This new doctest should contain the
  ticket number, for example ``See :trac:`12345```.

- **Speedup**: Can the ticket make any existing code slower? if the ticket
  claims to speed up some computation, does the ticket contain code examples to
  illustrate the claim? The ticket should explain how the speedup is achieved.

- **Build the manuals**: Does the reference manual build without
  errors (check both html and pdf)? See :ref:`chapter-sage_manuals` to
  learn how to build the manuals.

- **Look at the manuals**: Does the reference manual look okay? The
  changes may have typos that allow the documentation to build without
  apparent errors but that may cause badly formatted output or broken
  hyperlinks.

- **Run the tests**: Do all doctests pass without errors? Unrelated components
  of Sage may be affected by the change. Check all tests in the whole library,
  including "long" doctests (this can be done with ``make ptestlong``) and any
  optional doctests related to the functionality. See :ref:`chapter-doctesting`
  for more information.

You are now ready to change the ticket's status (see
:ref:`section-trac-ticket-status`):

- **positive review**: If the answers to the questions above and other
  reasonable questions are *"yes"*, you can set the ticket to
  ``positive_review``. Add your full name to the "reviewer" field (see
  :ref:`section-trac-fields`).

- **needs_work**: If something is not as it should, write a list of all points
  that need to be addressed in a comment and change the ticket's status to
  ``needs_work``.

- **needs_info**: If something is not clear to you and prevents you from going
  further with the review, ask your question and set the ticket's status to
  ``needs_info``.

- If you **do not know what to do**, for instance if you don't feel experienced
  enough to take a final decision, explain what you already did in a comment and
  ask if someone else could take a look.

**Reviewer's commit**: if you can fix the issues yourself, you may make a commit
in your own name and mark the commit as a reviewer's patch. To learn how
:ref:`click here <section-git_trac-push>` (git trac) or :ref:`here
<section-git-push>` (git only). This contribution must also be reviewed, for
example by the author of the original patch.

For more advice on reviewing, see [WSblog]_.

.. NOTE::

    "The perfect is the enemy of the good"

    The point of the review is to ensure that the Sage code guidelines
    are followed and that the implementation is mathematically
    correct. Please refrain from additional feature requests or
    open-ended discussion about alternative implementations. If you
    want the patch written differently, your suggestion should be a
    clear and actionable request.

REFERENCES:

.. [WSblog] William Stein, How to Referee Sage Trac Tickets,
   http://sagemath.blogspot.com/2010/10/how-to-referee-sage-trac-tickets.html
   (Caveat: mercurial was replaced with git)
