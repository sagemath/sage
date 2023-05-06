.. nodoctest

.. _chapter-review:

==============
Reviewing Code
==============

All code that goes into Sage is peer-reviewed. Two reasons for this are:

- Because a developer cannot think of everything at once
- Because a fresh pair of eyes may spot a mathematical error,
  a corner-case in the code, insufficient documentation, a missing
  consistency check, etc.

Anybody (e.g. you) can do this job for somebody else's PR. This document
lists things that the reviewer must check before deciding that a PR is
ready for inclusion into Sage.

You can now begin the review by reading the diff code.

**Read the diff:**

**Build the code:** while you read the code, you can :ref:`rebuild Sage with the
new code <section-walkthrough-make>`. If you do not know how to **download the
code**, :ref:`click here <section-workflows-pr-checkout>`.

The following should generally be checked while reading and testing the code:

- **The purpose**: Does the code address the PR's stated aim? Can it
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

- **Bugfixes**: If the PR contains a bugfix, does it add a doctest
  illustrating that the bug has been fixed? This new doctest should contain the
  issue or PR number, for example ``See :issue:`12345```.

- **Speedup**: Can the PR make any existing code slower? if the PR
  claims to speed up some computation, does the PR contain code examples to
  illustrate the claim? The PR should explain how the speedup is achieved.

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
:ref:`section-github-pr-status`):

- **positive review**: If the answers to the questions above and other
  reasonable questions are *"yes"*, you can set the PR to
  ``positive_review``.

- **needs_work**: If something is not as it should, write a list of all points
  that need to be addressed in a comment and change the PR's status to
  ``needs_work``.

- **needs_info**: If something is not clear to you and prevents you from going
  further with the review, ask your question and set the PR's status to
  ``needs_info``.

- If you **do not know what to do**, for instance if you don't feel experienced
  enough to take a final decision, explain what you already did in a comment and
  ask if someone else could take a look.

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

.. _section-workflows-review:

Reviewing a change
==================

Use the checks on GitHub Actions.

Use `pull request reviews
<https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/reviewing-changes-in-pull-requests/about-pull-request-reviews>`_.
You can add comments directly to changed lines.

Smaller suggestions can be made `through the github web interface
<https://egghead.io/lessons/github-add-suggestions-in-a-github-pr-review>`_ as
part of the `pull request review
<https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/reviewing-changes-in-pull-requests/incorporating-feedback-in-your-pull-request>`_.

If you want to be able to make changes directly to others' PRs** (when the
author selects to allow edits from maintainers), please contact one of the
`Sagemath Github admins
<https://github.com/orgs/sagemath/people?query=role%3Aowner>`_, who can give
you the relevant permissions.

Change the status of the PR (e.g., "positive review" or "needs work"), choose
the `correct type of your pull request review
<https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/reviewing-changes-in-pull-requests/reviewing-proposed-changes-in-a-pull-request#submitting-your-review>`_
(i.e. "approve" vs "request changes")

For trying the branch of a PR locally, use::

    git fetch upstream pull/PULL-REQUEST-ID/head:LOCAL-BRANCH-NAME

Consult https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/reviewing-changes-in-pull-requests/checking-out-pull-requests-locally.

Alternatively,::

    gh pr checkout PULL-REQUEST-ID

Consult https://cli.github.com/manual/gh_pr_checkout.


Reviewing and Closing PRs
=========================

PRs can be closed when they have positive review or for other reasons. To
learn how to review, please see :ref:`chapter-review`.

Use one of the *resolution labels* ``r: duplicate``, ``r: invalid``, ``r: wontfix``. Add a comment explaining why the issue has been closed if that's not already clear from the discussion.

Users with the necessary permissions can then directly `close the issue <https://docs.github.com/en/issues/tracking-your-work-with-issues/closing-an-issue>`_: In the dropdown menu on the "Close issue" button, select "Close as not planned"

Otherwise, use the labels "needs review" or "positive review", and someone else with the necessary rights will take care of closing the issue.

If you think an issue has been prematurely be closed, feel free to reopen it.


Reasons to Invalidate PRs
=========================

**One Issue Per One Issue**: An issue must cover only one issue
and should not be a laundry list of unrelated issues. If an issue
covers more than one issue, we cannot close it and while some of
the patches have been applied to a given release, the issue would
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

The Release Process
===================

The Sage Release Manager uses the following procedure to make releases, as of
2022.

**Beta Release Stage**: For preparing a new beta release or the first release
candidate, all positively reviewed tickets with the forthcoming release
milestone are considered. Tickets that have unmerged dependencies are ignored.
The Release Manager merges tickets in batches of 10 to 20 tickets, taking the
ticket priority into account. If a merge conflict of a ticket to the Release
Manager's branch occurs, the ticket is set back to "needs work" status by the
Release Manager, and the list of the tickets already merged to the Release
Manager's branch is posted. The author of the ticket needs to identify
conflicting tickets in the list, make merge commits and declare them as
dependencies, before setting back to "positive review" status. Each batch of
merged tickets then undergoes integration testing. If problems are detected, a
ticket will be set back to "needs work" status and unmerged. When a batch of
tickets is ready, the Release Manager closes these tickets and proceeds to the
next batch. After a few batches, a new beta release is tagged, pushed to the
``develop`` branch on the main git repository, and announced on
``sage-release``.

**Release Candidate Stage**: After the first release candidate has been made,
the project is in the release candidate stage, and a modified procedure is
used. Now **only tickets with a priority set to "blocker" are considered**.
Tickets with all other priorities, including "critical", are ignored. Hence if
a ticket is important enough to merit inclusion in this stage, it should be set
to "blocker".

**Blocker Tickets**: The goal of the release process is to make a stable
release of high quality. Be aware that there is a risk/benefit trade-off in
merging a ticket. The benefit of merging a ticket is the improvement that the
ticket brings, such as fixing a bug. However, any code change has a risk of
introducing unforeseen new problems and thus delaying the release: If a new
issue triggers another release candidate, it delays the release by 1-2 weeks.
Hence developers should use "blocker" priority sparingly and should indicate
the rationale on the ticket. Though there is no one fixed rule or authority
that determines what is appropriate for "blocker" status,

- Tickets introducing new features are usually not blockers -- unless perhaps
  they round out a set of features that were the focus of development of this
  release cycle.

- Tickets that make big changes to the code, for example refactoring tickets,
  are usually not blockers.

**Final Release**: If there is no blocker ticket for the last release
candidate, the Release Manager turns it to the final release. It is tagged with
the release milestone, and announced on ``sage-release``.

