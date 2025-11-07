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

**Check the GitHub checks:** We require all checks have passed.

**Read the diff:** Click "Files changed" tab of the PR. Read through the
changes of all modified files. We use `pull request reviews
<https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/reviewing-changes-in-pull-requests/about-pull-request-reviews>`_.
You can add comments directly to changed lines.

**Test the code:** (This is optional if the **build and test** check has passed.)
Checkout the **code of the PR** following :ref:`see here <section-workflows-pr-checkout>`
and run the relevant tests.

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
  <chapter-code-basics>`? :ref:`Python's conventions <chapter-python>`?
  :ref:`Cython's conventions <chapter-cython>`?

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

For changes that affect the **user interface**, in particular, upgrades to
IPython and Jupyter component packages, manual testing is crucial because
our automatic tests do not cover the user interface. We recommend to use
a `Jupyter notebook with comprehensive tests of graphics and typesetting
<https://github.com/egourgoulhon/SageMathTest/blob/master/Notebooks/test_display.ipynb>`_,
some of which is Sage-specific.

You are now ready to change the PR's status (see
:ref:`section-github-pr-status`):

- **positive review**: If the answers to the questions above and other
  reasonable questions are *"yes"*, you can set the PR to
  ``positive review`` status.

- **needs work**: If something is not as it should, write a list of all points
  that need to be addressed in a comment and change the PR's status to
  ``needs work`` status.

- **needs info**: If something is not clear to you and prevents you from going
  further with the review, ask your question and set the PR's status to
  ``needs info`` status.

- If you **do not know what to do**, for instance if you don't feel experienced
  enough to take a final decision, explain what you already did in a comment and
  ask if someone else could take a look.

For more advice on reviewing, see `How to Referee Sage Trac Tickets
<http://sagemath.blogspot.com/2010/10/how-to-referee-sage-trac-tickets.html>`_
(caveat: mercurial was replaced with Git and Trac with GitHub).

.. NOTE::

    "The perfect is the enemy of the good"

    The point of the review is to ensure that the Sage code guidelines
    are followed and that the implementation is mathematically
    correct. Please refrain from additional feature requests or
    open-ended discussion about alternative implementations. If you
    want the code written differently, your suggestion should be a
    clear and actionable request.


Reviewing and closing PRs
=========================

PRs can be closed when they have positive review or for other reasons.

If a PR is closed for a reason other than positive review, use one of the
**resolution** labels ``r: duplicate``, ``r: invalid``, ``r: wontfix``, and
``r: worksforme``. Add a comment explaining why the issue has been closed if
that's not already clear from the discussion.

If you think an issue has been prematurely be closed, feel free to reopen it.


Reasons to invalidate PRs
=========================

**One Issue Per One Issue**: An issue must cover only one issue
and should not be a laundry list of unrelated issues. If an issue
covers more than one issue, we cannot close it and while some of
the patches have been applied to a given release, the issue would
remain in limbo.

**No Patch Bombs**: Code that goes into Sage is peer-reviewed. If
you show up with 80,000 lines of code bundle that completely
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
out. Looking at the open number of PRs that are Sage specific,
you hopefully will understand why.

**No Support Discussions**: GitHub is not meant to
be a system to track down problems when using Sage. An issue should
be clearly a bug and not "I tried to do X and I couldn't get it to
work. How do I do this?" That is usually not a bug in Sage and it
is likely that ``sage-support`` can answer that question for you. If
it turns out that you did hit a bug, somebody will open a concise
and to-the-point PR.

**Solution Must Be Achievable**: Issues must be achievable. Many
times, issues that fall into this category usually ran afoul to
some of the other rules listed above. An example would be to
"Make Sage the best CAS in the world". There is no metric to
measure this properly and it is highly subjective.


The release process
===================

It is good for developers and reviewers to be aware of the procedure that the
Sage Release Manager uses to make releases. Here it is as of 2024:

**Beta Release Stage**: For preparing a new beta release or the first release
candidate, all positively reviewed PRs with the forthcoming release
milestone are considered. The Release Manager merges PRs in batches of
10 to 20 PRs. If a merge conflict of a PR to the Release Manager's
branch occurs, the PR is set back to "needs work" status by the
Release Manager. (The author of the PR can try to guess which other
PRs may be causing the conflict, make merge commits and declare them as
dependencies, before setting back to "positive review" status.
Alternatively, the PR author can wait until the next beta release and
resolve the conflict then.) Each batch of
merged PRs undergoes integration testing. If problems are detected, a
PR will be set back to "needs work" status and unmerged. When a batch of
PRs is ready, the Release Manager closes these PRs and proceeds to the
next batch. After a few batches, a new beta release is tagged, pushed to the
``develop`` branch on the Sage repository on GitHub, and announced on
``sage-release``.

**Release Candidate Stage**: After the first release candidate has been made,
the project is in the release candidate stage, and a modified procedure is
used. Now **only PRs with a priority set to "blocker" are considered**.  PRs
with all other priorities, including "critical", are ignored. Hence if a ticket
is important enough to merit inclusion in this stage, it should be set to
"blocker" by adding ``p: blocker / 1`` label.

**Blocker PRs**: The goal of the release process is to make a stable
release of high quality. Be aware that there is a risk/benefit trade-off in
merging a PR. The benefit of merging a PR is the improvement that the
PR brings, such as fixing a bug. However, any code change has a risk of
introducing unforeseen new problems and thus delaying the release: If a new
issue triggers another release candidate, it can delay the release by up to
2 weeks.
Hence developers should use "blocker" priority sparingly and should indicate
the rationale on the PR. Though there is no one fixed rule or authority
that determines what is appropriate for "blocker" status,

- PRs introducing new features are usually not blockers -- unless perhaps
  they round out a set of features that were the focus of development of this
  release cycle.

- PRs that make big changes to the code, for example refactoring PRs,
  are usually not blockers.

**Final Release**: If there is no blocker PR for the last release candidate,
the Release Manager turns it to the final release. It is tagged with the
release milestone, and announced on ``sage-release``.

Release management scripts are maintained in the repository
`sagemath/sage-release-management <https://github.com/sagemath/sage-release-management>`_.
