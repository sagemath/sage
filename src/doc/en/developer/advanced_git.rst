.. _chapter-advanced-git:

============
Advanced Git
============

This chapter covers some advanced uses of git that go beyond what is
required to work with branches. These features can be used in Sage
development, but are not really necessary to contribute to Sage. If
you are just getting started with Sage development, you should read
:ref:`chapter-walk-through` instead. If you are new to git, please see
:ref:`chapter-manual-git`.


Detached Heads and Reviewing Tickets
====================================

Each commit is a snapshot of the Sage source tree at a certain
point. So far, we always used commits organized in branches. But
secretly the branch is just a shortcut for a particular commit, the
head commit of the branch. But you can just go to a particular commit
without a branch, this is called "detached head". If you have the
commit already in your local history, you can directly check it
out without requiring internet access::

    [user@localhost sage]$ git checkout a63227d0636e29a8212c32eb9ca84e9588bbf80b
    Note: checking out 'a63227d0636e29a8212c32eb9ca84e9588bbf80b'.

    You are in 'detached HEAD' state. You can look around, make experimental
    changes and commit them, and you can discard any commits you make in this
    state without impacting any branches by performing another checkout.

    If you want to create a new branch to retain commits you create, you may
    do so (now or later) by using -b with the checkout command again. Example:

      git checkout -b new_branch_name

    HEAD is now at a63227d... Szekeres Snark Graph constructor

If it is not stored in your local git repository, you need to download
it from the trac server first::

    [user@localhost sage]$ git fetch trac a63227d0636e29a8212c32eb9ca84e9588bbf80b
    From ssh://trac/sage
     * branch            a63227d0636e29a8212c32eb9ca84e9588bbf80b -> FETCH_HEAD
    [user@localhost sage]$ git checkout FETCH_HEAD
    HEAD is now at a63227d... Szekeres Snark Graph constructor

Either way, you end up with your current HEAD and working directory
that is not associated to any local branch::

    [user@localhost sage]$ git status
    # HEAD detached at a63227d
    nothing to commit, working directory clean

This is perfectly fine. You can switch to an existing branch (with the
usual ``git checkout my_branch``) and back to your detached head.

Detached heads can be used to your advantage when reviewing
tickets. Just check out the commit (look at the "Commit:" field on the
trac ticket) that you are reviewing as a detached head. Then you can
look at the changes and run tests in the detached head. When you are
finished with the review, you just abandon the detached head. That way
you never create a new local branch, so you don't have to type ``git
branch -D my_branch`` at the end to delete the local branch that you
created only to review the ticket.


.. _section-git-recovery:

Reset and Recovery
==================

Git makes it very hard to truly mess up. Here is a short way to get
back onto your feet, no matter what. First, if you just want to go
back to a working Sage installation you can always abandon your
working branch by switching to your local copy of the ``master``
branch::

    [user@localhost sage]$ git checkout master

As long as you did not make any changes to the ``master`` branch
directly, this will give you back a working Sage.

If you want to keep your branch but go back to a previous commit you
can use the *reset* command. For this, look up the commit in the log
which is some 40-digit hexadecimal number (the SHA1 hash). Then use
``git reset --hard`` to revert your files back to the previous state::

    [user@localhost sage]$ git log
    ...
    commit eafaedad5b0ae2013f8ae1091d2f1df58b72bae3
    Author: First Last <user@email.com>
    Date:   Sat Jul 20 21:57:33 2013 -0400

        Commit message
    ...
    [user@localhost sage]$ git reset --hard eafae

.. warning::

    Any *uncommitted* changes will be lost!

You only need to type the first couple of hex digits, git will
complain if this does not uniquely specify a commit. Also, there is
the useful abbreviation ``HEAD~`` for the previous commit and
``HEAD~n``, with some integer ``n``, for the n-th previous commit.

Finally, perhaps the ultimate human error recovery tool is the
reflog. This is a chronological history of git operations that you can
undo if needed. For example, let us assume we messed up the *git
reset* command and went back too far (say, 5 commits back). And, on
top of that, deleted a file and committed that::

    [user@localhost sage]$ git reset --hard HEAD~5
    [user@localhost sage]$ git rm sage
    [user@localhost sage]$ git commit -m "I shot myself into my foot"

Now we cannot just checkout the repository from before the reset,
because it is no longer in the history. However, here is the reflog::

    [user@localhost sage]$ git reflog
    2eca2a2 HEAD@{0}: commit: I shot myself into my foot
    b4d86b9 HEAD@{1}: reset: moving to HEAD~5
    af353bb HEAD@{2}: checkout: moving from some_branch to master
    1142feb HEAD@{3}: checkout: moving from other_branch to some_branch
    ...

The ``HEAD@{n}`` revisions are shortcuts for the history of git
operations. Since we want to rewind to before the erroneous *git
reset* command, we just have to reset back into the future::

    [user@localhost sage]$ git reset --hard HEAD@{2}
    


.. _section-git-rewriting-history:

Rewriting History
=================

Git allows you to rewrite history, but be careful: the SHA1 hash of a
commit includes the parent's hash. This means that the hash really
depends on the entire content of the working directory; every source
file is in exactly the same state as when the hash was computed. This
also means that you can't change history without modifying the
hash. If others branched off your code and then you rewrite history,
then the others are thoroughly screwed. So, ideally, you would only
rewrite history on branches that you have not yet pushed to trac.

As an advanced example, consider three commits A, B, C that were made
on top of each other. For simplicity, we'll assume they just added a
file named ``file_A.py``, ``file_B.py``, and ``file_C.py`` ::

    [user@localhost]$ git log --oneline
    9621dae added file C
    7873447 added file B
    bf817a5 added file A
    5b5588e base commit

Now, let's assume that the commit B was really independent and ought
to be on a separate ticket. So we want to move it to a new branch,
which we'll call ``second_branch``. First, branch off at the base
commit before we added A::

    [user@localhost]$ git checkout 5b5588e
    Note: checking out '5b5588e'.

    You are in 'detached HEAD' state. You can look around, make experimental
    changes and commit them, and you can discard any commits you make in this
    state without impacting any branches by performing another checkout.

    If you want to create a new branch to retain commits you create, you may
    do so (now or later) by using -b with the checkout command again. Example:

      git checkout -b new_branch_name

    HEAD is now at 5b5588e... base commit
    [user@localhost]$ git checkout -b second_branch
    Switched to a new branch 'second_branch'
    [user@localhost]$ git branch
      first_branch
    * second_branch
    [user@localhost]$ git log --oneline
    5b5588e base commit

Now, we make a copy of commit B in the current branch::

    [user@localhost]$ git cherry-pick 7873447
    [second_branch 758522b] added file B
     1 file changed, 1 insertion(+)
     create mode 100644 file_B.py
    [user@localhost]$ git log --oneline
    758522b added file B
    5b5588e base commit

Note that this changes the SHA1 of the commit B, since its parent
changed! Also, cherry-picking *copies* commits, it does not remove
them from the source branch. So we now have to modify the first branch
to exclude commit B, otherwise there will be two commits adding
``file_B.py`` and our two branches would conflict later when they are
being merged into Sage. Hence, we first reset the first branch back to
before B was added::

    [user@localhost]$ git checkout first_branch 
    Switched to branch 'first_branch'
    [user@localhost]$ git reset --hard bf817a5
    HEAD is now at bf817a5 added file A

Now we still want commit C, so we cherry-pick it again. Note that this
works even though commit C is, at this point, not included in any
branch::

    [user@localhost]$ git cherry-pick 9621dae
    [first_branch 5844535] added file C
     1 file changed, 1 insertion(+)
     create mode 100644 file_C.py
    [user@localhost]$ git log --oneline
    5844535 added file C
    bf817a5 added file A
    5b5588e base commit

And, again, we note that the SHA1 of commit C changed because its
parent changed. Voila, now you have two branches where the first
contains commits A, C and the second contains commit B.


.. _section-git-interactive-rebase:

Interactively Rebasing
======================

An alternative approach to :ref:`section-git-rewriting-history` is to
use the interactive rebase feature. This will open an editor where you
can modify the most recent commits. Again, this will naturally modify
the hash of all changed commits and all of their children.

Now we start by making an identical branch to the first branch::

    [user@localhost]$ git log --oneline
    9621dae added file C
    7873447 added file B
    bf817a5 added file A
    5b5588e base commit
    [user@localhost]$ git checkout -b second_branch
    Switched to a new branch 'second_branch'
    [user@localhost]$ git rebase -i HEAD~3
    
This will open an editor with the last 3 (corresponding to ``HEAD~3``)
commits and instuctions for how to modify them::

    pick bf817a5 added file A
    pick 7873447 added file B
    pick 9621dae added file C
    
    # Rebase 5b5588e..9621dae onto 5b5588e
    #
    # Commands:
    #  p, pick = use commit
    #  r, reword = use commit, but edit the commit message
    #  e, edit = use commit, but stop for amending
    #  s, squash = use commit, but meld into previous commit
    #  f, fixup = like "squash", but discard this commit's log message
    #  x, exec = run command (the rest of the line) using shell
    #
    # These lines can be re-ordered; they are executed from top to bottom.
    #
    # If you remove a line here THAT COMMIT WILL BE LOST.
    #
    # However, if you remove everything, the rebase will be aborted.
    #
    # Note that empty commits are commented out
   
To only use commit B, we delete the first and third line. Then save
and quit your editor, and your branch now consists only of the B commit.

You still have to delete the B commit from the first branch, so you
would go back (``git checkout first_branch``) and then run the same
``git rebase -i`` command and delete the B commit.
 
