.. highlight:: shell-session

.. _chapter-git_trac:

=======================================
Collaborative Development with Git-Trac
=======================================

Sometimes you will only want to work on local changes to Sage, for
your own private needs.  However, typically it is beneficial to
share code and ideas with others; the manner in which the
`Sage project <https://www.sagemath.org>`_ does this (as well as fixing
bugs and upgrading components) is in a very collaborative and
public setting on `the Sage Trac server <https://trac.sagemath.org>`_
(the Sage bug and enhancement tracker).

One can use ``git`` :ref:`the hard way <chapter-manual-git>` for this,
but this section explains how to use the helper ``git trac`` command, which
simplifies many of the most common actions in collaboration on Sage. Some
of the :ref:`tutorials <section-git-tutorials>` we suggest may be helpful
in navigating what they are for.

Most of the commands in the following section will not work unless
you have an account on Trac. If you want to contribute to Sage, it
is a good idea to get an account now (see :ref:`section-trac-account`).


.. _section-git_trac-install:

Installing the Git-Trac Command
===============================

Git is a separate project from trac, and the two do not know how to
talk to each other. To simplify the development, we have a special
``git trac`` subcommand for the git suite. Note that this really is
only to simplify interaction with our trac issue management, you can
perform every development task with just git and a web browser. See
:ref:`chapter-manual-git` instead if you prefer to do everything by
hand::

    [user@localhost]$ git clone https://github.com/sagemath/git-trac-command.git
    Cloning into 'git-trac-command'...
    [...]
    Checking connectivity... done.
    [user@localhost]$ source git-trac-command/enable.sh
    Prepending the git-trac command to your search PATH

This creates a directory ``git-trac-command``.

Sourcing the ``enable.sh`` script in there is just a quick and dirty
way to enable it temporarily. For a more permanent installation on
your system later, make sure to put the ``git-trac`` command in your
``PATH``. Assuming that ``~/bin`` is already in your ``PATH``, you can
do this by symlinking::

    [user@localhost]$ echo $PATH
    /home/user/bin:/usr/local/bin:/usr/bin:/bin:/usr/local/sbin:/usr/sbin
    [user@localhost]$ cd git-trac-command
    [user@localhost git-trac-command]$ ln -s `pwd`/git-trac ~/bin/

See the `git-trac README <https://github.com/sagemath/git-trac-command>`_ for
more details. At this point you leave ``git-trac-command`` subdirectory, and only go 
there whenever you need to update the ``git-trac`` command.



.. _section-git_trac-setup:

Git and Trac Configuration
==========================

.. NOTE::

    * `trac <https://trac.sagemath.org>`_ uses username/password for
      authentication.

    * Our `git repository server <https://git.sagemath.org>`_ uses SSH
      public key authentication for write access.

You need to set up both authentication mechanisms to be able to upload
your changes with "git trac". For read-only access neither
authentication mechanism is needed. To set up ``git trac``, first go
to the Sage directory and tell ``git trac`` about your trac account::

    [user@localhost sage]$ git trac config --user USERNAME --pass 'PASSWORD'
    Trac xmlrpc URL:
        https://trac.sagemath.org/xmlrpc (anonymous)
        https://trac.sagemath.org/login/xmlrpc (authenticated)
        realm sage.math.washington.edu
    Username: USERNAME
    Password: PASSWORD
    Retrieving SSH keys...
        1024 ab:1b:7c:c9:9b:48:fe:dd:59:56:1e:9d:a4:a6:51:9d  My SSH Key

where you have to replace USERNAME with your trac user name and
PASSWORD with your trac password. If you don't have a trac account,
use ``git trac config`` without any arguments. The single quotes in
``'PASSWORD'`` escape special characters that you might have in your
password. The password is stored in plain-text in ``.git/config``, so
make sure that it is not readable by other users on your system. For
example, by running ``chmod 0600 .git/config`` if your home directory
is not already private.

Instead of a username and password you may also configure authentication via
a generated token by passing ``--token=<token>`` instead of ``--pass``::

    [user@localhost sage]$ git trac config --user=<username> --token=<token>

This is required if you authenticate to Trac with your GitHub account, as
you do not have a Trac password.  Logged in users can find their token
under `the token tab in preferences on the trac site <https://trac.sagemath.org/prefs/token>`_ .

.. NOTE::

   The username to be entered here is NOT the GitHub username, but rather the trac username which is gh-<GitHub-username>
   as given on the top right corner of the trac server.

If both a token and a username/password are configured, the token-based
authentication takes precedence.

If you do not want to store your trac username/password/token on disk you
can temporarily override it with the environment variables
``TRAC_USERNAME``,  ``TRAC_PASSWORD``, and ``TRAC_TOKEN`` respectively.
These take precedence over any other configuration.

If there is no SSH key listed then you haven't uploaded your SSH
public key to the trac server. You should do that now following the
instructions to :ref:`section-trac-ssh-key`, if you want to upload
any changes. You may have to add your private key to your authentication agent::

    [user@localhost sage]$ ssh-add

.. NOTE::

   The ``git trac config`` command will automatically add a ``trac``
   remote git repository to your list of remotes if necessary.

If you followed the above instructions then you will have two remote
repositories set up::

    [user@localhost sage]$ git remote -v
    origin      https://github.com/sagemath/sage.git (fetch)
    origin      https://github.com/sagemath/sage.git (push)
    trac        git@trac.sagemath.org:sage.git (fetch)
    trac        git@trac.sagemath.org:sage.git (push)

The ``git@...`` part of the push url means that write access is
secured with SSH keys, which you must have set up as in
:ref:`section-trac-ssh-key`. Read-only access happens through the
fetch url and does not require SSH.

Finally, if you do not want to use the ``git trac`` subcommand at all
then you can set up the remote by hand as described in the section on
:ref:`section-git-trac`.


Trac Tickets and Git Branches
=============================

Now let's start adding code to Sage!

.. _section-git_trac-create:

Create a Ticket
---------------

Suppose you have written an algorithm for calculating the last twin prime, and
want to add it to Sage. You would first open a ticket for that::

    [user@localhost sage]$ git trac create 'Last Twin Prime'
    Remote branch: u/user/last_twin_prime
    Newly-created ticket number: 12345
    Ticket URL: https://trac.sagemath.org/12345
    Local branch: t/12345/last_twin_prime

This will create a new trac ticket titled "Last Twin Prime" with a
*remote branch* ``u/user/last_twin_prime`` attached to it. The remote
branch name is automatically derived from the ticket title; If you
don't like this then you can use the ``-b`` switch to specify it
explicitly. See ``git trac create -h`` for details. This new branch is
automatically checked out for you with the *local branch* name
``t/12345/last_twin_prime``.

.. NOTE::

    Only some trac fields are filled in automatically. See
    :ref:`section-trac-fields` for what trac fields are available and
    how we use them.



.. _section-git_trac-checkout:

Check out an Existing Ticket
----------------------------

Alternatively, you can use the `web interface to the Sage trac
development server <https://trac.sagemath.org>`_ to open a new ticket.
Just log in and click on "Create Ticket".

Or maybe somebody else already opened a ticket. Then, to get a suitable
local branch to make your edits, you would just run::

    [user@localhost sage]$ git trac checkout 12345
    Loading ticket #12345...
    Checking out Trac #13744 remote branch u/user/last_twin_prime -> local branch t/12345/last_twin_prime...

The ``git trac checkout`` command downloads an existing branch (as
specified in the "Branch:" field on the trac ticket) or creates a new
one if there is none yet. Just like the create command, you can
specify the remote branch name explicitly using the ``-b`` switch if
you want.

.. _section-git_trac-branch-names:

Note on Branch Names
--------------------

The "Branch:" field of a trac ticket (see :ref:`section-trac-fields`) indicates
the git branch containing its code. Our git server implements the following
access restrictions for **remote branch names**:

* You can read/write/create a branch named
  ``u/your_username/whatever_you_like``. Everybody else can read.

* Everybody can read/write/create a branch named ``public/whatever_you_like``.

Depending on your style of collaboration, you can use one or the
other. The ``git trac`` subcommands defaults to the former.

As a convention, the ``git trac`` subcommand uses **local branch
names** of the form ``t/12345/description``, where the number is the
trac ticket number. The script uses this number to figure out the
ticket from the local branch name. You can rename the local branches
if you want, but if they don't contain the ticket number then you will
have to specify the ticket number manually when you are uploading your
changes.

.. _section-git_trac-editing:

Making Changes
--------------

Once you have checked out a ticket, edit the appropriate files and
commit your changes to the branch as described in
:ref:`section-walkthrough-add-edit` and
:ref:`section-walkthrough-commit`.

.. _section-git_trac-push:

Uploading Changes to Trac
=========================

.. _section-git_trac-push-auto:

Automatic Push
--------------

At some point, you may wish to share your changes with the rest of us:
maybe it is ready for review, or maybe you are collaborating with
someone and want to share your changes "up until now". This is simply
done by::

    [user@localhost sage]$ git trac push
    Pushing to Trac #12345...
    Guessed remote branch: u/user/last_twin_prime

    To git@trac.sagemath.org:sage.git
     * [new branch]      HEAD -> u/user/last_twin_prime

    Changing the trac "Branch:" field...

This uploads your changes to a remote branch on the `Sage git server
<https://git.sagemath.org/sage.git>`_. The ``git trac`` command uses
the following logic to find out the remote branch name:

* By default, the remote branch name will be whatever is already on
  the trac ticket.

* If there is no remote branch yet, the branch will be called
  ``u/user/description`` (``u/user/last_twin_prime`` in the example).

* You can use the ``--branch`` option to specify the remote branch
  name explicitly, but it needs to follow the naming convention from
  :ref:`section-git_trac-branch-names` for you to have write
  permission.


.. _section-git_trac-push-with-ticket-number:

Specifying the Ticket Number
----------------------------

You can upload any local branch to an existing ticket, whether or not
you created the local branch with ``git trac``. This works exactly
like in the case where you started with a ticket, except that you have
to specify the ticket number (since there is no way to tell which
ticket you have in mind). That is::

    [user@localhost sage]$ git trac push TICKETNUM

where you have to replace ``TICKETNUM`` with the number of the trac
ticket.


.. _section-git_trac-push-finish:

Finishing It Up
---------------

It is common to go through a few iterations of commits before you
upload, and you will probably also have pushed your changes a few
times before your changes are ready for review.

Once you are happy with the changes you uploaded, they must be
reviewed by somebody else before they can be included in the next
version of Sage. To mark your ticket as ready for review, you should
set it to ``needs_review`` on the trac server. Also, add yourself as
the (or one of the) author(s) for that ticket by inserting the
following as the first line:

.. CODE-BLOCK:: text

    Authors: Your Real Name


.. _section-git_trac-pull:

Downloading Changes from Trac
=============================

If somebody else worked on a ticket, or if you just switched
computers, you'll want to get the latest version of the branch from a
ticket into your local branch. This is done with::

    [user@localhost sage]$ git trac pull

Technically, this does a *merge* (just like the standard ``git pull``)
command. See :ref:`section-git-merge` for more background information.


.. _section-git_trac-merge:

Merging
=======

As soon as you are working on a bigger project that spans multiple
tickets you will want to base your work on branches that have not been
merged into Sage yet. This is natural in collaborative development,
and in fact you are very much encouraged to split your work into
logically different parts. Ideally, each part that is useful on its
own and can be reviewed independently should be a different ticket
instead of a huge patch bomb.

For this purpose, you can incorporate branches from other tickets (or
just other local branches) into your current branch. This is called
merging, and all it does is include commits from other branches into
your current branch. In particular, this is done when a new Sage
release is made: the finished tickets are merged with the Sage master
and the result is the next Sage version. Git is smart enough to not
merge commits twice. In particular, it is possible to merge two
branches, one of which had already merged the other branch. The syntax
for merging is easy::

    [user@localhost sage]$ git merge other_branch

This creates a new "merge" commit, joining your current branch and
``other_branch``.

.. WARNING::

    You should avoid merging branches both ways. Once A merged B and B
    merged A, there is no way to distinguish commits that were
    originally made in A or B. Effectively, merging both ways combines
    the branches and makes individual review impossible.

    In practice, you should only merge when one of the following holds:

    * Either two tickets conflict, then you have to merge one into the
      other in order to resolve the merge conflict.

    * Or you definitely need a feature that has been developed as part
      of another branch.

A special case of merging is merging in the ``develop`` branch. This
brings your local branch up to date with the newest Sage version. The
above warning against unnecessary merges still applies, though. Try to
do all of your development with the Sage version that you originally
started with. The only reason for merging in the ``develop`` branch is if
you need a new feature or if your branch conflicts. See
:ref:`section-git-update-latest` for details.


.. _section-git_trac-collaborate:

Collaboration and conflict resolution
=====================================

Exchanging Branches
-------------------

It is very easy to collaborate by just going through the above steps
any number of times. For example, Alice starts a ticket and adds some
initial code::

    [alice@laptop sage]$ git trac create "A and B Ticket"
    ... EDIT EDIT ...
    [alice@laptop sage]$ git add .
    [alice@laptop sage]$ git commit
    [alice@laptop sage]$ git trac push

The trac ticket now has "Branch:" set to
``u/alice/a_and_b_ticket``. Bob downloads the branch and works some
more on it::

    [bob@home sage]$ git trac checkout TICKET_NUMBER
    ... EDIT EDIT ...
    [bob@home sage]$ git add .
    [bob@home sage]$ git commit
    [bob@home sage]$ git trac push

The trac ticket now has "Branch:" set to ``u/bob/a_and_b_ticket``,
since Bob cannot write to ``u/alice/...``. Now the two authors just
pull/push in their collaboration::

    [alice@laptop sage]$ git trac pull
    ... EDIT EDIT ...
    [alice@laptop sage]$ git add .
    [alice@laptop sage]$ git commit
    [alice@laptop sage]$ git trac push

    [bob@home sage]$ git trac pull
    ... EDIT EDIT ...
    [bob@home sage]$ git add .
    [bob@home sage]$ git commit
    [bob@home sage]$ git trac push

Alice and Bob need not alternate, they can also add further commits on
top of their own remote branch.  As long as their changes do not
conflict (edit the same lines simultaneously), this is fine.


.. _section-git_trac-conflict:

Conflict Resolution
-------------------

Merge conflicts happen if there are overlapping edits, and they are an
unavoidable consequence of distributed development. Fortunately,
resolving them is common and easy with git. As a hypothetical example,
consider the following code snippet:

.. CODE-BLOCK:: python

    def fibonacci(i):
        """
        Return the `i`-th Fibonacci number
        """
        return fibonacci(i-1) * fibonacci(i-2)

This is clearly wrong; Two developers, namely Alice and Bob, decide to
fix it. First, in a cabin in the woods far away from any internet
connection, Alice corrects the seed values:

.. CODE-BLOCK:: python

    def fibonacci(i):
       """
       Return the `i`-th Fibonacci number
       """
       if i > 1:
           return fibonacci(i-1) * fibonacci(i-2)
       return [0, 1][i]

and turns those changes into a new commit::

    [alice@laptop sage]$ git add fibonacci.py
    [alice@laptop sage]$ git commit -m 'return correct seed values'

However, not having an internet connection, she cannot immediately
send her changes to the trac server. Meanwhile, Bob changes the
multiplication to an addition since that is the correct recursion
formula:

.. CODE-BLOCK:: python

    def fibonacci(i):
        """
        Return the `i`-th Fibonacci number
        """
        return fibonacci(i-1) + fibonacci(i-2)

and immediately uploads his change::

    [bob@home sage]$ git add fibonacci.py
    [bob@home sage]$ git commit -m 'corrected recursion formula, must be + instead of *'
    [bob@home sage]$ git trac push

Eventually, Alice returns to civilization. In her mailbox, she finds a
trac notification email that Bob has uploaded further changes to their
joint project. Hence, she starts out by getting his changes into her
own local branch::

    [alice@laptop sage]$ git trac pull
    ...
    CONFLICT (content): Merge conflict in fibonacci.py
    Automatic merge failed; fix conflicts and then commit the result.

The file now looks like this:

.. skip    # doctester confuses >>> with input marker

.. CODE-BLOCK:: python

    def fibonacci(i):
        """
        Return the `i`-th Fibonacci number
        """
    <<<<<<< HEAD
        if i > 1:
            return fibonacci(i-1) * fibonacci(i-2)
        return [0, 1][i]
    =======
        return fibonacci(i-1) + fibonacci(i-2)
    >>>>>>> 41675dfaedbfb89dcff0a47e520be4aa2b6c5d1b

The conflict is shown between the conflict markers ``<<<<<<<`` and
``>>>>>>>``. The first half (up to the ``=======`` marker) is Alice's
current version, the second half is Bob's version. The 40-digit hex
number after the second conflict marker is the SHA1 hash of the most
recent common parent of both.

It is now Alice's job to resolve the conflict by reconciling their
changes, for example by editing the file. Her result is:

.. CODE-BLOCK:: python

    def fibonacci(i):
        """
        Return the `i`-th Fibonacci number
        """
        if i > 1:
            return fibonacci(i-1) + fibonacci(i-2)
        return [0, 1][i]

And then upload both her original change *and* her merge commit to trac::

    [alice@laptop sage]$ git add fibonacci.py
    [alice@laptop sage]$ git commit -m "merged Bob's changes with mine"

The resulting commit graph now has a loop::

    [alice@laptop sage]$ git log --graph --oneline
    *   6316447 merged Bob's changes with mine
    |\
    | * 41675df corrected recursion formula, must be + instead of *
    * | 14ae1d3 return correct seed values
    |/
    * 14afe53 initial commit

If Bob decides to do further work on the ticket then he will have to
pull Alice's changes. However, this time there is no conflict on his
end: git downloads both Alice's conflicting commit and her resolution.


.. _section-git_trac-review:

Reviewing
=========

For an explanation of what should be checked by the reviewer, see
:ref:`chapter-review`.

If you go to the `web interface to the Sage trac development server
<https://trac.sagemath.org>`_ then you can click on the "Branch:" field and see
the code that is added by combining all commits of the ticket. This is what
needs to be reviewed.

The ``git trac`` command gives you two commands that might be handy
(replace ``12345`` with the actual ticket number) if you do not want
to use the web interface:

* ``git trac print 12345`` displays the trac ticket directly in your
  terminal.

* ``git trac review 12345`` downloads the branch from the ticket and
  shows you what is being added, analogous to clicking on the
  "Branch:" field.

To review tickets with minimal recompiling, start by building the "develop"
branch, that is, the latest beta. Just checking out an older ticket would most
likely reset the Sage tree to an older version, so you would have to compile
older versions of packages to make it work. Instead, you can create an anonymous
("detached HEAD") merge of the ticket and the develop branch using ::

    $ git trac try 12345

This will only touch files that are really modified by the ticket. In particular,
if only Python files are changed by the ticket (which is true for most tickets)
then you just have to run ``sage -b`` to rebuild the Sage library. If files other
than Python have been changed, you must run ``make``. When you are finished
reviewing, just check out a named branch, for example ::

    $ git checkout develop

If you want to edit the ticket branch (that is, add additional commits) you cannot
use ``git trac try``. You must :ref:`section-git_trac-checkout` to get the actual ticket
branch as a starting point.
