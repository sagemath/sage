r"""
SageDev

This module provides :class:`SageDev`, the central object of the developer
scripts for sage.

AUTHORS:

- David Roe, Frej Drejhammar, Julian Rueth, Martin Raum, Nicolas M. Thiery,
  R. Andrew Ohana, Robert Bradshaw, Timo Kluck: initial version

"""
#*****************************************************************************
#       Copyright (C) 2013 David Roe <roed.math@gmail.com>
#                          Frej Drejhammar <frej.drejhammar@gmail.com>
#                          Julian Rueth <julian.rueth@fsfe.org>
#                          Martin Raum <martin@raum-brothers.eu>
#                          Nicolas M. Thiery <Nicolas.Thiery@u-psud.fr>
#                          R. Andrew Ohana <andrew.ohana@gmail.com>
#                          Robert Bradshaw <robertwb@gmail.com>
#                          Timo Kluck <tkluck@infty.nl>
#                          Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import os
import urllib, urlparse
import re

from user_interface_error import OperationCancelledError
from trac_error import TracConnectionError, TracInternalError, TracError
from git_error import GitError
from patch import MercurialPatchMixin

from sage.env import TRAC_SERVER_URI

# regular expression to check validity of git options
# http://stackoverflow.com/questions/12093748/how-do-i-check-for-valid-git-branch-names
GIT_BRANCH_REGEX = re.compile(
    r'^(?!.*/\.)(?!.*\.\.)(?!/)(?!.*//)(?!.*@\{)(?!.*\\)'
    r'[^\040\177 ~^:?*[]+(?<!\.lock)(?<!/)(?<!\.)$')

# the name of the branch which holds the vanilla clone of sage
MASTER_BRANCH = "master"
USER_BRANCH = re.compile(r"^u/([^/]+)/")

COMMIT_GUIDE=r"""


# Please type your commit message above.
#
# The first line should contain a short summary of your changes, the
# following lines should contain a more detailed description. Lines
# starting with '#' are ignored.
#
# An empty file aborts the commit.
"""

class SageDev(MercurialPatchMixin):
    r"""
    The developer interface for sage.

    This class facilitates access to git and trac.

    INPUT:

    - ``config`` -- a :class:`~sage.dev.config.Config` or ``None``
      (default: ``None``), the configuration of this object; the
      defaults uses the configuration stored in ``DOT_SAGE/devrc``.

    - ``UI`` -- a :class:`~sage.dev.user_interface.UserInterface` or ``None`` (default:
      ``None``), the default creates a
      :class:`cmd_line_interface.CmdLineInterface` from ``config['UI']``.

    - ``trac`` -- a :class:`trac_interface.TracInterface` or ``None`` (default:
      ``None``), the default creates a :class:`trac_interface.TracInterface`
      from ``config['trac']``.

    - ``git`` -- a :class:`git_interface.GitInterface` or ``None`` (default:
      ``None``), the default creates a :class:`git_interface.GitInterface` from
      ``config['git']``.

    EXAMPLES::

        sage: dev._sagedev
        SageDev()
    """
    def __init__(self, config=None, UI=None, trac=None, git=None):
        r"""
        Initialization.

        TESTS::

            sage: type(dev._sagedev)
            <class 'sage.dev.test.sagedev.DoctestSageDev'>
        """
        self.config = config
        if self.config is None:
            from config import Config
            self.config = Config()

        # create some empty config sections if they do not yet exist
        for section in ['UI','trac','git','sagedev']:
            if section not in self.config:
                self.config[section] = {}

        self._UI = UI
        if self._UI is None:
            from cmd_line_interface import CmdLineInterface
            self._UI = CmdLineInterface(self.config['UI'])

        self.trac = trac
        if self.trac is None:
            from trac_interface import TracInterface
            self.trac = TracInterface(self.config['trac'], self._UI)

        self.git = git
        if self.git is None:
            from git_interface import GitInterface
            self.git = GitInterface(self.config['git'], self._UI)

        # create some SavingDicts to store the relations between branches and tickets
        from sage.env import DOT_SAGE
        import os
        def move_legacy_saving_dict(key, old_file, new_file):
            '''
            We used to have these files in DOT_SAGE - this is not a good idea
            because a user might have multiple copies of sage which should each
            have their own set of files.

            This method moves an existing file mentioned in the config to its
            new position to support repositories created earlier.
            '''
            import sage.doctest
            if sage.doctest.DOCTEST_MODE:
                return
            import shutil
            if not os.path.exists(new_file) and os.path.exists(old_file):
                shutil.move(old_file, new_file)
                self._UI.show('The developer scripts used to store some of their data in "{0}".'
                              ' This file has now moved to "{1}". I moved "{0}" to "{1}". This might'
                              ' cause trouble if this is a fresh clone of the repository in which'
                              ' you never used the developer scripts before. In that case you'
                              ' should manually delete "{1}" now.', old_file, new_file)
            if key in self.config['sagedev']:
                del self.config['sagedev'][key]

        ticket_file = os.path.join(self.git._dot_git, 'branch_to_ticket')
        move_legacy_saving_dict('ticketfile', self.config['sagedev'].get(
            'ticketfile', os.path.join(DOT_SAGE, 'branch_to_ticket')), ticket_file)
        branch_file = os.path.join(self.git._dot_git, 'ticket_to_branch')
        move_legacy_saving_dict('branchfile', self.config['sagedev'].get(
            'branchfile', os.path.join(DOT_SAGE, 'ticket_to_branch')), branch_file)
        dependencies_file = os.path.join(self.git._dot_git, 'dependencies')
        move_legacy_saving_dict('dependenciesfile', self.config['sagedev'].get(
            'dependenciesfile', os.path.join(DOT_SAGE, 'dependencies')), dependencies_file)
        remote_branches_file = os.path.join(self.git._dot_git, 'remote_branches')
        move_legacy_saving_dict('remotebranchesfile', self.config['sagedev'].get(
            'remotebranchesfile', os.path.join(DOT_SAGE, 'remote_branches')), remote_branches_file)

        # some people dislike double underscore fields; here you can very
        # seriously screw up your setup if you put something invalid into
        # these. Ideally these fields should only be touched by single
        # underscore methods such as _set_remote_branch which do some checking
        # on the parameters
        from saving_dict import SavingDict
        self.__branch_to_ticket = SavingDict(ticket_file)
        self.__ticket_to_branch = SavingDict(branch_file, paired=self.__branch_to_ticket)
        self.__ticket_dependencies = SavingDict(dependencies_file, default=tuple)
        self.__branch_to_remote_branch = SavingDict(remote_branches_file)

    @property
    def tmp_dir(self):
        r"""
        A lazy property to provide a temporary directory

        TESTS::

            sage: import os
            sage: os.path.isdir(dev._sagedev.tmp_dir)
            True
        """
        try:
            return self._tmp_dir
        except AttributeError:
            from sage.dev.misc import tmp_dir
            self._tmp_dir = tmp_dir()
            import atexit, shutil
            atexit.register(shutil.rmtree, self._tmp_dir)
            return self._tmp_dir

    def __repr__(self):
        r"""
        Return a printable representation of this object.

        TESTS::

            sage: dev # indirect doctest
            SageDev()
        """
        return "SageDev()"

    def create_ticket(self):
        r"""
        Create a new ticket on trac.

        OUTPUT:

        Returns the number of the newly created ticket as an int.

        .. SEEALSO::

            :meth:`checkout`, :meth:`pull`, :meth:`edit_ticket`

        TESTS:

        Set up a single user environment::

            sage: from sage.dev.test.sagedev import single_user_setup
            sage: dev, config, UI, server = single_user_setup()
            sage: dev._wrap("_dependencies_for_ticket")

        Create some tickets::

            sage: UI.append("Summary: ticket1\ndescription")
            sage: dev.create_ticket()
            Created ticket #1 at https://trac.sagemath.org/1.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=1" to create a new local branch)
            1

            sage: UI.append("Summary: ticket2\ndescription")
            sage: dev.create_ticket()
            Created ticket #2 at https://trac.sagemath.org/2.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=2" to create a new local branch)
            2

        This fails if the internet connection is broken::

            sage: dev.trac._connected = False
            sage: UI.append("Summary: ticket7\ndescription")
            sage: dev.create_ticket()
            A network error ocurred, ticket creation aborted.
            Your command failed because no connection to trac could be established.
            sage: dev.trac._connected = True
        """
        try:
            ticket = self.trac.create_ticket_interactive()
        except OperationCancelledError:
            self._UI.debug("Ticket creation aborted.")
            raise
        except TracConnectionError as e:
            self._UI.error("A network error ocurred, ticket creation aborted.")
            raise
        ticket_url = urlparse.urljoin(self.trac._config.get('server', TRAC_SERVER_URI), str(ticket))
        self._UI.show("Created ticket #{0} at {1}.".format(ticket, ticket_url))
        self._UI.info(['',
                       '(use "{0}" to create a new local branch)'
                       .format(self._format_command("checkout", ticket=ticket))])
        return ticket

    def checkout(self, ticket=None, branch=None, base=''):
        r"""
        Checkout another branch.

        If ``ticket`` is specified, and ``branch`` is an existing local branch,
        then ``ticket`` will be associated to it, and ``branch`` will be
        checked out into the working directory.
        Otherwise, if there is no local branch for ``ticket``, the branch
        specified on trac will be pulled to ``branch`` unless ``base`` is
        set to something other than the empty string ``''``. If the trac ticket
        does not specify a branch yet or if ``base`` is not the empty string,
        then a new one will be created from ``base`` (per default, the master
        branch).

        If ``ticket`` is not specified, then checkout the local branch
        ``branch`` into the working directory.

        INPUT:

        - ``ticket`` -- a string or an integer identifying a ticket or ``None``
          (default: ``None``)

        - ``branch`` -- a string, the name of a local branch; if ``ticket`` is
          specified, then this defaults to ticket/``ticket``.

        - ``base`` -- a string or ``None``, a branch on which to base a new
          branch if one is going to be created (default: the empty string
          ``''`` to create the new branch from the master branch), or a ticket;
          if ``base`` is set to ``None``, then the current ticket is used. If
          ``base`` is a ticket, then the corresponding dependency will be
          added. Must be ``''`` if ``ticket`` is not specified.

        .. SEEALSO::

            :meth:`pull`, :meth:`create_ticket`, :meth:`vanilla`

        TESTS:

        Set up a single user for doctesting::

            sage: from sage.dev.test.sagedev import single_user_setup
            sage: dev, config, UI, server = single_user_setup()

        Create a few branches::

            sage: dev.git.silent.branch("branch1")
            sage: dev.git.silent.branch("branch2")

        Checking out a branch::

            sage: dev.checkout(branch="branch1")
            On local branch "branch1" without associated ticket.
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.
            sage: dev.git.current_branch()
            'branch1'

        Create a ticket and checkout a branch for it::

            sage: UI.append("Summary: summary\ndescription")
            sage: dev.create_ticket()
            Created ticket #1 at https://trac.sagemath.org/1.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=1" to create a new local branch)
            1
            sage: dev.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.
            sage: dev.git.current_branch()
            'ticket/1'
        """
        if ticket is not None:
            self.checkout_ticket(ticket=ticket, branch=branch, base=base)
        elif branch is not None:
            if base != '':
                raise SageDevValueError("base must not be specified if no ticket is specified.")
            self.checkout_branch(branch=branch)
        else:
            raise SageDevValueError("at least one of ticket or branch must be specified.")

        ticket = self._current_ticket()
        branch = self.git.current_branch()
        if ticket:
            self._UI.show(['On ticket #{0} with associated local branch "{1}".'], ticket, branch)
        else:
            self._UI.show(['On local branch "{0}" without associated ticket.'], branch)
        self._UI.info(['',
                       'Use "{0}" to include another ticket/branch.',
                       'Use "{1}" to save changes into a new commit.'],
                      self._format_command("merge"),
                      self._format_command("commit"))


    def checkout_ticket(self, ticket, branch=None, base=''):
        r"""
        Checkout the branch associated to ``ticket``.

        If ``branch`` is an existing local branch, then ``ticket`` will be
        associated to it, and ``branch`` will be checked out into the working directory.

        Otherwise, if there is no local branch for ``ticket``, the branch
        specified on trac will be pulled to ``branch`` unless ``base`` is
        set to something other than the empty string ``''``. If the trac ticket
        does not specify a branch yet or if ``base`` is not the empty string,
        then a new one will be created from ``base`` (per default, the master
        branch).

        INPUT:

        - ``ticket`` -- a string or an integer identifying a ticket

        - ``branch`` -- a string, the name of the local branch that stores
          changes for ``ticket`` (default: ticket/``ticket``)

        - ``base`` -- a string or ``None``, a branch on which to base a new
          branch if one is going to be created (default: the empty string
          ``''`` to create the new branch from the master branch), or a ticket;
          if ``base`` is set to ``None``, then the current ticket is used. If
          ``base`` is a ticket, then the corresponding dependency will be
          added.

        .. SEEALSO::

            :meth:`pull`, :meth:`create_ticket`, :meth:`vanilla`

        TESTS:

        Create a doctest setup with two users::

            sage: from sage.dev.test.sagedev import two_user_setup
            sage: alice, config_alice, bob, config_bob, server = two_user_setup()

        Alice tries to checkout ticket #1 which does not exist yet::

            sage: alice._chdir()
            sage: alice.checkout(ticket=1)
            Ticket name "1" is not valid or ticket does not exist on trac.

        Bob creates that ticket::

            sage: bob._chdir()
            sage: bob._UI.append("Summary: summary1\ndescription")
            sage: bob.create_ticket()
            Created ticket #1 at https://trac.sagemath.org/1.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=1" to create a new local branch)
            1
            sage: bob.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.

        Now alice can check it out, even though there is no branch on the
        ticket description::

            sage: alice._chdir()
            sage: alice.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.

        If Bob commits something to the ticket, a ``checkout`` by Alice
        does not take his changes into account::

            sage: bob._chdir()
            sage: bob.git.super_silent.commit(allow_empty=True,message="empty commit")
            sage: bob._UI.append("y")
            sage: bob.push()
            The branch "u/bob/ticket/1" does not exist on the remote server.
            Create new remote branch? [Yes/no] y

            sage: alice._chdir()
            sage: alice.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.
            sage: alice.git.echo.log('--pretty=%s')
            initial commit

        If Alice had not checked that ticket out before, she would of course
        see Bob's changes (this also checks that we can handle a corrupt ticket
        database and a detached HEAD)::

            sage: alice.git.super_silent.checkout('HEAD', detach=True)
            sage: alice.git.super_silent.branch('-d','ticket/1')
            sage: alice.checkout(ticket=1) # ticket #1 refers to the non-existant branch 'ticket/1'
            Ticket #1 refers to the non-existant local branch "ticket/1". If you have not
            manually interacted with git, then this is a bug in sagedev. Removing the
            association from ticket #1 to branch "ticket/1".
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.
            sage: alice.git.current_branch()
            'ticket/1'
            sage: alice.git.echo.log('--pretty=%s')
            empty commit
            initial commit

        Checking out a ticket with untracked files::

            sage: alice._UI.append("Summary: summary2\ndescription")
            sage: alice.create_ticket()
            Created ticket #2 at https://trac.sagemath.org/2.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=2" to create a new local branch)
            2
            sage: alice.checkout(ticket=2)
            On ticket #2 with associated local branch "ticket/2".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.
            sage: alice.git.echo.log('--pretty=%s')
            initial commit
            sage: open("untracked","w").close()
            sage: alice.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.
            sage: alice.git.echo.log('--pretty=%s')
            empty commit
            initial commit

        Checking out a ticket with untracked files which make a checkout
        impossible::

            sage: alice.git.super_silent.add("untracked")
            sage: alice.git.super_silent.commit(message="added untracked")
            sage: alice.checkout(ticket=2)
            On ticket #2 with associated local branch "ticket/2".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.
            sage: open("untracked","w").close()
            sage: alice.checkout(ticket=1)
            GitError: git exited with a non-zero exit code (1).
            This happened while executing "git -c user.email=doc@test.test -c
            user.name=alice checkout ticket/1".
            git printed nothing to STDOUT.
            git printed the following to STDERR:
            error: The following untracked working tree files would be overwritten by checkout:
                untracked
            Please move or remove them before you can switch branches.
            Aborting

        Checking out a ticket with uncommited changes::

            sage: open("tracked", "w").close()
            sage: alice.git.super_silent.add("tracked")
            sage: alice._UI.append('d')
            sage: alice.checkout(ticket=2)
            The following files in your working directory contain uncommitted changes:
            <BLANKLINE>
                 tracked
            <BLANKLINE>
            Discard changes? [discard/Keep/stash] d
            On ticket #2 with associated local branch "ticket/2".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.

        Now follow some single user tests to check that the parameters are interpreted correctly::

            sage: from sage.dev.test.sagedev import single_user_setup
            sage: dev, config, UI, server = single_user_setup()
            sage: dev._wrap("_dependencies_for_ticket")

        First, create some tickets::

            sage: UI.append("Summary: ticket1\ndescription")
            sage: dev.create_ticket()
            Created ticket #1 at https://trac.sagemath.org/1.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=1" to create a new local branch)
            1
            sage: dev.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.
            sage: UI.append("Summary: ticket2\ndescription")
            sage: dev.create_ticket()
            Created ticket #2 at https://trac.sagemath.org/2.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=2" to create a new local branch)
            2
            sage: dev.checkout(ticket=2)
            On ticket #2 with associated local branch "ticket/2".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.
            sage: dev.git.silent.commit(allow_empty=True, message="second commit")
            sage: dev.git.commit_for_branch('ticket/2') != dev.git.commit_for_branch('ticket/1')
            True

        Check that ``base`` works::

            sage: UI.append("Summary: ticket3\ndescription")
            sage: dev.create_ticket()
            Created ticket #3 at https://trac.sagemath.org/3.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=3" to create a new local branch)
            3
            sage: dev.checkout(ticket=3, base=2)
            On ticket #3 with associated local branch "ticket/3".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.
            sage: dev.git.commit_for_branch('ticket/3') == dev.git.commit_for_branch('ticket/2')
            True
            sage: dev._dependencies_for_ticket(3)
            (2,)
            sage: UI.append("Summary: ticket4\ndescription")
            sage: dev.create_ticket()
            Created ticket #4 at https://trac.sagemath.org/4.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=4" to create a new local branch)
            4
            sage: dev.checkout(ticket=4, base='ticket/2')
            On ticket #4 with associated local branch "ticket/4".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.
            sage: dev.git.commit_for_branch('ticket/4') == dev.git.commit_for_branch('ticket/2')
            True
            sage: dev._dependencies_for_ticket(4)
            ()

        In this example ``base`` does not exist::

            sage: UI.append("Summary: ticket5\ndescription")
            sage: dev.create_ticket()
            Created ticket #5 at https://trac.sagemath.org/5.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=5" to create a new local branch)
            5
            sage: dev.checkout(ticket=5, base=1000)
            Ticket name "1000" is not valid or ticket does not exist on trac.

        In this example ``base`` does not exist locally::

            sage: UI.append("Summary: ticket6\ndescription")
            sage: dev.create_ticket()
            Created ticket #6 at https://trac.sagemath.org/6.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=6" to create a new local branch)
            6
            sage: dev.checkout(ticket=6, base=5)
            Branch field is not set for ticket #5 on trac.

        Creating a ticket when in detached HEAD state::

            sage: dev.git.super_silent.checkout('HEAD', detach=True)
            sage: UI.append("Summary: ticket detached\ndescription")
            sage: dev.create_ticket()
            Created ticket #7 at https://trac.sagemath.org/7.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=7" to create a new local branch)
            7
            sage: dev.checkout(ticket=7)
            On ticket #7 with associated local branch "ticket/7".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.
            sage: dev.git.current_branch()
            'ticket/7'

        Creating a ticket when in the middle of a merge::

            sage: dev.git.super_silent.checkout('-b','merge_branch')
            sage: with open('merge', 'w') as f: f.write("version 0")
            sage: dev.git.silent.add('merge')
            sage: dev.git.silent.commit('-m','some change')
            sage: dev.git.super_silent.checkout('ticket/7')
            sage: with open('merge', 'w') as f: f.write("version 1")
            sage: dev.git.silent.add('merge')
            sage: dev.git.silent.commit('-m','conflicting change')
            sage: from sage.dev.git_error import GitError
            sage: try:
            ....:     dev.git.silent.merge('merge_branch')
            ....: except GitError: pass
            sage: UI.append("Summary: ticket merge\ndescription")
            sage: dev.create_ticket()
            Created ticket #8 at https://trac.sagemath.org/8.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=8" to create a new local branch)
            8
            sage: UI.append("cancel")
            sage: dev.checkout(ticket=8)
            Repository is in an unclean state (merge). Resetting the state will discard any
            uncommited changes.
            Reset repository? [reset/Cancel] cancel
            Aborting checkout of branch "ticket/8".
            <BLANKLINE>
            #  (use "sage --dev commit" to save changes in a new commit)
            sage: dev.git.reset_to_clean_state()

        Creating a ticket with uncommitted changes::

            sage: open('tracked', 'w').close()
            sage: dev.git.silent.add('tracked')
            sage: UI.append("Summary: ticket merge\ndescription")
            sage: dev.create_ticket()
            Created ticket #9 at https://trac.sagemath.org/9.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=9" to create a new local branch)
            9

        The new branch is based on master which is not the same commit
        as the current branch ``ticket/7``, so it is not a valid
        option to ``'keep'`` changes::

            sage: UI.append("cancel")
            sage: dev.checkout(ticket=9)
            The following files in your working directory contain uncommitted changes:
            <BLANKLINE>
                 tracked
            <BLANKLINE>
            Discard changes? [discard/Cancel/stash] cancel
            Aborting checkout of branch "ticket/9".
            <BLANKLINE>
            #  (use "sage --dev commit" to save changes in a new commit)

        Finally, in this case we can keep changes because the base is
        the same commit as the current branch::

            sage: UI.append("Summary: ticket merge\ndescription")
            sage: dev.create_ticket()
            Created ticket #10 at https://trac.sagemath.org/10.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=10" to create a new local branch)
            10
            sage: UI.append("keep")
            sage: dev.checkout(ticket=10, base='ticket/7')
            The following files in your working directory contain uncommitted changes:
            <BLANKLINE>
                 tracked
            <BLANKLINE>
            Discard changes? [discard/Keep/stash] keep
            On ticket #10 with associated local branch "ticket/10".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.
        """
        self._check_ticket_name(ticket, exists=True)
        ticket = self._ticket_from_ticket_name(ticket)

        # if branch points to an existing branch make it the ticket's branch and check it out
        if branch is not None and self._is_local_branch_name(branch, exists=True):
            if base != MASTER_BRANCH:
                raise SageDevValueError("base must not be specified if branch is an existing branch")
            if branch == MASTER_BRANCH:
                raise SageDevValueError("branch must not be the master branch")

            self._set_local_branch_for_ticket(ticket, branch)
            self._UI.debug('The branch for ticket #{0} is now "{1}".', ticket, branch)
            self._UI.debug('Now checking out branch "{0}".', branch)
            self.checkout_branch(branch)
            return

        # if there is a branch for ticket locally, check it out
        if branch is None:
            if self._has_local_branch_for_ticket(ticket):
                branch = self._local_branch_for_ticket(ticket)
                self._UI.debug('Checking out branch "{0}".', branch)
                self.checkout_branch(branch)
                return
            else:
                branch = self._new_local_branch_for_ticket(ticket)

        # branch does not exist, so we have to create a new branch for ticket
        # depending on the value of base, this will either be base or a copy of
        # the branch mentioned on trac if any
        dependencies = self.trac.dependencies(ticket)
        if base is None:
            base = self._current_ticket()
        if base is None:
            raise SageDevValueError('currently on no ticket, "base" must not be None')
        if self._is_ticket_name(base):
            base = self._ticket_from_ticket_name(base)
            dependencies = [base] # we create a new branch for this ticket - ignore the dependencies which are on trac
            base = self._local_branch_for_ticket(base, pull_if_not_found=True)

        remote_branch = self.trac._branch_for_ticket(ticket)
        try:
            if base == '':
                base = MASTER_BRANCH
                if remote_branch is None: # branch field is not set on ticket
                    # create a new branch off master
                    self._UI.debug('The branch field on ticket #{0} is not set. Creating a new branch'
                                   ' "{1}" off the master branch "{2}".', ticket, branch, MASTER_BRANCH)
                    self.git.silent.branch(branch, MASTER_BRANCH)
                else:
                    # pull the branch mentioned on trac
                    if not self._is_remote_branch_name(remote_branch, exists=True):
                        self._UI.error('The branch field on ticket #{0} is set to the non-existent "{1}".'
                                       ' Please set the field on trac to a field value.',
                                       ticket, remote_branch)
                        self._UI.info(['', '(use "{0}" to edit the ticket description)'],
                                       self._format_command("edit-ticket", ticket=ticket))
                        raise OperationCancelledError("remote branch does not exist")

                    self.git.super_silent.fetch(self.git._repository_anonymous, remote_branch)
                    self.git.super_silent.branch(branch, 'FETCH_HEAD')
            else:
                self._check_local_branch_name(base, exists=True)
                if remote_branch is not None:
                    self._UI.show('About to create a new branch for #{0} based on "{1}". However, the trac'
                                  ' ticket for #{0} already refers to the branch "{2}". The new branch will'
                                  ' not contain any work that has already been done on "{2}".',
                                  ticket, base, remote_branch)
                    if not self._UI.confirm('Create fresh branch?', default=False):
                        command = ""
                        if self._has_local_branch_for_ticket(ticket):
                            command += self._format_command("abandon", self._local_branch_for_ticket(ticket)) + "; "
                        command += self._format_command("checkout", ticket=ticket)
                        self._UI.info(['', 'Use "{1}" to work on a local copy of the existing remote branch "{0}".'],
                                      remote_branch, command)
                        raise OperationCancelledError("user requested")

                self._UI.debug('Creating a new branch for #{0} based on "{1}".', ticket, base)
                self.git.silent.branch(branch, base)
        except:
            if self._is_local_branch_name(branch, exists=True):
                self._UI.debug('Deleting local branch "{0}".', branch)
                self.git.super_silent.branch(branch, D=True)
            raise

        self._set_local_branch_for_ticket(ticket, branch)
        if dependencies:
            self._UI.debug("Locally recording dependency on {0} for #{1}.",
                           ", ".join(["#"+str(dep) for dep in dependencies]), ticket)
            self._set_dependencies_for_ticket(ticket, dependencies)
        self._set_remote_branch_for_branch(branch, self._remote_branch_for_ticket(ticket))
        self._UI.debug('Checking out to newly created branch "{0}".'.format(branch))
        self.checkout_branch(branch)

    def checkout_branch(self, branch, helpful=True):
        r"""
        Checkout to the local branch ``branch``.

        INPUT:

        - ``branch`` -- a string, the name of a local branch

        TESTS:

        Set up a single user for doctesting::

            sage: from sage.dev.test.sagedev import single_user_setup
            sage: dev, config, UI, server = single_user_setup()

        Create a few branches::

            sage: dev.git.silent.branch("branch1")
            sage: dev.git.silent.branch("branch2")

        Checking out a branch::

            sage: dev.checkout(branch="branch1")
            On local branch "branch1" without associated ticket.
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.
            sage: dev.git.current_branch()
            'branch1'

        The branch must exist::

            sage: dev.checkout(branch="branch3")
            Branch "branch3" does not exist locally.
            <BLANKLINE>
            #  (use "sage --dev tickets" to list local branches)

        Checking out branches with untracked files::

            sage: open("untracked", "w").close()
            sage: dev.checkout(branch="branch2")
            On local branch "branch2" without associated ticket.
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.

        Checking out a branch with uncommitted changes::

            sage: open("tracked", "w").close()
            sage: dev.git.silent.add("tracked")
            sage: dev.git.silent.commit(message="added tracked")
            sage: with open("tracked", "w") as f: f.write("foo")
            sage: UI.append("cancel")
            sage: dev.checkout(branch="branch1")
            The following files in your working directory contain uncommitted changes:
            <BLANKLINE>
                 tracked
            <BLANKLINE>
            Discard changes? [discard/Cancel/stash] cancel
            Aborting checkout of branch "branch1".
            <BLANKLINE>
            #  (use "sage --dev commit" to save changes in a new commit)

        We can stash uncommitted changes::

            sage: UI.append("s")
            sage: dev.checkout(branch="branch1")
            The following files in your working directory contain uncommitted changes:
            <BLANKLINE>
                 tracked
            <BLANKLINE>
            Discard changes? [discard/Cancel/stash] s
            Your changes have been moved to the git stash stack. To re-apply your changes
            later use "git stash apply".
            On local branch "branch1" without associated ticket.
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.

        And retrieve the stashed changes later::

            sage: dev.checkout(branch='branch2')
            On local branch "branch2" without associated ticket.
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.
            sage: dev.git.echo.stash('apply')
            # On branch branch2
            # Changes not staged for commit:
            #   (use "git add <file>..." to update what will be committed)
            #   (use "git checkout -- <file>..." to discard changes in working directory)
            #
            #   modified:   tracked
            #
            # Untracked files:
            #   (use "git add <file>..." to include in what will be committed)
            #
            #   untracked
            no changes added to commit (use "git add" and/or "git commit -a")

        Or we can just discard the changes::

            sage: UI.append("discard")
            sage: dev.checkout(branch="branch1")
            The following files in your working directory contain uncommitted changes:
            <BLANKLINE>
                 tracked
            <BLANKLINE>
            Discard changes? [discard/Cancel/stash] discard
            On local branch "branch1" without associated ticket.
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.

        Checking out a branch when in the middle of a merge::

            sage: dev.git.super_silent.checkout('-b','merge_branch')
            sage: with open('merge', 'w') as f: f.write("version 0")
            sage: dev.git.silent.add('merge')
            sage: dev.git.silent.commit('-m','some change')
            sage: dev.git.super_silent.checkout('branch1')
            sage: with open('merge', 'w') as f: f.write("version 1")
            sage: dev.git.silent.add('merge')
            sage: dev.git.silent.commit('-m','conflicting change')
            sage: from sage.dev.git_error import GitError
            sage: try:
            ....:     dev.git.silent.merge('merge_branch')
            ....: except GitError: pass
            sage: UI.append('r')
            sage: dev.checkout(branch='merge_branch')
            Repository is in an unclean state (merge). Resetting the state will discard any
            uncommited changes.
            Reset repository? [reset/Cancel] r
            On local branch "merge_branch" without associated ticket.
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.

        Checking out a branch when in a detached HEAD::

            sage: dev.git.super_silent.checkout('branch2', detach=True)
            sage: dev.checkout(branch='branch1')
            On local branch "branch1" without associated ticket.
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.

        With uncommitted changes::

            sage: dev.git.super_silent.checkout('branch2', detach=True)
            sage: with open('tracked', 'w') as f: f.write("boo")
            sage: UI.append("discard")
            sage: dev.checkout(branch='branch1')
            The following files in your working directory contain uncommitted changes:
            <BLANKLINE>
                 tracked
            <BLANKLINE>
            Discard changes? [discard/Cancel/stash] discard
            On local branch "branch1" without associated ticket.
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.

        Checking out a branch with untracked files that would be overwritten by
        the checkout::

            sage: with open('tracked', 'w') as f: f.write("boo")
            sage: dev.checkout(branch='branch2')
            GitError: git exited with a non-zero exit code (1).
            This happened while executing "git -c user.email=doc@test.test -c
            user.name=doctest checkout branch2".
            git printed nothing to STDOUT.
            git printed the following to STDERR:
            error: The following untracked working tree files would be overwritten
            by checkout:
                tracked
            Please move or remove them before you can switch branches.
            Aborting
        """
        self._check_local_branch_name(branch, exists=True)

        try:
            self.reset_to_clean_state(helpful=False)
        except OperationCancelledError:
            if helpful:
                self._UI.show('Aborting checkout of branch "{0}".', branch)
                self._UI.info(['', '(use "{0}" to save changes in a new commit)'],
                              self._format_command("commit"))
            raise

        current_commit = self.git.commit_for_ref('HEAD')
        target_commit = self.git.commit_for_ref(branch)
        try:
            self.clean(error_unless_clean=(current_commit != target_commit))
        except OperationCancelledError:
            if helpful:
                self._UI.show('Aborting checkout of branch "{0}".', branch)
                self._UI.info(['', '(use "{0}" to save changes in a new commit)'],
                              self._format_command("commit"))
            raise

        try:
            # this leaves locally modified files intact (we only allow this to happen
            # if current_commit == target_commit
            self.git.super_silent.checkout(branch)
        except GitError as e:
            # the error message should be self explanatory
            raise

    def pull(self, ticket_or_remote_branch=None):
        r"""
        Pull ``ticket_or_remote_branch`` to ``branch``.

        INPUT:

        - ``ticket_or_remote_branch`` -- a string or an integer or ``None`` (default:
          ``None``), a ticket or a remote branch name; setting this to ``None``
          has the same effect as setting it to the :meth:`current_ticket`.

        TESTS:

        Create a doctest setup with two users::

            sage: from sage.dev.test.sagedev import two_user_setup
            sage: alice, config_alice, bob, config_bob, server = two_user_setup()

        Alice creates ticket 1::

            sage: alice._chdir()
            sage: alice._UI.append("Summary: summary1\ndescription")
            sage: alice.create_ticket()
            Created ticket #1 at https://trac.sagemath.org/1.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=1" to create a new local branch)
            1
            sage: alice.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.

        Bob attempts to pull for the ticket but fails because there is no
        branch for the ticket yet::

            sage: bob._chdir()
            sage: bob.pull(1)
            Branch field is not set for ticket #1 on trac.

        So, Bob starts to work on the ticket on a new branch::

            sage: bob.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.

        Alice pushes a commit::

            sage: alice._chdir()
            sage: alice.git.super_silent.commit(allow_empty=True, message="alice: empty commit")
            sage: alice._UI.append("y")
            sage: alice.push()
            The branch "u/alice/ticket/1" does not exist on the remote server.
            Create new remote branch? [Yes/no] y

        Bob pulls the changes for ticket 1::

            sage: bob._chdir()
            sage: bob.pull()
            Merging the remote branch "u/alice/ticket/1" into the local branch "ticket/1".
            Automatic merge successful.
            <BLANKLINE>
            #  (use "sage --dev commit" to commit your merge)
            sage: bob.git.echo.log('--pretty=%s')
            alice: empty commit
            initial commit

        Bob commits a change::

            sage: open("bobs_file","w").close()
            sage: bob.git.silent.add("bobs_file")
            sage: bob.git.super_silent.commit(message="bob: added bobs_file")
            sage: bob._UI.append("y")
            sage: bob._UI.append("y")
            sage: bob.push()
            The branch "u/bob/ticket/1" does not exist on the remote server.
            Create new remote branch? [Yes/no] y
            The branch field of ticket #1 needs to be updated from its current value
            "u/alice/ticket/1" to "u/bob/ticket/1"
            Change the "Branch:" field? [Yes/no] y

        Alice commits non-conflicting changes::

            sage: alice._chdir()
            sage: with open("alices_file","w") as f: f.write("1")
            sage: alice.git.silent.add("alices_file")
            sage: alice.git.super_silent.commit(message="alice: added alices_file")

        Alice can now pull the changes by Bob without the need to merge
        manually::

            sage: alice.pull()
            Merging the remote branch "u/bob/ticket/1" into the local branch "ticket/1".
            Automatic merge successful.
            <BLANKLINE>
            #  (use "sage --dev commit" to commit your merge)
            sage: alice.git.echo.log('--pretty=%s')
            Merge branch 'u/bob/ticket/1' of ... into ticket/1
            alice: added alices_file
            bob: added bobs_file
            alice: empty commit
            initial commit

        Now, Bob commits some conflicting changes::

            sage: bob._chdir()
            sage: with open("alices_file","w") as f: f.write("2")
            sage: bob.git.silent.add("alices_file")
            sage: bob.git.super_silent.commit(message="bob: added alices_file")
            sage: bob._UI.append('y')
            sage: bob.push()
            Local commits that are not on the remote branch "u/bob/ticket/1":
            <BLANKLINE>
                ...: bob: added alices_file
            <BLANKLINE>
            Push to remote branch? [Yes/no] y

        Now, the pull fails; one would have to use :meth:`merge`::

            sage: alice._chdir()
            sage: alice._UI.append("abort")
            sage: alice.pull()
            Merging the remote branch "u/bob/ticket/1" into the local branch "ticket/1".
            Automatic merge failed, there are conflicting commits.
            <BLANKLINE>
            Auto-merging alices_file
            CONFLICT (add/add): Merge conflict in alices_file
            <BLANKLINE>
            Please edit the affected files to resolve the conflicts. When you are finished,
            your resolution will be commited.
            Finished? [ok/Abort] abort

        Undo the latest commit by alice, so we can pull again::

            sage: alice.git.super_silent.reset('HEAD~~', hard=True)
            sage: alice.pull()
            Merging the remote branch "u/bob/ticket/1" into the local branch "ticket/1".
            Automatic merge successful.
            <BLANKLINE>
            #  (use "sage --dev commit" to commit your merge)
            sage: alice.git.echo.log('--pretty=%s')
            bob: added alices_file
            bob: added bobs_file
            alice: empty commit
            initial commit

        Now, Alice creates an untracked file which makes a trivial merge
        impossible::

            sage: alice._chdir()
            sage: open("bobs_other_file","w").close()

            sage: bob._chdir()
            sage: open("bobs_other_file","w").close()
            sage: bob.git.super_silent.add("bobs_other_file")
            sage: bob.git.super_silent.commit(message="bob: added bobs_other_file")
            sage: bob._UI.append('y')
            sage: bob.push()
            Local commits that are not on the remote branch "u/bob/ticket/1":
            <BLANKLINE>
                ...: bob: added bobs_other_file
            <BLANKLINE>
            Push to remote branch? [Yes/no] y

            sage: alice._chdir()
            sage: alice._UI.append("abort")
            sage: alice.pull()
            Merging the remote branch "u/bob/ticket/1" into the local branch "ticket/1".
            Automatic merge failed, there are conflicting commits.
            <BLANKLINE>
            Updating ...
            error: The following untracked working tree files would be overwritten by merge:
                    bobs_other_file
            Please move or remove them before you can merge.
            <BLANKLINE>
            Please edit the affected files to resolve the conflicts. When you are finished,
            your resolution will be commited.
            Finished? [ok/Abort] abort
        """
        if ticket_or_remote_branch is None:
            ticket_or_remote_branch = self._current_ticket()

        if self._is_ticket_name(ticket_or_remote_branch):
            ticket = self._ticket_from_ticket_name(ticket_or_remote_branch)
            self._check_ticket_name(ticket, exists=True)

            remote_branch = self.trac._branch_for_ticket(ticket)
            if remote_branch is None:
                raise SageDevValueError("Branch field is not set for ticket #{0} on trac.".format(ticket))
        else:
            remote_branch = ticket_or_remote_branch

        self.merge(remote_branch, pull=True)

    def commit(self, message=None, interactive=False):
        r"""
        Create a commit from the pending changes on the current branch.

        This is most akin to mercurial's commit command, not git's,
        since we do not require users to add files.

        INPUT:

        - ``message`` -- the message of the commit (default: ``None``), if
          ``None``, prompt for a message.

        - ``interactive`` -- if set, interactively select which part of the
          changes should be part of the commit

        .. SEEALSO::

            - :meth:`push` -- Push changes to the remote server.  This
              is the next step once you've committed some changes.

            - :meth:`diff` -- Show changes that will be committed.

        TESTS:

        Set up a single user for doctesting::

            sage: from sage.dev.test.sagedev import single_user_setup
            sage: dev, config, UI, server = single_user_setup()

        Commit an untracked file::

            sage: dev.git.super_silent.checkout('-b', 'branch1')
            sage: open("tracked","w").close()
            sage: dev._UI.extend(["y", "added tracked", "y", "y"])
            sage: dev.commit()
            The following files in your working directory are not tracked by git:
            <BLANKLINE>
                tracked
            <BLANKLINE>
            Start tracking any of these files? [yes/No] y
            Start tracking "tracked"? [yes/No] y
            Commit your changes to branch "branch1"? [Yes/no] y
            <BLANKLINE>
            #  Use "sage --dev push" to push your commits to the trac server once you are done.

        Commit a tracked file::

            sage: with open("tracked", "w") as F: F.write("foo")
            sage: dev._UI.append('y')
            sage: dev.commit(message='modified tracked')
            Commit your changes to branch "branch1"? [Yes/no] y
            <BLANKLINE>
            #  Use "sage --dev push" to push your commits to the trac server once you are done.
        """
        from git_error import DetachedHeadError
        try:
            branch = self.git.current_branch()
        except DetachedHeadError:
            self._UI.error("Cannot commit changes when not on any branch.")
            self._UI.info(['',
                           '(use "{0}" to checkout a branch)'
                           .format(self._format_command("checkout"))])
            raise OperationCancelledError("cannot proceed in detached HEAD mode")

        # make sure the index is clean
        self.git.super_silent.reset()

        try:
            self._UI.debug('Committing pending changes to branch "{0}".'.format(branch))

            try:
                untracked_files = self.git.untracked_files()
                if untracked_files:
                    self._UI.show(['The following files in your working directory are not tracked by git:', ''] +
                                  ['    ' + f for f in untracked_files ] +
                                  [''])
                    if self._UI.confirm('Start tracking any of these files?', default=False):
                        for file in untracked_files:
                            if self._UI.confirm('Start tracking "{0}"?'.format(file), default=False):
                                self.git.add(file)

                if interactive:
                    self.git.echo.add(patch=True)
                else:
                    self.git.echo.add(self.git._src, update=True)

                if message is None:
                    from sage.dev.misc import tmp_filename
                    commit_message = tmp_filename()
                    with open(commit_message, 'w') as f:
                        f.write(COMMIT_GUIDE)
                    self._UI.edit(commit_message)
                    message = "\n".join([line for line in open(commit_message).read().splitlines()
                                         if not line.startswith("#")]).strip()
                if not message:
                    raise OperationCancelledError("empty commit message")

                if not self._UI.confirm('Commit your changes to branch "{0}"?'.format(branch), default=True):
                    self._UI.info(['', 'Run "{0}" first if you want to commit to a different branch or ticket.'],
                                  self._format_command("checkout"))
                    raise OperationCancelledError("user does not want to create a commit")
                self.git.commit(message=message)
                self._UI.debug("A commit has been created.")
                self._UI.info(['', 'Use "{0}" to push your commits to the trac server once you are done.'],
                              self._format_command("push"))
            except OperationCancelledError:
                self._UI.debug("Not creating a commit.")
                raise
            except:
                self._UI.error("No commit has been created.")
                raise

        finally:
            # do not leave a non-clean index behind
            self.git.super_silent.reset()

    def set_remote(self, branch_or_ticket, remote_branch):
        r"""
        Set the remote branch to push to for ``branch_or_ticket`` to
        ``remote_branch``.

        INPUT:

        - ``branch_or_ticket`` -- a string, the name of a local branch, or a
          string or an integer identifying a ticket or ``None``; if ``None``,
          the current branch is used.

        - ``remote_branch`` -- a string, the name of a remote branch (this
          branch may not exist yet)

        .. SEEALSO::

            - :meth:`push` -- To push changes after setting the remote
              branch

        TESTS:

        Set up a single user for doctesting::

            sage: from sage.dev.test.sagedev import single_user_setup
            sage: dev, config, UI, server = single_user_setup()
            sage: dev._wrap("_remote_branch_for_ticket")

        Create a new branch::

            sage: UI.append("Summary: ticket1\ndescription")
            sage: dev.create_ticket()
            Created ticket #1 at https://trac.sagemath.org/1.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=1" to create a new local branch)
            1
            sage: dev.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.

        Modify the remote branch for this ticket's branch::

            sage: dev._remote_branch_for_ticket(1)
            'u/doctest/ticket/1'
            sage: dev.set_remote('ticket/1', 'u/doctest/foo')
            sage: dev._remote_branch_for_ticket(1)
            'u/doctest/foo'
            sage: dev.set_remote('ticket/1', 'foo')
            sage: dev._remote_branch_for_ticket(1)
            'foo'
            sage: dev.set_remote('#1', 'u/doctest/foo')
            sage: dev._remote_branch_for_ticket(1)
            'u/doctest/foo'
        """
        if branch_or_ticket is None:
            from git_error import DetachedHeadError
            try:
                branch = self.git.current_branch()
            except DetachedHeadError:
                self._UI.error('You must specify "branch" in detached HEAD state.')
                self._UI.info(['', 'Use "{0}" to checkout a branch'],
                              self._format_command('checkout'))
                raise OperationCancelledError("detached head state")
        elif self._is_ticket_name(branch_or_ticket):
            ticket = self._ticket_from_ticket_name(branch_or_ticket)
            if not self._has_local_branch_for_ticket(ticket):
                self._UI.error('no local branch for ticket #{0} found. Cannot set remote branch'
                               ' for that ticket.', ticket)
                raise OperationCancelledError("no such ticket")
            branch = self._local_branch_for_ticket(ticket)
        else:
            branch = branch_or_ticket

        self._check_local_branch_name(branch, exists=True)
        self._check_remote_branch_name(remote_branch)

        # If we add restrictions on which branches users may push to, we should append them here.
        m = USER_BRANCH.match(remote_branch)
        if remote_branch == 'master' or m and m.groups()[0] != self.trac._username:
            self._UI.warning('The remote branch "{0}" is not in your user scope. You probably'
                             ' do not have permission to push to that branch.', remote_branch)
            self._UI.info(['', 'You can always use "u/{1}/{0}" as the remote branch name.'],
                          remote_branch, self.trac._username)

        self._set_remote_branch_for_branch(branch, remote_branch)

    def push(self, ticket=None, remote_branch=None, force=False):
        r"""
        Push the current branch to the Sage repository.

        INPUT:

        - ``ticket`` -- an integer or string identifying a ticket or ``None``
          (default: ``None``), if ``None`` and currently working on a ticket or
          if ``ticket`` specifies a ticket, then the branch on that ticket is
          set to ``remote_branch`` after the current branch has been pushed there.

        - ``remote_branch`` -- a string or ``None`` (default: ``None``), the remote
          branch to push to; if ``None``, then a default is chosen

        - ``force`` -- a boolean (default: ``False``), whether to push if
          this is not a fast-forward.

        .. SEEALSO::

            - :meth:`commit` -- Save changes to the local repository.

            - :meth:`pull` -- Update a ticket with changes from the remote
              repository.

        TESTS:

        Create a doctest setup with two users::

            sage: from sage.dev.test.sagedev import two_user_setup
            sage: alice, config_alice, bob, config_bob, server = two_user_setup()

        Alice tries to push to ticket 1 which does not exist yet::

            sage: alice._chdir()
            sage: alice.push(ticket=1)
            Ticket name "1" is not valid or ticket does not exist on trac.

        Alice creates ticket 1 and pushes some changes to it::

            sage: alice._UI.append("Summary: summary1\ndescription")
            sage: alice.create_ticket()
            Created ticket #1 at https://trac.sagemath.org/1.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=1" to create a new local branch)
            1
            sage: alice.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.
            sage: open("tracked", "w").close()
            sage: alice.git.super_silent.add("tracked")
            sage: alice.git.super_silent.commit(message="alice: added tracked")
            sage: alice._UI.append("y")
            sage: alice.push()
            The branch "u/alice/ticket/1" does not exist on the remote server.
            Create new remote branch? [Yes/no] y

        Now Bob can check that ticket out and push changes himself::

            sage: bob._chdir()
            sage: bob.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.
            sage: with open("tracked", "w") as f: f.write("bob")
            sage: bob.git.super_silent.add("tracked")
            sage: bob.git.super_silent.commit(message="bob: modified tracked")
            sage: bob._UI.append("y")
            sage: bob._UI.append("y")
            sage: bob.push()
            The branch "u/bob/ticket/1" does not exist on the remote server.
            Create new remote branch? [Yes/no] y
            The branch field of ticket #1 needs to be updated from its current value
            "u/alice/ticket/1" to "u/bob/ticket/1"
            Change the "Branch:" field? [Yes/no] y

        Now Alice can pull these changes::

            sage: alice._chdir()
            sage: alice.pull()
            Merging the remote branch "u/bob/ticket/1" into the local branch "ticket/1".
            Automatic merge successful.
            <BLANKLINE>
            #  (use "sage --dev commit" to commit your merge)

        Alice and Bob make non-conflicting changes simultaneously::

            sage: with open("tracked", "w") as f: f.write("alice")
            sage: alice.git.super_silent.add("tracked")
            sage: alice.git.super_silent.commit(message="alice: modified tracked")

            sage: bob._chdir()
            sage: open("tracked2", "w").close()
            sage: bob.git.super_silent.add("tracked2")
            sage: bob.git.super_silent.commit(message="bob: added tracked2")

        After Alice pushed her changes, Bob can not set the branch field anymore::

            sage: alice._chdir()
            sage: alice._UI.append("y")
            sage: alice._UI.append("y")
            sage: alice.push()
            Local commits that are not on the remote branch "u/alice/ticket/1":
            <BLANKLINE>
                ...: alice: modified tracked
                ...: bob: modified tracked
            <BLANKLINE>
            Push to remote branch? [Yes/no] y
            The branch field of ticket #1 needs to be updated from its current value
            "u/bob/ticket/1" to "u/alice/ticket/1"
            Change the "Branch:" field? [Yes/no] y

            sage: bob._chdir()
            sage: bob._UI.append("y")
            sage: bob.push()
            Local commits that are not on the remote branch "u/bob/ticket/1":
            <BLANKLINE>
                ....: bob: added tracked2
            <BLANKLINE>
            Push to remote branch? [Yes/no] y
            Not setting the branch field for ticket #1 to "u/bob/ticket/1" because
            "u/bob/ticket/1" and the current value of the branch field "u/alice/ticket/1"
            have diverged.
            <BLANKLINE>
            #  Use "sage --dev pull --ticket=1" to merge the changes introduced by the remote "u/alice/ticket/1" into your local branch.

        After merging the changes, this works again::

            sage: bob.pull()
            Merging the remote branch "u/alice/ticket/1" into the local branch "ticket/1".
            Automatic merge successful.
            <BLANKLINE>
            #  (use "sage --dev commit" to commit your merge)
            sage: bob._UI.append("y")
            sage: bob._UI.append("y")
            sage: bob.push()
            Local commits that are not on the remote branch "u/bob/ticket/1":
            <BLANKLINE>
                ...: Merge branch 'u/alice/ticket/1' of ... into ticket/1
                ...: alice: modified tracked
            <BLANKLINE>
            Push to remote branch? [Yes/no] y
            The branch field of ticket #1 needs to be updated from its current value
            "u/alice/ticket/1" to "u/bob/ticket/1"
            Change the "Branch:" field? [Yes/no] y

        Check that ``ticket`` works::

            sage: bob.push(2)
            Ticket name "2" is not valid or ticket does not exist on trac.

        After creating the ticket, this works with a warning::

            sage: bob._UI.append("Summary: summary2\ndescription")
            sage: bob.create_ticket()
            Created ticket #2 at https://trac.sagemath.org/2.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=2" to create a new local branch)
            2
            sage: bob.checkout(ticket=2)
            On ticket #2 with associated local branch "ticket/2".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.
            sage: bob.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.
            sage: bob._UI.append("y")
            sage: bob._UI.append("y")
            sage: bob.push(2)
            About to push the branch "ticket/1" to "u/bob/ticket/2" for ticket #2. However,
            your local branch for ticket #2 seems to be "ticket/2".
             Do you really want to proceed? [yes/No] y
            <BLANKLINE>
            #  Use "sage --dev checkout --ticket=2 --branch=ticket/1" to permanently set "ticket/1" as the branch associated to ticket #2.
            The branch "u/bob/ticket/2" does not exist on the remote server.
            Create new remote branch? [Yes/no] y

        Check that ``remote_branch`` works::

            sage: bob._UI.append("y")
            sage: bob._UI.append("y")
            sage: bob.push(remote_branch="u/bob/branch1")
            The branch "u/bob/branch1" does not exist on the remote server.
            Create new remote branch? [Yes/no] y
            The branch field of ticket #1 needs to be updated from its current value
            "u/bob/ticket/1" to "u/bob/branch1"
            Change the "Branch:" field? [Yes/no] y

        Check that dependencies are pushed correctly::

            sage: bob.merge(2)
            Merging the remote branch "u/bob/ticket/2" into the local branch "ticket/1".
            Automatic merge successful.
            <BLANKLINE>
            #  (use "sage --dev commit" to commit your merge)
            <BLANKLINE>
            Added dependency on #2 to #1.
            sage: with open("another_file", "w") as f: f.write("bob after merge(2)")
            sage: bob._UI.append('n')
            sage: bob.push()
            The branch field of ticket #1 needs to be updated from its current value
            "u/bob/branch1" to "u/bob/ticket/1"
            Change the "Branch:" field? [Yes/no] n
            sage: bob._UI.extend(['y', 'y', 'y'])
            sage: bob.commit(message="Bob's merge")  # oops
            The following files in your working directory are not tracked by git:
            <BLANKLINE>
                another_file
            <BLANKLINE>
            Start tracking any of these files? [yes/No] y
            Start tracking "another_file"? [yes/No] y
            Commit your changes to branch "ticket/1"? [Yes/no] y
            <BLANKLINE>
            #  Use "sage --dev push" to push your commits to the trac server once you are done.
            sage: bob._UI.extend(['y', 'y'])
            sage: bob.push()
            Local commits that are not on the remote branch "u/bob/ticket/1":
            <BLANKLINE>
                ...: Bob's merge
            <BLANKLINE>
            Push to remote branch? [Yes/no] y
            The branch field of ticket #1 needs to be updated from its current value
            "u/bob/branch1" to "u/bob/ticket/1"
            Change the "Branch:" field? [Yes/no] y
            Uploading your dependencies for ticket #1: "" => "#2"
            sage: bob._sagedev._set_dependencies_for_ticket(1,())
            sage: with open("another_file", "w") as f: f.write("bob after push")
            sage: bob._UI.extend(['y', 'y', 'y'])
            sage: bob.commit(message='another commit')
            Commit your changes to branch "ticket/1"? [Yes/no] y
            <BLANKLINE>
            #  Use "sage --dev push" to push your commits to the trac server once you are done.
            sage: bob._UI.extend(['y', "keep", 'y'])
            sage: bob.push()
            Local commits that are not on the remote branch "u/bob/ticket/1":
            <BLANKLINE>
                ...: another commit
            <BLANKLINE>
            Push to remote branch? [Yes/no] y
            Trac ticket #1 depends on #2 while your local branch depends on no tickets.
            Updating dependencies is recommended but optional.
            Action for dependencies? [upload/download/keep] keep
            sage: with open("another_file", "w") as f: f.write("bob after 2nd push")
            sage: bob._UI.append('y')
            sage: bob.commit(message='final commit')
            Commit your changes to branch "ticket/1"? [Yes/no] y
            <BLANKLINE>
            #  Use "sage --dev push" to push your commits to the trac server once you are done.

            sage: bob._UI.extend(['y', 'download', 'y'])
            sage: bob.push()
            Local commits that are not on the remote branch "u/bob/ticket/1":
            <BLANKLINE>
                ...: final commit
            <BLANKLINE>
            Push to remote branch? [Yes/no] y
            Trac ticket #1 depends on #2 while your local branch depends on no tickets.
            Updating dependencies is recommended but optional.
            Action for dependencies? [upload/download/keep] download
        """
        if ticket is None:
            ticket = self._current_ticket()
        if ticket is not None:
            ticket = self._ticket_from_ticket_name(ticket)
            self._check_ticket_name(ticket, exists=True)

        from git_error import DetachedHeadError
        try:
            branch = self.git.current_branch()
        except DetachedHeadError:
            self._UI.error("Cannot push while in detached HEAD state.")
            raise OperationCancelledError("cannot push while in detached HEAD state")

        if remote_branch is None:
            if ticket:
                remote_branch = self._remote_branch_for_ticket(ticket)
                if remote_branch is None:
                    raise SageDevValueError("remote_branch must be specified since #{0}"
                                            " has no remote branch set.".format(ticket))
            else:
                remote_branch = self._remote_branch_for_branch(branch)
                if remote_branch is None:
                    raise SageDevValueError("remote_branch must be specified since the"
                                            " current branch has no remote branch set.")
        self._check_remote_branch_name(remote_branch)

        # whether the user already confirmed that he really wants to push and set the branch field
        user_confirmation = False

        if ticket is not None:
            if self._has_local_branch_for_ticket(ticket) and self._local_branch_for_ticket(ticket) == branch:
                pass
            elif self._has_local_branch_for_ticket(ticket) and self._local_branch_for_ticket(ticket) != branch:
                self._UI.show('About to push the branch "{0}" to "{1}" for ticket #{2}.'
                              ' However, your local branch for ticket #{2} seems to be "{3}".',
                              branch, remote_branch, ticket, self._local_branch_for_ticket(ticket))
                user_confirmation = self._UI.confirm(' Do you really want to proceed?', default=False)
                if user_confirmation:
                    self._UI.info(['',
                                   'Use "{2}" to permanently set "{1}" as the branch'
                                   ' associated to ticket #{0}.'],
                                  ticket, branch, self._format_command("checkout",ticket=ticket,branch=branch))
                else:
                    raise OperationCancelledError("user requsted")
            elif self._has_ticket_for_local_branch(branch) and self._ticket_for_local_branch(branch) != ticket:
                self._UI.show('About to push the local branch "{0}" to remote branch "{1}" for'
                              ' ticket #{2}. However, that branch is already associated to ticket #{3}.',
                              branch, remote_branch, ticket, self._ticket_for_local_branch(branch))
                user_confirmation = self._UI.confirm(' Do you really want to proceed?', default=False)
                if user_confirmation:
                    self._UI.info(['', 'Use "{2}" to permanently set the branch associated to'
                                   ' ticket #{0} to "{1}". To create a new branch from "{1}" for'
                                   ' #{0}, use "{3}" and "{4}".'],
                                  ticket, branch,
                                  self._format_command("checkout",ticket=ticket,branch=branch),
                                  self._format_command("checkout",ticket=ticket),
                                  self._format_command("merge", branch=branch))

        self._UI.debug('Pushing your changes in "{0}" to "{1}".'.format(branch, remote_branch))
        try:
            remote_branch_exists = self._is_remote_branch_name(remote_branch, exists=True)
            if not remote_branch_exists:
                self._UI.show('The branch "{0}" does not exist on the remote server.', remote_branch)
                if not self._UI.confirm('Create new remote branch?', default=True):
                    raise OperationCancelledError("User did not want to create remote branch.")
            else:
                self.git.super_silent.fetch(self.git._repository_anonymous, remote_branch)

            # check whether force is necessary
            if remote_branch_exists and not self.git.is_child_of(branch, 'FETCH_HEAD'):
                if not force:
                    self._UI.error('Not pushing your changes because they would discard some of'
                                   ' the commits on the remote branch "{0}".', remote_branch)
                    raise OperationCancelledError("not a fast-forward")

            # check whether this is a nop
            if remote_branch_exists and not force and \
               self.git.commit_for_branch(branch) == self.git.commit_for_ref('FETCH_HEAD'):
                self._UI.debug('Remote branch "{0}" is idential to your local branch "{1}',
                              remote_branch, branch)
                self._UI.debug(['', '(use "{0}" to commit changes before pushing)'],
                               self._format_command("commit"))
            else:
                try:
                    if not force:
                        if remote_branch_exists:
                            commits = self.git.log("{0}..{1}".format('FETCH_HEAD', branch), '--pretty=%h: %s')
                            self._UI.show(['Local commits that are not on the remote branch "{0}":', ''] +
                                          ['    ' + c for c in commits.splitlines()] +
                                          [''], remote_branch)
                            if not self._UI.confirm('Push to remote branch?', default=True):
                                raise OperationCancelledError("user requested")

                    self._upload_ssh_key() # make sure that we have access to the repository
                    self.git.super_silent.push(self.git._repository,
                                               "{0}:{1}".format(branch, remote_branch),
                                               force=force)
                except GitError as e:
                    # can we give any advice if this fails?
                    raise
                self._UI.debug('Changes in "{0}" have been pushed to "{1}".'.format(branch, remote_branch))
        except OperationCancelledError:
            self._UI.debug("Did not push any changes.")
            raise

        if ticket:
            current_remote_branch = self.trac._branch_for_ticket(ticket)
            if current_remote_branch == remote_branch:
                self._UI.debug('Not setting the branch field for ticket #{0} because it already'
                               ' points to your branch "{1}".'.format(ticket, remote_branch))
            else:
                self._UI.debug('Setting the branch field of ticket #{0} to "{1}".'.format(ticket, remote_branch))
                if current_remote_branch is not None:
                    self.git.super_silent.fetch(self.git._repository_anonymous, current_remote_branch)
                    if force or self.git.is_ancestor_of('FETCH_HEAD', branch):
                        pass
                    else:
                        self._UI.error('Not setting the branch field for ticket #{0} to "{1}" because'
                                       ' "{1}" and the current value of the branch field "{2}" have diverged.'
                                       .format(ticket, remote_branch, current_remote_branch))
                        self._UI.info(['',
                                       'Use "{0}" to merge the changes introduced by'
                                       ' the remote "{1}" into your local branch.'],
                                      self._format_command("pull", ticket=ticket),
                                      current_remote_branch)
                        raise OperationCancelledError("not a fast-forward")

                if current_remote_branch is not None and not force and not user_confirmation:
                    self._UI.show('The branch field of ticket #{0} needs to be'
                                  ' updated from its current value "{1}" to "{2}"'
                                  ,ticket, current_remote_branch, remote_branch)
                    if not self._UI.confirm('Change the "Branch:" field?', default=True):
                        raise OperationCancelledError("user requested")

                attributes = self.trac._get_attributes(ticket)
                attributes['branch'] = remote_branch
                self.trac._authenticated_server_proxy.ticket.update(ticket, "", attributes)

        if ticket and self._has_ticket_for_local_branch(branch):
            old_dependencies_ = self.trac.dependencies(ticket)
            old_dependencies = ", ".join(["#"+str(dep) for dep in old_dependencies_])
            new_dependencies_ = self._dependencies_for_ticket(self._ticket_for_local_branch(branch))
            new_dependencies = ", ".join(["#"+str(dep) for dep in new_dependencies_])

            upload = True
            if old_dependencies != new_dependencies:
                if old_dependencies:
                    self._UI.show('Trac ticket #{0} depends on {1} while your local branch depends'
                                  ' on {2}. Updating dependencies is recommended but optional.',
                                  ticket, old_dependencies, new_dependencies or "no tickets"),
                    sel = self._UI.select('Action for dependencies?', options=("upload", "download", "keep"))
                    if sel == "keep":
                        upload = False
                    elif sel == "download":
                        self._set_dependencies_for_ticket(ticket, old_dependencies_)
                        self._UI.debug("Setting dependencies for #{0} to {1}.", ticket, old_dependencies)
                        upload = False
                    elif sel == "upload":
                        pass
                    else:
                        raise NotImplementedError
            else:
                self._UI.debug("Not uploading your dependencies for ticket #{0} because the"
                               " dependencies on trac are already up-to-date.", ticket)
                upload = False

            if upload:
                self._UI.show('Uploading your dependencies for ticket #{0}: "{1}" => "{2}"',
                              ticket, old_dependencies, new_dependencies)
                attributes = self.trac._get_attributes(ticket)
                attributes['dependencies'] = new_dependencies
                # Don't send an e-mail notification
                self.trac._authenticated_server_proxy.ticket.update(ticket, "", attributes)

    def reset_to_clean_state(self, error_unless_clean=True, helpful=True):
        r"""
        Reset the current working directory to a clean state.

        INPUT:

        - ``error_unless_clean`` -- a boolean (default: ``True``),
          whether to raise an
          :class:`user_interface_error.OperationCancelledError` if the
          directory remains in an unclean state; used internally.

        TESTS:

        Set up a single user for doctesting::

            sage: from sage.dev.test.sagedev import single_user_setup
            sage: dev, config, UI, server = single_user_setup()
            sage: dev._wrap("reset_to_clean_state")

        Nothing happens if the directory is already clean::

            sage: dev.reset_to_clean_state()

        Bring the directory into a non-clean state::

            sage: dev.git.super_silent.checkout(b="branch1")
            sage: with open("tracked", "w") as f: f.write("boo")
            sage: dev.git.silent.add("tracked")
            sage: dev.git.silent.commit(message="added tracked")

            sage: dev.git.super_silent.checkout('HEAD~')
            sage: dev.git.super_silent.checkout(b="branch2")
            sage: with open("tracked", "w") as f: f.write("foo")
            sage: dev.git.silent.add("tracked")
            sage: dev.git.silent.commit(message="added tracked")
            sage: from sage.dev.git_error import GitError
            sage: try:
            ....:     dev.git.silent.merge("branch1")
            ....: except GitError: pass
            sage: UI.append("cancel")
            sage: dev.reset_to_clean_state()
            Repository is in an unclean state (merge). Resetting the state will discard any
            uncommited changes.
            Reset repository? [reset/Cancel] cancel
            <BLANKLINE>
            #  (use "sage --dev commit" to save changes in a new commit)
            sage: UI.append("reset")
            sage: dev.reset_to_clean_state()
            Repository is in an unclean state (merge). Resetting the state will discard any
            uncommited changes.
            Reset repository? [reset/Cancel] reset
            sage: dev.reset_to_clean_state()

        A detached HEAD does not count as a non-clean state::

            sage: dev.git.super_silent.checkout('HEAD', detach=True)
            sage: dev.reset_to_clean_state()
        """
        states = self.git.get_state()
        if not states:
            return
        self._UI.show('Repository is in an unclean state ({0}).'
                      ' Resetting the state will discard any uncommited changes.',
                      ', '.join(states))
        sel = self._UI.select('Reset repository?',
                              options=('reset', 'cancel'), default=1)
        if sel == 'cancel':
            if not error_unless_clean:
                return
            if helpful:
                self._UI.info(['', '(use "{0}" to save changes in a new commit)'],
                              self._format_command("commit"))
            raise OperationCancelledError("User requested not to clean the current state.")
        elif sel == 'reset':
            self.git.reset_to_clean_state()
        else:
            assert False

    def clean(self, error_unless_clean=True):
        r"""
        Restore the working directory to the most recent commit.

        INPUT:

        - ``error_unless_clean`` -- a boolean (default: ``True``),
          whether to raise an
          :class:`user_interface_error.OperationCancelledError` if the
          directory remains in an unclean state; used internally.

        TESTS:

        Set up a single user for doctesting::

            sage: from sage.dev.test.sagedev import single_user_setup
            sage: dev, config, UI, server = single_user_setup()

        Check that nothing happens if there no changes::

            sage: dev.clean()

        Check that nothing happens if there are only untracked files::

            sage: open("untracked","w").close()
            sage: dev.clean()

        Uncommitted changes can simply be dropped::

            sage: open("tracked","w").close()
            sage: dev.git.silent.add("tracked")
            sage: dev.git.silent.commit(message="added tracked")
            sage: with open("tracked", "w") as f: f.write("foo")
            sage: UI.append("discard")
            sage: dev.clean()
            The following files in your working directory contain uncommitted changes:
            <BLANKLINE>
                 tracked
            <BLANKLINE>
            Discard changes? [discard/Cancel/stash] discard
            sage: dev.clean()

        Uncommitted changes can be kept::

            sage: with open("tracked", "w") as f: f.write("foo")
            sage: UI.append("cancel")
            sage: dev.clean()
            The following files in your working directory contain uncommitted changes:
            <BLANKLINE>
                 tracked
            <BLANKLINE>
            Discard changes? [discard/Cancel/stash] cancel

        Or stashed::

            sage: UI.append("stash")
            sage: dev.clean()
            The following files in your working directory contain uncommitted changes:
            <BLANKLINE>
                 tracked
            <BLANKLINE>
            Discard changes? [discard/Cancel/stash] stash
            Your changes have been moved to the git stash stack. To re-apply your changes
            later use "git stash apply".
            sage: dev.clean()
        """
        try:
            self.reset_to_clean_state(error_unless_clean)
        except OperationCancelledError:
            self._UI.error("Can not clean the working directory unless in a clean state.")
            raise

        if not self.git.has_uncommitted_changes():
            return

        files = [line[2:] for line in self.git.status(porcelain=True).splitlines()
                 if not line.startswith('?')]

        self._UI.show(
            ['The following files in your working directory contain uncommitted changes:'] +
            [''] +
            ['    ' + f for f in files ] +
            [''])
        cancel = 'cancel' if error_unless_clean else 'keep'
        sel = self._UI.select('Discard changes?',
                              options=('discard', cancel, 'stash'), default=1)
        if sel == 'discard':
            self.git.clean_wrapper()
        elif sel == cancel:
            if error_unless_clean:
                raise OperationCancelledError("User requested not to clean the working directory.")
        elif sel == 'stash':
            self.git.super_silent.stash()
            self._UI.show('Your changes have been moved to the git stash stack. '
                          'To re-apply your changes later use "git stash apply".')
        else:
            assert False

    def edit_ticket(self, ticket=None):
        r"""
        Edit the description of ``ticket`` on trac.

        INPUT:

        - ``ticket`` -- an integer or string identifying a ticket or ``None``
          (default: ``None``), the number of the ticket to edit. If ``None``,
          edit the :meth:`_current_ticket`.

        .. SEEALSO::

            :meth:`create_ticket`, :meth:`comment`,
            :meth:`set_needs_review`, :meth:`set_needs_work`,
            :meth:`set_positive_review`, :meth:`set_needs_info`

        TESTS:

        Set up a single user for doctesting::

            sage: from sage.dev.test.sagedev import single_user_setup
            sage: dev, config, UI, server = single_user_setup()

        Create a ticket and edit it::

            sage: UI.append("Summary: summary1\ndescription")
            sage: dev.create_ticket()
            Created ticket #1 at https://trac.sagemath.org/1.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=1" to create a new local branch)
            1
            sage: dev.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.
            sage: UI.append("Summary: summary1\ndescription...")
            sage: dev.edit_ticket()
            sage: dev.trac._get_attributes(1)
            {'description': 'description...', 'summary': 'summary1'}
        """
        if ticket is None:
            ticket = self._current_ticket()

        if ticket is None:
            raise SageDevValueError("ticket must be specified if not currently on a ticket.")

        self._check_ticket_name(ticket, exists=True)
        ticket = self._ticket_from_ticket_name(ticket)
        self.trac.edit_ticket_interactive(ticket)

    def needs_review(self, ticket=None, comment=''):
        r"""
        Set a ticket on trac to ``needs_review``.

        INPUT:

        - ``ticket`` -- an integer or string identifying a ticket or
          ``None`` (default: ``None``), the number of the ticket to
          edit.  If ``None``, edit the :meth:`_current_ticket`.

        - ``comment`` -- a comment to go with the status change.

        .. SEEALSO::

            :meth:`edit_ticket`, :meth:`set_needs_work`,
            :meth:`set_positive_review`, :meth:`comment`,
            :meth:`set_needs_info`

        TESTS:

        Set up a single user for doctesting::

            sage: from sage.dev.test.sagedev import single_user_setup
            sage: dev, config, UI, server = single_user_setup()

        Create a ticket and set it to needs_review::

            sage: UI.append("Summary: summary1\ndescription")
            sage: dev.create_ticket()
            Created ticket #1 at https://trac.sagemath.org/1.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=1" to create a new local branch)
            1
            sage: dev.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.
            sage: open("tracked", "w").close()
            sage: dev.git.super_silent.add("tracked")
            sage: dev.git.super_silent.commit(message="alice: added tracked")
            sage: dev._UI.append("y")
            sage: dev.push()
            The branch "u/doctest/ticket/1" does not exist on the remote server.
            Create new remote branch? [Yes/no] y
            sage: dev.needs_review(comment='Review my ticket!')
            sage: dev.trac._get_attributes(1)['status']
            'needs_review'
        """
        if ticket is None:
            ticket = self._current_ticket()
        if ticket is None:
            raise SageDevValueError("ticket must be specified if not currently on a ticket.")
        self._check_ticket_name(ticket, exists=True)
        self.trac.set_attributes(ticket, comment, notify=True, status='needs_review')
        self._UI.debug("Ticket #%s marked as needing review"%ticket)

    def needs_work(self, ticket=None, comment=''):
        r"""
        Set a ticket on trac to ``needs_work``.

        INPUT:

        - ``ticket`` -- an integer or string identifying a ticket or
          ``None`` (default: ``None``), the number of the ticket to
          edit.  If ``None``, edit the :meth:`_current_ticket`.

        - ``comment`` -- a comment to go with the status change.

        .. SEEALSO::

            :meth:`edit_ticket`, :meth:`set_needs_review`,
            :meth:`set_positive_review`, :meth:`comment`,
            :meth:`set_needs_info`

        TESTS:

        Create a doctest setup with two users::

            sage: from sage.dev.test.sagedev import two_user_setup
            sage: alice, config_alice, bob, config_bob, server = two_user_setup()

        Alice creates a ticket and set it to needs_review::

            sage: alice._chdir()
            sage: alice._UI.append("Summary: summary1\ndescription")
            sage: alice.create_ticket()
            Created ticket #1 at https://trac.sagemath.org/1.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=1" to create a new local branch)
            1
            sage: alice.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.
            sage: open("tracked", "w").close()
            sage: alice.git.super_silent.add("tracked")
            sage: alice.git.super_silent.commit(message="alice: added tracked")
            sage: alice._UI.append("y")
            sage: alice.push()
            The branch "u/alice/ticket/1" does not exist on the remote server.
            Create new remote branch? [Yes/no] y
            sage: alice.needs_review(comment='Review my ticket!')

        Bob reviews the ticket and finds it lacking::

            sage: bob._chdir()
            sage: bob.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.
            sage: bob.needs_work(comment='Need to add an untracked file!')
            sage: bob.trac._get_attributes(1)['status']
            'needs_work'
        """
        if ticket is None:
            ticket = self._current_ticket()
        if ticket is None:
            raise SageDevValueError("ticket must be specified if not currently on a ticket.")
        self._check_ticket_name(ticket, exists=True)
        if not comment:
            comment = self._UI.get_input("Please add a comment for the author:")
        self.trac.set_attributes(ticket, comment, notify=True, status='needs_work')
        self._UI.debug("Ticket #%s marked as needing work"%ticket)

    def needs_info(self, ticket=None, comment=''):
        r"""
        Set a ticket on trac to ``needs_info``.

        INPUT:

        - ``ticket`` -- an integer or string identifying a ticket or
          ``None`` (default: ``None``), the number of the ticket to
          edit.  If ``None``, edit the :meth:`_current_ticket`.

        - ``comment`` -- a comment to go with the status change.

        .. SEEALSO::

            :meth:`edit_ticket`, :meth:`needs_review`,
            :meth:`positive_review`, :meth:`comment`,
            :meth:`needs_work`

        TESTS:

        Create a doctest setup with two users::

            sage: from sage.dev.test.sagedev import two_user_setup
            sage: alice, config_alice, bob, config_bob, server = two_user_setup()

        Alice creates a ticket and set it to needs_review::

            sage: alice._chdir()
            sage: alice._UI.append("Summary: summary1\ndescription")
            sage: alice.create_ticket()
            Created ticket #1 at https://trac.sagemath.org/1.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=1" to create a new local branch)
            1
            sage: alice.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.
            sage: open("tracked", "w").close()
            sage: alice.git.super_silent.add("tracked")
            sage: alice.git.super_silent.commit(message="alice: added tracked")
            sage: alice._UI.append("y")
            sage: alice.push()
            The branch "u/alice/ticket/1" does not exist on the remote server.
            Create new remote branch? [Yes/no] y
            sage: alice.needs_review(comment='Review my ticket!')

        Bob reviews the ticket and finds it lacking::

            sage: bob._chdir()
            sage: bob.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.
            sage: bob.needs_info(comment='Why is a tracked file enough?')
            sage: bob.trac._get_attributes(1)['status']
            'needs_info'
        """
        if ticket is None:
            ticket = self._current_ticket()
        if ticket is None:
            raise SageDevValueError("ticket must be specified if not currently on a ticket.")
        self._check_ticket_name(ticket, exists=True)
        if not comment:
            comment = self._UI.get_input("Please specify what information is required from the author:")
        self.trac.set_attributes(ticket, comment, notify=True, status='needs_info')
        self._UI.debug("Ticket #%s marked as needing info"%ticket)

    def positive_review(self, ticket=None, comment=''):
        r"""
        Set a ticket on trac to ``positive_review``.

        INPUT:

        - ``ticket`` -- an integer or string identifying a ticket or
          ``None`` (default: ``None``), the number of the ticket to
          edit.  If ``None``, edit the :meth:`_current_ticket`.

        - ``comment`` -- a comment to go with the status change.

        .. SEEALSO::

            :meth:`edit_ticket`, :meth:`needs_review`,
            :meth:`needs_info`, :meth:`comment`,
            :meth:`needs_work`

        TESTS:

        Create a doctest setup with two users::

            sage: from sage.dev.test.sagedev import two_user_setup
            sage: alice, config_alice, bob, config_bob, server = two_user_setup()

        Alice creates a ticket and set it to needs_review::

            sage: alice._chdir()
            sage: alice._UI.append("Summary: summary1\ndescription")
            sage: alice.create_ticket()
            Created ticket #1 at https://trac.sagemath.org/1.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=1" to create a new local branch)
            1
            sage: alice.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.
            sage: open("tracked", "w").close()
            sage: alice.git.super_silent.add("tracked")
            sage: alice.git.super_silent.commit(message="alice: added tracked")
            sage: alice._UI.append("y")
            sage: alice.push()
            The branch "u/alice/ticket/1" does not exist on the remote server.
            Create new remote branch? [Yes/no] y
            sage: alice.needs_review(comment='Review my ticket!')

        Bob reviews the ticket and finds it good::

            sage: bob._chdir()
            sage: bob.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.
            sage: bob.positive_review()
            sage: bob.trac._get_attributes(1)['status']
            'positive_review'
        """
        if ticket is None:
            ticket = self._current_ticket()
        if ticket is None:
            raise SageDevValueError("ticket must be specified if not currently on a ticket.")
        self._check_ticket_name(ticket, exists=True)
        self.trac.set_attributes(ticket, comment, notify=True, status='positive_review')
        self._UI.debug("Ticket #%s reviewed!"%ticket)

    def comment(self, ticket=None):
        r"""
        Add a comment to ``ticket`` on trac.

        INPUT:

        - ``ticket`` -- an integer or string identifying a ticket or ``None``
          (default: ``None``), the number of the ticket to edit. If ``None``,
          edit the :meth:`_current_ticket`.

        .. SEEALSO::

            :meth:`create_ticket`, :meth:`edit_ticket`

        TESTS:

        Set up a single user for doctesting::

            sage: from sage.dev.test.sagedev import single_user_setup
            sage: dev, config, UI, server = single_user_setup()

        Create a ticket and add a comment::

            sage: UI.append("Summary: summary1\ndescription")
            sage: dev.create_ticket()
            Created ticket #1 at https://trac.sagemath.org/1.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=1" to create a new local branch)
            1
            sage: dev.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.
            sage: UI.append("comment")
            sage: dev.comment()
            sage: server.tickets[1].comments
            ['comment']
        """
        if ticket is None:
            ticket = self._current_ticket()

        if ticket is None:
            raise SageDevValueError("ticket must be specified if not currently on a ticket.")

        self._check_ticket_name(ticket, exists=True)
        ticket = self._ticket_from_ticket_name(ticket)
        self.trac.add_comment_interactive(ticket)

    def browse_ticket(self, ticket=None):
        r"""
        Start a webbrowser at the ticket page on trac.

        INPUT:

        - ``ticket`` -- an integer or string identifying a ticket or ``None``
          (default: ``None``), the number of the ticket to edit. If ``None``,
          browse the :meth:`_current_ticket`.

        .. SEEALSO::

            :meth:`edit_ticket`, :meth:`comment`,
            :meth:`sage.dev.trac_interface.TracInterface.show_ticket`,
            :meth:`sage.dev.trac_interface.TracInterface.show_comments`

        EXAMPLES::

            sage: dev.browse_ticket(10000) # not tested
        """
        if ticket is None:
            ticket = self._current_ticket()

        if ticket is None:
            raise SageDevValueError("ticket must be specified if not currently on a ticket.")

        self._check_ticket_name(ticket, exists=True)
        ticket = self._ticket_from_ticket_name(ticket)

        from sage.misc.viewer import browser
        from sage.env import TRAC_SERVER_URI
        browser_cmdline = browser() + ' ' + TRAC_SERVER_URI + '/ticket/' + str(ticket)
        import os
        os.system(browser_cmdline)

    def remote_status(self, ticket=None):
        r"""
        Show information about the status of ``ticket``.

        INPUT:

        - ``ticket`` -- an integer or string identifying a ticket or ``None``
          (default: ``None``), the number of the ticket to edit.  If ``None``,
          show information for the :meth:`_current_ticket`.

        TESTS:

        Set up a single user for doctesting::

            sage: from sage.dev.test.sagedev import single_user_setup
            sage: dev, config, UI, server = single_user_setup()

        It is an error to call this without parameters if not on a ticket::

            sage: dev.remote_status()
            ticket must be specified if not currently on a ticket.

        Create a ticket and show its remote status::

            sage: UI.append("Summary: ticket1\ndescription")
            sage: dev.create_ticket()
            Created ticket #1 at https://trac.sagemath.org/1.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=1" to create a new local branch)
            1
            sage: dev.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.
            sage: dev.remote_status()
            Ticket #1 (https://trac.sagemath.org/ticket/1)
            ==============================================
            Your branch "ticket/1" has 0 commits.
            No branch has been set on the trac ticket yet.
            You have not created a remote branch yet.

        After pushing the local branch::

            sage: UI.append("y")
            sage: dev.push()
            The branch "u/doctest/ticket/1" does not exist on the remote server.
            Create new remote branch? [Yes/no] y
            sage: dev.remote_status()
            Ticket #1 (https://trac.sagemath.org/ticket/1)
            ==============================================
            Your branch "ticket/1" has 0 commits.
            The trac ticket points to the branch "u/doctest/ticket/1" which has 0 commits. It does not differ from "ticket/1".

        Making local changes::

            sage: open("tracked", "w").close()
            sage: dev.git.silent.add("tracked")
            sage: dev.git.silent.commit(message="added tracked")
            sage: dev.remote_status()
            Ticket #1 (https://trac.sagemath.org/ticket/1)
            ==============================================
            Your branch "ticket/1" has 1 commits.
            The trac ticket points to the branch "u/doctest/ticket/1" which has 0 commits. "ticket/1" is ahead of "u/doctest/ticket/1" by 1 commits:
            ...: added tracked

        Pushing them::

            sage: UI.append("y")
            sage: dev.push()
            Local commits that are not on the remote branch "u/doctest/ticket/1":
            <BLANKLINE>
                ...: added tracked
            <BLANKLINE>
            Push to remote branch? [Yes/no] y
            sage: dev.remote_status()
            Ticket #1 (https://trac.sagemath.org/ticket/1)
            ==============================================
            Your branch "ticket/1" has 1 commits.
            The trac ticket points to the branch "u/doctest/ticket/1" which has 1 commits. It does not differ from "ticket/1".

        The branch on the ticket is ahead of the local branch::

            sage: dev.git.silent.reset('HEAD~', hard=True)
            sage: dev.remote_status()
            Ticket #1 (https://trac.sagemath.org/ticket/1)
            ==============================================
            Your branch "ticket/1" has 0 commits.
            The trac ticket points to the branch "u/doctest/ticket/1" which has 1 commits. "u/doctest/ticket/1" is ahead of "ticket/1" by 1 commits:
            ...: added tracked

        A mixed case::

            sage: open("tracked2", "w").close()
            sage: dev.git.silent.add("tracked2")
            sage: dev.git.silent.commit(message="added tracked2")
            sage: open("tracked3", "w").close()
            sage: dev.git.silent.add("tracked3")
            sage: dev.git.silent.commit(message="added tracked3")
            sage: open("tracked4", "w").close()
            sage: dev.git.silent.add("tracked4")
            sage: dev.git.silent.commit(message="added tracked4")
            sage: dev._UI.append("y")
            sage: dev.push(remote_branch="u/doctest/branch1", force=True)
            The branch "u/doctest/branch1" does not exist on the remote server.
            Create new remote branch? [Yes/no] y
            sage: dev.git.silent.reset('HEAD~', hard=True)
            sage: dev.remote_status()
            Ticket #1 (https://trac.sagemath.org/ticket/1)
            ==============================================
            Your branch "ticket/1" has 2 commits.
            The trac ticket points to the branch "u/doctest/branch1" which has 3 commits. "u/doctest/branch1" is ahead of "ticket/1" by 1 commits:
            ...: added tracked4
            Your remote branch "u/doctest/ticket/1" has 1 commits. The branches "u/doctest/ticket/1" and "ticket/1" have diverged.
            "u/doctest/ticket/1" is ahead of "ticket/1" by 1 commits:
            ...: added tracked
            "ticket/1" is ahead of "u/doctest/ticket/1" by 2 commits:
            ...: added tracked2
            ...: added tracked3
        """
        if ticket is None:
            ticket = self._current_ticket()

        if ticket is None:
            raise SageDevValueError("ticket must be specified if not currently on a ticket.")

        self._check_ticket_name(ticket, exists=True)
        ticket = self._ticket_from_ticket_name(ticket)

        self._is_master_uptodate(action_if_not="warning")

        from sage.env import TRAC_SERVER_URI
        header = "Ticket #{0} ({1})".format(ticket, TRAC_SERVER_URI + '/ticket/' + str(ticket))
        underline = "="*len(header)

        commits = lambda a, b: list(reversed(
            self.git.log("{0}..{1}".format(a,b), "--pretty=%an <%ae>: %s").splitlines()))

        def detail(a, b, a_to_b, b_to_a):
            if not a_to_b and not b_to_a:
                return 'It does not differ from "{0}".'.format(b)
            elif not a_to_b:
                return '"{0}" is ahead of "{1}" by {2} commits:\n{3}'.format(a,b,len(b_to_a), "\n".join(b_to_a))
            elif not b_to_a:
                return '"{0}" is ahead of "{1}" by {2} commits:\n{3}'.format(b,a,len(a_to_b),"\n".join(a_to_b))
            else:
                return ('The branches "{0}" and "{1}" have diverged.\n"{0}" is ahead of'
                        ' "{1}" by {2} commits:\n{3}\n"{1}" is ahead of "{0}" by {4}'
                        ' commits:\n{5}'.format(a, b, len(b_to_a), "\n".join(b_to_a),
                                                len(a_to_b), "\n".join(a_to_b)))

        branch = None
        merge_base_local = None
        if self._has_local_branch_for_ticket(ticket):
            branch = self._local_branch_for_ticket(ticket)
            merge_base_local = self.git.merge_base(MASTER_BRANCH, branch).splitlines()[0]
            master_to_branch = commits(merge_base_local, branch)
            local_summary = 'Your branch "{0}" has {1} commits.'.format(branch, len(master_to_branch))
        else:
            local_summary = "You have no local branch for this ticket"

        ticket_branch = self.trac._branch_for_ticket(ticket)
        if ticket_branch:
            ticket_to_local = None
            local_to_ticket = None
            if not self._is_remote_branch_name(ticket_branch, exists=True):
                ticket_summary = 'The trac ticket points to the branch "{0}" which does not exist.'
            else:
                self.git.super_silent.fetch(self.git._repository_anonymous, ticket_branch)
                merge_base_ticket = self.git.merge_base(MASTER_BRANCH, 'FETCH_HEAD').splitlines()[0]
                master_to_ticket = commits(merge_base_ticket, 'FETCH_HEAD')
                ticket_summary = 'The trac ticket points to the' \
                    ' branch "{0}" which has {1} commits.'.format(ticket_branch, len(master_to_ticket))
                if branch is not None:
                    if merge_base_local != merge_base_ticket:
                        ticket_summary += ' The branch can not be compared to your local' \
                            ' branch "{0}" because the branches are based on different versions' \
                            ' of sage (i.e. the "master" branch).'
                    else:
                        ticket_to_local = commits('FETCH_HEAD', branch)
                        local_to_ticket = commits(branch, 'FETCH_HEAD')
                        ticket_summary += " "+detail(ticket_branch, branch, ticket_to_local, local_to_ticket)
        else:
            ticket_summary = "No branch has been set on the trac ticket yet."

        remote_branch = self._remote_branch_for_ticket(ticket)
        if self._is_remote_branch_name(remote_branch, exists=True):
            remote_to_local = None
            local_to_remote = None
            self.git.super_silent.fetch(self.git._repository_anonymous, remote_branch)
            merge_base_remote = self.git.merge_base(MASTER_BRANCH, 'FETCH_HEAD').splitlines()[0]
            master_to_remote = commits(merge_base_remote, 'FETCH_HEAD')
            remote_summary = 'Your remote branch "{0}" has {1} commits.'.format(
                remote_branch, len(master_to_remote))
            if branch is not None:
                if merge_base_remote != merge_base_local:
                    remote_summary += ' The branch can not be compared to your local' \
                        ' branch "{0}" because the branches are based on different version' \
                        ' of sage (i.e. the "master" branch).'
                else:
                    remote_to_local = commits('FETCH_HEAD', branch)
                    local_to_remote = commits(branch, 'FETCH_HEAD')
                    remote_summary += " "+detail(remote_branch, branch, remote_to_local, local_to_remote)
        else:
            remote_summary = "You have not created a remote branch yet."

        show = [header, underline, local_summary, ticket_summary]
        if not self._is_remote_branch_name(remote_branch, exists=True) or remote_branch != ticket_branch:
            show.append(remote_summary)
        self._UI.show("\n".join(show))

    def prune_tickets(self):
        r"""
        Remove branches for tickets that are already merged into master.

        .. SEEALSO::

            :meth:`abandon` -- Abandon a single ticket or branch.

        TESTS:

        Create a single user for doctesting::

            sage: from sage.dev.test.sagedev import single_user_setup
            sage: dev, config, UI, server = single_user_setup()

        Create a ticket branch::

            sage: UI.append("Summary: summary\ndescription")
            sage: dev.create_ticket()
            Created ticket #1 at https://trac.sagemath.org/1.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=1" to create a new local branch)
            1
            sage: dev.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.
            sage: dev.tickets()
                : master
            * #1: ticket/1 summary

        With a commit on it, the branch is not abandoned::

            sage: open("tracked","w").close()
            sage: dev.git.silent.add("tracked")
            sage: dev.git.super_silent.commit(message="added tracked")
            sage: dev.prune_tickets()
            sage: dev.tickets()
                : master
            * #1: ticket/1 summary

        After merging it to the master branch, it is abandoned. This does not
        work, because we cannot move the current branch::

            sage: dev.git.super_silent.checkout("master")
            sage: dev.git.super_silent.merge("ticket/1")

            sage: dev.git.super_silent.checkout("ticket/1")
            sage: dev.prune_tickets()
            Abandoning #1.
            Cannot delete "ticket/1": is the current branch.
            <BLANKLINE>
            #  (use "sage --dev vanilla" to switch to the master branch)

        Now, the branch is abandoned::

            sage: dev.vanilla()
            sage: dev.prune_tickets()
            Abandoning #1.
            Moved your branch "ticket/1" to "trash/ticket/1".
            sage: dev.tickets()
            : master
            sage: dev.prune_tickets()
        """
        for branch in self.git.local_branches():
            if self._has_ticket_for_local_branch(branch):
                ticket = self._ticket_for_local_branch(branch)
                if self.git.is_ancestor_of(branch, MASTER_BRANCH):
                    self._UI.show("Abandoning #{0}.".format(ticket))
                    self.abandon(ticket, helpful=False)

    def abandon(self, ticket_or_branch=None, helpful=True):
        r"""
        Abandon a ticket or branch.

        INPUT:

        - ``ticket_or_branch`` -- an integer or string identifying a ticket or
          the name of a local branch or ``None`` (default: ``None``), remove
          the branch ``ticket_or_branch`` or the branch for the ticket
          ``ticket_or_branch`` (or the current branch if ``None``). Also
          removes the users remote tracking branch.

        - ``helpful`` -- boolean (default: ``True``). Whether to print
          informational messages to guide new users.

        .. SEEALSO::

            - :meth:`prune_tickets` -- abandon tickets that have
              been closed.

            - :meth:`tickets` -- list local non-abandoned tickets.

        TESTS:

        Create a single user for doctesting::

            sage: from sage.dev.test.sagedev import single_user_setup
            sage: dev, config, UI, server = single_user_setup()

        Create a ticket branch and abandon it::

            sage: UI.append("Summary: summary\ndescription")
            sage: dev.create_ticket()
            Created ticket #1 at https://trac.sagemath.org/1.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=1" to create a new local branch)
            1
            sage: dev.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.
            sage: UI.append("y")
            sage: dev.push()
            The branch "u/doctest/ticket/1" does not exist on the remote server.
            Create new remote branch? [Yes/no] y
            sage: dev.abandon(1)
            Cannot delete "ticket/1": is the current branch.
            <BLANKLINE>
            #  (use "sage --dev vanilla" to switch to the master branch)
            sage: dev.vanilla()
            sage: dev.abandon(1)
            Moved your branch "ticket/1" to "trash/ticket/1".
            <BLANKLINE>
            #  Use "sage --dev checkout --ticket=1 --base=master" to restart working on #1 with a clean copy of the master branch.

        Start to work on a new branch for this ticket::

            sage: from sage.dev.sagedev import MASTER_BRANCH
            sage: UI.append("y")
            sage: dev.checkout(ticket=1, base=MASTER_BRANCH)
            About to create a new branch for #1 based on "master". However, the trac ticket
            for #1 already refers to the branch "u/doctest/ticket/1". The new branch will
            not contain any work that has already been done on "u/doctest/ticket/1".
            Create fresh branch? [yes/No] y
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.
        """
        ticket = None

        if self._is_ticket_name(ticket_or_branch):
            ticket = self._ticket_from_ticket_name(ticket_or_branch)

            if not self._has_local_branch_for_ticket(ticket):
                raise SageDevValueError("Cannot abandon #{0}: no local branch for this ticket.", ticket)
            ticket_or_branch = self._local_branch_for_ticket(ticket)

        if self._has_ticket_for_local_branch(ticket_or_branch):
            ticket = self._ticket_for_local_branch(ticket_or_branch)

        if self._is_local_branch_name(ticket_or_branch):
            branch = ticket_or_branch
            self._check_local_branch_name(branch, exists=True)

            if branch == MASTER_BRANCH:
                self._UI.error("Cannot delete the master branch.")
                raise OperationCancelledError("protecting the user")

            from git_error import DetachedHeadError
            try:
                if self.git.current_branch() == branch:
                    self._UI.error('Cannot delete "{0}": is the current branch.', branch)
                    self._UI.info(['', '(use "{0}" to switch to the master branch)'],
                                  self._format_command("vanilla"))
                    raise OperationCancelledError("can not delete current branch")
            except DetachedHeadError:
                pass

            new_branch = self._new_local_branch_for_trash(branch)
            self.git.super_silent.branch("-m", branch, new_branch)
            self._UI.show('Moved your branch "{0}" to "{1}".', branch, new_branch)
        else:
            raise SageDevValueError("ticket_or_branch must be the name of a ticket or a local branch")

        if ticket:
            self._set_local_branch_for_ticket(ticket, None)
            self._set_dependencies_for_ticket(ticket, None)
            if helpful:
                self._UI.info(['',
                               'Use "{0}" to restart working on #{1} with a clean copy of the master branch.'],
                               self._format_command("checkout", ticket=ticket, base=MASTER_BRANCH), ticket)

    def gather(self, branch, *tickets_or_branches):
        r"""
        Create a new branch ``branch`` with ``tickets_or_remote_branches``
        applied.

        This method is not wrapped in the commandline dev scripts. It
        does nothing that cannot be done with ``checkout`` and
        ``merge``, it just steepens the learning curve by introducing
        yet another command. Unless a clear use case emerges, it
        should be removed.

        INPUT:

        - ``branch`` -- a string, the name of the new branch

        - ``tickets_or_branches`` -- a list of integers and strings; for an
          integer or string identifying a ticket, the branch on the trac ticket
          gets merged, for the name of a local or remote branch, that branch
          gets merged.

        .. SEEALSO::

            - :meth:`merge` -- merge into the current branch rather
              than creating a new one

        TESTS:

        Create a doctest setup with a single user::

            sage: from sage.dev.test.sagedev import single_user_setup
            sage: dev, config, UI, server = single_user_setup()

        Create tickets and branches::

            sage: dev._UI.append("Summary: summary1\ndescription")
            sage: dev.create_ticket()
            Created ticket #1 at https://trac.sagemath.org/1.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=1" to create a new local branch)
            1
            sage: dev.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.
            sage: open("tracked","w").close()
            sage: dev.git.silent.add("tracked")
            sage: dev.git.super_silent.commit(message="added tracked")
            sage: dev._UI.append("y")
            sage: dev._UI.append("y")
            sage: dev.push()
            The branch "u/doctest/ticket/1" does not exist on the remote server.
            Create new remote branch? [Yes/no] y

        Gather all these branches::

            sage: dev._sagedev.gather("gather_branch", "#1", "ticket/1", "u/doctest/ticket/1")
        """
        try:
            self.reset_to_clean_state()
            self.clean()
        except OperationCancelledError:
            self._UI.error("Cannot gather branches because working directory is not in a clean state.")
            raise OperationCancelledError("working directory not clean")

        self._check_local_branch_name(branch, exists=False)

        branches = []
        for ticket_or_branch in tickets_or_branches:
            local_branch = None
            remote_branch = None
            if self._is_ticket_name(ticket_or_branch):
                ticket = self._ticket_from_ticket_name(ticket_or_branch)
                remote_branch = self.trac._branch_for_ticket(ticket)
                if remote_branch is None:
                    raise SageDevValueError("Ticket #{0} does not have a branch set yet.".format(ticket))
            elif self._is_local_branch_name(ticket_or_branch, exists=True):
                local_branch = ticket_or_branch
            else:
                remote_branch = ticket_or_branch

            if local_branch:
                self._check_local_branch_name(local_branch, exists=True)
                branches.append(("local",local_branch))
            if remote_branch:
                self._check_remote_branch_name(remote_branch, exists=True)
                branches.append(("remote",remote_branch))

        self._UI.debug('Creating a new branch "{0}".'.format(branch))
        self.git.super_silent.branch(branch, MASTER_BRANCH)
        self.git.super_silent.checkout(branch)

        try:
            for local_remote,branch_name in branches:
                self._UI.debug('Merging {2} branch "{0}" into "{1}".'
                              .format(branch_name, branch, local_remote))
                self.merge(branch, pull=local_remote=="remote")
        except:
            self.git.reset_to_clean_state()
            self.git.clean_wrapper()
            self.vanilla()
            self.git.super_silent.branch("-D", branch)
            self._UI.debug('Deleted branch "{0}".'.format(branch))

    def merge(self, ticket_or_branch=MASTER_BRANCH, pull=None, create_dependency=None):
        r"""
        Merge changes from ``ticket_or_branch`` into the current branch.

        Incorporate commits from other tickets/branches into the
        current branch.

        Optionally, you can add the merged ticket to the trac
        "Dependency:" field. Note that the merged commits become part
        of the current branch, regardless of whether they are noted on
        trac. Adding a dependency implies the following:

        - the other ticket must be positively reviewed and merged
          before this ticket may be merged into the official release
          of sage.  The commits included from a dependency don't need
          to be reviewed in this ticket, whereas commits reviewed in
          this ticket from a non-dependency may make reviewing the
          other ticket easier.

        - you can more easily merge in future changes to dependencies.
          So if you need a feature from another ticket it may be
          appropriate to create a dependency to that you may more
          easily benefit from others' work on that ticket.

        - if you depend on another ticket then you need to worry about
          the progress on that ticket.  If that ticket is still being
          actively developed then you may need to make further merges
          in the future if conflicts arise.

        INPUT:

        - ``ticket_or_branch`` -- an integer or strings (default:
          ``'master'``); for an integer or string identifying a ticket, the
          branch on the trac ticket gets merged (or the local branch for the
          ticket, if ``pull`` is ``False``), for the name of a local or
          remote branch, that branch gets merged. If ``'dependencies'``, the
          dependencies are merged in one by one.

        - ``pull`` -- a boolean or ``None`` (default: ``None``); if
          ``ticket_or_branch`` identifies a ticket, whether to pull the
          latest branch on the trac ticket (the default); if
          ``ticket_or_branch`` is a branch name, then ``pull`` controls
          whether it should be interpreted as a remote branch (``True``) or as
          a local branch (``False``). If it is set to ``None``, then it will
          take ``ticket_or_branch`` as a remote branch if it exists, and as a
          local branch otherwise.

        - ``create_dependency`` -- a boolean or ``None`` (default: ``None``),
          whether to create a dependency to ``ticket_or_branch``. If ``None``,
          then a dependency is created if ``ticket_or_branch`` identifies a
          ticket and if the current branch is associated to a ticket.

        .. NOTE::

            Dependencies are stored locally and only updated with respect to
            the remote server during :meth:`push` and :meth:`pull`.

        .. SEEALSO::

            - :meth:`show_dependencies` -- see the current
              dependencies.

            - :meth:`GitInterface.merge` -- git's merge command has
              more options and can merge multiple branches at once.

        TESTS:

        Create a doctest setup with two users::

            sage: from sage.dev.test.sagedev import two_user_setup
            sage: alice, config_alice, bob, config_bob, server = two_user_setup()

        Create tickets and branches::

            sage: alice._chdir()
            sage: alice._UI.append("Summary: summary1\ndescription")
            sage: alice.create_ticket()
            Created ticket #1 at https://trac.sagemath.org/1.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=1" to create a new local branch)
            1
            sage: alice._UI.append("Summary: summary2\ndescription")
            sage: alice.create_ticket()
            Created ticket #2 at https://trac.sagemath.org/2.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=2" to create a new local branch)
            2

        Alice creates two branches and merges them::

            sage: alice.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.
            sage: open("alice1","w").close()
            sage: alice.git.silent.add("alice1")
            sage: alice.git.super_silent.commit(message="added alice1")
            sage: alice.checkout(ticket=2)
            On ticket #2 with associated local branch "ticket/2".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.
            sage: with open("alice2","w") as f: f.write("alice")
            sage: alice.git.silent.add("alice2")
            sage: alice.git.super_silent.commit(message="added alice2")

        When merging for a ticket, the branch on the trac ticket matters::

            sage: alice.merge("#1")
            Cannot merge remote branch for #1 because no branch has been set on the trac
            ticket.
            sage: alice.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.
            sage: alice._UI.append("y")
            sage: alice.push()
            The branch "u/alice/ticket/1" does not exist on the remote server.
            Create new remote branch? [Yes/no] y
            sage: alice.checkout(ticket=2)
            On ticket #2 with associated local branch "ticket/2".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.
            sage: alice.merge("#1", pull=False)
            Merging the local branch "ticket/1" into the local branch "ticket/2".
            Automatic merge successful.
            <BLANKLINE>
            #  (use "sage --dev commit" to commit your merge)
            <BLANKLINE>
            Added dependency on #1 to #2.

        Check that merging dependencies works::

            sage: alice.merge("dependencies")
            Merging the remote branch "u/alice/ticket/1" into the local branch "ticket/2".
            Automatic merge successful.
            <BLANKLINE>
            #  (use "sage --dev commit" to commit your merge)

        Merging local branches::

            sage: alice.merge("ticket/1")
            Merging the local branch "ticket/1" into the local branch "ticket/2".
            Automatic merge successful.
            <BLANKLINE>
            #  (use "sage --dev commit" to commit your merge)

        A remote branch for a local branch is only merged in if ``pull`` is set::

            sage: alice._sagedev._set_remote_branch_for_branch("ticket/1", "nonexistant")
            sage: alice.merge("ticket/1")
            Merging the local branch "ticket/1" into the local branch "ticket/2".
            Automatic merge successful.
            <BLANKLINE>
            #  (use "sage --dev commit" to commit your merge)
            sage: alice.merge("ticket/1", pull=True)
            Branch "ticket/1" does not exist on the remote system.

        Bob creates a conflicting commit::

            sage: bob._chdir()
            sage: bob.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.
            sage: with open("alice2","w") as f: f.write("bob")
            sage: bob.git.silent.add("alice2")
            sage: bob.git.super_silent.commit(message="added alice2")
            sage: bob._UI.append("y")
            sage: bob._UI.append("y")
            sage: bob.push()
            The branch "u/bob/ticket/1" does not exist on the remote server.
            Create new remote branch? [Yes/no] y
            The branch field of ticket #1 needs to be updated from its current value
            "u/alice/ticket/1" to "u/bob/ticket/1"
            Change the "Branch:" field? [Yes/no] y

        The merge now requires manual conflict resolution::

            sage: alice._chdir()
            sage: alice._UI.append("abort")
            sage: alice.merge("#1")
            Merging the remote branch "u/bob/ticket/1" into the local branch "ticket/2".
            Automatic merge failed, there are conflicting commits.
            <BLANKLINE>
            Auto-merging alice2
            CONFLICT (add/add): Merge conflict in alice2
            <BLANKLINE>
            Please edit the affected files to resolve the conflicts. When you are finished,
            your resolution will be commited.
            Finished? [ok/Abort] abort
            sage: alice._UI.append("ok")
            sage: alice.merge("#1")
            Merging the remote branch "u/bob/ticket/1" into the local branch "ticket/2".
            Automatic merge failed, there are conflicting commits.
            <BLANKLINE>
            Auto-merging alice2
            CONFLICT (add/add): Merge conflict in alice2
            <BLANKLINE>
            Please edit the affected files to resolve the conflicts. When you are finished,
            your resolution will be commited.
            Finished? [ok/Abort] ok
            Created a commit from your conflict resolution.

        We cannot merge a ticket into itself::

            sage: alice.merge(2)
            cannot merge a ticket into itself

        We also cannot merge if the working directory has uncommited changes::

            sage: alice._UI.append("cancel")
            sage: with open("alice2","w") as f: f.write("uncommited change")
            sage: alice.merge(1)
            The following files in your working directory contain uncommitted changes:
            <BLANKLINE>
                 alice2
            <BLANKLINE>
            Discard changes? [discard/Cancel/stash] cancel
            Cannot merge because working directory is not in a clean state.
            <BLANKLINE>
            #  (use "sage --dev commit" to commit your changes)
        """
        try:
            self.reset_to_clean_state()
            self.clean()
        except OperationCancelledError:
            self._UI.error("Cannot merge because working directory is not in a clean state.")
            self._UI.info(['', '(use "{0}" to commit your changes)'],
                          self._format_command('commit'))
            raise OperationCancelledError("working directory not clean")

        from git_error import DetachedHeadError
        try:
            current_branch = self.git.current_branch()
        except DetachedHeadError:
            self._UI.error('Not on any branch.')
            self._UI.info(['', '(use "{0}" to checkout a branch)'],
                           self._format_command("checkout"))
            raise OperationCancelledError("detached head")

        current_ticket = self._current_ticket()

        ticket = None
        branch = None
        remote_branch = None

        if ticket_or_branch == 'dependencies':
            if current_ticket == None:
                raise SageDevValueError("dependencies can only be merged if currently on a ticket.")
            if pull == False:
                raise SageDevValueError('"pull" must not be "False" when merging dependencies.')
            if create_dependency != None:
                raise SageDevValueError('"create_dependency" must not be set when merging dependencies.')
            for dependency in self._dependencies_for_ticket(current_ticket):
                self._UI.debug("Merging dependency #{0}.".format(dependency))
                self.merge(ticket_or_branch=dependency, pull=True)
            return
        elif self._is_ticket_name(ticket_or_branch):
            ticket = self._ticket_from_ticket_name(ticket_or_branch)
            if ticket == current_ticket:
                raise SageDevValueError("cannot merge a ticket into itself")
            self._check_ticket_name(ticket, exists=True)
            if pull is None:
                pull = True
            if create_dependency is None:
                create_dependency = True
            if self._has_local_branch_for_ticket(ticket):
                branch = self._local_branch_for_ticket(ticket)
            if pull:
                remote_branch = self.trac._branch_for_ticket(ticket)
                if remote_branch is None:
                    self._UI.error("Cannot merge remote branch for #{0} because no branch has"
                                   " been set on the trac ticket.", ticket)
                    raise OperationCancelledError("remote branch not set on trac")
        elif pull == False or (pull is None and not
                               self._is_remote_branch_name(ticket_or_branch, exists=True)):
            # ticket_or_branch should be interpreted as a local branch name
            branch = ticket_or_branch
            self._check_local_branch_name(branch, exists=True)
            pull = False
            if create_dependency == True:
                if self._has_ticket_for_local_branch(branch):
                    ticket = self._ticket_for_local_branch(branch)
                else:
                    raise SageDevValueError('"create_dependency" must not be "True" if'
                                            ' "ticket_or_branch" is a local branch which'
                                            ' is not associated to a ticket.')
            else:
                create_dependency = False
        else:
            # ticket_or_branch should be interpreted as a remote branch name
            remote_branch = ticket_or_branch
            self._check_remote_branch_name(remote_branch, exists=True)
            pull = True
            if create_dependency == True:
                raise SageDevValueError('"create_dependency" must not be "True" if'
                                        ' "ticket_or_branch" is a local branch.')
            create_dependency = False

        if pull:
            assert remote_branch
            if not self._is_remote_branch_name(remote_branch, exists=True):
                self._UI.error('Can not merge remote branch "{0}". It does not exist.',
                               remote_branch)
                raise OperationCancelledError("no such branch")
            self._UI.show('Merging the remote branch "{0}" into the local branch "{1}".',
                          remote_branch, current_branch)
            self.git.super_silent.fetch(self.git._repository_anonymous, remote_branch)
            local_merge_branch = 'FETCH_HEAD'
        else:
            assert branch
            self._UI.show('Merging the local branch "{0}" into the local branch "{1}".',
                          branch, current_branch)
            local_merge_branch = branch

        from git_error import GitError
        try:
            self.git.super_silent.merge(local_merge_branch)
            self._UI.show('Automatic merge successful.')
            self._UI.info(['', '(use "{0}" to commit your merge)'],
                          self._format_command('commit'))
        except GitError as e:
            try:
                self._UI.show('Automatic merge failed, there are conflicting commits.')
                excluded = ['Aborting',
                    "Automatic merge failed; fix conflicts and then commit the result."]
                lines = e.stdout.splitlines() + e.stderr.splitlines()
                lines = [line for line in lines if line not in excluded]
                self._UI.show([''] + lines + [''])
                self._UI.show('Please edit the affected files to resolve the conflicts.'
                              ' When you are finished, your resolution will be commited.')
                sel = self._UI.select("Finished?", ['ok', 'abort'], default=1)
                if sel == 'ok':
                    self.git.silent.commit(a=True, no_edit=True)
                    self._UI.show("Created a commit from your conflict resolution.")
                elif sel == 'abort':
                    raise OperationCancelledError("user requested")
                else:
                    assert False
            except Exception as e:
                self.git.reset_to_clean_state()
                self.git.clean_wrapper()
                raise

        if create_dependency:
            assert ticket and current_ticket
            dependencies = list(self._dependencies_for_ticket(current_ticket))
            if ticket in dependencies:
                self._UI.debug("Not recording dependency on #{0} because #{1} already depends on #{0}.",
                               ticket, current_ticket)
            else:
                self._UI.show(['', "Added dependency on #{0} to #{1}."], ticket, current_ticket)
                self._set_dependencies_for_ticket(current_ticket, dependencies+[ticket])

    def tickets(self, include_abandoned=False, cached=True):
        r"""
        Print the tickets currently being worked on in your local
        repository.

        This function shows the branch names as well as the ticket numbers for
        all active tickets.  It also shows local branches that are not
        associated to ticket numbers.

        INPUT:

        - ``include_abandoned`` -- boolean (default: ``False``), whether to
          include abandoned branches.

        - ``cached`` -- boolean (default: ``True``), whether to try to pull the
          summaries from the ticket cache; if ``True``, then the summaries
          might not be accurate if they changed since they were last updated.
          To update the summaries, set this to ``False``.

        .. SEEALSO::

            - :meth:`abandon_ticket` -- hide tickets from this method.

            - :meth:`remote_status` -- also show status compared to
              the trac server.

            - :meth:`current_ticket` -- get the current ticket.

        TESTS:

        Create a doctest setup with a single user::

            sage: from sage.dev.test.sagedev import single_user_setup
            sage: dev, config, UI, server = single_user_setup()

        Create some tickets::

            sage: dev.tickets()
            * : master

            sage: UI.append("Summary: summary\ndescription")
            sage: dev.create_ticket()
            Created ticket #1 at https://trac.sagemath.org/1.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=1" to create a new local branch)
            1
            sage: dev.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.
            sage: UI.append("Summary: summary\ndescription")
            sage: dev.create_ticket()
            Created ticket #2 at https://trac.sagemath.org/2.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=2" to create a new local branch)
            2
            sage: dev.checkout(ticket=2)
            On ticket #2 with associated local branch "ticket/2".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.
            sage: dev.tickets()
                : master
              #1: ticket/1 summary
            * #2: ticket/2 summary
        """
        branches = self.git.local_branches()
        from git_error import DetachedHeadError
        try:
            current_branch = self.git.current_branch()
        except DetachedHeadError:
            current_branch = None
        branches = [ branch for branch in branches if include_abandoned or not self._is_trash_name(branch) ]
        if not branches:
            return
        ret = []
        for branch in branches:
            ticket = None
            ticket_summary = ""
            extra = " "
            if self._has_ticket_for_local_branch(branch):
                ticket = self._ticket_for_local_branch(branch)
                try:
                    try:
                        ticket_summary = self.trac._get_attributes(ticket, cached=cached)['summary']
                    except KeyError:
                        ticket_summary = self.trac._get_attributes(ticket, cached=False)['summary']
                except TracConnectionError:
                    ticket_summary = ""
            if current_branch == branch:
                extra = "*"
            ticket_str = "#"+str(ticket) if ticket else ""
            ret.append(("{0:>7}: {1} {2}".format(ticket_str, branch, ticket_summary), extra))
        while all([info.startswith(' ') for (info, extra) in ret]):
            ret = [(info[1:],extra) for (info, extra) in ret]
        ret = sorted(ret)
        ret = ["{0} {1}".format(extra,info) for (info,extra) in ret]
        self._UI.show("\n".join(ret))

    def vanilla(self, release=MASTER_BRANCH):
        r"""
        Return to a clean version of Sage.

        INPUT:

        - ``release`` -- a string or decimal giving the release name (default:
          ``'master'``).  In fact, any tag, commit or branch will work.  If the
          tag does not exist locally an attempt to fetch it from the server
          will be made.

        Git equivalent::

            Checks out a given tag, commit or branch in detached head mode.

        .. SEEALSO::

            - :meth:`checkout` -- checkout another branch, ready to
              develop on it.

            - :meth:`pull` -- pull a branch from the server and merge
              it.

        TESTS:

        Create a doctest setup with a single user::

            sage: from sage.dev.test.sagedev import single_user_setup
            sage: dev, config, UI, server = single_user_setup()

        Go to a sage release::

            sage: dev.git.current_branch()
            'master'
            sage: dev.vanilla()
            sage: dev.git.current_branch()
            Traceback (most recent call last):
            ...
            DetachedHeadError: unexpectedly, git is in a detached HEAD state
        """
        if hasattr(release, 'literal'):
            release = release.literal
        release = str(release)

        try:
            self.reset_to_clean_state()
            self.clean()
        except OperationCancelledError:
            self._UI.error("Cannot checkout a release while your working directory is not clean.")
            raise OperationCancelledError("working directory not clean")

        # we do not do any checking on the argument here, trying to be liberal
        # about what are valid inputs
        try:
            self.git.super_silent.checkout(release, detach=True)
        except GitError as e:
            try:
                self.git.super_silent.fetch(self.git._repository_anonymous, release)
            except GitError as e:
                self._UI.error('"{0}" does not exist locally or on the remote server.'.format(release))
                raise OperationCancelledError("no such tag/branch/...")
            self.git.super_silent.checkout('FETCH_HEAD', detach=True)

    def diff(self, base='commit'):
        r"""
        Show how the current file system differs from ``base``.

        INPUT:

        - ``base`` -- a string; show the differences against the latest
          ``'commit'`` (the default), against the branch ``'master'`` (or any
          other branch name), or the merge of the ``'dependencies'`` of the
          current ticket (if the dependencies merge cleanly)

        .. SEEALSO::

            - :meth:`commit` -- record changes into the repository.

            - :meth:`tickets` -- list local tickets (you may
              want to commit your changes to a branch other than the
              current one).

        TESTS:

        Create a doctest setup with a single user::

            sage: from sage.dev.test.sagedev import single_user_setup
            sage: dev, config, UI, server = single_user_setup()

        Create some tickets and make one depend on the others::

            sage: UI.append("Summary: summary\ndescription")
            sage: dev.create_ticket()
            Created ticket #1 at https://trac.sagemath.org/1.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=1" to create a new local branch)
            1
            sage: dev.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.
            sage: UI.append("y")
            sage: dev.push()
            The branch "u/doctest/ticket/1" does not exist on the remote server.
            Create new remote branch? [Yes/no] y
            sage: UI.append("Summary: summary\ndescription")
            sage: dev.create_ticket()
            Created ticket #2 at https://trac.sagemath.org/2.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=2" to create a new local branch)
            2
            sage: dev.checkout(ticket=2)
            On ticket #2 with associated local branch "ticket/2".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.
            sage: UI.append("y")
            sage: dev.push()
            The branch "u/doctest/ticket/2" does not exist on the remote server.
            Create new remote branch? [Yes/no] y
            sage: UI.append("Summary: summary\ndescription")
            sage: dev.create_ticket()
            Created ticket #3 at https://trac.sagemath.org/3.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=3" to create a new local branch)
            3
            sage: dev.checkout(ticket=3)
            On ticket #3 with associated local branch "ticket/3".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.
            sage: UI.append("y")
            sage: dev.push()
            The branch "u/doctest/ticket/3" does not exist on the remote server.
            Create new remote branch? [Yes/no] y
            sage: dev.merge("#1")
            Merging the remote branch "u/doctest/ticket/1" into the local branch "ticket/3".
            Automatic merge successful.
            <BLANKLINE>
            #  (use "sage --dev commit" to commit your merge)
            <BLANKLINE>
            Added dependency on #1 to #3.
            sage: dev.merge("#2")
            Merging the remote branch "u/doctest/ticket/2" into the local branch "ticket/3".
            Automatic merge successful.
            <BLANKLINE>
            #  (use "sage --dev commit" to commit your merge)
            <BLANKLINE>
            Added dependency on #2 to #3.

        Make some non-conflicting changes on the tickets::

            sage: dev.checkout(ticket="#1")
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.
            sage: with open("ticket1","w") as f: f.write("ticket1")
            sage: dev.git.silent.add("ticket1")
            sage: dev.git.super_silent.commit(message="added ticket1")

            sage: dev.checkout(ticket="#2")
            On ticket #2 with associated local branch "ticket/2".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.
            sage: with open("ticket2","w") as f: f.write("ticket2")
            sage: dev.git.silent.add("ticket2")
            sage: dev.git.super_silent.commit(message="added ticket2")
            sage: UI.append("y")
            sage: dev.push()
            Local commits that are not on the remote branch "u/doctest/ticket/2":
            <BLANKLINE>
                ...: added ticket2
            <BLANKLINE>
            Push to remote branch? [Yes/no] y

            sage: dev.checkout(ticket="#3")
            On ticket #3 with associated local branch "ticket/3".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.
            sage: open("ticket3","w").close()
            sage: dev.git.silent.add("ticket3")
            sage: dev.git.super_silent.commit(message="added ticket3")
            sage: UI.append("y")
            sage: dev.push()
            Local commits that are not on the remote branch "u/doctest/ticket/3":
            <BLANKLINE>
                ...: added ticket3
            <BLANKLINE>
            Push to remote branch? [Yes/no] y
            Uploading your dependencies for ticket #3: "" => "#1, #2"

        A diff against the previous commit::

            sage: dev.diff()

        A diff against a ticket will always take the branch on trac::

            sage: dev.diff("#1")
            diff --git a/ticket3 b/ticket3
            new file mode ...
            index ...
            sage: dev.diff("ticket/1")
            diff --git a/ticket1 b/ticket1
            deleted file mode ...
            index ...
            diff --git a/ticket3 b/ticket3
            new file mode ...
            index ...
            sage: dev.checkout(ticket="#1")
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.
            sage: UI.append("y")
            sage: dev.push()
            Local commits that are not on the remote branch "u/doctest/ticket/1":
            <BLANKLINE>
                ...: added ticket1
            <BLANKLINE>
            Push to remote branch? [Yes/no] y
            sage: dev.checkout(ticket="#3")
            On ticket #3 with associated local branch "ticket/3".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.
            sage: dev.diff("#1")
            diff --git a/ticket1 b/ticket1
            deleted file mode ...
            index ...
            diff --git a/ticket3 b/ticket3
            new file mode ...
            index ...

        A diff against the dependencies::

            sage: dev.diff("dependencies")
            Dependency #1 has not been merged into "ticket/3" (at least not its latest
            version).
            #  (use "sage --dev merge --ticket=1" to merge it)
            <BLANKLINE>
            Dependency #2 has not been merged into "ticket/3" (at least not its latest
            version).
            #  (use "sage --dev merge --ticket=2" to merge it)
            <BLANKLINE>
            diff --git a/ticket1 b/ticket1
            deleted file mode ...
            index ...
            diff --git a/ticket2 b/ticket2
            deleted file mode ...
            index ...
            diff --git a/ticket3 b/ticket3
            new file mode ...
            index ...
            sage: dev.merge("#1")
            Merging the remote branch "u/doctest/ticket/1" into the local branch "ticket/3".
            Automatic merge successful.
            <BLANKLINE>
            #  (use "sage --dev commit" to commit your merge)
            sage: dev.merge("#2")
            Merging the remote branch "u/doctest/ticket/2" into the local branch "ticket/3".
            Automatic merge successful.
            <BLANKLINE>
            #  (use "sage --dev commit" to commit your merge)
            sage: dev.diff("dependencies")
            diff --git a/ticket3 b/ticket3
            new file mode ...
            index ...

        This does not work if the dependencies do not merge::

            sage: dev.checkout(ticket="#1")
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.
            sage: with open("ticket2","w") as f: f.write("foo")
            sage: dev.git.silent.add("ticket2")
            sage: dev.git.super_silent.commit(message="added ticket2")
            sage: UI.append("y")
            sage: dev.push()
            Local commits that are not on the remote branch "u/doctest/ticket/1":
            <BLANKLINE>
                ...: added ticket2
            <BLANKLINE>
            Push to remote branch? [Yes/no] y

            sage: dev.checkout(ticket="#3")
            On ticket #3 with associated local branch "ticket/3".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.
            sage: dev.diff("dependencies")
            Dependency #1 has not been merged into "ticket/3" (at least not its latest
            version).
            #  (use "sage --dev merge --ticket=1" to merge it)
            <BLANKLINE>
            Dependency #2 does not merge cleanly with the other dependencies. Your diff
            could not be computed.
        """
        if base == "dependencies":
            current_ticket = self._current_ticket()
            if current_ticket is None:
                raise SageDevValueError("'dependencies' are only supported if currently on a ticket.")
            try:
                self.reset_to_clean_state()
                self.clean()
            except OperationCancelledError:
                self._UI.error("Cannot create merge of dependencies because working directory is not clean.")
                raise

            self._is_master_uptodate(action_if_not="warning")

            branch = self.git.current_branch()
            merge_base = self.git.merge_base(branch, MASTER_BRANCH).splitlines()[0]
            temporary_branch = self._new_local_branch_for_trash("diff")
            self.git.super_silent.branch(temporary_branch, merge_base)
            try:
                self.git.super_silent.checkout(temporary_branch)
                try:
                    self._UI.debug("Merging dependencies of #{0}.".format(current_ticket))
                    for dependency in self._dependencies_for_ticket(current_ticket):
                        self._check_ticket_name(dependency, exists=True)
                        remote_branch = self.trac._branch_for_ticket(dependency)
                        if remote_branch is None:
                            self._UI.warning("Dependency #{0} has no branch field set.".format(dependency))
                        self._check_remote_branch_name(remote_branch, exists=True)
                        self.git.super_silent.fetch(self.git._repository_anonymous, remote_branch)
                        merge_base_dependency = self.git.merge_base(MASTER_BRANCH, 'FETCH_HEAD').splitlines()[0]
                        if merge_base_dependency != merge_base and \
                           self.git.is_child_of(merge_base_dependency, merge_base):
                            self._UI.warning('The remote branch "{0}" is based on a later version of sage'
                                             ' compared to the local branch "{1}". The diff might therefore'
                                             ' contain unrelated changes.')
                            self._UI.info(['Use "{2}" to merge latest version of Sage into your branch.', ''],
                                          remote_branch, branch, self._format_command("merge"))
                        if self.git.is_child_of(merge_base, 'FETCH_HEAD'):
                            self._UI.debug('Dependency #{0} has already been merged into the master'
                                           ' branch of your version of sage.', dependency)
                        else:
                            if not self.git.is_child_of(branch, 'FETCH_HEAD'):
                                self._UI.warning('Dependency #{0} has not been merged into "{1}" (at'
                                                 ' least not its latest version).', dependency, branch)
                                self._UI.info(['(use "{0}" to merge it)', ''],
                                              self._format_command("merge", ticket_or_branch=str(dependency)))
                            from git_error import GitError
                            try:
                                self.git.super_silent.merge('FETCH_HEAD')
                            except GitError as e:
                                self._UI.error("Dependency #{0} does not merge cleanly with the other"
                                               " dependencies. Your diff could not be computed.", dependency)
                                raise OperationCancelledError("merge failed")

                    self.git.echo.diff("{0}..{1}".format(temporary_branch, branch))
                    return
                finally:
                    self.git.reset_to_clean_state()
                    self.git.clean_wrapper()
                    self.git.super_silent.checkout(branch)
            finally:
                self.git.super_silent.branch("-D", temporary_branch)

        if base == "commit":
            base = "HEAD"
        else:
            if self._is_ticket_name(base):
                ticket = self._ticket_from_ticket_name(base)
                self._check_ticket_name(ticket, exists=True)
                base = self.trac._branch_for_ticket(ticket)
                if base is None:
                    self._UI.error("Ticket #{0} has no branch set on trac.".format(ticket))

            if self._is_local_branch_name(base, exists=True):
                pass
            else:
                self._check_remote_branch_name(base, exists=True)
                self._is_master_uptodate(action_if_not="warning")
                self.git.super_silent.fetch(self.git._repository_anonymous, base)
                base = 'FETCH_HEAD'

        self.git.echo.diff(base)

    def show_dependencies(self, ticket=None, all=False, _seen=None): # all = recursive
        r"""
        Show the dependencies of ``ticket``.

        INPUT:

        - ``ticket`` -- a string or integer identifying a ticket or ``None``
          (default: ``None``), the ticket for which dependencies are displayed.
          If ``None``, then the dependencies for the current ticket are
          displayed.

        - ``all`` -- boolean (default: ``True``), whether to recursively list
          all tickets on which this ticket depends (in depth-first order), only
          including tickets that have a local branch.

        .. NOTE::

            Ticket dependencies are stored locally and only updated with
            respect to the remote server during :meth:`push` and
            :meth:`pull`.

        .. SEEALSO::

            - :meth:`TracInterface.dependencies` -- Query Trac to find
              dependencies.

            - :meth:`remote_status` -- will show the status of tickets
              with respect to the remote server.

            - :meth:`merge` -- Merge in changes from a dependency.

            - :meth:`diff` -- Show the changes in this branch over the
              dependencies.

        TESTS:

        Create a doctest setup with a single user::

            sage: from sage.dev.test.sagedev import single_user_setup
            sage: dev, config, UI, server = single_user_setup()

        Create some tickets and add dependencies::

            sage: UI.append("Summary: summary\ndescription")
            sage: dev.create_ticket()
            Created ticket #1 at https://trac.sagemath.org/1.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=1" to create a new local branch)
            1
            sage: dev.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.
            sage: UI.append("Summary: summary\ndescription")
            sage: dev.create_ticket()
            Created ticket #2 at https://trac.sagemath.org/2.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=2" to create a new local branch)
            2
            sage: dev.checkout(ticket=2)
            On ticket #2 with associated local branch "ticket/2".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.
            sage: UI.append("Summary: summary\ndescription")
            sage: dev.create_ticket()
            Created ticket #3 at https://trac.sagemath.org/3.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=3" to create a new local branch)
            3
            sage: dev.checkout(ticket=3)
            On ticket #3 with associated local branch "ticket/3".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.
            sage: UI.append("Summary: summary\ndescription")
            sage: dev.create_ticket()
            Created ticket #4 at https://trac.sagemath.org/4.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=4" to create a new local branch)
            4
            sage: dev.checkout(ticket=4)
            On ticket #4 with associated local branch "ticket/4".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.

            sage: dev.merge('ticket/2',create_dependency=True)
            Merging the local branch "ticket/2" into the local branch "ticket/4".
            Automatic merge successful.
            <BLANKLINE>
            #  (use "sage --dev commit" to commit your merge)
            <BLANKLINE>
            Added dependency on #2 to #4.
            sage: dev.merge('ticket/3',create_dependency=True)
            Merging the local branch "ticket/3" into the local branch "ticket/4".
            Automatic merge successful.
            <BLANKLINE>
            #  (use "sage --dev commit" to commit your merge)
            <BLANKLINE>
            Added dependency on #3 to #4.
            sage: dev.checkout(ticket='#2')
            On ticket #2 with associated local branch "ticket/2".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.
            sage: dev.merge('ticket/1', create_dependency=True)
            Merging the local branch "ticket/1" into the local branch "ticket/2".
            Automatic merge successful.
            <BLANKLINE>
            #  (use "sage --dev commit" to commit your merge)
            <BLANKLINE>
            Added dependency on #1 to #2.
            sage: dev.checkout(ticket='#3')
            On ticket #3 with associated local branch "ticket/3".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.
            sage: dev.merge('ticket/1', create_dependency=True)
            Merging the local branch "ticket/1" into the local branch "ticket/3".
            Automatic merge successful.
            <BLANKLINE>
            #  (use "sage --dev commit" to commit your merge)
            <BLANKLINE>
            Added dependency on #1 to #3.

        Check that the dependencies show correctly::

            sage: dev.checkout(ticket='#4')
            On ticket #4 with associated local branch "ticket/4".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.
            sage: dev.show_dependencies()
            Ticket #4 depends on #2, #3.
            sage: dev.show_dependencies('#4')
            Ticket #4 depends on #2, #3.
            sage: dev.show_dependencies('#3')
            Ticket #3 depends on #1.
            sage: dev.show_dependencies('#2')
            Ticket #2 depends on #1.
            sage: dev.show_dependencies('#1')
            Ticket #1 has no dependencies.
            sage: dev.show_dependencies('#4', all=True)
            Ticket #4 depends on #3, #1, #2.
        """
        if ticket is None:
            ticket = self._current_ticket()

        if ticket is None:
            raise SageDevValueError("ticket must be specified")

        self._check_ticket_name(ticket)
        ticket = self._ticket_from_ticket_name(ticket)

        if not self._has_local_branch_for_ticket(ticket):
            raise SageDevValueError('ticket must be a ticket with a local branch. Use "{0}" to checkout the ticket first.'.format(self._format_command("checkout",ticket=ticket)))

        branch = self._local_branch_for_ticket(ticket)
        if all:
            ret = []
            stack = [ticket]
            while stack:
                t = stack.pop()
                if t in ret: continue
                ret.append(t)
                if not self._has_local_branch_for_ticket(t):
                    self._UI.warning("no local branch for ticket #{0} present, some dependencies might be missing in the output.".format(t))
                    continue
                deps = self._dependencies_for_ticket(t)
                for d in deps:
                    if d not in stack and d not in ret:
                        stack.append(d)
            ret = ret[1:]
        else:
            ret = self._dependencies_for_ticket(ticket)

        if ret:
            self._UI.show("Ticket #{0} depends on {1}.".format(ticket,", ".join(["#{0}".format(d) for d in ret])))
        else:
            self._UI.show("Ticket #{0} has no dependencies.".format(ticket))

    def upload_ssh_key(self, public_key=None):
        r"""
        Upload ``public_key`` to gitolite through the trac interface.

        INPUT:

        - ``public_key`` -- a string or ``None`` (default: ``None``), the path
          of the key file, defaults to ``~/.ssh/id_rsa.pub`` (or
          ``~/.ssh/id_dsa.pub`` if it exists).

        TESTS:

        Create a doctest setup with a single user::

            sage: from sage.dev.test.sagedev import single_user_setup
            sage: dev, config, UI, server = single_user_setup()

        Create and upload a key file::

            sage: import os
            sage: public_key = os.path.join(dev._sagedev.tmp_dir, "id_rsa.pub")
            sage: UI.append("no")
            sage: UI.append("yes")
            sage: dev.upload_ssh_key(public_key=public_key)
            The trac git server requires your SSH public key to be able to identify you.
            Upload ".../id_rsa.pub" to trac? [Yes/no] yes
            File not found: ".../id_rsa.pub"
            Create new ssh key pair? [Yes/no] no
            <BLANKLINE>
            #  Use "sage --dev upload-ssh-key" to upload a public key. Or set your key manually at https://trac.sagemath.org/prefs/sshkeys.
            sage: UI.append("yes")
            sage: UI.append("yes")
            sage: dev.upload_ssh_key(public_key=public_key)
            The trac git server requires your SSH public key to be able to identify you.
            Upload ".../id_rsa.pub" to trac? [Yes/no] yes
            File not found: ".../id_rsa.pub"
            Create new ssh key pair? [Yes/no] yes
            Generating ssh key.
            Your key has been uploaded.
            sage: UI.append("yes")
            sage: dev.upload_ssh_key(public_key=public_key)
            The trac git server requires your SSH public key to be able to identify you.
            Upload ".../id_rsa.pub" to trac? [Yes/no] yes
            Your key has been uploaded.
        """
        try:
            import os
            if public_key is None:
                public_key = os.path.expanduser("~/.ssh/id_dsa.pub")
                if not os.path.exists(public_key):
                    public_key = os.path.expanduser("~/.ssh/id_rsa.pub")
            if not public_key.endswith(".pub"):
                raise SageDevValueError('public key must end with ".pub".')

            self._UI.show('The trac git server requires your SSH public key'
                          ' to be able to identify you.')
            if not self._UI.confirm('Upload "{0}" to trac?'
                                    .format(public_key), default=True):
                raise OperationCancelledError("do not upload key")

            if not os.path.exists(public_key):
                self._UI.warning('File not found: "{0}"'.format(public_key))
                if not self._UI.confirm('Create new ssh key pair?', default=True):
                    raise OperationCancelledError("no keyfile found")
                private_key = public_key[:-4]
                self._UI.show("Generating ssh key.")
                from subprocess import call
                success = call(['sage-native-execute', 'ssh-keygen', '-q',
                                '-f', private_key, '-P', '', '-t', 'rsa'])
                if success == 0:
                    self._UI.debug("Key generated.")
                else:
                    self._UI.error(["Key generation failed.",
                                    'Please create a key in "{0}" and retry.'.format(public_key)])
                    raise OperationCancelledError("ssh-keygen failed")

            with open(public_key, 'r') as F:
                public_key = F.read().strip()

            self.trac._authenticated_server_proxy.sshkeys.addkey(public_key)
            self._UI.show("Your key has been uploaded.")
        except OperationCancelledError:
            server = self.config.get('server', TRAC_SERVER_URI)

            url = urlparse.urljoin(server, urllib.pathname2url(os.path.join('prefs', 'sshkeys')))
            self._UI.info(['',
                           'Use "{0}" to upload a public key. Or set your key manually at {1}.'
                           .format(self._format_command("upload_ssh_key"), url)])
            raise

    def _upload_ssh_key(self):
        r"""
        Make sure that the public ssh key has been uploaded to the trac server.

        .. NOTE::

            This is a wrapper for :meth:`upload_ssh_key` which is only called
            one the user's first attempt to push to the repository, i.e., on
            the first attempt to acces ``SAGE_REPO_AUTHENTICATED``.

        TESTS:

        Create a doctest setup with a single user::

            sage: from sage.dev.test.sagedev import single_user_setup
            sage: dev, config, UI, server = single_user_setup()
            sage: del dev._sagedev.config['git']['ssh_key_set']

        We need to patch :meth:`upload_ssh_key` to get testable results since
        it depends on whether the user has an ssh key in ``.ssh/id_rsa.pub``::

            sage: from sage.dev.user_interface_error import OperationCancelledError
            sage: def upload_ssh_key():
            ....:     print "Uploading ssh key."
            ....:     raise OperationCancelledError("")
            sage: dev._sagedev.upload_ssh_key = upload_ssh_key

        The ssh key is only uploaded once::

            sage: dev._sagedev._upload_ssh_key()
            Uploading ssh key.
            sage: dev._sagedev._upload_ssh_key()
        """
        if self.config['git'].get('ssh_key_set', False):
            return

        from user_interface_error import OperationCancelledError
        try:
            self.upload_ssh_key()
        except OperationCancelledError:
            pass # do not bother the user again, probably the key has been uploaded manually already
        self.config['git']['ssh_key_set'] = "True"

    def _is_master_uptodate(self, action_if_not=None):
        r"""
        Check whether the master branch is up to date with respect to the
        remote master branch.

        INPUT:

        - ``action_if_not`` -- one of ``'error'``, ``'warning'``, or ``None``
          (default: ``None``), the action to perform if master is not up to
          date. If ``'error'``, then this raises a ``SageDevValueError``,
          otherwise return a boolean and print a warning if ``'warning'``.

        .. NOTE::

            In the transitional period from hg to git, this is a nop. This will
            change as soon as ``master`` is our actual master branch.

        TESTS:

        Create a doctest setup with a single user::

            sage: from sage.dev.test.sagedev import single_user_setup
            sage: dev, config, UI, server = single_user_setup()
            sage: dev._wrap("_is_master_uptodate")

        Initially ``master`` is up to date::

            sage: dev._is_master_uptodate()
            True

        When the remote ``master`` branches changes, this is not the case
        anymore::

            sage: server.git.super_silent.commit(allow_empty=True, message="a commit")
            sage: dev._is_master_uptodate()
            False
            sage: dev._is_master_uptodate(action_if_not="warning")
            Your version of sage, i.e., your "master" branch, is out of date. Your command might fail or produce unexpected results.
            False
            sage: dev._is_master_uptodate(action_if_not="error")
            Your version of sage, i.e., your "master" branch, is out of date.

        We upgrade the local master::

            sage: dev.pull(ticket_or_remote_branch="master")
            Merging the remote branch "master" into the local branch "master".
            Automatic merge successful.
            <BLANKLINE>
            #  (use "sage --dev commit" to commit your merge)
            sage: dev._is_master_uptodate()
            True
            sage: dev._is_master_uptodate(action_if_not="warning")
            True
            sage: dev._is_master_uptodate(action_if_not="error")
            True
        """
        remote_master = self._remote_branch_for_branch(MASTER_BRANCH)
        if remote_master is not None:
            self.git.fetch(self.git._repository_anonymous, remote_master)
            # In the transition from hg to git we are using
            # public/sage-git/master instead of master on the remote end.
            # This check makes sure that we are not printing any confusing
            # messages unless master is actually the latest (development)
            # version of sage.
            if self.git.is_child_of('FETCH_HEAD', MASTER_BRANCH):
                if self.git.commit_for_ref('FETCH_HEAD') != self.git.commit_for_branch(MASTER_BRANCH):
                    msg = ('To upgrade your "{0}" branch to the latest version, use "{1}".',
                           MASTER_BRANCH, self._format_command("pull", ticket_or_branch=remote_master,
                                                               branch=MASTER_BRANCH))
                    if action_if_not is None:
                        pass
                    elif action_if_not == "error":
                        self._UI.debug(*msg)
                        raise SageDevValueError('Your version of sage, i.e., your "{0}" branch, is out'
                                                ' of date.', MASTER_BRANCH)
                    elif action_if_not == "warning":
                        self._UI.warning('Your version of sage, i.e., your "{0}" branch, is out of date.'
                                         ' Your command might fail or produce unexpected results.',
                                         MASTER_BRANCH)
                        self._UI.debug(*msg)
                    else:
                        raise ValueError
                    return False

        return True

    def _is_ticket_name(self, name, exists=False):
        r"""
        Return whether ``name`` is a valid ticket name, i.e., an integer.

        INPUT:

        - ``name`` -- a string or an int

        - ``exists`` -- a boolean (default: ``False``), if ``True``, return
          whether ``name`` is the name of an existing ticket

        EXAMPLES::

            sage: dev = dev._sagedev
            sage: dev._is_ticket_name(1000)
            True
            sage: dev._is_ticket_name("1000")
            True
            sage: dev._is_ticket_name("1 000")
            False
            sage: dev._is_ticket_name("#1000")
            True
            sage: dev._is_ticket_name("master")
            False
            sage: dev._is_ticket_name(1000, exists=True) # optional: internet
            True
            sage: dev._is_ticket_name(2^30, exists=True) # optional: internet
            False
            sage: dev._is_ticket_name('')
            False
        """
        if name is None:
            return False
        if not isinstance(name, int):
            try:
                name = self._ticket_from_ticket_name(name)
            except SageDevValueError:
                return False
        if exists:
            try:
                self.trac._anonymous_server_proxy.ticket.get(name)
            except TracInternalError as e:
                if e.faultCode == 404: # ticket does not exist
                    return False
                raise
            except TracConnectionError as e:
                # if we cannot connect to trac, we assume that the ticket
                # exists; this makes more of the dev scripts usable in offline
                # scenarios
                pass
        return True

    def _check_ticket_name(self, name, exists=False):
        r"""
        Check that ``name`` is a valid ticket name.

        INPUT:

        - ``name`` -- a string or int

        - ``exists`` -- a boolean (default: ``False``), whether to check that
          the ticket exists on trac

        TESTS::

            sage: dev = dev._sagedev
            sage: dev._check_ticket_name(1000)
            sage: dev._check_ticket_name("1000")
            sage: dev._check_ticket_name("1 000")
            Traceback (most recent call last):
            ...
            SageDevValueError: Invalid ticket name "1 000".
            sage: dev._check_ticket_name("#1000")
            sage: dev._check_ticket_name("master")
            Traceback (most recent call last):
            ...
            SageDevValueError: Invalid ticket name "master".
            sage: dev._check_ticket_name(1000, exists=True) # optional: internet
            sage: dev._check_ticket_name(2^30, exists=True) # optional: internet
            Traceback (most recent call last):
            ...
            SageDevValueError: Ticket name "1073741824" is not valid or ticket does not exist on trac.
        """
        if not self._is_ticket_name(name, exists=exists):
            if exists:
                raise SageDevValueError('Ticket name "{0}" is not valid or ticket'
                                        ' does not exist on trac.', name)
            else:
                raise SageDevValueError('Invalid ticket name "{0}".', name)

    def _ticket_from_ticket_name(self, name):
        r"""
        Return the ticket number for the ticket ``name``.

        EXAMPLES::

            sage: from sage.dev.test.sagedev import single_user_setup
            sage: dev, config, UI, server = single_user_setup()
            sage: dev = dev._sagedev

            sage: dev._ticket_from_ticket_name("1000")
            1000
            sage: dev._ticket_from_ticket_name("#1000")
            1000
            sage: dev._ticket_from_ticket_name(1000)
            1000
            sage: dev._ticket_from_ticket_name(int(1000))
            1000
            sage: dev._ticket_from_ticket_name("1 000")
            Traceback (most recent call last):
            ...
            SageDevValueError: "1 000" is not a valid ticket name.
        """
        ticket = name
        if not isinstance(ticket, int):
            if isinstance(ticket, str) and ticket and ticket[0] == "#":
                ticket = ticket[1:]
            try:
                ticket = int(ticket)
            except ValueError:
                raise SageDevValueError('"{0}" is not a valid ticket name.'.format(name))

        if ticket < 0:
            raise SageDevValueError('"{0}" is not a valid ticket name.'.format(name))

        return ticket

    def _is_local_branch_name(self, name, exists=any):
        r"""
        Return whether ``name`` is a valid name for a local branch.

        INPUT:

        - ``name`` -- a string

        - ``exists`` -- a boolean or ``any`` (default: ``any``), if ``True``,
          check whether ``name`` is the name of an existing local branch; if
          ``False``, check whether ``name`` is the name of a branch that does
          not exist yet.

        TESTS::

            sage: from sage.dev.test.sagedev import single_user_setup
            sage: dev, config, UI, server = single_user_setup()
            sage: dev = dev._sagedev

            sage: dev._is_local_branch_name('')
            False
            sage: dev._is_local_branch_name('ticket/1')
            True
            sage: dev._is_local_branch_name('ticket/1', exists=True)
            False
            sage: dev._is_local_branch_name('ticket/1', exists=False)
            True
            sage: dev.git.silent.branch('ticket/1')
            sage: dev._is_local_branch_name('ticket/1', exists=True)
            True
            sage: dev._is_local_branch_name('ticket/1', exists=False)
            False
        """
        if not isinstance(name, str):
            raise ValueError("name must be a string")

        if not GIT_BRANCH_REGEX.match(name):
            return False
        # branches which could be tickets are calling for trouble - cowardly refuse to accept them
        if self._is_ticket_name(name):
            return False
        if name in ["None", "True", "False", "dependencies"]:
            return False

        if exists == True:
            return self.git.commit_for_branch(name) is not None
        elif exists == False:
            return self.git.commit_for_branch(name) is None
        elif exists is any:
            return True
        else:
            raise ValueError("exists")

    def _is_trash_name(self, name, exists=any):
        r"""
        Return whether ``name`` is a valid name for an abandoned branch.

        INPUT:

        - ``name`` -- a string

        - ``exists`` -- a boolean or ``any`` (default: ``any``), if ``True``,
          check whether ``name`` is the name of an existing branch; if
          ``False``, check whether ``name`` is the name of a branch that does
          not exist yet.

        TESTS::

            sage: from sage.dev.test.sagedev import single_user_setup
            sage: dev, config, UI, server = single_user_setup()
            sage: dev = dev._sagedev

            sage: dev._is_trash_name("branch1")
            False
            sage: dev._is_trash_name("trash")
            False
            sage: dev._is_trash_name("trash/")
            False
            sage: dev._is_trash_name("trash/1")
            True
            sage: dev._is_trash_name("trash/1", exists=True)
            False
        """
        if not isinstance(name, str):
            raise ValueError("name must be a string")
        if not name.startswith("trash/"):
            return False
        return self._is_local_branch_name(name, exists)

    def _is_remote_branch_name(self, name, exists=any):
        r"""
        Return whether ``name`` is a valid name for a remote branch.

        INPUT:

        - ``name`` -- a string

        - ``exists`` -- a boolean or ``any`` (default: ``any``), if ``True``,
          check whether ``name`` is the name of an existing remote branch; if
          ``False``, check whether ``name`` is the name of a branch that does
          not exist yet.

        .. NOTE::

            Currently, this does not check whether name is in accordance with
            naming scheme configured on gitolite.

        TESTS::

            sage: from sage.dev.test.sagedev import single_user_setup
            sage: dev, config, UI, server = single_user_setup()
            sage: dev = dev._sagedev

            sage: dev._is_remote_branch_name('')
            False
            sage: dev._is_remote_branch_name('ticket/1')
            True

            sage: dev._is_remote_branch_name('ticket/1', exists=True)
            False
            sage: dev._is_remote_branch_name('ticket/1', exists=False)
            True
        """
        if not isinstance(name, str):
            raise ValueError("name must be a string")

        if not GIT_BRANCH_REGEX.match(name):
            return False
        # branches which could be tickets are calling for trouble - cowardly refuse to accept them
        if self._is_ticket_name(name):
            return False

        if exists is any:
            return True

        from git_error import GitError
        try:
            self.git.super_silent.ls_remote(self.git._repository_anonymous, "refs/heads/"+name, exit_code=True)
            remote_exists = True
        except GitError as e:
            if e.exit_code == 2:
                remote_exists = False
            else:
                raise

        if exists == True or exists == False:
            return remote_exists == exists
        else:
            raise ValueError("exists")

    def _check_local_branch_name(self, name, exists=any):
        r"""
        Check whether ``name`` is a valid name for a local branch, raise a
        ``SageDevValueError`` if it is not.

        INPUT:

        same as for :meth:`_is_local_branch_name`

        TESTS::

            sage: from sage.dev.test.sagedev import single_user_setup
            sage: dev, config, UI, server = single_user_setup()
            sage: dev = dev._sagedev

            sage: dev._check_local_branch_name('')
            Traceback (most recent call last):
            ...
            SageDevValueError: Invalid branch name "".
            sage: dev._check_local_branch_name('ticket/1')
            sage: dev._check_local_branch_name('ticket/1', exists=True)
            Traceback (most recent call last):
            ...
            SageDevValueError: Branch "ticket/1" does not exist locally.
            sage: dev._check_local_branch_name('ticket/1', exists=False)
            sage: dev.git.silent.branch('ticket/1')
            sage: dev._check_local_branch_name('ticket/1', exists=True)
            sage: dev._check_local_branch_name('ticket/1', exists=False)
            Traceback (most recent call last):
            ...
            SageDevValueError: Branch "ticket/1" already exists, use a different name.
        """
        try:
            if not self._is_local_branch_name(name, exists=any):
                raise SageDevValueError("caught below")
        except SageDevValueError:
            raise SageDevValueError('Invalid branch name "{0}".'.format(name))

        if exists == any:
            return
        elif exists == True:
            if not self._is_local_branch_name(name, exists=exists):
                raise SageDevValueError('Branch "{0}" does not exist locally.', name).info(
                    ['', '(use "{0}" to list local branches)'], self._format_command('tickets'))
        elif exists == False:
            if not self._is_local_branch_name(name, exists=exists):
                raise SageDevValueError('Branch "{0}" already exists, use a different name.'.format(name))
        else:
            assert False

    def _check_remote_branch_name(self, name, exists=any):
        r"""
        Check whether ``name`` is a valid name for a remote branch, raise a
        ``SageDevValueError`` if it is not.

        INPUT:

        same as for :meth:`_is_remote_branch_name`

        TESTS::

            sage: from sage.dev.test.sagedev import single_user_setup
            sage: dev, config, UI, server = single_user_setup()
            sage: dev = dev._sagedev

            sage: dev._check_remote_branch_name('')
            Traceback (most recent call last):
            ...
            SageDevValueError: Invalid name "" for a remote branch.
            sage: dev._check_remote_branch_name('ticket/1')

            sage: dev._check_remote_branch_name('ticket/1', exists=True)
            Traceback (most recent call last):
            ...
            SageDevValueError: Branch "ticket/1" does not exist on the remote system.
            sage: dev._check_remote_branch_name('ticket/1', exists=False)
        """
        try:
            if not self._is_remote_branch_name(name, exists=any):
                raise SageDevValueError("caught below")
        except SageDevValueError:
            raise SageDevValueError('Invalid name "{0}" for a remote branch.'.format(name))

        if exists == any:
            return
        elif exists == True:
            if not self._is_remote_branch_name(name, exists=exists):
                raise SageDevValueError('Branch "{0}" does not exist on the remote system.'.format(name))
        elif exists == False:
            if not self._is_remote_branch_name(name, exists=exists):
                raise SageDevValueError('Branch "{0}" already exists, use a different name.'.format(name))
        else:
            assert False

    def _remote_branch_for_ticket(self, ticket):
        r"""
        Return the name of the remote branch for ``ticket``.

        INPUT:

        - ``ticket`` -- an int or a string identifying a ticket

        .. NOTE:

            This does not take into account the ``branch`` field of the ticket
            on trac.

        TESTS::

            sage: from sage.dev.test.sagedev import single_user_setup
            sage: dev, config, UI, server = single_user_setup()
            sage: dev = dev._sagedev

            sage: dev._remote_branch_for_ticket(1)
            'u/doctest/ticket/1'
            sage: dev._remote_branch_for_ticket("#1")
            'u/doctest/ticket/1'
            sage: dev._remote_branch_for_ticket("1")
            'u/doctest/ticket/1'
            sage: dev._remote_branch_for_ticket("master")
            Traceback (most recent call last):
            ...
            SageDevValueError: "master" is not a valid ticket name.

            sage: UI.append("Summary: summary1\ndescription")
            sage: dev.create_ticket()
            Created ticket #1 at https://trac.sagemath.org/1.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=1" to create a new local branch)
            1
            sage: dev.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.

            sage: dev._set_remote_branch_for_branch("ticket/1", "public/1")
            sage: dev._remote_branch_for_ticket(1)
            'public/1'
            sage: dev._set_remote_branch_for_branch("ticket/1", None)
            sage: dev._remote_branch_for_ticket(1)
            'u/doctest/ticket/1'
        """
        ticket = self._ticket_from_ticket_name(ticket)

        default = "u/{0}/ticket/{1}".format(self.trac._username, ticket)

        try:
            branch = self._local_branch_for_ticket(ticket)
        except KeyError: # ticket has no branch yet
            return default

        ret = self._remote_branch_for_branch(branch)
        if ret is None:
            return default
        return ret

    def _ticket_for_local_branch(self, branch):
        r"""
        Return the ticket associated to the local ``branch``.

        INPUT:

        - ``branch`` -- a string, the name of a local branch

        TESTS::

            sage: from sage.dev.test.sagedev import single_user_setup
            sage: dev, config, UI, server = single_user_setup()
            sage: UI.append("Summary: summary\ndescription")
            sage: dev.create_ticket()
            Created ticket #1 at https://trac.sagemath.org/1.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=1" to create a new local branch)
            1
            sage: dev.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.
            sage: dev._sagedev._ticket_for_local_branch("ticket/1")
            1
        """
        self._check_local_branch_name(branch, exists=True)
        if not self._has_ticket_for_local_branch(branch):
            raise SageDevValueError("branch must be associated to a ticket")
        return self.__branch_to_ticket[branch]

    def _has_ticket_for_local_branch(self, branch):
        r"""
        Return whether ``branch`` is associated to a ticket.

        INPUT:

        - ``branch`` -- a string, the name of a local branch

        TESTS::

            sage: from sage.dev.test.sagedev import single_user_setup
            sage: dev, config, UI, server = single_user_setup()
            sage: UI.append("Summary: summary\ndescription")
            sage: dev.create_ticket()
            Created ticket #1 at https://trac.sagemath.org/1.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=1" to create a new local branch)
            1
            sage: dev.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.
            sage: dev._sagedev._has_ticket_for_local_branch("ticket/1")
            True
        """
        self._check_local_branch_name(branch, exists=True)

        return branch in self.__branch_to_ticket

    def _has_local_branch_for_ticket(self, ticket):
        r"""
        Return whether there is a local branch for ``ticket``.

        INPUT:

        - ``ticket`` -- an int or a string identifying a ticket

        TESTS::

            sage: from sage.dev.test.sagedev import single_user_setup
            sage: dev, config, UI, server = single_user_setup()
            sage: dev._sagedev._has_local_branch_for_ticket(1)
            False
        """
        ticket = self._ticket_from_ticket_name(ticket)
        if ticket not in self.__ticket_to_branch:
            return False

        branch = self.__ticket_to_branch[ticket]
        if not self._is_local_branch_name(branch, exists=True):
            self._UI.warning('Ticket #{0} refers to the non-existant local branch "{1}".'
                             ' If you have not manually interacted with git, then this is'
                             ' a bug in sagedev. Removing the association from ticket #{0}'
                             ' to branch "{1}".', ticket, branch)
            del self.__ticket_to_branch[ticket]
            return False
        return True

    def _local_branch_for_ticket(self, ticket, pull_if_not_found=False):
        r"""
        Return the name of the local branch for ``ticket``.

        INPUT:

        - ``ticket`` -- an int or a string identifying a ticket

        - ``pull_if_not_found`` -- a boolean (default: ``False``), whether
          to attempt to pull a branch for ``ticket`` from trac if it does
          not exist locally

        TESTS:

        Create a doctest setup with two users::

            sage: from sage.dev.test.sagedev import two_user_setup
            sage: alice, config_alice, bob, config_bob, server = two_user_setup()

        If a local branch for the ticket exists, its name is returned::

            sage: alice._chdir()
            sage: alice._UI.append("Summary: ticket1\ndescription")
            sage: alice.create_ticket()
            Created ticket #1 at https://trac.sagemath.org/1.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=1" to create a new local branch)
            1
            sage: alice.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.
            sage: alice._sagedev._local_branch_for_ticket(1)
            'ticket/1'

        If no local branch exists, the behaviour depends on ``pull_if_not_found``::

            sage: bob._chdir()
            sage: bob._sagedev._local_branch_for_ticket(1)
            Traceback (most recent call last):
            ...
            KeyError: 'No branch for ticket #1 in your repository.'
            sage: bob._sagedev._local_branch_for_ticket(1, pull_if_not_found=True)
            Traceback (most recent call last):
            ...
            SageDevValueError: Branch field is not set for ticket #1 on trac.
            sage: attributes = alice.trac._get_attributes(1)
            sage: attributes['branch'] = 'public/ticket/1'
            sage: alice.trac._authenticated_server_proxy.ticket.update(1, "", attributes)
            'https://trac.sagemath.org/ticket/1#comment:1'
            sage: bob._sagedev._local_branch_for_ticket(1, pull_if_not_found=True)
            Traceback (most recent call last):
            ...
            SageDevValueError: Branch "public/ticket/1" does not exist on the remote server.

            sage: import os
            sage: os.chdir(server.git._config['src'])
            sage: server.git.silent.branch('public/ticket/1')
            sage: bob._chdir()
            sage: bob._sagedev._local_branch_for_ticket(1, pull_if_not_found=True)
            'ticket/1'
            sage: bob._sagedev._local_branch_for_ticket(1)
            'ticket/1'
        """
        ticket = self._ticket_from_ticket_name(ticket)

        if self._has_local_branch_for_ticket(ticket):
            return self.__ticket_to_branch[ticket]

        if not pull_if_not_found:
            raise KeyError("No branch for ticket #{0} in your repository.".format(ticket))

        branch = self._new_local_branch_for_ticket(ticket)
        self._check_ticket_name(ticket, exists=True)

        remote_branch = self.trac._branch_for_ticket(ticket)
        if remote_branch is None:
            raise SageDevValueError("Branch field is not set for ticket #{0} on trac.".format(ticket))

        try:
            self.git.super_silent.fetch(self.git._repository_anonymous, remote_branch)
        except GitError as e:
            raise SageDevValueError('Branch "%s" does not exist on the remote server.'%remote_branch)

        self.git.super_silent.branch(branch, 'FETCH_HEAD')

        self._set_local_branch_for_ticket(ticket, branch)

        return self._local_branch_for_ticket(ticket, pull_if_not_found=False)

    def _new_local_branch_for_trash(self, branch):
        r"""
        Return a new local branch name to trash ``branch``.

        TESTS::

            sage: from sage.dev.test.sagedev import single_user_setup
            sage: dev, config, UI, server = single_user_setup()
            sage: dev = dev._sagedev

            sage: dev._new_local_branch_for_trash('branch')
            'trash/branch'
            sage: dev.git.silent.branch('trash/branch')
            sage: dev._new_local_branch_for_trash('branch')
            'trash/branch_'
        """
        while True:
            trash_branch = 'trash/{0}'.format(branch)
            if self._is_trash_name(trash_branch, exists=False):
                return trash_branch
            branch = branch + "_"

    def _new_local_branch_for_ticket(self, ticket):
        r"""
        Return a local branch name for ``ticket`` which does not exist yet.

        INPUT:

        - ``ticket`` -- a string or an int identifying a ticket

        TESTS::

            sage: from sage.dev.test.sagedev import single_user_setup
            sage: dev, config, UI, server = single_user_setup()
            sage: dev = dev._sagedev

            sage: dev._new_local_branch_for_ticket(1)
            'ticket/1'
            sage: dev.git.silent.branch('ticket/1')
            sage: dev._new_local_branch_for_ticket(1)
            'ticket/1_'
        """
        ticket = self._ticket_from_ticket_name(ticket)
        branch = 'ticket/{0}'.format(ticket)

        while self._is_local_branch_name(branch, exists=True):
            branch = branch + "_"
        assert self._is_local_branch_name(branch, exists=False)
        return branch

    def _set_dependencies_for_ticket(self, ticket, dependencies):
        r"""
        Locally record ``dependencies`` for ``ticket``.

        INPUT:

        - ``ticket`` -- an int or string identifying a ticket

        - ``dependencies`` -- an iterable of ticket numbers or ``None`` for no
          dependencies

        TESTS::

            sage: from sage.dev.test.sagedev import single_user_setup
            sage: dev, config, UI, server = single_user_setup()
            sage: dev = dev._sagedev

            sage: UI.append("Summary: ticket1\ndescription")
            sage: dev.create_ticket()
            Created ticket #1 at https://trac.sagemath.org/1.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=1" to create a new local branch)
            1
            sage: dev.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.
            sage: dev._set_dependencies_for_ticket(1, [2, 3])
            sage: dev._dependencies_for_ticket(1)
            (2, 3)
            sage: dev._set_dependencies_for_ticket(1, None)
            sage: dev._dependencies_for_ticket(1)
            ()
        """
        ticket = self._ticket_from_ticket_name(ticket)
        if dependencies is None:
            dependencies = []
        dependencies = [self._ticket_from_ticket_name(dep) for dep in dependencies]
        if not(dependencies):
            if ticket in self.__ticket_dependencies:
                del self.__ticket_dependencies[ticket]
            return

        if not self._has_local_branch_for_ticket(ticket):
            raise KeyError("no local branch for ticket #{0} found.".format(ticket))
        self.__ticket_dependencies[ticket] = tuple(sorted(dependencies))

    def _dependencies_for_ticket(self, ticket, download_if_not_found=False):
        r"""
        Return the locally recorded dependencies for ``ticket``.

        INPUT:

        - ``ticket`` -- an int or string identifying a ticket

        - ``download_if_not_found`` -- a boolean (default: ``False``), whether
          to take the information from trac if the ticket does not exist
          locally

        TESTS::

            sage: from sage.dev.test.sagedev import single_user_setup
            sage: dev, config, UI, server = single_user_setup()
            sage: dev = dev._sagedev

            sage: UI.append("Summary: ticket1\ndescription")
            sage: dev.create_ticket()
            Created ticket #1 at https://trac.sagemath.org/1.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=1" to create a new local branch)
            1
            sage: dev.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.

            sage: dev._set_dependencies_for_ticket(1, [2, 3])
            sage: dev._dependencies_for_ticket(1)
            (2, 3)
            sage: dev._set_dependencies_for_ticket(1, None)
            sage: dev._dependencies_for_ticket(1)
            ()

            sage: dev._dependencies_for_ticket(2, download_if_not_found=True)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        ticket = self._ticket_from_ticket_name(ticket)

        if not self._has_local_branch_for_ticket(ticket):
            if download_if_not_found:
                raise NotImplementedError
            else:
                raise KeyError("no local branch for ticket #{0} found.".format(ticket))
        else:
            ret = self.__ticket_dependencies[ticket]

        return tuple(sorted([self._ticket_from_ticket_name(dep) for dep in ret]))

    def _set_remote_branch_for_branch(self, branch, remote_branch):
        r"""
        Set the remote branch of ``branch`` to ``remote_branch``.

        INPUT:

        - ``branch`` -- a string, a name of a local branch

        - ``remote_branch`` -- a string or ``None``, unset the remote branch if
          ``None``

        TESTS::

            sage: from sage.dev.test.sagedev import single_user_setup
            sage: dev, config, UI, server = single_user_setup()
            sage: dev = dev._sagedev

            sage: dev.git.silent.branch('ticket/1')

            sage: dev._remote_branch_for_ticket(1)
            'u/doctest/ticket/1'
            sage: dev._set_remote_branch_for_branch("ticket/1", "public/1")
            sage: dev._remote_branch_for_ticket(1) # ticket/1 has not been set to be the branch for ticket #1
            'u/doctest/ticket/1'
            sage: dev._set_local_branch_for_ticket(1, 'ticket/1')
            sage: dev._remote_branch_for_ticket(1)
            'public/1'
            sage: dev._set_remote_branch_for_branch("ticket/1", None)
            sage: dev._remote_branch_for_ticket(1)
            'u/doctest/ticket/1'
        """
        self._check_local_branch_name(branch, exists=any)

        if remote_branch is None:
            if branch in self.__branch_to_remote_branch:
                del self.__branch_to_remote_branch[branch]
            return

        self._check_local_branch_name(branch, exists=True)
        self._check_remote_branch_name(remote_branch)

        self.__branch_to_remote_branch[branch] = remote_branch

    def _remote_branch_for_branch(self, branch):
        r"""
        Return the remote branch of ``branch`` or ``None`` if no remote branch is set.

        INPUT:

        - ``branch`` -- a string, the name of a local branch

        TESTS::

            sage: from sage.dev.test.sagedev import single_user_setup
            sage: dev, config, UI, server = single_user_setup()
            sage: dev = dev._sagedev

            sage: dev.git.silent.branch('ticket/1')

            sage: dev._remote_branch_for_branch('ticket/1') is None
            True
            sage: dev._set_remote_branch_for_branch("ticket/1", "public/1")
            sage: dev._remote_branch_for_branch('ticket/1')
            'public/1'
            sage: dev._set_remote_branch_for_branch("ticket/1", None)
            sage: dev._remote_branch_for_branch('ticket/1') is None
            True
        """
        self._check_local_branch_name(branch, exists=True)

        if branch in self.__branch_to_remote_branch:
            return self.__branch_to_remote_branch[branch]
        if branch == MASTER_BRANCH:
            return MASTER_BRANCH

        return None

    def _set_local_branch_for_ticket(self, ticket, branch):
        r"""
        Record that ``branch`` is the local branch associated to ``ticket``.

        INPUT:

        - ``ticket`` -- a string or int identifying a ticket

        - ``branch`` -- a string, the name of a local branch, or ``None`` to
          delete the association

        TESTS::

            sage: from sage.dev.test.sagedev import single_user_setup
            sage: dev, config, UI, server = single_user_setup()
            sage: dev = dev._sagedev

            sage: dev._local_branch_for_ticket(1)
            Traceback (most recent call last):
            ...
            KeyError: 'No branch for ticket #1 in your repository.'

            sage: dev._set_local_branch_for_ticket(1, 'ticket/1')
            Traceback (most recent call last):
            ...
            SageDevValueError: Branch "ticket/1" does not exist locally.
            sage: dev.git.silent.branch('ticket/1')
            sage: dev._set_local_branch_for_ticket(1, 'ticket/1')
            sage: dev._local_branch_for_ticket(1)
            'ticket/1'
        """
        ticket = self._ticket_from_ticket_name(ticket)
        if branch is None:
            if ticket in self.__ticket_to_branch:
                del self.__ticket_to_branch[ticket]
            return

        self._check_local_branch_name(branch, exists=True)
        self.__ticket_to_branch[ticket] = branch

    def _format_command(self, command, *args, **kwargs):
        r"""
        Helper method for informational messages.

        OUTPUT:

        A command which the user can run from the command line/sage interactive
        shell to execute ``command`` with ``args`` and ``kwargs``.

        TESTS::

            sage: dev._sagedev._format_command('checkout')
            'sage --dev checkout'

            sage: dev._sagedev._format_command('checkout', ticket=int(1))
            'sage --dev checkout --ticket=1'
        """
        try:
            __IPYTHON__
        except NameError:
            args = [str(arg) for arg in args]
            kwargs = [ "--{0}{1}".format(str(key.split("_or_")[0]).replace("_","-"),"="+str(kwargs[key]) if kwargs[key] is not True else "") for key in kwargs ]
            return "sage --dev {0} {1}".format(command.replace("_","-"), " ".join(args+kwargs)).rstrip()
        else:
            args = [str(arg) for arg in args]
            kwargs = [ "{0}={1}".format(str(key).replace("-","_"),kwargs[key]) for key in kwargs ]
            return "dev.{0}({1})".format(command.replace("-","_"), ", ".join(args+kwargs))

    def _current_ticket(self):
        r"""
        Return the ticket corresponding to the current branch or ``None`` if
        there is no ticket associated to that branch.

        TESTS::

            sage: from sage.dev.test.sagedev import single_user_setup
            sage: dev, config, UI, server = single_user_setup()
            sage: dev = dev._sagedev

            sage: dev._current_ticket() is None
            True

            sage: UI.append("Summary: ticket1\ndescription")
            sage: dev.create_ticket()
            Created ticket #1 at https://trac.sagemath.org/1.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=1" to create a new local branch)
            1
            sage: dev._current_ticket()
            sage: dev.checkout(ticket=1)
            On ticket #1 with associated local branch "ticket/1".
            <BLANKLINE>
            #  Use "sage --dev merge" to include another ticket/branch.
            #  Use "sage --dev commit" to save changes into a new commit.
            sage: dev._current_ticket()
            1
        """
        from git_error import DetachedHeadError
        try:
            branch = self.git.current_branch()
        except DetachedHeadError:
            return None

        if branch in self.__branch_to_ticket:
            return self.__branch_to_ticket[branch]
        return None


class SageDevValueError(ValueError):
    r"""
    A ``ValueError`` to indicate that the user supplied an invaid value.

    EXAMPLES::

        sage: from sage.dev.test.sagedev import single_user_setup
        sage: dev, config, UI, server = single_user_setup()

        sage: dev.checkout(ticket=-1)
        Ticket name "-1" is not valid or ticket does not exist on trac.
    """
    def __init__(self, message, *args):
        r"""
        Initialization.

        TESTS::

            sage: from sage.dev.sagedev import SageDevValueError
            sage: type(SageDevValueError("message"))
            <class 'sage.dev.sagedev.SageDevValueError'>
        """
        ValueError.__init__(self, message.format(*args))
        self._error = (message,) + args
        self._info = None

    def show_error(self, user_interface):
        """
        Display helpful message if available.

        INPUT:

        - ``user_interface`` -- an instance of
          :class:`~sage.dev.user_interface.UserInterface`.

        TESTS::

            sage: from sage.dev.sagedev import SageDevValueError
            sage: e = SageDevValueError("message >{0}<", 123).info('1{0}3', 2)
            sage: e.show_error(dev._sagedev._UI)
            message >123<
        """
        user_interface.error(*self._error)

    def info(self, *args):
        """
        Store helpful message to be displayed if the exception is not
        caught.

        INPUT:

        - ``*args`` -- arguments to be passed to
          :meth:`~sage.dev.user_interface.UserInterface.info`.

        OUTPUT:

        Returns the exception.

        TESTS::

            sage: from sage.dev.sagedev import SageDevValueError
            sage: e = SageDevValueError("message").info('1{0}3', 2)
            sage: e.show_info(dev._sagedev._UI)
            #  123
        """
        self._info = args
        return self

    def show_info(self, user_interface):
        """
        Display helpful message if available.

        INPUT:

        - ``user_interface`` -- an instance of
          :class:`~sage.dev.user_interface.UserInterface`.

        TESTS::

            sage: from sage.dev.sagedev import SageDevValueError
            sage: e = SageDevValueError("message").info('1{0}3', 2)
            sage: e.show_info(dev._sagedev._UI)
            #  123
        """
        if self._info:
            user_interface.info(*self._info)


