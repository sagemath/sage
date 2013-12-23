r"""
Trac Interface

This module provides an interface to access sage's issue tracker 'trac' through
its RPC interface.

AUTHORS:

- David Roe, Julian Rueth, R. Andrew Ohana, Robert Bradshaw, Timo Kluck:
  initial version

"""
#*****************************************************************************
#       Copyright (C) 2013 David Roe <roed.math@gmail.com>
#                          Julian Rueth <julian.rueth@fsfe.org>
#                          R. Andrew Ohana <andrew.ohana@gmail.com>
#                          Robert Bradshaw <robertwb@gmail.com>
#                          Timo Kluck <tkluck@infty.nl>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import os, re, time, datetime
FIELD_REGEX = re.compile("^([A-Za-z ]+):(.*)$")
ALLOWED_FIELDS = {
        "author":           "Authors",
        "branch":           "Branch",
        "cc":               "Cc",
        "component":        "Component",
        "dependencies":     "Dependencies",
        "keywords":         "Keywords",
        "merged":           "Merged in",
        "milestone":        "Milestone",
        "owner":            "Owned by",
        "priority":         "Priority",
        "upstream":         "Report Upstream",
        "reviewer":         "Reviewers",
        "stopgaps":         "Stopgaps",
        "status":           "Status",
        "type":             "Type",
        "work_issues":      "Work issues",
        }
FIELDS_LOOKUP = {"summary":"summary"}
for _k, _v in ALLOWED_FIELDS.iteritems():
    FIELDS_LOOKUP[_v.lower()] = _k
TICKET_FILE_GUIDE = r"""
# Lines starting with `#` are ignored.
# Lines at the beginning of this file starting with `Field: ` correspond to
# fields of the trac ticket, and can be followed by text on the same line.
# They will be assigned to the corresponding field on the trac ticket.
#
# Lines not following this format will be put into the ticket description. Trac
# markup is supported.
#
# An empty file aborts ticket creation/editing.
#
# The following trac fields only allow certain values.
# priority:  blocker, critical, major, minor, trivial
# status:    closed, needs_info, needs_review, needs_work, new,
#            positive_review
# type:      defect, enhancement, task
# milestone: sage-duplicate/invalid/wontfix, sage-feature,
#            sage-pending, sage-wishlist, or sage-VERSION_NUMBER
#            (e.g. sage-6.0)
# component: algebra, algebraic geometry, algebraic topology, basic
#            arithmetic, build, calculus, categories, c_lib, coding
#            theory, coercion, combinatorics, commutative algebra,
#            cryptography, cython, distribution, doctest coverage,
#            doctest framework, documentation, elliptic curves,
#            factorization, finance, finite rings, fractals, geometry,
#            graphics, graph theory, group theory, interact,
#            interfaces, linear algebra, linear programming, matroid
#            theory, memleak, misc, modular forms, notebook, number
#            fields, number theory, numerical, packages: experimental,
#            packages: huge, packages: optional, packages: standard,
#            padics, performance, pickling, PLEASE CHANGE, porting,
#            porting: AIX or HP-UX, porting: BSD, porting: Cygwin,
#            porting: Solaris, quadratic forms, relocation, scripts,
#            spkg-check, statistics, symbolics, user interface,
#            website/wiki
"""
RESTRICTED_FIELDS = ["priority", "status", "type", "milestone", "component"]
ALLOWED_VALUES = {}
_a = 0
for _i, field in enumerate(RESTRICTED_FIELDS):
    _a = TICKET_FILE_GUIDE.find(field + ":", _a) + len(field) + 2
    if _i + 1 == len(RESTRICTED_FIELDS):
        _b = -1
    else:
        _b = TICKET_FILE_GUIDE.find(RESTRICTED_FIELDS[_i+1], _a)
    if field != "milestone":
        ALLOWED_VALUES[field] = re.compile(TICKET_FILE_GUIDE[_a:_b].replace("\n#           ","").replace("\n# ","").strip().replace(', ','|'))
ALLOWED_VALUES['milestone'] = re.compile("sage-(duplicate/invalid/wontfix|feature|pending|wishlist|\d+\.\d+)")
COMMENT_FILE_GUIDE = r"""
# Lines starting with `#` are ignored.
# An empty file aborts the comment.
"""

def _timerep(tm, curtime=None):
    """
    Prints the given time in terms of a rounded-down number of time-units ago.

    INPUT:

    - ``tm`` -- a datetime.datetime instance
    - ``curtime`` -- the current time as a datetime.datetime instance;
      if ``None`` then it will be computed using time.gmtime()

    EXAMPLES::

        sage: from sage.dev.trac_interface import _timerep
        sage: import datetime
        sage: curtime = datetime.datetime(2013, 9, 5, 16)
        sage: T = []
        sage: T.append(datetime.datetime(2013, 9, 5, 15, 54, 22))
        sage: T.append(datetime.datetime(2013, 9, 5, 10, 17, 17))
        sage: T.append(datetime.datetime(2013, 8, 25, 9))
        sage: T.append(datetime.datetime(2012, 11, 19))
        sage: T.append(datetime.datetime(2010, 1, 2))
        sage: for t in T: print _timerep(t, curtime)
        5 minutes ago
        5 hours ago
        11 days ago
        9 months ago
        3 years ago
    """
    year = datetime.timedelta(365)
    month = datetime.timedelta(30)
    day = datetime.timedelta(1)
    hour = datetime.timedelta(0,3600)
    minute = datetime.timedelta(0,60)
    second = datetime.timedelta(0,1)
    timelist = [(year, "year"), (month, "month"), (day, "day"), (hour, "hour"), (minute, "minute"), (second, "second")]
    def timediv(a, b):
        # Our version of datetime doesn't implement //
        x = a.total_seconds() / b.total_seconds()
        if x >= 0: return int(x)
        raise NotImplementedError
    if curtime is None:
        curtime = datetime.datetime(*(time.gmtime()[:6]))
    diff = curtime - tm
    if diff.total_seconds() < 0:
        return "in the future"
    for period, name in timelist:
        n = timediv(diff, period)
        if n > 0:
            return "%s %s%s ago"%(n, name, "s" if n > 1 else "")

class TicketSyntaxError(SyntaxError): # we don't want to catch normal syntax errors
    r"""
    A syntax error when parsing a ticket description modified by the user.

    EXAMPLES::

        sage: from sage.dev.trac_interface import TicketSyntaxError
        sage: raise TicketSyntaxError()
        Traceback (most recent call last):
        ...
        TicketSyntaxError: None
    """

class TracInterface(object):
    r"""
    Wrapper around the XML-RPC interface of trac.

    EXAMPLES::

        sage: from sage.dev.test.config import DoctestConfig
        sage: from sage.dev.test.user_interface import DoctestUserInterface
        sage: from sage.dev.trac_interface import TracInterface
        sage: config = DoctestConfig()
        sage: trac = TracInterface(config['trac'], DoctestUserInterface(config['UI']))
        sage: trac
        <sage.dev.trac_interface.TracInterface object at 0x...>
    """
    def __init__(self, config, UI):
        r"""
        Initialization.

        TESTS::

            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: from sage.dev.trac_interface import TracInterface
            sage: config = DoctestConfig()
            sage: trac = TracInterface(config['trac'], DoctestUserInterface(config['UI']))
            sage: type(trac)
            <class 'sage.dev.trac_interface.TracInterface'>
        """
        self._UI = UI
        self._config = config

        self.__anonymous_server_proxy = None
        self.__authenticated_server_proxy = None

        self.__username = None
        self.__password = None
        self.__auth_timeout = None

        # cache data of tickets locally here
        # this cache is only used to speed
        # up SageDev's local_tickets(), it is not doing any good invalidation
        # and should never be relied on for something other than informational
        # messages
        import os
        from saving_dict import SavingDict
        from sage.env import DOT_SAGE
        ticket_cache_file = self._config.get('ticket_cache', os.path.join(DOT_SAGE,'ticket_cache'))
        self.__ticket_cache = SavingDict(ticket_cache_file)

    @property
    def _username(self):
        r"""
        A lazy property to get the username on trac.

        EXAMPLES::

            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: from sage.dev.trac_interface import TracInterface
            sage: config = DoctestConfig()
            sage: UI = DoctestUserInterface(config['UI'])
            sage: trac = TracInterface(config['trac'], UI)
            sage: trac._username # username is read from config
            'doctest'
            sage: trac.reset_username()
            sage: UI.append('doctest2')
            sage: trac._username # user is prompted for a username
            Trac username: doctest2
            #  Your trac username has been written to a configuration file for future
            sessions. To reset your username, use "dev.trac.reset_username()".
            'doctest2'
            sage: config['trac']['username']
            'doctest2'
        """
        if self.__username is None:
            self.__username = self._config.get('username', None)

        if self.__username is None:
            self.__username = self._config['username'] = self._UI.get_input('Trac username:')
            self._UI.info('Your trac username has been written to a configuration file for'
                          ' future sessions. To reset your username, use "dev.trac.reset_username()".')
        return self.__username

    def reset_username(self):
        r"""
        Reset username and password stored in this object and in the
        configuration.

        EXAMPLES::

            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: from sage.dev.trac_interface import TracInterface
            sage: config = DoctestConfig()
            sage: UI = DoctestUserInterface(config['UI'])
            sage: trac = TracInterface(config['trac'], UI)
            sage: trac.reset_username()
            sage: UI.append("doctest2")
            sage: trac._username
            Trac username: doctest2
            #  Your trac username has been written to a configuration file for future
            sessions. To reset your username, use "dev.trac.reset_username()".
            'doctest2'
        """
        self.__username = None
        if 'username' in self._config:
            del self._config['username']

        self.reset_password()

    @property
    def _password(self):
        r"""
        A lazy property to get the password for trac.

        EXAMPLES::

            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: from sage.dev.trac_interface import TracInterface
            sage: config = DoctestConfig()
            sage: UI = DoctestUserInterface(config['UI'])
            sage: trac = TracInterface(config['trac'], UI)
            sage: UI.append('')
            sage: UI.append('secret')
            sage: trac._password
            Trac password:
            You can save your password in a configuration file. However, this file might be
            readable by privileged users on this system.
            Save password in file? [yes/No]
            'secret'
            sage: trac._password # password is stored for some time, so there is no need to type it immediately afterwards
            'secret'
            sage: config['trac']['password']
            Traceback (most recent call last):
            ...
            KeyError: 'password'
            sage: trac.reset_password()

            sage: UI.append('y')
            sage: UI.append('secret')
            sage: trac._password
            Trac password:
            You can save your password in a configuration file. However, this file might be
            readable by privileged users on this system.
            Save password in file? [yes/No] y
            #  Your trac password has been written to a configuration file. To reset your
            password, use "dev.trac.reset_password()".
            'secret'
            sage: config['trac']['password']
            'secret'
            sage: trac._password
            'secret'
        """
        self._check_password_timeout()

        if self.__password is None:
            self.__password = self._config.get('password', None)

        if self.__password is None:
            self.__password = self._UI.get_password('Trac password:')
            store_password = self._config.get('store_password', None)
            if store_password is None:
                self._UI.show("You can save your password in a configuration file."
                              " However, this file might be readable by privileged"
                              " users on this system.")
                store_password = "yes" if self._UI.confirm(
                    "Save password in file?", default=False) else "no"
                if store_password == "no":
                    # remember the user's decision (if negative) and do not ask every time
                    self._config['store_password'] = store_password
            if store_password == "yes":
                self._config['password'] = self.__password
                self._UI.info('Your trac password has been written to a configuration file. To reset'
                              ' your password, use "dev.trac.reset_password()".')
        self._postpone_password_timeout()
        return self.__password

    def reset_password(self):
        r"""
        Reset password stored in this object and in the configuration.

        EXAMPLES::

            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: from sage.dev.trac_interface import TracInterface
            sage: config = DoctestConfig()
            sage: UI = DoctestUserInterface(config['UI'])
            sage: trac = TracInterface(config['trac'], UI)
            sage: UI.append('y')
            sage: UI.append('secret')
            sage: trac._password
            Trac password:
            You can save your password in a configuration file. However, this file might be
            readable by privileged users on this system.
            Save password in file? [yes/No] y
            #  Your trac password has been written to a configuration file. To reset your
            password, use "dev.trac.reset_password()".
            'secret'
            sage: config['trac']['password']
            'secret'
            sage: trac.reset_password()
            sage: config['trac']['password']
            Traceback (most recent call last):
            ...
            KeyError: 'password'
        """
        self.__password = None
        self.__authenticated_server_proxy = None
        self.__auth_timeout = None
        if 'password' in self._config:
            del self._config['password']
        if 'store_password' in self._config:
            del self._config['store_password']

    def _check_password_timeout(self):
        r"""
        Reset all attributes that depend on the saved password if it has timed
        out (usually after 5 minutes without using it).

        EXAMPLES::

            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: from sage.dev.trac_interface import TracInterface
            sage: config = DoctestConfig()
            sage: UI = DoctestUserInterface(config['UI'])
            sage: trac = TracInterface(config['trac'], UI)
            sage: UI.append('')
            sage: UI.append('secret')
            sage: config['trac']['password_timeout'] = 0
            sage: trac._password
            Trac password:
            You can save your password in a configuration file. However, this file might be
            readable by privileged users on this system.
            Save password in file? [yes/No]
            'secret'
            sage: UI.append('secret')
            sage: trac._password # indirect doctest
            Trac password:
            'secret'
            sage: trac.reset_password()
            sage: UI.append('y')
            sage: UI.append('secret')
            sage: trac._password
            Trac password:
            You can save your password in a configuration file. However, this file might be
            readable by privileged users on this system.
            Save password in file? [yes/No] y
            #  Your trac password has been written to a configuration file. To reset your
            password, use "dev.trac.reset_password()".
            'secret'

        The timeout has no effect if the password can be read from the
        configuration file::

            sage: trac._password
            'secret'
        """
        import time
        if self.__auth_timeout is None or time.time() >= self.__auth_timeout:
            self.__password = None
            self.__authenticated_server_proxy = None
            self.__auth_timeout = None

    def _postpone_password_timeout(self):
        r"""
        Postpone the password timeout.

        EXAMPLES::

            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: from sage.dev.trac_interface import TracInterface
            sage: config = DoctestConfig()
            sage: UI = DoctestUserInterface(config['UI'])
            sage: trac = TracInterface(config['trac'], UI)
            sage: UI.append('')
            sage: UI.append('secret')
            sage: trac._password
            Trac password:
            You can save your password in a configuration file. However, this file might be
            readable by privileged users on this system.
            Save password in file? [yes/No]
            'secret'
            sage: trac._password # indirect doctest
            'secret'
        """
        import time
        new_timeout = time.time() + float(self._config.get('password_timeout', 300))

        if self.__auth_timeout is None or new_timeout > self.__auth_timeout:
            self.__auth_timeout = new_timeout

    @property
    def _anonymous_server_proxy(self):
        """
        Return a non-authenticated XML-RPC interface to trac.

        .. NOTE::

            Unlike the authenticated server proxy, this can be used in
            doctesting. However, all doctests relying on it talking to the
            actual trac server should be marked as ``optional: internet``.

        EXAMPLES::

            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: from sage.dev.trac_interface import TracInterface
            sage: config = DoctestConfig()
            sage: trac = TracInterface(config['trac'], DoctestUserInterface(config['UI']))
            sage: repr(trac._anonymous_server_proxy)
            '<ServerProxy for trac.sagemath.org/xmlrpc>'
        """
        if self.__anonymous_server_proxy is None:
            from sage.env import TRAC_SERVER_URI
            server = self._config.get('server', TRAC_SERVER_URI)
            import urlparse
            url = urlparse.urljoin(server, 'xmlrpc')
            from digest_transport import DigestTransport
            transport = DigestTransport()
            from xmlrpclib import ServerProxy
            self.__anonymous_server_proxy = ServerProxy(url, transport=transport)

        return self.__anonymous_server_proxy

    @property
    def _authenticated_server_proxy(self):
        r"""
        Get an XML-RPC proxy object that is authenticated using the users
        username and password.

        EXAMPLES::

            sage: dev.trac._authenticated_server_proxy  # not tested
            Trac username: username
            Trac password:
            Should I store your password in a configuration file for future sessions? (This configuration file might be readable by privileged users on this system.) [yes/No]
            <ServerProxy for trac.sagemath.org/login/xmlrpc>

        TESTS:

        To make sure that doctests do not tamper with the live trac
        server, it is an error to access this property during a
        doctest (The ``dev`` object during doctests is also modified
        to prevent this)::

            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: from sage.dev.trac_interface import TracInterface
            sage: config = DoctestConfig()
            sage: trac = TracInterface(config['trac'], DoctestUserInterface(config['UI']))
            sage: trac._authenticated_server_proxy
            Traceback (most recent call last):
            ...
            AssertionError: doctest tried to access an authenticated session to trac
        """
        import sage.doctest
        assert not sage.doctest.DOCTEST_MODE, \
            "doctest tried to access an authenticated session to trac"

        self._check_password_timeout()

        if self.__authenticated_server_proxy is None:
            from sage.env import REALM
            realm = self._config.get('realm', REALM)
            from sage.env import TRAC_SERVER_URI
            server = self._config.get('server', TRAC_SERVER_URI)
            import os, urllib, urllib2, urlparse
            url = urlparse.urljoin(server, urllib.pathname2url(os.path.join('login', 'xmlrpc')))
            while True:
                from xmlrpclib import ServerProxy
                from digest_transport import DigestTransport
                from trac_error import TracAuthenticationError
                transport = DigestTransport()
                transport.add_authentication(realm=realm, url=server, username=self._username, password=self._password)
                proxy = ServerProxy(url, transport=transport)
                try:
                    proxy.system.listMethods()
                    break
                except TracAuthenticationError:
                    self._UI.error("Invalid username/password")
                    self.reset_username()
            self.__authenticated_server_proxy = proxy
        self._postpone_password_timeout()
        return self.__authenticated_server_proxy

    def create_ticket(self, summary, description, attributes={}):
        r"""
        Create a ticket on trac and return the new ticket number.

        .. SEEALSO::

            :meth:`create_ticket_interactive`

        EXAMPLES::

            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: from sage.dev.test.trac_interface import DoctestTracInterface
            sage: from sage.dev.test.trac_server import DoctestTracServer
            sage: config = DoctestConfig()
            sage: config['trac']['password'] = 'secret'
            sage: UI = DoctestUserInterface(config['UI'])
            sage: trac = DoctestTracInterface(config['trac'], UI, DoctestTracServer())
            sage: trac.create_ticket('Summary', 'Description', {'type':'defect', 'component':'algebra'})
            1
        """
        return self._authenticated_server_proxy.ticket.create(
            summary, description, attributes, True) # notification e-mail sent.

    def add_comment(self, ticket, comment):
        r"""
        Add ``comment`` to ``ticket`` on trac.

        .. SEEALSO::

            :meth:`add_comment_interactive`

        EXAMPLES::

            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: from sage.dev.test.trac_interface import DoctestTracInterface
            sage: from sage.dev.test.trac_server import DoctestTracServer
            sage: config = DoctestConfig()
            sage: config['trac']['password'] = 'secret'
            sage: UI = DoctestUserInterface(config['UI'])
            sage: trac = DoctestTracInterface(config['trac'], UI, DoctestTracServer())
            sage: ticket = trac.create_ticket('Summary', 'Description', {'type':'defect', 'component':'algebra'})
            sage: trac.add_comment(ticket, "a comment")
        """
        ticket = int(ticket)
        attributes = self._get_attributes(ticket)
        self._authenticated_server_proxy.ticket.update(ticket, comment, attributes, True) # notification e-mail sent

    def _get_attributes(self, ticket, cached=False):
        r"""
        Retrieve the properties of ``ticket``.

        INPUT:

        - ``ticket`` -- an integer, the number of a ticket

        - ``cached`` -- a boolean (default: ``False``), whether to take the
          attributes from a local cache; used, e.g., by
          :meth:`sagedev.SageDev.local_tickets` to speedup display of ticket
          summaries.

        OUTPUT:

        Raises a ``KeyError`` if ``cached`` is ``True`` and the ticket could
        not be found in the cache.

        EXAMPLES::

            sage: from sage.dev.test.sagedev import single_user_setup_with_internet
            sage: dev = single_user_setup_with_internet()[0]
            sage: dev.trac._get_attributes(1000) # optional: internet
            {'status': 'closed',
             'changetime': <DateTime '...' at ...>,
             'description': '...',
             'reporter': 'was',
             'cc': '',
             'type': 'defect',
             'milestone': 'sage-2.10',
             '_ts': '...',
             'component': 'distribution',
             'summary': 'Sage does not have 10000 users yet.',
             'priority': 'major',
             'owner': 'was',
             'time': <DateTime '20071025T16:48:05' at ...>,
             'keywords': '',
             'resolution': 'fixed'}
        """
        ticket = int(ticket)   # must not pickle Sage integers in the ticket cache!
        if not cached:
            self.__ticket_cache[ticket] = self._anonymous_server_proxy.ticket.get(ticket)
        if ticket not in self.__ticket_cache:
            raise KeyError(ticket)
        return self.__ticket_cache[ticket][3]

    def show_ticket(self, ticket):
        r"""
        Show the important fields of the given ticket.

        .. SEEALSO::

            :meth:`_get_attributes`
            :meth:`show_comments`

        EXAMPLES::

            sage: from sage.dev.test.sagedev import single_user_setup_with_internet
            sage: dev = single_user_setup_with_internet()[0]
            sage: dev.trac.show_ticket(101) # optional: internet
            #101: closed enhancement
            == graph theory -- create a graph theory package for SAGE ==
            Opened: ... years ago
            Closed: ... years ago
            Priority: major
            Milestone: sage-2.8.5
            Component: combinatorics
            ----------------------------------------
            See http://sage.math.washington.edu:9001/graph for
            initial research that Robert Miller and Emily Kirkman are doing on this.
        """
        a = self._get_attributes(ticket)
        opentime = datetime.datetime(*(a['time'].timetuple()[:6]))
        changetime = datetime.datetime(*(a['changetime'].timetuple()[:6]))
        fields = ['Opened: %s'%(_timerep(opentime))]
        if a['status'] == 'closed':
            fields.append('Closed: %s'%(_timerep(changetime)))
        else:
            fields.append('Last modified: %s'%(_timerep(changetime)))
        def cap(fieldname):
            if fieldname == 'reporter': return "Reported by"
            elif fieldname == 'merged': return "Merged in"
            elif fieldname == 'work_issues': return "Work issues"
            elif fieldname == 'upstream': return "Report upstream"
            return fieldname.capitalize()
        def add_field(fieldname):
            try:
                fieldval = a[fieldname]
                if fieldval not in ['', 'N/A']:
                    fields.append(ALLOWED_FIELDS[fieldname] + ': ' + fieldval)
            except KeyError:
                pass
        for field in ['priority','milestone','component','cc','merged','author',
                      'reviewer','upstream','work_issues','branch','dependencies','stopgaps']:
            add_field(field)
        self._UI.show("#%s: %s %s\n== %s ==\n%s\n%s\n%s"%(ticket, a['status'], a['type'], a['summary'],
                                                     '\n'.join(fields), '-'*40, a['description']))

    def show_comments(self, ticket, ignore_git_user=True):
        """
        Shows the comments on a given ticket.

        INPUT:

        - ``ticket`` -- the ticket number
        - ``ignore_git_user`` -- whether to remove comments
          automatically added when the branch is updated.

        EXAMPLES::

            sage: dev.trac.show_comments(100) # optional: internet
            ====================
            was (6 years ago)
            fixed
        """
        comments = []
        changelog = self._anonymous_server_proxy.ticket.changeLog(int(ticket))
        for time, author, field, oldvalue, newvalue, permanent in changelog:
            if field == 'comment' and newvalue and not (ignore_git_user and author == 'git'):
                comments.append((_timerep(datetime.datetime(*(time.timetuple()[:6]))), author, newvalue))
        self._UI.show('\n'.join(['====================\n{0} ({1})\n{2}'.format(author, time, comment)
                                 for time, author, comment in reversed(comments)]))

    def query(self, qstr):
        """
        Return a list of ticket ids that match the given query string.

        INPUT:

        - ``qstr`` -- a query string.  All queries will use stored
          settings for maximum number of results per page and paging
          options. Use max=n to define number of results to receive,
          and use page=n to page through larger result sets. Using
          max=0 will turn off paging and return all results.

        EXAMPLES::

            sage: dev.trac.query('status!=closed&(component=padics||component=misc)&max=3') # random, optional: internet
            [329, 15130, 21]
        """
        return self._anonymous_server_proxy.ticket.query(qstr)

    def _branch_for_ticket(self, ticket):
        r"""
        Return the branch field for ``ticket`` or ``None`` if it is not set.

        INPUT:

        - ``ticket`` -- an int

        EXAMPLES::

            sage: from sage.dev.test.sagedev import single_user_setup_with_internet
            sage: dev = single_user_setup_with_internet()[0]
            sage: dev.trac._branch_for_ticket(1000) is None # optional: internet
            True
        """
        attributes = self._get_attributes(ticket)
        if 'branch' in attributes:
            return attributes['branch'] or None
        else:
            return None

    def dependencies(self, ticket, recurse=False, seen=None):
        r"""
        Retrieve dependencies of ``ticket``, sorted by ticket number.

        INPUT:

        - ``ticket`` -- an integer, the number of the ticket

        - ``recurse`` -- a boolean (default: ``False``), whether to get
          indirect dependencies of ``ticket``

        - ``seen`` -- a list (default: ``[]``), used internally to implement
          ``recurse``

       EXAMPLES::

            sage: from sage.dev.test.sagedev import single_user_setup_with_internet
            sage: dev = single_user_setup_with_internet()[0]
            sage: dev.trac.dependencies(1000)            # optional: internet (an old ticket with no dependency field)
            []
            sage: dev.trac.dependencies(13147)           # optional: internet
            [13579, 13681]
            sage: dev.trac.dependencies(13147, recurse=True) # long time, optional: internet
            [13579, 13631, 13681]
        """
        ticket = int(ticket)
        if seen is None:
            seen = []
        if ticket in seen:
            return []
        seen.append(ticket)
        dependencies = self._get_attributes(ticket).get('dependencies','').strip()
        dependencies = dependencies.split(',')
        dependencies = [dep.strip() for dep in dependencies]
        dependencies = [dep for dep in dependencies if dep]
        if not all(dep[0]=="#" for dep in dependencies):
            raise RuntimeError("malformatted dependency on ticket `%s`"%ticket)
        dependencies = [dep[1:] for dep in dependencies]
        try:
            dependencies = [int(dep) for dep in dependencies]
        except ValueError:
            raise RuntimeError("malformatted dependency on ticket `%s`"%ticket)

        if recurse:
            for dep in dependencies:
                self.dependencies(dep, recurse, seen)
        else:
            seen.extend(dependencies)

        ret = sorted(seen)
        ret.remove(ticket)
        return ret

    def attachment_names(self, ticket):
        """
        Retrieve the names of the attachments for ``ticket``.

        EXAMPLES::

            sage: dev.trac.attachment_names(1000) # optional: internet
            ()
            sage: dev.trac.attachment_names(13147) # optional: internet
            ('13147_move.patch',
             '13147_lazy.patch',
             '13147_lazy_spkg.patch',
             '13147_new.patch',
             '13147_over_13579.patch',
             'trac_13147-ref.patch',
             'trac_13147-rebased-to-13681.patch',
             'trac_13681_root.patch')
        """
        ticket = int(ticket)
        return tuple(a[0] for a in self._anonymous_server_proxy.ticket.listAttachments(ticket))

    def add_comment_interactive(self, ticket, comment=''):
        r"""
        Add a comment to ``ticket`` on trac.

        INPUT:

        - ``comment`` -- a string (default: ``''``), the default value for the
          comment to add.

        EXAMPLES::

            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: from sage.dev.test.trac_interface import DoctestTracInterface
            sage: from sage.dev.test.trac_server import DoctestTracServer
            sage: config = DoctestConfig()
            sage: config['trac']['password'] = 'secret'
            sage: UI = DoctestUserInterface(config['UI'])
            sage: trac = DoctestTracInterface(config['trac'], UI, DoctestTracServer())
            sage: ticket = trac.create_ticket('Summary', 'Description', {'type':'defect', 'component':'algebra'})

            sage: UI.append("# empty comment")
            sage: trac.add_comment_interactive(ticket)
            Traceback (most recent call last):
            ...
            OperationCancelledError: comment creation aborted

            sage: UI.append("a comment")
            sage: trac.add_comment_interactive(ticket)
        """
        ticket = int(ticket)

        attributes = self._get_attributes(ticket)

        from sage.dev.misc import tmp_filename
        filename = tmp_filename()
        with open(filename, "w") as f:
            f.write(comment)
            f.write("\n")
            f.write(COMMENT_FILE_GUIDE)
        self._UI.edit(filename)

        comment = list(open(filename).read().splitlines())
        comment = [line for line in comment if not line.startswith("#")]
        if all([line.strip()=="" for line in comment]):
            from user_interface_error import OperationCancelledError
            raise OperationCancelledError("comment creation aborted")
        comment = "\n".join(comment)

        url = self._authenticated_server_proxy.ticket.update(ticket, comment, attributes, True) # notification e-mail sent
        self._UI.debug("Your comment has been recorded: %s"%url)

    def edit_ticket_interactive(self, ticket):
        r"""
        Edit ``ticket`` on trac.

        EXAMPLES::

            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: from sage.dev.test.trac_interface import DoctestTracInterface
            sage: from sage.dev.test.trac_server import DoctestTracServer
            sage: config = DoctestConfig()
            sage: config['trac']['password'] = 'secret'
            sage: UI = DoctestUserInterface(config['UI'])
            sage: trac = DoctestTracInterface(config['trac'], UI, DoctestTracServer())
            sage: ticket = trac.create_ticket('Summary', 'Description', {'type':'defect', 'component':'algebra'})

            sage: UI.append("# empty")
            sage: trac.edit_ticket_interactive(ticket)
            Traceback (most recent call last):
            ...
            OperationCancelledError: ticket edit aborted

            sage: UI.append("Summary: summary\ndescription\n")
            sage: trac.edit_ticket_interactive(ticket)
        """
        ticket = int(ticket)
        attributes = self._get_attributes(ticket)
        summary = attributes.get('summary', 'No Summary')
        description = attributes.get('description', 'No Description')

        ret = self._edit_ticket_interactive(summary, description, attributes)
        if ret is None:
            from user_interface_error import OperationCancelledError
            raise OperationCancelledError("edit aborted")

        attributes['summary'] = ret[0]
        attributes['description'] = ret[1]
        attributes.update(ret[2])

        url = self._authenticated_server_proxy.ticket.update(ticket, "", attributes, True) # notification e-mail sent.
        self._UI.debug("Ticket modified: %s"%url)

    def _edit_ticket_interactive(self, summary, description, attributes):
        r"""
        Helper method for :meth:`edit_ticket_interactive` and
        :meth:`create_ticket_interactive`.

        INPUT:

        - ``summary`` -- a string, summary of ticket

        - ``description`` -- a string, description of ticket

        - ``attributes`` -- dictionary containing field, value pairs

        OUTPUT:

        A tuple ``(summary, description, attributes)``, the updated version of
        input after user has edited the ticket.

        TESTS::

            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: from sage.dev.trac_interface import TracInterface
            sage: config = DoctestConfig()
            sage: UI = DoctestUserInterface(config['UI'])
            sage: trac = TracInterface(config['trac'], UI)
            sage: UI.append("# abort")
            sage: trac._edit_ticket_interactive('summary', 'description', {'branch':'branch1'})
            Traceback (most recent call last):
            ...
            OperationCancelledError: ticket edit aborted

            sage: UI.append("Summary: new summary\nBranch: branch2\nnew description")
            sage: trac._edit_ticket_interactive('summary', 'description', {'branch':'branch1'})
            ('new summary', 'new description', {'branch': 'branch2'})

            sage: UI.append("Summary: new summary\nBranch: branch2\nnew description")
            sage: UI.append("")
            sage: UI.append("Summary: new summary\nInvalid: branch2\nnew description")
            sage: trac._edit_ticket_interactive('summary', 'description', {'branch':'branch1'})
            Syntax error: field "Invalid" not supported on line 2
            Edit ticket file again? [Yes/no]
            ('new summary', 'new description', {'branch': 'branch2'})
        """
        from sage.dev.misc import tmp_filename
        filename = tmp_filename()
        try:
            with open(filename, "w") as F:
                F.write("Summary: %s\n"%summary.encode('utf-8'))
                for k,v in attributes.items():
                    k = ALLOWED_FIELDS.get(k.lower())
                    if k is not None:
                        F.write("%s: %s\n"%(k.encode('utf-8'),v.encode('utf-8')))

                if description is None or not description.strip():
                    description = "\nADD DESCRIPTION\n"
                F.write("\n" + description.encode('utf-8') + "\n")
                F.write(TICKET_FILE_GUIDE)

            while True:
                try:
                    self._UI.edit(filename)
                    ret = self._parse_ticket_file(filename)
                    break
                except (RuntimeError, TicketSyntaxError) as error:
                    pass

                self._UI.error("Syntax error: " + error.message)

                if not self._UI.confirm("Edit ticket file again?", default=True):
                    ret = None
                    break

            if ret is None:
                from user_interface_error import OperationCancelledError
                raise OperationCancelledError("ticket edit aborted")

        finally:
            os.unlink(filename)
        return ret

    def create_ticket_interactive(self):
        r"""
        Drop user into an editor for creating a ticket.

        EXAMPLE::

            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: from sage.dev.test.trac_interface import DoctestTracInterface
            sage: from sage.dev.test.trac_server import DoctestTracServer
            sage: config = DoctestConfig()
            sage: config['trac']['password'] = 'secret'
            sage: UI = DoctestUserInterface(config['UI'])
            sage: trac = DoctestTracInterface(config['trac'], UI, DoctestTracServer())
            sage: UI.append("Summary: summary\nType: defect\nPriority: minor\nComponent: algebra\ndescription")
            sage: trac.create_ticket_interactive()
            1
        """
        attributes = {
                "Type":         "defect",
                "Priority":     "major",
                "Component":    "PLEASE CHANGE",
                "Reporter":     self._username
                }

        ret = self._edit_ticket_interactive("", None, attributes)

        if ret is None:
            from user_interface_error import OperationCancelledError
            raise OperationCancelledError("ticket creation aborted")

        ticket = self.create_ticket(*ret)
        import urlparse
        from sage.env import TRAC_SERVER_URI
        ticket_url = urlparse.urljoin(self._config.get('server', TRAC_SERVER_URI), str(ticket))
        self._UI.debug("Created ticket #%s (%s)."%(ticket, ticket_url))
        return ticket

    @classmethod
    def _parse_ticket_file(cls, filename):
        r"""
        Parse ticket file ``filename``, helper for
        :meth:`create_ticket_interactive` and :meth:`edit_ticket_interactive`.

        OUTPUT:

        ``None`` if the filename contains only comments; otherwise a triple
        ``(summary, description, attributes)``, where ``summary`` is a string
        consisting of the ticket's summary, ``description`` is a string
        containing the ticket's description, and ``attributes`` is a dictionary
        with additional fields of the ticket.

        TESTS::

            sage: from sage.dev.trac_interface import TracInterface
            sage: tmp = tmp_filename()
            sage: with open(tmp, 'w') as f:
            ....:     f.write("no summary\n")
            sage: TracInterface._parse_ticket_file(tmp)
            Traceback (most recent call last):
            ...
            TicketSyntaxError: no valid summary found
            sage: with open(tmp, 'w') as f:
            ....:     f.write("summary:no description\n")
            sage: TracInterface._parse_ticket_file(tmp)
            Traceback (most recent call last):
            ...
            TicketSyntaxError: no description found
            sage: with open(tmp, 'w') as f:
            ....:     f.write("summary:double summary\n")
            ....:     f.write("summary:double summary\n")
            sage: TracInterface._parse_ticket_file(tmp)
            Traceback (most recent call last):
            ...
            TicketSyntaxError: only one value for "summary" allowed on line 2
            sage: with open(tmp, 'w') as f:
            ....:     f.write("bad field:bad field entry\n")
            sage: TracInterface._parse_ticket_file(tmp)
            Traceback (most recent call last):
            ...
            TicketSyntaxError: field "bad field" not supported on line 1
            sage: with open(tmp, 'w') as f:
            ....:     f.write("summary:a summary\n")
            ....:     f.write("branch:a branch\n")
            ....:     f.write("some description\n")
            ....:     f.write("#an ignored line\n")
            ....:     f.write("more description\n")
            ....:     f.write("\n")
            sage: TracInterface._parse_ticket_file(tmp)
            ('a summary', 'some description\nmore description', {'branch': 'a branch'})
            sage: with open(tmp, 'w') as f:
            ....:     f.write("summary:a summary\n")
            ....:     f.write("some description\n")
            ....:     f.write("branch:a branch\n")
            sage: TracInterface._parse_ticket_file(tmp)
            ('a summary', 'some description\nbranch:a branch', {})
            sage: os.unlink(tmp)
        """
        lines = list(open(filename).read().splitlines())

        if all(l.rstrip().startswith('#') for l in lines if l.rstrip()):
            return

        fields = {}
        for i, line in enumerate(lines):
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            m = FIELD_REGEX.match(line)
            if m:
                display_field = m.groups()[0]
                try:
                    field = FIELDS_LOOKUP[display_field.lower()]
                except KeyError:
                    raise TicketSyntaxError('field "{0}" not supported on line {1}'
                                            .format(display_field, i+1))
                if field in fields:
                    raise TicketSyntaxError('only one value for "{0}" allowed on line {1}'
                                            .format(display_field, i+1))
                else:
                    value = m.groups()[1].strip()
                    if field in RESTRICTED_FIELDS and not ALLOWED_VALUES[field].match(value):
                        raise TicketSyntaxError('"{0}" is not a valid value for the field "{1}"'
                                                .format(value, field))
                    fields[field] = value
                    continue
            else:
                break
        else: # no description
            i += 1

        # separate summary from other fields
        try:
            summary = fields.pop('summary')
        except KeyError:
            summary = None

        description = [line.rstrip() for line in lines[i:]
                if not line.startswith('#')]

        # remove leading and trailing empty newlines
        while description and not description[0]:
            description.pop(0)
        while description and not description[-1]:
            description.pop()

        if not summary:
            raise TicketSyntaxError("no valid summary found")
        elif not description:
            raise TicketSyntaxError("no description found")
        else:
            return summary, "\n".join(description), fields

    def set_attributes(self, ticket, comment='', notify=False, **kwds):
        """
        Set attributes on a track ticket.

        INPUT:

        - ``ticket`` -- the ticket id
        - ``comment`` -- a comment when changing these attributes
        - ``kwds`` -- a dictionary of field:value pairs

        .. SEEALSO::

            :meth:`_get_attributes`

        EXAMPLES::

            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: from sage.dev.test.trac_interface import DoctestTracInterface
            sage: from sage.dev.test.trac_server import DoctestTracServer
            sage: config = DoctestConfig()
            sage: config['trac']['password'] = 'secret'
            sage: UI = DoctestUserInterface(config['UI'])
            sage: trac = DoctestTracInterface(config['trac'], UI, DoctestTracServer())
            sage: n = trac.create_ticket('Summary', 'Description', {'type':'defect', 'component':'algebra', 'status':'new'})
            sage: trac._get_attributes(n)['status']
            'new'
            sage: trac.set_attributes(n, status='needs_review')
            sage: trac._get_attributes(n)['status']
            'needs_review'

        Some error checking is done:

            sage: trac.set_attributes(n, status='invalid_status')
            Traceback (most recent call last):
            ...
            TicketSyntaxError: "invalid_status" is not a valid value for the field "status"
        """
        ticket = int(ticket)
        for field, value in kwds.iteritems():
            if field not in ALLOWED_FIELDS:
                raise TicketSyntaxError('field "{0}" not supported'.format(field))
            if field in ALLOWED_VALUES and not ALLOWED_VALUES[field].match(value):
                raise TicketSyntaxError('"{0}" is not a valid value for the field "{1}"'.format(value, field))
        attributes = self._get_attributes(ticket)
        attributes.update(**kwds)
        self._authenticated_server_proxy.ticket.update(ticket, comment, attributes, notify)
        if ticket in self.__ticket_cache:
            del self.__ticket_cache[ticket]
