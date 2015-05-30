r"""
Configuration wrapper

This module provides a wrapper for the ``devrc`` file in ``DOT_SAGE`` directory
which stores the configuration for :class:`SageDev`.

AUTHORS:

- David Roe, Frej Drejhammar, Julian Rueth, Martin Raum, Nicolas M. Thiery, R.
  Andrew Ohana, Robert Bradshaw, Timo Kluck: initial version

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
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import collections
import os
import ConfigParser as configparser

from sage.env import DOT_SAGE

DEVRC = os.path.join(DOT_SAGE, 'devrc')

class Config(collections.MutableMapping):
    r"""
    Wrapper for the ``devrc`` file storing the configuration for
    :class:`SageDev`.

    INPUT:

    - ``devrc`` -- a string (default: the absolute path of the ``devrc`` file
    in ``DOT_SAGE``)

    EXAMPLES::

        sage: from sage.dev.config import Config
        sage: Config() # random output, the output depends on the user's config
        Config('''
        [trac]
        username = doctest
        ''')
    """
    def __init__(self, devrc = DEVRC):
        r"""
        Initialization.

        TESTS::

            sage: from sage.dev.config import Config
            sage: type(Config())
            <class 'sage.dev.config.Config'>
        """
        self._config = configparser.ConfigParser()
        self._devrc = devrc
        self._read_config()
        self._doctest_config = False

    def __repr__(self):
        r"""
        Return a printable representation of this object.

        EXAMPLES::

            sage: from sage.dev.test.config import DoctestConfig
            sage: c = DoctestConfig()
            sage: repr(c)
            "Config('''\n[trac]\n...\n[UI]\n...\n[git]...\n[sagedev]\n''')"
        """
        config = "".join([ "[%s]\n"%s + 
                           "".join(["%s = %s\n"%(o,self[s][o]) for o in self[s] ]) 
                           for s in self ])
        return "Config('''\n"+ config +"''')"

    def _read_config(self):
        r"""
        Read the configuration from disk.

        EXAMPLES::

            sage: from sage.dev.test.config import DoctestConfig
            sage: c = DoctestConfig()

            sage: from sage.dev.config import Config
            sage: c2 = Config(c._devrc)
            sage: c["trac"]["username"] = "foo"
            sage: c2
            Config('''
            [trac]
            username = doctest
            ticket_cache = ...
            [UI]
            ...
            [git]
            ...
            [sagedev]
            ''')
            sage: c._write_config()
            sage: c2._read_config()
            sage: c2
            Config('''
            [trac]
            username = foo
            ticket_cache = ...
            [UI]
            ...
            [git]
            ...
            [sagedev]
            ''')
        """
        if os.path.exists(self._devrc):
            self._config.read(self._devrc)

    def _write_config(self):
        r"""
        Write the configuration to disk.

        EXAMPLES::

            sage: from sage.dev.test.config import DoctestConfig
            sage: c = DoctestConfig()
            sage: os.unlink(c._devrc)
            sage: os.path.exists(c._devrc)
            False
            sage: c._write_config()
            sage: os.path.exists(c._devrc)
            True
        """
        import sage.doctest
        assert not sage.doctest.DOCTEST_MODE or self._devrc != DEVRC, "attempt to overwrite devrc in doctest"

        with open(self._devrc, 'w') as F:
            self._config.write(F)
        # set the configuration file to read only by this user,
        # because it may contain the trac password
        os.chmod(self._devrc, 0o600)

    def __getitem__(self, section):
        r"""
        Return the configurations in ``section``.

        EXAMPLES::

            sage: from sage.dev.test.config import DoctestConfig
            sage: c = DoctestConfig()
            sage: c['trac']
            IndexableForSection('''
            username = doctest
            ticket_cache = ...
            ''')
            sage: c['tig']
            Traceback (most recent call last):
            ...
            KeyError: 'tig'
        """
        if not section in self:
            raise KeyError(section)

        class IndexableForSection(collections.MutableMapping):
            r"""
            A section of a :class:`Config`.

            TESTS::

                sage: from sage.dev.test.config import DoctestConfig
                sage: c = DoctestConfig()
                sage: c["trac"]
                IndexableForSection('''
                username = doctest
                ticket_cache = ...
                ''')
            """
            def __init__(this, section):
                r"""
                Initialization.

                TESTS::

                    sage: from sage.dev.test.config import DoctestConfig
                    sage: c = DoctestConfig()
                    sage: type(c['trac'])
                    <class 'sage.dev.config.IndexableForSection'>
                """
                this._section = section

            def __repr__(this):
                r"""
                Return a printable representation of this section.

                TESTS::

                    sage: from sage.dev.test.config import DoctestConfig
                    sage: c = DoctestConfig()
                    sage: c["trac"]
                    IndexableForSection('''
                    username = doctest
                    ticket_cache = ...
                    ''')
                """
                return "IndexableForSection('''\n"+"\n".join(["%s = %s"%(o,this[o]) for o in this])+"\n''')"

            def __getitem__(this, option):
                r"""
                Return the value for ``option`` or raise a ``KeyError`` if
                this option is not set.

                TESTS::

                    sage: from sage.dev.test.config import DoctestConfig
                    sage: c = DoctestConfig()
                    sage: c["trac"]["username"]
                    'doctest'
                    sage: c["trac"]["nousername"]
                    Traceback (most recent call last):
                    ...
                    KeyError: 'nousername'
                """
                try:
                    return self._config.get(this._section, option)
                except configparser.NoOptionError:
                    raise KeyError(option)

            def __iter__(this):
                r"""
                Return an iterator over the options in this section.

                TESTS::

                    sage: from sage.dev.test.config import DoctestConfig
                    sage: c = DoctestConfig()
                    sage: list(c["trac"])
                    ['username', 'ticket_cache']
                """
                return iter(self._config.options(this._section))

            def __setitem__(this, option, value):
                r"""
                Set ``option`` to ``value``.

                TESTS::

                    sage: from sage.dev.test.config import DoctestConfig
                    sage: c = DoctestConfig()
                    sage: c["trac"]["username"] = "foo"
                    sage: c["trac"]["username"]
                    'foo'
                """
                self._config.set(this._section, option, value)
                self._write_config()

            def __delitem__(this, option):
                r"""
                Remove ``option`` from this section.

                TESTS::

                    sage: from sage.dev.test.config import DoctestConfig
                    sage: c = DoctestConfig()
                    sage: del(c["trac"]["username"])
                    sage: c["trac"]["username"]
                    Traceback (most recent call last):
                    ...
                    KeyError: 'username'
                """
                self._config.remove_option(this._section, option)
                self._write_config()

            def __len__(this):
                r"""
                Return the number of options in this section.

                TESTS::

                    sage: from sage.dev.test.config import DoctestConfig
                    sage: c = DoctestConfig()
                    sage: len(c["trac"])
                    2
                """
                return len(self._config.options(this._section))

            def __contains__(this, option):
                r"""
                Return whether this section contains ``options``.

                TESTS::

                    sage: from sage.dev.test.config import DoctestConfig
                    sage: c = DoctestConfig()
                    sage: "username" in c["trac"]
                    True
                """
                return option in self._config.options(this._section)

        return IndexableForSection(section)

    def __contains__(self, section):
        r"""
        Returns ``True`` if ``section`` is in the configuration.

        EXAMPLES::

            sage: from sage.dev.test.config import DoctestConfig
            sage: c = DoctestConfig()
            sage: 'trac' in c
            True
            sage: 'nottrac' in c
            False
        """
        return section in self._config.sections()

    def __iter__(self):
        r"""
        Return an iterator over the section names.

        EXAMPLES::

            sage: from sage.dev.test.config import DoctestConfig
            sage: c = DoctestConfig()
            sage: list(c)
            ['trac', 'UI', 'git', 'sagedev']
        """
        return iter(self._config.sections())

    def __setitem__(self, section, dictionary):
        r"""
        Copy the options from ``dictionary`` to ``section``.

        EXAMPLES::

            sage: from sage.dev.test.config import DoctestConfig
            sage: c = DoctestConfig()
            sage: c["foo"] = {"foo":"foo"}
            sage: c["foo"]["foo"]
            'foo'
        """
        was_empty = True

        if self._config.has_section(section):
            was_empty = len(self[section]) == 0
            self._config.remove_section(section)
        self._config.add_section(section)
        for option, value in dictionary.iteritems():
            self._config.set(section, option, value)

        if not was_empty or len(dictionary) != 0:
            self._write_config()

    def __len__(self):
        r"""
        Get the number of sections in the configuration.

        EXAMPLES::

            sage: from sage.dev.test.config import DoctestConfig
            sage: c = DoctestConfig()
            sage: len(c)
            4
        """
        return len(self._config.sections())

    def __delitem__(self, section):
        r"""
        Remove ``section`` from the configuration

        EXAMPLES::

            sage: from sage.dev.test.config import DoctestConfig
            sage: c = DoctestConfig()
            sage: c
            Config('''
            [trac]
            username = doctest
            ticket_cache = ...
            [UI]
            log_level = 1
            [git]
            ssh_key_set = True
            repository_anonymous = remote_repository_undefined
            repository = remote_repository_undefined
            src = ...
            dot_git = ...
            user_email_set = True
            [sagedev]
            ''')
            sage: del c['git']
            sage: c
            Config('''
            [trac]
            username = doctest
            ticket_cache = ...
            [UI]
            log_level = 1
            [sagedev]
            ''')
        """
        self._config.remove_section(section)
        self._write_config()
