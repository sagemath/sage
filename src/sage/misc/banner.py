# -*- coding: utf-8 -*-
r"""
SageMath version and banner info
"""

#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function

from sage.env import SAGE_VERSION, SAGE_DATE, SAGE_SRC, SAGE_BANNER

def version():
    """
    Return the version of Sage.

    OUTPUT:

    str

    EXAMPLES::

       sage: version()
       'SageMath version ..., Release Date: ...'
    """
    return 'SageMath version %s, Release Date: %s' % (SAGE_VERSION, SAGE_DATE)


def banner_text(full=None):
    """
    Text for the Sage banner.

    INPUT:

    - full -- boolean (optional, default = None)

    OUTPUT:

    A string containing the banner message.

    If option full is False, a simplified plain ASCII banner is displayed; if
    True the full banner with box art is displayed.  By default this is
    determined from the SAGE_BANNER environment variable-- if its value is
    "bare" this implies full=False.  Otherwise full=True by default.

    EXAMPLES::

        sage: print(sage.misc.banner.banner_text(full=True))
        ┌────────────────────────────────────────────────────────────────────┐
        │ SageMath version ...
        sage: print(sage.misc.banner.banner_text(full=False))
        SageMath version ..., Release Date: ...
    """

    if full is None:
        full = (SAGE_BANNER.lower() != 'bare')

    if not full:
        return version()

    bars = u"─"*68
    s = []
    a = s.append
    a(u'┌' + bars + u'┐')
    a(u"\n│ %-66s │\n" % version())
    a(u"│ %-66s │\n" % 'Type "notebook()" for the browser-based notebook interface.')
    a(u"│ %-66s │\n" % 'Type "help()" for help.')
    #s += u"│ %-66s │\n" % 'Distributed under the GNU General Public License V2.'
    a(u'└' + bars + u'┘')
    pre = version_dict()['prerelease']
    if pre:
        red_in = '\033[31m'
        red_out = '\033[0m'
        bars2 = bars.replace(u'─', u'━')
        a('\n')
        a(red_in + u'┏' + bars2 + u'┓' + '\n')
        a(u"┃ %-66s ┃\n" % 'Warning: this is a prerelease version, and it may be unstable.')
        a(u'┗' + bars2 + u'┛' + red_out)
    return ''.join(s).encode('utf8')


def banner(full=None):
    """
    Print the Sage banner.

    INPUT:

    - full -- boolean (optional, default = None)

    OUTPUT:

    None

    If option full is False, a simplified plain ASCII banner is displayed; if
    True the full banner with box art is displayed.  By default this is
    determined from the SAGE_BANNER environment variable-- if its value is
    "bare" this implies full=False.  Otherwise full=True by default.

    EXAMPLES::

        sage: banner(full=True)
        ┌────────────────────────────────────────────────────────────────────┐
        │ SageMath version ..., Release Date: ...
        │ Type "notebook()" for the browser-based notebook interface.        │
        │ Type "help()" for help.                                            │
        ...
    """
    print(banner_text(full=full))


def version_dict():
    """
    A dictionary describing the version of Sage.

    INPUT:

    nothing

    OUTPUT:

    dictionary with keys 'major', 'minor', 'tiny', 'prerelease'

    This process the Sage version string and produces a dictionary.
    It expects the Sage version to be in one of these forms::

       N.N
       N.N.N
       N.N.N.N
       N.N.str
       N.N.N.str
       N.N.N.N.str

    where 'N' stands for an integer and 'str' stands for a string.
    The first integer is stored under the 'major' key and the second
    integer under 'minor'.  If there is one more integer, it is stored
    under 'tiny'; if there are two more integers, then they are stored
    together as a float N.N under 'tiny'.  If there is a string, then
    the key 'prerelease' returns True.

    For example, if the Sage version is '3.2.1', then the dictionary
    is {'major': 3, 'minor': 2, 'tiny': 1, 'prerelease': False}.
    If the Sage version is '3.2.1.2', then the dictionary is
    {'major': 3, 'minor': 2, 'tiny': 1.200..., 'prerelease': False}.
    If the Sage version is '3.2.alpha0', then the dictionary is
    {'major': 3, 'minor': 2, 'tiny': 0, 'prerelease': True}.

    EXAMPLES::

        sage: from sage.misc.banner import version_dict
        sage: print("SageMath major version is %s" % version_dict()['major'])
        SageMath major version is ...
        sage: version_dict()['major'] == int(sage.version.version.split('.')[0])
        True
    """
    v = SAGE_VERSION.split('.')
    dict = {}
    dict['major'] = int(v[0])
    dict['minor'] = int(v[1])
    dict['tiny'] = 0
    dict['prerelease'] = False
    try:
        dummy = int(v[-1])
    except ValueError:  # when last entry is not an integer
        dict['prerelease'] = True
    if (len(v) == 3 and not dict['prerelease']) or len(v) > 3:
        dict['tiny'] = int(v[2])
    try:
        teeny = int(v[3])
        dict['tiny'] += 0.1 * teeny
    except (ValueError, IndexError):
        pass
    return dict

def require_version(major, minor=0, tiny=0, prerelease=False,
                    print_message=False):
    """
    True if Sage version is at least major.minor.tiny.

    INPUT:

    - major -- integer
    - minor -- integer (optional, default = 0)
    - tiny -- float (optional, default = 0)
    - prerelease -- boolean (optional, default = False)
    - print_message -- boolean (optional, default = False)

    OUTPUT:

    True if major.minor.tiny is <= version of Sage, False otherwise

    For example, if the Sage version number is 3.1.2, then
    require_version(3, 1, 3) will return False, while
    require_version(3, 1, 2) will return True.
    If the Sage version is 3.1.2.alpha0, then
    require_version(3, 1, 1) will return True, while, by default,
    require_version(3, 1, 2) will return False.  Note, though, that
    require_version(3, 1, 2, prerelease=True) will return True:
    if the optional argument prerelease is True, then a prerelease
    version of Sage counts as if it were the released version.

    If optional argument print_message is True and this function
    is returning False, print a warning message.

    EXAMPLES::

        sage: from sage.misc.banner import require_version
        sage: require_version(2, 1, 3)
        True
        sage: require_version(821, 4)
        False
        sage: require_version(821, 4, print_message=True)
        This code requires at least version 821.4 of SageMath to run correctly.
        You are running version ...
        False
    """
    vers = version_dict()
    prerelease_checked = (prerelease if vers['prerelease'] else True)
    if (vers['major'] > major
        or (vers['major'] == major and vers['minor'] > minor)
        or (vers['major'] == major and vers['minor'] == minor
            and vers['tiny'] > tiny)
        or (vers['major'] == major and vers['minor'] == minor
            and vers['tiny'] == tiny and prerelease_checked)):
        return True
    else:
        if print_message:
            print("This code requires at least version {} of SageMath to run correctly.".
                   format(major + 0.1 * minor + 0.01 * tiny))
            print("You are running version {}.".format(SAGE_VERSION))
        return False
