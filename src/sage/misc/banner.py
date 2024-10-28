# sage_setup: distribution = sagemath-repl
r"""
SageMath version and banner info
"""

# ****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
import sys

from sage.env import (SAGE_VERSION, SAGE_VERSION_BANNER, SAGE_BANNER)


def version():
    """
    Return the version of Sage.

    OUTPUT: string

    EXAMPLES::

       sage: version()
       'SageMath version ..., Release Date: ...'
    """
    return SAGE_VERSION_BANNER


def banner_text(full=True):
    """
    Text for the Sage banner.

    INPUT:

    - ``full`` -- boolean (default: ``True``)

    OUTPUT:

    A string containing the banner message.

    If option full is ``False``, a simplified plain ASCII banner is
    displayed; if ``True`` the full banner with box art is displayed.

    EXAMPLES::

        sage: print(sage.misc.banner.banner_text(full=True))
        ┌────────────────────────────────────────────────────────────────────┐
        │ SageMath version ...
        sage: print(sage.misc.banner.banner_text(full=False))
        SageMath version ..., Release Date: ...
    """
    if not full:
        return version()

    bars = "─" * 68
    s = []
    a = s.append
    a('┌' + bars + '┐')
    a("\n│ %-66s │\n" % version())
    python_version = sys.version_info[:3]
    a("│ %-66s │\n" % 'Using Python {}.{}.{}. Type "help()" for help.'.format(*python_version))
    a('└' + bars + '┘')
    pre = version_dict()['prerelease']
    try:
        import sage.all
        have_sage_all = True
    except ImportError:
        have_sage_all = False
    if pre or not have_sage_all:
        red_in = '\033[31m'
        red_out = '\033[0m'
        bars2 = bars.replace('─', '━')
        a('\n')
        a(red_in + '┏' + bars2 + '┓' + '\n')
        if pre:
            a("┃ %-66s ┃\n" % 'Warning: this is a prerelease version, and it may be unstable.')
        if not have_sage_all:
            a("┃ %-66s ┃\n" % 'Warning: sage.all is not available; this is a limited REPL.')
        a('┗' + bars2 + '┛' + red_out)
    return ''.join(s)


def banner():
    """
    Print the Sage banner.

    OUTPUT: none

    If the environment variable ``SAGE_BANNER`` is set to ``no``, no
    banner is displayed. If ``SAGE_BANNER`` is set to ``bare``, a
    simplified plain ASCII banner is displayed. Otherwise, the full
    banner with box art is displayed.

    EXAMPLES::

        sage: import sage.misc.banner; sage.misc.banner.SAGE_BANNER = ''
        sage: sage.misc.banner.banner()
        ┌────────────────────────────────────────────────────────────────────┐
        │ SageMath version ..., Release Date: ...                            │
        │ Using Python .... Type "help()" for help.                          │
        ...
    """
    typ = SAGE_BANNER.lower()

    if typ == "no":
        return

    if typ != "bare":
        try:
            print(banner_text(full=True))
            return
        except UnicodeEncodeError:
            pass

    print(banner_text(full=False))


def version_dict():
    """
    A dictionary describing the version of Sage.

    OUTPUT: dictionary with keys 'major', 'minor', 'tiny', 'prerelease'

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
        int(v[-1])
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
    Return ``True`` if Sage version is at least ``major.minor.tiny``.

    INPUT:

    - ``major`` -- integer
    - ``minor`` -- integer (default: 0)
    - ``tiny`` -- float (default: 0)
    - ``prerelease`` -- boolean (default: ``False``)
    - ``print_message`` -- boolean (default: ``False``)

    OUTPUT: ``True`` if ``major.minor.tiny`` is <= version of Sage, ``False``
    otherwise

    For example, if the Sage version number is 3.1.2, then
    require_version(3, 1, 3) will return False, while
    require_version(3, 1, 2) will return True.
    If the Sage version is 3.1.2.alpha0, then
    require_version(3, 1, 1) will return True, while, by default,
    require_version(3, 1, 2) will return False.  Note, though, that
    require_version(3, 1, 2, prerelease=True) will return True:
    if the optional argument prerelease is True, then a prerelease
    version of Sage counts as if it were the released version.

    If optional argument print_message is ``True`` and this function
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
            txt = "This code requires at least version {} of SageMath to run correctly."
            print(txt.format(major + 0.1 * minor + 0.01 * tiny))
            print("You are running version {}.".format(SAGE_VERSION))
        return False
