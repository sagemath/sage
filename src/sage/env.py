"""
Sage Runtime Environment

AUTHORS:

- \R. Andrew Ohana (2012): Initial version.

"""

########################################################################
#       Copyright (C) 2013 R. Andrew Ohana <andrew.ohana@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  http://www.gnu.org/licenses/
########################################################################

import os, socket
import version

opj = os.path.join

# set default values for sage environment variables
# every variable can be overwritten by os.environ
SAGE_ENV = dict()

# Helper to build the SAGE_ENV dictionary
def _add_variable_or_fallback(key, fallback, force=False):
    """
    Set ``SAGE_ENV[key]``.

    If ``key`` is an environment variable, this is the
    value. Otherwise, the ``fallback`` is used.

    INPUT:

    - ``key`` -- string.

    - ``fallback`` -- anything.

    - ``force`` -- boolean (optional, default is ``False``). Whether
      to always use the fallback, regardless of environment variables.

    EXAMPLES::

        sage: import os, sage.env
        sage: sage.env.SAGE_ENV = dict()
        sage: os.environ['SAGE_FOO'] = 'foo'
        sage: sage.env._add_variable_or_fallback('SAGE_FOO', '---$SAGE_URL---')
        sage: sage.env.SAGE_FOO
        'foo'
        sage: sage.env.SAGE_ENV['SAGE_FOO']
        'foo'

    If the environment variable does not exist, the fallback is
    used. Previously-declared variables are replaced if they are
    prefixed with a dollar sign::

        sage: _ = os.environ.pop('SAGE_BAR', None)  # ensure that SAGE_BAR does not exist
        sage: sage.env._add_variable_or_fallback('SAGE_BAR', '---$SAGE_FOO---')
        sage: sage.env.SAGE_BAR
        '---foo---'
        sage: sage.env.SAGE_ENV['SAGE_BAR']
        '---foo---'
    """
    global SAGE_ENV
    try:
        import os
        value = os.environ[key]
    except KeyError:
        value = fallback
    if force:
        value = fallback
    for k,v in SAGE_ENV.iteritems():
        if isinstance(k, basestring):
            value = value.replace('$'+k, v)
    SAGE_ENV[key] = value
    globals()[key] = value

# system info
_add_variable_or_fallback('UNAME',           os.uname()[0])
_add_variable_or_fallback('HOSTNAME',        socket.gethostname())
_add_variable_or_fallback('LOCAL_IDENTIFIER','$HOSTNAME.%s'%os.getpid())

# bunch of sage directories and files
_add_variable_or_fallback('SAGE_ROOT',       None)
_add_variable_or_fallback('SAGE_LOCAL',      opj('$SAGE_ROOT', 'local'))
_add_variable_or_fallback('SAGE_ETC',        opj('$SAGE_LOCAL', 'etc'))
_add_variable_or_fallback('SAGE_SHARE',      opj('$SAGE_LOCAL', 'share'))

# SAGE_LIB is the site-packages directory if the sage library
# has been installed, otherwise it is the same of SAGE_SRC
_add_variable_or_fallback('SAGE_SRC',        opj('$SAGE_ROOT', 'src'))
_add_variable_or_fallback('SAGE_LIB',        os.path.dirname(os.path.dirname(__file__)))

_add_variable_or_fallback('SAGE_EXTCODE',    opj('$SAGE_SHARE', 'sage', 'ext'))
_add_variable_or_fallback('SAGE_LOGS',       opj('$SAGE_ROOT', 'logs', 'pkgs'))
_add_variable_or_fallback('SAGE_SPKG_INST',  opj('$SAGE_LOCAL', 'var', 'lib', 'sage', 'installed'))
_add_variable_or_fallback('SAGE_DOC',        opj('$SAGE_SRC', 'doc'))
_add_variable_or_fallback('DOT_SAGE',        opj(os.environ.get('HOME','$SAGE_ROOT'), '.sage'))
_add_variable_or_fallback('SAGE_DOT_GIT',    opj('$SAGE_ROOT', '.git'))

# misc
_add_variable_or_fallback('SAGE_URL',                'http://sage.math.washington.edu/sage/')
_add_variable_or_fallback('REALM',                   'sage.math.washington.edu')
_add_variable_or_fallback('TRAC_SERVER_URI',         'https://trac.sagemath.org')
_add_variable_or_fallback('SAGE_REPO_AUTHENTICATED', 'ssh://git@trac.sagemath.org:2222/sage.git')
_add_variable_or_fallback('SAGE_REPO_ANONYMOUS',     'git://trac.sagemath.org/sage.git')
_add_variable_or_fallback('SAGE_VERSION',            version.version)
_add_variable_or_fallback('SAGE_DATE',               version.date)

# post process
if ' ' in DOT_SAGE:
    if UNAME[:6] == 'CYGWIN':
        # on windows/cygwin it is typical for the home directory
        # to have a space in it.  Fortunately, users also have
        # write privileges to c:\cygwin\home, so we just put
        # .sage there.
        _add_variable_or_fallback('DOT_SAGE', "/home/.sage", force=True)
    else:
        print("Your home directory has a space in it.  This")
        print("will probably break some functionality of Sage.  E.g.,")
        print("the GAP interface will not work. A workaround")
        print("is to set the environment variable HOME to a")
        print("directory with no spaces that you have write")
        print("permissions to before you start sage.")

# delete temporary variables used for setting up sage.env
del opj, os, socket, version
