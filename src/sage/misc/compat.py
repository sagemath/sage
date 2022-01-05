"""Cross-platform compatibility routines and wrappers."""

#*****************************************************************************
#       Copyright (C) 2017 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import os
import subprocess
import sys

from sage.env import SAGE_LOCAL


#################################################################
# Replacements (as needed) for Python stdlib functions to provide
# better platform compatibility
#################################################################
if sys.platform == 'cygwin':
    # find_library that works in cygwin adapted from
    # http://cygwin-ports.svn.sourceforge.net/viewvc/cygwin-ports/ports/trunk/lang/python/2.5.2-ctypes-util-find_library.patch?revision=8245&view=markup
    def _find_library(name):
        libdirs = []
        if SAGE_LOCAL:
            libdirs.append(os.path.join(SAGE_LOCAL, 'lib'))
        libdirs.extend(['/usr/local/lib', '/usr/lib'])
        for libdir in libdirs:
            for libext in ['dll.a', 'a']:
                implib = os.path.join(libdir,
                                      'lib{0}.{1}'.format(name, libext))
                if not os.path.exists(implib):
                    continue

                cmd = ['dlltool', '-I', implib]

                p = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                          stderr=subprocess.PIPE,
                                          universal_newlines=True)

                stdout, stderr = p.communicate()

                if p.returncode == 0:
                    return stdout.strip()
elif sys.platform == 'darwin':
    # On OSX non-standard library paths are not automatically found by the
    # find_library implementation without setting DYLD_LIBRARY_PATH; see
    # https://trac.sagemath.org/ticket/21399#comment:25
    from ctypes.util import find_library as _orig_find_library

    def _find_library(name):
        libdirs = []
        if SAGE_LOCAL:
            libdirs.append(os.path.join(SAGE_LOCAL, 'lib'))
        orig_dyld_library_path = os.environ.get('DYLD_LIBRARY_PATH')
        try:
            if libdirs:
                colon_sep_path = ':'.join(libdirs)
                if orig_dyld_library_path:
                    colon_sep_path += ":" + orig_dyld_library_path
                os.environ['DYLD_LIBRARY_PATH'] = colon_sep_path

            return _orig_find_library(name)
        finally:
            if libdirs:
                # Set os.environ back to what it was
                if orig_dyld_library_path is not None:
                    os.environ['DYLD_LIBRARY_PATH'] = orig_dyld_library_path
                else:
                    os.environ.pop('DYLD_LIBRARY_PATH', None)
else:
    # On other Unix-like platforms, at least where gcc is available,
    # ctypes.util.find_library works, because it takes into account where gcc
    # finds the library as if were being linked with
    from ctypes.util import find_library as _find_library


def find_library(name):
    """
    Returns the shared library filename for a given library.

    The library name is given without any prefixes or suffixes--(e.g.
    just "Singular", not "libSingular", as shared library naming is
    platform-specific.

    This does ''not'' currently return the absolute path of the file on most
    platforms; see https://bugs.python.org/issue21042

    EXAMPLES::

        sage: from sage.misc.compat import find_library
        sage: find_library('Singular')
        '...Singular...'

    """

    result = _find_library(name)
    if result:
        return result

    # if no result found, check DYLD_LIBRARY_PATH
    LDPATH_STR = os.environ.get('DYLD_LIBRARY_PATH')
    lib_dirs = (LDPATH_STR.split(':') if LDPATH_STR else []) + [os.path.join(SAGE_LOCAL, 'lib')]
    for libdir in lib_dirs:
        for libext in ['so', 'a']:
            implib = os.path.join(libdir,
                                  'lib{0}.{1}'.format(name, libext))
            if os.path.exists(implib):
                return implib
