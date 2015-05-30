"""
Auxiliary functions for the dev scripts

These are implemented in the Sage library outside of
:mod:`sage.env`. To reduce the risk of things breaking during
development, a fallback is provided.
"""

#*****************************************************************************
#       Copyright (C) 2013 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


_sage_tmp = None

def get_sage_tmp():
    """
    Get the Sage temporary directory.

    OUTPUT:

    String.

    EXAMPLES::

        sage: from sage.dev.misc import get_sage_tmp
        sage: get_sage_tmp() == SAGE_TMP
        True
    """
    global _sage_tmp
    if _sage_tmp is not None:
        return _sage_tmp
    try:
        from sage.misc.misc import SAGE_TMP
        _sage_tmp = SAGE_TMP
    except ImportError:
        from tempfile import mkdtemp
        _sage_tmp = mkdtemp()
    return _sage_tmp
        

def tmp_filename():
    """
    Return a temporary file.

    OUTPUT:

    String. The absolute filename of the temporary file.

    EXAMPLES::

        sage: from sage.dev.misc import tmp_filename
        sage: tmp_filename().startswith(str(SAGE_TMP))
        True
    """
    try:
        from sage.misc.temporary_file import tmp_filename
        return tmp_filename()
    except ImportError:
        from tempfile import NamedTemporaryFile
        f = NamedTemporaryFile(dir=get_sage_tmp())
        f.close()
        return f.name


def tmp_dir():
    """
    Return a temporary directory.

    OUTPUT:

    String. The absolute filename of the directory.

    EXAMPLES::

        sage: from sage.dev.misc import tmp_dir
        sage: tmp_dir().startswith(str(SAGE_TMP))
        True
    """
    try:
        from sage.misc.all import tmp_dir
        return tmp_dir()
    except ImportError:
        from tempfile import mkdtemp
        return mkdtemp(dir=get_sage_tmp())
    
