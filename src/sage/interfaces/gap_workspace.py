# sage_setup: distribution = sagemath-gap
r"""
Support for (lib)GAP workspace files
"""

# ****************************************************************************
#       Copyright (C) 2017 Jeroen Demeyer <J.Demeyer@UGent.be>
#                     2019 Vincent Delecroix <vincent.delecroix@u-bordeaux.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import os
import time
import hashlib
import subprocess
from sage.env import DOT_SAGE, HOSTNAME, GAP_ROOT_PATHS


def gap_workspace_file(system='gap', name='workspace', dir=None):
    r"""
    Return the filename for the GAP workspace.

    INPUT:

    - ``system`` -- the name of the system, either ``'gap'`` or
      ``'libgap'``

    - ``name`` -- the kind of workspace, usually ``'workspace'`` but
      the library interface also uses other files

    - ``dir`` -- the directory where the workspaces should be stored
      By default, this is ``DOT_SAGE/gap``

    EXAMPLES::

        sage: from sage.interfaces.gap_workspace import gap_workspace_file
        sage: gap_workspace_file("foo", "bar", "/somewhere")
        '/somewhere/foo-bar-...'

    TESTS::

        sage: from sage.env import DOT_SAGE
        sage: D = gap_workspace_file()
        sage: D.startswith(os.path.join(DOT_SAGE, "gap", "gap-workspace-"))
        True

    Check that the name generated is independent of the session::

        sage: from subprocess import Popen, PIPE
        sage: import sys
        sage: cmd = 'import sage.all, sage.interfaces.gap_workspace; print(sage.interfaces.gap_workspace.gap_workspace_file())'
        sage: name1 = Popen([sys.executable, '-c', cmd], stdout=PIPE).communicate()[0]
        sage: name2 = Popen([sys.executable, '-c', cmd], stdout=PIPE).communicate()[0]
        sage: assert name1 == name2
    """
    if dir is None:
        dir = os.path.join(DOT_SAGE, 'gap')

    data = f'{GAP_ROOT_PATHS}'
    for path in GAP_ROOT_PATHS.split(";"):
        if not path:
            # If GAP_ROOT_PATHS begins or ends with a semicolon,
            # we'll get one empty path.
            continue
        sysinfo = os.path.join(path, "sysinfo.gap")
        if os.path.exists(sysinfo):
            data += subprocess.getoutput(f'. "{sysinfo}" && echo ":$GAP_VERSION:$GAParch"')
    h = hashlib.sha1(data.encode('utf-8')).hexdigest()
    return os.path.join(dir, f'{system}-{name}-{HOSTNAME}-{h}')


def prepare_workspace_dir(dir=None):
    r"""
    Create and clean up the directory for GAP workspaces.

    INPUT:

    - ``dir`` -- the directory where the workspaces should be stored
      By default, this is ``DOT_SAGE/gap``

    OUTPUT: the actual workspace directory

    EXAMPLES::

        sage: from sage.interfaces.gap_workspace import prepare_workspace_dir
        sage: prepare_workspace_dir()
        '.../gap'

    TESTS::

        sage: import tempfile
        sage: with tempfile.TemporaryDirectory() as d:
        ....:     prepare_workspace_dir(os.path.join(d, "new"))
        '.../new'
    """
    if dir is None:
        dir = os.path.join(DOT_SAGE, 'gap')

    # Make sure that the workspace directory exists
    try:
        os.makedirs(dir)
    except OSError:
        if not os.path.isdir(dir):
            raise
    else:
        # Directory was created, add a README file
        with open(os.path.join(dir, 'README.txt'), 'w') as f:
            f.write("It is OK to delete all these cache files.  They will be recreated as needed.\n")

    # Delete all gap workspaces that haven't been used in the last
    # week, to avoid needless cruft.  I had an install on sage.math
    # with 90 of these, since I run a lot of different versions of
    # Sage, and it totalled 1.3GB of wasted space!  See trac #4936.
    # We only do this after creating a new workspace, since this cruft
    # issue is only a problem if workspaces get created every so
    # often.  We don't want to have to do this on every startup.
    now = time.time()
    for F in os.listdir(dir):
        if "workspace" in F:
            W = os.path.join(dir, F)
            try:
                age = now - os.path.getatime(W)
                if age >= 604800:    # 1 week in seconds
                    os.unlink(W)
            except OSError:
                # It's not a problem if W doesn't exist, everything
                # else is an error.
                if os.path.exists(W):
                    raise

    return dir
