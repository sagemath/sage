# sage_setup: distribution = sagemath-repl
"""
Interface to the Sage cleaner

Triva Note: For the name "sage-cleaner", think of the
"The Cleaner" from Pulp Fiction:
http://www.frankjankowski.de/quiz/illus/keitel.jpg
"""
# ****************************************************************************
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
import os
import atexit
import tempfile

_spd = tempfile.TemporaryDirectory()
SAGE_SPAWNED_PROCESS_FILE = os.path.join(_spd.name, "spawned_processes")
atexit.register(lambda: _spd.cleanup())


def start_cleaner():
    """
    Start ``sage-cleaner`` in a new process group.
    """
    if not os.fork():
        os.setpgid(os.getpid(), os.getpid())

        # Redirect stdin, stdout, stderr to /dev/null
        with open(os.devnull, "r+") as f:
            os.dup2(f.fileno(), 0)
            os.dup2(f.fileno(), 1)
            os.dup2(f.fileno(), 2)

        # Close all other file descriptors
        try:
            maxopenfiles = os.sysconf("SC_OPEN_MAX")
            if maxopenfiles <= 0:
                raise ValueError
        except ValueError:
            maxopenfiles = 1024
        os.closerange(3, maxopenfiles)

        os.execlp("sage-cleaner", "sage-cleaner")
