# Fix build of subprocess32
#
# subprocess32 is a backport from Python 3 to Python 2. But only the C
# sources were backported, not the configure script. This way, macros
# like HAVE_DIRFD which come from the Python 3 configure script are
# always undefined. This causes breakage on certain platforms.
# See upstream bug
# https://github.com/google/python-subprocess32/issues/40
#
# In Sage, we fix this by using the actual pyconfig.h file from our
# Python 3 installation.
#
# This Python script should be run with the Python version where
# subprocess32 will eventually be installed.


import os
from sysconfig import get_path
from subprocess import check_output


# Path to the Python 3 includes
cmd = "from sysconfig import get_path; print(get_path('include'), end='')"
py3incdir = check_output(["python3", "-c", cmd])

# Path to the includes of the Python installation where subprocess32
# will be installed
incdir = get_path("include")


# Create a fake "Python.h" file which includes "pyconfig.h" from
# Python 3 and then includes the real Python.h header
header = '''
/* Include pyconfig.h from Python 3 */
#include "{}/pyconfig.h"

/* Make sure that the Python 2 version of pyconfig.h can also be included */
#undef Py_PYCONFIG_H

/* Include the real Python.h file */
#include "{}/Python.h"
'''.format(py3incdir, incdir)


print("NOTE: Using Python 3 configuration to build subprocess32 for Python 2")

with open("Python.h", "w") as f:
    f.write(header)
