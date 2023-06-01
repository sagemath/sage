#from .all__sagemath_objects import *
from .all__sagemath_environment import *
from .all__sagemath_modules import *
from .all__sagemath_repl import *

from .dev_tools import runsnake, import_statements

from .edit_module import edit, set_edit_template

from .session import load_session, save_session, show_identifiers

from .remote_file import get_remote_file

from .profiler import Profiler

from .dist import install_scripts

lazy_import('sage.misc.package', ('installed_packages', 'is_package_installed',
                                  'standard_packages', 'optional_packages',
                                  'experimental_packages', 'package_versions'))

lazy_import('sage.misc.pager', 'pager')

from .classgraph import class_graph

from .reset import reset, restore

lazy_import("sage.misc.cython", "cython_lambda")
lazy_import("sage.misc.cython", "cython_compile", "cython")

from .trace import trace

from .explain_pickle import explain_pickle, unpickle_newobj, unpickle_build, unpickle_instantiate, unpickle_persistent, unpickle_extension, unpickle_appends

lazy_import('sage.misc.inline_fortran', 'fortran')


##########################################################################
def benchmark(n=-1):
    """
    Run a well-chosen range of Sage commands and record the time it
    takes for each to run.

    INPUT:

    - ``n`` -- int (default: -1); the benchmark number. The default
      of -1 runs all the benchmarks.

    OUTPUT:

    - ``list`` -- summary of timings for each benchmark
    """
    import sage.misc.benchmark
    return sage.misc.benchmark.benchmark(n)


class logstr(str):
    def __repr__(self):
        return self

    def _latex_(self):
        # return "\\begin{verbatim}%s\\end{verbatim}"%self
        if '#' not in self:
            delim = '#'
        elif '@' not in self:
            delim = '@'
        elif '~' not in self:
            delim = '~'
        return r"""\verb%s%s%s""" % (delim, self.replace('\n\n', '\n').replace('\n', '; '), delim)
