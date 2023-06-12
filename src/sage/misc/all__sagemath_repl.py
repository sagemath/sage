from .all__sagemath_objects import *

from .sage_eval import sage_eval, sageobj

from .sage_input import sage_input

from .banner import version, banner

lazy_import('sage.misc.sagedoc', ['browse_sage_doc',
        'search_src', 'search_def', 'search_doc',
        'tutorial', 'reference', 'manual', 'developer',
        'constructions', 'help'])

lazy_import('pydoc', 'help', 'python_help')

from .explain_pickle import explain_pickle, unpickle_newobj, unpickle_build, unpickle_instantiate, unpickle_persistent, unpickle_extension, unpickle_appends

from .trace import trace

from .profiler import Profiler

from .dev_tools import runsnake, import_statements

from .edit_module import edit, set_edit_template

lazy_import('sage.misc.pager', 'pager')


lazy_import("sage.misc.cython", "cython_lambda")
lazy_import("sage.misc.cython", "cython_compile", "cython")
lazy_import('sage.misc.inline_fortran', 'fortran')

lazy_import('sage.misc.package', ('installed_packages', 'is_package_installed',
                                  'standard_packages', 'optional_packages',
                                  'experimental_packages', 'package_versions'))


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
