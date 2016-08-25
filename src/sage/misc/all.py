from __future__ import absolute_import
from .lazy_attribute import lazy_attribute, lazy_class_attribute
from .lazy_import import lazy_import

from .misc import (alarm, cancel_alarm,
                  ellipsis_range, ellipsis_iter, xsrange, sxrange,
                  BackslashOperator, getitem,
                  cputime, verbose, set_verbose, set_verbose_files,
                  get_verbose_files, unset_verbose_files, get_verbose,
                  union, uniq, powerset, subsets,
                  exists, forall, is_iterator,
                  random_sublist, walltime, generic_cmp,
                  repr_lincomb,
                  pad_zeros, attrcall,
                  SAGE_DB, SAGE_TMP,
                  newton_method_sizes, compose,
                  self_compose, nest)

lazy_import('sage.arith.srange', 'srange', deprecation=20334)

from .banner import version, banner

from .temporary_file import tmp_dir, tmp_filename

from .misc_c import prod, running_total, balanced_sum
lazy_import('sage.misc.misc_c', ['is_32_bit', 'is_64_bit'], deprecation=17460)
mul = prod
add = sum

from .dev_tools import runsnake, import_statements

from .html import html

from .table import table

from .sage_timeit_class import timeit

from .edit_module import edit, set_edit_template

from .flatten import flatten

from .map_threaded import map_threaded

from .session import load_session, save_session, show_identifiers

from .remote_file import get_remote_file

from .profiler import Profiler

from .mrange import xmrange, mrange, xmrange_iter, mrange_iter, cartesian_product_iterator

from .fpickle import pickle_function, unpickle_function

from .dist import install_scripts

from .package import (install_package,
        installed_packages, is_package_installed,
        standard_packages, optional_packages, experimental_packages,
        upgrade, package_versions)

from .pager import pager

lazy_import('sage.misc.sagedoc', ['browse_sage_doc',
        'search_src', 'search_def', 'search_doc',
        'tutorial', 'reference', 'manual', 'developer',
        'constructions', 'python_help', 'help'])

from .classgraph import class_graph

from .reset import reset, restore

from .getusage import top, get_memory_usage

from .mathml import mathml

from .defaults import (set_default_variable_name,
                        series_precision, set_series_precision)

from .sage_eval import sage_eval, sageobj

from .sage_input import sage_input

lazy_import("sage.misc.cython", ["cython_lambda", "cython_create_local_so"]) 
lazy_import("sage.misc.cython_c", "cython_compile", "cython")
lazy_import("sage.misc.cython_c", "cython_compile", "pyrex", deprecation=9552)
lazy_import("sage.misc.cython_c", "cython_compile", "sagex", deprecation=9552)

from .persist import save, load, dumps, loads, db, db_save

from .func_persist import func_persist

from .functional import (additive_order,
                        base_ring,
                        base_field,
                        basis,
                        category,
                        charpoly,
                        characteristic_polynomial,
                        coerce,
                        cyclotomic_polynomial,
                        decomposition,
                        denominator,
                        det,
                        dimension,
                        dim,
                        discriminant,
                        disc,
                        eta,
                        fcp,
                        gen,
                        gens,
                        hecke_operator,
                        image,
                        integral, integrate,
                        integral_closure,
                        interval,
                        xinterval,
                        is_commutative,
                        is_even,
                        is_integrally_closed,
                        is_field,
                        is_odd,
                        kernel,
                        krull_dimension,
                        lift,
                        log as log_b,
                        minimal_polynomial,
                        minpoly,
                        multiplicative_order,
                        ngens,
                        norm,
                        numerator,
                        numerical_approx,
                        n, N,
                        objgens,
                        objgen,
                        one,
                        order,
                        rank,
                        regulator,
                        round,
                        quotient,
                        quo,
                        isqrt,
                        squarefree_part,
                        symbolic_sum as sum,
                        transpose,
                        zero)


from .latex import LatexExpr, latex, view, pretty_print_default

from .trace import trace

from .constant_function import ConstantFunction

from .cachefunc import CachedFunction, cached_function, cached_method, cached_in_parent_method, disk_cached_function

from .abstract_method import abstract_method

from .randstate import seed, set_random_seed, initial_seed, current_randstate

from .prandom import *

from .sage_unittest import TestSuite

from .explain_pickle import explain_pickle, unpickle_newobj, unpickle_global, unpickle_build, unpickle_instantiate, unpickle_persistent, unpickle_extension, unpickle_appends

from .decorators import specialize, sage_wraps, infix_operator

from .unknown import Unknown

lazy_import('sage.misc.inline_fortran', 'fortran')

##########################################################################
def benchmark(n=-1):
    """
    Run a well-chosen range of Sage commands and record the time it
    takes for each to run.

    INPUT:
        n -- int (default: -1) the benchmark number; the default
             of -1 runs all the benchmarks.
    OUTPUT:
        list -- summary of timings for each benchmark.
    """
    import sage.misc.benchmark
    return sage.misc.benchmark.benchmark(n)


class logstr(str):
    def __repr__(self):
        return self

    def _latex_(self):
        #return "\\begin{verbatim}%s\\end{verbatim}"%self
        if not '#' in self:
         delim = '#'
        elif not '@' in self:
         delim = '@'
        elif not '~' in self:
         delim = '~'
        return r"""\verb%s%s%s"""%(delim, self.replace('\n\n','\n').replace('\n','; '), delim)


lazy_import("sage.misc", "messaging", deprecation=18140)

