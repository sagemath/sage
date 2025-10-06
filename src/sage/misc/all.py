from sage.misc.lazy_attribute import lazy_attribute, lazy_class_attribute
from sage.misc.lazy_import import lazy_import

import sage.structure.all   # to break a cyclic import

from sage.misc.verbose import (set_verbose, set_verbose_files,
                               get_verbose_files, unset_verbose_files, get_verbose)
lazy_import('sage.misc.verbose', 'verbose',
            deprecation=17815)
from sage.misc.call import attrcall

from sage.misc.misc_c import prod, running_total, balanced_sum
mul = prod
add = sum

from sage.misc.repr import repr_lincomb

from sage.misc.flatten import flatten

from sage.misc.persist import save, load, dumps, loads, db, db_save

from sage.misc.constant_function import ConstantFunction

from sage.misc.sage_unittest import TestSuite

from sage.misc.decorators import specialize, sage_wraps, infix_operator

from sage.misc.unknown import Unknown, UnknownError

from sage.misc.cachefunc import CachedFunction, cached_function, cached_method, cached_in_parent_method, disk_cached_function

from sage.misc.abstract_method import abstract_method

from sage.misc.timing import walltime, cputime

from sage.misc.temporary_file import tmp_dir, tmp_filename

from sage.misc.sage_eval import sage_eval, sageobj

from sage.misc.sage_input import sage_input

from sage.misc.misc import (
    exists, forall, is_iterator, random_sublist, pad_zeros,
    newton_method_sizes, compose, nest
)

from sage.misc.banner import version

from sage.misc.dev_tools import import_statements

from sage.misc.html import html, pretty_print_default

from sage.misc.table import table

from sage.misc.sage_timeit_class import timeit

from sage.misc.edit_module import edit

from sage.misc.map_threaded import map_threaded

from sage.misc.session import load_session, save_session, show_identifiers

from sage.misc.remote_file import get_remote_file

from sage.misc.mrange import xmrange, mrange, xmrange_iter, mrange_iter, cartesian_product_iterator

from sage.misc.fpickle import pickle_function, unpickle_function

lazy_import('sage.misc.pager', 'pager')

lazy_import('sage.misc.sagedoc', ['browse_sage_doc',
        'search_src', 'search_def', 'search_doc',
        'tutorial', 'reference', 'manual', 'developer',
        'constructions', 'help'])
lazy_import('pydoc', 'help', 'python_help')

from sage.misc.classgraph import class_graph

from sage.misc.reset import reset, restore

from sage.misc.mathml import mathml

from sage.misc.defaults import (set_default_variable_name,
                       series_precision, set_series_precision)

lazy_import("sage.misc.cython", "cython_lambda")
lazy_import("sage.misc.cython", "cython_compile", "cython")

from sage.misc.func_persist import func_persist

from sage.misc.functional import (additive_order,
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
                        is_even,
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
                        order,
                        rank,
                        regulator,
                        round,
                        quotient,
                        quo,
                        isqrt,
                        squarefree_part,
                        sqrt,
                        symbolic_sum as sum,
                        symbolic_prod as product,
                        transpose)


from sage.misc.latex import LatexExpr, latex, view

from sage.misc.randstate import seed, set_random_seed, initial_seed, current_randstate

from sage.misc.prandom import *

from sage.misc.explain_pickle import explain_pickle, unpickle_newobj, unpickle_global, unpickle_build, unpickle_instantiate, unpickle_persistent, unpickle_extension, unpickle_appends

lazy_import('sage.misc.inline_fortran', 'fortran')

lazy_import('sage.misc.banner', 'banner', deprecation=34259)
lazy_import('sage.misc.edit_module', 'set_edit_template', deprecation=34259)
lazy_import('sage.misc.profiler', 'Profiler', deprecation=34259)
lazy_import('sage.misc.trace', 'trace', deprecation=34259)
lazy_import('sage.misc.package', ('installed_packages', 'is_package_installed',
                                  'package_versions'),
            deprecation=34259)
lazy_import('sage.misc.benchmark', 'benchmark', deprecation=34259)
lazy_import('sage.repl.interpreter', 'logstr', deprecation=34259)
