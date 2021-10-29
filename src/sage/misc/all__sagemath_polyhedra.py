# All of sage.misc.all except for development tools, session management,
# and deprecated functionality

from .lazy_attribute import lazy_attribute, lazy_class_attribute
from .lazy_import import lazy_import

from .all__sagemath_objects import *

from .misc import (BackslashOperator,
                  cputime,
                  union, uniq, powerset, subsets,
                  exists, forall, is_iterator,
                  random_sublist, walltime,
                  pad_zeros,
                  SAGE_DB, SAGE_TMP,
                   newton_method_sizes, compose,
                  nest)

from .temporary_file import tmp_dir, tmp_filename

from .html import html, pretty_print_default

from .table import table

from .map_threaded import map_threaded

from .mrange import xmrange, mrange, xmrange_iter, mrange_iter, cartesian_product_iterator

from .fpickle import pickle_function, unpickle_function

from .mathml import mathml

from .defaults import (set_default_variable_name,
                       series_precision, set_series_precision)

from .sage_eval import sage_eval, sageobj

from .sage_input import sage_input

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


from .latex import LatexExpr, latex, view

from .randstate import seed, set_random_seed, initial_seed, current_randstate

from .prandom import *
