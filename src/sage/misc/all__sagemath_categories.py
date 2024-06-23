# sage_setup: distribution = sagemath-categories

from sage.misc.all__sagemath_objects import *

from sage.misc.html import html, pretty_print_default

from sage.misc.mathml import mathml

from sage.misc.table import table

from sage.misc.map_threaded import map_threaded

from sage.misc.mrange import xmrange, mrange, xmrange_iter, mrange_iter, cartesian_product_iterator

from sage.misc.defaults import (set_default_variable_name,
                                series_precision, set_series_precision)


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

from sage.misc.fpickle import pickle_function, unpickle_function

from sage.misc.persist import unpickle_global
