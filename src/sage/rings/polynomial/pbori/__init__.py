"""The PolyBoRi package implements a framework for computations with Polynomials in Boolean Ring.

The core of PolyBoRi is a C++ library, which provides high-level data types for Boolean polynomials and monomials,
exponent vectors, as well as for the underlying polynomial rings and subsets of the powerset of the Boolean variables.
The description of the latter can be found in the description of the 'dynamic' submodule, as well as in the doxygen-based documentation.

As a unique approach, binary decision diagrams are used as internal storage type for polynomial structures.
On top of this C++-library we provide a Python interface. This allows parsing of complex polynomial systems,
as well as sophisticated and extendable strategies for Groebner base computation.
PolyBoRi features a powerful reference implementation for Groebner basis computation.

AUTHOR:

The PolyBoRi Team, 2007-2011

REFERENCES:
M. Brickenstein, A. Dreyer, G. Greuel, M. Wedler, O. Wienand,
New developments in the theory of Groebner bases and applications
to formal Verification,  Preprint at http://arxiv.org/abs/0801.1177

M. Brickenstein, A. Dreyer, PolyBoRi:
A Groebner Basis Framework for Boolean Polynomials,
Reports of Fraunhofer ITWM, No. 122, Kaiserslautern, Germany, 2007.
http://www.itwm.fraunhofer.de/zentral/download/berichte/bericht122.pdf

M. Brickenstein, A. Dreyer, PolyBoRi:
A framework for Groebner basis computations with Boolean polynomials,
Electronic Proceedings of the MEGA 2007 - Effective Methods in Algebraic Geometry, Strobl, Austria, June 2007.
http://www.ricam.oeaw.ac.at/mega2007/electronic/electronic.html
"""
from sage.misc.lazy_import import lazy_import
from .PyPolyBoRi import *

# Get all-inclusive groebner routine
from .gbcore import groebner_basis
from .nf import normal_form

# Import some high-level modelling functionality
from .blocks import declare_ring
from .blocks import HigherOrderBlock, AlternatingBlock, Block
from .gbrefs import load_file
from .specialsets import all_monomials_of_degree_d, power_set

# Advertised reimports
# ... any from below? ...

# Deprecated reimports

lazy_import('sage.rings.polynomial.pbori.pbori',
            ['BooleConstant',
             'BooleSet',
             'BooleSetIterator',
             'BooleanMonomial',
             'BooleanMonomialIterator',
             'BooleanMonomialMonoid',
             'BooleanMonomialVariableIterator',
             'BooleanMulAction',
             'BooleanPolynomial',
             'BooleanPolynomialEntry',
             'BooleanPolynomialIdeal',
             'BooleanPolynomialIterator',
             'BooleanPolynomialRing',
             'BooleanPolynomialVector',
             'BooleanPolynomialVectorIterator',
             'CCuddNavigator',
             'FGLMStrategy',
             'GroebnerStrategy',
             'MonomialConstruct',
             'MonomialFactory',
             'PolynomialConstruct',
             'PolynomialFactory',
             'ReductionStrategy',
             'TermOrder_from_pb_order',
             'VariableBlock',
             'VariableConstruct',
             'add_up_polynomials',
             'block_dlex',
             'block_dp_asc',
             'contained_vars',
             'dlex',
             'dp',
             'dp_asc',
             'easy_linear_factors',
             'gauss_on_polys',
             'get_var_mapping',
             'if_then_else',
             'interpolate',
             'interpolate_smallest_lex',
             'inv_order_dict',
             'll_red_nf_noredsb',
             'll_red_nf_noredsb_single_recursive_call',
             'll_red_nf_redsb',
             'lp',
             'map_every_x_to_x_plus_one',
             'mod_mon_set',
             'mod_var_set',
             'mult_fact_sim_C',
             'nf3',
             'order_dict',
             'order_mapping',
             'parallel_reduce',
             'random_set',
             'recursively_insert',
             'red_tail',
             'rings',
             'set_random_seed',
             'singular_default',
             'substitute_variables',
             'top_index',
             'unpickle_BooleanPolynomial',
             'unpickle_BooleanPolynomial0',
             'unpickle_BooleanPolynomialRing',
             'zeros'],
            deprecation=30332)
