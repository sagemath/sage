r"""
Cython wrapper for the Parma Polyhedra Library (PPL)

The Parma Polyhedra Library (PPL) is a library for polyhedral
computations over `\QQ`. This interface tries to reproduce the C++ API
as faithfully as possible in Cython/Sage. For example, the following
C++ excerpt:

.. code-block:: c++

    Variable x(0);
    Variable y(1);
    Constraint_System cs;
    cs.insert(x >= 0);
    cs.insert(x <= 3);
    cs.insert(y >= 0);
    cs.insert(y <= 3);
    C_Polyhedron poly_from_constraints(cs);

translates into::

    sage: from sage.libs.ppl import Variable, Constraint_System, C_Polyhedron
    sage: x = Variable(0)
    sage: y = Variable(1)
    sage: cs = Constraint_System()
    sage: cs.insert(x >= 0)
    sage: cs.insert(x <= 3)
    sage: cs.insert(y >= 0)
    sage: cs.insert(y <= 3)
    sage: poly_from_constraints = C_Polyhedron(cs)

The same polyhedron constructed from generators::

    sage: from sage.libs.ppl import Variable, Generator_System, C_Polyhedron, point
    sage: gs = Generator_System()
    sage: gs.insert(point(0*x + 0*y))
    sage: gs.insert(point(0*x + 3*y))
    sage: gs.insert(point(3*x + 0*y))
    sage: gs.insert(point(3*x + 3*y))
    sage: poly_from_generators = C_Polyhedron(gs)

Rich comparisons test equality/inequality and strict/non-strict
containment::

    sage: poly_from_generators == poly_from_constraints
    True
    sage: poly_from_generators >= poly_from_constraints
    True
    sage: poly_from_generators <  poly_from_constraints
    False
    sage: poly_from_constraints.minimized_generators()
    Generator_System {point(0/1, 0/1), point(0/1, 3/1), point(3/1, 0/1), point(3/1, 3/1)}
    sage: poly_from_constraints.minimized_constraints()
    Constraint_System {-x0+3>=0, -x1+3>=0, x0>=0, x1>=0}

As we see above, the library is generally easy to use. There are a few
pitfalls that are not entirely obvious without consulting the
documentation, in particular:

* There are no vectors used to describe :class:`Generator` (points,
  closure points, rays, lines) or :class:`Constraint` (strict
  inequalities, non-strict inequalities, or equations). Coordinates
  are always specified via linear polynomials in :class:`Variable`

* All coordinates of rays and lines as well as all coefficients of
  constraint relations are (arbitrary precision) integers. Only the
  generators :func:`point` and :func:`closure_point` allow one to
  specify an overall divisor of the otherwise integral
  coordinates. For example::

      sage: from sage.libs.ppl import Variable, point
      sage: x = Variable(0); y = Variable(1)
      sage: p = point( 2*x+3*y, 5 ); p
      point(2/5, 3/5)
      sage: p.coefficient(x)
      2
      sage: p.coefficient(y)
      3
      sage: p.divisor()
      5

* PPL supports (topologically) closed polyhedra
  (:class:`C_Polyhedron`) as well as not neccesarily closed polyhedra
  (:class:`NNC_Polyhedron`). Only the latter allows closure points
  (=points of the closure but not of the actual polyhedron) and strict
  inequalities (``>`` and ``<``)

The naming convention for the C++ classes is that they start with
``PPL_``, for example, the original ``Linear_Expression`` becomes
``PPL_Linear_Expression``. The Python wrapper has the same name as the
original library class, that is, just ``Linear_Expression``. In short:

* If you are using the Python wrapper (if in doubt: thats you), then
  you use the same names as the PPL C++ class library.

* If you are writing your own Cython code, you can access the
  underlying C++ classes by adding the prefix ``PPL_``.

Finally, PPL is fast. For example, here is the permutahedron of 5
basis vectors::

    sage: from sage.libs.ppl import Variable, Generator_System, point, C_Polyhedron
    sage: basis = range(0,5)
    sage: x = [ Variable(i) for i in basis ]
    sage: gs = Generator_System();
    sage: for coeff in Permutations(basis):
    ....:    gs.insert(point( sum( (coeff[i]+1)*x[i] for i in basis ) ))
    sage: C_Polyhedron(gs)
    A 4-dimensional polyhedron in QQ^5 defined as the convex hull of 120 points

The above computation (using PPL) finishes without noticeable delay (timeit
measures it to be 90 microseconds on sage.math). Below we do the same
computation with cddlib, which needs more than 3 seconds on the same
hardware::

    sage: basis = range(0,5)
    sage: gs = [ tuple(coeff) for coeff in Permutations(basis) ]
    sage: Polyhedron(vertices=gs, backend='cdd')  # long time (3s on sage.math, 2011)
    A 4-dimensional polyhedron in QQ^5 defined as the convex hull of 120 vertices

DIFFERENCES VS. C++

Since Python and C++ syntax are not always compatible, there are
necessarily some differences. The main ones are:

* The :class:`Linear_Expression` also accepts an iterable as input for
  the homogeneous cooefficients.

* :class:`Polyhedron` and its subclasses as well as
  :class:`Generator_System` and :class:`Constraint_System` can be set
  immutable via a ``set_immutable()`` method. This is the analog of
  declaring a C++ instance ``const``. All other classes are immutable
  by themselves.

AUTHORS:

- Volker Braun (2010-10-08): initial version.
- Risan (2012-02-19): extension for MIP_Problem class
"""

#*****************************************************************************
#       Copyright (C) 2010 Volker Braun  <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at youroption) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.sage_object cimport SageObject
from sage.libs.gmp.mpz cimport mpz_t, mpz_set
from sage.rings.integer cimport Integer
from sage.rings.rational import Rational

include 'sage/ext/interrupt.pxi'
include "sage/ext/stdsage.pxi"
include "sage/ext/cdefs.pxi"

from libcpp cimport bool as cppbool

####################################################
# Potentially expensive operations:
#  - compute dual description
#  - solve linear program
# These can only be triggered by methods in the Polyhedron class
# they need to be wrapped in sig_on() / sig_off()
####################################################
cdef extern from "gmpxx.h":
    cdef cppclass mpz_class:
        mpz_class()
        mpz_class(int i)
        mpz_class(mpz_t z)
        mpz_class(mpz_class)
        mpz_t* get_mpz_t()


####################################################
# PPL can use floating-point arithmetic to compute integers
cdef extern from "ppl.hh" namespace "Parma_Polyhedra_Library":
    cdef void set_rounding_for_PPL()
    cdef void restore_pre_PPL_rounding()

# but with PPL's rounding the gsl will be very unhappy; must turn off!
restore_pre_PPL_rounding()


####################################################
# Cython does not support ctypedef within cppclass; Hack around this restriction:
cdef extern from "ppl.hh" namespace "Parma_Polyhedra_Library::Generator":
    ctypedef enum PPL_GeneratorType:
        LINE, RAY, POINT, CLOSURE_POINT

cdef extern from "ppl.hh" namespace "Parma_Polyhedra_Library::Constraint":
    ctypedef enum PPL_ConstraintType:
        EQUALITY, NONSTRICT_INEQUALITY, STRICT_INEQUALITY

cdef extern from "ppl.hh" namespace "Parma_Polyhedra_Library::MIP_Problem":
    ctypedef enum PPL_MIP_Problem_Control_Parameter_Name:
        PRICING
    ctypedef enum PPL_MIP_Problem_Control_Parameter_Value:
        PRICING_STEEPEST_EDGE_FLOAT, PRICING_STEEPEST_EDGE_EXACT, PRICING_TEXTBOOK


####################################################
cdef extern from "ppl.hh" namespace "Parma_Polyhedra_Library":

    ctypedef size_t PPL_dimension_type  "Parma_Polyhedra_Library::dimension_type"
    ctypedef mpz_class PPL_Coefficient  "Parma_Polyhedra_Library::Coefficient"
    cdef cppclass PPL_Variable          "Parma_Polyhedra_Library::Variable"
    cdef cppclass PPL_Linear_Expression "Parma_Polyhedra_Library::Linear_Expression"
    cdef cppclass PPL_Generator         "Parma_Polyhedra_Library::Generator"
    cdef cppclass PPL_Generator_System  "Parma_Polyhedra_Library::Generator_System"
    cdef cppclass PPL_Constraint        "Parma_Polyhedra_Library::Constraint"
    cdef cppclass PPL_Constraint_System "Parma_Polyhedra_Library::Constraint_System"
    cdef cppclass PPL_Polyhedron        "Parma_Polyhedra_Library::Polyhedron"
    cdef cppclass PPL_C_Polyhedron      "Parma_Polyhedra_Library::C_Polyhedron"   (PPL_Polyhedron)
    cdef cppclass PPL_NNC_Polyhedron    "Parma_Polyhedra_Library::NNC_Polyhedron" (PPL_Polyhedron)
    cdef cppclass PPL_Poly_Gen_Relation "Parma_Polyhedra_Library::Poly_Gen_Relation"
    cdef cppclass PPL_Poly_Con_Relation "Parma_Polyhedra_Library::Poly_Con_Relation"
    cdef cppclass PPL_MIP_Problem       "Parma_Polyhedra_Library::MIP_Problem"

    cdef cppclass PPL_Variable:
        PPL_Variable(PPL_dimension_type i)
        PPL_dimension_type id()
        bint OK()
        PPL_dimension_type space_dimension()

    cdef cppclass PPL_Linear_Expression:
        PPL_Linear_Expression()
        PPL_Linear_Expression(PPL_Linear_Expression &e)
        PPL_Linear_Expression(PPL_Coefficient n)
        PPL_Linear_Expression(PPL_Variable v)
        PPL_dimension_type space_dimension()
        PPL_Coefficient coefficient(PPL_Variable v)
        PPL_Coefficient inhomogeneous_term()
        bint is_zero()
        bint all_homogeneous_terms_are_zero()
        void ascii_dump()
        bint OK()
        PPL_Linear_Expression operator+(PPL_Linear_Expression& e)
        PPL_Linear_Expression operator-(PPL_Linear_Expression& e)
        PPL_Linear_Expression operator*(PPL_Coefficient n)
        PPL_Constraint operator> (PPL_Linear_Expression& e)
        PPL_Constraint operator>=(PPL_Linear_Expression& e)
        PPL_Constraint operator==(PPL_Linear_Expression& e)
        PPL_Constraint operator<=(PPL_Linear_Expression& e)
        PPL_Constraint operator< (PPL_Linear_Expression& e)

    cdef cppclass PPL_Generator:
        PPL_Generator(PPL_Generator &g)
        # Cython does not support static members
        #PPL_Generator line(PPL_Linear_Expression &e)
        #PPL_Generator ray(PPL_Linear_Expression &e)
        #PPL_Generator point(PPL_Linear_Expression &e, PPL_Coefficient d)
        #PPL_Generator closure_point(PPL_Linear_Expression &e)
        PPL_dimension_type space_dimension()
        PPL_GeneratorType type()
        bint is_line()
        bint is_ray()
        bint is_line_or_ray()
        bint is_point()
        bint is_closure_point()
        PPL_Coefficient coefficient(PPL_Variable v)
        PPL_Coefficient divisor() except +
        bint is_equivalent_to(PPL_Generator &y)
        void ascii_dump()
        bint OK()

    cdef cppclass PPL_Constraint:
        PPL_Constraint(PPL_Constraint &g)
        PPL_dimension_type space_dimension()
        PPL_ConstraintType type()
        bint is_equality()
        bint is_inequality()
        bint is_nonstrict_inequality()
        bint is_strict_inequality()
        PPL_Coefficient coefficient(PPL_Variable v)
        PPL_Coefficient inhomogeneous_term()
        bint is_tautological()
        bint is_inconsistent()
        bint is_equivalent_to(PPL_Constraint &y)
        void ascii_dump()
        bint OK()

    cdef cppclass PPL_Generator_System:
        # This seems to not work in cython
        #cppclass PPL_const_iterator "const_iterator":
        #    PPL_Generator operator*()
        #    PPL_const_iterator operator++()
        #    bint operator==(PPL_const_iterator&)
        #    bint operator!=(PPL_const_iterator&)
        #PPL_const_iterator begin()
        #PPL_const_iterator end()
        PPL_Generator_System()
        PPL_Generator_System(PPL_Generator &g)
        PPL_Generator_System(PPL_Generator_System &gs)
        PPL_dimension_type space_dimension()
        void clear()
        void insert(PPL_Generator &g)
        bint empty()
        void ascii_dump()
        bint OK()

    cdef cppclass PPL_Constraint_System:
        # This seems to not work in cython
        #cppclass PPL_const_iterator "const_iterator":
        #    PPL_Constraint operator*()
        #    PPL_const_iterator operator++()
        #    bint operator==(PPL_const_iterator&)
        #    bint operator!=(PPL_const_iterator&)
        #PPL_const_iterator begin()
        #PPL_const_iterator end()
        PPL_Constraint_System()
        PPL_Constraint_System(PPL_Constraint &g)
        PPL_Constraint_System(PPL_Constraint_System &gs)
        PPL_dimension_type space_dimension()
        bint has_equalities()
        bint has_strict_inequalities()
        void clear()
        void insert(PPL_Constraint &g)
        bint empty()
        void ascii_dump()
        bint OK()

    cdef enum PPL_Degenerate_Element:
        UNIVERSE, EMPTY

    cdef enum PPL_Optimization_Mode:
        MINIMIZATION, MAXIMIZATION

    cdef enum PPL_MIP_Problem_Status:
        UNFEASIBLE_MIP_PROBLEM, UNBOUNDED_MIP_PROBLEM, OPTIMIZED_MIP_PROBLEM

    cdef cppclass PPL_Polyhedron:
        PPL_dimension_type space_dimension()
        PPL_dimension_type affine_dimension()
        PPL_Constraint_System& constraints()
        PPL_Constraint_System& minimized_constraints()
        PPL_Generator_System& generators()
        PPL_Generator_System& minimized_generators()
        PPL_Poly_Con_Relation relation_with(PPL_Constraint &c) except +ValueError
        PPL_Poly_Gen_Relation relation_with(PPL_Generator &g) except +ValueError
        bint is_empty()
        bint is_universe()
        bint is_topologically_closed()
        bint is_disjoint_from(PPL_Polyhedron &y) except +ValueError
        bint is_discrete()
        bint is_bounded()
        bint contains_integer_point()
        bint constrains(PPL_Variable var) except +ValueError
        bint bounds_from_above(PPL_Linear_Expression &expr) except +ValueError
        bint bounds_from_below(PPL_Linear_Expression &expr) except +ValueError
        bint maximize(PPL_Linear_Expression &expr, PPL_Coefficient &sup_n, PPL_Coefficient &sup_d,
                      cppbool &maximum)
        bint maximize(PPL_Linear_Expression &expr, PPL_Coefficient &sup_n, PPL_Coefficient &sup_d,
                      cppbool &maximum, PPL_Generator &g)
        bint minimize(PPL_Linear_Expression &expr, PPL_Coefficient &inf_n, PPL_Coefficient &inf_d,
                      cppbool &minimum)
        bint minimize(PPL_Linear_Expression &expr, PPL_Coefficient &inf_n, PPL_Coefficient &inf_d,
                      cppbool &minimum, PPL_Generator &g)
        bint frequency(PPL_Linear_Expression &expr, PPL_Coefficient &freq_n, PPL_Coefficient &freq_d,
                       PPL_Coefficient &val_n, PPL_Coefficient &val_d)
        bint contains(PPL_Polyhedron &y) except +ValueError
        bint strictly_contains(PPL_Polyhedron &y) except +ValueError
        void add_constraint(PPL_Constraint &c) except +ValueError
        void add_generator(PPL_Generator &g) except +ValueError
        void add_constraints(PPL_Constraint_System &cs) except +ValueError
        void add_generators(PPL_Generator_System &gs) except +ValueError
        void refine_with_constraint(PPL_Constraint &c) except +ValueError
        void refine_with_constraints(PPL_Constraint_System &cs) except +ValueError
        void unconstrain(PPL_Variable var) except +ValueError
        void intersection_assign(PPL_Polyhedron &y) except +ValueError
        void poly_hull_assign(PPL_Polyhedron &y) except +ValueError
        void upper_bound_assign(PPL_Polyhedron &y) except +ValueError
        void poly_difference_assign(PPL_Polyhedron &y) except +ValueError
        void difference_assign(PPL_Polyhedron &y) except +ValueError
        void drop_some_non_integer_points()
        void topological_closure_assign()
        void add_space_dimensions_and_embed(PPL_dimension_type m) except +ValueError
        void add_space_dimensions_and_project(PPL_dimension_type m) except +ValueError
        void concatenate_assign(PPL_Polyhedron &y) except +ValueError
        void remove_higher_space_dimensions(PPL_dimension_type new_dimension) except +ValueError
        void ascii_dump()
        int hash_code()
        PPL_dimension_type max_space_dimension()
        bint OK(cppbool check_not_empty=false)
        bint operator!=(PPL_Polyhedron &y)
        bint operator==(PPL_Polyhedron &y)

    cdef cppclass PPL_C_Polyhedron(PPL_Polyhedron):
        PPL_C_Polyhedron(PPL_dimension_type num_dimensions, PPL_Degenerate_Element)
        PPL_C_Polyhedron(PPL_Constraint_System &cs) except +ValueError
        PPL_C_Polyhedron(PPL_Generator_System &gs) except +ValueError
        PPL_C_Polyhedron(PPL_C_Polyhedron &y)

    cdef cppclass PPL_NNC_Polyhedron(PPL_Polyhedron):
        PPL_NNC_Polyhedron(PPL_dimension_type num_dimensions, PPL_Degenerate_Element kind)
        PPL_NNC_Polyhedron(PPL_Constraint_System &cs) except +ValueError
        PPL_NNC_Polyhedron(PPL_Generator_System &gs) except +ValueError
        PPL_NNC_Polyhedron(PPL_NNC_Polyhedron &y)
        PPL_NNC_Polyhedron(PPL_C_Polyhedron &y)

    cdef cppclass PPL_Poly_Gen_Relation:
        PPL_Poly_Gen_Relation(PPL_Poly_Gen_Relation &cpy_from)
        bint implies(PPL_Poly_Gen_Relation &y)
        void ascii_dump()
        bint OK()

    cdef cppclass PPL_Poly_Con_Relation:
        PPL_Poly_Con_Relation(PPL_Poly_Con_Relation &cpy_from)
        bint implies(PPL_Poly_Con_Relation &y)
        void ascii_dump()
        bint OK()

    cdef cppclass PPL_MIP_Problem:
        PPL_MIP_Problem(PPL_MIP_Problem &cpy_from)
        PPL_MIP_Problem(PPL_dimension_type dim) except +ValueError
        PPL_MIP_Problem(PPL_dimension_type dim, PPL_Constraint_System &cs, PPL_Linear_Expression &obj, PPL_Optimization_Mode) except +ValueError
        PPL_dimension_type space_dimension()
        PPL_Linear_Expression& objective_function()
        void clear()
        void add_space_dimensions_and_embed(PPL_dimension_type m) except +ValueError
        void add_constraint(PPL_Constraint &c) except +ValueError
        void add_constraints(PPL_Constraint_System &cs) except +ValueError
        void set_objective_function(PPL_Linear_Expression &obj) except +ValueError
        void set_optimization_mode(PPL_Optimization_Mode mode)
        PPL_Optimization_Mode optimization_mode()
        bint is_satisfiable()
        PPL_MIP_Problem_Status solve()
        void evaluate_objective_function(PPL_Generator evaluating_point, PPL_Coefficient &num, PPL_Coefficient &den) except +ValueError
        PPL_Generator& feasible_point()
        PPL_Generator optimizing_point() except +ValueError
        void optimal_value(PPL_Coefficient &num, PPL_Coefficient &den) except +ValueError
        bint OK()
        PPL_MIP_Problem_Control_Parameter_Value get_control_parameter(PPL_MIP_Problem_Control_Parameter_Name name)
        void set_control_parameter(PPL_MIP_Problem_Control_Parameter_Value value)

cdef extern from "ppl.hh":
    PPL_Generator PPL_line          "Parma_Polyhedra_Library::line"          (PPL_Linear_Expression &e) except +ValueError
    PPL_Generator PPL_ray           "Parma_Polyhedra_Library::ray"           (PPL_Linear_Expression &e) except +ValueError
    PPL_Generator PPL_point         "Parma_Polyhedra_Library::point"         (PPL_Linear_Expression &e, PPL_Coefficient &d) except +ValueError
    PPL_Generator PPL_closure_point "Parma_Polyhedra_Library::closure_point" (PPL_Linear_Expression &e, PPL_Coefficient &d) except +ValueError


####################################################
# Cython does not support static methods; hack around
cdef extern from "ppl.hh":

    PPL_Poly_Gen_Relation PPL_Poly_Gen_Relation_nothing  "Parma_Polyhedra_Library::Poly_Gen_Relation::nothing"  ()
    PPL_Poly_Gen_Relation PPL_Poly_Gen_Relation_subsumes "Parma_Polyhedra_Library::Poly_Gen_Relation::subsumes" ()

    PPL_Poly_Con_Relation PPL_Poly_Con_Relation_nothing  "Parma_Polyhedra_Library::Poly_Con_Relation::nothing"  ()
    PPL_Poly_Con_Relation PPL_Poly_Con_Relation_is_disjoint "Parma_Polyhedra_Library::Poly_Con_Relation::is_disjoint" ()
    PPL_Poly_Con_Relation PPL_Poly_Con_Relation_strictly_intersects "Parma_Polyhedra_Library::Poly_Con_Relation::strictly_intersects" ()
    PPL_Poly_Con_Relation PPL_Poly_Con_Relation_is_included "Parma_Polyhedra_Library::Poly_Con_Relation::is_included" ()
    PPL_Poly_Con_Relation PPL_Poly_Con_Relation_saturates "Parma_Polyhedra_Library::Poly_Con_Relation::saturates" ()



####################################################
# Workaround for private constructors
cdef extern from "ppl_shim.hh":
    PPL_Poly_Gen_Relation* new_relation_with(PPL_Polyhedron &p, PPL_Generator &g) except +ValueError
    PPL_Poly_Con_Relation* new_relation_with(PPL_Polyhedron &p, PPL_Constraint &c) except +ValueError


### Forward declarations ###########################
cdef class _mutable_or_immutable(SageObject)
cdef class Variable(object)
cdef class Linear_Expression(object)
cdef class Generator(object)
cdef class Generator_System(_mutable_or_immutable)
cdef class Generator_System_iterator(object)
cdef class Constraint(object)
cdef class Constraint_System(_mutable_or_immutable)
cdef class Constraint_System_iterator(object)
cdef class Polyhedron(_mutable_or_immutable)
cdef class C_Polyhedron(Polyhedron)
cdef class NNC_Polyhedron(Polyhedron)
cdef class Poly_Gen_Relation(object)
cdef class Poly_Con_Relation(object)
cdef class MIP_Problem(_mutable_or_immutable)


####################################################
### _mutable_or_immutable ##########################
####################################################
cdef class _mutable_or_immutable(SageObject):
    r"""
    A base class for mutable or immutable objects.

    By default, any object is mutable. It can then be
    :meth:`set_immutable`.

    EXAMPLES::

        sage: from sage.libs.ppl import _mutable_or_immutable as ExampleObj
        sage: x = ExampleObj()
        sage: x.is_mutable()
        True
        sage: x.is_immutable()
        False
        sage: x.set_immutable()
        sage: x.is_mutable()
        False
    """

    cdef bint _is_mutable

    def __cinit__(self):
        """
        The Cython constructor.

        TESTS::

            sage: from sage.libs.ppl import _mutable_or_immutable as ExampleObj
            sage: x = ExampleObj()    # indirect doctest
            sage: x.is_mutable()
            True
        """
        self._is_mutable = True


    def set_immutable(self):
        """
        Make this object immutable.

        This operation cannot be undone.

        EXAMPLES::

            sage: from sage.libs.ppl import _mutable_or_immutable as ExampleObj
            sage: x = ExampleObj()
            sage: x.is_mutable()
            True
            sage: x.set_immutable()
            sage: x.is_mutable()
            False
        """
        self._is_mutable = False


    def is_mutable(self):
        """
        Return whether this object is mutable.

        The data members of the object can only be modified if the
        object is mutable.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: from sage.libs.ppl import _mutable_or_immutable as ExampleObj
            sage: x = ExampleObj()
            sage: x.is_mutable()
            True
        """
        return self._is_mutable


    def is_immutable(self):
        """
        Return whether this object is immutable.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: from sage.libs.ppl import _mutable_or_immutable as ExampleObj
            sage: x = ExampleObj()
            sage: x.is_immutable()
            False
        """
        return not self._is_mutable


    def assert_mutable(self, msg):
        r"""
        Raise ``ValueError`` if the object is not mutable.

        INPUT:

        - ``msg`` -- a string. The message to be returned together
          with the ``ValueError``

        OUTPUT:

        This method returns no output. A ``ValueError``` is raised if
        the object is not mutable.

        EXAMPLES::

            sage: from sage.libs.ppl import _mutable_or_immutable as ExampleObj
            sage: x = ExampleObj()
            sage: x.assert_mutable("this will not trigger")
            sage: x.set_immutable()
            sage: x.assert_mutable("this will trigger")
            Traceback (most recent call last):
            ...
            ValueError: this will trigger
        """
        if not self._is_mutable:
            raise ValueError(msg)

####################################################
### MIP_Problem ####################################
####################################################
cdef class MIP_Problem(_mutable_or_immutable):
    r"""
    wrapper for PPL's MIP_Problem class

    An object of the class MIP_Problem represents a Mixed Integer
    (Linear) Program problem.

    INPUT:

    - ``dim`` -- integer
    - ``args`` -- an array of the defining data of the MIP_Problem.
      For each element, any one of the following is accepted:

      * A :class:`Constraint_System`.

      * A :class:`Linear_Expression`.

    OUTPUT:

    A :class:`MIP_Problem`.

    EXAMPLES::

        sage: from sage.libs.ppl import Variable, Constraint_System, MIP_Problem
        sage: x = Variable(0)
        sage: y = Variable(1)
        sage: cs = Constraint_System()
        sage: cs.insert( x >= 0)
        sage: cs.insert( y >= 0 )
        sage: cs.insert( 3 * x + 5 * y <= 10 )
        sage: m = MIP_Problem(2, cs, x + y)
        sage: m.optimal_value()
        10/3
        sage: m.optimizing_point()
        point(10/3, 0/3)
    """
    cdef PPL_MIP_Problem *thisptr

    def __repr__(self):
        """
        String representation of MIP Problem.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, Constraint_System, MIP_Problem
            sage: x = Variable(0)
            sage: y = Variable(1)
            sage: cs = Constraint_System()
            sage: cs.insert( x >= 0 )
            sage: cs.insert( y >= 0 )
            sage: cs.insert( 3 * x + 5 * y <= 10 )
            sage: m = MIP_Problem(2, cs, x + y)
            sage: m
            A MIP_Problem
            Maximize: x0+x1
            Subject to constraints
        """
        ret = 'A MIP_Problem\n'
        if self.optimization_mode() == 'maximization':
            ret += 'Maximize'
        else:
            ret += 'Minimize'
        ret += ': ' + str(self.objective_function()) + '\n'
        ret += 'Subject to constraints\n'

        return ret

    def __cinit__(self, PPL_dimension_type dim = 0, *args):
        """
        The Cython constructor.

        TESTS::

            sage: from sage.libs.ppl import Variable, Constraint_System, MIP_Problem
            sage: MIP_Problem(0)
            A MIP_Problem
            Maximize: 0
            Subject to constraints
        """
        if len(args) == 0:
            self.thisptr = new PPL_MIP_Problem(dim)
        elif len(args) == 2:
            cs = <Constraint_System>args[0]
            obj = <Linear_Expression>args[1]
            self.thisptr = new PPL_MIP_Problem(dim, cs.thisptr[0], obj.thisptr[0], MAXIMIZATION)
        elif len(args) == 3:
            cs = <Constraint_System>args[0]
            obj = <Linear_Expression>args[1]

            mode = str(args[2])
            if mode == 'maximization':
                self.thisptr = new PPL_MIP_Problem(dim, cs.thisptr[0], obj.thisptr[0], MAXIMIZATION)
            elif mode == 'minimization':
                self.thisptr = new PPL_MIP_Problem(dim, cs.thisptr[0], obj.thisptr[0], MINIMIZATION)
            else:
                raise ValueError, 'Unknown value: mode='+str(mode)+'.'
        else:
            raise ValueError, 'Cannot initialize with '+str(args)+'.'

    def __dealloc__(self):
        """
        The Cython destructor
        """
        del self.thisptr

    def optimization_mode(self):
        """
        Return the optimization mode used in the MIP_Problem.

        It will return "maximization" if the MIP_Problem was set
        to MAXIMIZATION mode, and "minimization" otherwise.

        EXAMPLES::

            sage: from sage.libs.ppl import MIP_Problem
            sage: m = MIP_Problem()
            sage: m.optimization_mode()
            'maximization'
        """
        if self.thisptr.optimization_mode() == MAXIMIZATION:
            return "maximization"
        elif self.thisptr.optimization_mode() == MINIMIZATION:
            return "minimization"

    def optimal_value(self):
        """
        Return the optimal value of the MIP_Problem. ValueError thrown if self does not
        have an optimizing point, i.e., if the MIP problem is unbounded or not satisfiable.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, Constraint_System, MIP_Problem
            sage: x = Variable(0)
            sage: y = Variable(1)
            sage: cs = Constraint_System()
            sage: cs.insert( x >= 0 )
            sage: cs.insert( y >= 0 )
            sage: cs.insert( 3 * x + 5 * y <= 10 )
            sage: m = MIP_Problem(2, cs, x + y)
            sage: m.optimal_value()
            10/3
            sage: cs = Constraint_System()
            sage: cs.insert( x >= 0 )
            sage: m = MIP_Problem(1, cs, x + x )
            sage: m.optimal_value()
            Traceback (most recent call last):
            ...
            ValueError: PPL::MIP_Problem::optimizing_point():
            *this does not have an optimizing point.
        """
        cdef PPL_Coefficient sup_n
        cdef PPL_Coefficient sup_d

        sig_on()
        try:
            self.thisptr.optimal_value(sup_n, sup_d)
        finally:
            sig_off()

        cdef Integer Int_sup_n = Integer(0)
        mpz_set(Int_sup_n.value, sup_n.get_mpz_t())
        cdef Integer Int_sup_d = Integer(0)
        mpz_set(Int_sup_d.value, sup_d.get_mpz_t())

        return Rational((Int_sup_n, Int_sup_d))

    def space_dimension(self):
        """
        Return the space dimension of the MIP_Problem.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, Constraint_System, MIP_Problem
            sage: x = Variable(0)
            sage: y = Variable(1)
            sage: cs = Constraint_System()
            sage: cs.insert( x >= 0)
            sage: cs.insert( y >= 0 )
            sage: cs.insert( 3 * x + 5 * y <= 10 )
            sage: m = MIP_Problem(2, cs, x + y)
            sage: m.space_dimension()
            2
        """
        return self.thisptr.space_dimension()

    def objective_function(self):
        """
        Return the optimal value of the MIP_Problem.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, Constraint_System, MIP_Problem
            sage: x = Variable(0)
            sage: y = Variable(1)
            sage: cs = Constraint_System()
            sage: cs.insert( x >= 0)
            sage: cs.insert( y >= 0 )
            sage: cs.insert( 3 * x + 5 * y <= 10 )
            sage: m = MIP_Problem(2, cs, x + y)
            sage: m.objective_function()
            x0+x1
        """
        rc = Linear_Expression()
        rc.thisptr[0] = self.thisptr.objective_function()
        return rc

    def clear(self):
        """
        Reset the MIP_Problem to be equal to the trivial MIP_Problem.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, Constraint_System, MIP_Problem
            sage: x = Variable(0)
            sage: y = Variable(1)
            sage: cs = Constraint_System()
            sage: cs.insert( x >= 0)
            sage: cs.insert( y >= 0 )
            sage: cs.insert( 3 * x + 5 * y <= 10 )
            sage: m = MIP_Problem(2, cs, x + y)
            sage: m.objective_function()
            x0+x1
            sage: m.clear()
            sage: m.objective_function()
            0
        """
        self.thisptr.clear()

    def add_space_dimensions_and_embed(self, PPL_dimension_type m):
        """
        Adds m new space dimensions and embeds the old MIP problem in the new vector space.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, Constraint_System, MIP_Problem
            sage: x = Variable(0)
            sage: y = Variable(1)
            sage: cs = Constraint_System()
            sage: cs.insert( x >= 0)
            sage: cs.insert( y >= 0 )
            sage: cs.insert( 3 * x + 5 * y <= 10 )
            sage: m = MIP_Problem(2, cs, x + y)
            sage: m.add_space_dimensions_and_embed(5)
            sage: m.space_dimension()
            7
        """
        self.assert_mutable("The MIP_Problem is not mutable!");
        sig_on()
        self.thisptr.add_space_dimensions_and_embed(m)
        sig_off()

    def add_constraint(self, Constraint c):
        """
        Adds a copy of constraint c to the MIP problem.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, Constraint_System, MIP_Problem
            sage: x = Variable(0)
            sage: y = Variable(1)
            sage: m = MIP_Problem()
            sage: m.add_space_dimensions_and_embed(2)
            sage: m.add_constraint(x >= 0)
            sage: m.add_constraint(y >= 0)
            sage: m.add_constraint(3 * x + 5 * y <= 10)
            sage: m.set_objective_function(x + y)
            sage: m.optimal_value()
            10/3

        TESTS::

            sage: z = Variable(2)
            sage: m.add_constraint(z >= -3)
            Traceback (most recent call last):
            ...
            ValueError: PPL::MIP_Problem::add_constraint(c):
            c.space_dimension() == 3 exceeds this->space_dimension == 2.
        """
        self.assert_mutable("The MIP_Problem is not mutable!");
        sig_on()
        try:
            self.thisptr.add_constraint(c.thisptr[0])
        finally:
            sig_off()

    def add_constraints(self, Constraint_System cs):
        """
        Adds a copy of the constraints in cs to the MIP problem.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, Constraint_System, MIP_Problem
            sage: x = Variable(0)
            sage: y = Variable(1)
            sage: cs = Constraint_System()
            sage: cs.insert( x >= 0)
            sage: cs.insert( y >= 0 )
            sage: cs.insert( 3 * x + 5 * y <= 10 )
            sage: m = MIP_Problem(2)
            sage: m.set_objective_function(x + y)
            sage: m.add_constraints(cs)
            sage: m.optimal_value()
            10/3

        TESTS::

            sage: p = Variable(9)
            sage: cs.insert(p >= -3)
            sage: m.add_constraints(cs)
            Traceback (most recent call last):
            ...
            ValueError: PPL::MIP_Problem::add_constraints(cs):
            cs.space_dimension() == 10 exceeds this->space_dimension() == 2.
        """
        self.assert_mutable("The MIP_Problem is not mutable!");
        sig_on()
        try:
            self.thisptr.add_constraints(cs.thisptr[0])
        finally:
            sig_off()

    def set_objective_function(self, Linear_Expression obj):
        """
        Sets the objective function to obj.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, Constraint_System, MIP_Problem
            sage: x = Variable(0)
            sage: y = Variable(1)
            sage: m = MIP_Problem()
            sage: m.add_space_dimensions_and_embed(2)
            sage: m.add_constraint(x >= 0)
            sage: m.add_constraint(y >= 0)
            sage: m.add_constraint(3 * x + 5 * y <= 10)
            sage: m.set_objective_function(x + y)
            sage: m.optimal_value()
            10/3

        TESTS::

            sage: z = Variable(2)
            sage: m.set_objective_function(x + y + z)
            Traceback (most recent call last):
            ...
            ValueError: PPL::MIP_Problem::set_objective_function(obj):
            obj.space_dimension() == 3 exceeds this->space_dimension == 2.
        """
        self.assert_mutable("The MIP_Problem is not mutable!");
        self.thisptr.set_objective_function(obj.thisptr[0])

    def set_optimization_mode(self, mode):
        """
        Sets the optimization mode to mode.

        EXAMPLES::

            sage: from sage.libs.ppl import MIP_Problem
            sage: m = MIP_Problem()
            sage: m.optimization_mode()
            'maximization'
            sage: m.set_optimization_mode('minimization')
            sage: m.optimization_mode()
            'minimization'

        TESTS::

            sage: m.set_optimization_mode('max')
            Traceback (most recent call last):
            ...
            ValueError: Unknown value: mode=max.
        """
        if mode == 'minimization':
            self.thisptr.set_optimization_mode(MINIMIZATION)
        elif mode == 'maximization':
            self.thisptr.set_optimization_mode(MAXIMIZATION)
        else:
            raise ValueError, 'Unknown value: mode='+str(mode)+'.'

    def is_satisfiable(self):
        """
        Check if the MIP_Problem is satisfiable

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, Constraint_System, MIP_Problem
            sage: x = Variable(0)
            sage: y = Variable(1)
            sage: m = MIP_Problem()
            sage: m.add_space_dimensions_and_embed(2)
            sage: m.add_constraint(x >= 0)
            sage: m.add_constraint(y >= 0)
            sage: m.add_constraint(3 * x + 5 * y <= 10)
            sage: m.is_satisfiable()
            True
        """
        ret = self.thisptr.is_satisfiable()

        return ret

    def evaluate_objective_function(self, Generator evaluating_point):
        """
        Return the result of evaluating the objective function on evaluating_point. ValueError thrown
        if self and evaluating_point are dimension-incompatible or if the generator evaluating_point is not a point.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, Constraint_System, MIP_Problem, Generator
            sage: x = Variable(0)
            sage: y = Variable(1)
            sage: m = MIP_Problem()
            sage: m.add_space_dimensions_and_embed(2)
            sage: m.add_constraint(x >= 0)
            sage: m.add_constraint(y >= 0)
            sage: m.add_constraint(3 * x + 5 * y <= 10)
            sage: m.set_objective_function(x + y)
            sage: g = Generator.point(5 * x - 2 * y, 7)
            sage: m.evaluate_objective_function(g)
            3/7
            sage: z = Variable(2)
            sage: g = Generator.point(5 * x - 2 * z, 7)
            sage: m.evaluate_objective_function(g)
            Traceback (most recent call last):
            ...
            ValueError: PPL::MIP_Problem::evaluate_objective_function(p, n, d):
            *this and p are dimension incompatible.
        """
        cdef PPL_Coefficient sup_n
        cdef PPL_Coefficient sup_d

        sig_on()
        try:
            self.thisptr.evaluate_objective_function(evaluating_point.thisptr[0], sup_n, sup_d)
        finally:
            sig_off()

        cdef Integer Int_sup_n = Integer(0)
        mpz_set(Int_sup_n.value, sup_n.get_mpz_t())
        cdef Integer Int_sup_d = Integer(0)
        mpz_set(Int_sup_d.value, sup_d.get_mpz_t())

        return Rational((Int_sup_n, Int_sup_d))

    def solve(self):
        """
        Optimizes the MIP_Problem

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, Constraint_System, MIP_Problem
            sage: x = Variable(0)
            sage: y = Variable(1)
            sage: m = MIP_Problem()
            sage: m.add_space_dimensions_and_embed(2)
            sage: m.add_constraint(x >= 0)
            sage: m.add_constraint(y >= 0)
            sage: m.add_constraint(3 * x + 5 * y <= 10)
            sage: m.set_objective_function(x + y)
            sage: m.solve()
            {'status': 'optimized'}
        """
        sig_on()
        try:
            tmp = self.thisptr.solve()
        finally:
            sig_off()
        if tmp == UNFEASIBLE_MIP_PROBLEM:
            return { 'status':'unfeasible' }
        elif tmp == UNBOUNDED_MIP_PROBLEM:
            return { 'status':'unbounded' }
        else:
            return { 'status':'optimized' }

    def optimizing_point(self):
        """
        Returns an optimal point for the MIP_Problem, if it exists.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, Constraint_System, MIP_Problem
            sage: x = Variable(0)
            sage: y = Variable(1)
            sage: m = MIP_Problem()
            sage: m.add_space_dimensions_and_embed(2)
            sage: m.add_constraint(x >= 0)
            sage: m.add_constraint(y >= 0)
            sage: m.add_constraint(3 * x + 5 * y <= 10)
            sage: m.set_objective_function(x + y)
            sage: m.optimizing_point()
            point(10/3, 0/3)
        """
        cdef PPL_Generator *g
        sig_on()
        try:
            g = new_MIP_optimizing_point(self.thisptr[0])
        finally:
            sig_off()
        return _wrap_Generator(g[0])

    def OK(self):
        """
        Check if all the invariants are satisfied.

        OUTPUT:

        ``True`` if and only if ``self`` satisfies all the invariants.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, Constraint_System, MIP_Problem
            sage: x = Variable(0)
            sage: y = Variable(1)
            sage: m = MIP_Problem()
            sage: m.add_space_dimensions_and_embed(2)
            sage: m.add_constraint(x >= 0)
            sage: m.OK()
            True
        """
        return self.thisptr.OK()

####################################################
### Polyhedron #####################################
####################################################
cdef class Polyhedron(_mutable_or_immutable):
    r"""
    Wrapper for PPL's ``Polyhedron`` class.

    An object of the class Polyhedron represents a convex polyhedron
    in the vector space.

    A polyhedron can be specified as either a finite system of
    constraints or a finite system of generators (see Section
    Representations of Convex Polyhedra) and it is always possible to
    obtain either representation. That is, if we know the system of
    constraints, we can obtain from this the system of generators that
    define the same polyhedron and vice versa. These systems can
    contain redundant members: in this case we say that they are not
    in the minimal form.

    INPUT/OUTPUT:

    This is an abstract base for :class:`C_Polyhedron` and
    :class:`NNC_Polyhedron`. You cannot instantiate this class.
    """

    cdef PPL_Polyhedron *thisptr


    def __init__(self):
        r"""
        The Python constructor.

        See also :class:`C_Polyhedron` and
        :class:`NNC_Polyhedron`. You must not instantiate
        :class:`Polyhedron` objects.

        TESTS::

            sage: from sage.libs.ppl import Polyhedron
            sage: Polyhedron()
            Traceback (most recent call last):
            ...
            NotImplementedError: The Polyhedron class is abstract, you must not instantiate it.
        """
        raise NotImplementedError, 'The Polyhedron class is abstract, you must not instantiate it.'


    def _repr_(self):
        """
        Return a string representation.

        OUTPUT:

        String.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, C_Polyhedron
            sage: x = Variable(0)
            sage: y = Variable(1)
            sage: C_Polyhedron( 5*x-2*y >=  x+y-1 )._repr_()
            'A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 point, 1 ray, 1 line'

        Special cases::

            sage: C_Polyhedron(3, 'empty')._repr_()
            'The empty polyhedron in QQ^3'
            sage: C_Polyhedron(3, 'universe')._repr_()
            'The space-filling polyhedron in QQ^3'
        """
        dim = self.affine_dimension()
        ambient_dim = self.space_dimension()
        gs = self.minimized_generators()
        n_points = 0
        n_closure_points = 0
        n_lines = 0
        n_rays = 0
        for g in gs:
            if g.is_line():
                n_lines += 1
            elif g.is_ray():
                n_rays += 1
            elif g.is_point():
                n_points += 1
            elif g.is_closure_point():
                n_closure_points += 1
            else:
                assert False
        if self.is_empty():
            return 'The empty polyhedron in QQ^'+str(ambient_dim)
        if self.is_universe():
            return 'The space-filling polyhedron in QQ^'+str(ambient_dim)
        desc = 'A ' + str(dim) + '-dimensional polyhedron'
        desc += ' in QQ'
        desc += '^' + str(ambient_dim)
        desc += ' defined as the convex hull of '
        first = True
        if n_points>0:
            if not first:
                desc += ", "
            first = False
            desc += str(n_points)
            if n_points==1:  desc += ' point'
            else:          desc += ' points'
        if n_closure_points>0:
            if not first:
                desc += ", "
            first = False
            desc += str(n_closure_points)
            if n_closure_points==1:  desc += ' closure_point'
            else:          desc += ' closure_points'
        if n_rays>0:
            if not first:
                desc += ", "
            first = False
            desc += str(n_rays)
            if n_rays==1:  desc += ' ray'
            else:          desc += ' rays'
        if n_lines>0:
            if not first:
                desc += ", "
            first = False
            desc += repr(n_lines)
            if n_lines==1: desc +=' line'
            else:          desc +=' lines'
        return desc;


    def space_dimension(self):
        r"""
        Return the dimension of the vector space enclosing ``self``.

        OUTPUT:

        Integer.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, C_Polyhedron
            sage: x = Variable(0)
            sage: y = Variable(1)
            sage: p = C_Polyhedron( 5*x-2*y >=  x+y-1 )
            sage: p.space_dimension()
            2
        """
        return self.thisptr.space_dimension()


    def affine_dimension(self):
        r"""
        Return the affine dimension of ``self``.

        OUTPUT:

        An integer. Returns 0 if ``self`` is empty. Otherwise, returns
        the affine dimension of ``self``.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, C_Polyhedron
            sage: x = Variable(0)
            sage: y = Variable(1)
            sage: p = C_Polyhedron( 5*x-2*y ==  x+y-1 )
            sage: p.affine_dimension()
            1
        """
        sig_on()
        cdef size_t dim = self.thisptr.affine_dimension()
        sig_off()
        return dim


    def constraints(self):
        r"""
        Returns the system of constraints.

        See also :meth:`minimized_constraints`.

        OUTPUT:

        A :class:`Constraint_System`.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, C_Polyhedron
            sage: x = Variable(0)
            sage: y = Variable(1)
            sage: p = C_Polyhedron( y>=0 )
            sage: p.add_constraint( x>=0 )
            sage: p.add_constraint( x+y>=0 )
            sage: p.constraints()
            Constraint_System {x1>=0, x0>=0, x0+x1>=0}
            sage: p.minimized_constraints()
            Constraint_System {x1>=0, x0>=0}
        """
        sig_on()
        cdef PPL_Constraint_System cs = self.thisptr.constraints()
        sig_off()
        return _wrap_Constraint_System(cs)


    def minimized_constraints(self):
        r"""
        Returns the minimized system of constraints.

        See also :meth:`constraints`.

        OUTPUT:

        A :class:`Constraint_System`.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, C_Polyhedron
            sage: x = Variable(0)
            sage: y = Variable(1)
            sage: p = C_Polyhedron( y>=0 )
            sage: p.add_constraint( x>=0 )
            sage: p.add_constraint( x+y>=0 )
            sage: p.constraints()
            Constraint_System {x1>=0, x0>=0, x0+x1>=0}
            sage: p.minimized_constraints()
            Constraint_System {x1>=0, x0>=0}
        """
        sig_on()
        cdef PPL_Constraint_System cs = self.thisptr.minimized_constraints()
        sig_off()
        return _wrap_Constraint_System(cs)


    def generators(self):
        r"""
        Returns the system of generators.

        See also :meth:`minimized_generators`.

        OUTPUT:

        A :class:`Generator_System`.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, C_Polyhedron, point
            sage: x = Variable(0)
            sage: y = Variable(1)
            sage: p = C_Polyhedron(3,'empty')
            sage: p.add_generator( point(-x-y) )
            sage: p.add_generator( point(0) )
            sage: p.add_generator( point(+x+y) )
            sage: p.generators()
            Generator_System {point(-1/1, -1/1, 0/1), point(0/1, 0/1, 0/1), point(1/1, 1/1, 0/1)}
            sage: p.minimized_generators()
            Generator_System {point(-1/1, -1/1, 0/1), point(1/1, 1/1, 0/1)}
        """
        sig_on()
        cdef PPL_Generator_System gs = self.thisptr.generators()
        sig_off()
        return _wrap_Generator_System(gs)


    def minimized_generators(self):
        r"""
        Returns the minimized system of generators.

        See also :meth:`generators`.

        OUTPUT:

        A :class:`Generator_System`.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, C_Polyhedron, point
            sage: x = Variable(0)
            sage: y = Variable(1)
            sage: p = C_Polyhedron(3,'empty')
            sage: p.add_generator( point(-x-y) )
            sage: p.add_generator( point(0) )
            sage: p.add_generator( point(+x+y) )
            sage: p.generators()
            Generator_System {point(-1/1, -1/1, 0/1), point(0/1, 0/1, 0/1), point(1/1, 1/1, 0/1)}
            sage: p.minimized_generators()
            Generator_System {point(-1/1, -1/1, 0/1), point(1/1, 1/1, 0/1)}
        """
        sig_on()
        cdef PPL_Generator_System gs = self.thisptr.minimized_generators()
        sig_off()
        return _wrap_Generator_System(gs)


    cdef _relation_with_generator(Polyhedron self, Generator g):
        r"""
        Helper method for :meth:`relation_with`.
        """
        rel = Poly_Gen_Relation(True)
        try:
            sig_on()
            try:
                rel.thisptr = new_relation_with(self.thisptr[0], g.thisptr[0])
            finally:
                sig_off()
        except BaseException:
            # rel.thisptr must be set to something valid or rel.__dealloc__() will segfault
            rel.thisptr = new PPL_Poly_Gen_Relation(PPL_Poly_Gen_Relation_nothing())
            raise
        return rel


    cdef _relation_with_constraint(Polyhedron self, Constraint c):
        r"""
        Helper method for :meth:`relation_with`.
        """
        rel = Poly_Con_Relation(True)
        try:
            sig_on()
            try:
                rel.thisptr = new_relation_with(self.thisptr[0], c.thisptr[0])
            finally:
                sig_off()
        except BaseException:
            # rel.thisptr must be set to something valid or rel.__dealloc__() will segfault
            rel.thisptr = new PPL_Poly_Con_Relation(PPL_Poly_Con_Relation_nothing())
            raise
        return rel


    def relation_with(self, arg):
        r"""
        Return the relations holding between the polyhedron ``self``
        and the generator or constraint ``arg``.

        INPUT:

        - ``arg`` -- a :class:`Generator` or a :class:`Constraint`.

        OUTPUT:

        A :class:`Poly_Gen_Relation` or a :class:`Poly_Con_Relation`
        according to the type of the input.

        Raises ``ValueError`` if ``self`` and the generator/constraint
        ``arg`` are dimension-incompatible.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, C_Polyhedron, point, ray, Poly_Con_Relation
            sage: x = Variable(0);  y = Variable(1)
            sage: p = C_Polyhedron(2, 'empty')
            sage: p.add_generator( point(1*x+0*y) )
            sage: p.add_generator( point(0*x+1*y) )
            sage: p.minimized_constraints()
            Constraint_System {x0+x1-1==0, -x1+1>=0, x1>=0}
            sage: p.relation_with( point(1*x+1*y) )
            nothing
            sage: p.relation_with( point(1*x+1*y, 2) )
            subsumes
            sage: p.relation_with( x+y==-1 )
            is_disjoint
            sage: p.relation_with( x==y )
            strictly_intersects
            sage: p.relation_with( x+y<=1 )
            is_included, saturates
            sage: p.relation_with( x+y<1 )
            is_disjoint, saturates

        In a Sage program you will usually use :meth:`relation_with`
        together with :meth:`~sage.libs.ppl.Poly_Gen_Relation.implies`
        or :meth:`~sage.libs.ppl.Poly_Con_Relation.implies`, for
        example::

            sage: p.relation_with( x+y<1 ).implies(Poly_Con_Relation.saturates())
            True

        You can only get relations with dimension-compatible
        generators or constraints::

            sage: z = Variable(2)
            sage: p.relation_with( point(x+y+z) )
            Traceback (most recent call last):
            ...
            ValueError: PPL::C_Polyhedron::relation_with(g):
            this->space_dimension() == 2, g.space_dimension() == 3.
            sage: p.relation_with( z>0 )
            Traceback (most recent call last):
            ...
            ValueError: PPL::C_Polyhedron::relation_with(c):
            this->space_dimension() == 2, c.space_dimension() == 3.
        """
        if isinstance(arg, Generator):
            return self._relation_with_generator(arg)
        if isinstance(arg, Constraint):
            return self._relation_with_constraint(arg)
        else:
            raise TypeError, 'Argument must be Generator or a Constraint'


    def is_empty(self):
        """
        Test if ``self`` is an empty polyhedron.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: from sage.libs.ppl import C_Polyhedron
            sage: C_Polyhedron(3, 'empty').is_empty()
            True
            sage: C_Polyhedron(3, 'universe').is_empty()
            False
        """
        sig_on()
        cdef bint result = self.thisptr.is_empty()
        sig_off()
        return result


    def is_universe(self):
        """
        Test if ``self`` is a universe (space-filling) polyhedron.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: from sage.libs.ppl import C_Polyhedron
            sage: C_Polyhedron(3, 'empty').is_universe()
            False
            sage: C_Polyhedron(3, 'universe').is_universe()
            True
        """
        sig_on()
        cdef bint result = self.thisptr.is_universe()
        sig_off()
        return result


    def is_topologically_closed(self):
        """
        Tests if ``self`` is topologically closed.

        OUTPUT:

        Returns ``True`` if and only if ``self`` is a topologically
        closed subset of the ambient vector space.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, C_Polyhedron, NNC_Polyhedron
            sage: x = Variable(0);  y = Variable(1)
            sage: C_Polyhedron(3, 'universe').is_topologically_closed()
            True
            sage: C_Polyhedron( x>=1 ).is_topologically_closed()
            True
            sage: NNC_Polyhedron( x>1 ).is_topologically_closed()
            False
        """
        sig_on()
        cdef bint result = self.thisptr.is_topologically_closed()
        sig_off()
        return result


    def is_disjoint_from(self, Polyhedron y):
        r"""
        Tests whether ``self`` and ``y`` are disjoint.

        INPUT:

        - ``y`` -- a :class:`Polyhedron`.

        OUTPUT:

        Boolean. Returns ``True`` if and only if ``self`` and ``y``
        are disjoint.

        Rayises a ``ValueError`` if ``self`` and ``y`` are
        topology-incompatible or dimension-incompatible.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, C_Polyhedron, NNC_Polyhedron
            sage: x = Variable(0);  y = Variable(1)
            sage: C_Polyhedron(x<=0).is_disjoint_from( C_Polyhedron(x>=1) )
            True

        This is not allowed::

            sage: x = Variable(0);  y = Variable(1)
            sage: poly_1d = C_Polyhedron(x<=0)
            sage: poly_2d = C_Polyhedron(x+0*y>=1)
            sage: poly_1d.is_disjoint_from(poly_2d)
            Traceback (most recent call last):
            ...
            ValueError: PPL::C_Polyhedron::intersection_assign(y):
            this->space_dimension() == 1, y.space_dimension() == 2.

        Nor is this::

            sage: x = Variable(0);  y = Variable(1)
            sage: c_poly   =   C_Polyhedron( x<=0 )
            sage: nnc_poly = NNC_Polyhedron( x >0 )
            sage: c_poly.is_disjoint_from(nnc_poly)
            Traceback (most recent call last):
            ...
            ValueError: PPL::C_Polyhedron::intersection_assign(y):
            y is a NNC_Polyhedron.
            sage: NNC_Polyhedron(c_poly).is_disjoint_from(nnc_poly)
            True
        """
        cdef bint result
        sig_on()
        try:
            result = self.thisptr.is_disjoint_from(y.thisptr[0])
        finally:
            sig_off()
        return result


    def is_discrete(self):
        r"""
        Test whether ``self`` is discrete.

        OUTPUT:

        Boolean. Returns ``True`` if and only if ``self`` is discrete.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, C_Polyhedron, point, ray
            sage: x = Variable(0);  y = Variable(1)
            sage: p = C_Polyhedron( point(1*x+2*y) )
            sage: p.is_discrete()
            True
            sage: p.add_generator( point(x) )
            sage: p.is_discrete()
            False
        """
        sig_on()
        cdef bint result = self.thisptr.is_discrete()
        sig_off()
        return result


    def is_bounded(self):
        r"""
        Test whether ``self`` is bounded.

        OUTPUT:

        Boolean. Returns ``True`` if and only if ``self`` is a bounded polyhedron.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, NNC_Polyhedron, point, closure_point, ray
            sage: x = Variable(0)
            sage: p = NNC_Polyhedron( point(0*x) )
            sage: p.add_generator( closure_point(1*x) )
            sage: p.is_bounded()
            True
            sage: p.add_generator( ray(1*x) )
            sage: p.is_bounded()
            False
        """
        sig_on()
        cdef bint result = self.thisptr.is_bounded()
        sig_off()
        return result


    def contains_integer_point(self):
        r"""
        Test whether ``self`` contains an integer point.

        OUTPUT:

        Boolean. Returns ``True`` if and only if ``self`` contains an
        integer point.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, NNC_Polyhedron
            sage: x = Variable(0)
            sage: p = NNC_Polyhedron(x>0)
            sage: p.add_constraint(x<1)
            sage: p.contains_integer_point()
            False
            sage: p.topological_closure_assign()
            sage: p.contains_integer_point()
            True
        """
        sig_on()
        cdef bint result = self.thisptr.contains_integer_point()
        sig_off()
        return result


    def constrains(self, Variable var):
        r"""
        Test whether ``var`` is constrained in ``self``.

        INPUT:

        - ``var`` -- a :class:`Variable`.

        OUTPUT:

        Boolean. Returns ``True`` if and only if ``var`` is
        constrained in ``self``.

        Raises a ``ValueError`` if ``var`` is not a space dimension of
        ``self``.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, C_Polyhedron
            sage: x = Variable(0)
            sage: p = C_Polyhedron(1, 'universe')
            sage: p.constrains(x)
            False
            sage: p = C_Polyhedron(x>=0)
            sage: p.constrains(x)
            True
            sage: y = Variable(1)
            sage: p.constrains(y)
            Traceback (most recent call last):
            ...
            ValueError: PPL::C_Polyhedron::constrains(v):
            this->space_dimension() == 1, v.space_dimension() == 2.
        """
        cdef bint result
        sig_on()
        try:
            result = self.thisptr.constrains(var.thisptr[0])
        finally:
            sig_off()
        return result


    def bounds_from_above(self, Linear_Expression expr):
        r"""
        Test whether the ``expr`` is bounded from above.

        INPUT:

        - ``expr`` -- a :class:`Linear_Expression`

        OUTPUT:

        Boolean. Returns ``True`` if and only if ``expr`` is bounded
        from above in ``self``.

        Raises a ``ValueError`` if ``expr`` and ``this`` are
        dimension-incompatible.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, C_Polyhedron, Linear_Expression
            sage: x = Variable(0);  y = Variable(1)
            sage: p = C_Polyhedron(y<=0)
            sage: p.bounds_from_above(x+1)
            False
            sage: p.bounds_from_above(Linear_Expression(y))
            True
            sage: p = C_Polyhedron(x<=0)
            sage: p.bounds_from_above(y+1)
            Traceback (most recent call last):
            ...
            ValueError: PPL::C_Polyhedron::bounds_from_above(e):
            this->space_dimension() == 1, e.space_dimension() == 2.
        """
        cdef bint result
        sig_on()
        try:
            result = self.thisptr.bounds_from_above(expr.thisptr[0])
        finally:
            sig_off()
        return result


    def bounds_from_below(self, Linear_Expression expr):
        r"""
        Test whether the ``expr`` is bounded from above.

        INPUT:

        - ``expr`` -- a :class:`Linear_Expression`

        OUTPUT:

        Boolean. Returns ``True`` if and only if ``expr`` is bounded
        from above in ``self``.

        Raises a ``ValueError`` if ``expr`` and ``this`` are
        dimension-incompatible.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, C_Polyhedron, Linear_Expression
            sage: x = Variable(0);  y = Variable(1)
            sage: p = C_Polyhedron(y>=0)
            sage: p.bounds_from_below(x+1)
            False
            sage: p.bounds_from_below(Linear_Expression(y))
            True
            sage: p = C_Polyhedron(x<=0)
            sage: p.bounds_from_below(y+1)
            Traceback (most recent call last):
            ...
            ValueError: PPL::C_Polyhedron::bounds_from_below(e):
            this->space_dimension() == 1, e.space_dimension() == 2.
        """
        cdef bint result
        sig_on()
        try:
            result = self.thisptr.bounds_from_below(expr.thisptr[0])
        finally:
            sig_off()
        return result


    def maximize(self, Linear_Expression expr):
        r"""
        Maximize ``expr``.

        INPUT:

        - ``expr`` -- a :class:`Linear_Expression`.

        OUTPUT:

        A dictionary with the following keyword:value pair:

        * ``'bounded'``: Boolean. Whether the linear expression
          ``expr`` is bounded from above on ``self``.

        If ``expr`` is bounded from above, the following additional
        keyword:value pairs are set to provide information about the
        supremum:

        * ``'sup_n'``: Integer. The numerator of the supremum value.

        * ``'sup_d'``: Non-zero integer. The denominator of the supremum
          value.

        * ``'maximum'``: Boolean. ``True`` if and only if the supremum
          is also the maximum value.

        * ``'generator'``: a :class:`Generator`. A point or closure
          point where expr reaches its supremum value.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, C_Polyhedron, NNC_Polyhedron, Constraint_System, Linear_Expression
            sage: x = Variable(0);  y = Variable(1)
            sage: cs = Constraint_System()
            sage: cs.insert( x>=0 )
            sage: cs.insert( y>=0 )
            sage: cs.insert( 3*x+5*y<=10 )
            sage: p = C_Polyhedron(cs)
            sage: p.maximize( x+y )
            {'sup_d': 3, 'sup_n': 10, 'bounded': True, 'maximum': True, 'generator': point(10/3, 0/3)}

        Unbounded case::

            sage: cs = Constraint_System()
            sage: cs.insert( x>0 )
            sage: p = NNC_Polyhedron(cs)
            sage: p.maximize( +x )
            {'bounded': False}
            sage: p.maximize( -x )
            {'sup_d': 1, 'sup_n': 0, 'bounded': True, 'maximum': False, 'generator': closure_point(0/1)}
        """
        cdef PPL_Coefficient sup_n
        cdef PPL_Coefficient sup_d
        cdef Generator g = point()
        cdef cppbool maximum
        sig_on()
        rc = self.thisptr.maximize(<PPL_Linear_Expression&>expr.thisptr[0], sup_n, sup_d, maximum, g.thisptr[0])
        sig_off()

        cdef Integer Int_sup_n = Integer(0)
        mpz_set(Int_sup_n.value, sup_n.get_mpz_t())
        cdef Integer Int_sup_d = Integer(0)
        mpz_set(Int_sup_d.value, sup_d.get_mpz_t())

        if rc:
            return { 'bounded':True, 'sup_n':Int_sup_n, 'sup_d':Int_sup_d, 'maximum':maximum, 'generator':g }
        else:
            return { 'bounded':False }


    def minimize(self, Linear_Expression expr):
        r"""
        Minimize ``expr``.

        INPUT:

        - ``expr`` -- a :class:`Linear_Expression`.

        OUTPUT:

        A dictionary with the following keyword:value pair:

        * ``'bounded'``: Boolean. Whether the linear expression
          ``expr`` is bounded from below on ``self``.

        If ``expr`` is bounded from below, the following additional
        keyword:value pairs are set to provide information about the
        infimum:

        * ``'inf_n'``: Integer. The numerator of the infimum value.

        * ``'inf_d'``: Non-zero integer. The denominator of the infimum
          value.

        * ``'minimum'``: Boolean. ``True`` if and only if the infimum
          is also the minimum value.

        * ``'generator'``: a :class:`Generator`. A point or closure
          point where expr reaches its infimum value.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, C_Polyhedron, NNC_Polyhedron, Constraint_System, Linear_Expression
            sage: x = Variable(0);  y = Variable(1)
            sage: cs = Constraint_System()
            sage: cs.insert( x>=0 )
            sage: cs.insert( y>=0 )
            sage: cs.insert( 3*x+5*y<=10 )
            sage: p = C_Polyhedron(cs)
            sage: p.minimize( x+y )
            {'minimum': True, 'bounded': True, 'inf_d': 1, 'generator': point(0/1, 0/1), 'inf_n': 0}

        Unbounded case::

            sage: cs = Constraint_System()
            sage: cs.insert( x>0 )
            sage: p = NNC_Polyhedron(cs)
            sage: p.minimize( +x )
            {'minimum': False, 'bounded': True, 'inf_d': 1, 'generator': closure_point(0/1), 'inf_n': 0}
            sage: p.minimize( -x )
            {'bounded': False}
        """
        cdef PPL_Coefficient inf_n
        cdef PPL_Coefficient inf_d
        cdef Generator g = point()
        cdef cppbool minimum
        sig_on()
        rc = self.thisptr.minimize(<PPL_Linear_Expression&>expr.thisptr[0], inf_n, inf_d, minimum, g.thisptr[0])
        sig_off()

        cdef Integer Int_inf_n = Integer(0)
        mpz_set(Int_inf_n.value, inf_n.get_mpz_t())
        cdef Integer Int_inf_d = Integer(0)
        mpz_set(Int_inf_d.value, inf_d.get_mpz_t())

        if rc:
            return { 'bounded':True, 'inf_n':Int_inf_n, 'inf_d':Int_inf_d, 'minimum':minimum, 'generator':g }
        else:
            return { 'bounded':False }


    def contains(self, Polyhedron y):
        r"""
        Test whether ``self`` contains ``y``.

        INPUT:

        - ``y`` -- a :class:`Polyhedron`.

        OUTPUT:

        Boolean. Returns ``True`` if and only if ``self`` contains ``y``.

        Raises a ``ValueError`` if ``self`` and ``y`` are
        topology-incompatible or dimension-incompatible.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, C_Polyhedron, NNC_Polyhedron
            sage: x = Variable(0)
            sage: y = Variable(1)
            sage: p0 = C_Polyhedron( x>=0 )
            sage: p1 = C_Polyhedron( x>=1 )
            sage: p0.contains(p1)
            True
            sage: p1.contains(p0)
            False

        Errors are raised if the dimension or topology is not compatible::

            sage: p0.contains(C_Polyhedron(y>=0))
            Traceback (most recent call last):
            ...
            ValueError: PPL::C_Polyhedron::contains(y):
            this->space_dimension() == 1, y.space_dimension() == 2.
            sage: p0.contains(NNC_Polyhedron(x>0))
            Traceback (most recent call last):
            ...
            ValueError: PPL::C_Polyhedron::contains(y):
            y is a NNC_Polyhedron.
        """
        cdef bint result
        sig_on()
        try:
            result = self.thisptr.contains(y.thisptr[0])
        finally:
            sig_off()
        return result


    def strictly_contains(self, Polyhedron y):
        r"""
        Test whether ``self`` strictly contains ``y``.

        INPUT:

        - ``y`` -- a :class:`Polyhedron`.

        OUTPUT:

        Boolean. Returns ``True`` if and only if ``self`` contains
        ``y`` and ``self`` does not equal ``y``.

        Raises a ``ValueError`` if ``self`` and ``y`` are
        topology-incompatible or dimension-incompatible.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, C_Polyhedron, NNC_Polyhedron
            sage: x = Variable(0)
            sage: y = Variable(1)
            sage: p0 = C_Polyhedron( x>=0 )
            sage: p1 = C_Polyhedron( x>=1 )
            sage: p0.strictly_contains(p1)
            True
            sage: p1.strictly_contains(p0)
            False

        Errors are raised if the dimension or topology is not compatible::

            sage: p0.strictly_contains(C_Polyhedron(y>=0))
            Traceback (most recent call last):
            ...
            ValueError: PPL::C_Polyhedron::contains(y):
            this->space_dimension() == 1, y.space_dimension() == 2.
            sage: p0.strictly_contains(NNC_Polyhedron(x>0))
            Traceback (most recent call last):
            ...
            ValueError: PPL::C_Polyhedron::contains(y):
            y is a NNC_Polyhedron.
        """
        cdef bint result
        sig_on()
        try:
            result = self.thisptr.strictly_contains(y.thisptr[0])
        finally:
            sig_off()
        return result


    def add_constraint(self, Constraint c):
        r"""
        Add a constraint to the polyhedron.

        Adds a copy of constraint ``c`` to the system of constraints
        of ``self``, without minimizing the result.

        See alse :meth:`add_constraints`.

        INPUT:

        - ``c`` -- the :class:`Constraint` that will be added to the
          system of constraints of ``self``.

        OUTPUT:

        This method modifies the polyhedron ``self`` and does not
        return anything.

        Raises a ``ValueError`` if ``self`` and the constraint ``c`` are
        topology-incompatible or dimension-incompatible.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, C_Polyhedron
            sage: x = Variable(0)
            sage: y = Variable(1)
            sage: p = C_Polyhedron( y>=0 )
            sage: p.add_constraint( x>=0 )

         We just added a 1-d constraint to a 2-d polyhedron, this is
         fine. The other way is not::

            sage: p = C_Polyhedron( x>=0 )
            sage: p.add_constraint( y>=0 )
            Traceback (most recent call last):
            ...
            ValueError: PPL::C_Polyhedron::add_constraint(c):
            this->space_dimension() == 1, c.space_dimension() == 2.

         The constraint must also be topology-compatible, that is,
         :class:`C_Polyhedron` only allows non-strict inequalities::

            sage: p = C_Polyhedron( x>=0 )
            sage: p.add_constraint( x< 1 )
            Traceback (most recent call last):
            ...
            ValueError: PPL::C_Polyhedron::add_constraint(c):
            c is a strict inequality.
        """
        self.assert_mutable('The Polyhedron is not mutable!')
        sig_on()
        try:
            self.thisptr.add_constraint(c.thisptr[0])
        finally:
            sig_off()


    def add_generator(self, Generator g):
        r"""
        Add a generator to the polyhedron.

        Adds a copy of constraint ``c`` to the system of generators
        of ``self``, without minimizing the result.

        INPUT:

        - ``g`` -- the :class:`Generator` that will be added to the
          system of Generators of ``self``.

        OUTPUT:

        This method modifies the polyhedron ``self`` and does not
        return anything.

        Raises a ``ValueError`` if ``self`` and the generator ``g``
        are topology-incompatible or dimension-incompatible, or if
        ``self`` is an empty polyhedron and ``g`` is not a point.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, C_Polyhedron, point, closure_point, ray
            sage: x = Variable(0)
            sage: y = Variable(1)
            sage: p = C_Polyhedron(1, 'empty')
            sage: p.add_generator( point(0*x) )

         We just added a 1-d generator to a 2-d polyhedron, this is
         fine. The other way is not::

            sage: p = C_Polyhedron(1, 'empty')
            sage: p.add_generator(  point(0*y) )
            Traceback (most recent call last):
            ...
            ValueError: PPL::C_Polyhedron::add_generator(g):
            this->space_dimension() == 1, g.space_dimension() == 2.

         The constraint must also be topology-compatible, that is,
         :class:`C_Polyhedron` does not allow :func:`closure_point`
         generators::

            sage: p = C_Polyhedron( point(0*x+0*y) )
            sage: p.add_generator( closure_point(0*x) )
            Traceback (most recent call last):
            ...
            ValueError: PPL::C_Polyhedron::add_generator(g):
            g is a closure point.

        Finally, ever non-empty polyhedron must have at least one
        point generator::

            sage: p = C_Polyhedron(3, 'empty')
            sage: p.add_generator( ray(x) )
            Traceback (most recent call last):
            ...
            ValueError: PPL::C_Polyhedron::add_generator(g):
            *this is an empty polyhedron and g is not a point.
        """
        self.assert_mutable('The Polyhedron is not mutable!')
        sig_on()
        try:
            self.thisptr.add_generator(g.thisptr[0])
        finally:
            sig_off()


    def add_constraints(self, Constraint_System cs):
        r"""
        Add constraints to the polyhedron.

        Adds a copy of constraints in ``cs`` to the system of constraints
        of ``self``, without minimizing the result.

        See alse :meth:`add_constraint`.

        INPUT:

        - ``cs`` -- the :class:`Constraint_System` that will be added
          to the system of constraints of ``self``.

        OUTPUT:

        This method modifies the polyhedron ``self`` and does not
        return anything.

        Raises a ``ValueError`` if ``self`` and the constraints in
        ``cs`` are topology-incompatible or dimension-incompatible.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, C_Polyhedron, Constraint_System
            sage: x = Variable(0)
            sage: y = Variable(1)
            sage: cs = Constraint_System()
            sage: cs.insert(x>=0)
            sage: cs.insert(y>=0)
            sage: p = C_Polyhedron( y<=1 )
            sage: p.add_constraints(cs)

         We just added a 1-d constraint to a 2-d polyhedron, this is
         fine. The other way is not::

            sage: p = C_Polyhedron( x<=1 )
            sage: p.add_constraints(cs)
            Traceback (most recent call last):
            ...
            ValueError: PPL::C_Polyhedron::add_recycled_constraints(cs):
            this->space_dimension() == 1, cs.space_dimension() == 2.

         The constraints must also be topology-compatible, that is,
         :class:`C_Polyhedron` only allows non-strict inequalities::

            sage: p = C_Polyhedron( x>=0 )
            sage: p.add_constraints( Constraint_System(x<0) )
            Traceback (most recent call last):
            ...
            ValueError: PPL::C_Polyhedron::add_recycled_constraints(cs):
            cs contains strict inequalities.
        """
        self.assert_mutable('The Polyhedron is not mutable!')
        sig_on()
        try:
            self.thisptr.add_constraints(cs.thisptr[0])
        finally:
            sig_off()


    def add_generators(self, Generator_System gs):
        r"""
        Add generators to the polyhedron.

        Adds a copy of the generators in ``gs`` to the system of
        generators of ``self``, without minimizing the result.

        See alse :meth:`add_generator`.

        INPUT:

        - ``gs`` -- the :class:`Generator_System` that will be added
          to the system of constraints of ``self``.

        OUTPUT:

        This method modifies the polyhedron ``self`` and does not
        return anything.

        Raises a ``ValueError`` if ``self`` and one of the generators
        in ``gs`` are topology-incompatible or dimension-incompatible,
        or if ``self`` is an empty polyhedron and ``gs`` does not
        contain a point.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, C_Polyhedron, Generator_System, point, ray, closure_point
            sage: x = Variable(0)
            sage: y = Variable(1)
            sage: gs = Generator_System()
            sage: gs.insert(point(0*x+0*y))
            sage: gs.insert(point(1*x+1*y))
            sage: p = C_Polyhedron(2, 'empty')
            sage: p.add_generators(gs)

         We just added a 1-d constraint to a 2-d polyhedron, this is
         fine. The other way is not::

            sage: p = C_Polyhedron(1, 'empty')
            sage: p.add_generators(gs)
            Traceback (most recent call last):
            ...
            ValueError: PPL::C_Polyhedron::add_recycled_generators(gs):
            this->space_dimension() == 1, gs.space_dimension() == 2.

         The constraints must also be topology-compatible, that is,
         :class:`C_Polyhedron` does not allow :func:`closure_point`
         generators::

            sage: p = C_Polyhedron( point(0*x+0*y) )
            sage: p.add_generators( Generator_System(closure_point(x) ))
            Traceback (most recent call last):
            ...
            ValueError: PPL::C_Polyhedron::add_recycled_generators(gs):
            gs contains closure points.
        """
        self.assert_mutable('The Polyhedron is not mutable!')
        sig_on()
        try:
            self.thisptr.add_generators(gs.thisptr[0])
        finally:
            sig_off()


    def unconstrain(self, Variable var):
        r"""
        Compute the cylindrification of ``self`` with respect to space
        dimension ``var``.

        INPUT:

        - ``var`` -- a :class:`Variable`. The space dimension that
          will be unconstrained.  Exceptions:

        OUTPUT:

        This method assigns the cylindrification to ``self`` and does
        not return anything.

        Raises a ``ValueError`` if ``var`` is not a space dimension of
        ``self``.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, C_Polyhedron, point
            sage: x = Variable(0)
            sage: y = Variable(1)
            sage: p = C_Polyhedron( point(x+y) ); p
            A 0-dimensional polyhedron in QQ^2 defined as the convex hull of 1 point
            sage: p.unconstrain(x); p
            A 1-dimensional polyhedron in QQ^2 defined as the convex hull of 1 point, 1 line
            sage: z = Variable(2)
            sage: p.unconstrain(z)
            Traceback (most recent call last):
            ...
            ValueError: PPL::C_Polyhedron::unconstrain(var):
            this->space_dimension() == 2, required space dimension == 3.
        """
        sig_on()
        try:
            self.thisptr.unconstrain(var.thisptr[0])
        finally:
            sig_off()


    def intersection_assign(self, Polyhedron y):
        r"""
        Assign to ``self`` the intersection of ``self`` and ``y``.

        INPUT:

        - ``y`` -- a :class:`Polyhedron`

        OUTPUT:

        This method assigns the intersection to ``self`` and does not
        return anything.

        Raises a ``ValueError`` if ``self`` and and ``y`` are
        topology-incompatible or dimension-incompatible.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, C_Polyhedron, NNC_Polyhedron
            sage: x = Variable(0)
            sage: y = Variable(1)
            sage: p = C_Polyhedron( 1*x+0*y >= 0 )
            sage: p.intersection_assign( C_Polyhedron(y>=0) )
            sage: p.constraints()
            Constraint_System {x0>=0, x1>=0}
            sage: z = Variable(2)
            sage: p.intersection_assign( C_Polyhedron(z>=0) )
            Traceback (most recent call last):
            ...
            ValueError: PPL::C_Polyhedron::intersection_assign(y):
            this->space_dimension() == 2, y.space_dimension() == 3.
            sage: p.intersection_assign( NNC_Polyhedron(x+y<1) )
            Traceback (most recent call last):
            ...
            ValueError: PPL::C_Polyhedron::intersection_assign(y):
            y is a NNC_Polyhedron.
        """
        self.assert_mutable('The Polyhedron is not mutable!')
        sig_on()
        try:
            self.thisptr.intersection_assign(y.thisptr[0])
        finally:
            sig_off()


    def poly_hull_assign(self, Polyhedron y):
        r"""
        Assign to ``self`` the poly-hull of ``self`` and ``y``.

        For any pair of NNC polyhedra `P_1` and `P_2`, the convex
        polyhedral hull (or poly-hull) of is the smallest NNC
        polyhedron that includes both `P_1` and `P_2`. The poly-hull
        of any pair of closed polyhedra in is also closed.

        INPUT:

        - ``y`` -- a :class:`Polyhedron`

        OUTPUT:

        This method assigns the poly-hull to ``self`` and does not
        return anything.

        Raises a ``ValueError`` if ``self`` and and ``y`` are
        topology-incompatible or dimension-incompatible.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, C_Polyhedron, point, NNC_Polyhedron
            sage: x = Variable(0)
            sage: y = Variable(1)
            sage: p = C_Polyhedron( point(1*x+0*y) )
            sage: p.poly_hull_assign(C_Polyhedron( point(0*x+1*y) ))
            sage: p.generators()
            Generator_System {point(0/1, 1/1), point(1/1, 0/1)}

        ``self`` and ``y`` must be dimension- and topology-compatible,
        or an exception is raised::

            sage: z = Variable(2)
            sage: p.poly_hull_assign( C_Polyhedron(z>=0) )
            Traceback (most recent call last):
            ...
            ValueError: PPL::C_Polyhedron::poly_hull_assign(y):
            this->space_dimension() == 2, y.space_dimension() == 3.
            sage: p.poly_hull_assign( NNC_Polyhedron(x+y<1) )
            Traceback (most recent call last):
            ...
            ValueError: PPL::C_Polyhedron::poly_hull_assign(y):
            y is a NNC_Polyhedron.
        """
        self.assert_mutable('The Polyhedron is not mutable!')
        sig_on()
        try:
            self.thisptr.poly_hull_assign(y.thisptr[0])
        finally:
            sig_off()


    upper_bound_assign = poly_hull_assign


    def poly_difference_assign(self, Polyhedron y):
        r"""
        Assign to ``self`` the poly-difference of ``self`` and ``y``.

        For any pair of NNC polyhedra `P_1` and `P_2` the convex
        polyhedral difference (or poly-difference) of `P_1` and `P_2`
        is defined as the smallest convex polyhedron containing the
        set-theoretic difference `P_1\setminus P_2` of `P_1` and
        `P_2`.

        In general, even if `P_1` and `P_2` are topologically closed
        polyhedra, their poly-difference may be a convex polyhedron
        that is not topologically closed. For this reason, when
        computing the poly-difference of two :class:`C_Polyhedron`,
        the library will enforce the topological closure of the
        result.

        INPUT:

        - ``y`` -- a :class:`Polyhedron`

        OUTPUT:

        This method assigns the poly-difference to ``self`` and does
        not return anything.

        Raises a ``ValueError`` if ``self`` and and ``y`` are
        topology-incompatible or dimension-incompatible.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, C_Polyhedron, point, closure_point, NNC_Polyhedron
            sage: x = Variable(0)
            sage: p = NNC_Polyhedron( point(0*x) )
            sage: p.add_generator( point(1*x) )
            sage: p.poly_difference_assign(NNC_Polyhedron( point(0*x) ))
            sage: p.minimized_constraints()
            Constraint_System {-x0+1>=0, x0>0}

        The poly-difference of :class:`C_polyhedron` is really its closure::

            sage: p = C_Polyhedron( point(0*x) )
            sage: p.add_generator( point(1*x) )
            sage: p.poly_difference_assign(C_Polyhedron( point(0*x) ))
            sage: p.minimized_constraints()
            Constraint_System {x0>=0, -x0+1>=0}

        ``self`` and ``y`` must be dimension- and topology-compatible,
        or an exception is raised::

            sage: y = Variable(1)
            sage: p.poly_difference_assign( C_Polyhedron(y>=0) )
            Traceback (most recent call last):
            ...
            ValueError: PPL::C_Polyhedron::poly_difference_assign(y):
            this->space_dimension() == 1, y.space_dimension() == 2.
            sage: p.poly_difference_assign( NNC_Polyhedron(x+y<1) )
            Traceback (most recent call last):
            ...
            ValueError: PPL::C_Polyhedron::poly_difference_assign(y):
            y is a NNC_Polyhedron.
        """
        self.assert_mutable('The Polyhedron is not mutable!')
        sig_on()
        try:
            self.thisptr.poly_difference_assign(y.thisptr[0])
        finally:
            sig_off()


    difference_assign = poly_difference_assign


    def drop_some_non_integer_points(self):
        r"""
        Possibly tighten ``self`` by dropping some points with
        non-integer coordinates.

        The modified polyhedron satisfies:

        * it is (not necessarily strictly) contained in the original
          polyhedron.

        * integral vertices (generating points with integer
          coordinates) of the original polyhedron are not removed.

        .. NOTE::

            The modified polyhedron is not neccessarily a lattice
            polyhedron; Some vertices will, in general, still be
            rational. Lattice points interior to the polyhedron may be
            lost in the process.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, NNC_Polyhedron, Constraint_System
            sage: x = Variable(0)
            sage: y = Variable(1)
            sage: cs = Constraint_System()
            sage: cs.insert( x>=0 )
            sage: cs.insert( y>=0 )
            sage: cs.insert( 3*x+2*y<5 )
            sage: p = NNC_Polyhedron(cs)
            sage: p.minimized_generators()
            Generator_System {point(0/1, 0/1), closure_point(0/2, 5/2), closure_point(5/3, 0/3)}
            sage: p.drop_some_non_integer_points()
            sage: p.minimized_generators()
            Generator_System {point(0/1, 0/1), point(0/1, 2/1), point(4/3, 0/3)}
        """
        self.assert_mutable('The Polyhedron is not mutable!')
        sig_on()
        self.thisptr.drop_some_non_integer_points()
        sig_off()


    def topological_closure_assign(self):
        r"""
        Assign to ``self`` its topological closure.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, NNC_Polyhedron
            sage: x = Variable(0)
            sage: p = NNC_Polyhedron(x>0)
            sage: p.is_topologically_closed()
            False
            sage: p.topological_closure_assign()
            sage: p.is_topologically_closed()
            True
            sage: p.minimized_constraints()
            Constraint_System {x0>=0}
        """
        self.assert_mutable('The Polyhedron is not mutable!')
        sig_on()
        self.thisptr.topological_closure_assign()
        sig_off()


    def add_space_dimensions_and_embed(self, m):
        r"""
        Add ``m`` new space dimensions and embed ``self`` in the new
        vector space.

        The new space dimensions will be those having the highest
        indexes in the new polyhedron, which is characterized by a
        system of constraints in which the variables running through
        the new dimensions are not constrained. For instance, when
        starting from the polyhedron `P` and adding a third space
        dimension, the result will be the polyhedron

        .. MATH::

            \Big\{
            (x,y,z)^T \in \RR^3
            \Big|
            (x,y)^T \in P
            \Big\}

        INPUT:

        - ``m`` -- integer.

        OUTPUT:

        This method assigns the embedded polyhedron to ``self`` and
        does not return anything.

        Raises a ``ValueError`` if adding ``m`` new space dimensions
        would cause the vector space to exceed dimension
        ``self.max_space_dimension()``.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, C_Polyhedron, point
            sage: x = Variable(0)
            sage: p = C_Polyhedron( point(3*x) )
            sage: p.add_space_dimensions_and_embed(1)
            sage: p.minimized_generators()
            Generator_System {line(0, 1), point(3/1, 0/1)}
            sage: p.add_space_dimensions_and_embed( p.max_space_dimension() )
            Traceback (most recent call last):
            ...
            ValueError: PPL::C_Polyhedron::add_space_dimensions_and_embed(m):
            adding m new space dimensions exceeds the maximum allowed space dimension.
        """
        self.assert_mutable('The Polyhedron is not mutable!')
        m = int(m)
        sig_on()
        try:
            self.thisptr.add_space_dimensions_and_embed(m)
        finally:
            sig_off()


    def add_space_dimensions_and_project(self, m):
        r"""
        Add ``m`` new space dimensions and embed ``self`` in the new
        vector space.

        The new space dimensions will be those having the highest
        indexes in the new polyhedron, which is characterized by a
        system of constraints in which the variables running through
        the new dimensions are all constrained to be equal to `0`.
        For instance, when starting from the polyhedron `P` and adding
        a third space dimension, the result will be the polyhedron

        .. MATH::

            \Big\{
            (x,y,0)^T \in \RR^3
            \Big|
            (x,y)^T \in P
            \Big\}

        INPUT:

        - ``m`` -- integer.

        OUTPUT:

        This method assigns the projected polyhedron to ``self`` and
        does not return anything.

        Raises a ``ValueError`` if adding ``m`` new space dimensions
        would cause the vector space to exceed dimension
        ``self.max_space_dimension()``.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, C_Polyhedron, point
            sage: x = Variable(0)
            sage: p = C_Polyhedron( point(3*x) )
            sage: p.add_space_dimensions_and_project(1)
            sage: p.minimized_generators()
            Generator_System {point(3/1, 0/1)}
            sage: p.add_space_dimensions_and_project( p.max_space_dimension() )
            Traceback (most recent call last):
            ...
            ValueError: PPL::C_Polyhedron::add_space_dimensions_and_project(m):
            adding m new space dimensions exceeds the maximum allowed space dimension.
        """
        self.assert_mutable('The Polyhedron is not mutable!')
        m = int(m)
        sig_on()
        try:
            self.thisptr.add_space_dimensions_and_project(m)
        finally:
            sig_off()


    def concatenate_assign(self, Polyhedron y):
        r"""
        Assign to ``self`` the concatenation of ``self`` and ``y``.

        This functions returns the Cartiesian product of ``self`` and
        ``y``.

        Viewing a polyhedron as a set of tuples (its points), it is
        sometimes useful to consider the set of tuples obtained by
        concatenating an ordered pair of polyhedra. Formally, the
        concatenation of the polyhedra `P` and `Q` (taken in this
        order) is the polyhedron such that

        .. MATH::

            R =
            \Big\{
            (x_0,\dots,x_{n-1},y_0,\dots,y_{m-1})^T \in \RR^{n+m}
            \Big|
            (x_0,\dots,x_{n-1})^T \in P
            ,~
            (y_0,\dots,y_{m-1})^T \in Q
            \Big\}

        Another way of seeing it is as follows: first embed polyhedron
        `P` into a vector space of dimension `n+m` and then add a
        suitably renamed-apart version of the constraints defining
        `Q`.

        INPUT:

        - ``m`` -- integer.

        OUTPUT:

        This method assigns the concatenated polyhedron to ``self`` and
        does not return anything.

        Raises a ``ValueError`` if ``self`` and ``y`` are
        topology-incompatible or if adding ``y.space_dimension()`` new
        space dimensions would cause the vector space to exceed
        dimension ``self.max_space_dimension()``.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, C_Polyhedron, NNC_Polyhedron, point
            sage: x = Variable(0)
            sage: p1 = C_Polyhedron( point(1*x) )
            sage: p2 = C_Polyhedron( point(2*x) )
            sage: p1.concatenate_assign(p2)
            sage: p1.minimized_generators()
            Generator_System {point(1/1, 2/1)}

        The polyhedra must be topology-compatible and not exceed the
        maximum space dimension::

            sage: p1.concatenate_assign( NNC_Polyhedron(1, 'universe') )
            Traceback (most recent call last):
            ...
            ValueError: PPL::C_Polyhedron::concatenate_assign(y):
            y is a NNC_Polyhedron.
            sage: p1.concatenate_assign( C_Polyhedron(p1.max_space_dimension(), 'empty') )
            Traceback (most recent call last):
            ...
            ValueError: PPL::C_Polyhedron::concatenate_assign(y):
            concatenation exceeds the maximum allowed space dimension.
        """
        self.assert_mutable('The Polyhedron is not mutable!')
        sig_on()
        try:
            self.thisptr.concatenate_assign(y.thisptr[0])
        finally:
            sig_off()


    def remove_higher_space_dimensions(self, new_dimension):
        r"""
        Remove the higher dimensions of the vector space so that the
        resulting space will have dimension ``new_dimension``.

        OUTPUT:

        This method modifies ``self`` and does not return anything.

        Raises a ``ValueError`` if ``new_dimensions`` is greater than
        the space dimension of ``self``.

        EXAMPLES::

            sage: from sage.libs.ppl import C_Polyhedron, Variable
            sage: x = Variable(0)
            sage: y = Variable(1)
            sage: p = C_Polyhedron(3*x+0*y==2)
            sage: p.remove_higher_space_dimensions(1)
            sage: p.minimized_constraints()
            Constraint_System {3*x0-2==0}
            sage: p.remove_higher_space_dimensions(2)
            Traceback (most recent call last):
            ...
            ValueError: PPL::C_Polyhedron::remove_higher_space_dimensions(nd):
            this->space_dimension() == 1, required space dimension == 2.
        """
        self.assert_mutable('The Polyhedron is not mutable!')
        new_dimension = int(new_dimension)
        sig_on()
        try:
            self.thisptr.remove_higher_space_dimensions(new_dimension)
        finally:
            sig_off()


    def ascii_dump(self):
        r"""
        Write an ASCII dump to stderr.

        EXAMPLES::

            sage: sage_cmd  = 'from sage.libs.ppl import C_Polyhedron, Variable\n'
            sage: sage_cmd += 'x = Variable(0)\n'
            sage: sage_cmd += 'y = Variable(1)\n'
            sage: sage_cmd += 'p = C_Polyhedron(3*x+2*y==1)\n'
            sage: sage_cmd += 'p.minimized_generators()\n'
            sage: sage_cmd += 'p.ascii_dump()\n'
            sage: from sage.tests.cmdline import test_executable
            sage: (out, err, ret) = test_executable(['sage', '-c', sage_cmd], timeout=100);  # long time, indirect doctest
            sage: print err  # long time
            space_dim 2
            -ZE -EM  +CM +GM  +CS +GS  -CP -GP  -SC +SG
            con_sys (up-to-date)
            topology NECESSARILY_CLOSED
            2 x 2 SPARSE (sorted)
            index_first_pending 2
            size 3 -1 3 2 = (C)
            size 3 1 0 0 >= (C)
            <BLANKLINE>
            gen_sys (up-to-date)
            topology NECESSARILY_CLOSED
            2 x 2 DENSE (not_sorted)
            index_first_pending 2
            size 3 0 2 -3 L (C)
            size 3 2 0 1 P (C)
            <BLANKLINE>
            sat_c
            0 x 0
            <BLANKLINE>
            sat_g
            2 x 2
            0 0
            0 1
        """
        sig_on()
        self.thisptr.ascii_dump()
        sig_off()


    def max_space_dimension(self):
        r"""
        Return the maximum space dimension all kinds of Polyhedron can handle.

        OUTPUT:

        Integer.

        EXAMPLES::

            sage: from sage.libs.ppl import C_Polyhedron
            sage: C_Polyhedron(1, 'empty').max_space_dimension()   # random output
            1152921504606846974
            sage: C_Polyhedron(1, 'empty').max_space_dimension()
            357913940            # 32-bit
            1152921504606846974  # 64-bit
        """
        return self.thisptr.max_space_dimension()


    def OK(self, check_non_empty=False):
        """
        Check if all the invariants are satisfied.

        The check is performed so as to intrude as little as
        possible. If the library has been compiled with run-time
        assertions enabled, error messages are written on std::cerr in
        case invariants are violated. This is useful for the purpose
        of debugging the library.

        INPUT:

        - ``check_not_empty`` -- boolean. ``True`` if and only if, in
          addition to checking the invariants, ``self`` must be
          checked to be not empty.

        OUTPUT:

        ``True`` if and only if ``self`` satisfies all the invariants
        and either ``check_not_empty`` is ``False`` or ``self`` is not
        empty.

        EXAMPLES::

            sage: from sage.libs.ppl import Linear_Expression, Variable
            sage: x = Variable(0)
            sage: y = Variable(1)
            sage: e = 3*x+2*y+1
            sage: e.OK()
            True
        """
        sig_on()
        cdef bint result = self.thisptr.OK()
        sig_off()
        return result


    def __richcmp__(Polyhedron lhs, Polyhedron rhs, op):
        r"""
        Comparison for polyhedra.

        INPUT:

        - ``lhs``, ``rhs`` -- :class:`Polyhedron`.

        - ``op`` -- integer. The comparison operation to be performed.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, C_Polyhedron
            sage: x = Variable(0)
            sage: C_Polyhedron(x>=0) > C_Polyhedron(x>=1)    # indirect doctest
            True
        """
        cdef result
        sig_on()
        if op==0:      # <   0
            result = rhs.strictly_contains(lhs)
        elif op==1:    # <=  1
            result = rhs.contains(lhs)
        elif op==2:    # ==  2
            result = (lhs.thisptr[0] == rhs.thisptr[0])
        elif op==4:    # >   4
            result = lhs.strictly_contains(rhs)
        elif op==5:    # >=  5
            result = lhs.contains(rhs)
        elif op==3:    # !=  3
            result = (lhs.thisptr[0] != rhs.thisptr[0])
        else:
            assert False  # unreachable
        sig_off()
        return result



####################################################
### C_Polyhedron ###################################
####################################################
cdef class C_Polyhedron(Polyhedron):
    r"""
    Wrapper for PPL's ``C_Polyhedron`` class.

    An object of the class :class:`C_Polyhedron` represents a
    topologically closed convex polyhedron in the vector space. See
    :class:`NNC_Polyhedron` for more general (not necessarily closed)
    polyhedra.

    When building a closed polyhedron starting from a system of
    constraints, an exception is thrown if the system contains a
    strict inequality constraint. Similarly, an exception is thrown
    when building a closed polyhedron starting from a system of
    generators containing a closure point.

    INPUT:

    - ``arg`` -- the defining data of the polyhedron. Any one of the
      following is accepted:

      * A non-negative integer. Depending on ``degenerate_element``,
        either the space-filling or the empty polytope in the given
        dimension ``arg`` is constructed.

      * A :class:`Constraint_System`.

      * A :class:`Generator_System`.

      * A single :class:`Constraint`.

      * A single :class:`Generator`.

      * A :class:`C_Polyhedron`.

    - ``degenerate_element`` -- string, either ``'universe'`` or
      ``'empty'``. Only used if ``arg`` is an integer.

    OUTPUT:

    A :class:`C_Polyhedron`.

    EXAMPLES::

        sage: from sage.libs.ppl import Constraint, Constraint_System, Generator, Generator_System, Variable, C_Polyhedron, point, ray
        sage: x = Variable(0)
        sage: y = Variable(1)
        sage: C_Polyhedron( 5*x-2*y >=  x+y-1 )
        A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 point, 1 ray, 1 line
        sage: cs = Constraint_System()
        sage: cs.insert( x >= 0 )
        sage: cs.insert( y >= 0 )
        sage: C_Polyhedron(cs)
        A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 point, 2 rays
        sage: C_Polyhedron( point(x+y) )
        A 0-dimensional polyhedron in QQ^2 defined as the convex hull of 1 point
        sage: gs = Generator_System()
        sage: gs.insert( point(-x-y) )
        sage: gs.insert( ray(x) )
        sage: C_Polyhedron(gs)
        A 1-dimensional polyhedron in QQ^2 defined as the convex hull of 1 point, 1 ray

    The empty and universe polyhedra are constructed like this::

        sage: C_Polyhedron(3, 'empty')
        The empty polyhedron in QQ^3
        sage: C_Polyhedron(3, 'empty').constraints()
        Constraint_System {-1==0}
        sage: C_Polyhedron(3, 'universe')
        The space-filling polyhedron in QQ^3
        sage: C_Polyhedron(3, 'universe').constraints()
        Constraint_System {}

    Note that, by convention, the generator system of a polyhedron is
    either empty or contains at least one point. In particular, if you
    define a polyhedron via a non-empty :class:`Generator_System` it
    must contain a point (at any position). If you start with a single
    generator, this generator must be a point::

        sage: C_Polyhedron( ray(x) )
        Traceback (most recent call last):
        ...
        ValueError: PPL::C_Polyhedron::C_Polyhedron(gs):
        *this is an empty polyhedron and
        the non-empty generator system gs contains no points.
    """


    def __cinit__(self, arg, degenerate_element='universe'):
        """
        The Cython constructor.

        See :class:`C_Polyhedron` for documentation.

        TESTS::

            sage: from sage.libs.ppl import C_Polyhedron
            sage: C_Polyhedron(3, 'empty')   # indirect doctest
            The empty polyhedron in QQ^3
        """
        if isinstance(arg, C_Polyhedron):
            ph = <C_Polyhedron>arg
            self.thisptr = new PPL_C_Polyhedron(<PPL_C_Polyhedron&>ph.thisptr[0])
            return
        if isinstance(arg, Generator):
            arg = Generator_System(arg)
        if isinstance(arg, Constraint):
            arg = Constraint_System(arg)
        if isinstance(arg, Generator_System):
            gs = <Generator_System>arg
            self.thisptr = new PPL_C_Polyhedron(gs.thisptr[0])
            return
        if isinstance(arg, Constraint_System):
            cs = <Constraint_System>arg
            self.thisptr = new PPL_C_Polyhedron(cs.thisptr[0])
            return
        try:
            dim = int(arg)
            assert dim>=0
        except ValueError:
            raise ValueError, 'Cannot initialize C_Polyhedron with '+str(arg)+'.'
        degenerate_element = degenerate_element.lower()
        if degenerate_element=='universe':
            self.thisptr = new PPL_C_Polyhedron(<PPL_dimension_type>dim, UNIVERSE)
            return
        elif degenerate_element=='empty':
            self.thisptr = new PPL_C_Polyhedron(<PPL_dimension_type>dim, EMPTY)
            return
        else:
            raise ValueError, 'Unknown value: degenerate_element='+str(degenerate_element)+'.'


    def __init__(self, *args):
        """
        The Python destructor.

        See :class:`C_Polyhedron` for documentation.

        TESTS::

            sage: from sage.libs.ppl import C_Polyhedron
            sage: C_Polyhedron(3, 'empty')   # indirect doctest
            The empty polyhedron in QQ^3
        """
        # override Polyhedron.__init__
        pass


    def __dealloc__(self):
        """
        The Cython destructor.
        """
        del self.thisptr


    def __reduce__(self):
        """
        Pickle object

        TESTS::

            sage: from sage.libs.ppl import C_Polyhedron, Variable
            sage: P = C_Polyhedron(3, 'empty')
            sage: loads(dumps(P))
            The empty polyhedron in QQ^3

            sage: Q = C_Polyhedron(5, 'universe')
            sage: loads(dumps(Q))
            The space-filling polyhedron in QQ^5

            sage: x = Variable(0)
            sage: y = Variable(1)
            sage: H = C_Polyhedron( 5*x-2*y >=  x+y-1 )
            sage: loads(dumps(H))
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 point, 1 ray, 1 line
        """
        if self.is_empty():
            return (C_Polyhedron, (self.space_dimension(), 'empty'))
        elif self.is_universe():
            return (C_Polyhedron, (self.space_dimension(), 'universe'))
        else:
            return (C_Polyhedron, (self.generators(),))


####################################################
### NNC_Polyhedron #################################
####################################################
cdef class NNC_Polyhedron(Polyhedron):
    r"""
    Wrapper for PPL's ``NNC_Polyhedron`` class.

    An object of the class ``NNC_Polyhedron`` represents a not
    necessarily closed (NNC) convex polyhedron in the vector space.

    Note: Since NNC polyhedra are a generalization of closed
    polyhedra, any object of the class :class:`C_Polyhedron` can be
    (explicitly) converted into an object of the class
    :class:`NNC_Polyhedron`. The reason for defining two different
    classes is that objects of the class :class:`C_Polyhedron` are
    characterized by a more efficient implementation, requiring less
    time and memory resources.

    INPUT:

    - ``arg`` -- the defining data of the polyhedron. Any one of the
      following is accepted:

      * An non-negative integer. Depending on ``degenerate_element``,
        either the space-filling or the empty polytope in the given
        dimension ``arg`` is constructed.

      * A :class:`Constraint_System`.

      * A :class:`Generator_System`.

      * A single :class:`Constraint`.

      * A single :class:`Generator`.

      * A :class:`NNC_Polyhedron`.

      * A :class:`C_Polyhedron`.

    - ``degenerate_element`` -- string, either ``'universe'`` or
      ``'empty'``. Only used if ``arg`` is an integer.

    OUTPUT:

    A :class:`C_Polyhedron`.

    EXAMPLES::

        sage: from sage.libs.ppl import Constraint, Constraint_System, Generator, Generator_System, Variable, NNC_Polyhedron, point, ray, closure_point
        sage: x = Variable(0)
        sage: y = Variable(1)
        sage: NNC_Polyhedron( 5*x-2*y >  x+y-1 )
        A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 point, 1 closure_point, 1 ray, 1 line
        sage: cs = Constraint_System()
        sage: cs.insert( x > 0 )
        sage: cs.insert( y > 0 )
        sage: NNC_Polyhedron(cs)
        A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 point, 1 closure_point, 2 rays
        sage: NNC_Polyhedron( point(x+y) )
        A 0-dimensional polyhedron in QQ^2 defined as the convex hull of 1 point
        sage: gs = Generator_System()
        sage: gs.insert( point(-y) )
        sage: gs.insert( closure_point(-x-y) )
        sage: gs.insert( ray(x) )
        sage: p = NNC_Polyhedron(gs); p
        A 1-dimensional polyhedron in QQ^2 defined as the convex hull of 1 point, 1 closure_point, 1 ray
        sage: p.minimized_constraints()
        Constraint_System {x1+1==0, x0+1>0}

    Note that, by convention, every polyhedron must contain a point::

        sage: NNC_Polyhedron( closure_point(x+y) )
        Traceback (most recent call last):
        ...
        ValueError: PPL::NNC_Polyhedron::NNC_Polyhedron(gs):
        *this is an empty polyhedron and
        the non-empty generator system gs contains no points.
    """


    def __cinit__(self, arg, degenerate_element='universe'):
        """
        The Cython constructor.

        See :class:`NNC_Polyhedron` for documentation.

        TESTS::

            sage: from sage.libs.ppl import NNC_Polyhedron
            sage: NNC_Polyhedron(3, 'empty')   # indirect doctest
            The empty polyhedron in QQ^3
        """
        if isinstance(arg, NNC_Polyhedron):
            p_nnc = <NNC_Polyhedron>arg
            self.thisptr = new PPL_NNC_Polyhedron(<PPL_NNC_Polyhedron&>p_nnc.thisptr[0])
            return
        if isinstance(arg, C_Polyhedron):
            p_c = <C_Polyhedron>arg
            self.thisptr = new PPL_NNC_Polyhedron(<PPL_C_Polyhedron&>p_c.thisptr[0])
            return
        if isinstance(arg, Generator):
            arg = Generator_System(arg)
        if isinstance(arg, Constraint):
            arg = Constraint_System(arg)
        if isinstance(arg, Generator_System):
            gs = <Generator_System>arg
            self.thisptr = new PPL_NNC_Polyhedron(gs.thisptr[0])
            return
        if isinstance(arg, Constraint_System):
            cs = <Constraint_System>arg
            self.thisptr = new PPL_NNC_Polyhedron(cs.thisptr[0])
            return
        try:
            dim = int(arg)
            assert dim>=0
        except ValueError:
            raise ValueError, 'Cannot initialize NNC_Polyhedron with '+str(arg)+'.'
        degenerate_element = degenerate_element.lower()
        if degenerate_element=='universe':
            self.thisptr = new PPL_NNC_Polyhedron(<PPL_dimension_type>dim, UNIVERSE)
            return
        elif degenerate_element=='empty':
            self.thisptr = new PPL_NNC_Polyhedron(<PPL_dimension_type>dim, EMPTY)
            return
        else:
            raise ValueError, 'Unknown value: degenerate_element='+str(degenerate_element)+'.'


    def __init__(self, *args):
        """
        The Python destructor.

        See :class:`NNC_Polyhedron` for documentation.

        TESTS::

            sage: from sage.libs.ppl import NNC_Polyhedron
            sage: NNC_Polyhedron(3, 'empty')   # indirect doctest
            The empty polyhedron in QQ^3
        """
        # override Polyhedron.__init__
        pass


    def __dealloc__(self):
        """
        The Cython destructor.
        """
        del self.thisptr


    def __reduce__(self):
        """
        Pickle object

        TESTS::

            sage: from sage.libs.ppl import NNC_Polyhedron, Variable
            sage: P = NNC_Polyhedron(3, 'empty')
            sage: loads(dumps(P))
            The empty polyhedron in QQ^3

            sage: Q = NNC_Polyhedron(5, 'universe')
            sage: loads(dumps(Q))
            The space-filling polyhedron in QQ^5

            sage: x = Variable(0)
            sage: y = Variable(1)
            sage: H = NNC_Polyhedron( 5*x-2*y >  x+y-1 )
            sage: loads(dumps(H))
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 point,
            1 closure_point, 1 ray, 1 line
        """
        if self.is_empty():
            return (NNC_Polyhedron, (self.space_dimension(), 'empty'))
        elif self.is_universe():
            return (NNC_Polyhedron, (self.space_dimension(), 'universe'))
        else:
            return (NNC_Polyhedron, (self.generators(),))


####################################################
### Variable #######################################
####################################################
cdef class Variable(object):
    r"""
    Wrapper for PPL's ``Variable`` class.

    A dimension of the vector space.

    An object of the class Variable represents a dimension of the
    space, that is one of the Cartesian axes. Variables are used as
    basic blocks in order to build more complex linear
    expressions. Each variable is identified by a non-negative
    integer, representing the index of the corresponding Cartesian
    axis (the first axis has index 0). The space dimension of a
    variable is the dimension of the vector space made by all the
    Cartesian axes having an index less than or equal to that of the
    considered variable; thus, if a variable has index `i`, its space
    dimension is `i+1`.

    INPUT:

    - ``i`` -- integer. The index of the axis.

    OUTPUT:

    A :class:`Variable`

    EXAMPLES::

        sage: from sage.libs.ppl import Variable
        sage: x = Variable(123)
        sage: x.id()
        123
        sage: x
        x123

    Note that the "meaning" of an object of the class Variable is
    completely specified by the integer index provided to its
    constructor: be careful not to be mislead by C++ language variable
    names. For instance, in the following example the linear
    expressions ``e1`` and ``e2`` are equivalent, since the two
    variables ``x`` and ``z`` denote the same Cartesian axis::

        sage: x = Variable(0)
        sage: y = Variable(1)
        sage: z = Variable(0)
        sage: e1 = x + y; e1
        x0+x1
        sage: e2 = y + z; e2
        x0+x1
        sage: e1 - e2
        0
    """

    cdef PPL_Variable *thisptr


    def __cinit__(self, PPL_dimension_type i):
        """
        The Cython constructor.

        See :class:`Variable` for documentation.

        TESTS::

            sage: from sage.libs.ppl import Variable
            sage: Variable(123)   # indirect doctest
            x123
        """
        self.thisptr = new PPL_Variable(i)


    def __dealloc__(self):
        """
        The Cython destructor.
        """
        del self.thisptr


    def id(self):
        """
        Return the index of the Cartesian axis associated to the variable.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable
            sage: x = Variable(123)
            sage: x.id()
            123
        """
        return self.thisptr.id()


    def OK(self):
        """
        Checks if all the invariants are satisfied.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable
            sage: x = Variable(0)
            sage: x.OK()
            True
        """
        return self.thisptr.OK()


    def space_dimension(self):
        r"""
        Return the dimension of the vector space enclosing ``self``.

        OUPUT:

        Integer. The returned value is ``self.id()+1``.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable
            sage: x = Variable(0)
            sage: x.space_dimension()
            1
        """
        return self.thisptr.space_dimension()


    def __repr__(self):
        """
        Return a string representation.

        OUTPUT:

        String.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable
            sage: x = Variable(0)
            sage: x.__repr__()
            'x0'
        """
        return 'x{0}'.format(self.id())


    def __add__(self, other):
        r"""
        Return the sum ``self`` + ``other``.

        INPUT:

        - ``self``, ``other`` -- anything convertible to
          ``Linear_Expression``: An integer, a :class:`Variable`, or a
          :class:`Linear_Expression`.

        OUTPUT:

        A :class:`Linear_Expression` representing ``self`` + ``other``.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable
            sage: x = Variable(0); y = Variable(1)
            sage: x + 15
            x0+15
            sage: 15 + y
            x1+15
        """
        return Linear_Expression(self)+Linear_Expression(other)


    def __sub__(self, other):
        r"""
        Return the difference ``self`` - ``other``.

        INPUT:

        - ``self``, ``other`` -- anything convertible to
          ``Linear_Expression``: An integer, a :class:`Variable`, or a
          :class:`Linear_Expression`.

        OUTPUT:

        A :class:`Linear_Expression` representing ``self`` - ``other``.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable
            sage: x = Variable(0); y = Variable(1)
            sage: x - 15
            x0-15
            sage: 15 - y
            -x1+15
        """
        return Linear_Expression(self)-Linear_Expression(other)


    def __mul__(self, other):
        r"""
        Return the product ``self`` * ``other``.

        INPUT:

        - ``self``, ``other`` -- One must be an integer, the other a
          :class:`Variable`.

        OUTPUT:

        A :class:`Linear_Expression` representing ``self`` * ``other``.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable
            sage: x = Variable(0); y = Variable(1)
            sage: x * 15
            15*x0
            sage: 15 * y
            15*x1
        """
        if isinstance(self, Variable):
            return Linear_Expression(self) * other
        else:
            return Linear_Expression(other) * self


    def __pos__(self):
        r"""
        Return ``self`` as :class:`Linear_Expression`

        OUTPUT:

        The :class:`Linear_Expression` ``+self``

        EXAMPLES::

            sage: from sage.libs.ppl import Variable
            sage: x = Variable(0); x
            x0
            sage: +x
            x0
        """
        return Linear_Expression(self)


    def __neg__(self):
        r"""
        Return -``self`` as :class:`Linear_Expression`

        OUTPUT:

        The :class:`Linear_Expression` ``-self``

        EXAMPLES::

            sage: from sage.libs.ppl import Variable
            sage: x = Variable(0); x
            x0
            sage: -x
            -x0
        """
        return Linear_Expression(self)*(-1)


    def __richcmp__(self, other, op):
        """
        Construct :class:`Constraint` from equalities or inequalities.

        INPUT:

        - ``self``, ``other`` -- anything convertible to a
          :class:`Linear_Expression`

        - ``op`` -- the operation.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable
            sage: x = Variable(0)
            sage: y = Variable(1)
            sage: x <  y
            -x0+x1>0
            sage: x <= 0
            -x0>=0
            sage: x == y-y
            x0==0
            sage: x >= -2
            x0+2>=0
            sage: x >  0
            x0>0
            sage: 0 == 1    # watch out!
            False
            sage: 0*x == 1
            -1==0
        """
        return _make_Constraint_from_richcmp(self, other, op)


####################################################
### Linear_Expression ##############################
####################################################
cdef class Linear_Expression(object):
    r"""
    Wrapper for PPL's ``PPL_Linear_Expression`` class.

    INPUT:

    The constructor accepts zero, one, or two arguments.

    If there are two arguments ``Linear_Expression(a,b)``, they are
    interpreted as

    - ``a`` -- an iterable of integer coefficients, for example a
      list.

    - ``b`` -- an integer. The inhomogeneous term.

    A single argument ``Linear_Expression(arg)`` is interpreted as

    - ``arg`` -- something that determines a linear
      expression. Possibilities are:

      * a :class:`Variable`: The linear expression given by that
        variable.

      * a :class:`Linear_Expression`: The copy constructor.

      * an integer: Constructs the constant linear expression.

    No argument is the default constructor and returns the zero linear
    expression.

    OUTPUT:

    A :class:`Linear_Expression`

    EXAMPLES::

        sage: from sage.libs.ppl import Variable, Linear_Expression
        sage: Linear_Expression([1,2,3,4],5)
        x0+2*x1+3*x2+4*x3+5
        sage: Linear_Expression(10)
        10
        sage: Linear_Expression()
        0
        sage: Linear_Expression(10).inhomogeneous_term()
        10
        sage: x = Variable(123)
        sage: expr = x+1; expr
        x123+1
        sage: expr.OK()
        True
        sage: expr.coefficient(x)
        1
        sage: expr.coefficient( Variable(124) )
        0
    """

    cdef PPL_Linear_Expression *thisptr


    def __cinit__(self, *args):
        """
        The Cython constructor.

        See :class:`Linear_Expression` for documentation.

        TESTS::

            sage: from sage.libs.ppl import Linear_Expression
            sage: Linear_Expression(10)   # indirect doctest
            10
        """
        if len(args)==2:
            a = args[0]
            b = args[1]
            ex = Linear_Expression(0)
            for i in range(0,len(a)):
                ex = ex + Variable(i).__mul__(Integer(a[i]))
            arg = ex + b
        elif len(args)==1:
            arg = args[0]
        elif len(args)==0:
            self.thisptr = new PPL_Linear_Expression()
            return
        else:
            assert False, 'Cannot initialize with more than 2 arguments.'

        if isinstance(arg, Variable):
            v = <Variable>arg
            self.thisptr = new PPL_Linear_Expression(v.thisptr[0])
            return
        if isinstance(arg, Linear_Expression):
            e = <Linear_Expression>arg
            self.thisptr = new PPL_Linear_Expression(e.thisptr[0])
            return
        try:
            c = Integer(arg)
            self.thisptr = new PPL_Linear_Expression(PPL_Coefficient(c.value))
            return
        except ValueError:
            raise ValueError, 'Cannot initialize with '+str(args)+'.'


    def __dealloc__(self):
        """
        The Cython destructor.
        """
        del self.thisptr


    def space_dimension(self):
        """
        Return the dimension of the vector space necessary for the
        linear expression.

        OUTPUT:

        Integer.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable
            sage: x = Variable(0)
            sage: y = Variable(1)
            sage: ( x+y+1 ).space_dimension()
            2
            sage: ( x+y   ).space_dimension()
            2
            sage: (   y+1 ).space_dimension()
            2
            sage: ( x  +1 ).space_dimension()
            1
            sage: ( y+1-y ).space_dimension()
            2
        """
        return self.thisptr.space_dimension()


    def coefficient(self, Variable v):
        """
        Return the coefficient of the variable ``v``.

        INPUT:

        - ``v`` -- a :class:`Variable`.

        OUTPUT:

        An integer.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable
            sage: x = Variable(0)
            sage: e = 3*x+1
            sage: e.coefficient(x)
            3
        """
        cdef Integer c = Integer(0)
        mpz_set(c.value, self.thisptr.coefficient(v.thisptr[0]).get_mpz_t())
        return c


    def coefficients(self):
        """
        Return the coefficients of the linear expression.

        OUTPUT:

        A tuple of integers of length :meth:`space_dimension`.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable
            sage: x = Variable(0);  y = Variable(1)
            sage: e = 3*x+5*y+1
            sage: e.coefficients()
            (3, 5)
        """
        cdef int d = self.space_dimension()
        cdef int i
        cdef Integer c = Integer(0)
        coeffs = []
        for i in range(0,d):
            mpz_set(c.value, self.thisptr.coefficient(PPL_Variable(i)).get_mpz_t())
            coeffs.append(Integer(c))
        return tuple(coeffs)


    def inhomogeneous_term(self):
        """
        Return the inhomogeneous term of the linear expression.

        OUTPUT:

        Integer.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, Linear_Expression
            sage: Linear_Expression(10).inhomogeneous_term()
            10
        """
        cdef Integer c = Integer(0)
        mpz_set(c.value, self.thisptr.inhomogeneous_term().get_mpz_t())
        return c


    def is_zero(self):
        """
        Test if ``self`` is the zero linear expression.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, Linear_Expression
            sage: Linear_Expression(0).is_zero()
            True
            sage: Linear_Expression(10).is_zero()
            False
        """
        return self.thisptr.is_zero()


    def all_homogeneous_terms_are_zero(self):
        """
        Test if ``self`` is a constant linear expression.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, Linear_Expression
            sage: Linear_Expression(10).all_homogeneous_terms_are_zero()
            True
        """
        return self.thisptr.all_homogeneous_terms_are_zero()


    def ascii_dump(self):
        r"""
        Write an ASCII dump to stderr.

        EXAMPLES::

            sage: sage_cmd  = 'from sage.libs.ppl import Linear_Expression, Variable\n'
            sage: sage_cmd += 'x = Variable(0)\n'
            sage: sage_cmd += 'y = Variable(1)\n'
            sage: sage_cmd += 'e = 3*x+2*y+1\n'
            sage: sage_cmd += 'e.ascii_dump()\n'
            sage: from sage.tests.cmdline import test_executable
            sage: (out, err, ret) = test_executable(['sage', '-c', sage_cmd], timeout=100);  # long time, indirect doctest
            sage: print err  # long time
            size 3 1 3 2
        """
        self.thisptr.ascii_dump()


    def OK(self):
        """
        Check if all the invariants are satisfied.

        EXAMPLES::

            sage: from sage.libs.ppl import Linear_Expression, Variable
            sage: x = Variable(0)
            sage: y = Variable(1)
            sage: e = 3*x+2*y+1
            sage: e.OK()
            True
        """
        return self.thisptr.OK()


    def __repr__(self):
        r"""
        Return a string representation of the linear expression.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: from sage.libs.ppl import Linear_Expression, Variable
            sage: x = Variable(0)
            sage: y = Variable(1)
            sage: x+1
            x0+1
            sage: x+1-x
            1
            sage: 2*x
            2*x0
            sage: x-x-1
            -1
            sage: x-x
            0
        """
        s = ''
        first = True
        for i in range(0,self.space_dimension()):
            x = Variable(i)
            coeff = self.coefficient(x)
            if coeff==0: continue
            if first and coeff==1:
                s += '%r' % x
                first = False
            elif first and coeff==-1:
                s += '-%r' % x
                first = False
            elif first and coeff!=1:
                s += '%d*%r' % (coeff, x)
                first = False
            elif coeff==1:
                s += '+%r' % x
            elif coeff==-1:
                s += '-%r' % x
            else:
                s += '%+d*%r' % (coeff, x)
        inhomog = self.inhomogeneous_term()
        if inhomog!=0:
            if first:
                s += '%d' % inhomog
                first = False
            else:
                s += '%+d' % inhomog
        if first:
            s = '0'
        return s


    def __add__(self, other):
        r"""
        Add ``self`` and ``other``.

        INPUT:

        - ``self``, ``other`` -- anything that can be used to
          construct a :class:`Linear_Expression`. One of them, not
          necessarily ``self``, is guaranteed to be a
          :class:``Linear_Expression``, otherwise Python would not
          have called this method.

        OUTPUT:

        The sum as a :class:`Linear_Expression`

        EXAMPLES::

            sage: from sage.libs.ppl import Linear_Expression, Variable
            sage: x = Variable(0)
            sage: y = Variable(1)
            sage: 9+x+y+(1+x)+y+y
            2*x0+3*x1+10
        """
        cdef Linear_Expression lhs = Linear_Expression(self)
        cdef Linear_Expression rhs = Linear_Expression(other)
        cdef Linear_Expression result = Linear_Expression()
        result.thisptr[0] = lhs.thisptr[0] + rhs.thisptr[0]
        return result


    def __sub__(self, other):
        r"""
        Subtract ``other`` from ``self``.

        INPUT:

        - ``self``, ``other`` -- anything that can be used to
          construct a :class:`Linear_Expression`. One of them, not
          necessarily ``self``, is guaranteed to be a
          :class:``Linear_Expression``, otherwise Python would not
          have called this method.

        OUTPUT:

        The difference as a :class:`Linear_Expression`

        EXAMPLES::

            sage: from sage.libs.ppl import Linear_Expression, Variable
            sage: x = Variable(0)
            sage: y = Variable(1)
            sage: 9-x-y-(1-x)-y-y
            -3*x1+8
        """
        cdef Linear_Expression lhs = Linear_Expression(self)
        cdef Linear_Expression rhs = Linear_Expression(other)
        cdef Linear_Expression result = Linear_Expression()
        result.thisptr[0] = lhs.thisptr[0] - rhs.thisptr[0]
        return result


    def __mul__(self, other):
        r"""
        Multiply ``self`` with ``other``.

        INPUT:

        - ``self``, ``other`` -- anything that can be used to
          construct a :class:`Linear_Expression`. One of them, not
          necessarily ``self``, is guaranteed to be a
          :class:``Linear_Expression``, otherwise Python would not
          have called this method.

        OUTPUT:

        The product as a :class:`Linear_Expression`

        EXAMPLES::

            sage: from sage.libs.ppl import Linear_Expression, Variable
            sage: x = Variable(0)
            sage: y = Variable(1)
            sage: 8*(x+1)
            8*x0+8
            sage: y*8
            8*x1
        """
        cdef Linear_Expression e
        cdef Integer c
        if isinstance(self, Linear_Expression):
            e = <Linear_Expression>self
            c = Integer(other)
        else:
            e = <Linear_Expression>other
            c = Integer(self)

        cdef Linear_Expression result = Linear_Expression()
        result.thisptr[0] = e.thisptr[0] * PPL_Coefficient(c.value)
        return result


    def __pos__(self):
        """
        Return ``self``.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, Linear_Expression
            sage: +Linear_Expression(1)
            1
            sage: x = Variable(0)
            sage: +(x+1)
            x0+1
        """
        return self


    def __neg__(self):
        """
        Return the negative of ``self``.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, Linear_Expression
            sage: -Linear_Expression(1)
            -1
            sage: x = Variable(0)
            sage: -(x+1)
            -x0-1
        """
        return self*(-1)


    def __richcmp__(self, other, int op):
        """
        Construct :class:`Constraint`s

        EXAMPLES::

            sage: from sage.libs.ppl import Variable
            sage: x = Variable(0)
            sage: y = Variable(1)
            sage: x+1 <  y-2
            -x0+x1-3>0
            sage: x+1 <= y-2
            -x0+x1-3>=0
            sage: x+1 == y-2
            x0-x1+3==0
            sage: x+1 >= y-2
            x0-x1+3>=0
            sage: x+1 >  y-2
            x0-x1+3>0
        """
        return _make_Constraint_from_richcmp(self, other, op)


    def __reduce__(self):
        """
        Pickle object

        EXAMPLES::

            sage: from sage.libs.ppl import Linear_Expression
            sage: le = loads(dumps(Linear_Expression([1,2,3],4)))
            sage: le.coefficients() == (1,2,3)
            True
            sage: le.inhomogeneous_term() == 4
            True
        """
        return (Linear_Expression, (self.coefficients(), self.inhomogeneous_term()))


####################################################
### Generator ######################################
####################################################
cdef _wrap_Generator(PPL_Generator generator):
    """
    Wrap a C++ ``PPL_Generator`` into a Cython ``Generator``.
    """
    cdef Generator g = Generator(True)
    g.thisptr = new PPL_Generator(generator)
    return g


####################################################
# C++ static methods not supported
# Note that the PPL_Generator default constructor is private, hence we must return pointers
cdef extern from "ppl_shim.hh":
    PPL_Generator* new_line(PPL_Linear_Expression &e) except +ValueError
    PPL_Generator* new_ray(PPL_Linear_Expression &e) except +ValueError
    PPL_Generator* new_point(PPL_Linear_Expression &e, PPL_Coefficient d) except +ValueError
    PPL_Generator* new_closure_point(PPL_Linear_Expression &e, PPL_Coefficient d) except +ValueError
    PPL_Generator* new_MIP_optimizing_point(PPL_MIP_Problem &problem) except +ValueError


####################################################
cdef class Generator(object):
    r"""
    Wrapper for PPL's ``Generator`` class.

    An object of the class Generator is one of the following:

      * a line `\ell = (a_0, \dots, a_{n-1})^T`

      * a ray `r = (a_0, \dots, a_{n-1})^T`

      * a point `p = (\tfrac{a_0}{d}, \dots, \tfrac{a_{n-1}}{d})^T`

      * a closure point `c = (\tfrac{a_0}{d}, \dots, \tfrac{a_{n-1}}{d})^T`

    where `n` is the dimension of the space and, for points and
    closure points, `d` is the divisor.

    INPUT/OUTPUT:

    Use the helper functions :func:`line`, :func:`ray`, :func:`point`,
    and :func:`closure_point` to construct generators. Analogous class
    methods are also available, see :meth:`Generator.line`,
    :meth:`Generator.ray`, :meth:`Generator.point`,
    :meth:`Generator.closure_point`. Do not attempt to construct
    generators manually.

    .. NOTE::

        The generators are constructed from linear expressions. The
        inhomogeneous term is always silently discarded.

    EXAMPLES::

        sage: from sage.libs.ppl import Generator, Variable
        sage: x = Variable(0)
        sage: y = Variable(1)
        sage: Generator.line(5*x-2*y)
        line(5, -2)
        sage: Generator.ray(5*x-2*y)
        ray(5, -2)
        sage: Generator.point(5*x-2*y, 7)
        point(5/7, -2/7)
        sage: Generator.closure_point(5*x-2*y, 7)
        closure_point(5/7, -2/7)
    """

    cdef PPL_Generator *thisptr


    def __cinit__(self, do_not_construct_manually=False):
        """
        The Cython constructor.

        See :class:`Variable` for documentation.

        TESTS::

            sage: from sage.libs.ppl import Variable, line
            sage: x = Variable(0)
            sage: line(x)   # indirect doctest
            line(1)
        """
        assert(do_not_construct_manually)
        self.thisptr = NULL


    def __dealloc__(self):
        """
        The Cython destructor.
        """
        assert self.thisptr!=NULL, 'Do not construct Generators manually!'
        del self.thisptr


    @classmethod
    def line(cls, expression):
        """
        Construct a line.

        INPUT:

        - ``expression`` -- a :class:`Linear_Expression` or something
          convertible to it (:class:`Variable` or integer).

        OUTPUT:

        A new :class:`Generator` representing the line.

        Raises a ``ValueError` if the homogeneous part of
        ``expression`` represents the origin of the vector space.

        EXAMPLES::

            sage: from sage.libs.ppl import Generator, Variable
            sage: y = Variable(1)
            sage: Generator.line(2*y)
            line(0, 1)
            sage: Generator.line(y)
            line(0, 1)
            sage: Generator.line(1)
            Traceback (most recent call last):
            ...
            ValueError: PPL::line(e):
            e == 0, but the origin cannot be a line.
        """
        cdef Linear_Expression e = Linear_Expression(expression)
        # This does not work as Cython gets confused by the private default ctor
        #   return _wrap_Generator(PPL_line(e.thisptr[0]))
        # workaround follows
        cdef Generator g = Generator(True)
        try:
            g.thisptr = new_line(e.thisptr[0])
        except BaseException:
            # g.thisptr must be set to something valid or g.__dealloc__() will segfault
            g.thisptr = new_point(e.thisptr[0],PPL_Coefficient(1))
            raise
        return g


    @classmethod
    def ray(cls, expression):
        """
        Construct a ray.

        INPUT:

        - ``expression`` -- a :class:`Linear_Expression` or something
          convertible to it (:class:`Variable` or integer).

        OUTPUT:

        A new :class:`Generator` representing the ray.

        Raises a ``ValueError` if the homogeneous part of
        ``expression`` represents the origin of the vector space.

        EXAMPLES::

            sage: from sage.libs.ppl import Generator, Variable
            sage: y = Variable(1)
            sage: Generator.ray(2*y)
            ray(0, 1)
            sage: Generator.ray(y)
            ray(0, 1)
            sage: Generator.ray(1)
            Traceback (most recent call last):
            ...
            ValueError: PPL::ray(e):
            e == 0, but the origin cannot be a ray.
        """
        cdef Linear_Expression e = Linear_Expression(expression)
        # This does not work as Cython gets confused by the private default ctor
        #   return _wrap_Generator(PPL_ray(e.thisptr[0]))
        # workaround follows
        cdef Generator g = Generator(True)
        try:
            g.thisptr = new_ray(e.thisptr[0])
        except BaseException:
            # g.thisptr must be set to something valid or g.__dealloc__() will segfault
            g.thisptr = new_point(e.thisptr[0],PPL_Coefficient(1))
            raise
        return g


    @classmethod
    def point(cls, expression=0, divisor=1):
        """
        Construct a point.

        INPUT:

        - ``expression`` -- a :class:`Linear_Expression` or something
          convertible to it (:class:`Variable` or integer).

        - ``divisor`` -- an integer.

        OUTPUT:

        A new :class:`Generator` representing the point.

        Raises a ``ValueError` if ``divisor==0``.

        EXAMPLES::

            sage: from sage.libs.ppl import Generator, Variable
            sage: y = Variable(1)
            sage: Generator.point(2*y+7, 3)
            point(0/3, 2/3)
            sage: Generator.point(y+7, 3)
            point(0/3, 1/3)
            sage: Generator.point(7, 3)
            point()
            sage: Generator.point(0, 0)
            Traceback (most recent call last):
            ...
            ValueError: PPL::point(e, d):
            d == 0.
        """
        cdef Linear_Expression e = Linear_Expression(expression)
        cdef Integer d = Integer(divisor)
        # This does not work as Cython gets confused by the private default ctor
        #   return _wrap_Generator(PPL_point(e.thisptr[0], PPL_Coefficient(d.value)))
        # workaround follows
        cdef Generator g = Generator(True)
        try:
            g.thisptr = new_point(e.thisptr[0], PPL_Coefficient(d.value))
        except BaseException:
            # g.thisptr must be set to something valid or g.__dealloc__() will segfault
            g.thisptr = new_point(e.thisptr[0],PPL_Coefficient(1))
            raise
        return g


    @classmethod
    def closure_point(cls, expression=0, divisor=1):
        """
        Construct a closure point.

        A closure point is a point of the topological closure of a
        polyhedron that is not a point of the polyhedron itself.

        INPUT:

        - ``expression`` -- a :class:`Linear_Expression` or something
          convertible to it (:class:`Variable` or integer).

        - ``divisor`` -- an integer.

        OUTPUT:

        A new :class:`Generator` representing the point.

        Raises a ``ValueError` if ``divisor==0``.

        EXAMPLES::

            sage: from sage.libs.ppl import Generator, Variable
            sage: y = Variable(1)
            sage: Generator.closure_point(2*y+7, 3)
            closure_point(0/3, 2/3)
            sage: Generator.closure_point(y+7, 3)
            closure_point(0/3, 1/3)
            sage: Generator.closure_point(7, 3)
            closure_point()
            sage: Generator.closure_point(0, 0)
            Traceback (most recent call last):
            ...
            ValueError: PPL::closure_point(e, d):
            d == 0.
        """
        cdef Linear_Expression e = Linear_Expression(expression)
        cdef Integer d = Integer(divisor)
        # This does not work as Cython gets confused by the private default ctor
        #   return _wrap_Generator(PPL_closure_point(e.thisptr[0], PPL_Coefficient(d.value)))
        # workaround follows
        cdef Generator g = Generator(True)
        try:
            g.thisptr = new_closure_point(e.thisptr[0], PPL_Coefficient(d.value))
        except BaseException:
            # g.thisptr must be set to something valid or g.__dealloc__() will segfault
            g.thisptr = new_point(e.thisptr[0],PPL_Coefficient(1))
            raise
        return g


    def __repr__(self):
        """
        Return a string representation of the generator.

        OUTPUT:

        String.

        EXAMPLES::

            sage: from sage.libs.ppl import Generator, Variable
            sage: x = Variable(0)
            sage: y = Variable(1)
            sage: e = 2*x-y+5
            sage: Generator.line(e)
            line(2, -1)
            sage: Generator.ray(e)
            ray(2, -1)
            sage: Generator.point(e, 3)
            point(2/3, -1/3)
            sage: Generator.closure_point(e, 3)
            closure_point(2/3, -1/3)
        """
        t = self.type()
        if t=='line':
            s = 'line('
            div = ''
        elif t=='ray':
            s = 'ray('
            div = ''
        elif t=='point':
            s = 'point('
            div = '/'+str(self.divisor())
        elif t=='closure_point':
            s = 'closure_point('
            div = '/'+str(self.divisor())
        else:
            assert(False)

        for i in range(0,self.space_dimension()):
            if i>0:
                s += ', '
            s += str(self.coefficient(Variable(i))) + div

        s += ')'
        return s


    def space_dimension(self):
        r"""
        Return the dimension of the vector space enclosing ``self``.

        OUTPUT:

        Integer.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, point
            sage: x = Variable(0)
            sage: y = Variable(1)
            sage: point(x).space_dimension()
            1
            sage: point(y).space_dimension()
            2
        """
        return self.thisptr.space_dimension()


    def type(self):
        r"""
        Return the generator type of ``self``.

        OUTPUT:

        String. One of ``'line'``, ``'ray'``, ``'point'``, or
        ``'closure_point'``.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, point, closure_point, ray, line
            sage: x = Variable(0)
            sage: line(x).type()
            'line'
            sage: ray(x).type()
            'ray'
            sage: point(x,2).type()
            'point'
            sage: closure_point(x,2).type()
            'closure_point'
        """
        t = self.thisptr.type()
        if t==LINE:
            return 'line'
        elif t==RAY:
            return 'ray'
        elif t==POINT:
            return 'point'
        elif t==CLOSURE_POINT:
            return 'closure_point'
        assert False


    def is_line(self):
        r"""
        Test whether ``self`` is a line.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, point, closure_point, ray, line
            sage: x = Variable(0)
            sage: line(x).is_line()
            True
            sage: ray(x).is_line()
            False
            sage: point(x,2).is_line()
            False
            sage: closure_point(x,2).is_line()
            False
        """
        return self.thisptr.is_line()


    def is_ray(self):
        r"""
        Test whether ``self`` is a ray.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, point, closure_point, ray, line
            sage: x = Variable(0)
            sage: line(x).is_ray()
            False
            sage: ray(x).is_ray()
            True
            sage: point(x,2).is_ray()
            False
            sage: closure_point(x,2).is_ray()
            False
        """
        return self.thisptr.is_ray()


    def is_line_or_ray(self):
        r"""
        Test whether ``self`` is a line or a ray.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, point, closure_point, ray, line
            sage: x = Variable(0)
            sage: line(x).is_line_or_ray()
            True
            sage: ray(x).is_line_or_ray()
            True
            sage: point(x,2).is_line_or_ray()
            False
            sage: closure_point(x,2).is_line_or_ray()
            False
        """
        return self.thisptr.is_line_or_ray()


    def is_point(self):
        r"""
        Test whether ``self`` is a point.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, point, closure_point, ray, line
            sage: x = Variable(0)
            sage: line(x).is_point()
            False
            sage: ray(x).is_point()
            False
            sage: point(x,2).is_point()
            True
            sage: closure_point(x,2).is_point()
            False
        """
        return self.thisptr.is_point()


    def is_closure_point(self):
        r"""
        Test whether ``self`` is a closure point.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, point, closure_point, ray, line
            sage: x = Variable(0)
            sage: line(x).is_closure_point()
            False
            sage: ray(x).is_closure_point()
            False
            sage: point(x,2).is_closure_point()
            False
            sage: closure_point(x,2).is_closure_point()
            True
        """
        return self.thisptr.is_closure_point()


    def coefficient(self, Variable v):
        """
        Return the coefficient of the variable ``v``.

        INPUT:

        - ``v`` -- a :class:`Variable`.

        OUTPUT:

        An integer.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, line
            sage: x = Variable(0)
            sage: line = line(3*x+1)
            sage: line
            line(1)
            sage: line.coefficient(x)
            1
        """
        cdef Integer c = Integer(0)
        mpz_set(c.value, self.thisptr.coefficient(v.thisptr[0]).get_mpz_t())
        return c


    def coefficients(self):
        """
        Return the coefficients of the generator.

        See also :meth:`coefficient`.

        OUTPUT:

        A tuple of integers of length :meth:`space_dimension`.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, point
            sage: x = Variable(0);  y = Variable(1)
            sage: p = point(3*x+5*y+1, 2); p
            point(3/2, 5/2)
            sage: p.coefficients()
            (3, 5)
        """
        cdef int d = self.space_dimension()
        cdef int i
        cdef Integer c = Integer(0)
        coeffs = []
        for i in range(0,d):
            mpz_set(c.value, self.thisptr.coefficient(PPL_Variable(i)).get_mpz_t())
            coeffs.append(Integer(c))
        return tuple(coeffs)


    def divisor(self):
        """
        If ``self`` is either a point or a closure point, return its
        divisor.

        OUTPUT:

        An integer. If ``self`` is a ray or a line, raises
        ``ValueError``.

        EXAMPLES::

            sage: from sage.libs.ppl import Generator, Variable
            sage: x = Variable(0)
            sage: y = Variable(1)
            sage: point = Generator.point(2*x-y+5)
            sage: point.divisor()
            1
            sage: line = Generator.line(2*x-y+5)
            sage: line.divisor()
            Traceback (most recent call last):
            ...
            ValueError: PPL::Generator::divisor():
            *this is neither a point nor a closure point.
        """
        cdef Integer c = Integer(0)
        mpz_set(c.value, self.thisptr.divisor().get_mpz_t())
        return c


    def is_equivalent_to(self, Generator g):
        r"""
        Test whether ``self`` and ``g`` are equivalent.

        INPUT:

        - ``g`` -- a :class:`Generator`.

        OUTPUT:

        Boolean. Returns ``True`` if and only if ``self`` and ``g``
        are equivalent generators.

        Note that generators having different space dimensions are not
        equivalent.

        EXAMPLES::

            sage: from sage.libs.ppl import Generator, Variable, point, line
            sage: x = Variable(0)
            sage: y = Variable(1)
            sage: point(2*x    , 2).is_equivalent_to( point(x) )
            True
            sage: point(2*x+0*y, 2).is_equivalent_to( point(x) )
            False
            sage: line(4*x).is_equivalent_to(line(x))
            True
        """
        return self.thisptr.is_equivalent_to(g.thisptr[0])


    def ascii_dump(self):
        r"""
        Write an ASCII dump to stderr.

        EXAMPLES::

            sage: sage_cmd  = 'from sage.libs.ppl import Linear_Expression, Variable, point\n'
            sage: sage_cmd += 'x = Variable(0)\n'
            sage: sage_cmd += 'y = Variable(1)\n'
            sage: sage_cmd += 'p = point(3*x+2*y)\n'
            sage: sage_cmd += 'p.ascii_dump()\n'
            sage: from sage.tests.cmdline import test_executable
            sage: (out, err, ret) = test_executable(['sage', '-c', sage_cmd], timeout=100);  # long time, indirect doctest
            sage: print err  # long time
            size 3 1 3 2 P (C)
        """
        self.thisptr.ascii_dump()


    def OK(self):
        """
        Check if all the invariants are satisfied.

        EXAMPLES::

            sage: from sage.libs.ppl import Linear_Expression, Variable
            sage: x = Variable(0)
            sage: y = Variable(1)
            sage: e = 3*x+2*y+1
            sage: e.OK()
            True
        """
        return self.thisptr.OK()


    def __reduce__(self):
        """
        Pickle object.

        TESTS::

            sage: from sage.libs.ppl import Generator, Variable, line, ray, point, closure_point
            sage: x = Variable(0); y = Variable(1);
            sage: loads(dumps(Generator.point(2*x+7*y, 3)))
            point(2/3, 7/3)
            sage: loads(dumps(Generator.closure_point(2*x+7*y, 3)))
            closure_point(2/3, 7/3)
            sage: loads(dumps(Generator.line(2*x+7*y)))
            line(2, 7)
            sage: loads(dumps(Generator.ray(2*x+7*y)))
            ray(2, 7)
        """
        t = self.thisptr.type()
        le = Linear_Expression(self.coefficients(), 0)
        if t==LINE:
            return (line, (le,))
        elif t==RAY:
            return (ray, (le,))
        elif t==POINT:
            return (point, (le, self.divisor()))
        elif t==CLOSURE_POINT:
            return (closure_point, (le, self.divisor()))
        assert False



####################################################
def line(expression):
    """
    Constuct a line.

    See :meth:`Generator.line` for documentation.

    EXAMPLES::

        sage: from sage.libs.ppl import Variable, line
        sage: y = Variable(1)
        sage: line(2*y)
        line(0, 1)
    """
    return Generator.line(expression)


####################################################
def ray(expression):
    """
    Constuct a ray.

    See :meth:`Generator.ray` for documentation.

    EXAMPLES::

        sage: from sage.libs.ppl import Variable, ray
        sage: y = Variable(1)
        sage: ray(2*y)
        ray(0, 1)
    """
    return Generator.ray(expression)


####################################################
def point(expression=0, divisor=1):
    """
    Constuct a point.

    See :meth:`Generator.point` for documentation.

    EXAMPLES::

        sage: from sage.libs.ppl import Variable, point
        sage: y = Variable(1)
        sage: point(2*y, 5)
        point(0/5, 2/5)
    """
    return Generator.point(expression, divisor)


####################################################
def closure_point(expression=0, divisor=1):
    """
    Constuct a closure point.

    See :meth:`Generator.closure_point` for documentation.

    EXAMPLES::

        sage: from sage.libs.ppl import Variable, closure_point
        sage: y = Variable(1)
        sage: closure_point(2*y, 5)
        closure_point(0/5, 2/5)
    """
    return Generator.closure_point(expression, divisor)



####################################################
### Generator_System  ##############################
####################################################
cdef _wrap_Generator_System(PPL_Generator_System generator_system):
    """
    Wrap a C++ ``PPL_Generator_System`` into a Cython ``Generator_System``.
    """
    cdef Generator_System gs = Generator_System()
    del gs.thisptr
    gs.thisptr = new PPL_Generator_System(generator_system)
    return gs


####################################################
cdef class Generator_System(_mutable_or_immutable):
    """
    Wrapper for PPL's ``Generator_System`` class.

    An object of the class Generator_System is a system of generators,
    i.e., a multiset of objects of the class Generator (lines, rays,
    points and closure points). When inserting generators in a system,
    space dimensions are automatically adjusted so that all the
    generators in the system are defined on the same vector space. A
    system of generators which is meant to define a non-empty
    polyhedron must include at least one point: the reason is that
    lines, rays and closure points need a supporting point (lines and
    rays only specify directions while closure points only specify
    points in the topological closure of the NNC polyhedron).

    EXAMPLES::

        sage: from sage.libs.ppl import Generator_System, Variable, line, ray, point, closure_point
        sage: x = Variable(0)
        sage: y = Variable(1)
        sage: gs = Generator_System( line(5*x-2*y) )
        sage: gs.insert( ray(6*x-3*y) )
        sage: gs.insert( point(2*x-7*y, 5) )
        sage: gs.insert( closure_point(9*x-1*y, 2) )
        sage: gs
        Generator_System {line(5, -2), ray(2, -1), point(2/5, -7/5), closure_point(9/2, -1/2)}
    """

    cdef PPL_Generator_System *thisptr


    def __cinit__(self, arg=None):
        """
        The Cython constructor.

        See :class:`Generator_System` for documentation.

        TESTS::

            sage: from sage.libs.ppl import Generator_System
            sage: Generator_System()   # indirect doctest
            Generator_System {}
        """
        if arg is None:
            self.thisptr = new PPL_Generator_System()
            return
        if isinstance(arg, Generator):
            g = <Generator>arg
            self.thisptr = new PPL_Generator_System(g.thisptr[0])
            return
        if isinstance(arg, Generator_System):
            gs = <Generator_System>arg
            self.thisptr = new PPL_Generator_System(gs.thisptr[0])
            return
        if isinstance(arg, (list,tuple)):
            self.thisptr = new PPL_Generator_System()
            for generator in arg:
                self.insert(generator)
            return
        raise ValueError, 'Cannot initialize with '+str(arg)+'.'


    def __dealloc__(self):
        """
        The Cython destructor.
        """
        del self.thisptr


    def space_dimension(self):
        r"""
        Return the dimension of the vector space enclosing ``self``.

        OUTPUT:

        Integer.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, Generator_System, point
            sage: x = Variable(0)
            sage: gs = Generator_System( point(3*x) )
            sage: gs.space_dimension()
            1
        """
        return self.thisptr.space_dimension()


    def clear(self):
        r"""
        Removes all generators from the generator system and sets its
        space dimension to 0.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, Generator_System, point
            sage: x = Variable(0)
            sage: gs = Generator_System( point(3*x) ); gs
            Generator_System {point(3/1)}
            sage: gs.clear()
            sage: gs
            Generator_System {}
        """
        self.assert_mutable('The Generator_System is not mutable!')
        self.thisptr.clear()


    def insert(self, Generator g):
        """
        Insert ``g`` into the generator system.

        The number of space dimensions of ``self`` is increased, if needed.

        INPUT:

        - ``g`` -- a :class:`Generator`.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, Generator_System, point
            sage: x = Variable(0)
            sage: gs = Generator_System( point(3*x) )
            sage: gs.insert( point(-3*x) )
            sage: gs
            Generator_System {point(3/1), point(-3/1)}
        """
        self.assert_mutable('The Generator_System is not mutable!')
        self.thisptr.insert(g.thisptr[0])


    def empty(self):
        """
        Return ``True`` if and only if ``self`` has no generators.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, Generator_System, point
            sage: x = Variable(0)
            sage: gs = Generator_System()
            sage: gs.empty()
            True
            sage: gs.insert( point(-3*x) )
            sage: gs.empty()
            False
        """
        return self.thisptr.empty()


    def ascii_dump(self):
        r"""
        Write an ASCII dump to stderr.

        EXAMPLES::

            sage: sage_cmd  = 'from sage.libs.ppl import Generator_System, point, Variable\n'
            sage: sage_cmd += 'x = Variable(0)\n'
            sage: sage_cmd += 'y = Variable(1)\n'
            sage: sage_cmd += 'gs = Generator_System( point(3*x+2*y+1) )\n'
            sage: sage_cmd += 'gs.ascii_dump()\n'
            sage: from sage.tests.cmdline import test_executable
            sage: (out, err, ret) = test_executable(['sage', '-c', sage_cmd], timeout=100);  # long time, indirect doctest
            sage: print err  # long time
            topology NECESSARILY_CLOSED
            1 x 2 SPARSE (sorted)
            index_first_pending 1
            size 3 1 3 2 P (C)
        """
        self.thisptr.ascii_dump()


    def OK(self):
        """
        Check if all the invariants are satisfied.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, Generator_System, point
            sage: x = Variable(0)
            sage: y = Variable(1)
            sage: gs = Generator_System( point(3*x+2*y+1) )
            sage: gs.OK()
            True
        """
        return self.thisptr.OK()


    def __len__(self):
        """
        Return the number of generators in the system.

            sage: from sage.libs.ppl import Variable, Generator_System, point
            sage: x = Variable(0)
            sage: y = Variable(1)
            sage: gs = Generator_System()
            sage: gs.insert(point(3*x+2*y))
            sage: gs.insert(point(x))
            sage: gs.insert(point(y))
            sage: len(gs)
            3
        """
        return sum([1 for g in self])


    def __iter__(self):
        """
        Iterate through the generators of the system.

        EXAMPLES::

            sage: from sage.libs.ppl import Generator_System, Variable, point
            sage: x = Variable(0)
            sage: gs = Generator_System(point(3*x))
            sage: iter = gs.__iter__()
            sage: iter.next()
            point(3/1)
        """
        return Generator_System_iterator(self)


    def __getitem__(self, int k):
        """
        Return the ``k``-th generator.

        The correct way to read the individual generators is to
        iterate over the generator system. This method is for
        convenience only.

        INPUT:

        - ``k`` -- integer. The index of the generator.

        OUTPUT:

        The ``k``-th constraint of the generator system.

        EXAMPLES::

            sage: from sage.libs.ppl import Generator_System, Variable, point
            sage: x = Variable(0)
            sage: gs = Generator_System()
            sage: gs.insert(point(3*x))
            sage: gs.insert(point(-2*x))
            sage: gs
            Generator_System {point(3/1), point(-2/1)}
            sage: gs[0]
            point(3/1)
            sage: gs[1]
            point(-2/1)
        """
        if k < 0:
            raise IndexError('index must be nonnegative')
        iterator = self.__iter__()
        try:
            for i in range(k):
                iterator.next()
        except StopIteration:
            raise IndexError('index is past-the-end')
        return iterator.next()


    def __repr__(self):
        r"""
        Return a string representation of the generator system.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: from sage.libs.ppl import Generator_System, Variable, point, ray
            sage: x = Variable(0)
            sage: y = Variable(1)
            sage: gs = Generator_System(point(3*x+2*y+1))
            sage: gs.insert(ray(x))
            sage: gs.__repr__()
            'Generator_System {point(3/1, 2/1), ray(1, 0)}'
        """
        s = 'Generator_System {'
        s += ', '.join([ g.__repr__() for g in self ])
        s += '}'
        return s


    def __reduce__(self):
        """
        Pickle object.

        TESTS::

            sage: from sage.libs.ppl import Generator_System, Variable, point, ray
            sage: x = Variable(0)
            sage: y = Variable(1)
            sage: gs = Generator_System((point(3*x+2*y+1), ray(x)));  gs
            Generator_System {point(3/1, 2/1), ray(1, 0)}
            sage: loads(dumps(gs))
            Generator_System {point(3/1, 2/1), ray(1, 0)}
        """
        return (Generator_System, (tuple(self), ))



####################################################
### Generator_System_iterator ######################
####################################################
cdef extern from "ppl_shim.hh":
    ctypedef void* gs_iterator_ptr
    cdef gs_iterator_ptr init_gs_iterator(PPL_Generator_System &gs)
    cdef PPL_Generator next_gs_iterator(gs_iterator_ptr)
    cdef bint is_end_gs_iterator(PPL_Generator_System &gs, gs_iterator_ptr gsi_ptr)
    cdef void delete_gs_iterator(gs_iterator_ptr)


####################################################
cdef class Generator_System_iterator(object):
    """
    Wrapper for PPL's ``Generator_System::const_iterator`` class.

    EXAMPLES::

        sage: from sage.libs.ppl import Generator_System, Variable, line, ray, point, closure_point, Generator_System_iterator
        sage: x = Variable(0)
        sage: y = Variable(1)
        sage: gs = Generator_System( line(5*x-2*y) )
        sage: gs.insert( ray(6*x-3*y) )
        sage: gs.insert( point(2*x-7*y, 5) )
        sage: gs.insert( closure_point(9*x-1*y, 2) )
        sage: Generator_System_iterator(gs).next()
        line(5, -2)
        sage: list(gs)
        [line(5, -2), ray(2, -1), point(2/5, -7/5), closure_point(9/2, -1/2)]
    """

    cdef Generator_System gs
    cdef gs_iterator_ptr gsi_ptr


    def __cinit__(self, Generator_System gs):
        r"""
        The Cython constructor.

        TESTS::

            sage: from sage.libs.ppl import Generator_System, Generator_System_iterator
            sage: iter = Generator_System_iterator(Generator_System())   # indirect doctest
        """
        self.gs = gs
        self.gsi_ptr = init_gs_iterator(gs.thisptr[0])


    def __dealloc__(self):
        """
        The Cython destructor.
        """
        delete_gs_iterator(self.gsi_ptr)


    def __next__(Generator_System_iterator self):
        r"""
        The next iteration.

        OUTPUT:

        A :class:`Generator`.

        EXAMPLES::

            sage: from sage.libs.ppl import Generator_System, Variable, point, Generator_System_iterator
            sage: x = Variable(0)
            sage: y = Variable(1)
            sage: gs = Generator_System( point(5*x-2*y) )
            sage: Generator_System_iterator(gs).next()
            point(5/1, -2/1)
        """
        if is_end_gs_iterator((<Generator_System>self.gs).thisptr[0], self.gsi_ptr):
            raise StopIteration
        return _wrap_Generator(next_gs_iterator(self.gsi_ptr))


####################################################
### Constraint ######################################
####################################################
cdef _wrap_Constraint(PPL_Constraint constraint):
    """
    Wrap a C++ ``PPL_Constraint`` into a Cython ``Constraint``.
    """
    cdef Constraint c = Constraint(True)
    c.thisptr = new PPL_Constraint(constraint)
    return c


####################################################
cdef _make_Constraint_from_richcmp(lhs_, rhs_, op):
    cdef Linear_Expression lhs = Linear_Expression(lhs_)
    cdef Linear_Expression rhs = Linear_Expression(rhs_)
    if op==0:      # <   0
        return _wrap_Constraint(lhs.thisptr[0] <  rhs.thisptr[0])
    elif op==1:    # <=  1
        return _wrap_Constraint(lhs.thisptr[0] <= rhs.thisptr[0])
    elif op==2:    # ==  2
        return _wrap_Constraint(lhs.thisptr[0] == rhs.thisptr[0])
    elif op==4:    # >   4
        return _wrap_Constraint(lhs.thisptr[0] >  rhs.thisptr[0])
    elif op==5:    # >=  5
        return _wrap_Constraint(lhs.thisptr[0] >= rhs.thisptr[0])
    elif op==3:    # !=  3
        raise NotImplementedError
    else:
        assert(False)


####################################################
cdef class Constraint(object):
    """
    Wrapper for PPL's ``Constraint`` class.

    An object of the class ``Constraint`` is either:

    * an equality `\sum_{i=0}^{n-1} a_i x_i + b = 0`

    * a non-strict inequality `\sum_{i=0}^{n-1} a_i x_i + b \geq 0`

    * a strict inequality `\sum_{i=0}^{n-1} a_i x_i + b > 0`

    where `n` is the dimension of the space, `a_i` is the integer
    coefficient of variable `x_i`, and `b_i` is the integer
    inhomogeneous term.

    INPUT/OUTPUT:

    You construct constraints by writing inequalities in
    :class:`Linear_Expression`. Do not attempt to manually construct
    constraints.

    EXAMPLES::

        sage: from sage.libs.ppl import Constraint, Variable, Linear_Expression
        sage: x = Variable(0)
        sage: y = Variable(1)
        sage: 5*x-2*y >  x+y-1
        4*x0-3*x1+1>0
        sage: 5*x-2*y >= x+y-1
        4*x0-3*x1+1>=0
        sage: 5*x-2*y == x+y-1
        4*x0-3*x1+1==0
        sage: 5*x-2*y <= x+y-1
        -4*x0+3*x1-1>=0
        sage: 5*x-2*y <  x+y-1
        -4*x0+3*x1-1>0
        sage: x > 0
        x0>0

    Special care is needed if the left hand side is a constant::

        sage: 0 == 1    # watch out!
        False
        sage: Linear_Expression(0) == 1
        -1==0
    """

    cdef PPL_Constraint *thisptr


    def __cinit__(self, do_not_construct_manually=False):
        """
        The Cython constructor.

        See :class:`Constraint` for documentation.

        TESTS::

            sage: from sage.libs.ppl import Constraint, Variable, Linear_Expression
            sage: x = Variable(0)
            sage: x>0   # indirect doctest
            x0>0
        """
        assert(do_not_construct_manually)
        self.thisptr = NULL


    def __dealloc__(self):
        """
        The Cython destructor.
        """
        assert self.thisptr!=NULL, 'Do not construct Constraints manually!'
        del self.thisptr


    def __repr__(self):
        """
        Return a string representation of the constraint.

        OUTPUT:

        String.

        EXAMPLES::

            sage: from sage.libs.ppl import Constraint, Variable
            sage: x = Variable(0)
            sage: y = Variable(1)
            sage: (2*x-y+5 >  x).__repr__()
            'x0-x1+5>0'
            sage: (2*x-y+5 == x).__repr__()
            'x0-x1+5==0'
            sage: (2*x-y+5 >= x).__repr__()
            'x0-x1+5>=0'
        """
        e = sum([ self.coefficient(x)*x
                  for x in [Variable(i)
                            for i in range(0,self.space_dimension())] ])
        e += self.inhomogeneous_term()
        s = e.__repr__()
        t = self.type()
        if t=='equality':
            s += '==0'
        elif t=='nonstrict_inequality':
            s += '>=0'
        elif t=='strict_inequality':
            s += '>0'
        else:
            assert(False)
        return s


    def space_dimension(self):
        r"""
        Return the dimension of the vector space enclosing ``self``.

        OUTPUT:

        Integer.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable
            sage: x = Variable(0)
            sage: y = Variable(1)
            sage: (x>=0).space_dimension()
            1
            sage: (y==1).space_dimension()
            2
        """
        return self.thisptr.space_dimension()


    def type(self):
        r"""
        Return the constraint type of ``self``.

        OUTPUT:

        String. One of ``'equality'``, ``'nonstrict_inequality'``, or
        ``'strict_inequality'``.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable
            sage: x = Variable(0)
            sage: (x==0).type()
            'equality'
            sage: (x>=0).type()
            'nonstrict_inequality'
            sage: (x>0).type()
            'strict_inequality'
        """
        t = self.thisptr.type()
        if t==EQUALITY:
            return 'equality'
        elif t==NONSTRICT_INEQUALITY:
            return 'nonstrict_inequality'
        elif t==STRICT_INEQUALITY:
            return 'strict_inequality'


    def is_equality(self):
        r"""
        Test whether ``self`` is an equality.

        OUTPUT:

        Boolean. Returns ``True`` if and only if ``self`` is an
        equality constraint.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable
            sage: x = Variable(0)
            sage: (x==0).is_equality()
            True
            sage: (x>=0).is_equality()
            False
            sage: (x>0).is_equality()
            False
        """
        return self.thisptr.is_equality()


    def is_inequality(self):
        r"""
        Test whether ``self`` is an inequality.

        OUTPUT:

        Boolean. Returns ``True`` if and only if ``self`` is an
        inequality constraint, either strict or non-strict.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable
            sage: x = Variable(0)
            sage: (x==0).is_inequality()
            False
            sage: (x>=0).is_inequality()
            True
            sage: (x>0).is_inequality()
            True
        """
        return self.thisptr.is_inequality()


    def is_nonstrict_inequality(self):
        r"""
        Test whether ``self`` is a non-strict inequality.

        OUTPUT:

        Boolean. Returns ``True`` if and only if ``self`` is an
        non-strict inequality constraint.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable
            sage: x = Variable(0)
            sage: (x==0).is_nonstrict_inequality()
            False
            sage: (x>=0).is_nonstrict_inequality()
            True
            sage: (x>0).is_nonstrict_inequality()
            False
        """
        return self.thisptr.is_nonstrict_inequality()


    def is_strict_inequality(self):
        r"""
        Test whether ``self`` is a strict inequality.

        OUTPUT:

        Boolean. Returns ``True`` if and only if ``self`` is an
        strict inequality constraint.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable
            sage: x = Variable(0)
            sage: (x==0).is_strict_inequality()
            False
            sage: (x>=0).is_strict_inequality()
            False
            sage: (x>0).is_strict_inequality()
            True
        """
        return self.thisptr.is_strict_inequality()


    def coefficient(self, Variable v):
        """
        Return the coefficient of the variable ``v``.

        INPUT:

        - ``v`` -- a :class:`Variable`.

        OUTPUT:

        An integer.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable
            sage: x = Variable(0)
            sage: ineq = (3*x+1 > 0)
            sage: ineq.coefficient(x)
            3
        """
        cdef Integer c = Integer(0)
        mpz_set(c.value, self.thisptr.coefficient(v.thisptr[0]).get_mpz_t())
        return c


    def coefficients(self):
        """
        Return the coefficients of the constraint.

        See also :meth:`coefficient`.

        OUTPUT:

        A tuple of integers of length :meth:`space_dimension`.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable
            sage: x = Variable(0);  y = Variable(1)
            sage: ineq = ( 3*x+5*y+1 ==  2);  ineq
            3*x0+5*x1-1==0
            sage: ineq.coefficients()
            (3, 5)
        """
        cdef int d = self.space_dimension()
        cdef int i
        cdef Integer c = Integer(0)
        coeffs = []
        for i in range(0,d):
            mpz_set(c.value, self.thisptr.coefficient(PPL_Variable(i)).get_mpz_t())
            coeffs.append(Integer(c))
        return tuple(coeffs)


    def inhomogeneous_term(self):
        """
        Return the inhomogeneous term of the constraint.

        OUTPUT:

        Integer.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable
            sage: y = Variable(1)
            sage: ineq = ( 10+y > 9 )
            sage: ineq
            x1+1>0
            sage: ineq.inhomogeneous_term()
            1
        """
        cdef Integer c = Integer(0)
        mpz_set(c.value, self.thisptr.inhomogeneous_term().get_mpz_t())
        return c


    def is_tautological(self):
        r"""
        Test whether ``self`` is a tautological constraint.

        A tautology can have either one of the following forms:

        * an equality: `\sum 0 x_i + 0 = 0`,

        * a non-strict inequality: `\sum 0 x_i + b \geq 0` with `b\geq 0`, or

        * a strict inequality: `\sum 0 x_i + b > 0` with `b> 0`.

        OUTPUT:

        Boolean. Returns ``True`` if and only if ``self`` is a
        tautological constraint.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable
            sage: x = Variable(0)
            sage: (x==0).is_tautological()
            False
            sage: (0*x>=0).is_tautological()
            True
        """
        return self.thisptr.is_tautological()


    def is_inconsistent(self):
        r"""
        Test whether ``self`` is an inconsistent constraint, that is, always false.

        An inconsistent constraint can have either one of the
        following forms:

        * an equality: `\sum 0 x_i + b = 0` with `b\not=0`,

        * a non-strict inequality: `\sum 0 x_i + b \geq 0` with `b< 0`, or

        * a strict inequality: `\sum 0 x_i + b > 0` with `b\leq 0`.

        OUTPUT:

        Boolean. Returns ``True`` if and only if ``self`` is an
        inconsistent constraint.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable
            sage: x = Variable(0)
            sage: (x==1).is_inconsistent()
            False
            sage: (0*x>=1).is_inconsistent()
            True
        """
        return self.thisptr.is_inconsistent()


    def is_equivalent_to(self, Constraint c):
        r"""
        Test whether ``self`` and ``c`` are equivalent.

        INPUT:

        - ``c`` -- a :class:`Constraint`.

        OUTPUT:

        Boolean. Returns ``True`` if and only if ``self`` and ``c``
        are equivalent constraints.

        Note that constraints having different space dimensions are
        not equivalent. However, constraints having different types
        may nonetheless be equivalent, if they both are tautologies or
        inconsistent.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, Linear_Expression
            sage: x = Variable(0)
            sage: y = Variable(1)
            sage: ( x>0 ).is_equivalent_to( Linear_Expression(0)<x )
            True
            sage: ( x>0 ).is_equivalent_to( 0*y<x )
            False
            sage: ( 0*x>1 ).is_equivalent_to( 0*x==-2 )
            True
        """
        return self.thisptr.is_equivalent_to(c.thisptr[0])


    def ascii_dump(self):
        r"""
        Write an ASCII dump to stderr.

        EXAMPLES::

            sage: sage_cmd  = 'from sage.libs.ppl import Linear_Expression, Variable\n'
            sage: sage_cmd += 'x = Variable(0)\n'
            sage: sage_cmd += 'y = Variable(1)\n'
            sage: sage_cmd += 'e = (3*x+2*y+1 > 0)\n'
            sage: sage_cmd += 'e.ascii_dump()\n'
            sage: from sage.tests.cmdline import test_executable
            sage: (out, err, ret) = test_executable(['sage', '-c', sage_cmd], timeout=100);  # long time, indirect doctest
            sage: print err  # long time
            size 4 1 3 2 -1 > (NNC)
        """
        self.thisptr.ascii_dump()


    def OK(self):
        """
        Check if all the invariants are satisfied.

        EXAMPLES::

            sage: from sage.libs.ppl import Linear_Expression, Variable
            sage: x = Variable(0)
            sage: y = Variable(1)
            sage: ineq = (3*x+2*y+1>=0)
            sage: ineq.OK()
            True
        """
        return self.thisptr.OK()


    def __reduce__(self):
        """
        Pickle object.

        EXAMPLES::

            sage: from sage.libs.ppl import Linear_Expression, Variable
            sage: x = Variable(0)
            sage: y = Variable(1)
            sage: loads(dumps(3*x+2*y+1>=5))
            3*x0+2*x1-4>=0
            sage: loads(dumps(3*x+2*y+1>5))
            3*x0+2*x1-4>0
            sage: loads(dumps(3*x+2*y+1==5))
            3*x0+2*x1-4==0
        """
        le = Linear_Expression(self.coefficients(), self.inhomogeneous_term())
        if self.is_nonstrict_inequality():
            return (inequality, (le, ))
        elif self.is_strict_inequality():
            return (strict_inequality, (le, ))
        elif self.is_equality():
            return (equation, (le, ))
        else:
            assert False



####################################################
def inequality(expression):
    """
    Constuct an inequality.

    INPUT:

    - ``expression`` -- a :class:`Linear_Expression`.

    OUTPUT:

    The inequality ``expression`` >= 0.

    EXAMPLES::

        sage: from sage.libs.ppl import Variable, inequality
        sage: y = Variable(1)
        sage: 2*y+1 >= 0
        2*x1+1>=0
        sage: inequality(2*y+1)
        2*x1+1>=0
    """
    return expression >= 0


####################################################
def strict_inequality(expression):
    """
    Constuct a strict inequality.

    INPUT:

    - ``expression`` -- a :class:`Linear_Expression`.

    OUTPUT:

    The inequality ``expression`` > 0.

    EXAMPLES::

        sage: from sage.libs.ppl import Variable, strict_inequality
        sage: y = Variable(1)
        sage: 2*y+1 > 0
        2*x1+1>0
        sage: strict_inequality(2*y+1)
        2*x1+1>0
    """
    return expression > 0


####################################################
def equation(expression):
    """
    Constuct an equation.

    INPUT:

    - ``expression`` -- a :class:`Linear_Expression`.

    OUTPUT:

    The equation ``expression`` == 0.

    EXAMPLES::

        sage: from sage.libs.ppl import Variable, equation
        sage: y = Variable(1)
        sage: 2*y+1 == 0
        2*x1+1==0
        sage: equation(2*y+1)
        2*x1+1==0
    """
    return expression == 0



####################################################
### Constraint_System  ##############################
####################################################
cdef _wrap_Constraint_System(PPL_Constraint_System constraint_system):
    """
    Wrap a C++ ``PPL_Constraint_System`` into a Cython ``Constraint_System``.
    """
    cdef Constraint_System cs = Constraint_System()
    del cs.thisptr
    cs.thisptr = new PPL_Constraint_System(constraint_system)
    return cs


####################################################
cdef class Constraint_System(object):
    """
    Wrapper for PPL's ``Constraint_System`` class.

    An object of the class Constraint_System is a system of
    constraints, i.e., a multiset of objects of the class
    Constraint. When inserting constraints in a system, space
    dimensions are automatically adjusted so that all the constraints
    in the system are defined on the same vector space.

    EXAMPLES::

        sage: from sage.libs.ppl import Constraint_System, Variable
        sage: x = Variable(0)
        sage: y = Variable(1)
        sage: cs = Constraint_System( 5*x-2*y > 0 )
        sage: cs.insert( 6*x<3*y )
        sage: cs.insert( x >= 2*x-7*y )
        sage: cs
        Constraint_System {5*x0-2*x1>0, -2*x0+x1>0, -x0+7*x1>=0}
    """

    cdef PPL_Constraint_System *thisptr


    def __cinit__(self, arg=None):
        """
        The Cython constructor.

        See :class:`Constraint_System` for documentation.

        TESTS::

            sage: from sage.libs.ppl import Constraint_System
            sage: Constraint_System()
            Constraint_System {}
        """
        if arg is None:
            self.thisptr = new PPL_Constraint_System()
            return
        if isinstance(arg, Constraint):
            g = <Constraint>arg
            self.thisptr = new PPL_Constraint_System(g.thisptr[0])
            return
        if isinstance(arg, Constraint_System):
            gs = <Constraint_System>arg
            self.thisptr = new PPL_Constraint_System(gs.thisptr[0])
            return
        if isinstance(arg, (list,tuple)):
            self.thisptr = new PPL_Constraint_System()
            for constraint in arg:
                self.insert(constraint)
            return
        raise ValueError, 'Cannot initialize with '+str(arg)+'.'


    def __dealloc__(self):
        """
        The Cython destructor.
        """
        del self.thisptr


    def space_dimension(self):
        r"""
        Return the dimension of the vector space enclosing ``self``.

        OUTPUT:

        Integer.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, Constraint_System
            sage: x = Variable(0)
            sage: cs = Constraint_System( x>0 )
            sage: cs.space_dimension()
            1
        """
        return self.thisptr.space_dimension()


    def has_equalities(self):
        r"""
        Tests whether ``self`` contains one or more equality constraints.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, Constraint_System
            sage: x = Variable(0)
            sage: cs = Constraint_System()
            sage: cs.insert( x>0 )
            sage: cs.insert( x<0 )
            sage: cs.has_equalities()
            False
            sage: cs.insert( x==0 )
            sage: cs.has_equalities()
            True
        """
        return self.thisptr.has_equalities()


    def has_strict_inequalities(self):
        r"""
        Tests whether ``self`` contains one or more strict inequality constraints.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, Constraint_System
            sage: x = Variable(0)
            sage: cs = Constraint_System()
            sage: cs.insert( x>=0 )
            sage: cs.insert( x==-1 )
            sage: cs.has_strict_inequalities()
            False
            sage: cs.insert( x>0 )
            sage: cs.has_strict_inequalities()
            True
        """
        return self.thisptr.has_strict_inequalities()


    def clear(self):
        r"""
        Removes all constraints from the constraint system and sets its
        space dimension to 0.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, Constraint_System
            sage: x = Variable(0)
            sage: cs = Constraint_System(x>0)
            sage: cs
            Constraint_System {x0>0}
            sage: cs.clear()
            sage: cs
            Constraint_System {}
        """
        self.assert_mutable('The Constraint_System is not mutable!')
        self.thisptr.clear()


    def insert(self, Constraint c):
        """
        Insert ``c`` into the constraint system.

        INPUT:

        - ``c`` -- a :class:`Constraint`.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, Constraint_System
            sage: x = Variable(0)
            sage: cs = Constraint_System()
            sage: cs.insert( x>0 )
            sage: cs
            Constraint_System {x0>0}
        """
        self.assert_mutable('The Constraint_System is not mutable!')
        self.thisptr.insert(c.thisptr[0])


    def empty(self):
        """
        Return ``True`` if and only if ``self`` has no constraints.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, Constraint_System, point
            sage: x = Variable(0)
            sage: cs = Constraint_System()
            sage: cs.empty()
            True
            sage: cs.insert( x>0 )
            sage: cs.empty()
            False
        """
        return self.thisptr.empty()


    def ascii_dump(self):
        r"""
        Write an ASCII dump to stderr.

        EXAMPLES::

            sage: sage_cmd  = 'from sage.libs.ppl import Constraint_System, Variable\n'
            sage: sage_cmd += 'x = Variable(0)\n'
            sage: sage_cmd += 'y = Variable(1)\n'
            sage: sage_cmd += 'cs = Constraint_System( 3*x > 2*y+1 )\n'
            sage: sage_cmd += 'cs.ascii_dump()\n'
            sage: from sage.tests.cmdline import test_executable
            sage: (out, err, ret) = test_executable(['sage', '-c', sage_cmd], timeout=100);  # long time, indirect doctest
            sage: print err  # long time
            topology NOT_NECESSARILY_CLOSED
            1 x 2 SPARSE (sorted)
            index_first_pending 1
            size 4 -1 3 -2 -1 > (NNC)
        """
        self.thisptr.ascii_dump()


    def OK(self):
        """
        Check if all the invariants are satisfied.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, Constraint_System
            sage: x = Variable(0)
            sage: y = Variable(1)
            sage: cs = Constraint_System( 3*x+2*y+1 <= 10 )
            sage: cs.OK()
            True
        """
        return self.thisptr.OK()


    def __len__(self):
        """
        Return the number of constraints in the system.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, Constraint_System
            sage: x = Variable(0)
            sage: cs = Constraint_System( x>0 )
            sage: cs.insert( x<1 )
            sage: len(cs)
            2
        """
        return sum([1 for c in self])


    def __iter__(self):
        """
        Iterate through the constraints of the system.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, Constraint_System
            sage: x = Variable(0)
            sage: cs = Constraint_System( x>0 )
            sage: iter = cs.__iter__()
            sage: iter.next()
            x0>0
            sage: list(cs)   # uses __iter__() internally
            [x0>0]
        """
        return Constraint_System_iterator(self)


    def __getitem__(self, int k):
        """
        Return the k-th constraint.

        The correct way to read the individual constraints is to
        iterate over the constraint system. This method is for
        convenience only.

        INPUT:

        - ``k`` -- integer. The index of the constraint.

        OUTPUT:

        The `k`-th constraint of the constraint system.

        EXAMPLES::

            sage: from sage.libs.ppl import Variable, Constraint_System
            sage: x = Variable(0)
            sage: cs = Constraint_System( x>0 )
            sage: cs.insert( x<1 )
            sage: cs
            Constraint_System {x0>0, -x0+1>0}
            sage: cs[0]
            x0>0
            sage: cs[1]
            -x0+1>0
        """
        if k < 0:
            raise IndexError('index must be nonnegative')
        iterator = self.__iter__()
        try:
            for i in range(k):
                iterator.next()
        except StopIteration:
            raise IndexError('index is past-the-end')
        return iterator.next()


    def __repr__(self):
        r"""
        Return a string representation of the constraint system.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: from sage.libs.ppl import Constraint_System, Variable
            sage: x = Variable(0)
            sage: y = Variable(1)
            sage: cs = Constraint_System([3*x+2*y+1 < 3, 0*x>x+1])
            sage: cs.__repr__()
            'Constraint_System {-3*x0-2*x1+2>0, -x0-1>0}'
        """
        s = 'Constraint_System {'
        s += ', '.join([ c.__repr__() for c in self ])
        s += '}'
        return s


    def __reduce__(self):
        """
        Pickle object.

        TESTS::

            sage: from sage.libs.ppl import Constraint_System, Variable
            sage: x = Variable(0)
            sage: y = Variable(1)
            sage: cs = Constraint_System([3*x+2*y+1 < 3, 0*x>x+1]);  cs
            Constraint_System {-3*x0-2*x1+2>0, -x0-1>0}
            sage: loads(dumps(cs))
            Constraint_System {-3*x0-2*x1+2>0, -x0-1>0}
        """
        return (Constraint_System, (tuple(self), ))


####################################################
### Constraint_System_iterator #####################
####################################################
cdef extern from "ppl_shim.hh":
    ctypedef void* cs_iterator_ptr
    cdef cs_iterator_ptr init_cs_iterator(PPL_Constraint_System &cs)
    cdef PPL_Constraint next_cs_iterator(cs_iterator_ptr)
    cdef bint is_end_cs_iterator(PPL_Constraint_System &cs, cs_iterator_ptr csi_ptr)
    cdef void delete_cs_iterator(cs_iterator_ptr)


####################################################
cdef class Constraint_System_iterator(object):
    """
    Wrapper for PPL's ``Constraint_System::const_iterator`` class.

    EXAMPLES::

        sage: from sage.libs.ppl import Constraint_System, Variable, Constraint_System_iterator
        sage: x = Variable(0)
        sage: y = Variable(1)
        sage: cs = Constraint_System( 5*x < 2*y )
        sage: cs.insert( 6*x-3*y==0 )
        sage: cs.insert( x >= 2*x-7*y )
        sage: Constraint_System_iterator(cs).next()
        -5*x0+2*x1>0
        sage: list(cs)
        [-5*x0+2*x1>0, 2*x0-x1==0, -x0+7*x1>=0]
    """

    cdef Constraint_System cs
    cdef cs_iterator_ptr csi_ptr


    def __cinit__(self, Constraint_System cs):
        """
        The Cython constructor.

        See :class:`Constraint_System_iterator` for documentation.

        TESTS::

            sage: from sage.libs.ppl import Constraint_System, Constraint_System_iterator
            sage: iter = Constraint_System_iterator( Constraint_System() )   # indirect doctest
        """
        self.cs = cs
        self.csi_ptr = init_cs_iterator(cs.thisptr[0])


    def __dealloc__(self):
        """
        The Cython destructor.
        """
        delete_cs_iterator(self.csi_ptr)


    def __next__(Constraint_System_iterator self):
        r"""
        The next iteration.

        OUTPUT:

        A :class:`Generator`.

        EXAMPLES::

            sage: from sage.libs.ppl import Constraint_System, Variable, Constraint_System_iterator
            sage: x = Variable(0)
            sage: cs = Constraint_System( 5*x > 0 )
            sage: Constraint_System_iterator(cs).next()
            x0>0
        """
        if is_end_cs_iterator((<Constraint_System>self.cs).thisptr[0], self.csi_ptr):
            raise StopIteration
        return _wrap_Constraint(next_cs_iterator(self.csi_ptr))



####################################################
### Poly_Gen_Relation ##############################
####################################################
cdef _wrap_Poly_Gen_Relation(PPL_Poly_Gen_Relation relation):
    """
    Wrap a C++ ``PPL_Poly_Gen_Relation`` into a Cython ``Poly_Gen_Relation``.
    """
    cdef Poly_Gen_Relation rel = Poly_Gen_Relation(True)
    rel.thisptr = new PPL_Poly_Gen_Relation(relation)
    return rel


####################################################
cdef class Poly_Gen_Relation(object):
    r"""
    Wrapper for PPL's ``Poly_Con_Relation`` class.

    INPUT/OUTPUT:

    You must not construct :class:`Poly_Gen_Relation` objects
    manually. You will usually get them from
    :meth:`~sage.libs.ppl.Polyhedron.relation_with`. You can also get
    pre-defined relations from the class methods :meth:`nothing` and
    :meth:`subsumes`.

    EXAMPLES::

        sage: from sage.libs.ppl import Poly_Gen_Relation
        sage: nothing = Poly_Gen_Relation.nothing(); nothing
        nothing
        sage: subsumes = Poly_Gen_Relation.subsumes(); subsumes
        subsumes
        sage: nothing.implies( subsumes )
        False
        sage: subsumes.implies( nothing )
        True
    """

    cdef PPL_Poly_Gen_Relation *thisptr


    def __cinit__(self, do_not_construct_manually=False):
        """
        The Cython constructor.

        See :class:`Poly_Gen_Relation` for documentation.

        TESTS::

            sage: from sage.libs.ppl import Poly_Gen_Relation
            sage: Poly_Gen_Relation.nothing()
            nothing
        """
        assert(do_not_construct_manually)
        self.thisptr = NULL


    def __dealloc__(self):
        """
        The Cython destructor.
        """
        assert self.thisptr!=NULL, 'Do not construct Poly_Gen_Relation objects manually!'
        del self.thisptr


    def implies(self, Poly_Gen_Relation y):
        r"""
        Test whether ``self`` implies ``y``.

        INPUT:

        - ``y`` -- a :class:`Poly_Gen_Relation`.

        OUTPUT:

        Boolean. ``True`` if and only if ``self`` implies ``y``.

        EXAMPLES::

            sage: from sage.libs.ppl import Poly_Gen_Relation
            sage: nothing = Poly_Gen_Relation.nothing()
            sage: nothing.implies( nothing )
            True
        """
        return self.thisptr.implies(y.thisptr[0])


    @classmethod
    def nothing(cls):
        r"""
        Return the assertion that says nothing.

        OUTPUT:

        A :class:`Poly_Gen_Relation`.

        EXAMPLES::

            sage: from sage.libs.ppl import Poly_Gen_Relation
            sage: Poly_Gen_Relation.nothing()
            nothing
        """
        return _wrap_Poly_Gen_Relation(PPL_Poly_Gen_Relation_nothing())


    @classmethod
    def subsumes(cls):
        r"""
        Return the assertion "Adding the generator would not change
        the polyhedron".

        OUTPUT:

        A :class:`Poly_Gen_Relation`.

        EXAMPLES::

            sage: from sage.libs.ppl import Poly_Gen_Relation
            sage: Poly_Gen_Relation.subsumes()
            subsumes
        """
        return _wrap_Poly_Gen_Relation(PPL_Poly_Gen_Relation_subsumes())


    def ascii_dump(self):
        r"""
        Write an ASCII dump to stderr.

        EXAMPLES::

            sage: sage_cmd  = 'from sage.libs.ppl import Poly_Gen_Relation\n'
            sage: sage_cmd += 'Poly_Gen_Relation.nothing().ascii_dump()\n'
            sage: from sage.tests.cmdline import test_executable
            sage: (out, err, ret) = test_executable(['sage', '-c', sage_cmd], timeout=100);  # long time, indirect doctest
            sage: print err  # long time
            NOTHING
        """
        self.thisptr.ascii_dump()


    def OK(self, check_non_empty=False):
        """
        Check if all the invariants are satisfied.

        EXAMPLES::

            sage: from sage.libs.ppl import Poly_Gen_Relation
            sage: Poly_Gen_Relation.nothing().OK()
            True
        """
        return self.thisptr.OK()


    def __repr__(self):
        r"""
        Return a string representation.

        OUTPUT:

        String.

        EXAMPLES::

            sage: from sage.libs.ppl import Poly_Gen_Relation
            sage: Poly_Gen_Relation.nothing().__repr__()
            'nothing'
        """
        if self.implies(Poly_Gen_Relation.subsumes()):
            return 'subsumes'
        else:
            return 'nothing'


####################################################
### Poly_Con_Relation ##############################
####################################################
cdef _wrap_Poly_Con_Relation(PPL_Poly_Con_Relation relation):
    """
    Wrap a C++ ``PPL_Poly_Con_Relation`` into a Cython ``Poly_Con_Relation``.
    """
    cdef Poly_Con_Relation rel = Poly_Con_Relation(True)
    rel.thisptr = new PPL_Poly_Con_Relation(relation)
    return rel


####################################################
cdef class Poly_Con_Relation(object):
    r"""
    Wrapper for PPL's ``Poly_Con_Relation`` class.

    INPUT/OUTPUT:

    You must not construct :class:`Poly_Con_Relation` objects
    manually. You will usually get them from
    :meth:`~sage.libs.ppl.Polyhedron.relation_with`. You can also get
    pre-defined relations from the class methods :meth:`nothing`,
    :meth:`is_disjoint`, :meth:`strictly_intersects`,
    :meth:`is_included`, and :meth:`saturates`.

    EXAMPLES::

        sage: from sage.libs.ppl import Poly_Con_Relation
        sage: saturates     = Poly_Con_Relation.saturates();  saturates
        saturates
        sage: is_included   = Poly_Con_Relation.is_included(); is_included
        is_included
        sage: is_included.implies(saturates)
        False
        sage: saturates.implies(is_included)
        False
        sage: rels = []
        sage: rels.append( Poly_Con_Relation.nothing() )
        sage: rels.append( Poly_Con_Relation.is_disjoint() )
        sage: rels.append( Poly_Con_Relation.strictly_intersects() )
        sage: rels.append( Poly_Con_Relation.is_included() )
        sage: rels.append( Poly_Con_Relation.saturates() )
        sage: rels
        [nothing, is_disjoint, strictly_intersects, is_included, saturates]
        sage: from sage.matrix.constructor import matrix
        sage: m = matrix(5,5)
        sage: for i, rel_i in enumerate(rels):
        ...       for j, rel_j in enumerate(rels):
        ...           m[i,j] = rel_i.implies(rel_j)
        sage: m
        [1 0 0 0 0]
        [1 1 0 0 0]
        [1 0 1 0 0]
        [1 0 0 1 0]
        [1 0 0 0 1]
    """

    cdef PPL_Poly_Con_Relation *thisptr


    def __cinit__(self, do_not_construct_manually=False):
        """
        The Cython constructor.

        See :class:`Poly_Con_Relation` for documentation.

        TESTS::

            sage: from sage.libs.ppl import Poly_Con_Relation
            sage: Poly_Con_Relation.nothing()
            nothing
        """
        assert(do_not_construct_manually)
        self.thisptr = NULL


    def __dealloc__(self):
        """
        The Cython destructor.
        """
        assert self.thisptr!=NULL, 'Do not construct Poly_Con_Relation objects manually!'
        del self.thisptr


    def implies(self, Poly_Con_Relation y):
        r"""
        Test whether ``self`` implies ``y``.

        INPUT:

        - ``y`` -- a :class:`Poly_Con_Relation`.

        OUTPUT:

        Boolean. ``True`` if and only if ``self`` implies ``y``.

        EXAMPLES::

            sage: from sage.libs.ppl import Poly_Con_Relation
            sage: nothing = Poly_Con_Relation.nothing()
            sage: nothing.implies( nothing )
            True
        """
        return self.thisptr.implies(y.thisptr[0])


    @classmethod
    def nothing(cls):
        r"""
        Return the assertion that says nothing.

        OUTPUT:

        A :class:`Poly_Con_Relation`.

        EXAMPLES::

            sage: from sage.libs.ppl import Poly_Con_Relation
            sage: Poly_Con_Relation.nothing()
            nothing
        """
        return _wrap_Poly_Con_Relation(PPL_Poly_Con_Relation_nothing())


    @classmethod
    def is_disjoint(cls):
        r"""
        Return the assertion "The polyhedron and the set of points
        satisfying the constraint are disjoint".

        OUTPUT:

        A :class:`Poly_Con_Relation`.

        EXAMPLES::

            sage: from sage.libs.ppl import Poly_Con_Relation
            sage: Poly_Con_Relation.is_disjoint()
            is_disjoint
        """
        return _wrap_Poly_Con_Relation(PPL_Poly_Con_Relation_is_disjoint())


    @classmethod
    def strictly_intersects(cls):
        r"""
        Return the assertion "The polyhedron intersects the set of
        points satisfying the constraint, but it is not included in
        it".

        OUTPUT:

        A :class:`Poly_Con_Relation`.

        EXAMPLES::

            sage: from sage.libs.ppl import Poly_Con_Relation
            sage: Poly_Con_Relation.strictly_intersects()
            strictly_intersects
        """
        return _wrap_Poly_Con_Relation(PPL_Poly_Con_Relation_strictly_intersects())


    @classmethod
    def is_included(cls):
        r"""
        Return the assertion "The polyhedron is included in the set of
        points satisfying the constraint".

        OUTPUT:

        A :class:`Poly_Con_Relation`.

        EXAMPLES::

            sage: from sage.libs.ppl import Poly_Con_Relation
            sage: Poly_Con_Relation.is_included()
            is_included
        """
        return _wrap_Poly_Con_Relation(PPL_Poly_Con_Relation_is_included())


    @classmethod
    def saturates(cls):
        r"""
        Return the assertion "".

        OUTPUT:

        A :class:`Poly_Con_Relation`.

        EXAMPLES::

            sage: from sage.libs.ppl import Poly_Con_Relation
            sage: Poly_Con_Relation.saturates()
            saturates
        """
        return _wrap_Poly_Con_Relation(PPL_Poly_Con_Relation_saturates())


    def ascii_dump(self):
        r"""
        Write an ASCII dump to stderr.

        EXAMPLES::

            sage: sage_cmd  = 'from sage.libs.ppl import Poly_Con_Relation\n'
            sage: sage_cmd += 'Poly_Con_Relation.nothing().ascii_dump()\n'
            sage: from sage.tests.cmdline import test_executable
            sage: (out, err, ret) = test_executable(['sage', '-c', sage_cmd], timeout=100);  # long time, indirect doctest
            sage: print err  # long time
            NOTHING
        """
        self.thisptr.ascii_dump()


    def OK(self, check_non_empty=False):
        """
        Check if all the invariants are satisfied.

        EXAMPLES::

            sage: from sage.libs.ppl import Poly_Con_Relation
            sage: Poly_Con_Relation.nothing().OK()
            True
        """
        return self.thisptr.OK()


    def __repr__(self):
        r"""
        Return a string representation.

        OUTPUT:

        String.

        EXAMPLES::

            sage: from sage.libs.ppl import Poly_Con_Relation
            sage: Poly_Con_Relation.nothing().__repr__()
            'nothing'
        """
        rel = []
        if self.implies(Poly_Con_Relation.is_disjoint()):
            rel.append('is_disjoint')
        if self.implies(Poly_Con_Relation.strictly_intersects()):
            rel.append('strictly_intersects')
        if self.implies(Poly_Con_Relation.is_included()):
            rel.append('is_included')
        if self.implies(Poly_Con_Relation.saturates()):
            rel.append('saturates')

        if len(rel)>0:
            return ', '.join(rel)
        else:
            return 'nothing'


####################################################
####################################################
####################################################
####################################################
