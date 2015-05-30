r"""
Interactive Simplex Method

This module, meant for **educational purposes only**, supports learning and
exploring of the simplex method.

Do you want to solve Linear Programs efficiently? use
:class:`MixedIntegerLinearProgram` instead.

The methods implemented here allow solving Linear Programming Problems (LPPs) in
a number of ways, may require explicit (and correct!) description of steps and
are likely to be much slower than "regular" LP solvers. If, however, you want to
learn how the simplex method works and see what happens in different situations
using different strategies, but don't want to deal with tedious arithmetic, this
module is for you!

Historically it was created to complement the Math 373 course on Mathematical
Programming and Optimization at the University of Alberta, Edmonton, Canada.

AUTHORS:

- Andrey Novoseltsev (2013-03-16): initial version.

EXAMPLES:

Most of the module functionality is demonstrated on the following problem.

.. admonition:: Corn & Barley

    A farmer has 1000 acres available to grow corn and barley.
    Corn has a net profit of 10 dollars per acre while barley has a net
    profit of 5 dollars per acre.
    The farmer has 1500 kg of fertilizer available with 3 kg per acre
    needed for corn and 1 kg per acre needed for barley.
    The farmer wants to maximize profit.
    (Sometimes we also add one more constraint to make the initial dictionary
    infeasible: the farmer has to use at least 40% of the available land.)

Using variables `C` and `B` for land used to grow corn and barley respectively,
in acres, we can construct the following LP problem::

    sage: A = ([1, 1], [3, 1])
    sage: b = (1000, 1500)
    sage: c = (10, 5)
    sage: P = InteractiveLPProblem(A, b, c, ["C", "B"], variable_type=">=")
    sage: P
    LP problem (use typeset mode to see details)

It is recommended to copy-paste such examples into your own worksheet, so that
you can run these commands with typeset mode on and get

.. MATH::

    \begin{array}{l}
    \begin{array}{lcrcrcl}
     \max \mspace{-6mu}&\mspace{-6mu}  \mspace{-6mu}&\mspace{-6mu} 10 C \mspace{-6mu}&\mspace{-6mu} + \mspace{-6mu}&\mspace{-6mu} 5 B \mspace{-6mu} \\
     \mspace{-6mu}&\mspace{-6mu}  \mspace{-6mu}&\mspace{-6mu} C \mspace{-6mu}&\mspace{-6mu} + \mspace{-6mu}&\mspace{-6mu} B \mspace{-6mu}&\mspace{-6mu} \leq \mspace{-6mu}&\mspace{-6mu} 1000 \\
     \mspace{-6mu}&\mspace{-6mu}  \mspace{-6mu}&\mspace{-6mu} 3 C \mspace{-6mu}&\mspace{-6mu} + \mspace{-6mu}&\mspace{-6mu} B \mspace{-6mu}&\mspace{-6mu} \leq \mspace{-6mu}&\mspace{-6mu} 1500 \\
    \end{array} \\
    C, B \geq 0
    \end{array}

Since it has only two variables, we can solve it graphically::

    sage: P.plot()
    Graphics object consisting of 19 graphics primitives


The simplex method can be applied only to :class:`problems in standard form
<InteractiveLPProblemStandardForm>`, which can be created either directly ::

    sage: InteractiveLPProblemStandardForm(A, b, c, ["C", "B"])
    LP problem (use typeset mode to see details)

or from an already constructed problem of "general type"::

    sage: P = P.standard_form()

In this case the problem does not require any modifications to be written in
standard form, but this step is still necessary to enable methods related to
the simplex method.

The simplest way to use the simplex method is::

    sage: P.run_simplex_method()
    %notruncate
    ...

(This method produces quite long formulas which have been omitted here.)
But, of course, it is much more fun to do most of the steps by hand. Let's start
by creating the initial dictionary::

    sage: D = P.initial_dictionary()
    sage: D
    LP problem dictionary (use typeset mode to see details)

Using typeset mode as recommended, you'll see

.. MATH::

    \renewcommand{\arraystretch}{1.5}
    \begin{array}{|rcrcrcr|}
    \hline
    x_{3} \mspace{-6mu}&\mspace{-6mu} = \mspace{-6mu}&\mspace{-6mu} 1000 \mspace{-6mu}&\mspace{-6mu} - \mspace{-6mu}&\mspace{-6mu} C \mspace{-6mu}&\mspace{-6mu} - \mspace{-6mu}&\mspace{-6mu} B\\
    x_{4} \mspace{-6mu}&\mspace{-6mu} = \mspace{-6mu}&\mspace{-6mu} 1500 \mspace{-6mu}&\mspace{-6mu} - \mspace{-6mu}&\mspace{-6mu} 3 C \mspace{-6mu}&\mspace{-6mu} - \mspace{-6mu}&\mspace{-6mu} B\\
    \hline
    z \mspace{-6mu}&\mspace{-6mu} = \mspace{-6mu}&\mspace{-6mu} 0 \mspace{-6mu}&\mspace{-6mu} + \mspace{-6mu}&\mspace{-6mu} 10 C \mspace{-6mu}&\mspace{-6mu} + \mspace{-6mu}&\mspace{-6mu} 5 B\\
    \hline
    \end{array}

With the initial or any other dictionary you can perform a number of checks::

    sage: D.is_feasible()
    True
    sage: D.is_optimal()
    False

You can look at many of its pieces and associated data::

    sage: D.basic_variables()
    (x3, x4)
    sage: D.basic_solution()
    (0, 0)
    sage: D.objective_value()
    0

Most importantly, you can perform steps of the simplex method by picking an
entering variable, a leaving variable, and updating the dictionary::

    sage: D.enter("C")
    sage: D.leave(4)
    sage: D.update()

If everything was done correctly, the new dictionary is still feasible and the
objective value did not decrease::

    sage: D.is_feasible()
    True
    sage: D.objective_value()
    5000

If you are unsure about picking entering and leaving variables, you can use
helper methods that will try their best to tell you what are your next
options::

    sage: D.possible_entering()
    [B]
    sage: D.possible_leaving()
    Traceback (most recent call last):
    ...
    ValueError: leaving variables can be determined
     for feasible dictionaries with a set entering variable
     or for dual feasible dictionaries

It is also possible to obtain :meth:`feasible sets <InteractiveLPProblem.feasible_set>`
and :meth:`final dictionaries <InteractiveLPProblemStandardForm.final_dictionary>` of
problems, work with :class:`revised dictionaries <LPRevisedDictionary>`,
and use the dual simplex method!

.. NOTE::

    Currently this does not have a display format for the terminal.

Classes and functions
---------------------
"""


#*****************************************************************************
#       Copyright (C) 2013 Andrey Novoseltsev <novoselt@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


import operator, re


from copy import copy


from sage.geometry.all import Polyhedron
from sage.matrix.all import (column_matrix,
                             identity_matrix,
                             matrix,
                             random_matrix)
from sage.misc.all import (LatexExpr,
                           cached_function,
                           cached_method,
                           latex,
                           randint,
                           random)
from sage.misc.misc import get_main_globals
from sage.modules.all import random_vector, vector
from sage.plot.all import Graphics, arrow, line, point, rainbow, text
from sage.rings.all import Infinity, PolynomialRing, QQ, RDF, ZZ
from sage.structure.all import SageObject
from sage.symbolic.all import SR


# We produce rather complicated LaTeX code which needs some tweaks to be
# displayed nicely by MathJax, which make it look slightly different from real
# LaTeX. We use our own variable as it may be convenient to override it.
# Hopefully, some day there will be no need in it at all and only "if" parts
# will have to be left.
generate_real_LaTeX = False

def _assemble_arrayl(lines, stretch=None):
    r"""
    Return ``lines`` assembled in a left-justified array.

    INPUT:

    - ``lines`` -- a list of strings suitable for math mode typesetting

    - ``stretch`` -- (default: None) if given, a command setting
      ``\arraystretch`` to this value will be added before the array

    OUTPUT:

    - a :class:`LatexExpr`

    EXAMPLES::

        sage: from sage.numerical.interactive_simplex_method \
        ....:     import _assemble_arrayl
        sage: lines = ["1 + 1", "2"]
        sage: print _assemble_arrayl(lines)
        %notruncate
        \begin{array}{l}
        1 + 1\\
        2
        \end{array}
        sage: print _assemble_arrayl(lines, 1.5)
        %notruncate
        \renewcommand{\arraystretch}{1.500000}
        \begin{array}{l}
        1 + 1\\
        2
        \end{array}
    """
    # Even simple LP problems tend to generate long output, so we prohibit
    # truncation in the notebook cells and hope for the best!
    return LatexExpr(("" if generate_real_LaTeX else "%notruncate\n") +
                     ("" if stretch is None else
                      "\\renewcommand{\\arraystretch}{%f}\n" % stretch) +
                     "\\begin{array}{l}\n" +
                     "\\\\\n".join(lines) +
                     "\n\\end{array}")


def _latex_product(coefficients, variables,
                   separator=None, head=None, tail=None,
                   drop_plus=True, allow_empty=False):
    r"""
    Generate LaTeX code for a linear function.

    This function is intended for internal use by LaTeX methods of LP problems
    and their dictionaries.

    INPUT:

    - ``coefficients`` -- a list of coefficients

    - ``variables`` -- a list of variables

    - ``separator`` -- (default: ``"&"`` with some extra space adjustment) a
      string to be inserted between elements of the generated expression

    - ``head`` -- either ``None`` (default) or a list of entries to be
      added to the beginning of the output

    - ``tail`` -- either ``None`` (default) or a list of entries to be
      added to the end of the output

    - ``drop_plus`` -- (default: ``True``) whether to drop the leading plus
      sign or not

    - ``allow_empty`` -- (default: ``False``) whether to allow empty output or
      produce at least "0"

    OUTPUT:

    - A string joining ``head``, each sign and coefficient-variable product,
      and ``tail`` using ``separator``. Strings in ``head`` and ``tail`` are
      used as is except for "<=", "==", and ">=", which are replaced by LaTeX
      commands. Other elements in ``head`` in ``tail`` are processed by
      :func:`latex`.

    TESTS::

        sage: from sage.numerical.interactive_simplex_method import \
        ....:       _latex_product
        sage: var("x, y")
        (x, y)
        sage: print _latex_product([-1, 3], [x, y])
        - \mspace{-6mu}&\mspace{-6mu} x \mspace{-6mu}&\mspace{-6mu} + \mspace{-6mu}&\mspace{-6mu} 3 y
    """
    entries = []
    for c, v in zip(coefficients, variables):
        if c == 0:
            entries.extend(["", ""])
            continue
        sign = "+"
        if latex(c).strip().startswith("-"):
            sign = "-"
            c = - c
        if c == 1:
            t = latex(v)
        else:
            t = latex(c)
            if SR(c).operator() in [operator.add, operator.sub]:
                t = r"\left( " + t + r" \right)"
            t += " " + latex(v)
        entries.extend([sign, t])
    if drop_plus:   # Don't start with +
        for i, e in enumerate(entries):
            if e:   # The first non-empty
                if e == "+":
                    entries[i] = ""
                break
    if not (allow_empty or any(entries)):   # Return at least 0
        entries[-1] = "0"
    latex_relations = {"<=": r"\leq", "==": "=", ">=": r"\geq"}
    if head is not None:
        for e in reversed(head):
            if not isinstance(e, str):
                e = latex(e)
            elif e in latex_relations:
                e = latex_relations[e]
            entries.insert(0, e)
    if tail is not None:
        for e in tail:
            if not isinstance(e, str):
                e = latex(e)
            elif e in latex_relations:
                e = latex_relations[e]
            entries.append(e)
    if separator is None:
        if generate_real_LaTeX:
            separator = " & "
        else:
            separator = r" \mspace{-6mu}&\mspace{-6mu} "
    return separator.join(entries)


@cached_function
def variable(R, v):
    r"""
    Interpret ``v`` as a variable of ``R``.

    INPUT:

    - ``R`` -- a polynomial ring

    - ``v`` -- a variable of ``R`` or convertible into ``R``, a string
      with the name of a variable of ``R`` or an index of a variable in ``R``

    OUTPUT:

    - a variable of ``R``

    EXAMPLES::

        sage: from sage.numerical.interactive_simplex_method \
        ....:     import variable
        sage: R = PolynomialRing(QQ, "x3, y5, x5, y")
        sage: R.inject_variables()
        Defining x3, y5, x5, y
        sage: variable(R, "x3")
        x3
        sage: variable(R, x3)
        x3
        sage: variable(R, 3)
        x3
        sage: variable(R, 0)
        Traceback (most recent call last):
        ...
        ValueError: there is no variable with the given index
        sage: variable(R, 5)
        Traceback (most recent call last):
        ...
        ValueError: the given index is ambiguous
        sage: variable(R, 2 * x3)
        Traceback (most recent call last):
        ...
        ValueError: cannot interpret given data as a variable
        sage: variable(R, "z")
        Traceback (most recent call last):
        ...
        ValueError: cannot interpret given data as a variable
    """
    if v in ZZ:
        v = str(v)
        tail = re.compile(r"\d+$")
        matches = []
        for g in R.gens():
            match = tail.search(str(g))
            if match is not None and match.group() == v:
                matches.append(g)
        if not matches:
            raise ValueError("there is no variable with the given index")
        if len(matches) > 1:
            raise ValueError("the given index is ambiguous")
        return matches[0]
    else:
        try:
            v = R(v)
            if v in R.gens():
                return v
        except TypeError:
            pass
        raise ValueError("cannot interpret given data as a variable")


class InteractiveLPProblem(SageObject):
    r"""
    Construct an LP (Linear Programming) problem.

    .. NOTE::

        This class is for **educational purposes only**: if you want to solve
        Linear Programs efficiently, use :class:`MixedIntegerLinearProgram`
        instead.

    This class supports LP problems with "variables on the left" constraints.

    INPUT:

    - ``A`` -- a matrix of constraint coefficients

    - ``b`` -- a vector of constraint constant terms

    - ``c`` -- a vector of objective coefficients

    - ``x`` -- (default: ``"x"``) a vector of decision variables or a
      string giving the base name

    - ``constraint_type`` -- (default: ``"<="``) a string specifying constraint
      type(s): either ``"<="``, ``">="``, ``"=="``, or a list of them

    - ``variable_type`` -- (default: ``""``) a string specifying variable
      type(s): either ``">="``, ``"<="``, ``""`` (the empty string), or a
      list of them, corresponding, respectively, to non-negative,
      non-positive, and free variables

    - ``problem_type`` -- (default: ``"max"``) a string specifying the
      problem type: ``"max"``, ``"min"``, ``"-max"``, or ``"-min"``

    - ``prefix`` -- (default: parameter ``x`` if it was given as a string,
      string ``"x"`` otherwise) a string giving the base name of automatically
      created variables in :meth:`standard_form`

    - ``base_ring`` -- (default: the fraction field of a common ring for all
      input coefficients) a field to which all input coefficients will be
      converted

    EXAMPLES:

    We will construct the following problem:

    .. MATH::

        \begin{array}{l}
        \begin{array}{lcrcrcl}
         \max \mspace{-6mu}&\mspace{-6mu}  \mspace{-6mu}&\mspace{-6mu} 10 C \mspace{-6mu}&\mspace{-6mu} + \mspace{-6mu}&\mspace{-6mu} 5 B \mspace{-6mu} \\
         \mspace{-6mu}&\mspace{-6mu}  \mspace{-6mu}&\mspace{-6mu} C \mspace{-6mu}&\mspace{-6mu} + \mspace{-6mu}&\mspace{-6mu} B \mspace{-6mu}&\mspace{-6mu} \leq \mspace{-6mu}&\mspace{-6mu} 1000 \\
         \mspace{-6mu}&\mspace{-6mu}  \mspace{-6mu}&\mspace{-6mu} 3 C \mspace{-6mu}&\mspace{-6mu} + \mspace{-6mu}&\mspace{-6mu} B \mspace{-6mu}&\mspace{-6mu} \leq \mspace{-6mu}&\mspace{-6mu} 1500 \\
        \end{array} \\
        C, B \geq 0
        \end{array}

    ::

        sage: A = ([1, 1], [3, 1])
        sage: b = (1000, 1500)
        sage: c = (10, 5)
        sage: P = InteractiveLPProblem(A, b, c, ["C", "B"], variable_type=">=")

    Same problem, but more explicitly::

        sage: P = InteractiveLPProblem(A, b, c, ["C", "B"],
        ....:     constraint_type="<=", variable_type=">=")

    Even more explicitly::

        sage: P = InteractiveLPProblem(A, b, c, ["C", "B"], problem_type="max",
        ....:     constraint_type=["<=", "<="], variable_type=[">=", ">="])

    Using the last form you should be able to represent any LP problem, as long
    as all like terms are collected and in constraints variables and constants
    are on different sides.
    """

    def __init__(self, A, b, c, x="x",
                 constraint_type="<=", variable_type="", problem_type="max",
                 prefix="x", base_ring=None):
        r"""
        See :class:`InteractiveLPProblem` for documentation.

        TESTS::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblem(A, b, c, ["C", "B"], variable_type=">=")
            sage: TestSuite(P).run()
        """
        super(InteractiveLPProblem, self).__init__()
        A = matrix(A)
        b = vector(b)
        c = vector(c)
        if base_ring is None:
            base_ring = vector(A.list() + list(b) + list(c)).base_ring()
        base_ring = base_ring.fraction_field()
        A = A.change_ring(base_ring)
        A.set_immutable()
        m, n = A.nrows(), A.ncols()
        b = b.change_ring(base_ring)
        b.set_immutable()
        if b.degree() != m:
            raise ValueError("A and b have incompatible dimensions")
        c = c.change_ring(base_ring)
        c.set_immutable()
        if c.degree() != n:
            raise ValueError("A and c have incompatible dimensions")
        if isinstance(x, str):
            prefix = x
            x = ["{}{:d}".format(x, i) for i in range(1, n+1)]
        else:
            x = map(str, x)
            if len(x) != n:
                raise ValueError("A and x have incompatible dimensions")
        R = PolynomialRing(base_ring, x, order="neglex")
        x = vector(R, R.gens()) # All variables as a vector
        self._Abcx = A, b, c, x

        if constraint_type in ["<=", ">=", "=="]:
            constraint_type = (constraint_type, ) * m
        else:
            constraint_type = tuple(constraint_type)
            if any(ct not in ["<=", ">=", "=="] for ct in constraint_type):
                raise ValueError("unknown constraint type")
            if len(constraint_type) != m:
                raise ValueError("wrong number of constraint types")
        self._constraint_types = constraint_type

        if variable_type in ["<=", ">=", ""]:
            variable_type = (variable_type, ) * n
        else:
            variable_type = tuple(variable_type)
            if any(vt not in ["<=", ">=", ""] for vt in variable_type):
                raise ValueError("unknown variable type")
            if len(variable_type) != n:
                raise ValueError("wrong number of variable types")
        self._variable_types = variable_type

        if problem_type.startswith("-"):
            self._is_negative = True
            problem_type = problem_type[1:]
        else:
            self._is_negative = False
        if problem_type not in ["max", "min"]:
            raise ValueError("unknown problem type")
        self._problem_type = problem_type

        self._prefix = prefix

    def __eq__(self, other):
        r"""
        Check if two LP problems are equal.

        INPUT:

        - ``other`` -- anything

        OUTPUT:

        - ``True`` if ``other`` is an :class:`InteractiveLPProblem` with all details the
          same as ``self``, ``False`` otherwise.

        TESTS::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblem(A, b, c, ["C", "B"], variable_type=">=")
            sage: P2 = InteractiveLPProblem(A, b, c, ["C", "B"], variable_type=">=")
            sage: P == P2
            True
            sage: P3 = InteractiveLPProblem(A, c, b, ["C", "B"], variable_type=">=")
            sage: P == P3
            False
        """
        return (isinstance(other, InteractiveLPProblem) and
                self.Abcx() == other.Abcx() and
                self._problem_type == other._problem_type and
                self._is_negative == other._is_negative and
                self._constraint_types == other._constraint_types and
                self._variable_types == other._variable_types and
                self._prefix == other._prefix)

    def _latex_(self):
        r"""
        Return a LaTeX representation of ``self``.

        OUTPUT:

        - a string

        TESTS::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblem(A, b, c, ["C", "B"], variable_type=">=")
            sage: print P._latex_()
            \begin{array}{l}
            \begin{array}{lcrcrcl}
             \max \mspace{-6mu}&\mspace{-6mu}  \mspace{-6mu}&\mspace{-6mu} 10 C \mspace{-6mu}&\mspace{-6mu} + \mspace{-6mu}&\mspace{-6mu} 5 B \mspace{-6mu} \\
             \mspace{-6mu}&\mspace{-6mu}  \mspace{-6mu}&\mspace{-6mu} C \mspace{-6mu}&\mspace{-6mu} + \mspace{-6mu}&\mspace{-6mu} B \mspace{-6mu}&\mspace{-6mu} \leq \mspace{-6mu}&\mspace{-6mu} 1000 \\
             \mspace{-6mu}&\mspace{-6mu}  \mspace{-6mu}&\mspace{-6mu} 3 C \mspace{-6mu}&\mspace{-6mu} + \mspace{-6mu}&\mspace{-6mu} B \mspace{-6mu}&\mspace{-6mu} \leq \mspace{-6mu}&\mspace{-6mu} 1500 \\
            \end{array} \\
            C, B \geq 0
            \end{array}
        """
        A, b, c, x = self.Abcx()
        lines = []
        lines.append(r"\begin{array}{l}")
        if generate_real_LaTeX:
            lines[-1] += r" \setlength{\arraycolsep}{0.125em}"
        lines.append(r"\begin{array}{l" + "cr" * len(x) + "cl}")
        head = [r"{} \{}".format("- " if self._is_negative else "",
                                 self._problem_type)]
        lines.append(_latex_product(c, x, head=head) +
                     (r"\\" if generate_real_LaTeX else r" \mspace{-6mu} \\"))
        for Ai, ri, bi in zip(A.rows(), self._constraint_types, b):
            lines.append(_latex_product(Ai, x, head=[""], tail=[ri, bi]) +
                         r" \\")
        lines.append(r"\end{array} \\")
        if set(self._variable_types) == set([">="]):
            lines.append(r"{} \geq 0".format(", ".join(map(latex, x))))
        else:
            lines.append(r",\ ".join(r"{} {} 0".format(
                                latex(xj), r"\geq" if vt == ">=" else r"\leq")
                            for xj, vt in zip(x, self._variable_types) if vt))
        lines.append(r"\end{array}")
        return  "\n".join(lines)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        OUTPUT:

        - a string

        TESTS::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblem(A, b, c, ["C", "B"], variable_type=">=")
            sage: print P._repr_()
            LP problem (use typeset mode to see details)
        """
        return "LP problem (use typeset mode to see details)"

    @cached_method
    def _solve(self):
        r"""
        Return an optimal solution and the optimal value of ``self``.

        OUTPUT:

        - A pair consisting of a vector and a number. If the problem is
          infeasible, both components are ``None``. If the problem is
          unbounded, the first component is ``None`` and the second is
          `\pm \infty`.

        This function uses "brute force" solution technique of evaluating the
        objective at all vertices of the feasible set and taking into account
        its rays and lines.

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblem(A, b, c, ["C", "B"], variable_type=">=")
            sage: P._solve()
            ((250, 750), 6250)
        """
        F = self.feasible_set()
        R = self.base_ring()
        A, b, c, x = self._Abcx
        if F.n_vertices() == 0:
            return (None, None)
        elif c.is_zero():
            M, S = 0, F.vertices()[0]
        elif self._problem_type == "max":
            if any(c * vector(R, ray) > 0 for ray in F.rays()) or \
               any(c * vector(R, line) != 0 for line in F.lines()):
                M, S = Infinity, None
            else:
                M, S = max((c * vector(R, v), v) for v in F.vertices())
        elif self._problem_type == "min":
            if any(c * vector(R, ray) < 0 for ray in F.rays()) or \
               any(c * vector(R, line) != 0 for line in F.lines()):
                M, S = -Infinity, None
            else:
                M, S = min((c * vector(R, v), v) for v in F.vertices())
        if self._is_negative:
            M = - M
        if S is not None:
            S = vector(R, S)
            S.set_immutable()
        return S, M

    def Abcx(self):
        r"""
        Return `A`, `b`, `c`, and `x` of ``self`` as a tuple.

        OUTPUT:

        - a tuple

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblem(A, b, c, ["C", "B"], variable_type=">=")
            sage: P.Abcx()
            (
            [1 1]
            [3 1], (1000, 1500), (10, 5), (C, B)
            )
        """
        return self._Abcx

    def base_ring(self):
        r"""
        Return the base ring of ``self``.

        .. NOTE::

            The base ring of LP problems is always a field.

        OUTPUT:

        - a ring

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblem(A, b, c, ["C", "B"], variable_type=">=")
            sage: P.base_ring()
            Rational Field

            sage: c = (10, 5.)
            sage: P = InteractiveLPProblem(A, b, c, ["C", "B"], variable_type=">=")
            sage: P.base_ring()
            Real Field with 53 bits of precision
        """
        return self._Abcx[0].base_ring()

    def constant_terms(self):
        r"""
        Return constant terms of constraints of ``self``, i.e. `b`.

        OUTPUT:

        - a vector

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblem(A, b, c, ["C", "B"], variable_type=">=")
            sage: P.constant_terms()
            (1000, 1500)
            sage: P.b()
            (1000, 1500)
        """
        return self._Abcx[1]

    def constraint_coefficients(self):
        r"""
        Return coefficients of constraints of ``self``, i.e. `A`.

        OUTPUT:

        - a matrix

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblem(A, b, c, ["C", "B"], variable_type=">=")
            sage: P.constraint_coefficients()
            [1 1]
            [3 1]
            sage: P.A()
            [1 1]
            [3 1]
        """
        return self._Abcx[0]

    def decision_variables(self):
        r"""
        Return decision variables of ``self``, i.e. `x`.

        OUTPUT:

        - a vector

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblem(A, b, c, ["C", "B"], variable_type=">=")
            sage: P.decision_variables()
            (C, B)
            sage: P.x()
            (C, B)
        """
        return self._Abcx[3]

    def dual(self, y=None):
        r"""
        Construct the dual LP problem for ``self``.

        INPUT:

        - ``y`` -- (default: ``"x"`` if the prefix of ``self`` is ``"y"``,
          ``"y"`` otherwise) a vector of dual decision variables or a string
          giving the base name

        OUTPUT:

        - an :class:`InteractiveLPProblem`

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblem(A, b, c, ["C", "B"], variable_type=">=")
            sage: DP = P.dual()
            sage: DP.b() == P.c()
            True
            sage: DP.dual(["C", "B"]) == P
            True
        """
        A, c, b, x = self.Abcx()
        A = A.transpose()
        if y is None:
            y = "x" if self._prefix == "y" else "y"
        problem_type = "min" if self._problem_type == "max" else "max"
        constraint_type = []
        for vt in self._variable_types:
            if (vt == ">=" and problem_type == "min" or
                vt == "<=" and problem_type == "max"):
                constraint_type.append(">=")
            elif (vt == "<=" and problem_type == "min" or
                vt == ">=" and problem_type == "max"):
                constraint_type.append("<=")
            else:
                constraint_type.append("==")
        variable_type = []
        for ct in self._constraint_types:
            if (ct == ">=" and problem_type == "min" or
                ct == "<=" and problem_type == "max"):
                variable_type.append("<=")
            elif (ct == "<=" and problem_type == "min" or
                ct == ">=" and problem_type == "max"):
                variable_type.append(">=")
            else:
                variable_type.append("")
        if self._is_negative:
            problem_type = "-" + problem_type
        return InteractiveLPProblem(A, b, c, y,
                         constraint_type, variable_type, problem_type)

    @cached_method
    def feasible_set(self):
        r"""
        Return the feasible set of ``self``.

        OUTPUT:

        - a :mod:`Polyhedron <sage.geometry.polyhedron.constructor>`

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblem(A, b, c, ["C", "B"], variable_type=">=")
            sage: P.feasible_set()
            A 2-dimensional polyhedron in QQ^2
            defined as the convex hull of 4 vertices
        """
        ieqs = []
        eqns = []
        for a, r, b in zip(self.A().rows(), self._constraint_types, self.b()):
            if r == "<=":
                ieqs.append([b] + list(-a))
            elif r == ">=":
                ieqs.append([-b] + list(a))
            else:
                eqns.append([-b] + list(a))
        for n, r in zip(identity_matrix(self.n()).rows(), self._variable_types):
            if r == "<=":
                ieqs.append([0] + list(-n))
            elif r == ">=":
                ieqs.append([0] + list(n))
        if self.base_ring() is QQ:
            R = QQ
        else:
            R = RDF
            ieqs = [map(R, ieq) for ieq in ieqs]
            eqns = [map(R, eqn) for eqn in eqns]
        return Polyhedron(ieqs=ieqs, eqns=eqns, base_ring=R)

    def is_bounded(self):
        r"""
        Check if ``self`` is bounded.

        OUTPUT:

        - ``True`` is ``self`` is bounded, ``False`` otherwise

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblem(A, b, c, ["C", "B"], variable_type=">=")
            sage: P.is_bounded()
            True
        """
        return self._solve()[0] is not None

    def is_feasible(self):
        r"""
        Check if ``self`` is feasible.

        OUTPUT:

        - ``True`` is ``self`` is feasible, ``False`` otherwise

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblem(A, b, c, ["C", "B"], variable_type=">=")
            sage: P.is_feasible()
            True
        """
        return self._solve()[1] is not None

    def n_constraints(self):
        r"""
        Return the number of constraints of ``self``, i.e. `m`.

        OUTPUT:

        - an integer

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblem(A, b, c, ["C", "B"], variable_type=">=")
            sage: P.n_constraints()
            2
            sage: P.m()
            2
        """
        return self._Abcx[0].nrows()

    def n_variables(self):
        r"""
        Return the number of decision variables of ``self``, i.e. `n`.

        OUTPUT:

        - an integer

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblem(A, b, c, ["C", "B"], variable_type=">=")
            sage: P.n_variables()
            2
            sage: P.n()
            2
        """
        return self._Abcx[0].ncols()

    def objective_coefficients(self):
        r"""
        Return coefficients of the objective of ``self``, i.e. `c`.

        OUTPUT:

        - a vector

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblem(A, b, c, ["C", "B"], variable_type=">=")
            sage: P.objective_coefficients()
            (10, 5)
            sage: P.c()
            (10, 5)
        """
        return self._Abcx[2]

    def optimal_solution(self):
        r"""
        Return **an** optimal solution of ``self``.

        OUTPUT:

        - a vector or ``None`` if there are no optimal solutions

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblem(A, b, c, ["C", "B"], variable_type=">=")
            sage: P.optimal_solution()
            (250, 750)
        """
        return self._solve()[0]

    def optimal_value(self):
        r"""
        Return the optimal value for ``self``.

        OUTPUT:

        - a number if the problem is bounded, `\pm \infty` if it is unbounded,
          or ``None`` if it is infeasible

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblem(A, b, c, ["C", "B"], variable_type=">=")
            sage: P.optimal_value()
            6250
        """
        return self._solve()[1]

    def plot(self, *args, **kwds):
        r"""
        Return a plot for solving ``self`` graphically.

        INPUT:

        - ``xmin``, ``xmax``, ``ymin``, ``ymax`` -- bounds for the axes, if
          not given, an attempt will be made to pick reasonable values

        - ``alpha`` -- (default: 0.2) determines how opaque are shadows

        OUTPUT:

        - a plot

        This only works for problems with two decision variables. On the plot
        the black arrow indicates the direction of growth of the objective. The
        lines perpendicular to it are level curves of the objective. If there
        are optimal solutions, the arrow originates in one of them and the
        corresponding level curve is solid: all points of the feasible set
        on it are optimal solutions. Otherwise the arrow is placed in the
        center. If the problem is infeasible or the objective is zero, a plot
        of the feasible set only is returned.

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblem(A, b, c, ["C", "B"], variable_type=">=")
            sage: p = P.plot()
            sage: p.show()

        In this case the plot works better with the following axes ranges::

            sage: p = P.plot(0, 1000, 0, 1500)
            sage: p.show()

        TESTS:

        We check that zero objective can be dealt with::

            sage: InteractiveLPProblem(A, b, (0, 0), ["C", "B"], variable_type=">=").plot()
            Graphics object consisting of 8 graphics primitives
        """
        FP = self.plot_feasible_set(*args, **kwds)
        c = self.c().n().change_ring(QQ)
        if c.is_zero():
            return FP
        xmin = FP.xmin()
        xmax = FP.xmax()
        ymin = FP.ymin()
        ymax = FP.ymax()
        xmin, xmax, ymin, ymax = map(QQ, [xmin, xmax, ymin, ymax])
        start = self.optimal_solution()
        start = vector(QQ, start.n() if start is not None
                            else [xmin + (xmax-xmin)/2, ymin + (ymax-ymin)/2])
        length = min(xmax - xmin, ymax - ymin) / 5
        end = start + (c * length / c.norm()).n().change_ring(QQ)
        result = FP + point(start, color="black", size=50, zorder=10)
        result += arrow(start, end, color="black", zorder=10)
        ieqs = [(xmax, -1, 0), (- xmin, 1, 0),
                (ymax, 0, -1), (- ymin, 0, 1)]
        box = Polyhedron(ieqs=ieqs)
        d = vector([c[1], -c[0]])
        for i in range(-10, 11):
            level = Polyhedron(vertices=[start + i*(end-start)], lines=[d])
            level = box.intersection(level)
            if level.vertices():
                if i == 0 and self.is_bounded():
                    result += line(level.vertices(), color="black",
                                   thickness=2)
                else:
                    result += line(level.vertices(), color="black",
                                   linestyle="--")
        result.set_axes_range(xmin, xmax, ymin, ymax)
        result.axes_labels(FP.axes_labels())    #FIXME: should be preserved!
        return result

    def plot_feasible_set(self, xmin=None, xmax=None, ymin=None, ymax=None,
                          alpha=0.2):
        r"""
        Return a plot of the feasible set of ``self``.

        INPUT:

        - ``xmin``, ``xmax``, ``ymin``, ``ymax`` -- bounds for the axes, if
          not given, an attempt will be made to pick reasonable values

        - ``alpha`` -- (default: 0.2) determines how opaque are shadows

        OUTPUT:

        - a plot

        This only works for a problem with two decision variables. The plot
        shows boundaries of constraints with a shadow on one side for
        inequalities. If the :meth:`feasible_set` is not empty and at least
        part of it is in the given boundaries, it will be shaded gray and `F`
        will be placed in its middle.

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblem(A, b, c, ["C", "B"], variable_type=">=")
            sage: p = P.plot_feasible_set()
            sage: p.show()

        In this case the plot works better with the following axes ranges::

            sage: p = P.plot_feasible_set(0, 1000, 0, 1500)
            sage: p.show()
        """
        if self.n() != 2:
            raise ValueError("only problems with 2 variables can be plotted")
        A, b, c, x = self.Abcx()
        if self.base_ring() is not QQ:
            # Either we use QQ or crash
            A = A.n().change_ring(QQ)
            b = b.n().change_ring(QQ)
        F = self.feasible_set()
        if ymax is None:
            ymax = max(map(abs, b) + [v[1] for v in F.vertices()])
        if ymin is None:
            ymin = min([-ymax/4.0] + [v[1] for v in F.vertices()])
        if xmax is None:
            xmax = max([1.5*ymax] + [v[0] for v in F.vertices()])
        if xmin is None:
            xmin = min([-xmax/4.0] + [v[0] for v in F.vertices()])
        xmin, xmax, ymin, ymax = map(QQ, [xmin, xmax, ymin, ymax])
        pad = max(xmax - xmin, ymax - ymin) / 20
        ieqs = [(xmax, -1, 0), (- xmin, 1, 0),
                (ymax, 0, -1), (- ymin, 0, 1)]
        box = Polyhedron(ieqs=ieqs)
        F = box.intersection(F)
        result = Graphics()
        colors = rainbow(self.m() + 2)
        for Ai, ri, bi, color in zip(A.rows(), self._constraint_types,
                                           b, colors[:-2]):
            border = box.intersection(Polyhedron(eqns=[[-bi] + list(Ai)]))
            vertices = border.vertices()
            if not vertices:
                continue
            label = r"${}$".format(_latex_product(Ai, x, " ", tail=[ri, bi]))
            result += line(vertices, color=color, legend_label=label)
            if ri == "<=":
                ieqs = [[bi] + list(-Ai), [-bi+pad*Ai.norm().n()] + list(Ai)]
            elif ri == ">=":
                ieqs = [[-bi] + list(Ai), [bi+pad*Ai.norm().n()] + list(-Ai)]
            else:
                continue
            ieqs = map(lambda ieq: map(QQ, ieq), ieqs)
            halfplane = box.intersection(Polyhedron(ieqs=ieqs))
            result += halfplane.render_solid(alpha=alpha, color=color)
        # Same for variables, but no legend
        for ni, ri, color in zip((QQ**2).gens(), self._variable_types,
                                 colors[-2:]):
            border = box.intersection(Polyhedron(eqns=[[0] + list(ni)]))
            if not border.vertices():
                continue
            if ri == "<=":
                ieqs = [[0] + list(-ni), [pad] + list(ni)]
            elif ri == ">=":
                ieqs = [[0] + list(ni), [pad] + list(-ni)]
            else:
                continue
            ieqs = map(lambda ieq: map(QQ, ieq), ieqs)
            halfplane = box.intersection(Polyhedron(ieqs=ieqs))
            result += halfplane.render_solid(alpha=alpha, color=color)
        if F.vertices():
            result += F.render_solid(alpha=alpha, color="gray")
            result += text("$F$", F.center(),
                           fontsize=20, color="black", zorder=5)
        result.set_axes_range(xmin, xmax, ymin, ymax)
        result.axes_labels(map(lambda xi: "${}$".format(latex(xi)), x))
        result.legend(True)
        result.set_legend_options(fancybox=True, handlelength=1.5, loc=1,
                                  shadow=True)
        result._extra_kwds["aspect_ratio"] = 1
        result.set_aspect_ratio(1)
        return result

    def standard_form(self):
        r"""
        Construct the LP problem in standard form equivalent to ``self``.

        OUTPUT:

        - an :class:`InteractiveLPProblemStandardForm`

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblem(A, b, c, ["C", "B"], variable_type=">=")
            sage: DP = P.dual()
            sage: DPSF = DP.standard_form()
            sage: DPSF.b()
            (-10, -5)
        """
        A, b, c, x = self.Abcx()
        if not all(ct == "<=" for ct in self._constraint_types):
            newA = []
            newb = []
            for ct, Ai, bi in zip(self._constraint_types, A, b):
                if ct in ["<=", "=="]:
                    newA.append(Ai)
                    newb.append(bi)
                if ct in [">=", "=="]:
                    newA.append(-Ai)
                    newb.append(-bi)
            A = matrix(newA)
            b = vector(newb)
        if not all(vt == ">=" for vt in self._variable_types):
            newA = []
            newc = []
            newx = []
            for vt, Aj, cj, xj in zip(self._variable_types, A.columns(), c, x):
                xj = str(xj)
                if vt in [">=", ""]:
                    newA.append(Aj)
                    newc.append(cj)
                if vt == ">=":
                    newx.append(xj)
                if vt == "":
                    newx.append(xj + "_p")
                if vt in ["<=", ""]:
                    newA.append(-Aj)
                    newc.append(-cj)
                    newx.append(xj + "_n")
            A = column_matrix(newA)
            c = vector(newc)
            x = newx
        is_negative = self._is_negative
        if self._problem_type == "min":
            is_negative = not is_negative
            c = - c
        problem_type = "-max" if is_negative else "max"
        return InteractiveLPProblemStandardForm(A, b, c, x, problem_type,
                                     self._prefix, self._prefix + "0")

    # Aliases for the standard notation
    A = constraint_coefficients
    b = constant_terms
    c = objective_coefficients
    x = decision_variables
    m = n_constraints
    n = n_variables

class InteractiveLPProblemStandardForm(InteractiveLPProblem):
    r"""
    Construct an LP (Linear Programming) problem in standard form.

    .. NOTE::

        This class is for **educational purposes only**: if you want to solve
        Linear Programs efficiently, use :class:`MixedIntegerLinearProgram`
        instead.

    The used standard form is:

    .. MATH::

        \begin{array}{l}
        \pm \max cx \\
        Ax \leq b \\
        x \geq 0
        \end{array}

    INPUT:

    - ``A`` -- a matrix of constraint coefficients

    - ``b`` -- a vector of constraint constant terms

    - ``c`` -- a vector of objective coefficients

    - ``x`` -- (default: ``"x"``) a vector of decision variables or a string
      the base name giving

    - ``problem_type`` -- (default: ``"max"``) a string specifying the
      problem type: either ``"max"`` or ``"-max"``

    - ``slack_variables`` -- (default: same as ``x`` parameter, if it was given
      as a string, otherwise string ``"x"``) a vector of slack variables or
      a sting giving the base name

    - ``auxiliary_variable`` -- (default: same as ``x`` parameter with adjoined
      ``"0"`` if it was given as a string, otherwise ``"x0"``) the auxiliary
      name, expected to be the same as the first decision variable for
      auxiliary problems

    - ``objective`` -- (default: ``"z"``) the objective variable (used for the
      initial dictionary)

    - ``base_ring`` -- (default: the fraction field of a common ring for all
      input coefficients) a field to which all input coefficients will be
      converted

    EXAMPLES::

        sage: A = ([1, 1], [3, 1])
        sage: b = (1000, 1500)
        sage: c = (10, 5)
        sage: P = InteractiveLPProblemStandardForm(A, b, c)

    Unlike :class:`InteractiveLPProblem`, this class does not allow you to adjust types of
    constraints (they are always ``"<="``) and variables (they are always
    ``">="``), and the problem type may only be ``"max"`` or ``"-max"``.
    You may give custom names to slack and auxiliary variables, but in
    most cases defaults should work::

        sage: P.decision_variables()
        (x1, x2)
        sage: P.slack_variables()
        (x3, x4)
    """

    def __init__(self, A, b, c, x="x", problem_type="max",
                 slack_variables=None, auxiliary_variable=None, objective="z",
                 base_ring=None):
        r"""
        See :class:`StandardFormLPP` for documentation.

        TESTS::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: TestSuite(P).run()
        """
        if problem_type not in ("max", "-max"):
            raise ValueError("problems in standard form must be of (negative) "
                             "maximization type")
        super(InteractiveLPProblemStandardForm, self).__init__(A, b, c, x,
                                                    problem_type=problem_type,
                                                    constraint_type="<=",
                                                    variable_type=">=",
                                                    base_ring=base_ring)
        n, m = self.n(), self.m()
        if slack_variables is None:
           slack_variables = self._prefix
        if isinstance(slack_variables, str):
            slack_variables = ["{}{:d}".format(slack_variables, i)
                               for i in range(n + 1, n + m + 1)]
        else:
            slack_variables = map(str, slack_variables)
            if len(slack_variables) != m:
                raise ValueError("wrong number of slack variables")
        if auxiliary_variable is None:
           auxiliary_variable = self._prefix + "0"
        names = [str(auxiliary_variable)] + map(str, self.x()) + slack_variables
        if names[0] == names[1]:
            names.pop(0)
        R = PolynomialRing(self.base_ring(), names, order="neglex")
        self._R = R
        x = vector(R.gens()[-n-m:-m])
        x.set_immutable()
        self._Abcx = self._Abcx[:-1] + (x, )
        self._objective = objective

    def auxiliary_problem(self):
        r"""
        Construct the auxiliary problem for ``self``.

        OUTPUT:

        - an :class:`LP problem in standard form <InteractiveLPProblemStandardForm>`

        The auxiliary problem with the auxiliary variable `x_0` is

        .. MATH::

            \begin{array}{l}
            \max - x_0 \\
            - x_0 + A_i x \leq b_i \text{ for all } i \\
            x \geq 0
            \end{array}\ .

        Such problems are used when the :meth:`initial_dictionary` is
        infeasible.

        EXAMPLES::

            sage: A = ([1, 1], [3, 1], [-1, -1])
            sage: b = (1000, 1500, -400)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: AP = P.auxiliary_problem()
        """
        X = self.coordinate_ring().gens()
        m, n = self.m(), self.n()
        if len(X) == m + n:
            raise ValueError("auxiliary variable is already among decision "
                             "ones")
        F = self.base_ring()
        A = column_matrix(F, [-1] * m).augment(self.A())
        c = vector(F, [-1] + [0] * n)
        return InteractiveLPProblemStandardForm(A, self.b(), c, X[:-m],
                                     slack_variables=X[-m:],
                                     auxiliary_variable=X[0],
                                     objective="w")

    def auxiliary_variable(self):
        r"""
        Return the auxiliary variable of ``self``.

        Note that the auxiliary variable may or may not be among
        :meth:`~InteractiveLPProblem.decision_variables`.

        OUTPUT:

        - a variable of the :meth:`coordinate_ring` of ``self``

        EXAMPLES::

            sage: A = ([1, 1], [3, 1], [-1, -1])
            sage: b = (1000, 1500, -400)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: P.auxiliary_variable()
            x0
            sage: P.decision_variables()
            (x1, x2)
            sage: AP = P.auxiliary_problem()
            sage: AP.auxiliary_variable()
            x0
            sage: AP.decision_variables()
            (x0, x1, x2)
        """
        return self._R.gen(0)

    def coordinate_ring(self):
        r"""
        Return the coordinate ring of ``self``.

        OUTPUT:

        - a polynomial ring over the :meth:`~InteractiveLPProblem.base_ring` of ``self`` in
          the :meth:`auxiliary_variable`, :meth:`~InteractiveLPProblem.decision_variables`,
          and :meth:`slack_variables` with "neglex" order

        EXAMPLES::

            sage: A = ([1, 1], [3, 1], [-1, -1])
            sage: b = (1000, 1500, -400)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: P.coordinate_ring()
            Multivariate Polynomial Ring in x0, x1, x2, x3, x4, x5
            over Rational Field
            sage: P.base_ring()
            Rational Field
            sage: P.auxiliary_variable()
            x0
            sage: P.decision_variables()
            (x1, x2)
            sage: P.slack_variables()
            (x3, x4, x5)
        """
        return self._R

    def dictionary(self, *x_B):
        r"""
        Construct a dictionary for ``self`` with given basic variables.

        INPUT:

        - basic variables for the dictionary to be constructed

        OUTPUT:

        - a :class:`dictionary <LPDictionary>`

        .. NOTE::

            This is a synonym for ``self.revised_dictionary(x_B).dictionary()``,
            but basic variables are mandatory.

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: D = P.dictionary("x1", "x2")
            sage: D.basic_variables()
            (x1, x2)
        """
        if not x_B:
            raise ValueError("basic variables must be given explicitly")
        return self.revised_dictionary(*x_B).dictionary()

    def feasible_dictionary(self, auxiliary_dictionary):
        r"""
        Construct a feasible dictionary for ``self``.

        INPUT:

        - ``auxiliary_dictionary`` -- an optimal dictionary for the
          :meth:`auxiliary_problem` of ``self`` with the optimal value `0` and
          a non-basic auxiliary variable

        OUTPUT:

        - a feasible :class:`dictionary <LPDictionary>` for ``self``

        EXAMPLES::

            sage: A = ([1, 1], [3, 1], [-1, -1])
            sage: b = (1000, 1500, -400)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: AP = P.auxiliary_problem()
            sage: D = AP.initial_dictionary()
            sage: D.enter(0)
            sage: D.leave(5)
            sage: D.update()
            sage: D.enter(1)
            sage: D.leave(0)
            sage: D.update()
            sage: D.is_optimal()
            True
            sage: D.objective_value()
            0
            sage: D.basic_solution()
            (0, 400, 0)
            sage: D = P.feasible_dictionary(D)
            sage: D.is_optimal()
            False
            sage: D.is_feasible()
            True
            sage: D.objective_value()
            4000
            sage: D.basic_solution()
            (400, 0)
        """
        # It is good to have sanity checks in this function, but they are a bit
        # problematic with numerical dictionaries, so we do only few.
        x0 = self.auxiliary_variable()
        if x0 not in auxiliary_dictionary.nonbasic_variables():
            raise ValueError("the auxiliary variable must be non-basic")
        if not auxiliary_dictionary.is_feasible():
            raise ValueError("the auxiliary dictionary must be feasible")
        A, b, c, v, B, N, z = auxiliary_dictionary._AbcvBNz
        B = tuple(B)
        N = tuple(N)
        k = N.index(x0)
        N = N[:k] + N[k+1:]
        n = len(c)
        A = A.matrix_from_columns(range(k) + range(k + 1, n))
        b = copy(b)
        c = vector(self.base_ring(), n - 1)
        for cj, xj in zip(*self.Abcx()[-2:]):
            if xj in N:
                c[N.index(xj)] += cj
            else:
                i = B.index(xj)
                c -= cj * A[i]
                v += cj * b[i]
        B = map(self._R, B)
        N = map(self._R, N)
        return LPDictionary(A, b, c, v, B, N, self._objective)

    def final_dictionary(self):
        r"""
        Return the final dictionary of the simplex method applied to ``self``.

        See :meth:`run_simplex_method` for the description of possibilities.

        OUTPUT:

        - a :class:`dictionary <LPDictionary>`

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: D = P.final_dictionary()
            sage: D.is_optimal()
            True

        TESTS::

            sage: P.final_dictionary() is P.final_dictionary()
            False
        """
        try:
            D = self._final_dictionary
            # Since dictionaries are "highly mutable", forget the returned one.
            del self._final_dictionary
            return D
        except AttributeError:
            self.run_simplex_method()
            return self.final_dictionary()

    def final_revised_dictionary(self):
        r"""
        Return the final dictionary of the revised simplex method applied
        to ``self``.

        See :meth:`run_revised_simplex_method` for the description of
        possibilities.

        OUTPUT:

        - a :class:`revised dictionary <LPRevisedDictionary>`

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: D = P.final_revised_dictionary()
            sage: D.is_optimal()
            True

        TESTS::

            sage: P.final_revised_dictionary() is P.final_revised_dictionary()
            False
        """
        try:
            D = self._final_revised_dictionary
            # Since dictionaries are "highly mutable", forget the returned one.
            del self._final_revised_dictionary
            return D
        except AttributeError:
            self.run_revised_simplex_method()
            return self.final_revised_dictionary()

    def initial_dictionary(self):
        r"""
        Construct the initial dictionary of ``self``.

        The initial dictionary "defines" :meth:`slack_variables` in terms
        of the :meth:`~InteractiveLPProblem.decision_variables`, i.e. it has slack
        variables as basic ones.

        OUTPUT:

        - a :class:`dictionary <LPDictionary>`

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: D = P.initial_dictionary()
        """
        A, b, c, x = self.Abcx()
        x = self._R.gens()
        m, n = self.m(), self.n()
        return LPDictionary(A, b, c, 0, x[-m:], x[-m-n:-m], self._objective)

    def inject_variables(self, scope=None, verbose=True):
        r"""
        Inject variables of ``self`` into ``scope``.

        INPUT:

        - ``scope`` -- namespace (default: global)

        - ``verbose`` -- if ``True`` (default), names of injected variables
          will be printed

        OUTPUT:

        - none

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: P.inject_variables()
            Defining x0, x1, x2, x3, x4
            sage: 3*x1 + x2
            x2 + 3*x1
        """
        if scope is None:
            scope = get_main_globals()
        try:
            self._R.inject_variables(scope, verbose)
        except AttributeError:
            pass

    def revised_dictionary(self, *x_B):
        r"""
        Construct a revised dictionary for ``self``.

        INPUT:

        - basic variables for the dictionary to be constructed; if not given,
          :meth:`slack_variables` will be used, perhaps with the
          :meth:`auxiliary_variable` to give a feasible dictionary

        OUTPUT:

        - a :class:`revised dictionary <LPRevisedDictionary>`

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: D = P.revised_dictionary("x1", "x2")
            sage: D.basic_variables()
            (x1, x2)

        If basic variables are not given the initial dictionary is
        constructed::

            sage: P.revised_dictionary().basic_variables()
            (x3, x4)
            sage: P.initial_dictionary().basic_variables()
            (x3, x4)

        Unless it is infeasible, in which case a feasible dictionary for the
        auxiliary problem is constructed::

            sage: A = ([1, 1], [3, 1], [-1,-1])
            sage: b = (1000, 1500, -400)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: P.initial_dictionary().is_feasible()
            False
            sage: P.revised_dictionary().basic_variables()
            (x3, x4, x0)
        """
        if not x_B:
            x_B = list(self.slack_variables())
            bm = min(self.b())
            if bm < 0:
                x_B[self.b().list().index(bm)] = self.auxiliary_variable()
        return LPRevisedDictionary(self, x_B)

    def run_revised_simplex_method(self):
        r"""
        Apply the revised simplex method to solve ``self`` and show the steps.

        OUTPUT:

        - a string with `\LaTeX` code of intermediate dictionaries

        .. NOTE::

            You can access the :meth:`final_revised_dictionary`, which can be
            one of the following:

            - an optimal dictionary with the :meth:`auxiliary_variable` among
              :meth:`~LPRevisedDictionary.basic_variables` and a non-zero
              optimal value indicating
              that ``self`` is infeasible;

            - a non-optimal dictionary that has marked entering
              variable for which there is no choice of the leaving variable,
              indicating that ``self`` is unbounded;

            - an optimal dictionary.

        EXAMPLES::

            sage: A = ([1, 1], [3, 1], [-1, -1])
            sage: b = (1000, 1500, -400)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: P.run_revised_simplex_method()
            %notruncate
            \renewcommand{\arraystretch}{1.500000}
            \begin{array}{l}
            ...
            \text{Entering: $x_{1}$. Leaving: $x_{0}$.}\\
            ...
            \text{Entering: $x_{5}$. Leaving: $x_{4}$.}\\
            ...
            \text{Entering: $x_{2}$. Leaving: $x_{3}$.}\\
            ...
            \text{The optimal value: $6250$. An optimal solution: $\left(250,\,750\right)$.}
            \end{array}
        """
        result = []
        d = self.revised_dictionary()
        while not d.is_optimal():
            entering, leaving = min(d.possible_simplex_method_steps())
            d.enter(entering)
            if leaving:
                leaving = min(leaving)
                d.leave(leaving)
                result.append(latex(d))
                result.append(r"\text{{Entering: ${}$. Leaving: ${}$.}}"
                              .format(latex(entering), latex(leaving)))
                result.append("")
                result.append(
                        r"B_\mathrm{new}^{-1} = E^{-1} B_\mathrm{old}^{-1} = " +
                        latex(d.E_inverse()) + latex(d.B_inverse()))
                result.append("")
                d.update()
            else:
                result.append(latex(d))
                result.append(r"\text{{The problem is unbounded in the "
                              r"${}$ direction.}}".format(latex(entering)))
                break
        if d.is_optimal():
            result.append(latex(d))
            if self.auxiliary_variable() in d.basic_variables():
                result.append(r"\text{The problem is infeasible.}")
            else:
                v = d.objective_value()
                if self._is_negative:
                    v = - v
                result.append((r"\text{{The optimal value: ${}$. "
                               "An optimal solution: ${}$.}}").format(
                               latex(v), latex(d.basic_solution())))
        self._final_revised_dictionary = d
        return _assemble_arrayl(result, 1.5)

    def run_simplex_method(self):
        r"""
        Apply the simplex method to solve ``self`` and show the steps.

        OUTPUT:

        - a string with `\LaTeX` code of intermediate dictionaries

        .. NOTE::

            You can access the :meth:`final_dictionary`, which can be one
            of the following:

            - an optimal dictionary for the :meth:`auxiliary_problem` with a
              non-zero optimal value indicating that ``self`` is infeasible;

            - a non-optimal dictionary for ``self`` that has marked entering
              variable for which there is no choice of the leaving variable,
              indicating that ``self`` is unbounded;

            - an optimal dictionary for ``self``.

        EXAMPLES::

            sage: A = ([1, 1], [3, 1], [-1, -1])
            sage: b = (1000, 1500, -400)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: P.run_simplex_method()  # not tested

        You should use the typeset mode as the command above generates long
        `\LaTeX` code::

            sage: print P.run_simplex_method()
            %notruncate
            ...
            \text{The initial dictionary is infeasible, solving auxiliary problem.}\\
            ...
            \text{Entering: $x_{0}$. Leaving: $x_{5}$.}\\
            ...
            \text{Entering: $x_{1}$. Leaving: $x_{0}$.}\\
            ...
            \text{Back to the original problem.}\\
            ...
            \text{Entering: $x_{5}$. Leaving: $x_{4}$.}\\
            ...
            \text{Entering: $x_{2}$. Leaving: $x_{3}$.}\\
            ...
            \text{The optimal value: $6250$. An optimal solution: $\left(250,\,750\right)$.}
            ...
        """
        result = []
        d = self.initial_dictionary()
        result.append(latex(d))

        def step(entering, leaving):
            result.append(r"\text{{Entering: ${}$. Leaving: ${}$.}}"
                          .format(latex(entering), latex(leaving)))
            result.append(d.ELLUL(entering, leaving))

        if d.is_feasible():
            is_feasible = True
        else:
            result.append(r"\text{The initial dictionary is infeasible, "
                          "solving auxiliary problem.}")
            d = self.auxiliary_problem().initial_dictionary()
            result.append(latex(d))
            x0 = self.auxiliary_variable()
            _, leaving = min(zip(d.constant_terms(), d.basic_variables()))
            step(x0, leaving)
            # while not d.is_optimal():
            # either optimality check should handle rounding errors or
            while not (d.is_optimal() or x0 in d.nonbasic_variables()):
                entering, leaving = min(d.possible_simplex_method_steps())
                step(entering, min(leaving))
            is_feasible = x0 in d.nonbasic_variables()
            if is_feasible:
                result.append(r"\text{Back to the original problem.}")
                d = self.feasible_dictionary(d)
                result.append(latex(d))
            else:
                result.append(r"\text{The original problem is infeasible.}")
        if is_feasible:
            while not d.is_optimal():
                entering, leaving = min(d.possible_simplex_method_steps())
                if leaving:
                    step(entering, min(leaving))
                else:
                    d.enter(entering)
                    result.append(r"\text{{The problem is unbounded in the "
                                  r"${}$ direction.}}".format(latex(entering)))
                    result.append(latex(d))
                    break
            if d.is_optimal():
                v = d.objective_value()
                if self._is_negative:
                    v = - v
                result.append((r"\text{{The optimal value: ${}$. "
                               "An optimal solution: ${}$.}}").format(
                               latex(v), latex(d.basic_solution())))
        self._final_dictionary = d
        if generate_real_LaTeX:
            # This will have to be used via \sagestr, as line&display breaks
            # don't get along with reference substitution without wrapping.
            return ("\\begin{gather*}\n\\allowdisplaybreaks\n" +
                    "\\displaybreak[0]\\\\\n".join(result) +
                    "\n\\end{gather*}")
        return _assemble_arrayl(result, 1.5)

    def slack_variables(self):
        r"""
        Return slack variables of ``self``.

        Slack variables are differences between the constant terms and
        left hand sides of the constraints.

        If you want to give custom names to slack variables, you have to do so
        during construction of the problem.

        OUTPUT:

        - a tuple

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: P.slack_variables()
            (x3, x4)
            sage: P = InteractiveLPProblemStandardForm(A, b, c, ["C", "B"],
            ....:     slack_variables=["L", "F"])
            sage: P.slack_variables()
            (L, F)
        """
        return self._R.gens()[-self.m():]


class LPAbstractDictionary(SageObject):
    r"""
    Abstract base class for dictionaries for LP problems.

    Instantiating this class directly is meaningless, see :class:`LPDictionary`
    and :class:`LPRevisedDictionary` for useful extensions.
    """

    def __init__(self):
        r"""
        Initialize internal fields for entering and leaving variables.

        TESTS::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: D = P.initial_dictionary()    # indirect doctest
        """
        super(LPAbstractDictionary, self).__init__()
        self._entering = None
        self._leaving = None

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        OUTPUT:

        - a string

        TESTS::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: D = P.initial_dictionary()
            sage: print D._repr_()
            LP problem dictionary (use typeset mode to see details)
            sage: D = P.revised_dictionary()
            sage: print D._repr_()
            LP problem dictionary (use typeset mode to see details)
        """
        return "LP problem dictionary (use typeset mode to see details)"

    def base_ring(self):
        r"""
        Return the base ring of ``self``, i.e. the ring of coefficients.

        OUTPUT:

        - a ring

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: D = P.initial_dictionary()
            sage: D.base_ring()
            Rational Field
            sage: D = P.revised_dictionary()
            sage: D.base_ring()
            Rational Field
        """
        return self.coordinate_ring().base_ring()

    def basic_solution(self, include_slack_variables=False):
        r"""
        Return the basic solution of ``self``.

        The basic solution associated to a dictionary is obtained by setting to
        zero all :meth:`~LPDictionary.nonbasic_variables`, in which case
        :meth:`~LPDictionary.basic_variables` have to be equal to
        :meth:`~LPDictionary.constant_terms` in equations.
        It may refer to values of :meth:`~InteractiveLPProblem.decision_variables` only or
        include :meth:`~InteractiveLPProblemStandardForm.slack_variables` as well.

        INPUT:

        - ``include_slack_variables`` -- (default: ``False``) if ``True``,
          values of slack variables will be appended at the end

        OUTPUT:

        - a vector

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: D = P.initial_dictionary()
            sage: D.basic_solution()
            (0, 0)
            sage: D.basic_solution(True)
            (0, 0, 1000, 1500)
            sage: D = P.revised_dictionary()
            sage: D.basic_solution()
            (0, 0)
            sage: D.basic_solution(True)
            (0, 0, 1000, 1500)
        """
        vv = zip(self.basic_variables(), self.constant_terms())
        N = self.nonbasic_variables()
        vv += [(v, 0) for v in N]
        vv.sort()   # We use neglex order
        v = [value for _, value in vv]
        return vector(self.base_ring(),
                      v if include_slack_variables else v[:len(N)])

    def coordinate_ring(self):
        r"""
        Return the coordinate ring of ``self``.

        OUTPUT:

        - a polynomial ring in
          :meth:`~InteractiveLPProblemStandardForm.auxiliary_variable`,
          :meth:`~InteractiveLPProblem.decision_variables`, and
          :meth:`~InteractiveLPProblemStandardForm.slack_variables` of ``self`` over the
          :meth:`base_ring`

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: D = P.initial_dictionary()
            sage: D.coordinate_ring()
            Multivariate Polynomial Ring in x0, x1, x2, x3, x4
            over Rational Field
            sage: D = P.revised_dictionary()
            sage: D.coordinate_ring()
            Multivariate Polynomial Ring in x0, x1, x2, x3, x4
            over Rational Field
        """
        return self.basic_variables()[0].parent()

    def dual_ratios(self):
        r"""
        Return ratios used to determine the entering variable based on leaving.

        OUTPUT:

        - A list of pairs `(r_j, x_j)` where `x_j` is a non-basic variable and
          `r_j = c_j / a_{ij}` is the ratio of the objective coefficient `c_j`
          to the coefficient `a_{ij}` of `x_j` in the relation for the leaving
          variable `x_i`:

          .. MATH::

              x_i = b_i - \cdots - a_{ij} x_j - \cdots.

          The order of pairs matches the order of
          :meth:`~LPDictionary.nonbasic_variables`,
          but only `x_j` with negative `a_{ij}` are considered.

        EXAMPLES::

            sage: A = ([1, 1], [3, 1], [-1, -1])
            sage: b = (1000, 1500, -400)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: D = P.dictionary(2, 3, 5)
            sage: D.leave(3)
            sage: D.dual_ratios()
            [(5/2, x1), (5, x4)]
            sage: D = P.revised_dictionary(2, 3, 5)
            sage: D.leave(3)
            sage: D.dual_ratios()
            [(5/2, x1), (5, x4)]
        """
        return [(c / a, x) for c, a, x in zip(self.objective_coefficients(),
                                              self.leaving_coefficients(),
                                            self.nonbasic_variables()) if a < 0]

    def enter(self, v):
        r"""
        Set ``v`` as the entering variable of ``self``.

        INPUT:

        - ``v`` -- a non-basic variable of ``self``, can be given as a string,
          an actual variable, or an integer interpreted as the index of a
          variable

        OUTPUT:

        - none, but the selected variable will be used as entering by methods
          that require an entering variable and the corresponding column
          will be typeset in green

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: D = P.initial_dictionary()
            sage: D.enter("x1")

        We can also use indices of variables::

            sage: D.enter(1)

        Or variable names without quotes after injecting them::

            sage: P.inject_variables()
            Defining x0, x1, x2, x3, x4
            sage: D.enter(x1)

        The same works for revised dictionaries as well::

            sage: D = P.revised_dictionary()
            sage: D.enter(x1)
        """
        v = variable(self.coordinate_ring(), v)
        if v not in self.nonbasic_variables():
            raise ValueError("entering variable must be non-basic")
        self._entering = v

    def is_dual_feasible(self):
        r"""
        Check if ``self`` is dual feasible.

        OUTPUT:

        - ``True`` if all :meth:`~LPDictionary.objective_coefficients` are
          non-positive, ``False`` otherwise

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: D = P.initial_dictionary()
            sage: D.is_dual_feasible()
            False
            sage: D = P.revised_dictionary()
            sage: D.is_dual_feasible()
            False
        """
        return all(ci <= 0 for ci in self.objective_coefficients())

    def is_feasible(self):
        r"""
        Check if ``self`` is feasible.

        OUTPUT:

        - ``True`` if all :meth:`~LPDictionary.constant_terms` are
          non-negative,``False`` otherwise

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: D = P.initial_dictionary()
            sage: D.is_feasible()
            True
            sage: D = P.revised_dictionary()
            sage: D.is_feasible()
            True
        """
        return all(bi >= 0 for bi in self.constant_terms())

    def is_optimal(self):
        r"""
        Check if ``self`` is optimal.

        OUTPUT:

        - ``True`` if ``self`` :meth:`is_feasible` and :meth:`is_dual_feasible`
          (i.e. all :meth:`~LPDictionary.constant_terms` are non-negative and
          all :meth:`~LPDictionary.objective_coefficients` are non-positive),
          ``False`` otherwise.

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: D = P.initial_dictionary()
            sage: D.is_optimal()
            False
            sage: D = P.revised_dictionary()
            sage: D.is_optimal()
            False
            sage: D = P.revised_dictionary(1, 2)
            sage: D.is_optimal()
            True
        """
        return self.is_feasible() and self.is_dual_feasible()

    def leave(self, v):
        r"""
        Set ``v`` as the leaving variable of ``self``.

        INPUT:

        - ``v`` -- a basic variable of ``self``, can be given as a string, an
          actual variable, or an integer interpreted as the index of a variable

        OUTPUT:

        - none, but the selected variable will be used as leaving by methods
          that require a leaving variable and the corresponding row will be
          typeset in red

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: D = P.initial_dictionary()
            sage: D.leave("x4")

        We can also use indices of variables::

            sage: D.leave(4)

        Or variable names without quotes after injecting them::

            sage: P.inject_variables()
            Defining x0, x1, x2, x3, x4
            sage: D.leave(x4)

        The same works for revised dictionaries as well::

            sage: D = P.revised_dictionary()
            sage: D.leave(x4)
        """
        v = variable(self.coordinate_ring(), v)
        if v not in self.basic_variables():
            raise ValueError("leaving variable must be basic")
        self._leaving = v

    def possible_dual_simplex_method_steps(self):
        r"""
        Return possible dual simplex method steps for ``self``.

        OUTPUT:

        - A list of pairs ``(leaving, entering)``, where ``leaving`` is a
          basic variable that may :meth:`leave` and ``entering`` is a list of
          non-basic variables that may :meth:`enter` when ``leaving`` leaves.
          Note that ``entering`` may be empty, indicating that the problem is
          infeasible (since the dual one is unbounded).

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: D = P.dictionary(2, 3)
            sage: D.possible_dual_simplex_method_steps()
            [(x3, [x1])]
            sage: D = P.revised_dictionary(2, 3)
            sage: D.possible_dual_simplex_method_steps()
            [(x3, [x1])]
        """
        if not self.is_dual_feasible():
            raise ValueError("dual simplex method steps are applicable to "
                             "dual feasible dictionaries only")
        steps = []
        old_entering = self._entering
        self._entering = None
        old_leaving = self._leaving
        for l in self.possible_leaving():
            self.leave(l)
            steps.append((l, self.possible_entering()))
        self._entering = old_entering
        self._leaving = old_leaving
        return steps

    def possible_entering(self):
        r"""
        Return possible entering variables for ``self``.

        OUTPUT:

        - a list of non-basic variables of ``self`` that can :meth:`enter` on
          the next step of the (dual) simplex method

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: D = P.initial_dictionary()
            sage: D.possible_entering()
            [x1, x2]
            sage: D = P.revised_dictionary()
            sage: D.possible_entering()
            [x1, x2]
        """
        if self.is_dual_feasible() and self._leaving is not None:
            ratios = self.dual_ratios()
            if not ratios:
                return []
            min_ratio = min(ratios)[0]
            return [v for r, v in ratios if r == min_ratio]
        if self.is_feasible():
            return [v for c, v in zip(self.objective_coefficients(),
                                      self.nonbasic_variables()) if c > 0]
        raise ValueError("entering variables can be determined for feasible "
                         "dictionaries or for dual feasible dictionaries "
                         "with a set leaving variable")

    def possible_leaving(self):
        r"""
        Return possible leaving variables for ``self``.

        OUTPUT:

        - a list of basic variables of ``self`` that can :meth:`leave` on
          the next step of the (dual) simplex method

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: D = P.initial_dictionary()
            sage: D.enter(1)
            sage: D.possible_leaving()
            [x4]
            sage: D = P.revised_dictionary()
            sage: D.enter(1)
            sage: D.possible_leaving()
            [x4]
        """
        if self.is_feasible() and self._entering is not None:
            ratios = self.ratios()
            if not ratios:
                return []
            min_ratio = min(ratios)[0]
            return [v for r, v in ratios if r == min_ratio]
        if self.is_dual_feasible():
            return [v for b, v in zip(self.constant_terms(),
                                      self.basic_variables()) if b < 0]
        raise ValueError("leaving variables can be determined for feasible "
                         "dictionaries with a set entering variable "
                         "or for dual feasible dictionaries")

    def possible_simplex_method_steps(self):
        r"""
        Return possible simplex method steps for ``self``.

        OUTPUT:

        - A list of pairs ``(entering, leaving)``, where ``entering`` is a
          non-basic variable that may :meth:`enter` and ``leaving`` is a list
          of basic variables that may :meth:`leave` when ``entering`` enters.
          Note that ``leaving`` may be empty, indicating that the problem is
          unbounded.

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: D = P.initial_dictionary()
            sage: D.possible_simplex_method_steps()
            [(x1, [x4]), (x2, [x3])]
            sage: D = P.revised_dictionary()
            sage: D.possible_simplex_method_steps()
            [(x1, [x4]), (x2, [x3])]
        """
        if not self.is_feasible():
            raise ValueError("simplex method steps are applicable to feasible "
                             "dictionaries only")
        steps = []
        old_entering = self._entering
        old_leaving = self._leaving
        self._leaving = None
        for e in self.possible_entering():
            self.enter(e)
            steps.append((e, self.possible_leaving()))
        self._entering = old_entering
        self._leaving = old_leaving
        return steps

    def ratios(self):
        r"""
        Return ratios used to determine the leaving variable based on entering.

        OUTPUT:

        - A list of pairs `(r_i, x_i)` where `x_i` is a basic variable and
          `r_i = b_i / a_{ik}` is the ratio of the constant term `b_i` to the
          coefficient `a_{ik}` of the entering variable `x_k` in the relation
          for `x_i`:

          .. MATH::

              x_i = b_i - \cdots - a_{ik} x_k - \cdots.

          The order of pairs matches the order of
          :meth:`~LPDictionary.basic_variables`,
          but only `x_i` with positive `a_{ik}` are considered.

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: D = P.initial_dictionary()
            sage: D.enter(1)
            sage: D.ratios()
            [(1000, x3), (500, x4)]
            sage: D = P.revised_dictionary()
            sage: D.enter(1)
            sage: D.ratios()
            [(1000, x3), (500, x4)]
        """
        return [(b / a, x) for b, a, x in zip(self.constant_terms(),
                                              self.entering_coefficients(),
                                              self.basic_variables()) if a > 0]


class LPDictionary(LPAbstractDictionary):
    r"""
    Construct a dictionary for an LP problem.

    A dictionary consists of the following data:

    .. MATH::

        \begin{array}{|l|}
        \hline
        x_B = b - A x_N\\
        \hline
        z = z^* + c x_N\\
        \hline
        \end{array}

    INPUT:

    - ``A`` -- a matrix of relation coefficients

    - ``b`` -- a vector of relation constant terms

    - ``c`` -- a vector of objective coefficients

    - ``objective_value`` -- current value of the objective `z^*`

    - ``basic_variables`` -- a list of basic variables `x_B`

    - ``nonbasic_variables`` -- a list of non-basic variables `x_N`

    - ``objective_variable`` -- an objective variable `z`

    OUTPUT:

    - a :class:`dictionary for an LP problem <LPDictionary>`

    .. NOTE::

        This constructor does not check correctness of input, as it is
        intended to be used internally by :class:`InteractiveLPProblemStandardForm`.

    EXAMPLES:

    The intended way to use this class is indirect::

        sage: A = ([1, 1], [3, 1])
        sage: b = (1000, 1500)
        sage: c = (10, 5)
        sage: P = InteractiveLPProblemStandardForm(A, b, c)
        sage: D = P.initial_dictionary()
        sage: D
        LP problem dictionary (use typeset mode to see details)

    But if you want you can create a dictionary without starting with an LP
    problem, here is construction of the same dictionary as above::

        sage: A = matrix(QQ, ([1, 1], [3, 1]))
        sage: b = vector(QQ, (1000, 1500))
        sage: c = vector(QQ, (10, 5))
        sage: R = PolynomialRing(QQ, "x1, x2, x3, x4", order="neglex")
        sage: from sage.numerical.interactive_simplex_method \
        ....:     import LPDictionary
        sage: D2 = LPDictionary(A, b, c, 0, R.gens()[2:], R.gens()[:2], "z")
        sage: D2 == D
        True
    """

    def __init__(self, A, b, c, objective_value,
                 basic_variables, nonbasic_variables, objective_variable):
        r"""
        See :class:`LPDictionary` for documentation.

        TESTS:::

            sage: A = matrix(QQ, ([1, 1], [3, 1]))
            sage: b = vector(QQ, (1000, 1500))
            sage: c = vector(QQ, (10, 5))
            sage: R = PolynomialRing(QQ, "x1, x2, x3, x4", order="neglex")
            sage: from sage.numerical.interactive_simplex_method \
            ....:     import LPDictionary
            sage: D = LPDictionary(A, b, c, 0, R.gens()[2:], R.gens()[:2], "z")
            sage: TestSuite(D).run()
        """
        super(LPDictionary, self).__init__()
        # We are going to change stuff while InteractiveLPProblem has immutable data.
        A = copy(A)
        b = copy(b)
        c = copy(c)
        B = vector(basic_variables)
        N = vector(nonbasic_variables)
        self._AbcvBNz = [A, b, c, objective_value, B, N, SR(objective_variable)]

    def __eq__(self, other):
        r"""
        Check if two LP problem dictionaries are equal.

        INPUT:

        - ``other`` -- anything

        OUTPUT:

        - ``True`` if ``other`` is an :class:`LPDictionary` with all
          details the same as ``self``, ``False`` otherwise.

        TESTS::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: D = P.initial_dictionary()

            sage: A = matrix(QQ, ([1, 1], [3, 1]))
            sage: b = vector(QQ, (1000, 1500))
            sage: c = vector(QQ, (10, 5))
            sage: R = PolynomialRing(QQ, "x1, x2, x3, x4", order="neglex")
            sage: from sage.numerical.interactive_simplex_method \
            ....:     import LPDictionary
            sage: D2 = LPDictionary(A, b, c, 0, R.gens()[2:], R.gens()[:2], "z")
            sage: D2 == D
            True

            sage: D3 = LPDictionary(A, b, c, 0, R.gens()[2:], R.gens()[:2], "w")
            sage: D2 == D3
            False
        """
        return (isinstance(other, LPDictionary) and
                self._AbcvBNz == other._AbcvBNz)

    def _latex_(self):
        r"""
        Return a LaTeX representation of ``self``.

        OUTPUT:

        - a string

        TESTS::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: D = P.initial_dictionary()
            sage: print D._latex_()
            \renewcommand{\arraystretch}{1.5}
            \begin{array}{|rcrcrcr|}
            \hline
            x_{3} \mspace{-6mu}&\mspace{-6mu} = \mspace{-6mu}&\mspace{-6mu} 1000 \mspace{-6mu}&\mspace{-6mu} - \mspace{-6mu}&\mspace{-6mu} x_{1} \mspace{-6mu}&\mspace{-6mu} - \mspace{-6mu}&\mspace{-6mu} x_{2}\\
            x_{4} \mspace{-6mu}&\mspace{-6mu} = \mspace{-6mu}&\mspace{-6mu} 1500 \mspace{-6mu}&\mspace{-6mu} - \mspace{-6mu}&\mspace{-6mu} 3 x_{1} \mspace{-6mu}&\mspace{-6mu} - \mspace{-6mu}&\mspace{-6mu} x_{2}\\
            \hline
            z \mspace{-6mu}&\mspace{-6mu} = \mspace{-6mu}&\mspace{-6mu} 0 \mspace{-6mu}&\mspace{-6mu} + \mspace{-6mu}&\mspace{-6mu} 10 x_{1} \mspace{-6mu}&\mspace{-6mu} + \mspace{-6mu}&\mspace{-6mu} 5 x_{2}\\
            \hline
            \end{array}
        """
        A, b, c, v, B, N, z = self._AbcvBNz
        lines = []
        lines.append(r"\renewcommand{\arraystretch}{1.5}")
        if generate_real_LaTeX:
            lines[-1] += r" \setlength{\arraycolsep}{0.125em}"
#        else:
#            lines[-1] += r"\require{color}"
        lines.append(r"\begin{array}{|rcr%s|}" % ("cr"*len(N)))
        lines.append(r"\hline")
        for xi, bi, Ai in zip(B, b, A.rows()):
            lines.append(_latex_product(-Ai,N, head=[xi, "=", bi],
                                    drop_plus=False, allow_empty=True) + r"\\")
        lines.append(r"\hline")
        lines.append(_latex_product(c, N, head=[z, "=", v],
                                    drop_plus=False, allow_empty=True) + r"\\")
        lines.append(r"\hline")
        lines.append(r"\end{array}")
        latex.add_package_to_preamble_if_available("color")
        if self._entering is not None:
            # Highlight the entering variable column
            e = 2 * tuple(N).index(self._entering) + 4
            for i, line in enumerate(lines):
                line = line.split("&")
                if len(line) > 1:
                    line[e] = r"\color{green}" + line[e]
                    lines[i] = "&".join(line)
        if self._leaving is not None:
            # Highlight the leaving variable row
            l = tuple(B).index(self._leaving) + 3
            line = lines[l].split("&")
            for i, term in enumerate(line):
                line[i] = r"\color{red}" + term
            line = "&".join(line)
            line = line.replace(r"\color{red}\color{green}", r"\color{blue}")
            lines[l] = line
        return  "\n".join(lines)

    def ELLUL(self, entering, leaving):
        r"""
        Perform the Enter-Leave-LaTeX-Update-LaTeX step sequence on ``self``.

        INPUT:

        - ``entering`` -- the entering variable

        - ``leaving`` -- the leaving variable

        OUTPUT:

        - a string with LaTeX code for ``self`` before and after update

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: D = P.initial_dictionary()
            sage: D.ELLUL("x1", "x4")
            \renewcommand{\arraystretch}{1.5}
            \begin{array}{|rcrcrcr|}
            \hline
            x_{3} \mspace{-6mu}&\mspace{-6mu} = \mspace{-6mu}&\mspace{-6mu} 1000 \mspace{-6mu}&\mspace{-6mu} - \mspace{-6mu}&\color{green}\mspace{-6mu} x_{1} \mspace{-6mu}&\mspace{-6mu} - \mspace{-6mu}&\mspace{-6mu} x_{2}\\
            \color{red}x_{4} \mspace{-6mu}&\color{red}\mspace{-6mu} = \mspace{-6mu}&\color{red}\mspace{-6mu} 1500 \mspace{-6mu}&\color{red}\mspace{-6mu} - \mspace{-6mu}&\color{blue}\mspace{-6mu} 3 x_{1} \mspace{-6mu}&\color{red}\mspace{-6mu} - \mspace{-6mu}&\color{red}\mspace{-6mu} x_{2}\\
            \hline
            z \mspace{-6mu}&\mspace{-6mu} = \mspace{-6mu}&\mspace{-6mu} 0 \mspace{-6mu}&\mspace{-6mu} + \mspace{-6mu}&\color{green}\mspace{-6mu} 10 x_{1} \mspace{-6mu}&\mspace{-6mu} + \mspace{-6mu}&\mspace{-6mu} 5 x_{2}\\
            \hline
            \\
            \hline
            x_{3} \mspace{-6mu}&\mspace{-6mu} = \mspace{-6mu}&\mspace{-6mu} 500 \mspace{-6mu}&\mspace{-6mu} + \mspace{-6mu}&\mspace{-6mu} \frac{1}{3} x_{4} \mspace{-6mu}&\mspace{-6mu} - \mspace{-6mu}&\mspace{-6mu} \frac{2}{3} x_{2}\\
            x_{1} \mspace{-6mu}&\mspace{-6mu} = \mspace{-6mu}&\mspace{-6mu} 500 \mspace{-6mu}&\mspace{-6mu} - \mspace{-6mu}&\mspace{-6mu} \frac{1}{3} x_{4} \mspace{-6mu}&\mspace{-6mu} - \mspace{-6mu}&\mspace{-6mu} \frac{1}{3} x_{2}\\
            \hline
            z \mspace{-6mu}&\mspace{-6mu} = \mspace{-6mu}&\mspace{-6mu} 5000 \mspace{-6mu}&\mspace{-6mu} - \mspace{-6mu}&\mspace{-6mu} \frac{10}{3} x_{4} \mspace{-6mu}&\mspace{-6mu} + \mspace{-6mu}&\mspace{-6mu} \frac{5}{3} x_{2}\\
            \hline
            \end{array}

        This is how the above output looks when rendered:

        .. MATH::

            \renewcommand{\arraystretch}{1.5}
            \begin{array}{|rcrcrcr|}
            \hline
            x_{3} \mspace{-6mu}&\mspace{-6mu} = \mspace{-6mu}&\mspace{-6mu} 1000 \mspace{-6mu}&\mspace{-6mu} - \mspace{-6mu}&\color{green}\mspace{-6mu} x_{1} \mspace{-6mu}&\mspace{-6mu} - \mspace{-6mu}&\mspace{-6mu} x_{2}\\
            \color{red}x_{4} \mspace{-6mu}&\color{red}\mspace{-6mu} = \mspace{-6mu}&\color{red}\mspace{-6mu} 1500 \mspace{-6mu}&\color{red}\mspace{-6mu} - \mspace{-6mu}&\color{blue}\mspace{-6mu} 3 x_{1} \mspace{-6mu}&\color{red}\mspace{-6mu} - \mspace{-6mu}&\color{red}\mspace{-6mu} x_{2}\\
            \hline
            z \mspace{-6mu}&\mspace{-6mu} = \mspace{-6mu}&\mspace{-6mu} 0 \mspace{-6mu}&\mspace{-6mu} + \mspace{-6mu}&\color{green}\mspace{-6mu} 10 x_{1} \mspace{-6mu}&\mspace{-6mu} + \mspace{-6mu}&\mspace{-6mu} 5 x_{2}\\
            \hline
            \\
            \hline
            x_{3} \mspace{-6mu}&\mspace{-6mu} = \mspace{-6mu}&\mspace{-6mu} 500 \mspace{-6mu}&\mspace{-6mu} + \mspace{-6mu}&\mspace{-6mu} \frac{1}{3} x_{4} \mspace{-6mu}&\mspace{-6mu} - \mspace{-6mu}&\mspace{-6mu} \frac{2}{3} x_{2}\\
            x_{1} \mspace{-6mu}&\mspace{-6mu} = \mspace{-6mu}&\mspace{-6mu} 500 \mspace{-6mu}&\mspace{-6mu} - \mspace{-6mu}&\mspace{-6mu} \frac{1}{3} x_{4} \mspace{-6mu}&\mspace{-6mu} - \mspace{-6mu}&\mspace{-6mu} \frac{1}{3} x_{2}\\
            \hline
            z \mspace{-6mu}&\mspace{-6mu} = \mspace{-6mu}&\mspace{-6mu} 5000 \mspace{-6mu}&\mspace{-6mu} - \mspace{-6mu}&\mspace{-6mu} \frac{10}{3} x_{4} \mspace{-6mu}&\mspace{-6mu} + \mspace{-6mu}&\mspace{-6mu} \frac{5}{3} x_{2}\\
            \hline
            \end{array}

        The column of the entering variable is green, while the row of the
        leaving variable is red in the original dictionary state on the top.
        The new state after the update step is shown on the bottom.
        """
        self.enter(entering)
        self.leave(leaving)
        result = latex(self).rsplit("\n", 1)[0] # Remove \end{array}
        # Make an empty line in the array
        if generate_real_LaTeX:
            result += "\n" r"\multicolumn{2}{c}{}\\[-3ex]" "\n"
        else:
            result += "\n\\\\\n"
        self.update()
        result += latex(self).split("\n", 2)[2] # Remove array header
        return LatexExpr(result)

    def basic_variables(self):
        r"""
        Return the basic variables of ``self``.

        OUTPUT:

        - a vector

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: D = P.initial_dictionary()
            sage: D.basic_variables()
            (x3, x4)
        """
        return self._AbcvBNz[4]

    def constant_terms(self):
        r"""
        Return the constant terms of relations of ``self``.

        OUTPUT:

        - a vector.

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: D = P.initial_dictionary()
            sage: D.constant_terms()
            (1000, 1500)
        """
        return self._AbcvBNz[1]

    def entering_coefficients(self):
        r"""
        Return coefficients of the entering variable.

        OUTPUT:

        - a vector

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: D = P.initial_dictionary()
            sage: D.enter(1)
            sage: D.entering_coefficients()
            (1, 3)
        """
        if self._entering is None:
            raise ValueError("entering variable must be chosen to compute "
                             "its coefficients")
        k = tuple(self.nonbasic_variables()).index(self._entering)
        return self._AbcvBNz[0].column(k)

    def leaving_coefficients(self):
        r"""
        Return coefficients of the leaving variable.

        OUTPUT:

        - a vector

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: D = P.dictionary(2, 3)
            sage: D.leave(3)
            sage: D.leaving_coefficients()
            (-2, -1)
        """
        if self._leaving is None:
            raise ValueError("leaving variable must be chosen to compute "
                             "its coefficients")
        i = tuple(self.basic_variables()).index(self._leaving)
        return self._AbcvBNz[0][i]

    def nonbasic_variables(self):
        r"""
        Return non-basic variables of ``self``.

        OUTPUT:

        - a vector

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: D = P.initial_dictionary()
            sage: D.nonbasic_variables()
            (x1, x2)
        """
        return self._AbcvBNz[5]

    def objective_coefficients(self):
        r"""
        Return coefficients of the objective of ``self``.

        OUTPUT:

        - a vector

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: D = P.initial_dictionary()
            sage: D.objective_coefficients()
            (10, 5)
        """
        return self._AbcvBNz[2]

    def objective_value(self):
        r"""
        Return the value of the objective at the
        :meth:`~LPAbstractDictionary.basic_solution` of ``self``.

        OUTPUT:

        - a number

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: D = P.initial_dictionary()
            sage: D.objective_value()
            0
        """
        return self._AbcvBNz[3]

    def update(self):
        r"""
        Update ``self`` using previously set entering and leaving variables.

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: D = P.initial_dictionary()
            sage: D.objective_value()
            0
            sage: D.enter("x1")
            sage: D.leave("x4")
            sage: D.update()
            sage: D.objective_value()
            5000
        """
        A, b, c, v, B, N, z = self._AbcvBNz
        entering = self._entering
        if entering is None:
            raise ValueError("entering variable must be set before updating")
        leaving = self._leaving
        if leaving is None:
            raise ValueError("leaving variable must be set before updating")
        l = tuple(B).index(leaving)
        e = tuple(N).index(entering)
        Ale = A[l, e]
        if Ale == 0:
            raise ValueError("incompatible choice of entering and leaving "
                             "variables")
        # Variables
        B[l] = entering
        N[e] = leaving
        # "The Changing Relation"
        b[l] /= Ale
        A[l] /= Ale
        A[l, e] = 1 / Ale
        # Other relations
        for i in range(A.nrows()):
            if i != l:
                Aie = A[i, e]
                A[i, e] = 0
                b[i] -= Aie * b[l]
                A[i] -= Aie * A[l]
        # Objective
        ce = c[e]
        c[e] = 0
        self._AbcvBNz[2] = c - ce * A[l]
        self._AbcvBNz[3] += ce * b[l]
        self._entering = None
        self._leaving = None


def random_dictionary(m, n, bound=5, special_probability=0.2):
    r"""
    Construct a random dictionary.

    INPUT:

    - ``m`` -- the number of constraints/basic variables

    - ``n`` -- the number of decision/non-basic variables

    - ``bound`` -- (default: 5) a bound on dictionary entries

    - ``special_probability`` -- (default: 0.2) probability of constructing a
      potentially infeasible or potentially optimal dictionary

    OUTPUT:

    - an :class:`LP problem dictionary <LPDictionary>`

    EXAMPLES::

        sage: from sage.numerical.interactive_simplex_method \
        ....:     import random_dictionary
        sage: random_dictionary(3, 4)
        LP problem dictionary (use typeset mode to see details)
    """
    A = random_matrix(ZZ, m, n, x=-bound, y=bound).change_ring(QQ)
    if special_probability < random():
        b = random_vector(ZZ, m, x=0, y=bound).change_ring(QQ)
    else:   # Allow infeasible dictionary
        b = random_vector(ZZ, m, x=-bound, y=bound).change_ring(QQ)
    if special_probability < random():
        c = random_vector(ZZ, n, x=-bound, y=bound).change_ring(QQ)
    else:   # Make dual feasible dictionary
        c = random_vector(ZZ, n, x=-bound, y=0).change_ring(QQ)
    x_N = list(PolynomialRing(QQ, "x", m + n + 1, order="neglex").gens())
    x_N.pop(0)
    x_B = []
    for i in range(m):
        x_B.append(x_N.pop(randint(0, n + m - i - 1)))
    return LPDictionary(A, b, c, randint(-bound, bound), x_B, x_N, "z")


class LPRevisedDictionary(LPAbstractDictionary):
    r"""
    Construct a revised dictionary for an LP problem.

    INPUT:

    - ``problem`` -- an :class:`LP problem in standard form
      <InteractiveLPProblemStandardForm>`

    - ``basic_variables`` -- a list of basic variables or their indices

    OUTPUT:

    - a :class:`revised dictionary for an LP problem <LPRevisedDictionary>`

    A revised dictionary encodes the same relations as a
    :class:`regular dictionary <LPDictionary>`, but stores only what is
    "necessary to efficiently compute data for the simplex method".

    Let the original problem be

    .. MATH::

        \begin{array}{l}
        \pm \max cx \\
        Ax \leq b \\
        x \geq 0
        \end{array}

    Let `\bar{x}` be the vector of :meth:`~InteractiveLPProblem.decision_variables` `x`
    followed by the :meth:`~InteractiveLPProblemStandardForm.slack_variables`.
    Let `\bar{c}` be the vector of :meth:`~InteractiveLPProblem.objective_coefficients` `c`
    followed by zeroes for all slack variables.
    Let `\bar{A} = (A | I)` be the matrix of
    :meth:`~InteractiveLPProblem.constraint_coefficients` `A` augmented by the identity
    matrix as columns corresponding to the slack variables. Then the problem
    above can be written as

    .. MATH::

        \begin{array}{l}
        \pm \max \bar{c} \bar{x} \\
        \bar{A} \bar{x} = b \\
        \bar{x} \geq 0
        \end{array}

    and any dictionary is a system of equations equivalent to
    `\bar{A} \bar{x} = b`, but resolved for :meth:`basic_variables` `x_B` in
    terms of :meth:`nonbasic_variables` `x_N` together with the expression for
    the objective in terms of `x_N`. Let :meth:`c_B` and :meth:`c_N` be vectors
    "splitting `\bar{c}` into basic and non-basic parts". Let :meth:`B` and
    :meth:`A_N` be the splitting of `\bar{A}`. Then the corresponding dictionary
    is

    .. MATH::

        \begin{array}{|l|}
        \hline
        x_B = B^{-1} b - B^{-1} A_N x_N\\
        \hline
        z = y b + \left(c_N - y^T A_N\right) x_N\\
        \hline
        \end{array}

    where `y = c_B^T B^{-1}`. To proceed with the simplex method, it is not
    necessary to compute all entries of this dictionary. On the other hand, any
    entry is easy to compute, if you know `B^{-1}`, so we keep track of it
    through the update steps.

    EXAMPLES::

        sage: A = ([1, 1], [3, 1])
        sage: b = (1000, 1500)
        sage: c = (10, 5)
        sage: P = InteractiveLPProblemStandardForm(A, b, c)
        sage: from sage.numerical.interactive_simplex_method \
        ....:     import LPRevisedDictionary
        sage: D = LPRevisedDictionary(P, [1, 2])
        sage: D.basic_variables()
        (x1, x2)
        sage: D
        LP problem dictionary (use typeset mode to see details)

    The same dictionary can be constructed through the problem::

        sage: P.revised_dictionary(1, 2) == D
        True

    When this dictionary is typeset, you will see two tables like these ones:

    .. MATH::

        \renewcommand{\arraystretch}{1.500000}
        \begin{array}{l}
        \begin{array}{l|r|rr||r||r}
        x_B & c_B &  & \mspace{-16mu} B^{-1} & y & B^{-1} b \\
        \hline
        x_{1} & 10 & -\frac{1}{2} & \frac{1}{2} & \frac{5}{2} & 250 \\
        x_{2} & 5 & \frac{3}{2} & -\frac{1}{2} & \frac{5}{2} & 750 \\
        \end{array}\\
        \\
        \begin{array}{r|rr}
        x_N & x_{3} & x_{4} \\
        \hline
        c_N^T & 0 & 0 \\
        \hline
        y^T A_N & \frac{5}{2} & \frac{5}{2} \\
        \hline
        c_N^T - y^T A_N & -\frac{5}{2} & -\frac{5}{2} \\
        \end{array}
        \end{array}

    More details will be shown if entering and leaving variables are set, but in
    any case the top table shows `B^{-1}` and a few extra columns, while the
    bottom one shows several rows: these are related to columns and rows of
    dictionary entries.
    """

    def __init__(self, problem, basic_variables):
        r"""
        See :class:`LPRevisedDictionary` for documentation.

        TESTS:::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: from sage.numerical.interactive_simplex_method \
            ....:     import LPRevisedDictionary
            sage: D = LPRevisedDictionary(P, [1, 2])
            sage: TestSuite(D).run()
        """
        if problem.auxiliary_variable() == problem.decision_variables()[0]:
            raise ValueError("revised dictionaries should not be constructed "
                             "for auxiliary problems")
        super(LPRevisedDictionary, self).__init__()
        self._problem = problem
        R =  problem.coordinate_ring()
        self._x_B = vector(R, [variable(R, v) for v in basic_variables])

    def __eq__(self, other):
        r"""
        Check if two revised LP problem dictionaries are equal.

        INPUT:

        - ``other`` -- anything

        OUTPUT:

        - ``True`` if ``other`` is an :class:`LPRevisedDictionary` for the same
          :class:`InteractiveLPProblemStandardForm` with the same :meth:`basic_variables`,
          ``False`` otherwise

        TESTS::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: from sage.numerical.interactive_simplex_method \
            ....:     import LPRevisedDictionary
            sage: D1 = LPRevisedDictionary(P, [1, 2])
            sage: D2 = LPRevisedDictionary(P, [1, 2])
            sage: D1 is D2
            False
            sage: D1 == D2
            True
            sage: D3 = LPRevisedDictionary(P, [2, 0])
            sage: D1 == D3
            False
        """
        return (isinstance(other, LPRevisedDictionary) and
                self._problem == other._problem and
                self._x_B == other._x_B)

    def _latex_(self):
        r"""
        Return a LaTeX representation of ``self``.

        OUTPUT:

        - a string

        TESTS::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: D = P.revised_dictionary()
            sage: D.enter(1)
            sage: D.leave(3)
            sage: print D._latex_()
            %notruncate
            \renewcommand{\arraystretch}{1.500000}
            \begin{array}{l}
            \begin{array}{l|r|rr||r||r|r|r}
            x_B & c_B &  & \mspace{-16mu} B^{-1} & y & B^{-1} b & B^{-1} A_{x_{1}} & \hbox{Ratio} \\
            \hline
            \color{red} x_{3} & \color{red} 0 & \color{red} 1 & \color{red} 0 & 0 & \color{red} 1000 & \color{red} 1 & \color{red} 1000 \\
            x_{4} & 0 & 0 & 1 & 0 & 1500 & 3 & 500 \\
            \end{array}\\
            \\
            \begin{array}{r|rr}
            x_N & \color{green} x_{1} & x_{2} \\
            \hline
            c_N^T & \color{green} 10 & 5 \\
            \hline
            y^T A_N & \color{green} 0 & 0 \\
            \hline
            c_N^T - y^T A_N & \color{green} 10 & 5 \\
            \end{array}
            \end{array}
        """
        latex.add_package_to_preamble_if_available("color")
        x_B = self._x_B
        m = len(x_B)
        entering = self._entering
        leaving = self._leaving
        show_ratios = entering is not None and self.is_feasible()
        if leaving is not None:
            l = x_B.list().index(leaving)
        lines = []
        lines.append(r"\begin{array}{l|r|%s||r||r%s%s}" % ("r"*m,
            "|r" if entering is not None else "", "|r" if show_ratios else ""))
        headers = ["x_B", "c_B"]
        if generate_real_LaTeX:
            headers.append(r"\multicolumn{%d}{c||}{B^{-1}}" % m)
        else:
            headers.extend([""] * (m//2))
            headers.append(r"\mspace{-16mu} B^{-1}")
            headers.extend([""] * ((m-1)//2))
        headers.extend(["y", "B^{-1} b"])
        if entering is not None:
            headers.append("B^{-1} A_{%s}" % latex(entering))
        if show_ratios:
            headers.append(r"\hbox{Ratio}")
        lines.append(" & ".join(headers) +  r" \\")
        lines.append(r"\hline")
        Bi = self.B_inverse()
        c_B = self.c_B()
        y = self.y()
        Bib = self.constant_terms()
        if entering is not None:
            Biae = self.entering_coefficients()
        if show_ratios:
            ratios = self.ratios()
        for i in range(m):
            entries = [x_B[i], c_B[i]]
            entries.extend(Bi.row(i))
            entries.extend([y[i], Bib[i]])
            if entering is not None:
                entries.append(Biae[i])
            if show_ratios:
                if ratios and ratios[0][1] == x_B[i]:
                    entries.append(ratios.pop(0)[0])
            terms = map(latex, entries)
            if leaving is not None and i == l:
                for j, t in enumerate(terms):
                    if j == m + 2:
                        continue
                    terms[j] = r"\color{red} " + t
            lines.append(" & ".join(terms) + r" \\")
        lines.append(r"\end{array}")
        top = "\n".join(lines)

        def make_line(header, terms):
            terms = map(latex, terms)
            if entering is not None:
                terms[k] = r"\color{green} " + terms[k]
            lines.append(" & ".join([header] + terms) + r" \\")

        lines = []
        x_N = self.x_N()
        if entering is not None:
            k = x_N.list().index(entering)
        lines.append(r"\begin{array}{r|" + "r" * len(x_N) + "}")
        make_line("x_N", x_N)
        lines.append(r"\hline")
        make_line("c_N^T", self.c_N())
        lines.append(r"\hline")
        make_line("y^T A_N", y * self.A_N())
        lines.append(r"\hline")
        make_line("c_N^T - y^T A_N", self.objective_coefficients())
        if leaving is not None and self.is_dual_feasible():
            lines.append(r"\hline")
            make_line("B^{-1}_{%s} A_N" % latex(leaving),
                      self.leaving_coefficients())
            lines.append(r"\hline")
            ratios = self.dual_ratios()
            make_line(r"\hbox{Ratio}", [ratios.pop(0)[0]
                                        if ratios and ratios[0][1] == x else ""
                                        for x in x_N])
        lines.append(r"\end{array}")
        bottom = "\n".join(lines)
        return _assemble_arrayl([top, "", bottom], 1.5)

    def A(self, v):
        r"""
        Return the column of constraint coefficients corresponding to ``v``.

        INPUT:

        - ``v`` -- a variable, its name, or its index

        OUTPUT:

        - a vector

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: D = P.revised_dictionary()
            sage: D.A(1)
            (1, 3)
            sage: D.A(0)
            (-1, -1)
            sage: D.A("x3")
            (1, 0)
        """
        P = self.problem()
        R = P.coordinate_ring()
        v = variable(R, v)
        k = R.gens().index(v)
        R = R.base_ring()
        m, n = P.m(), P.n()
        if k == 0:
            return vector(R, [-1] * m)
        elif k <= n:
            return P.A().column(k - 1)
        else:
            return identity_matrix(R, m).column(k - n - 1)

    def A_N(self):
        r"""
        Return the `A_N` matrix, constraint coefficients of
        non-basic variables.

        OUTPUT:

        - a matrix

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: D = P.revised_dictionary()
            sage: D.A_N()
            [1 1]
            [3 1]
        """
        return column_matrix(self.problem().base_ring(),
                             [self.A(x) for x in self.x_N()])

    def B(self):
        r"""
        Return the `B` matrix, i.e. constraint coefficients of
        basic variables.

        OUTPUT:

        - a matrix

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: D = P.revised_dictionary(1, 2)
            sage: D.B()
            [1 1]
            [3 1]
        """
        return column_matrix(self.problem().base_ring(),
                             [self.A(x) for x in self._x_B])

    def B_inverse(self):
        r"""
        Return the inverse of the :meth:`B` matrix.

        This inverse matrix is stored and computed during dictionary update in
        a more efficient way than generic inversion.

        OUTPUT:

        - a matrix

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: D = P.revised_dictionary(1, 2)
            sage: D.B_inverse()
            [-1/2  1/2]
            [ 3/2 -1/2]
        """
        try:
            return self._B_inverse
        except AttributeError:
            self._B_inverse = self.B().inverse()
            return self._B_inverse

    def E(self):
        r"""
        Return the eta matrix between ``self`` and the next dictionary.

        OUTPUT:

        - a matrix

        If `B_{\mathrm{old}}` is the current matrix `B` and `B_{\mathrm{new}}`
        is the `B` matrix of the next dictionary (after the update step), then
        `B_{\mathrm{new}} = B_{\mathrm{old}} E`.

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: D = P.revised_dictionary()
            sage: D.enter(1)
            sage: D.leave(4)
            sage: D.E()
            [1 1]
            [0 3]
        """
        if self._entering is None:
            raise ValueError("entering variable must be set to compute the "
                             "eta matrix")
        leaving = self._leaving
        if leaving is None:
            raise ValueError("leaving variable must be set to compute the "
                             "eta matrix")
        l = self._x_B.list().index(leaving)
        E = identity_matrix(self.base_ring(), self.problem().m())
        E.set_column(l, self.entering_coefficients())
        return E

    def E_inverse(self):
        r"""
        Return the inverse of the matrix :meth:`E`.

        This inverse matrix is computed in a more efficient way than generic
        inversion.

        OUTPUT:

        - a matrix

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: D = P.revised_dictionary()
            sage: D.enter(1)
            sage: D.leave(4)
            sage: D.E_inverse()
            [   1 -1/3]
            [   0  1/3]
        """
        E = self.E()
        l = self._x_B.list().index(self._leaving)
        d = E[l, l]
        if d == 0:
            raise ValueError("eta matrix is not invertible due to incompatible "
                             "choice of entering and leaving variables")
        E.set_col_to_multiple_of_col(l, l, -1/d)
        E[l, l] = 1 / d
        return E

    def basic_indices(self):
        r"""
        Return the basic indices of ``self``.

        .. NOTE::

            Basic indices are indices of :meth:`basic_variables` in the list of
            generators of the :meth:`~InteractiveLPProblemStandardForm.coordinate_ring` of
            the :meth:`problem` of ``self``, they may not coincide with the
            indices of variables which are parts of their names. (They will for
            the default indexed names.)

        OUTPUT:

        - a list.

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: D = P.revised_dictionary()
            sage: D.basic_indices()
            [3, 4]
        """
        gens = self.coordinate_ring().gens()
        return [gens.index(x) for x in self._x_B]

    def basic_variables(self):
        r"""
        Return the basic variables of ``self``.

        OUTPUT:

        - a vector

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: D = P.revised_dictionary()
            sage: D.basic_variables()
            (x3, x4)
        """
        return vector(self._x_B[0].parent(), self._x_B)

    def c_B(self):
        r"""
        Return the `c_B` vector, objective coefficients of basic variables.

        OUTPUT:

        - a vector

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: D = P.revised_dictionary(1, 2)
            sage: D.c_B()
            (10, 5)
        """
        P = self.problem()
        R = self.base_ring()
        BB = self.basic_indices()
        if 0 in BB:
            c_B = vector(R, P.m())
            c_B[BB.index(0)] = -1
            return c_B
        else:
            c_D = P.c()
            n = P.n()
            return vector(R, [c_D[k - 1] if k <= n else 0 for k in BB])

    def c_N(self):
        r"""
        Return the `c_N` vector, objective coefficients of non-basic variables.

        OUTPUT:

        - a vector

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: D = P.revised_dictionary()
            sage: D.c_N()
            (10, 5)
        """
        P = self.problem()
        n = P.n()
        R = P.base_ring()
        if 0 in self.basic_indices():
            return vector(R, n + 1)
        else:
            c_D = P.c()
            return vector(R, (c_D[k - 1] if k <= n else 0
                              for k in self.nonbasic_indices()))

    def constant_terms(self):
        r"""
        Return constant terms in the relations of ``self``.

        OUTPUT:

        - a vector

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: D = P.revised_dictionary()
            sage: D.constant_terms()
            (1000, 1500)
        """
        return self.B_inverse() * self.problem().b()

    def dictionary(self):
        r"""
        Return a regular LP dictionary matching ``self``.

        OUTPUT:

        - an :class:`LP dictionary <LPDictionary>`

        EXAMPLES::

            sage: A = ([1, 1], [3, 1], [-1, -1])
            sage: b = (1000, 1500, -400)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: D = P.revised_dictionary()
            sage: D.dictionary()
            LP problem dictionary (use typeset mode to see details)
        """
        D = LPDictionary(self.B_inverse() * self.A_N(),
                         self.constant_terms(),
                         self.objective_coefficients(),
                         self.objective_value(),
                         self.basic_variables(),
                         self.nonbasic_variables(),
                         "z")
        D._entering = self._entering
        D._leaving = self._leaving
        return D

    def entering_coefficients(self):
        r"""
        Return coefficients of the entering variable.

        OUTPUT:

        - a vector

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: D = P.revised_dictionary()
            sage: D.enter(1)
            sage: D.entering_coefficients()
            (1, 3)
        """
        if self._entering is None:
            raise ValueError("entering variable must be chosen to compute "
                             "its coefficients")
        return self.B_inverse() * self.A(self._entering)

    def leaving_coefficients(self):
        r"""
        Return coefficients of the leaving variable.

        OUTPUT:

        - a vector

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: D = P.revised_dictionary(2, 3)
            sage: D.leave(3)
            sage: D.leaving_coefficients()
            (-2, -1)
        """
        if self._leaving is None:
            raise ValueError("leaving variable must be chosen to compute "
                             "its coefficients")
        i = self.basic_variables().list().index(self._leaving)
        return self.B_inverse()[i] * self.A_N()

    def nonbasic_indices(self):
        r"""
        Return the non-basic indices of ``self``.

        .. NOTE::

            Non-basic indices are indices of :meth:`nonbasic_variables` in the
            list of generators of the
            :meth:`~InteractiveLPProblemStandardForm.coordinate_ring` of the
            :meth:`problem` of ``self``, they may not coincide with the indices
            of variables which are parts of their names. (They will for the
            default indexed names.)

        OUTPUT:

        - a list

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: D = P.revised_dictionary()
            sage: D.nonbasic_indices()
            [1, 2]
        """
        gens = self.coordinate_ring().gens()
        return [gens.index(x) for x in self.x_N()]

    def nonbasic_variables(self):
        r"""
        Return non-basic variables of ``self``.

        OUTPUT:

        - a vector

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: D = P.revised_dictionary()
            sage: D.nonbasic_variables()
            (x1, x2)
        """
        R = self.coordinate_ring()
        return vector(R, [xi for xi in R.gens()[1:] if xi not in self._x_B])

    def objective_coefficients(self):
        r"""
        Return coefficients of the objective of ``self``.

        OUTPUT:

        - a vector

        These are coefficients of non-basic variables when basic variables are
        eliminated.

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: D = P.revised_dictionary()
            sage: D.objective_coefficients()
            (10, 5)
        """
        return self.c_N() - self.y() * self.A_N()

    def objective_value(self):
        r"""
        Return the value of the objective at the basic solution of ``self``.

        OUTPUT:

        - a number

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: D = P.revised_dictionary()
            sage: D.objective_value()
            0
        """
        return self.y() * self.problem().b()

    def problem(self):
        r"""
        Return the original problem.

        OUTPUT:

        - an :class:`LP problem in standard form <InteractiveLPProblemStandardForm>`

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: D = P.revised_dictionary()
            sage: D.problem() is P
            True
        """
        return self._problem

    def update(self):
        r"""
        Update ``self`` using previously set entering and leaving variables.

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: D = P.revised_dictionary()
            sage: D.objective_value()
            0
            sage: D.enter("x1")
            sage: D.leave("x4")
            sage: D.update()
            sage: D.objective_value()
            5000
        """
        # Update the inverse of B first, in case it is impossible
        self._B_inverse = self.E_inverse() * self.B_inverse()
        # Now update the rest and clear settings
        self._x_B[self._x_B.list().index(self._leaving)] = self._entering
        self._entering = None
        self._leaving = None

    def y(self):
        r"""
        Return the `y` vector, the product of :meth:`c_B` and
        :meth:`B_inverse`.

        OUTPUT:

        - a vector

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: D = P.revised_dictionary()
            sage: D.y()
            (0, 0)
        """
        return self.c_B() * self.B_inverse()

    # Aliases for the standard notation
    x_B = basic_variables
    x_N = nonbasic_variables

# DEPRECATION (those two lines should be removed when cleaning #17867)
LPProblem = InteractiveLPProblem
LPProblemStandardForm = InteractiveLPProblemStandardForm
