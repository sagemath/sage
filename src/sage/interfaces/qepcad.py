r"""
Interface to QEPCAD
===================

The basic function of QEPCAD is to construct cylindrical algebraic
decompositions (CADs) of $\RR^k$, given a list of polynomials.  Using
this CAD, it is possible to perform quantifier elimination and formula
simplification.

A CAD for a set $A$ of $k$-variate polynomials decomposes $\RR^j$ into
disjoint cells, for each $j$ in $0 \leq j \leq k$.  The sign of each
polynomial in $A$ is constant in each cell of $\RR^k$, and for each
cell in $\RR^j$ ($j > 1$), the projection of that cell into
$\RR^{j-1}$ is a cell of $\RR^{j-1}$.  (This property makes the
decomposition 'cylindrical'.)

Given a formula $\exists x. P(a,b,x) = 0$ (for a polynomial $P$), and
a cylindrical algebraic decomposition for $P$, we can eliminate the
quantifier (find an equivalent formula in the two variables $a$, $b$
without the quantifier $\exists$) as follows.  For each cell $C$ in
$\RR^2$, find the cells of $\RR^3$ which project to $C$.  (This
collection is called the ``stack`` over $C$.)  Mark $C$ as true if
some member of the stack has sign $= 0$; otherwise, mark $C$ as false.
Then, construct a polynomial formula in $a$, $b$ which specifies
exactly the true cells (this is always possible).  The same technique
works if the body of the quantifier is any boolean combination of
polynomial equalities and inequalities.

Formula simplification is a similar technique.  Given a formula which
describes a simple set of $\RR^k$ in a complicated way as a boolean
combination of polynomial equalities and inequalities, QEPCAD can
construct a CAD for the polynomials and recover a simple equivalent
formula.

Note that while the following documentation is tutorial in nature, and
is written for people who may not be familiar with QEPCAD, it is
documentation for the \sage interface rather than for QEPCAD.  As
such, it does not cover several issues that are very important to use
QEPCAD efficiently, such as variable ordering, the efficient use of
the alternate quantifiers and ``_root_`` expressions, the
``measure-zero-error`` command, etc.  For more information on
QEPCAD, see the online documentation at
\url{http://www.cs.usna.edu/~qepcad/B/QEPCAD.html} and Chris Brown's
tutorial handout and slides from
\url{http://www.cs.usna.edu/~wcbrown/research/ISSAC04/Tutorial.html}.
(Several of the examples in this documentation came from these
sources.)

The examples below require that the optional qepcad package is installed.

QEPCAD can be run in a fully automatic fashion, or interactively.
We first demonstrate the automatic use of QEPCAD.

Since \sage has no built-in support for quantifiers, this interface
provides ``qepcad_formula`` which helps construct quantified
formulas in the syntax QEPCAD requires. ::

    sage: var('a,b,c,d,x,y,z')
    (a, b, c, d, x, y, z)
    sage: qf = qepcad_formula

We start with a simple example.  Consider an arbitrarily-selected
ellipse::

    sage: ellipse = 3*x^2 + 2*x*y + y^2 - x + y - 7

What is the projection onto the $x$ axis of this ellipse?  First we
construct a formula asking this question. ::

    sage: F = qf.exists(y, ellipse == 0); F
    (E y)[3 x^2 + 2 x y + y^2 - x + y - 7 = 0]

Then we run qepcad to get the answer::

    sage: qepcad(F)                            # optional - qepcad
    8 x^2 - 8 x - 29 <= 0

How about the projection onto the $y$ axis? ::

    sage: qepcad(qf.exists(x, ellipse == 0))   # optional - qepcad
    8 y^2 + 16 y - 85 <= 0

QEPCAD deals with more quantifiers than just 'exists', of course.
Besides the standard 'forall', there are also 'for infinitely
many', 'for all but finitely many', 'for a connected subset', and
'for exactly k'.  The :func:`qepcad` documentation has examples
of all of these; here we will just give one example.

First we construct a circle::

    sage: circle = x^2 + y^2 - 3

For what values $k$ does a vertical line $x=k$ intersect the combined
figure of the circle and ellipse exactly three times? ::

    sage: F = qf.exactly_k(3, y, circle * ellipse == 0); F
    (X3 y)[(3 x^2 + 2 x y + y^2 - x + y - 7) (x^2 + y^2 - 3) = 0]
    sage: qepcad(F)                                         # not tested (random order)
    x^2 - 3 <= 0 /\ 8 x^2 - 8 x - 29 <= 0 /\ 8 x^4 - 26 x^2 - 4 x + 13 >= 0 /\ [ 8 x^4 - 26 x^2 - 4 x + 13 = 0 \/ x^2 - 3 = 0 \/ 8 x^2 - 8 x - 29 = 0 ]

Here we see that the solutions are among the eight ($4 + 2 + 2$) roots
of the three polynomials inside the brackets, but not all of these
roots are solutions; the polynomial inequalities outside the brackets
are needed to select those roots that are solutions.

QEPCAD also supports an extended formula language, where
``_root_k`` `P(\bar x, y)` refers to a particular zero of
`P(\bar x, y)` (viewed as a polynomial in `y`).  If there are `n` roots,
then ``_root_1`` refers to the least root and ``_root_n`` refers to
the greatest. Also, ``_root_-n`` refers to the least root and
``_root_-1`` refers to the greatest.

This extended language is available both on input and output; see the
QEPCAD documentation for more information on how to use this syntax on
input.  We can request output that is intended to be easy to interpret
geometrically; then QEPCAD will use the extended language to produce
a solution formula without the selection polynomials. ::

    sage: qepcad(F, solution='geometric')                 # not tested (random order)
    x = _root_1 8 x^2 - 8 x - 29
    \/
    8 x^4 - 26 x^2 - 4 x + 13 = 0
    \/
    x = _root_-1 x^2 - 3

We then see that the 6 solutions correspond to the vertical tangent on
the left side of the ellipse, the four intersections between the
ellipse and the circle, and the vertical tangent on the right side of
the circle.

Let us do some basic formula simplification and visualization.
We will look at the region which is inside both the ellipse and the circle::

    sage: F = qf.and_(ellipse < 0, circle < 0); F
    [3 x^2 + 2 x y + y^2 - x + y - 7 < 0 /\ x^2 + y^2 - 3 < 0]
    sage: qepcad(F)                                       # not tested (random order)
    y^2 + 2 x y + y + 3 x^2 - x - 7 < 0 /\ y^2 + x^2 - 3 < 0

We get back the same formula we put in.  This is not surprising (we
started with a pretty simple formula, after all), but it is not very
enlightening either.  Again, if we ask for a 'geometric' output, then we see
an output that lets us understand something about the shape of the solution
set. ::

    sage: qepcad(F, solution='geometric')                 # not tested (random order)
    [
      [
        x = _root_-2 8 x^4 - 26 x^2 - 4 x + 13
        \/
        x = _root_2 8 x^4 - 26 x^2 - 4 x + 13
        \/
        8 x^4 - 26 x^2 - 4 x + 13 < 0
      ]
      /\
      y^2 + 2 x y + y + 3 x^2 - x - 7 < 0
      /\
      y^2 + x^2 - 3 < 0
    ]
    \/
    [
      x > _root_2 8 x^4 - 26 x^2 - 4 x + 13
      /\
      x < _root_-2 8 x^4 - 26 x^2 - 4 x + 13
      /\
      y^2 + x^2 - 3 < 0
    ]

There is another reason to prefer output using _root_ expressions; not
only does it sometimes give added insight into the geometric
structure, it also can be more efficient to construct.  Consider this
formula for the projection of a particular semicircle onto the $x$
axis::

    sage: F = qf.exists(y, qf.and_(circle == 0, x + y > 0)); F
    (E y)[x^2 + y^2 - 3 = 0 /\ x + y > 0]
    sage: qepcad(F)                                             # not tested (random order)
    x^2 - 3 <= 0 /\ [ x > 0 \/ 2 x^2 - 3 < 0 ]

Here, the formula $x > 0$ had to be introduced in order to get a
solution formula; the original CAD of F did not include the
polynomial $x$.  To avoid having QEPCAD do the extra work to come up
with a solution formula, we can tell it to use the extended language;
it is always possible to construct a solution formula in the extended
language without introducing new polynomials. ::

    sage: qepcad(F, solution='extended')                        # not tested (random order)
    x^2 - 3 <= 0 /\ x > _root_1 2 x^2 - 3

Up to this point, all the output we have seen has basically been in the
form of strings; there is no support (yet) for parsing these outputs
back into \sage polynomials (partly because \sage does not yet have
support for symbolic conjunctions and disjunctions).  The function
:func:`qepcad` supports three more output types that give numbers which
can be manipulated in \sage: any-point, all-points, and cell-points.

These output types give dictionaries mapping variable names to values.
With any-point, :func:`qepcad` either produces a single dictionary
specifying a point where the formula is true, or raises an exception
if the formula is false everywhere.  With all-points,
:func:`qepcad` either produces a list of dictionaries for all
points where the formula is true, or raises an exception if the
formula is true on infinitely many points.  With cell-points,
:func:`qepcad` produces a list of dictionaries with one point for
each cell where the formula is true.  (This means you will have at least
one point in each connected component of the solution, although you will
often have many more points than that.)

Let us revisit some of the above examples and get some points to play
with.  We will start by finding a point on our ellipse. ::

    sage: p = qepcad(ellipse == 0, solution='any-point'); p  # optional - qepcad
    {'x': -1.468501968502953?, 'y': 0.9685019685029527?}

(Note that despite the decimal printing and the question marks, these
are really exact numbers.)

We can verify that this point is a solution.  To do so, we create
a copy of ellipse as a polynomial over `\QQ` (instead of a symbolic
expression). ::

    sage: pellipse = QQ['x,y'](ellipse)
    sage: pellipse(**p) == 0                           # optional - qepcad
    True

For cell-points, let us look at points *not* on the ellipse. ::

    sage: pts = qepcad(ellipse != 0, solution='cell-points'); pts   # optional - qepcad
    [{'x': 4, 'y': 0},
     {'x': 2.468501968502953?, 'y': 1},
     {'x': 2.468501968502953?, 'y': -9},
     {'x': 1/2, 'y': 9},
     {'x': 1/2, 'y': -1},
     {'x': 1/2, 'y': -5},
     {'x': -1.468501968502953?, 'y': 3},
     {'x': -1.468501968502953?, 'y': -1},
     {'x': -3, 'y': 0}]

For the points here which are in full-dimensional cells, QEPCAD has the
freedom to choose rational sample points, and it does so.

And, of course, all these points really are not on the ellipse. ::

    sage: [pellipse(**p) != 0 for p in pts]               # optional - qepcad
    [True, True, True, True, True, True, True, True, True]

Finally, for all-points, let us look again at finding vertical lines that
intersect the union of the circle and the ellipse exactly three times. ::

    sage: F = qf.exactly_k(3, y, circle * ellipse == 0); F
    (X3 y)[(3 x^2 + 2 x y + y^2 - x + y - 7) (x^2 + y^2 - 3) = 0]
    sage: pts = qepcad(F, solution='all-points'); pts       # optional - qepcad
    [{'x': 1.732050807568878?}, {'x': 1.731054913462534?}, {'x': 0.678911384208004?}, {'x': -0.9417727377417167?}, {'x': -1.468193559928821?}, {'x': -1.468501968502953?}]

Since $y$ is bound by the quantifier, the solutions only refer to $x$.

We can substitute one of these solutions into the original equation::

    sage: pt = pts[0]                             # optional - qepcad
    sage: pcombo = QQ['x,y'](circle * ellipse)
    sage: intersections = pcombo(y=polygen(AA, 'y'), **pt); intersections       # optional - qepcad
    y^4 + 4.464101615137755?*y^3 + 0.2679491924311227?*y^2

and verify that it does have three roots::

    sage: intersections.roots()                                                 # optional - qepcad
    [(-4.403249005600958?, 1), (-0.06085260953679653?, 1), (0, 2)]

Let us check all six solutions. ::

    sage: [len(pcombo(y=polygen(AA, 'y'), **p).roots()) for p in pts]    # optional - qepcad
    [3, 3, 3, 3, 3, 3]

We said earlier that we can run QEPCAD either automatically or
interactively.  Now that we have discussed the automatic modes, let us turn
to interactive uses.

If the :func:`qepcad` function is passed ``interact=True``, then
instead of returning a result, it returns an object of class
:class:`Qepcad` representing a running instance of QEPCAD that you can
interact with.  For example::

    sage: qe = qepcad(qf.forall(x, x^2 + b*x + c > 0), interact=True); qe     # optional - qepcad
    QEPCAD object in phase 'Before Normalization'

This object is a fairly thin wrapper over QEPCAD; most QEPCAD commands
are available as methods on the :class:`Qepcad` object.  Given a
:class:`Qepcad` object ``qe``, you can type ``qe.[tab]`` to see
the available QEPCAD commands; to see the documentation for an
individual QEPCAD command, for example ``d_setting``, you can type
``qe.d_setting?``.  (In QEPCAD, this command is called
``d-setting``. We systematically replace hyphens with underscores
for this interface.)

The execution of QEPCAD is divided into four phases. Most commands are
not available during all phases.  We saw above that QEPCAD starts out
in phase ``'Before Normalization'``. We see that the ``d_cell`` command
is not available in this phase::

    sage: qe.d_cell()                                # optional - qepcad
    Error GETCID: This command is not active here.

We will focus here on the fourth (and last) phase, ``'Before Solution'``,
because this interface has special support for some operations in this
phase.  Consult the QEPCAD documentation for information on the other
phases.

We can tell QEPCAD to finish off the current phase and move to the
next with its ``go`` command.  (There is also the ``step``
command, which partially completes a phase for phases that have
multiple steps, and the ``finish`` command, which runs QEPCAD to
completion.) ::

    sage: qe.go()                                        # optional - qepcad
    QEPCAD object has moved to phase 'Before Projection (x)'
    sage: qe.go()                                        # optional - qepcad
    QEPCAD object has moved to phase 'Before Choice'
    sage: qe.go()                                        # optional - qepcad
    QEPCAD object has moved to phase 'Before Solution'

Note that the :class:`Qepcad` object returns the new phase whenever the
phase changes, as a convenience for interactive use; except that when
the new phase is ``'EXITED'``, the solution formula printed by QEPCAD is
returned instead. ::

    sage: qe.go()                                        # optional - qepcad
    4 c - b^2 > 0
    sage: qe                                             # optional - qepcad
    QEPCAD object in phase 'EXITED'

Let us pick a nice, simple example, return to phase 4, and explore the
resulting ``qe`` object. ::

    sage: qe = qepcad(circle == 0, interact=True); qe       # optional - qepcad
    QEPCAD object in phase 'Before Normalization'
    sage: qe.go(); qe.go(); qe.go()                         # optional - qepcad
    QEPCAD object has moved to phase 'Before Projection (y)'
    QEPCAD object has moved to phase 'Before Choice'
    QEPCAD object has moved to phase 'Before Solution'

We said before that QEPCAD creates 'cylindrical algebraic
decompositions'; since we have a bivariate polynomial, we get
decompositions of $\RR^0$, $\RR^1$, and $\RR^2$.  In this case, where
our example is a circle of radius $\sqrt{3}$ centered on the origin,
these decompositions are as follows:

The decomposition of $\RR^0$ is trivial (of course).  The
decomposition of $\RR^1$ has five cells: $x < -\sqrt{3}$, $x =
-\sqrt{3}$, $-\sqrt{3} < x < \sqrt{3}$, $x = \sqrt{3}$, and $x >
\sqrt{3}$.  These five cells comprise the ``stack`` over the single
cell in the trivial decomposition of $\RR^0$.

These five cells give rise to five stacks in $\RR^2$.  The first and
fifth stack have just one cell apiece.  The second and fourth stacks
have three cells: $y < 0$, $y = 0$, and $y > 0$.  The third stack has
five cells: below the circle, the lower semicircle, the interior of
the circle, the upper semicircle, and above the circle.

QEPCAD (and this QEPCAD interface) number the cells in a stack
starting with 1.  Each cell has an ``index``, which is a tuple of
integers describing the path to the cell in the tree of all cells.
For example, the cell 'below the circle' has index (3,1) (the first
cell in the stack over the third cell of $\RR^1$) and the interior of
the circle has index (3,3).

We can view these cells with the QEPCAD command ``d_cell``.  For
instance, let us look at the cell for the upper semicircle::

    sage: qe.d_cell(3, 4)                              # optional - qepcad
    ---------- Information about the cell (3,4) ----------
    Level                       : 2
    Dimension                   : 1
    Number of children          : 0
    Truth value                 : T    by trial evaluation.
    Degrees after substitution  : Not known yet or No polynomial.
    Multiplicities              : ((1,1))
    Signs of Projection Factors
    Level 1  : (-)
    Level 2  : (0)
    ----------   Sample point  ----------
    The sample point is in a PRIMITIVE representation.
    <BLANKLINE>
    alpha = the unique root of x^2 - 3 between 0 and 4
          = 1.7320508076-
    <BLANKLINE>
    Coordinate 1 = 0
                 = 0.0000000000
    Coordinate 2 = alpha
                 = 1.7320508076-
    ----------------------------------------------------

We see that, the level of this cell is 2, meaning that it is part of
the decomposition of $\RR^2$.  The dimension is 1, meaning that the
cell is homeomorphic to a line (rather than a plane or a point).  The
sample point gives the coordinates of one point in the cell, both
symbolically and numerically.

For programmatic access to cells, we have defined a \sage wrapper class
:class:`QepcadCell`.  These cells can be created with the
:meth:`cell` method; for example::

    sage: c = qe.cell(3, 4); c                       # optional - qepcad
    QEPCAD cell (3, 4)

A :class:`QepcadCell` has accessor methods for the important state held
within a cell.  For instance::

    sage: c.level()                                  # optional - qepcad
    2
    sage: c.index()                                  # optional - qepcad
    (3, 4)
    sage: qe.cell(3).number_of_children()            # optional - qepcad
    5
    sage: len(qe.cell(3))                            # optional - qepcad
    5

One particularly useful thing we can get from a cell is its sample point,
as \sage algebraic real numbers. ::

    sage: c.sample_point()                           # optional - qepcad
    (0, 1.732050807568878?)
    sage: c.sample_point_dict()                      # optional - qepcad
    {'x': 0, 'y': 1.732050807568878?}

We have seen that we can get cells using the :meth:`cell` method.
There are several QEPCAD commands that print lists of cells; we can
also get cells using the :meth:`make_cells` method, passing it the
output of one of these commands. ::

    sage: qe.make_cells(qe.d_true_cells())           # optional - qepcad
    [QEPCAD cell (4, 2), QEPCAD cell (3, 4), QEPCAD cell (3, 2),
    QEPCAD cell (2, 2)]

Also, the cells in the stack over a given cell can be accessed using
array subscripting or iteration.  (Remember that cells in a stack are
numbered starting with one; we preserve this convention in the
array-subscripting syntax.) ::

    sage: c = qe.cell(3)                             # optional - qepcad
    sage: c[1]                                       # optional - qepcad
    QEPCAD cell (3, 1)
    sage: [c2 for c2 in c]                           # optional - qepcad
    [QEPCAD cell (3, 1), QEPCAD cell (3, 2), QEPCAD cell (3, 3),
    QEPCAD cell (3, 4), QEPCAD cell (3, 5)]

We can do one more thing with a cell: we can set its truth value.
Once the truth values of the cells have been set, we can get QEPCAD to
produce a formula which is true in exactly the cells we have selected.
This is useful if QEPCAD's quantifier language is insufficient to
express your problem.

For example, consider again our combined figure of the circle and the
ellipse.  Suppose you want to find all vertical lines that intersect
the circle twice, and also intersect the ellipse twice.  The vertical
lines that intersect the circle twice can be found by simplifying::

    sage: F = qf.exactly_k(2, y, circle == 0); F
    (X2 y)[x^2 + y^2 - 3 = 0]

and the vertical lines that intersect the ellipse twice are expressed by::

    sage: G = qf.exactly_k(2, y, ellipse == 0); G
    (X2 y)[3 x^2 + 2 x y + y^2 - x + y - 7 = 0]

and the lines that intersect both figures would be::

    sage: qf.and_(F, G)
    Traceback (most recent call last):
    ...
    ValueError: QEPCAD formulas must be in prenex (quantifiers outermost) form

...except that QEPCAD does not support formulas like this; in QEPCAD
input, all logical connectives must be inside all quantifiers.

Instead, we can get QEPCAD to construct a CAD for our combined figure
and set the truth values ourselves.  (The exact formula we use doesn't
matter, since we're going to replace the truth values in the cells; we
just need to use a formula that uses both polynomials.) ::

    sage: qe = qepcad(qf.and_(ellipse == 0, circle == 0), interact=True)       # optional - qepcad
    sage: qe.go(); qe.go(); qe.go()                         # optional - qepcad
    QEPCAD object has moved to phase 'Before Projection (y)'
    QEPCAD object has moved to phase 'Before Choice'
    QEPCAD object has moved to phase 'Before Solution'

Now we want to find all cells $c$ in the decomposition of $\RR^1$ such
that the stack over $c$ contains exactly two cells on the ellipse, and
also contains exactly two cells on the circle.

Our input polynomials are 'level-2 projection factors', we see::

    sage: qe.d_proj_factors()                               # optional - qepcad
    P_1,1  = fac(J_1,1) = fac(dis(A_2,1))
           = 8 x^2 - 8 x - 29
    P_1,2  = fac(J_1,2) = fac(dis(A_2,2))
           = x^2 - 3
    P_1,3  = fac(J_1,3) = fac(res(A_2,1|A_2,2))
           = 8 x^4 - 26 x^2 - 4 x + 13
    A_2,1  = input
           = y^2 + 2 x y + y + 3 x^2 - x - 7
    A_2,2  = input
           = y^2 + x^2 - 3

so we can test whether a cell is on the ellipse by checking that the
sign of the corresponding projection factor is 0 in our cell.
For instance, the cell (12,2) is on the ellipse::

    sage: qe.cell(12,2).signs()[1][0]                       # optional - qepcad
    0

So we can update the truth values as desired like this::

    sage: for c in qe.cell():                               # optional - qepcad
    ....:     count_ellipse = 0
    ....:     count_circle = 0
    ....:     for c2 in c:
    ....:         count_ellipse += (c2.signs()[1][0] == 0)
    ....:         count_circle += (c2.signs()[1][1] == 0)
    ....:     c.set_truth(count_ellipse == 2 and count_circle == 2)

and then we can get our desired solution formula.  (The ``'G'`` stands for
``'geometric'``, and gives solutions using the same rules as
``solution='geometric'`` described above.) ::

    sage: qe.solution_extension('G')                        # not tested (random order)
    8 x^2 - 8 x - 29 < 0
    /\
    x^2 - 3 < 0

TESTS:

Check the qepcad configuration file::

    sage: with open(os.path.join(SAGE_LOCAL, 'etc', 'default.qepcadrc')) as f:  # optional - qepcad
    ....:     f.readlines()[-1]
    'SINGULAR yes\n'

Tests related to the not tested examples (nondeterministic order of atoms)::

    sage: from sage.interfaces.qepcad import _qepcad_atoms
    sage: var('a,b,c,d,x,y,z')
    (a, b, c, d, x, y, z)
    sage: qf = qepcad_formula
    sage: ellipse = 3*x^2 + 2*x*y + y^2 - x + y - 7
    sage: circle = x^2 + y^2 - 3

    sage: F = qf.exactly_k(3, y, circle * ellipse == 0)
    sage: _qepcad_atoms(qepcad(F))                                 # optional - qepcad
    {'8 x^2 - 8 x - 29 <= 0',
     '8 x^2 - 8 x - 29 = 0',
     '8 x^4 - 26 x^2 - 4 x + 13 = 0',
     '8 x^4 - 26 x^2 - 4 x + 13 >= 0',
     'x^2 - 3 <= 0',
     'x^2 - 3 = 0'}
    sage: _qepcad_atoms(qepcad(F, solution='geometric'))           # optional - qepcad
    {'8 x^4 - 26 x^2 - 4 x + 13 = 0',
     'x = _root_-1 x^2 - 3',
     'x = _root_1 8 x^2 - 8 x - 29'}

    sage: F = qf.and_(ellipse < 0, circle < 0)
    sage: _qepcad_atoms(qepcad(F))                                 # optional - qepcad
    {'y^2 + 2 x y + y + 3 x^2 - x - 7 < 0', 'y^2 + x^2 - 3 < 0'}
    sage: _qepcad_atoms(qepcad(F, solution='geometric'))           # optional - qepcad
    {'8 x^4 - 26 x^2 - 4 x + 13 < 0',
     'x < _root_-2 8 x^4 - 26 x^2 - 4 x + 13',
     'x = _root_-2 8 x^4 - 26 x^2 - 4 x + 13',
     'x = _root_2 8 x^4 - 26 x^2 - 4 x + 13',
     'x > _root_2 8 x^4 - 26 x^2 - 4 x + 13',
     'y^2 + 2 x y + y + 3 x^2 - x - 7 < 0',
     'y^2 + x^2 - 3 < 0'}

    sage: F = qf.exists(y, qf.and_(circle == 0, x + y > 0))
    sage: _qepcad_atoms(qepcad(F))                                 # optional - qepcad
    {'2 x^2 - 3 < 0', 'x > 0', 'x^2 - 3 <= 0'}
    sage: _qepcad_atoms(qepcad(F, solution='extended'))            # optional - qepcad
    {'x > _root_1 2 x^2 - 3', 'x^2 - 3 <= 0'}

    sage: qe = qepcad(qf.and_(ellipse == 0, circle == 0), interact=True)       # optional - qepcad
    sage: qe.go(); qe.go(); qe.go()                         # optional - qepcad
    QEPCAD object has moved to phase 'Before Projection (y)'
    QEPCAD object has moved to phase 'Before Choice'
    QEPCAD object has moved to phase 'Before Solution'
    sage: for c in qe.cell():                               # optional - qepcad
    ....:     count_ellipse = 0
    ....:     count_circle = 0
    ....:     for c2 in c:
    ....:         count_ellipse += (c2.signs()[1][0] == 0)
    ....:         count_circle += (c2.signs()[1][1] == 0)
    ....:     c.set_truth(count_ellipse == 2 and count_circle == 2)
    sage: _qepcad_atoms(qe.solution_extension('G'))                # optional - qepcad
    {'8 x^2 - 8 x - 29 < 0', 'x^2 - 3 < 0'}


AUTHORS:

- Carl Witty (2008-03): initial version
- Thierry Monteil (2015-07) repackaging + noncommutative doctests.
"""

#*****************************************************************************
#       Copyright (C) 2008 Carl Witty <Carl.Witty@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.env import SAGE_LOCAL
import pexpect
import re
import sys

from sage.cpython.string import bytes_to_str
from sage.misc.flatten import flatten
from sage.misc.sage_eval import sage_eval
from sage.repl.preparse import implicit_mul
from sage.interfaces.tab_completion import ExtraTabCompletion
from sage.docs.instancedoc import instancedoc
from .expect import Expect, ExpectFunction
from sage.interfaces.interface import AsciiArtString


def _qepcad_atoms(formula):
    r"""
    Return the atoms of a qepcad quantifier-free formula, as a set of strings.

    INPUT:

    - `formula` (string) - a quantifier-free formula.

    .. note:: this function is pis-aller used for doctesting, not a complete
    parser, which should be written in a further ticket.

    EXAMPLES::

    sage: from sage.interfaces.qepcad import _qepcad_atoms
    sage: _qepcad_atoms('y^5 + 4 y + 8 >= 0 /\\ y <= 0 /\\ [ y = 0 \\/ y^5 + 4 y + 8 = 0 ]')
    {'y <= 0', 'y = 0', 'y^5 + 4 y + 8 = 0', 'y^5 + 4 y + 8 >= 0'}
    """
    return set(i.strip() for i in flatten([i.split('\\/') for i in formula.replace('[','').replace(']','').split('/\\')]))

def _qepcad_cmd(memcells=None):
    r"""
    Construct a QEPCAD command line.

    Optionally set the number of memory cells to use.

    EXAMPLES::

        sage: from sage.interfaces.qepcad import _qepcad_cmd
        sage: s = _qepcad_cmd()
        sage: s == 'env qe=%s qepcad '%SAGE_LOCAL
        True
        sage: s = _qepcad_cmd(memcells=8000000)
        sage: s == 'env qe=%s qepcad +N8000000'%SAGE_LOCAL
        True
    """
    if memcells is not None:
        memcells_arg = '+N%s' % memcells
    else:
        memcells_arg = ''
    return "env qe=%s qepcad %s"%(SAGE_LOCAL, memcells_arg)

_command_info_cache = None


def _update_command_info():
    r"""
    Read the file ``qepcad.help`` to find the list of commands
    supported by QEPCAD.

    Used for tab-completion and documentation.

    EXAMPLES::

        sage: from sage.interfaces.qepcad import _update_command_info, _command_info_cache
        sage: _update_command_info()                # optional - qepcad
        sage: _command_info_cache['approx_precision']        # optional - qepcad
        ('46', 'abcde', 'm', 'approx-precision N\n\nApproximate algeraic numbers to N decimal places.\n', None)
    """
    global _command_info_cache
    if _command_info_cache is not None:
        return

    cache = {}

    with open(os.path.join(SAGE_LOCAL, 'share/qepcad', 'qepcad.help')) as help:
        assert(help.readline().strip() == '@')

        while True:
            cmd_line = help.readline()
            while not cmd_line.strip():
                cmd_line = help.readline()
            cmd_line = cmd_line.strip()
            if cmd_line == '@@@':
                break

            (cmd, id, phases, kind) = cmd_line.split()
            assert(help.readline().strip() == '@')

            help_text = ''
            help_line = help.readline()
            while help_line.strip() != '@':
                help_text += help_line
                help_line = help.readline()

            # I went through qepcad.help and picked out commands that
            # I thought might be worth a little extra tweaking...
            special = None

            # These commands have been tweaked.
            if cmd in ['d-all-cells-in-subtree', 'd-cell', 'd-pcad',
                       'd-pscad', 'd-stack', 'manual-choose-cell']:
                special = 'cell'
            if cmd in ['ipfzt', 'rational-sample', 'triv-convert', 'use-db',
                       'use-selected-cells-cond', 'verbose']:
                special = 'yn'

            # The tweaking for these commands has not been implemented yet.
            if cmd in ['selected-cells-cond']:
                special = 'formula'
            if cmd in ['ch-pivot', 'rem-pf', 'rem-pj']:
                special = 'i,j'
            if cmd in ['d-2d-cad', 'set-truth-value', 'solution-extension']:
                special = 'interactive'
            if cmd in ['p-2d-cad', 'trace-alg', 'trace-data']:
                special = 'optional'

            cmd = cmd.replace('-', '_')

            cache[cmd] = (id, phases, kind, help_text, special)

    _command_info_cache = cache

# QEPCAD does not have a typical "computer algebra system" interaction
# model.  Instead, you run QEPCAD once for each problem you wish to solve,
# then interact with it while you solve that problem.

# For this reason, most of the methods on the Expect class are
# inapplicable.  To avoid confusing the user with useless commands
# in tab completion, we split the control into two classes:
# Qepcad_expect is a subclass of Expect, and is never seen by the user;
# Qepcad is a wrapper for Qepcad_expect, and is what the user interacts with.


class Qepcad_expect(ExtraTabCompletion, Expect):
    r"""
    The low-level wrapper for QEPCAD.
    """
    def __init__(self, memcells=None,
                 maxread=None,
                 logfile=None,
                 server=None):
        r"""
        Initialize a low-level wrapper for QEPCAD.

        You can specify
        the number of memory cells that QEPCAD allocates on startup
        (which controls the maximum problem size QEPCAD can handle),
        and specify a logfile.  You can also specify a server, and the
        interface will run QEPCAD on that server, using ssh.  (UNTESTED)

        EXAMPLES::

            sage: from sage.interfaces.qepcad import Qepcad_expect
            sage: Qepcad_expect(memcells=100000, logfile=sys.stdout)
            Qepcad
        """
        Expect.__init__(self,
                        name="QEPCAD",
                        # yuck: when QEPCAD first starts,
                        # it doesn't give prompts
                        prompt="\nEnter an .*:\r",
                        command=_qepcad_cmd(memcells),
                        server=server,
                        restart_on_ctrlc=False,
                        verbose_start=False,
                        logfile=logfile)


class Qepcad:
    r"""
    The wrapper for QEPCAD.
    """

    def __init__(self, formula,
                 vars=None, logfile=None, verbose=False,
                 memcells=None, server=None):
        r"""
        Construct a QEPCAD wrapper object.

        Requires a formula, which
        may be a :class:`qformula` as returned by the methods of
        ``qepcad_formula``, a symbolic equality or inequality, a
        polynomial $p$ (meaning $p = 0$), or a string, which is passed
        straight to QEPCAD.

        ``vars`` specifies the variables to use; this gives the variable
        ordering, which may be very important.  If ``formula`` is
        given as a string, then ``vars`` is required; otherwise,
        if ``vars`` is omitted, then a default ordering is used
        (alphabetical ordering for the free variables).

        A logfile can be specified with ``logfile``.
        If ``verbose=True`` is given, then the logfile is automatically
        set to ``sys.stdout``, so all QEPCAD interaction is echoed to
        the terminal.

        You can set the amount of memory that QEPCAD allocates with
        ``memcells``, and you can use ``server`` to run QEPCAD on
        another machine using ssh.  (UNTESTED)

        Usually you will not call this directly, but use ``qepcad``
        to do so.  Check the ``qepcad`` documentation for more
        information.

        EXAMPLES::

            sage: from sage.interfaces.qepcad import Qepcad
            sage: Qepcad(x^2 - 1 == 0)            # optional - qepcad
            QEPCAD object in phase 'Before Normalization'

        To check that :trac:`20126` is fixed::

            sage: (x, y, z) = var('x, y, z')
            sage: conds = [-z < 0, -y + z < 0, x^2 + x*y + 2*x*z + 2*y*z - x < 0, \
                           x^2 + x*y + 3*x*z + 2*y*z + 2*z^2 - x - z < 0, \
                           -2*x + 1 < 0, -x*y - x*z - 2*y*z - 2*z^2 + z < 0, \
                           x + 3*y + 3*z - 1 < 0]
            sage: qepcad(conds, memcells=2000000) # optional - qepcad
            2 x - 1 > 0 /\ z > 0 /\ z - y < 0 /\ 3 z + 3 y + x - 1 < 0
        """
        self._cell_cache = {}

        if verbose:
            logfile=sys.stdout

        varlist = None
        if vars is not None:
            if isinstance(vars, str):
                varlist = vars.strip('()').split(',')
            else:
                varlist = [str(v) for v in vars]

        if isinstance(formula, str):
            if varlist is None:
                raise ValueError("vars must be specified if formula is a string")

            free_vars = len(varlist) - formula.count('(')
        else:
            formula = qepcad_formula.formula(formula)
            fvars = formula.vars
            fqvars = formula.qvars
            if varlist is None:
                varlist = sorted(fvars) + fqvars
            else:
                # Do some error-checking here.  Parse the given vars
                # and ensure they match up with the variables in the formula.
                if frozenset(varlist) != (fvars | frozenset(fqvars)):
                    raise ValueError("specified vars don't match vars in formula")
                if len(fqvars) and varlist[-len(fqvars):] != fqvars:
                    raise ValueError("specified vars don't match quantified vars")
            free_vars = len(fvars)
            formula = repr(formula)

        self._varlist = varlist
        self._free_vars = free_vars

        varlist = [v.replace('_', '') for v in varlist]
        if len(frozenset(varlist)) != len(varlist):
            raise ValueError("variables collide after stripping underscores")
        formula = formula.replace('_', '')

        qex = Qepcad_expect(logfile=logfile, memcells=memcells, server=server)
        qex._send('[ input from Sage ]')
        qex._send('(' + ','.join(varlist) + ')')
        qex._send(str(free_vars))
        # I hope this prompt is distinctive enough...
        # if not, we could list all the cases separately
        qex._change_prompt(['\r\n([^\r\n]*) >\r\n', pexpect.EOF])
        qex.eval(formula + '.')
        self._qex = qex

    def __repr__(self):
        r"""
        Return a string representation of this :class:`Qepcad` object.

        EXAMPLES::

            sage: qepcad(x - 1 == 0, interact=True) # optional - qepcad
            QEPCAD object in phase 'Before Normalization'
        """
        return "QEPCAD object in phase '{}'".format(self.phase())

    def assume(self, assume):
        r"""
        The following documentation is from ``qepcad.help``.

        Add an assumption to the problem.  These will not be
        included in the solution formula.

        For example, with input  (E x)[ a x^2 + b x + c = 0],
        if we issue the command

            assume [ a /= 0 ]

        we will get the solution formula b^2 - 4 a c >= 0.  Without
        the assumption we'd get something like [a = 0 /\ b /= 0] \/
        [a /= 0 /\ 4 a c - b^2 <= 0] \/ [a = 0 /\ b = 0 /\ c = 0].

        EXAMPLES::

            sage: var('a,b,c,x')
            (a, b, c, x)
            sage: qf = qepcad_formula
            sage: qe = qepcad(qf.exists(x, a*x^2 + b*x + c == 0), interact=True) # optional - qepcad
            sage: qe.assume(a != 0) # optional - qepcad
            sage: qe.finish() # optional - qepcad
            4 a c - b^2 <= 0
        """
        if not isinstance(assume, str):
            assume = qepcad_formula.formula(assume)
            if len(assume.qvars):
                raise ValueError("assumptions cannot be quantified")
            if not assume.vars.issubset(frozenset(self._varlist[:self._free_vars])):
                raise ValueError("assumption contains variables not "
                                 "present in formula")
            assume = repr(assume)
        assume = assume.replace('_', '')
        result = self._eval_line("assume [%s]" % assume)
        if len(result):
            return AsciiArtString(result)

    def solution_extension(self, kind):
        r"""
        The following documentation is modified from ``qepcad.help``:

        solution-extension x

        Use an alternative solution formula construction method.  The
        parameter x is allowed to be T,E, or G.  If x is T, then a
        formula in the usual language of Tarski formulas is produced.
        If x is E, a formula in the language of Extended Tarski formulas
        is produced.  If x is G, then a geometry-based formula is
        produced.

        EXAMPLES::

            sage: var('x,y')
            (x, y)
            sage: qf = qepcad_formula
            sage: qe = qepcad(qf.and_(x^2 + y^2 - 3 == 0, x + y > 0), interact=True) # optional - qepcad
            sage: qe.go(); qe.go(); qe.go() # optional - qepcad
            QEPCAD object has moved to phase 'Before Projection (y)'
            QEPCAD object has moved to phase 'Before Choice'
            QEPCAD object has moved to phase 'Before Solution'
            sage: qe.solution_extension('E')                    # not tested (random order)
            x > _root_1 2 x^2 - 3 /\ y^2 + x^2 - 3 = 0 /\ [ 2 x^2 - 3 > 0 \/ y = _root_-1 y^2 + x^2 - 3 ]
            sage: qe.solution_extension('G')                    # not tested (random order)
            [
              [
                2 x^2 - 3 < 0
                \/
                x = _root_-1 2 x^2 - 3
              ]
              /\
              y = _root_-1 y^2 + x^2 - 3
            ]
            \/
            [
              x^2 - 3 <= 0
              /\
              x > _root_-1 2 x^2 - 3
              /\
              y^2 + x^2 - 3 = 0
            ]
            sage: qe.solution_extension('T')                    # not tested (random order)
            y + x > 0 /\ y^2 + x^2 - 3 = 0

        TESTS:

        Tests related to the not tested examples (nondeterministic order of atoms)::

            sage: from sage.interfaces.qepcad import _qepcad_atoms
            sage: var('x,y')
            (x, y)
            sage: qf = qepcad_formula
            sage: qe = qepcad(qf.and_(x^2 + y^2 - 3 == 0, x + y > 0), interact=True) # optional - qepcad
            sage: qe.go(); qe.go(); qe.go() # optional - qepcad
            QEPCAD object has moved to phase 'Before Projection (y)'
            QEPCAD object has moved to phase 'Before Choice'
            QEPCAD object has moved to phase 'Before Solution'
            sage: _qepcad_atoms(qe.solution_extension('E'))        # optional - qepcad
            {'2 x^2 - 3 > 0',
             'x > _root_1 2 x^2 - 3',
             'y = _root_-1 y^2 + x^2 - 3',
             'y^2 + x^2 - 3 = 0'}
            sage: _qepcad_atoms(qe.solution_extension('G'))        # optional - qepcad
            {'2 x^2 - 3 < 0',
             'x = _root_-1 2 x^2 - 3',
             'x > _root_-1 2 x^2 - 3',
             'x^2 - 3 <= 0',
             'y = _root_-1 y^2 + x^2 - 3',
             'y^2 + x^2 - 3 = 0'}
            sage: _qepcad_atoms(qe.solution_extension('T'))        # optional - qepcad
            {'y + x > 0', 'y^2 + x^2 - 3 = 0'}
        """
        if kind == 'I':
            raise ValueError("Interactive solution construction not "
                             "handled by Sage interface")
        result = self._eval_line('solution-extension %s' % kind)
        tagline = 'An equivalent quantifier-free formula:'
        loc = result.find(tagline)
        if loc >= 0:
            result = result[loc + len(tagline):]
        result = result.strip()
        if len(result):
            return AsciiArtString(result)

    def set_truth_value(self, index, nv):
        r"""
        Given a cell index (or a cell) and an integer, set the truth value
        of the cell to that integer.

        Valid integers are 0 (false), 1 (true), and 2 (undetermined).

        EXAMPLES::

            sage: qe = qepcad(x == 1, interact=True) # optional - qepcad
            sage: qe.go(); qe.go(); qe.go() # optional - qepcad
            QEPCAD object has moved to phase 'At the end of projection phase'
            QEPCAD object has moved to phase 'Before Choice'
            QEPCAD object has moved to phase 'Before Solution'
            sage: qe.set_truth_value(1, 1) # optional - qepcad
        """
        index_str = _format_cell_index([index])
        self._eval_line('set-truth-value\n%s\n%s' % (index_str, nv))

    def phase(self):
        r"""
        Return the current phase of the QEPCAD program.

        EXAMPLES::

            sage: qe = qepcad(x > 2/3, interact=True) # optional - qepcad
            sage: qe.phase() # optional - qepcad
            'Before Normalization'
            sage: qe.go() # optional - qepcad
            QEPCAD object has moved to phase 'At the end of projection phase'
            sage: qe.phase() # optional  - qepcad
            'At the end of projection phase'
            sage: qe.go() # optional - qepcad
            QEPCAD object has moved to phase 'Before Choice'
            sage: qe.phase() # optional - qepcad
            'Before Choice'
            sage: qe.go() # optional - qepcad
            QEPCAD object has moved to phase 'Before Solution'
            sage: qe.phase() # optional - qepcad
            'Before Solution'
            sage: qe.go() # optional - qepcad
            3 x - 2 > 0
            sage: qe.phase() # optional - qepcad
            'EXITED'
        """
        match = self._qex.expect().match
        if match == pexpect.EOF:
            return 'EXITED'
        else:
            return bytes_to_str(match.group(1))

    def _parse_answer_stats(self):
        r"""
        Parse the final string printed by QEPCAD, which should be
        the simplified quantifier-free formula followed by some
        statistics.  Return a pair of the formula and the statistics
        (both as strings).

        EXAMPLES::

            sage: qe = qepcad(x^2 > 2, interact=True) # optional - qepcad
            sage: qe.finish() # optional - qepcad
            x^2 - 2 > 0
            sage: (ans, stats) = qe._parse_answer_stats() # optional - qepcad
            sage: ans # optional - qepcad
            'x^2 - 2 > 0'
            sage: stats # random, optional - qepcad
            '-----------------------------------------------------------------------------\r\n0 Garbage collections, 0 Cells and 0 Arrays reclaimed, in 0 milliseconds.\r\n492514 Cells in AVAIL, 500000 Cells in SPACE.\r\n\r\nSystem time: 16 milliseconds.\r\nSystem time after the initialization: 4 milliseconds.\r\n-----------------------------------------------------------------------------\r\n'
        """
        if self.phase() != 'EXITED':
            raise ValueError("QEPCAD is not finished yet")
        final = bytes_to_str(self._qex.expect().before)
        match = re.search('\nAn equivalent quantifier-free formula:(.*)\n=+  The End  =+\r\n\r\n(.*)$', final, re.DOTALL)

        if match:
            return (match.group(1).strip(), match.group(2))
        else:
            return (final, '')

    def answer(self):
        r"""
        For a QEPCAD instance which is finished, return the
        simplified quantifier-free formula that it printed just before
        exiting.

        EXAMPLES::

            sage: qe = qepcad(x^3 - x == 0, interact=True)  # optional - qepcad
            sage: qe.finish()                               # not tested (random order)
            x - 1 <= 0 /\ x + 1 >= 0 /\ [ x = 0 \/ x - 1 = 0 \/ x + 1 = 0 ]
            sage: qe.answer()                               # not tested (random order)
            x - 1 <= 0 /\ x + 1 >= 0 /\ [ x = 0 \/ x - 1 = 0 \/ x + 1 = 0 ]

        TESTS:

        Tests related to the not tested examples (nondeterministic order of atoms)::

            sage: from sage.interfaces.qepcad import _qepcad_atoms
            sage: qe = qepcad(x^3 - x == 0, interact=True)  # optional - qepcad
            sage: qe.finish()                               # random, optional - qepcad
            x - 1 <= 0 /\ x + 1 >= 0 /\ [ x = 0 \/ x - 1 = 0 \/ x + 1 = 0 ]
            sage: _qepcad_atoms(qe.answer())                                       # optional - qepcad
            {'x + 1 = 0', 'x + 1 >= 0', 'x - 1 <= 0', 'x - 1 = 0', 'x = 0'}
        """
        return AsciiArtString(self._parse_answer_stats()[0])

    def final_stats(self):
        r"""
        For a QEPCAD instance which is finished, return the
        statistics that it printed just before exiting.

        EXAMPLES::

            sage: qe = qepcad(x == 0, interact=True)  # optional - qepcad
            sage: qe.finish()  # optional - qepcad
            x = 0
            sage: qe.final_stats()  # random, optional - qepcad
            -----------------------------------------------------------------------------
            0 Garbage collections, 0 Cells and 0 Arrays reclaimed, in 0 milliseconds.
            492840 Cells in AVAIL, 500000 Cells in SPACE.
            System time: 8 milliseconds.
            System time after the initialization: 4 milliseconds.
            -----------------------------------------------------------------------------
        """
        return AsciiArtString(self._parse_answer_stats()[1])

    def _tab_completion(self):
        r"""
        Return a list of the QEPCAD commands which are available as
        extra methods on a :class:`Qepcad` object.

        EXAMPLES::

            sage: qe = qepcad(x^2 < 0, interact=True) # optional - qepcad
            sage: len(qe._tab_completion()) # random, optional - qepcad
            97
            sage: 'd_cell' in qe._tab_completion() # optional - qepcad
            True
        """
        _update_command_info()
        return _command_info_cache.keys()

    def cell(self, *index):
        r"""
        Given a cell index, returns a :class:`QepcadCell` wrapper for that
        cell.  Uses a cache for efficiency.

        EXAMPLES::

            sage: qe = qepcad(x + 3 == 42, interact=True) # optional - qepcad
            sage: qe.go(); qe.go(); qe.go() # optional - qepcad
            QEPCAD object has moved to phase 'At the end of projection phase'
            QEPCAD object has moved to phase 'Before Choice'
            QEPCAD object has moved to phase 'Before Solution'
            sage: qe.cell(2) # optional - qepcad
            QEPCAD cell (2)
            sage: qe.cell(2) is qe.cell(2) # optional - qepcad
            True
        """
        index_str = _format_cell_index(index)
        if index_str in self._cell_cache:
            return self._cell_cache[index_str]
        else:
            c = self.make_cells(self.d_cell(index))[0]
            self._cell_cache[index_str] = c
            return c

    def make_cells(self, text):
        r"""
        Given the result of some QEPCAD command that returns cells
        (such as :meth:`d_cell`, :meth:`d_witness_list`, etc.),
        return a list of cell objects.

        EXAMPLES::

            sage: var('x,y')
            (x, y)
            sage: qe = qepcad(x^2 + y^2 == 1, interact=True) # optional - qepcad
            sage: qe.go(); qe.go(); qe.go() # optional - qepcad
            QEPCAD object has moved to phase 'Before Projection (y)'
            QEPCAD object has moved to phase 'Before Choice'
            QEPCAD object has moved to phase 'Before Solution'
            sage: qe.make_cells(qe.d_false_cells()) # optional - qepcad
            [QEPCAD cell (5, 1), QEPCAD cell (4, 3), QEPCAD cell (4, 1), QEPCAD cell (3, 5), QEPCAD cell (3, 3), QEPCAD cell (3, 1), QEPCAD cell (2, 3), QEPCAD cell (2, 1), QEPCAD cell (1, 1)]
        """
        # get rid of AsciiArtString
        text = str(text)

        lines = text.strip().splitlines()

        cells = []
        cell_lines = []
        in_cell = False

        for line in lines:
            if 'Information about the cell' in line:
                in_cell = True
            if in_cell:
                cell_lines.append(line)
            if line == '----------------------------------------------------':
                cells.append(QepcadCell(self, cell_lines))
                cell_lines = []
                in_cell = False

        return cells

    def __getattr__(self, attrname):
        r"""
        Return a :class:`QepcadFunction` object for any QEPCAD command.

        EXAMPLES::

            sage: qe = qepcad(x^3 == 8, interact=True) # optional - qepcad
            sage: qe.d_cell # optional - qepcad
            d_cell
        """
        if attrname[:1] == "_":
            raise AttributeError
        if attrname not in self._tab_completion():
            raise AttributeError
        return QepcadFunction(self, attrname)

    def _eval_line(self, cmd, restart_if_needed=False):
        r"""
        Send a command to QEPCAD, wait for a prompt, and return the
        text printed by QEPCAD before the prompt.

        Not intended for direct use.

        EXAMPLES::

            sage: qe = qepcad(x^2 == x^3, interact=True) # optional - qepcad
            sage: qe._eval_line('d-formula') # optional - qepcad
            '-x^3 + x^2 = 0'
        """
        # Evidently when the QEPCAD command parser reads a command, it
        # parses from the beginning of a line until it has found a command,
        # ignoring newlines.  It then throws away the rest of the line
        # after the end of the command.  There is no command terminator;
        # each command does its own argument parsing and executes
        # as soon as it sees all the arguments.

        # Thus, if the user gives a function call with a missing argument,
        # QEPCAD will wait forever for the missing argument to show up
        # (with no further prompting), and Sage will wait forever for
        # a prompt.  Not good.

        # To avoid this, we add a trailing "&" to each command.
        # If there is a missing argument, then the QEPCAD argument parsing
        # will give a syntax error on the "&", instead of waiting forever;
        # if there is no missing argument then QEPCAD will silently
        # ignore the "&".  (Hopefully all QEPCAD commands will treat "&"
        # as a syntax error; as far as I know that's the case.)

        # This trailing "&" also helps with another problem; sometimes
        # we get a copy of the command we send to QEPCAD echoed back as
        # the first line of the result.  I'm guessing this is because
        # readline turns raw mode on when it's reading and off when it isn't;
        # if we manage to send a command with raw mode off, then it's
        # echoed twice (once by the normal UNIX tty handling and once by
        # readline).  To work around this problem, we check for an "&"
        # in the first line of the input, and remove the text up to the "&"
        # if it is present.  (Hopefully QEPCAD doesn't print "&" as part
        # of its normal operation.)

        result = self._qex._eval_line(cmd + ' &')

        nl = result.find('\n')
        if nl < 0:
            nl = len(result)

        amp = result.find('&', 0, nl)
        if amp > 0:
            result = result[amp+1:]

        result = result.strip()

        return result

    def _function_call(self, name, args):
        r"""
        Given a command name and a list of arguments, send the command
        to QEPCAD and return the result.

        Not intended for direct use.

        EXAMPLES::

            sage: qe = qepcad(x^2 - 1 == 0, interact=True) # optional - qepcad
            sage: qe.go(); qe.go() # optional - qepcad
            QEPCAD object has moved to phase 'At the end of projection phase'
            QEPCAD object has moved to phase 'Before Choice'
            sage: qe._function_call('d_level_factors', (1,)) # optional - qepcad
            A_1,1  = input
                   = x + 1
            A_1,2  = input
                   = x - 1
        """
        name = name.replace('_', '-')
        args = [str(_) for _ in args]
        pre_phase = self.phase()
        result = self._eval_line('%s %s'%(name, ' '.join(args)))
        post_phase = self.phase()
        if len(result) and post_phase != 'EXITED':
            return AsciiArtString(result)
        if pre_phase != post_phase:
            if post_phase == 'EXITED' and name != 'quit':
                return self.answer()
            return AsciiArtString("QEPCAD object has moved to phase '%s'"%post_phase)

def _format_cell_index(a):
    """
    Given a tuple (or list, etc.) containing a QEPCAD cell index, return a
    string with a properly-formatted index.

    The input is flattened, so
    extra levels of brackets are ignored.  Also, if the first item
    in the flattened list is a :class:`QepcadCell` object, then its index
    is used in place of the object.

    EXAMPLES::

        sage: from sage.interfaces.qepcad import _format_cell_index

        sage: qe = qepcad(x == 0, interact=True); qe.go(); qe.go(); qe.go() # optional - qepcad
        QEPCAD object has moved to phase 'At the end of projection phase'
        QEPCAD object has moved to phase 'Before Choice'
        QEPCAD object has moved to phase 'Before Solution'
        sage: _format_cell_index(qe.cell(1)) # optional - qepcad
        '(1)'
        sage: _format_cell_index((qe.cell(1), 2, 3)) # optional - qepcad
        '(1, 2, 3)'
        sage: _format_cell_index(())
        '()'
        sage: _format_cell_index(5)
        '(5)'
    """
    a = flatten([a])
    if len(a) and isinstance(a[0], QepcadCell):
        a[0:1] = a[0].index()
    if len(a) == 1:
        return '(%s)' % a[0]
    else:
        return str(tuple(a))


@instancedoc
class QepcadFunction(ExpectFunction):
    r"""
    A wrapper for a QEPCAD command.
    """
    def _instancedoc_(self):
        r"""
        Return the documentation for a QEPCAD command, from
        ``qepcad.help``.

        EXAMPLES::

            sage: qe = qepcad(x == 17, interact=True) # optional - qepcad
            sage: cmd = qe.approx_precision # optional - qepcad
            sage: cmd.__doc__  # optional - qepcad
            'approx-precision N\n\nApproximate algeraic numbers to N decimal places.\n'
        """
        _update_command_info()
        return _command_info_cache[self._name][3]

    def __call__(self, *args):
        r"""
        Call QEPCAD with the command this is a wrapper for.

        For commands which take cell indexes, :func:`_format_cell_index`
        is automatically called.  For commands which take 'y' or 'n',
        booleans are also allowed.  (These special commands were hand-selected
        after reading ``qepcad.help``.)

        EXAMPLES::

            sage: qe = qepcad(x^2 < 1, interact=True) # optional - qepcad
            sage: cmd = qe.d_formula # optional - qepcad
            sage: cmd.__call__() # optional - qepcad
            x^2 - 1 < 0
        """
        _update_command_info()
        info = _command_info_cache[self._name]
        special = info[4]
        # Tweak the argument processing for a few commands.
        args = list(args)

        if special == 'cell':
            args = [_format_cell_index(args)]

        if special == 'yn':
            if isinstance(args[0], bool):
                args[0] = 'y' if args[0] else 'n'

        if special == 'interactive':
            raise ValueError("Cannot call %s through Sage interface... "
                             "interactive commands not handled")

        return self._parent._function_call(self._name, args)


def qepcad(formula, assume=None, interact=False, solution=None,
           vars=None, **kwargs):
    r"""
    Quantifier elimination and formula simplification using QEPCAD B.

    If ``assume`` is specified, then the given formula is ``'assumed'``,
    which is taken into account during final solution formula construction.

    If ``interact=True`` is given, then a :class:`Qepcad` object is
    returned which can be interacted with either at the command line
    or programmatically.

    The type of solution returned can be adjusted with ``solution``.
    The options are ``'geometric'``, which tries to construct a
    solution formula with geometric meaning; ``'extended'``, which
    gives a solution formula in an extended language that may be more
    efficient to construct; ``'any-point'``, which returns any
    point where the formula is true; ``'all-points'``, which
    returns a list of all points where the formula is true (or raises
    an exception if there are infinitely many); and ``'cell-points'``,
    which returns one point in each cell where the formula is true.

    All other keyword arguments are passed through to the :class:`Qepcad`
    constructor.

    For much more documentation and many more examples, see the module
    docstring for this module (type ``sage.interfaces.qepcad?`` to
    read this docstring from the \sage command line).

    The examples below require that the optional qepcad package is installed.

    EXAMPLES::

        sage: qf = qepcad_formula

        sage: var('a,b,c,d,x,y,z,long_with_underscore_314159')
        (a, b, c, d, x, y, z, long_with_underscore_314159)
        sage: K.<q,r> = QQ[]

        sage: qepcad('(E x)[a x + b > 0]', vars='(a,b,x)')      # not tested (random order)
        a /= 0 \/ b > 0

        sage: qepcad(a > b)                            # optional - qepcad
        b - a < 0

        sage: qepcad(qf.exists(x, a*x^2 + b*x + c == 0))        # not tested (random order)
        4 a c - b^2 <= 0 /\ [ c = 0 \/ a /= 0 \/ 4 a c - b^2 < 0 ]

        sage: qepcad(qf.exists(x, a*x^2 + b*x + c == 0), assume=(a != 0))    # optional - qepcad
        4 a c - b^2 <= 0

    For which values of $a$, $b$, $c$ does $a x^2 + b x + c$ have
    2 real zeroes? ::

        sage: exact2 = qepcad(qf.exactly_k(2, x, a*x^2 + b*x + c == 0)); exact2   # not tested (random order)
        a /= 0 /\ 4 a c - b^2 < 0

    one real zero? ::

        sage: exact1 = qepcad(qf.exactly_k(1, x, a*x^2 + b*x + c == 0)); exact1   # not tested (random order)
        [ a > 0 /\ 4 a c - b^2 = 0 ] \/ [ a < 0 /\ 4 a c - b^2 = 0 ] \/ [ a = 0 /\ 4 a c - b^2 < 0 ]

    No real zeroes? ::

        sage: exact0 = qepcad(qf.forall(x, a*x^2 + b*x + c != 0)); exact0     # not tested (random order)
        4 a c - b^2 >= 0 /\ c /= 0 /\ [ b = 0 \/ 4 a c - b^2 > 0 ]

    $3^{75}$ real zeroes? ::

        sage: qepcad(qf.exactly_k(3^75, x, a*x^2 + b*x + c == 0))    # optional - qepcad
        FALSE

    We can check that the results don't overlap::

        sage: qepcad(r'[[%s] /\ [%s]]' % (exact0, exact1), vars='a,b,c')      # not tested (random order)
        FALSE
        sage: qepcad(r'[[%s] /\ [%s]]' % (exact0, exact2), vars='a,b,c')      # not tested (random order)
        FALSE
        sage: qepcad(r'[[%s] /\ [%s]]' % (exact1, exact2), vars='a,b,c')      # not tested (random order)
        FALSE

    and that the union of the results is as expected::

        sage: qepcad(r'[[%s] \/ [%s] \/ [%s]]' % (exact0, exact1, exact2), vars=(a,b,c))  # not tested (random order)
        b /= 0 \/ a /= 0 \/ c /= 0

    So we have finitely many zeroes if $a$, $b$, or $c$ is nonzero;
    which means we should have infinitely many zeroes if they are all
    zero. ::

        sage: qepcad(qf.infinitely_many(x, a*x^2 + b*x + c == 0))          # not tested (random order)
        a = 0 /\ b = 0 /\ c = 0

    The polynomial is nonzero almost everywhere iff it is not
    identically zero. ::

        sage: qepcad(qf.all_but_finitely_many(x, a*x^2 + b*x + c != 0))    # not tested (random order)
        b /= 0 \/ a /= 0 \/ c /= 0

    The non-zeroes are continuous iff there are no zeroes or if
    the polynomial is zero. ::

        sage: qepcad(qf.connected_subset(x, a*x^2 + b*x + c != 0))         # not tested (random order)
        4 a c - b^2 >= 0 /\ [ a = 0 \/ 4 a c - b^2 > 0 ]

    The zeroes are continuous iff there are no or one zeroes, or if the
    polynomial is zero::

        sage: qepcad(qf.connected_subset(x, a*x^2 + b*x + c == 0))         # not tested (random order)
        a = 0 \/ 4 a c - b^2 >= 0
        sage: qepcad(r'[[%s] \/ [%s] \/ [a = 0 /\ b = 0 /\ c = 0]]' % (exact0, exact1), vars='a,b,c')   # not tested (random order)
        a = 0 \/ 4 a c - b^2 >= 0

    Since polynomials are continuous and $y > 0$ is an open set,
    they are positive infinitely often iff they are positive at
    least once. ::

        sage: qepcad(qf.infinitely_many(x, a*x^2 + b*x + c > 0))        # not tested (random order)
        c > 0 \/ a > 0 \/ 4 a c - b^2 < 0
        sage: qepcad(qf.exists(x, a*x^2 + b*x + c > 0))                 # not tested (random order)
        c > 0 \/ a > 0 \/ 4 a c - b^2 < 0

    However, since $y >= 0$ is not open, the equivalence does not
    hold if you replace 'positive' with 'nonnegative'.
    (We assume $a \neq 0$ to get simpler formulas.) ::

        sage: qepcad(qf.infinitely_many(x, a*x^2 + b*x + c >= 0), assume=(a != 0))    # not tested (random order)
        a > 0 \/ 4 a c - b^2 < 0
        sage: qepcad(qf.exists(x, a*x^2 + b*x + c >= 0), assume=(a != 0))             # not tested (random order)
        a > 0 \/ 4 a c - b^2 <= 0

    TESTS:

    We verify that long variable names work.  (Note that QEPCAD
    does not support underscores, so they are stripped from the formula.) ::

        sage: qepcad(qf.exists(a, a*long_with_underscore_314159 == 1))                # optional - qepcad
        longwithunderscore314159 /= 0

    Tests related to the not tested examples (nondeterministic order of atoms)::

        sage: from sage.interfaces.qepcad import _qepcad_atoms
        sage: var('a,b,c,d,x,y,z,long_with_underscore_314159')
        (a, b, c, d, x, y, z, long_with_underscore_314159)
        sage: K.<q,r> = QQ[]

        sage: _qepcad_atoms(qepcad('(E x)[a x + b > 0]', vars='(a,b,x)'))                       # optional - qepcad
        {'a /= 0', 'b > 0'}

        sage: _qepcad_atoms(qepcad(qf.exists(x, a*x^2 + b*x + c == 0)))                         # optional - qepcad
        {'4 a c - b^2 < 0', '4 a c - b^2 <= 0', 'a /= 0', 'c = 0'}

        sage: exact2 = qepcad(qf.exactly_k(2, x, a*x^2 + b*x + c == 0)); _qepcad_atoms(exact2)  # optional - qepcad
        {'4 a c - b^2 < 0', 'a /= 0'}

        sage: exact1 = qepcad(qf.exactly_k(1, x, a*x^2 + b*x + c == 0)); _qepcad_atoms(exact1)  # optional - qepcad
        {'4 a c - b^2 < 0', '4 a c - b^2 = 0', 'a < 0', 'a = 0', 'a > 0'}

        sage: exact0 = qepcad(qf.forall(x, a*x^2 + b*x + c != 0)); _qepcad_atoms(exact0)        # optional - qepcad
        {'4 a c - b^2 > 0', '4 a c - b^2 >= 0', 'b = 0', 'c /= 0'}

        sage: qepcad(r'[[%s] /\ [%s]]' % (exact0, exact1), vars='a,b,c')                        # optional - qepcad
        FALSE
        sage: qepcad(r'[[%s] /\ [%s]]' % (exact0, exact2), vars='a,b,c')                        # optional - qepcad
        FALSE
        sage: qepcad(r'[[%s] /\ [%s]]' % (exact1, exact2), vars='a,b,c')                        # optional - qepcad
        FALSE

        sage: _qepcad_atoms(qepcad(r'[[%s] \/ [%s] \/ [%s]]' % (exact0, exact1, exact2), vars=(a,b,c)))    # optional - qepcad
        {'a /= 0', 'b /= 0', 'c /= 0'}

        sage: _qepcad_atoms(qepcad(qf.infinitely_many(x, a*x^2 + b*x + c == 0)))               # optional - qepcad
        {'a = 0', 'b = 0', 'c = 0'}

        sage: _qepcad_atoms(qepcad(qf.all_but_finitely_many(x, a*x^2 + b*x + c != 0)))         # optional - qepcad
        {'a /= 0', 'b /= 0', 'c /= 0'}

        sage: _qepcad_atoms(qepcad(qf.connected_subset(x, a*x^2 + b*x + c != 0)))              # optional - qepcad
        {'4 a c - b^2 > 0', '4 a c - b^2 >= 0', 'a = 0'}

        sage: _qepcad_atoms(qepcad(qf.connected_subset(x, a*x^2 + b*x + c == 0)))              # optional - qepcad
        {'4 a c - b^2 >= 0', 'a = 0'}

        sage: _qepcad_atoms(qepcad(r'[[%s] \/ [%s] \/ [a = 0 /\ b = 0 /\ c = 0]]' % (exact0, exact1), vars='a,b,c'))   # optional - qepcad
        {'4 a c - b^2 >= 0', 'a = 0'}

        sage: _qepcad_atoms(qepcad(qf.infinitely_many(x, a*x^2 + b*x + c > 0)))                # optional - qepcad
        {'4 a c - b^2 < 0', 'a > 0', 'c > 0'}

        sage: _qepcad_atoms(qepcad(qf.exists(x, a*x^2 + b*x + c > 0)))                         # optional - qepcad
        {'4 a c - b^2 < 0', 'a > 0', 'c > 0'}

        sage: _qepcad_atoms(qepcad(qf.infinitely_many(x, a*x^2 + b*x + c >= 0), assume=(a != 0)))  # optional - qepcad
        {'4 a c - b^2 < 0', 'a > 0'}

        sage: _qepcad_atoms(qepcad(qf.exists(x, a*x^2 + b*x + c >= 0), assume=(a != 0)))       # optional - qepcad
        {'4 a c - b^2 <= 0', 'a > 0'}
    """
    use_witness = False
    if solution == 'any-point':
        formula = qepcad_formula.formula(formula)
        if len(formula.qvars) == 0:
            if vars is None:
                vars = sorted(formula.vars)
            formula = qepcad_formula.exists(vars, formula)
            vars = None
            use_witness = True

    qe = Qepcad(formula, vars=vars, **kwargs)
    if assume is not None:
        qe.assume(assume)
    if interact:
        if solution is not None:
            print("WARNING: 'solution=' is ignored for interactive use")
        return qe
    else:
        qe.go()
        qe.go()
        qe.go()
        if solution is None:
            qe.finish()
            return qe.answer()
        elif solution == 'geometric':
            s = qe.solution_extension('G')
            qe.quit()
            return s
        elif solution == 'extended':
            s = qe.solution_extension('E')
            qe.quit()
            return s
        elif solution == 'any-point':
            if use_witness:
                cells = qe.make_cells(qe.d_witness_list())
            else:
                cells = qe.make_cells(qe.d_true_cells())
            qe.quit()
            if len(cells) == 0:
                raise ValueError("input formula is false everywhere")
            return cells[0].sample_point_dict()
        elif solution == 'cell-points':
            cells = qe.make_cells(qe.d_true_cells())
            qe.quit()
            return [c.sample_point_dict() for c in cells]
        elif solution == 'all-points':
            cells = qe.make_cells(qe.d_true_cells())
            qe.quit()
            for c in cells:
                if c._dimension > 0:
                    raise ValueError("input formula is true for "
                                     "infinitely many points")
            return [c.sample_point_dict() for c in cells]
        else:
            raise ValueError("Unknown solution type ({})".format(solution))


import os
def qepcad_console(memcells=None):
    r"""
    Run QEPCAD directly.  To exit early, press Control-C.

    EXAMPLES::

        sage: qepcad_console() # not tested
        ...
        Enter an informal description  between '[' and ']':
    """
    from sage.repl.rich_output.display_manager import get_display_manager
    if not get_display_manager().is_in_terminal():
        raise RuntimeError('Can use the console only in the terminal. Try %%qepcad magics instead.')
    # This will only spawn local processes
    os.system(_qepcad_cmd(memcells))


def qepcad_banner():
    """
    Return the QEPCAD startup banner.

    EXAMPLES::

        sage: from sage.interfaces.qepcad import qepcad_banner
        sage: qepcad_banner() # optional - qepcad
        =======================================================
                        Quantifier Elimination
                                  in
                    Elementary Algebra and Geometry
                                  by
              Partial Cylindrical Algebraic Decomposition
        ...
                                  by
                               Hoon Hong
                         (hhong@math.ncsu.edu)
        With contributions by: Christopher W. Brown, George E.
        Collins, Mark J. Encarnacion, Jeremy R. Johnson
        Werner Krandick, Richard Liska, Scott McCallum,
        Nicolas Robidoux, and Stanly Steinberg
        =======================================================
    """
    qex = Qepcad_expect()
    qex._start()
    banner = bytes_to_str(qex.expect().before)
    return AsciiArtString(banner)

def qepcad_version():
    """
    Return a string containing the current QEPCAD version number.

    EXAMPLES::

        sage: qepcad_version() # random, optional - qepcad
        'Version B 1.69, 16 Mar 2012'

    TESTS::

        sage: qepcad_version() # optional - qepcad
        'Version B ..., ...'
    """
    banner = str(qepcad_banner())
    lines = banner.split('\n')
    for line in lines:
        if 'Version' in line:
            return line.strip()


class qformula:
    """
    A qformula holds a string describing a formula in QEPCAD's syntax,
    and a set of variables used.
    """
    def __init__(self, formula, vars, qvars=[]):
        r"""
        Construct a qformula from a string, a frozenset of variable names,
        and (optionally) a list of ordered quantified names.

        EXAMPLES::

            sage: from sage.interfaces.qepcad import qformula
            sage: f = qformula('x + y = 0', frozenset(['x','y']))
            sage: f
            x + y = 0
            sage: f.formula
            'x + y = 0'
            sage: f.vars
            frozenset({'x', 'y'})
            sage: f.qvars
            []
        """
        self.formula = formula
        self.vars = vars
        self.qvars = qvars

    def __repr__(self):
        r"""
        Return a string representation of a qformula (this is just
        the formula it holds).

        EXAMPLES::

            sage: from sage.interfaces.qepcad import qformula
            sage: f = qformula('x + y = 0', frozenset(['x','y']))
            sage: f.__repr__()
            'x + y = 0'
        """
        return self.formula


class qepcad_formula_factory:
    r"""
    Contains routines to help construct formulas in QEPCAD syntax.
    """

    def _normalize_op(self, op):
        r"""
        Given a relational operator (either a string, or a function from
        the \module{operator} module) return the corresponding QEPCAD
        operator.

        EXAMPLES::

            sage: qf = qepcad_formula
            sage: qf._normalize_op(operator.ne)
            '/='
            sage: qf._normalize_op('==')
            '='
        """
        import operator
        if op == operator.eq:
            return '='
        if op == operator.ne:
            return '/='
        if op == operator.lt:
            return '<'
        if op == operator.gt:
            return '>'
        if op == operator.le:
            return '<='
        if op == operator.ge:
            return '>='

        if op == '==':
            return '='
        if op == '!=':
            return '/='

        return op

    def _varset(self, p):
        r"""
        Given a polynomial or a symbolic expression, return a frozenset
        of the variables involved.

        EXAMPLES::

            sage: qf = qepcad_formula
            sage: var('x,y')
            (x, y)
            sage: K.<p,q> = QQ[]
            sage: qf._varset(x)
            frozenset({'x'})
            sage: qf._varset(x*y)
            frozenset({'x', 'y'})
            sage: qf._varset(q)
            frozenset({'q'})
            sage: qf._varset(p*q)
            frozenset({'p', 'q'})
        """
        try:
            vars = p.variables()
        except AttributeError:
            vars = []
        return frozenset([str(v) for v in vars])

    def _combine_formulas(self, formulas):
        r"""
        Given a list of formulas, convert them to strings and return
        the free variables found in the formulas.

        Since this is used to build boolean expressions for QEPCAD,
        and QEPCAD boolean expressions cannot be built from quantified
        formulas, an exception is raised if an input formula is quantified.

        INPUT:

        - ``formulas`` -- a list of unquantified formulas

        OUTPUT:

        - form_strs -- a list of formulas as strings
        - vars -- a frozenset of all variables in the formulas

        EXAMPLES::

            sage: var('x,y')
            (x, y)
            sage: qf = qepcad_formula
            sage: qf._combine_formulas([x^2 == 0, y < 17])
            (['x^2 = 0', 'y < 17'], frozenset({'x', 'y'}))
        """
        formulas = [self.atomic(f) for f in formulas]
        formula_strs = [repr(f) for f in formulas]
        vars = frozenset()
        for f in formulas:
            vars = vars | f.vars
            if len(f.qvars):
                raise ValueError("QEPCAD formulas must be in prenex"
                                 " (quantifiers outermost) form")
        return formula_strs, vars

    def atomic(self, lhs, op='=', rhs=0):
        r"""
        Construct a QEPCAD formula from the given inputs.

        INPUT:

        - ``lhs`` -- a polynomial, or a symbolic equality or inequality
        - ``op`` -- a relational operator, default '='
        - ``rhs`` -- a polynomial, default 0

        If ``lhs`` is a symbolic equality or inequality, then ``op``
        and ``rhs`` are ignored.

        This method works by printing the given polynomials, so we do not
        care what ring they are in (as long as they print with integral
        or rational coefficients).

        EXAMPLES::

            sage: qf = qepcad_formula
            sage: var('a,b,c')
            (a, b, c)
            sage: K.<x,y> = QQ[]
            sage: def test_qf(qf):
            ....:     return qf, qf.vars
            sage: test_qf(qf.atomic(a^2 + 17))
            (a^2 + 17 = 0, frozenset({'a'}))
            sage: test_qf(qf.atomic(a*b*c <= c^3))
            (a b c <= c^3, frozenset({'a', 'b', 'c'}))
            sage: test_qf(qf.atomic(x+y^2, '!=', a+b))
            (y^2 + x /= a + b, frozenset({'a', 'b', 'x', 'y'}))
            sage: test_qf(qf.atomic(x, operator.lt))
            (x < 0, frozenset({'x'}))
        """
        if isinstance(lhs, qformula):
            return lhs

        from sage.structure.element import Expression
        if isinstance(lhs, Expression) and lhs.is_relational():
            lhs, op, rhs = lhs.lhs(), lhs.operator(), lhs.rhs()

        op = self._normalize_op(op)

        formula = ('%r %s %r' % (lhs, op, rhs))
        formula = formula.replace('*', ' ')
        vars = self._varset(lhs) | self._varset(rhs)

        return qformula(formula, vars)

    def formula(self, formula):
        r"""
        Constructs a QEPCAD formula from the given input.

        INPUT:

        - ``formula`` -- a polynomial, a symbolic equality or inequality,
          or a list of polynomials, equalities, or inequalities

        A polynomial $p$ is interpreted as the equation $p = 0$.
        A list is interpreted as the conjunction ('and') of the elements.

        EXAMPLES::

            sage: var('a,b,c,x')
            (a, b, c, x)
            sage: qf = qepcad_formula
            sage: qf.formula(a*x + b)
            a x + b = 0
            sage: qf.formula((a*x^2 + b*x + c, a != 0))
            [a x^2 + b x + c = 0 /\ a /= 0]
        """
        if isinstance(formula, (list, tuple)):
            return self.and_(formula)
        else:
            return self.atomic(formula)

    def and_(self, *formulas):
        r"""
        Return the conjunction of its input formulas.

        (This method would be named 'and' if that were not a Python
        keyword.)

        Each input formula may be a :class:`qformula` as returned by the
        methods of ``qepcad_formula``, a symbolic equality or
        inequality, or a polynomial $p$ (meaning $p = 0$).

        EXAMPLES::

            sage: var('a,b,c,x')
            (a, b, c, x)
            sage: qf = qepcad_formula
            sage: qf.and_(a*b, a*c, b*c != 0)
            [a b = 0 /\ a c = 0 /\ b c /= 0]
            sage: qf.and_(a*x^2 == 3, qf.or_(a > b, b > c))
            [a x^2 = 3 /\ [a > b \/ b > c]]
        """
        formulas = flatten(formulas)
        formula_strs, vars = self._combine_formulas(formulas)
        formula = '[' + r' /\ '.join(formula_strs) + ']'
        return qformula(formula, vars)

    def or_(self, *formulas):
        r"""
        Return the disjunction of its input formulas.

        (This method would be named 'or' if that were not a Python
        keyword.)

        Each input formula may be a :class:`qformula` as returned by the
        methods of ``qepcad_formula``, a symbolic equality or
        inequality, or a polynomial $p$ (meaning $p = 0$).

        EXAMPLES::

            sage: var('a,b,c,x')
            (a, b, c, x)
            sage: qf = qepcad_formula
            sage: qf.or_(a*b, a*c, b*c != 0)
            [a b = 0 \/ a c = 0 \/ b c /= 0]
            sage: qf.or_(a*x^2 == 3, qf.and_(a > b, b > c))
            [a x^2 = 3 \/ [a > b /\ b > c]]
        """
        formulas = flatten(formulas)
        formula_strs, vars = self._combine_formulas(formulas)
        formula = '[' + r' \/ '.join(formula_strs) + ']'
        return qformula(formula, vars)

    def not_(self, formula):
        r"""
        Return the negation of its input formula.

        (This method would be named 'not' if that were not a Python
        keyword.)

        The input formula may be a :class:`qformula` as returned by the
        methods of ``qepcad_formula``, a symbolic equality or
        inequality, or a polynomial $p$ (meaning $p = 0$).

        EXAMPLES::

            sage: var('a,b')
            (a, b)
            sage: qf = qepcad_formula
            sage: qf.not_(a > b)
            [~a > b]
            sage: qf.not_(a^2 + b^2)
            [~a^2 + b^2 = 0]
            sage: qf.not_(qf.and_(a > 0, b < 0))
            [~[a > 0 /\ b < 0]]
        """
        (formula_str,), vars = self._combine_formulas([formula])
        f = '[~' + formula_str + ']'
        return qformula(f, vars)

    def implies(self, f1, f2):
        r"""
        Return the implication of its input formulas (that is, given
        formulas $P$ and $Q$, returns '$P$ implies $Q$').

        The input formulas may be a :class:`qformula` as returned by the
        methods of ``qepcad_formula``, a symbolic equality or
        inequality, or a polynomial $p$ (meaning $p = 0$).

        EXAMPLES::

            sage: var('a,b')
            (a, b)
            sage: qf = qepcad_formula
            sage: qf.implies(a, b)
            [a = 0 ==> b = 0]
            sage: qf.implies(a^2 < b, b^2 < a)
            [a^2 < b ==> b^2 < a]
        """
        (f1s, f2s), vars = self._combine_formulas([f1, f2])
        f = '[' + f1s + ' ==> ' + f2s + ']'
        return qformula(f, vars)

    # QEPCAD also has a "reverse implication" symbol "<==", which
    # I'm not bothering to support.

    def iff(self, f1, f2):
        r"""
        Return the equivalence of its input formulas (that is, given
        formulas $P$ and $Q$, returns '$P$ iff $Q$').

        The input formulas may be a :class:`qformula` as returned by the
        methods of ``qepcad_formula``, a symbolic equality or
        inequality, or a polynomial $p$ (meaning $p = 0$).

        EXAMPLES::

            sage: var('a,b')
            (a, b)
            sage: qf = qepcad_formula
            sage: qf.iff(a, b)
            [a = 0 <==> b = 0]
            sage: qf.iff(a^2 < b, b^2 < a)
            [a^2 < b <==> b^2 < a]
        """
        (f1s, f2s), vars = self._combine_formulas([f1, f2])
        f = '[' + f1s + ' <==> ' + f2s + ']'
        return qformula(f, vars)

    def exists(self, v, formula):
        r"""
        Given a variable (or list of variables) and a formula, returns
        the existential quantification of the formula over the variables.

        This method is available both as :meth:`exists` and :meth:`E`
        (the QEPCAD name for this quantifier).

        The input formula may be a :class:`qformula` as returned by the
        methods of ``qepcad_formula``, a symbolic equality or
        inequality, or a polynomial $p$ (meaning $p = 0$).

        EXAMPLES::

            sage: var('a,b')
            (a, b)
            sage: qf = qepcad_formula
            sage: qf.exists(a, a^2 + b > b^2 + a)
            (E a)[a^2 + b > b^2 + a]
            sage: qf.exists((a, b), a^2 + b^2 < 0)
            (E a)(E b)[a^2 + b^2 < 0]
            sage: qf.E(b, b^2 == a)
            (E b)[b^2 = a]
        """
        return self.quantifier('E', v, formula)
    E = exists

    def forall(self, v, formula):
        r"""
        Given a variable (or list of variables) and a formula, returns
        the universal quantification of the formula over the variables.

        This method is available both as :meth:`forall` and :meth:`A`
        (the QEPCAD name for this quantifier).

        The input formula may be a :class:`qformula` as returned by the
        methods of ``qepcad_formula``, a symbolic equality or
        inequality, or a polynomial $p$ (meaning $p = 0$).

        EXAMPLES::

            sage: var('a,b')
            (a, b)
            sage: qf = qepcad_formula
            sage: qf.forall(a, a^2 + b > b^2 + a)
            (A a)[a^2 + b > b^2 + a]
            sage: qf.forall((a, b), a^2 + b^2 > 0)
            (A a)(A b)[a^2 + b^2 > 0]
            sage: qf.A(b, b^2 != a)
            (A b)[b^2 /= a]
        """
        return self.quantifier('A', v, formula)
    A = forall

    def infinitely_many(self, v, formula):
        r"""
        Given a variable and a formula, returns a new formula which is
        true iff the original formula was true for infinitely many
        values of the variable.

        This method is available both as :meth:`infinitely_many`
        and :meth:`F` (the QEPCAD name for this quantifier).

        The input formula may be a :class:`qformula` as returned by the
        methods of ``qepcad_formula``, a symbolic equality or
        inequality, or a polynomial $p$ (meaning $p = 0$).

        EXAMPLES::

            sage: var('a,b')
            (a, b)
            sage: qf = qepcad_formula
            sage: qf.infinitely_many(a, a^2 + b > b^2 + a)
            (F a)[a^2 + b > b^2 + a]
            sage: qf.F(b, b^2 != a)
            (F b)[b^2 /= a]
        """
        return self.quantifier('F', v, formula, allow_multi=False)
    F = infinitely_many

    def all_but_finitely_many(self, v, formula):
        r"""
        Given a variable and a formula, returns a new formula which is
        true iff the original formula was true for all but finitely many
        values of the variable.

        This method is available both as :meth:`all_but_finitely_many`
        and :meth:`G` (the QEPCAD name for this quantifier).

        The input formula may be a :class:`qformula` as returned by the
        methods of ``qepcad_formula``, a symbolic equality or
        inequality, or a polynomial $p$ (meaning $p = 0$).

        EXAMPLES::

            sage: var('a,b')
            (a, b)
            sage: qf = qepcad_formula
            sage: qf.all_but_finitely_many(a, a^2 + b > b^2 + a)
            (G a)[a^2 + b > b^2 + a]
            sage: qf.G(b, b^2 != a)
            (G b)[b^2 /= a]
        """
        return self.quantifier('G', v, formula, allow_multi=False)
    G = all_but_finitely_many

    def connected_subset(self, v, formula, allow_multi=False):
        r"""
        Given a variable and a formula, returns a new formula which is
        true iff the set of values for the variable at which the
        original formula was true is connected (including cases where
        this set is empty or is a single point).

        This method is available both as :meth:`connected_subset`
        and :meth:`C` (the QEPCAD name for this quantifier).

        The input formula may be a :class:`qformula` as returned by the
        methods of ``qepcad_formula``, a symbolic equality or
        inequality, or a polynomial $p$ (meaning $p = 0$).

        EXAMPLES::

            sage: var('a,b')
            (a, b)
            sage: qf = qepcad_formula
            sage: qf.connected_subset(a, a^2 + b > b^2 + a)
            (C a)[a^2 + b > b^2 + a]
            sage: qf.C(b, b^2 != a)
            (C b)[b^2 /= a]
        """
        return self.quantifier('C', v, formula)
    C = connected_subset

    def exactly_k(self, k, v, formula, allow_multi=False):
        r"""
        Given a nonnegative integer $k$, a variable, and a formula,
        returns a new formula which is true iff the original formula
        is true for exactly $k$ values of the variable.

        This method is available both as :meth:`exactly_k`
        and :meth:`X` (the QEPCAD name for this quantifier).

        (Note that QEPCAD does not support $k=0$ with this syntax, so if
        $k=0$ is requested we implement it with :meth:`forall` and
        :meth:`not_`.)

        The input formula may be a :class:`qformula` as returned by the
        methods of ``qepcad_formula``, a symbolic equality or
        inequality, or a polynomial $p$ (meaning $p = 0$).

        EXAMPLES::

            sage: var('a,b')
            (a, b)
            sage: qf = qepcad_formula
            sage: qf.exactly_k(1, x, x^2 + a*x + b == 0)
            (X1 x)[a x + x^2 + b = 0]
            sage: qf.exactly_k(0, b, a*b == 1)
            (A b)[~a b = 1]
        """
        from sage.rings.integer_ring import ZZ
        k = ZZ(k)
        if k < 0:
            raise ValueError("negative k in exactly_k quantifier")

        if k == 0:
            return self.forall(v, self.not_(formula))

        return self.quantifier('X%s' % k, v, formula)
    X = exactly_k

    def quantifier(self, kind, v, formula, allow_multi=True):
        r"""
        A helper method for building quantified QEPCAD formulas; not
        expected to be called directly.

        Takes the quantifier kind (the string label of this quantifier),
        a variable or list of variables, and a formula, and returns
        the quantified formula.

        EXAMPLES::

            sage: var('a,b')
            (a, b)
            sage: qf = qepcad_formula
            sage: qf.quantifier('NOT_A_REAL_QEPCAD_QUANTIFIER', a, a*b==0)
            (NOT_A_REAL_QEPCAD_QUANTIFIER a)[a b = 0]
            sage: qf.quantifier('FOO', (a, b), a*b)
            (FOO a)(FOO b)[a b = 0]
        """

        formula = self.formula(formula)

        if allow_multi and isinstance(v, (list, tuple)):
            if len(v) == 0:
                return formula
            else:
                return self.quantifier(kind, v[0],
                                       self.quantifier(kind, v[1:], formula))

        form_str = str(formula)
        if form_str[-1] != ']':
            form_str = '[' + form_str + ']'
        v = str(v)
        if not (v in formula.vars):
            raise ValueError("Attempting to quantify variable which "
                             "does not occur in formula")
        form_str = "(%s %s)%s" % (kind, v, form_str)
        return qformula(form_str, formula.vars - frozenset([v]),
                        [v] + formula.qvars)


qepcad_formula = qepcad_formula_factory()

_qepcad_algebraic_re = \
    re.compile(' ?the unique root of (.*) between (.*) and (.*)$')


def _eval_qepcad_algebraic(text):
    r"""
    Given a string of the form:
    'the unique root of 8 x^2 - 8 x - 29 between -47/32 and -1503/1024'
    (as produced by QEPCAD) this returns the corresponding \sage
    algebraic real number.

    Requires that the given rational bounds are exactly representable as
    arbitrary-precision floating-point numbers (that is, that the
    denominators are powers of two); this is true of the expressions
    given by QEPCAD.

    EXAMPLES::

        sage: from sage.interfaces.qepcad import _eval_qepcad_algebraic
        sage: x = _eval_qepcad_algebraic('the unique root of 8 x^2 - 8 x - 29 between -47/32 and -1503/1024'); x
        -1.468501968502953?
        sage: 8*x^2 - 8*x - 29 == 0
        True
    """
    from sage.rings.rational_field import QQ
    from sage.rings.polynomial.polynomial_ring import polygen
    from sage.rings.real_mpfi import RealIntervalField
    from sage.rings.qqbar import AA

    match = _qepcad_algebraic_re.match(text)

    p_text = match.group(1)
    lbound = QQ(match.group(2))
    ubound = QQ(match.group(3))

    p = sage_eval(implicit_mul(p_text), {'x': polygen(QQ)})

    for prec_scale in range(15):
        fld = RealIntervalField(53 << prec_scale)
        intv = fld(lbound, ubound)
        if intv.lower().exact_rational() == lbound and intv.upper().exact_rational() == ubound:
            return AA.polynomial_root(p, intv)

    raise ValueError("%s or %s not an exact floating-point number" % (lbound,
                                                                      ubound))


class QepcadCell:
    r"""
    A wrapper for a QEPCAD cell.
    """
    def __init__(self, parent, lines):
        r"""
        Construct a :class:`QepcadCell` wrapper for a QEPCAD cell, given
        a :class:`Qepcad` object and a list of lines holding QEPCAD's
        cell output.

        This is typically called by the :meth:`Qepcad.make_cells`
        method of :class:`Qepcad`.

        EXAMPLES::

            sage: from sage.interfaces.qepcad import QepcadCell
            sage: var('a,b')
            (a, b)
            sage: qe = qepcad(a^2 + b^2 == 5, interact=True) # optional - qepcad
            sage: qe.go(); qe.go(); qe.go() # optional - qepcad
            QEPCAD object has moved to phase 'Before Projection (b)'
            QEPCAD object has moved to phase 'Before Choice'
            QEPCAD object has moved to phase 'Before Solution'
            sage: qe.d_cell(4, 3) # optional - qepcad
            ---------- Information about the cell (4,3) ----------
            Level                       : 2
            Dimension                   : 1
            Number of children          : 0
            Truth value                 : F    by trial evaluation.
            Degrees after substitution  : Not known yet or No polynomial.
            Multiplicities              : ()
            Signs of Projection Factors
            Level 1  : (0)
            Level 2  : (+)
            ----------   Sample point  ----------
            The sample point is in a PRIMITIVE representation.
            <BLANKLINE>
            alpha = the unique root of x^2 - 5 between 2289/1024 and 1145/512
                  = 2.2360679775-
            <BLANKLINE>
            Coordinate 1 = alpha
                         = 2.2360679775-
            Coordinate 2 = 1
                         = 1.0000000000
            ----------------------------------------------------

            sage: QepcadCell(qe, str(qe.d_cell(4, 3)).splitlines()) # optional - qepcad
            QEPCAD cell (4, 3)
        """
        self._parent = parent
        self._lines = lines

        max_level = len(parent._varlist)

        all_signs = []
        saw_signs = False

        saw_primitive = False
        saw_extended = False

        grab_extended = False

        all_coordinates = []

        for line in lines:
            if 'Information about the cell' in line:
                tail = line.split('(')[1]
                index = tail.split(')')[0]
                if index == '':
                    index = ()
                else:
                    index = sage_eval(index)
                    if not isinstance(index, tuple):
                        index = (index,)
                self._index = index

                self._dimension = sum([r&1 for r in index])
            if 'Level        ' in line:
                self._level = int(line.split(':')[1].strip())
            if 'Number of children' in line:
                nkids = int(line.split(':')[1].strip())
                if nkids > 0 or self._level == max_level:
                    self._number_of_children = nkids
                else:
                    self._number_of_children = None
            if 'Truth value' in line:
                pass # might change
            if 'Degrees after substitution' in line:
                if self._level == max_level or self._level == 0:
                    self._degrees = None
                else:
                    self._degrees = sage_eval(line.split(':')[1].strip())
            if 'Multiplicities' in line:
                self._multiplicities = sage_eval(line.split(':')[1].strip())
            if 'Signs of Projection Factors' in line:
                saw_signs = True
            if saw_signs and 'Level' in line:
                (lev, n, colon, signs) = line.split()
                assert(lev == 'Level' and colon == ':')
                assert(int(n) == len(all_signs) + 1)
                signs = signs.replace('+','1').replace('-','-1').replace(')',',)')
                all_signs.append(sage_eval(signs))
            if 'PRIMITIVE' in line:
                saw_primitive = True
            if 'EXTENDED' in line:
                saw_extended = True

            if 'alpha = ' in line:
                alpha_text = line[8:]
            if 'Coordinate ' in line:
                (coord_n, val) = line.split('=')
                n = int(coord_n.split()[1])
                assert(n == len(all_coordinates) + 1)
                if n == self._level and saw_extended:
                    grab_extended = True
                else:
                    all_coordinates.append(val)
            elif grab_extended:
                assert('=' in line)
                grab_extended = False
                all_coordinates.append(line.split('=')[1])

        if saw_signs:
            self._signs = all_signs

        if saw_primitive or saw_extended:
            self._raw_sample_point = (saw_extended, alpha_text, all_coordinates)

    def __iter__(self):
        r"""
        Iterate through the stack over a QEPCAD cell.

        EXAMPLES::

            sage: var('x,y')
            (x, y)
            sage: qe = qepcad(x^2 + y^2 == 1, interact=True) # optional - qepcad
            sage: qe.go(); qe.go(); qe.go() # optional - qepcad
            QEPCAD object has moved to phase 'Before Projection (y)'
            QEPCAD object has moved to phase 'Before Choice'
            QEPCAD object has moved to phase 'Before Solution'
            sage: [c.sample_point() for c in qe.cell(3)] # optional - qepcad
            [(0, -3), (0, -1), (0, -1/2), (0, 1), (0, 3)]
        """
        for i in range(1, self._number_of_children + 1):
            yield self[i]

    def __getitem__(self, i):
        r"""
        Select an element from the stack over a QEPCAD cell.

        Note that cells are numbered starting with 1.

        EXAMPLES::

            sage: var('x,y')
            (x, y)
            sage: qe = qepcad(x^2 + y^2 == 1, interact=True)  # optional - qepcad
            sage: qe.go(); qe.go(); qe.go()  # optional - qepcad
            QEPCAD object has moved to phase 'Before Projection (y)'
            QEPCAD object has moved to phase 'Before Choice'
            QEPCAD object has moved to phase 'Before Solution'
            sage: qe.cell(2).__getitem__(3)  # optional - qepcad
            QEPCAD cell (2, 3)
        """
        return self._parent.cell(self.index() + (i,))

    def __len__(self):
        r"""
        Return the number of elements in the stack over a QEPCAD cell.

        This is always an odd number, if the stack has been constructed.

        EXAMPLES::

            sage: var('x,y')
            (x, y)
            sage: qe = qepcad(x^2 + y^2 == 1, interact=True) # optional - qepcad
            sage: qe.go(); qe.go(); qe.go() # optional - qepcad
            QEPCAD object has moved to phase 'Before Projection (y)'
            QEPCAD object has moved to phase 'Before Choice'
            QEPCAD object has moved to phase 'Before Solution'
            sage: len(qe.cell()) # optional - qepcad
            5
            sage: [len(c) for c in qe.cell()] # optional - qepcad
            [1, 3, 5, 3, 1]
        """
        return self._number_of_children

    def __repr__(self):
        r"""
        Return a string representation of a QEPCAD cell.

        Just gives the cell index.

        EXAMPLES::

            sage: var('x,y')
            (x, y)
            sage: qe = qepcad(x^2 + y^2 == 1, interact=True) # optional - qepcad
            sage: qe.go(); qe.go(); qe.go() # optional - qepcad
            QEPCAD object has moved to phase 'Before Projection (y)'
            QEPCAD object has moved to phase 'Before Choice'
            QEPCAD object has moved to phase 'Before Solution'
            sage: qe.cell() # optional - qepcad
            QEPCAD cell ()
            sage: qe.cell(1) # optional - qepcad
            QEPCAD cell (1)
            sage: qe.cell(2, 2) # optional - qepcad
            QEPCAD cell (2, 2)
        """
        ind = self.index()
        if len(ind) == 1:
            ind = '(%s)' % ind[0]
        else:
            ind = str(ind)
        return ('QEPCAD cell %s' % ind)

    def index(self):
        r"""
        Give the index of a QEPCAD cell.

        EXAMPLES::

            sage: var('x,y')
            (x, y)
            sage: qe = qepcad(x^2 + y^2 == 1, interact=True) # optional - qepcad
            sage: qe.go(); qe.go(); qe.go() # optional - qepcad
            QEPCAD object has moved to phase 'Before Projection (y)'
            QEPCAD object has moved to phase 'Before Choice'
            QEPCAD object has moved to phase 'Before Solution'
            sage: qe.cell().index() # optional - qepcad
            ()
            sage: qe.cell(1).index() # optional - qepcad
            (1,)
            sage: qe.cell(2, 2).index() # optional - qepcad
            (2, 2)
        """
        return self._index

    def level(self):
        r"""
        Return the level of a QEPCAD cell.

        EXAMPLES::

            sage: var('x,y')
            (x, y)
            sage: qe = qepcad(x^2 + y^2 == 1, interact=True) # optional - qepcad
            sage: qe.go(); qe.go(); qe.go() # optional - qepcad
            QEPCAD object has moved to phase 'Before Projection (y)'
            QEPCAD object has moved to phase 'Before Choice'
            QEPCAD object has moved to phase 'Before Solution'
            sage: qe.cell().level() # optional - qepcad
            0
            sage: qe.cell(1).level() # optional - qepcad
            1
            sage: qe.cell(2, 2).level() # optional - qepcad
            2
        """
        return self._level

    def signs(self):
        r"""
        Return the sign vector of a QEPCAD cell.

        This is a list of lists.
        The outer list contains one element for each level of the cell;
        the inner list contains one element for each projection factor at
        that level.  These elements are either -1, 0, or 1.

        EXAMPLES::

            sage: var('x,y')
            (x, y)
            sage: qe = qepcad(x^2 + y^2 == 1, interact=True) # optional - qepcad
            sage: qe.go(); qe.go(); qe.go() # optional - qepcad
            QEPCAD object has moved to phase 'Before Projection (y)'
            QEPCAD object has moved to phase 'Before Choice'
            QEPCAD object has moved to phase 'Before Solution'
            sage: from sage.interfaces.qepcad import QepcadCell
            sage: all_cells = flatten(qe.cell(), ltypes=QepcadCell, max_level=1) # optional - qepcad
            sage: [(c, c.signs()[1][0]) for c in all_cells] # optional - qepcad
            [(QEPCAD cell (1, 1), 1), (QEPCAD cell (2, 1), 1), (QEPCAD cell (2, 2), 0), (QEPCAD cell (2, 3), 1), (QEPCAD cell (3, 1), 1), (QEPCAD cell (3, 2), 0), (QEPCAD cell (3, 3), -1), (QEPCAD cell (3, 4), 0), (QEPCAD cell (3, 5), 1), (QEPCAD cell (4, 1), 1), (QEPCAD cell (4, 2), 0), (QEPCAD cell (4, 3), 1), (QEPCAD cell (5, 1), 1)]
        """
        return self._signs

    def number_of_children(self):
        r"""
        Return the number of elements in the stack over a QEPCAD cell.
        (This is always an odd number, if the stack has been constructed.)

        EXAMPLES::

            sage: var('x,y')
            (x, y)
            sage: qe = qepcad(x^2 + y^2 == 1, interact=True) # optional - qepcad
            sage: qe.go(); qe.go(); qe.go() # optional - qepcad
            QEPCAD object has moved to phase 'Before Projection (y)'
            QEPCAD object has moved to phase 'Before Choice'
            QEPCAD object has moved to phase 'Before Solution'
            sage: qe.cell().number_of_children() # optional - qepcad
            5
            sage: [c.number_of_children() for c in qe.cell()] # optional - qepcad
            [1, 3, 5, 3, 1]
        """
        return self._number_of_children

    def set_truth(self, v):
        r"""
        Set the truth value of this cell, as used by QEPCAD for solution
        formula construction.

        The argument ``v`` should be either a boolean or ``None``
        (which will set the truth value to ``'undetermined'``).

        EXAMPLES::

            sage: var('x,y')
            (x, y)
            sage: qe = qepcad(x^2 + y^2 == 1, interact=True) # optional - qepcad
            sage: qe.go(); qe.go(); qe.go() # optional - qepcad
            QEPCAD object has moved to phase 'Before Projection (y)'
            QEPCAD object has moved to phase 'Before Choice'
            QEPCAD object has moved to phase 'Before Solution'
            sage: qe.solution_extension('T') # optional - qepcad
            y^2 + x^2 - 1 = 0
            sage: qe.cell(3, 3).set_truth(True) # optional - qepcad
            sage: qe.solution_extension('T') # optional - qepcad
            y^2 + x^2 - 1 <= 0
        """
        if v is None:
            nv = 2
        else:
            nv = 1 if v else 0
        self._parent.set_truth_value(self, nv)

    def sample_point(self):
        r"""
        Return the coordinates of a point in the cell, as a tuple of
        \sage algebraic reals.

        EXAMPLES::

            sage: qe = qepcad(x^2 - x - 1 == 0, interact=True) # optional - qepcad
            sage: qe.go(); qe.go(); qe.go() # optional - qepcad
            QEPCAD object has moved to phase 'At the end of projection phase'
            QEPCAD object has moved to phase 'Before Choice'
            QEPCAD object has moved to phase 'Before Solution'
            sage: v1 = qe.cell(2).sample_point()[0]; v1 # optional - qepcad
            -0.618033988749895?
            sage: v2 = qe.cell(4).sample_point()[0]; v2 # optional - qepcad
            1.618033988749895?
            sage: v1 + v2 == 1 # optional - qepcad
            True
        """
        try:
            return self._sample_point
        except AttributeError:
            (extended, alpha, coordinates) = self._raw_sample_point

            need_alpha = False
            for c in coordinates:
                if 'alpha' in c:
                    need_alpha = True
                    break

            if need_alpha:
                locals = {'alpha': _eval_qepcad_algebraic(alpha)}
            else:
                locals = {}

            points = []

            for c in coordinates:
                if 'unique' in c:
                    points.append(_eval_qepcad_algebraic(c))
                else:
                    points.append(sage_eval(implicit_mul(c), locals))

            return tuple(points)

    def sample_point_dict(self):
        r"""
        Return the coordinates of a point in the cell, as a dictionary
        mapping variable names (as strings) to \sage algebraic reals.

        EXAMPLES::

            sage: qe = qepcad(x^2 - x - 1 == 0, interact=True) # optional - qepcad
            sage: qe.go(); qe.go(); qe.go() # optional - qepcad
            QEPCAD object has moved to phase 'At the end of projection phase'
            QEPCAD object has moved to phase 'Before Choice'
            QEPCAD object has moved to phase 'Before Solution'
            sage: qe.cell(4).sample_point_dict() # optional - qepcad
            {'x': 1.618033988749895?}
        """
        points = self.sample_point()
        vars = self._parent._varlist

        return dict([(vars[i], points[i]) for i in range(len(points))])

