.. -*- coding: utf-8 -*-
.. _group_theory:

=====================
Group Theory and Sage
=====================

.. MODULEAUTHOR:: Robert A. Beezer, University of Puget Sound

This compilation collects Sage commands that are useful for a student
in an introductory course on group theory.  It is not intended to
teach Sage or to teach group theory.  (There are many introductory
texts on group theory and more information on Sage can be found via
`<www.sagemath.org>`_) Rather, by presenting commands roughly in the
order a student would learn the corresponding mathematics they might
be encouraged to experiment and learn more about mathematics and learn
more about Sage.  Not coincidentally, when Sage was the acronym SAGE,
the "E" in Sage stood for "Experimentation."

This guide is also distributed in PDF format and as a Sage worksheet.
The worksheet version can be imported into the Sage notebook
environment running in a web browser, and then the displayed chunks of
code may be executed by Sage if one clicks on the small "evaluate"
link below each cell, for a fully interactive experience. A PDF and
Sage worksheet versions of this tutorial are available at
http://abstract.ups.edu/sage-aata.html.

.. contents:: Table of contents
   :depth: 2
   :class: this-will-duplicate-information-and-it-is-still-useful-here

Changelog:

* 2009/01/30  Version 1.0, first complete release
* 2009/03/03  Version 1.1, added cyclic group size interact
* 2010/03/10  Version 1.3, dropped US on license, some edits.

Basic properties of the integers
================================

Integer division
----------------

The command ``a % b`` will return the remainder upon division of `a`
by `b`.  In other words, the value is the unique integer `r` such that:

#. `0 \leq r < b`; and
#. `a = bq + r` for some integer `q` (the quotient).

Then `(a - r) / b` will equal `q`.  For example::

    sage: r = 14 % 3
    sage: q = (14 - r) / 3
    sage: r, q
    (2, 4)

will return ``2`` for the value of ``r`` and ``4`` for the value of
``q``.  Note that the "``/``" is *integer* division, where any
remainder is cast away and the result is always an integer.  So, for
example, ``14 / 3`` will again equal ``4``, not ``4.66666``.


Greatest common divisor
-----------------------

The greatest common divisor of `a` and `b` is obtained with the
command ``gcd(a,b)``, where in our first uses, `a` and `b` are
integers.  Later, `a` and `b` can be other objects with a notion
of divisibility and "greatness," such as polynomials.  For example::

    sage: gcd(2776, 2452)
    4


Extended greatest common divisor
--------------------------------

The command ``xgcd(a, b)`` ("eXtended GCD") returns a triple where the
first element is the greatest common divisor of `a` and `b` (as with
the ``gcd(a, b)`` command above), but the next two elements are the
values of `r` and `s` such that `ra + sb = \gcd(a, b)`. For example,
``xgcd(633, 331)`` returns ``(1, 194, -371)``.  Portions of the triple
can be extracted using ``[ ]`` to access the entries of the triple,
starting with the first as number ``0``.  For example, the following
should return the result ``True`` (even if you change the values of
``a`` and ``b``).  Studying this block of code will go a long way
towards helping you get the most out of Sage's output. (Note that
"``=``" is how a value is assigned to a variable, while as in the last
line, "``==``" is how we determine equality of two items.) ::

    sage: a = 633
    sage: b = 331
    sage: extended = xgcd(a, b)
    sage: g = extended[0]
    sage: r = extended[1]
    sage: s = extended[2]
    sage: g == r*a + s*b
    True


Divisibility
------------

A remainder of zero indicates divisibility.  So ``(a % b) == 0`` will
return ``True`` if `b` divides `a`, and will otherwise return
``False``.  For example, ``(9 % 3) == 0`` is ``True``, but
``(9 % 4) == 0`` is ``False``.  Try predicting the output of the
following before executing it in Sage. ::

    sage: answer1 = ((20 % 5) == 0)
    sage: answer2 = ((17 % 4) == 0)
    sage: answer1, answer2
    (True, False)


Factoring
---------

As promised by the Fundamental Theorem of Arithmetic, ``factor(a)``
will return a unique expression for `a` as a product of powers of
primes.  It will print in a nicely-readable form, but can also be
manipulated with Python as a list of pairs `(p_i, e_i)` containing
primes as bases, and their associated exponents. For example::

    sage: factor(2600)
    2^3 * 5^2 * 13

If you just want the prime divisors of an integer, then use the
``prime_divisors(a)`` command, which will return a list of all the
prime divisors of `a`.  For example::

    sage: prime_divisors(2600)
    [2, 5, 13]

We can strip off other pieces of the prime decomposition using two
levels of ``[ ]``.  This is another good example to study in order to
learn about how to drill down into Python lists. ::

    sage: n = 2600
    sage: decomposition = factor(n)
    sage: print("{} decomposes as {}".format(n, decomposition))
    2600 decomposes as 2^3 * 5^2 * 13
    sage: secondterm = decomposition[1]
    sage: print("Base and exponent (pair) for second prime: "+str(secondterm))
    Base and exponent (pair) for second prime: (5, 2)
    sage: base = secondterm[0]
    sage: exponent = secondterm[1]
    sage: print("Base is "+str(base))
    Base is 5
    sage: print("Exponent is "+str(exponent))
    Exponent is 2
    sage: thirdbase = decomposition[2][0]
    sage: thirdexponent = decomposition[2][1]
    sage: print("Base of third term is {} with exponent {}".format(thirdbase, thirdexponent))
    Base of third term is 13 with exponent  1

With a bit more work, the ``factor()`` command can be used to factor
more complicated items, such as polynomials.


Multiplicative inverse, modular arithmetic
------------------------------------------

The command ``inverse_mod(a, n)`` yields the multiplicative inverse of
`a` mod `n` (or an error if it doesn't exist).  For example::

    sage: inverse_mod(352, 917)
    508

(As a check, find the integer `m` such that ``352*508 == m*917+1``.)
Then try ::

    sage: inverse_mod(4, 24)
    Traceback (most recent call last):
    ...
    ZeroDivisionError: inverse of Mod(4, 24) does not exist

and explain the result.


Powers with modular arithmetic
------------------------------

The command ``power_mod(a, m, n)`` yields `a^m` mod `n`.  For example::

    sage: power_mod(15, 831, 23)
    10

If `m = -1`, then this command will duplicate the function of
``inverse_mod()``.


Euler `\phi`-function
---------------------

The command ``euler_phi(n)`` will return the number of positive
integers less than `n` and relatively prime to `n` (i.e. having
greatest common divisor with `n` equal to 1).  For example::

    sage: euler_phi(345)
    176

Experiment by running the following code several times::

    sage: m = random_prime(10000)
    sage: n = random_prime(10000)
    sage: euler_phi(m*n) == euler_phi(m) * euler_phi(n) or m == n
    True

Feel a conjecture coming on?  Can you generalize this result?


Primes
------

The command ``is_prime(a)`` returns ``True`` or ``False`` depending on
if `a` is prime or not.  For example, ::

    sage: is_prime(117371)
    True

while ::

    sage: is_prime(14547073)
    False

since `14547073 = 1597 * 9109` (as you could determine with the
``factor()`` command).

The command ``random_prime(a, True)`` will return a random prime
between 2 and `a`. Experiment with::

    sage: p = random_prime(10^21, True)
    sage: is_prime(p)
    True

(Replacing ``True`` by ``False`` will speed up the search, but there
will be a very small probability the result will not be prime.)

The command ``prime_range(a, b)`` returns an ordered list of all the
primes from `a` to `b - 1`, inclusive.  For example, ::

    sage: prime_range(500, 550)
    [503, 509, 521, 523, 541, 547]

The commands ``next_prime(a)`` and ``previous_prime(a)`` are other
ways to get a single prime number of a desired size.  Give them a try.


Permutation groups
==================

A good portion of Sage's support for group theory is based on routines
from GAP (Groups, Algorithms, and Programming at
https://www.gap-system.org.  Groups can be described in many different
ways, such as sets of matrices or sets of symbols subject to a few
defining relations.  A very concrete way to represent groups is via
permutations (one-to-one and onto functions of the integers 1
through `n`), using function composition as the operation in the
group.  Sage has many routines designed to work with groups of this
type and they are also a good way for those learning group theory to
gain experience with the basic ideas of group theory.  For both these
reasons, we will concentrate on these types of groups.


Writing permutations
--------------------

Sage uses "disjoint cycle notation" for permutations, see any
introductory text on group theory (such as Judson, Section 4.1) for
more on this.  Composition occurs *left to right*, which is not what
you might expect and is exactly the reverse of what Judson and many
others use.  (There are good reasons to support either direction, you
just need to be certain you know which one is in play.)  There are two
ways to write the permutation `\sigma = (1\,3) (2\,5\,4)`:

#. As a text string (include quotes): ``"(1,3) (2,5,4)"``
#. As a Python list of "tuples": ``[(1,3), (2,5,4)]``


Groups
------

Sage knows many popular groups as sets of permutations.  More are
listed below, but for starters, the full "symmetric group" of all
possible permutations of 1 through `n` can be built with the command
``SymmetricGroup(n)``.

**Permutation elements** Elements of a group can be created, and
composed, as follows ::

    sage: G = SymmetricGroup(5)
    sage: sigma = G("(1,3) (2,5,4)")
    sage: rho = G([(2,4), (1,5)])
    sage: rho^(-1) * sigma * rho
    (1,2,4)(3,5)

Available functions for elements of a permutation group include
finding the order of an element, i.e. for a permutation `\sigma` the
order is the smallest power of `k` such that `\sigma^k` equals the
identity element `()`.  For example::

    sage: G = SymmetricGroup(5)
    sage: sigma = G("(1,3) (2,5,4)")
    sage: sigma.order()
    6

The sign of the permutation `\sigma` is defined to be 1 for an even
permutation and `-1` for an odd permutation.  For example::

    sage: G = SymmetricGroup(5)
    sage: sigma = G("(1,3) (2,5,4)")
    sage: sigma.sign()
    -1

since `\sigma` is an odd permutation.

Many more available functions that can be applied to a permutation can
be found via "tab-completion."  With ``sigma`` defined as an element
of a permutation group, in a Sage cell, type ``sigma.`` (Note the
"``.``") and then press the :kbd:`Tab` key.  You will get a list of available
functions (you may need to scroll down to see the whole list).
Experiment and explore!  It is what Sage is all about.  You really
cannot break anything.


Creating groups
---------------

This is an annotated list of some small well-known permutation groups
that can be created simply in Sage. You can find more in the source
code file :sage_root:`src/sage/groups/perm_gps/permgroup_named.py`.

* ``SymmetricGroup(n)``: All `n!` permutations on `n` symbols.
* ``DihedralGroup(n)``: Symmetries of an `n`-gon.  Rotations and
  flips, `2n` in total.
* ``CyclicPermutationGroup(n)``: Rotations of an `n`-gon (no flips),
  `n` in total.
* ``AlternatingGroup(n)``:  Alternating group on `n` symbols having
  `n!/2` elements.
* ``KleinFourGroup()``:  The non-cyclic group of order 4.


Group functions
===============

Individual elements of permutation groups are important, but we
primarily wish to study groups as objects on their own.  So a
wide variety of computations are available for groups. Define a group,
for example ::

    sage: H = DihedralGroup(6)
    sage: H
    Dihedral group of order 12 as a permutation group

and then a variety of functions become available.

After trying the examples below, experiment with tab-completion.
Having defined ``H``, type ``H.`` (note the "``.``") and then press
the :kbd:`Tab` key.  You will get a list of available functions (you may need
to scroll down to see the whole list).  As before,
*experiment and explore*---it is really hard to break anything.

Here is another couple of ways to experiment and explore.  Find a
function that looks interesting, say ``is_abelian()``.  Type
``H.is_abelian?`` (note the question mark) followed by the :kbd:`Enter` key.
This will display a portion of the source code for the
``is_abelian()`` function, describing the inputs and output, possibly
illustrated with example uses.

If you want to learn more about how Sage works, or possibly extend its
functionality, then you can start by examining the complete Python
source code.  For example, try ``H.is_abelian??``, which will allow
you to determine that the ``is_abelian()`` function is basically
riding on GAP's ``IsAbelian()`` command and asking GAP do the
heavy-lifting for us.  (To get the maximum advantage of using Sage it
helps to know some basic Python programming, but it is not required.)

OK, on to some popular command for groups.  If you are using the
worksheet, be sure you have defined the group `H` as the dihedral
group `D_6`, since we will not keep repeating its definition below.


Abelian?
--------

The command ::

    sage: H = DihedralGroup(6)
    sage: H.is_abelian()
    False

will return ``False`` since `D_6` is a non-abelian group.


Order
-----

The command ::

    sage: H = DihedralGroup(6)
    sage: H.order()
    12

will return ``12`` since `D_6` is a group of with 12 elements.


All elements
------------

The command ::

    sage: H = DihedralGroup(6)
    sage: H.list()
    [(),
     (1,5,3)(2,6,4),
     (1,3,5)(2,4,6),
     (1,6,5,4,3,2),
     (1,4)(2,5)(3,6),
     (1,2,3,4,5,6),
     (2,6)(3,5),
     (1,5)(2,4),
     (1,3)(4,6),
     (1,6)(2,5)(3,4),
     (1,4)(2,3)(5,6),
     (1,2)(3,6)(4,5)]

will return all of the elements of `H` in a fixed order as a Python
list.  Indexing (``[ ]``) can be used to extract the individual
elements of the list, remembering that counting the elements of the
list begins at zero. ::

    sage: H = DihedralGroup(6)
    sage: elements = H.list()
    sage: elements[2]
    (1,3,5)(2,4,6)


Cayley table
------------

The command ::

    sage: H = DihedralGroup(6)
    sage: H.cayley_table()
    *  a b c d e f g h i j k l
     +------------------------
    a| a b c d e f g h i j k l
    b| b a d c f e h g j i l k
    c| c k a e d g f i h l b j
    d| d l b f c h e j g k a i
    e| e j k g a i d l f b c h
    f| f i l h b j c k e a d g
    g| g h j i k l a b d c e f
    h| h g i j l k b a c d f e
    i| i f h l j b k c a e g d
    j| j e g k i a l d b f h c
    k| k c e a g d i f l h j b
    l| l d f b h c j e k g i a


will construct the Cayley table (or "multiplication table") of `H`.
By default the table uses lowercase Latin letters to name the elements
of the group.  The actual elements used can be found using the
``row_keys()`` or ``column_keys()`` commands for the table.
For example to determine the fifth element in the table, the
element named ``e``::

    sage: H = DihedralGroup(6)
    sage: T = H.cayley_table()
    sage: headings = T.row_keys()
    sage: headings[4]
    (1,3)(4,6)

Center
------

The command ``H.center()`` will return a subgroup that is the center
of the group `H` (see Exercise 2.46 in Judson).  Try ::

    sage: H = DihedralGroup(6)
    sage: H.center().list()
    [(), (1,4)(2,5)(3,6)]

to see which elements of `H` commute with *every* element of `H`.


Cayley graph
------------

For fun, try ``show(H.cayley_graph())``.


Subgroups
=========


Cyclic subgroups
----------------

If ``G`` is a group and ``a`` is an element of the group (try
``a = G.random_element()``), then ::

    a = G.random_element()
    H = G.subgroup([a])

will create ``H`` as the cyclic subgroup of ``G`` with generator
``a``.

For example the code below will:

#. create ``G`` as the symmetric group on five symbols;
#. specify ``sigma`` as an element of ``G``;
#. use ``sigma`` as the generator of a cyclic subgroup ``H``;
#. list all the elements of ``H``.

In more mathematical notation, we might write
`\langle (1\,2\,3) (4\,5) \rangle = H \subseteq G = S_5`. ::

    sage: G = SymmetricGroup(5)
    sage: sigma = G("(1,2,3) (4,5)")
    sage: H = G.subgroup([sigma])
    sage: H.list()
    [(), (1,2,3)(4,5), (1,3,2), (4,5), (1,2,3), (1,3,2)(4,5)]

Experiment by trying different permutations for ``sigma`` and
observing the effect on ``H``.


Cyclic groups
-------------

Groups that are cyclic themselves are both important and rich in
structure.  The command ``CyclicPermutationGroup(n)`` will create a
permutation group that is cyclic with ``n`` elements.  Consider the
following example (note that the indentation of the third line is
critical) which will list the elements of a cyclic group of order 20,
preceded by the order of each element. ::

    sage: n = 20
    sage: CN = CyclicPermutationGroup(n)
    sage: for g in CN:
    ....:     print("{}   {}".format(g.order(), g))
    1   ()
    20   (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)
    10   (1,3,5,7,9,11,13,15,17,19)(2,4,6,8,10,12,14,16,18,20)
    20   (1,4,7,10,13,16,19,2,5,8,11,14,17,20,3,6,9,12,15,18)
    5   (1,5,9,13,17)(2,6,10,14,18)(3,7,11,15,19)(4,8,12,16,20)
    4   (1,6,11,16)(2,7,12,17)(3,8,13,18)(4,9,14,19)(5,10,15,20)
    10   (1,7,13,19,5,11,17,3,9,15)(2,8,14,20,6,12,18,4,10,16)
    20   (1,8,15,2,9,16,3,10,17,4,11,18,5,12,19,6,13,20,7,14)
    5   (1,9,17,5,13)(2,10,18,6,14)(3,11,19,7,15)(4,12,20,8,16)
    20   (1,10,19,8,17,6,15,4,13,2,11,20,9,18,7,16,5,14,3,12)
    2   (1,11)(2,12)(3,13)(4,14)(5,15)(6,16)(7,17)(8,18)(9,19)(10,20)
    20   (1,12,3,14,5,16,7,18,9,20,11,2,13,4,15,6,17,8,19,10)
    5   (1,13,5,17,9)(2,14,6,18,10)(3,15,7,19,11)(4,16,8,20,12)
    20   (1,14,7,20,13,6,19,12,5,18,11,4,17,10,3,16,9,2,15,8)
    10   (1,15,9,3,17,11,5,19,13,7)(2,16,10,4,18,12,6,20,14,8)
    4   (1,16,11,6)(2,17,12,7)(3,18,13,8)(4,19,14,9)(5,20,15,10)
    5   (1,17,13,9,5)(2,18,14,10,6)(3,19,15,11,7)(4,20,16,12,8)
    20   (1,18,15,12,9,6,3,20,17,14,11,8,5,2,19,16,13,10,7,4)
    10   (1,19,17,15,13,11,9,7,5,3)(2,20,18,16,14,12,10,8,6,4)
    20   (1,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2)

By varying the size of the group (change the value of ``n``) you can
begin to illustrate some of the structure of a cyclic group (for
example, try a prime).

We can cut/paste an element of order 5 from the output above (in the
case when the cyclic group has 20 elements) and quickly build a
subgroup::

    sage: C20 = CyclicPermutationGroup(20)
    sage: rho = C20("(1,17,13,9,5)(2,18,14,10,6)(3,19,15,11,7)(4,20,16,12,8)")
    sage: H = C20.subgroup([rho])
    sage: H.list()
    [(),
     (1,17,13,9,5)(2,18,14,10,6)(3,19,15,11,7)(4,20,16,12,8),
     (1,13,5,17,9)(2,14,6,18,10)(3,15,7,19,11)(4,16,8,20,12),
     (1,9,17,5,13)(2,10,18,6,14)(3,11,19,7,15)(4,12,20,8,16),
     (1,5,9,13,17)(2,6,10,14,18)(3,7,11,15,19)(4,8,12,16,20)]

For a cyclic group, the following command will list *all* of the
subgroups. ::

    sage: C20 = CyclicPermutationGroup(20)
    sage: C20.conjugacy_classes_subgroups()
    [Subgroup generated by [()] of (Cyclic group of order 20 as a permutation group),
     Subgroup generated by [(1,11)(2,12)(3,13)(4,14)(5,15)(6,16)(7,17)(8,18)(9,19)(10,20)] of (Cyclic group of order 20 as a permutation group),
     Subgroup generated by [(1,6,11,16)(2,7,12,17)(3,8,13,18)(4,9,14,19)(5,10,15,20), (1,11)(2,12)(3,13)(4,14)(5,15)(6,16)(7,17)(8,18)(9,19)(10,20)] of (Cyclic group of order 20 as a permutation group),
     Subgroup generated by [(1,5,9,13,17)(2,6,10,14,18)(3,7,11,15,19)(4,8,12,16,20)] of (Cyclic group of order 20 as a permutation group),
     Subgroup generated by [(1,5,9,13,17)(2,6,10,14,18)(3,7,11,15,19)(4,8,12,16,20), (1,3,5,7,9,11,13,15,17,19)(2,4,6,8,10,12,14,16,18,20)] of (Cyclic group of order 20 as a permutation group),
     Subgroup generated by [(1,5,9,13,17)(2,6,10,14,18)(3,7,11,15,19)(4,8,12,16,20), (1,3,5,7,9,11,13,15,17,19)(2,4,6,8,10,12,14,16,18,20), (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)] of (Cyclic group of order 20 as a permutation group)]

Be careful, this command uses some more advanced ideas and will not
usually list *all* of the subgroups of a group. Here we are relying on
special properties of cyclic groups (but see the next section).

If you are viewing this as a PDF, you can safely skip over the next
bit of code.  However, if you are viewing this as a worksheet in Sage,
then this is a place where you can experiment with the structure of
the subgroups of a cyclic group.  In the input box, enter the order of
a cyclic group (numbers between 1 and 40 are good initial choices) and
Sage will list each subgroup as a cyclic group with its generator.
The factorization at the bottom might help you formulate a
conjecture. ::

    %auto
    @interact
    def _(n = input_box(default=12, label = "Cyclic group of order:", type=Integer) ):
        cyclic = CyclicPermutationGroup(n)
        subgroups = cyclic.conjugacy_classes_subgroups()
        html( "All subgroups of a cyclic group of order $%s$\n" % latex(n) )
        table = "$\\begin{array}{ll}"
        for sg in subgroups:
          table = table + latex(sg.order()) + \
                  " & \\left\\langle" + latex(sg.gens()[0]) + \
                  "\\right\\rangle\\\\"
        table = table + "\\end{array}$"
        html(table)
        html("\nHint: $%s$ factors as $%s$" % ( latex(n), latex(factor(n)) ) )


All subgroups
-------------

If `H` is a subgroup of `G` and `g \in G`, then
`gHg^{-1} = \{ghg^{-1} \mid h \in G\}` will also be a subgroup of `G`.
If ``G`` is a group, then the command
``G.conjugacy_classes_subgroups()`` will return a list of subgroups of
``G``, but not all of the subgroups.  However, every subgroup can be
constructed from one on the list by the `gHg^{-1}` construction with a
suitable `g`. As an illustration, the code below:

#. creates ``K`` as the dihedral group of order 24, `D_{12}`;
#. stores the list of subgroups output by
   ``K.conjugacy_classes_subgroups()`` in the variable ``sg``;
#. prints the elements of the list;
#. selects the second subgroup in the list, and lists its elements.

::

    sage: K = DihedralGroup(12)
    sage: sg = K.conjugacy_classes_subgroups()
    sage: sg
    [Subgroup generated by [()] of (Dihedral group of order 24 as a permutation group),
     Subgroup generated by [(1,7)(2,8)(3,9)(4,10)(5,11)(6,12)] of (Dihedral group of order 24 as a permutation group),
     Subgroup generated by [(2,12)(3,11)(4,10)(5,9)(6,8)] of (Dihedral group of order 24 as a permutation group),
     Subgroup generated by [(1,2)(3,12)(4,11)(5,10)(6,9)(7,8)] of (Dihedral group of order 24 as a permutation group),
     Subgroup generated by [(1,5,9)(2,6,10)(3,7,11)(4,8,12)] of (Dihedral group of order 24 as a permutation group),
     Subgroup generated by [(2,12)(3,11)(4,10)(5,9)(6,8), (1,7)(2,8)(3,9)(4,10)(5,11)(6,12)] of (Dihedral group of order 24 as a permutation group),
     Subgroup generated by [(1,10,7,4)(2,11,8,5)(3,12,9,6), (1,7)(2,8)(3,9)(4,10)(5,11)(6,12)] of (Dihedral group of order 24 as a permutation group),
     Subgroup generated by [(1,2)(3,12)(4,11)(5,10)(6,9)(7,8), (1,7)(2,8)(3,9)(4,10)(5,11)(6,12)] of (Dihedral group of order 24 as a permutation group),
     Subgroup generated by [(1,5,9)(2,6,10)(3,7,11)(4,8,12), (1,3,5,7,9,11)(2,4,6,8,10,12)] of (Dihedral group of order 24 as a permutation group),
     Subgroup generated by [(1,5,9)(2,6,10)(3,7,11)(4,8,12), (2,12)(3,11)(4,10)(5,9)(6,8)] of (Dihedral group of order 24 as a permutation group),
     Subgroup generated by [(1,5,9)(2,6,10)(3,7,11)(4,8,12), (1,2)(3,12)(4,11)(5,10)(6,9)(7,8)] of (Dihedral group of order 24 as a permutation group),
     Subgroup generated by [(2,12)(3,11)(4,10)(5,9)(6,8), (1,10,7,4)(2,11,8,5)(3,12,9,6), (1,7)(2,8)(3,9)(4,10)(5,11)(6,12)] of (Dihedral group of order 24 as a permutation group),
     Subgroup generated by [(1,5,9)(2,6,10)(3,7,11)(4,8,12), (1,3,5,7,9,11)(2,4,6,8,10,12), (2,12)(3,11)(4,10)(5,9)(6,8)] of (Dihedral group of order 24 as a permutation group),
     Subgroup generated by [(1,5,9)(2,6,10)(3,7,11)(4,8,12), (1,3,5,7,9,11)(2,4,6,8,10,12), (1,2,3,4,5,6,7,8,9,10,11,12)] of (Dihedral group of order 24 as a permutation group),
     Subgroup generated by [(1,5,9)(2,6,10)(3,7,11)(4,8,12), (1,3,5,7,9,11)(2,4,6,8,10,12), (1,2)(3,12)(4,11)(5,10)(6,9)(7,8)] of (Dihedral group of order 24 as a permutation group),
     Subgroup generated by [(1,5,9)(2,6,10)(3,7,11)(4,8,12), (1,3,5,7,9,11)(2,4,6,8,10,12), (2,12)(3,11)(4,10)(5,9)(6,8), (1,2,3,4,5,6,7,8,9,10,11,12)] of (Dihedral group of order 24 as a permutation group)]

    sage: print("An order two subgroup:\n{}".format(sg[1].list()))
    An order two subgroup:
    [(), (1,7)(2,8)(3,9)(4,10)(5,11)(6,12)]

It is important to note that this is a nice long list of subgroups,
but will rarely create *every* such subgroup.  For example, the
code below:

#. creates ``rho`` as an element of the group ``K``;
#. creates ``L`` as a cyclic subgroup of ``K``;
#. prints the two elements of ``L``; and finally
#. tests to see if this subgroup is part of the output of the list
   ``sg`` created just above (it is not).

::

    sage: K = DihedralGroup(12)
    sage: sg = K.conjugacy_classes_subgroups()
    sage: rho = K("(1,4) (2,3) (5,12) (6,11) (7,10) (8,9)")
    sage: L = PermutationGroup([rho])
    sage: L.list()
    [(), (1,4)(2,3)(5,12)(6,11)(7,10)(8,9)]
    sage: L in sg
    False


Symmetry groups
===============

You can give Sage a short list of elements of a permutation group and
Sage will find the smallest subgroup that contains those elements.  We
say the list "generates" the subgroup.  We list a few interesting
subgroups you can create this way.


Symmetries of an equilateral triangle
-------------------------------------

Label the vertices of an equilateral triangle as 1, 2 and 3.  Then
*any* permutation of the vertices will be a symmetry of the triangle.
So either ``SymmetricGroup(3)`` or ``DihedralGroup(3)`` will create
the full symmetry group.


Symmetries of an `n`-gon
------------------------

A regular, `n`-sided figure in the plane (an `n`-gon) has `2n`
symmetries, comprised of `n` rotations (including the trivial one) and
`n` "flips" about various axes.  The dihedral group
``DihedralGroup(n)`` is frequently defined as exactly the symmetry
group of an `n`-gon.


Symmetries of a tetrahedron
---------------------------

Label the 4 vertices of a regular tetrahedron as 1, 2, 3 and 4.  Fix
the vertex labeled 4 and rotate the opposite face through 120 degrees.
This will create the permutation/symmetry `(1\,2\, 3)`.  Similarly,
fixing vertex 1, and rotating the opposite face will create the
permutation `(2\,3\,4)`.  These two permutations are enough to
generate the full group of the twelve symmetries of the tetrahedron.
Another symmetry can be visualized by running an axis through the
midpoint of an edge of the tetrahedron through to the midpoint of the
opposite edge, and then rotating by 180 degrees about this axis.  For
example, the 1--2 edge is opposite the 3--4 edge, and the symmetry is
described by the permutation `(1\,2) (3\,4)`.  This permutation, along
with either of the above permutations will also generate the group.
So here are two ways to create this group::

    sage: tetra_one = PermutationGroup(["(1,2,3)", "(2,3,4)"])
    sage: tetra_one
    Permutation Group with generators [(2,3,4), (1,2,3)]
    sage: tetra_two = PermutationGroup(["(1,2,3)", "(1,2)(3,4)"])
    sage: tetra_two
    Permutation Group with generators [(1,2)(3,4), (1,2,3)]

This group has a variety of interesting properties, so it is worth
experimenting with.  You may also know it as the "alternating group
on 4 symbols," which Sage will create with the command
``AlternatingGroup(4)``.


Symmetries of a cube
--------------------

Label vertices of one face of a cube with 1, 2, 3 and 4, and on the
opposite face label the vertices 5, 6, 7 and 8 (5 opposite 1, 6
opposite 2, etc.).  Consider three axes that run from the center of a
face to the center of the opposite face, and consider a quarter-turn
rotation about each axis.  These three rotations will construct the
entire symmetry group.  Use ::

    sage: cube = PermutationGroup(["(3,2,6,7)(4,1,5,8)",
    ....:     "(1,2,6,5)(4,3,7,8)", "(1,2,3,4)(5,6,7,8)"])
    sage: cube.list()
    [(),
     (1,3)(2,4)(5,7)(6,8),
     (1,6)(2,5)(3,8)(4,7),
     (1,8)(2,7)(3,6)(4,5),
     (1,4,3,2)(5,8,7,6),
     (1,2,3,4)(5,6,7,8),
     (1,5)(2,8)(3,7)(4,6),
     (1,7)(2,6)(3,5)(4,8),
     (2,5,4)(3,6,8),
     (1,3,8)(2,7,5),
     (1,6,3)(4,5,7),
     (1,8,6)(2,4,7),
     (1,4)(2,8)(3,5)(6,7),
     (1,2,6,5)(3,7,8,4),
     (1,5,6,2)(3,4,8,7),
     (1,7)(2,3)(4,6)(5,8),
     (2,4,5)(3,8,6),
     (1,3,6)(4,7,5),
     (1,6,8)(2,7,4),
     (1,8,3)(2,5,7),
     (1,4,8,5)(2,3,7,6),
     (1,2)(3,5)(4,6)(7,8),
     (1,5,8,4)(2,6,7,3),
     (1,7)(2,8)(3,4)(5,6)]

A cube has four distinct diagonals (joining opposite vertices through
the center of the cube).  Each symmetry of the cube will cause the
diagonals to arrange differently. In this way, we can view an element
of the symmetry group as a permutation of four "symbols"---the
diagonals.  It happens that *each* of the 24 permutations of the
diagonals is created by exactly one symmetry of the 8 vertices of the
cube.  So this subgroup of `S_8` is "the same as" `S_4`.  In Sage::

    sage: cube = PermutationGroup(["(3,2,6,7)(4,1,5,8)",
    ....:     "(1,2,6,5)(4,3,7,8)", "(1,2,3,4)(5,6,7,8)"])
    sage: cube.is_isomorphic(SymmetricGroup(4))
    True

will test to see if the group of symmetries of the cube are "the same
as" `S_4` and so will return ``True``.

Here is an another way to create the symmetries of a cube.  Number the
six *faces* of the cube as follows:  1 on top, 2 on the bottom, 3
in front, 4 on the right, 5 in back, 6 on the left.  Now the same
rotations as before (quarter-turns about axes through the centers of
two opposite faces) can be used as generators of the symmetry group::

    sage: cubeface = PermutationGroup(["(1,3,2,5)", "(1,4,2,6)", "(3,4,5,6)"])
    sage: cubeface.list()
    [(),
     (3,5)(4,6),
     (1,6,5)(2,4,3),
     (1,6,3)(2,4,5),
     (1,5,6)(2,3,4),
     (1,5,4)(2,3,6),
     (1,2)(4,6),
     (1,2)(3,5),
     (1,4,5)(2,6,3),
     (1,4,3)(2,6,5),
     (1,3,4)(2,5,6),
     (1,3,6)(2,5,4),
     (1,4,2,6),
     (1,6)(2,4)(3,5),
     (1,3,2,5),
     (1,5,2,3),
     (1,2)(3,4)(5,6),
     (3,6,5,4),
     (1,6,2,4),
     (1,4)(2,6)(3,5),
     (1,5)(2,3)(4,6),
     (1,3)(2,5)(4,6),
     (3,4,5,6),
     (1,2)(3,6)(4,5)]

Again, this subgroup of `S_6` is "same as" the full symmetric group, `S_4`::

    sage: cubeface = PermutationGroup(["(1,3,2,5)", "(1,4,2,6)", "(3,4,5,6)"])
    sage: cubeface.is_isomorphic(SymmetricGroup(4))
    True

It turns out that in each of the above constructions, it is sufficient
to use just two of the three generators (any two).  But one generator
is not enough.  Give it a try and use Sage to convince yourself that
a generator can be sacrificed in each case.


Normal subgroups
================


Checking normality
------------------

The code below:

#. begins with the alternating group `A_4`;
#. specifies three elements of the group (the three symmetries of the
   tetrahedron that are 180 degree rotations about axes through midpoints
   of opposite edges);
#. uses these three elements to generate a subgroup; and finally
#. illustrates the command for testing if the subgroup ``H`` is a
   normal subgroup of the group ``A4``.

::

    sage: A4 = AlternatingGroup(4)
    sage: r1 = A4("(1,2) (3,4)")
    sage: r2 = A4("(1,3) (2,4)")
    sage: r3 = A4("(1,4) (2,3)")
    sage: H = A4.subgroup([r1, r2, r3])
    sage: H.is_normal(A4)
    True


Quotient group
--------------

Extending the previous example, we can create the quotient (factor)
group of `A_4` by `H`. The commands ::

    sage: A4 = AlternatingGroup(4)
    sage: r1 = A4("(1,2) (3,4)")
    sage: r2 = A4("(1,3) (2,4)")
    sage: r3 = A4("(1,4) (2,3)")
    sage: H = A4.subgroup([r1, r2, r3])
    sage: A4.quotient(H)
    Permutation Group with generators [(1,2,3)]

returns a permutation group generated by ``(1,2,3)``.  As expected
this is a group of order 3.  Notice that we do not get back a group of
the actual cosets, but instead we get a group *isomorphic* to the
factor group.


Simple groups
-------------

It is easy to check to see if a group is void of any normal
subgroups. The commands ::

    sage: AlternatingGroup(5).is_simple()
    True
    sage: AlternatingGroup(4).is_simple()
    False

prints ``True`` and then ``False``.


Composition series
------------------

For any group, it is easy to obtain a composition series.  There is an
element of randomness in the algorithm, so you may not always get the
same results.  (But the list of factor groups is unique, according to
the Jordan-Hölder theorem.)  Also, the subgroups generated sometimes
have more generators than necessary, so you might want to
"study" each subgroup carefully by checking properties like its order.

An interesting example is::

    DihedralGroup(105).composition_series()

The output will be a list of 5 subgroups of `D_{105}`, each a normal
subgroup of its predecessor.

Several other series are possible, such as the derived series.  Use
tab-completion to see the possibilities.


Conjugacy
=========

Given a group `G`, we can define a relation `\sim` on `G` by:  for
`a,b \in G`, `a \sim b` if and only if there exists an element
`g \in G` such that `gag^{-1} = b`.

Since this is an equivalence relation, there is an associated
partition of the elements of `G` into equivalence classes.  For this
very important relation, the classes are known as "conjugacy
classes."  A representative of each of these equivalence classes can
be found as follows.  Suppose ``G`` is a permutation group, then
``G.conjugacy_classes_representatives()`` will return a list of
elements of $G$, one per conjugacy class.

Given an element `g \in G`, the "centralizer" of `g` is the set
`C(g) = \{h \in G \mid hgh^{-1} = g\}`, which is a subgroup of `G`.  A
theorem tells us that the size of each conjugacy class is the order
of the group divided by the order of the centralizer of an element of
the class.  With the following code we can determine the size of the
conjugacy classes of the full symmetric group on 5 symbols::

    sage: G = SymmetricGroup(5)
    sage: group_order = G.order()
    sage: reps = G.conjugacy_classes_representatives()
    sage: class_sizes = []
    sage: for g in reps:
    ....:     class_sizes.append(group_order / G.centralizer(g).order())
    ...
    sage: class_sizes
    [1, 10, 15, 20, 20, 30, 24]

This should produce the list ``[1, 10, 15, 20, 20, 30, 24]`` which you
can check sums to 120, the order of the group.  You might be able to
produce this list by counting elements of the group `S_5` with
identical cycle structure (which will require a few simple
combinatorial arguments).


Sylow subgroups
===============

Sylow's Theorems assert the existence of certain subgroups.  For
example, if `p` is a prime, and `p^r` divides the order of a group
`G`, then `G` must have a subgroup of order `p^r`.  Such a subgroup
could be found among the output of the
``conjugacy_classes_subgroups()`` command by checking the orders of
the subgroups produced.  The ``map()`` command is a quick way
to do this.  The symmetric group on 7 symbols, `S_7`, has order
`7! = 5040` and is divisible by `2^4 = 16`.  Let's find one example
of a subgroup of permutations on 4 symbols with order 16::

    sage: G = SymmetricGroup(7)
    sage: subgroups = G.conjugacy_classes_subgroups()
    sage: list(map(order, subgroups))
    [1, 2, 2, 2, 3, 3, 4, 4, 4, 4, 4, 4, 4, 5, 6, 6, 6, 6, 6, 6, 6, 6, 7, 8, 8, 8, 8, 8, 8, 8, 9, 10, 10, 10, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 14, 16, 18, 18, 18, 20, 20, 20, 21, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 36, 36, 36, 36, 40, 42, 48, 48, 48, 60, 60, 72, 72, 72, 72, 120, 120, 120, 120, 144, 168, 240, 360, 720, 2520, 5040]

The ``map(order, subgroups)`` command will apply the ``order()``
function to each of the subgroups in the list ``subgroups``.  The
output is thus a large list of the orders of many subgroups (96 to be
precise).  If you count carefully, you will see that the 49-th
subgroup has order 16.  You can retrieve this group for further study
by referencing it as ``subgroups[48]`` (remember that counting starts
at zero).

If `p^r` is the highest power of `p` to divide the order of `G`, then
a subgroup of order `p^r` is known as a "Sylow `p`-subgroup." Sylow's
Theorems also say any two Sylow `p`-subgroups are conjugate, so the
output of ``conjugacy_classes_subgroups()`` should only contain each
Sylow `p`-subgroup once.  But there is an easier way,
``sylow_subgroup(p)`` will return one. Notice that the argument of the
command is just the prime $p$, not the full power `p^r`.  Failure to
use a prime will generate an informative error message.


Groups of small order as permutation groups
============================================

We list here constructions, as permutation groups, for all of the
groups of order less than 16.

.. CODE-BLOCK:: text

    ---------------------------------------------------------------------------------------------
    Size  Construction                                Notes
    ---------------------------------------------------------------------------------------------
    1     SymmetricGroup(1)                           Trivial
    2     SymmetricGroup(2)                           Also CyclicPermutationGroup(2)
    3     CyclicPermutationGroup(3)                   Prime order
    4     CyclicPermutationGroup(4)                   Cyclic
    4     KleinFourGroup()                            Abelian, non-cyclic
    5     CyclicPermutationGroup(5)                   Prime order
    6     CyclicPermutationGroup(6)                   Cyclic
    6     SymmetricGroup(3)                           Non-abelian, also DihedralGroup(3)
    7     CyclicPermutationGroup(7)                   Prime order
    8     CyclicPermutationGroup(8)                   Cyclic
    8     D1 = CyclicPermutationGroup(4)
          D2 = CyclicPermutationGroup(2)
          G = direct_product_permgroups([D1,D2])      Abelian, non-cyclic
    8     D1 = CyclicPermutationGroup(2)
          D2 = CyclicPermutationGroup(2)
          D3 = CyclicPermutationGroup(2)
          G = direct_product_permgroups([D1,D2,D3])   Abelian, non-cyclic
    8     DihedralGroup(4)                            Non-abelian
    8     QuaternionGroup()                           Quaternions, also DiCyclicGroup(2)
    9     CyclicPermutationGroup(9)                   Cyclic
    9     D1 = CyclicPermutationGroup(3)
          D2 = CyclicPermutationGroup(3)
          G = direct_product_permgroups([D1,D2])      Abelian, non-cyclic
    10    CyclicPermutationGroup(10)                  Cyclic
    10    DihedralGroup(5)                            Non-abelian
    11    CyclicPermutationGroup(11)                  Prime order
    12    CyclicPermutationGroup(12)                  Cyclic
    12    D1 = CyclicPermutationGroup(6)
          D2 = CyclicPermutationGroup(2)
          G = direct_product_permgroups([D1,D2])      Abelian, non-cyclic
    12    DihedralGroup(6)                            Non-abelian
    12    AlternatingGroup(4)                         Non-abelian, symmetries of tetrahedron
    12    DiCyclicGroup(3)                            Non-abelian
                                                      Also semi-direct  product $Z_3 \rtimes Z_4$
    13    CyclicPermutationGroup(13)                  Prime order
    14    CyclicPermutationGroup(14)                  Cyclic
    14    DihedralGroup(7)                            Non-abelian
    15    CyclicPermutationGroup(15)                  Cyclic
    ----------------------------------------------------------------------------------------------


Acknowledgements
================

The construction of Sage is the work of many people, and the group
theory portion is made possible by the extensive work of the creators
of GAP.  However, we will single out three people from the Sage team
to thank for major contributions toward bringing you the group theory
portion of Sage:  David Joyner, William Stein, and Robert Bradshaw.
Thanks!
