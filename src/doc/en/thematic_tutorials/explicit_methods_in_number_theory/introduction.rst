Introduction
============

What is Sage?
-------------

Sage (see http://sagemath.org) is a comprehensive mathematical
software system for computations in many areas of pure and applied
mathematics. We program Sage using the mainstream programming
language Python (see http://python.org), or its compiled variant
Cython. It is also very easy to efficiently use code written in
C/C++ from Sage.

The author of this article started the Sage project in 2005.

Sage is free and open source, meaning you can change any part of
Sage and redistribute the result without having to pay any license
fees, and Sage can also leverage the power of commercial
mathematical software such as Magma and Mathematica, if you happen
to have access to those closed source commercial systems.

This document assumes no prior knowledge of either Python or Sage.
Our goal is to help number theorists do computations involving
number fields and modular forms using Sage.

TODO: Overview of Article

As you read this article, please try every example in Sage, and
make sure things works as I claim, and do all of the exercises.
Moreover, you should experiment by typing in similar examples and
checking that the output you get agrees with what you expect.

Using Sage
----------

To use Sage, install it on your computer, and use either the command
line or use the notebook by starting Sage as ``sage -n``.

We show Sage sessions as follows::

    sage: factor(123456)
    2^6 * 3 * 643


This means that if you type ``factor(123456)`` as input to Sage, then
you'll get ``2^6 * 3 * 643`` as output. If you're using the Sage
command line, you type ``factor(123456)`` and press enter; if you're
using the Sage notebook via your web browser, you type
``factor(123456)`` into an input cell and press shift-enter; in the
output cell you'll see ``2^6 * 3 * 643``.

After trying the ``factor`` command in the previous
paragraph (do this now!), you should try factoring some other
numbers.

.. note::

    What happens if you factor a negative number? a rational number?


You can also draw both 2d and 3d pictures using Sage. For example,
the following input plots the number of prime divisors of each
positive integer up to :math:`500`.

::

    sage: line([(n, len(factor(n))) for n in [1..500]])
    Graphics object consisting of 1 graphics primitive


And, this example draws a similar 3d plot::

    sage: import warnings
    sage: warnings.simplefilter('ignore', UserWarning)
    sage: v = [[len(factor(n*m)) for n in [1..15]] for m in [1..15]]
    sage: list_plot3d(v, interpolation_type='clough')
    Graphics3d Object


The Sage-Pari-Magma Ecosystem
-----------------------------

* *The main difference between Sage and Pari is that Sage is vastly
  larger than Pari with a much wider range of functionality, and has
  many more data types and much more structured objects.* Sage in fact
  includes Pari, and a typical Sage install takes nearly a gigabyte of
  disk space, whereas a typical Pari install is much more nimble, using
  only a few megabytes. There are many number-theoretic algorithms that
  are included in Sage, which have never been implemented in Pari, and
  Sage has 2d and 3d graphics which can be helpful for visualizing
  number theoretic ideas, and a graphical user interface. Both Pari and
  Sage are free and open source, which means anybody can read or change
  anything in either program, and the software is free.

* *The biggest difference between Sage and Magma is that Magma is
  closed source, not free, and difficult for users to extend.* This
  means that most of Magma cannot be changed except by the core Magma
  developers, since Magma itself is well over two million lines of
  compiled C code, combined with about a half million lines of
  interpreted Magma code (that anybody can read and modify). In
  designing Sage, we carried over some of the excellent design ideas
  from Magma, such as the parent, element, category hierarchy.

* *Any mathematician who is serious about doing extensive computational
  work in algebraic number theory and arithmetic geometry is strongly
  urged to become familiar with all three systems*, since they all have
  their pros and cons. Pari is sleek and small, Magma has much unique
  functionality for computations in arithmetic geometry, and Sage has a
  wide range of functionality in most areas of mathematics, a large
  developer community, and much unique new code.
