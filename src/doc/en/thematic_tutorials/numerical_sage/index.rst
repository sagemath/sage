.. Numerical Sage documentation master file, created by sphinx-quickstart on Sat Dec  6 11:08:04 2008.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. _numerical_computing:

Numerical Computing with Sage
=============================

.. WARNING::

    Beware that this document may be obsolete.

This document is designed to introduce the reader to the tools in Sage
that are useful for doing numerical computation. By numerical
computation we essentially mean machine precision floating point
computations.  In particular, things such as optimization, numerical
linear algebra, solving ODE's or PDE's numerically, etc.

In the first part of this document the reader is only assumed to be
familiar with Python/Sage. In the second section on using compiled
code, the computational prerequisites increase and I assume the reader
is comfortable with writing programs in C or Fortran. The third
section is on mpi and parallel programming and only requires knowledge
of Python, though familiarity with mpi would be helpful.

In the current version of this document the reader is assumed to be
familiar with the techniques of numerical analysis. The goal of this
document is not to teach you numerical analysis, but to explain how to
express your ideas in Sage and Python. Also this document is not meant
to be comprehensive. Instead the goal is to be a road map and orient
the reader to the packages relevant to numerical computation and to
where they can find more information.

.. toctree::
   :maxdepth: 2

   numerical_tools
   using_compiled_code_iteractively
   parallel_computation
