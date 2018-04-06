Using Compiled Code Interactively
=================================

This section is about using compiled code in Sage. However, since Sage
is built on top of Python most of this is valid for Python in
general. The exception is that these notes assume you are using Sage's
interface to f2py which makes it more convenient to work with f2py
interactively. You should look at the f2py website for information on
using the command line f2py tool. The ctypes example will work in any 
recent Python install. If you are using Sage, then ctypes and f2py are 
all there already.

Firstly why would we want to write compiled code? Obviously, because
its fast, far faster than interpreted Python code.  Sage has very
powerful facilities that allow one to interactively call compiled code
written in C or Fortran. In fact there 2-4 ways to do this depending
on exactly what you want to accomplish. One way is to use
Cython. Cython is a language that is a hybrid of C and Python based on
Pyrex. It has the ability to call external shared object libraries and
is very useful for writing Python extension modules. Cython/Pyrex is
covered in detail elsewhere in the Sage documentation.

Suppose that you really want to just write Python code, but there is
some particularly time intensive piece of your code that you would
like to either write in C/Fortran or simply call an external shared
library to accomplish. In this case you have three options with
varying strengths and weaknesses.

Note that before you try to use compiled code to speed up your
bottleneck make sure there isn't an easier way. In particular, first
try to vectorize, that is express your algorithm as arithmetic on
vectors or numpy arrays. These arithmetic operations are done directly
in C so will be very fast. If your problem does not lend itself to
being expressed in a vectorized form them read on.

Before we start let us note that this is in no way a complete
introduction to any of the programs we discuss. This is more meant to
orient you to what is possible and what the different options will
feel like.

.. toctree::
   :maxdepth: 2

   f2py
   f2py_examples
   ctypes
   ctypes_examples
   comparison_to_cython
