simplegeneric: Simple single-dispatch generic functions for Python
==================================================================

Description
-----------

Simple generic functions (similar to Python's own len(), pickle.dump(),
etc.)

The simplegeneric module lets you define simple single-dispatch generic
functions, akin to Python's built-in generic functions like len() iter()
and so on. However, instead of using specially-named methods, these
generic functions use simple lookup tables, akin to those used by e.g.
pickle.dump() and other generic functions found in the Python standard
library.

As you can see from the above examples, generic functions are actually
quite common in Python already, but there is no standard way to create
simple ones. This library attempts to fill that gap, as generic
functions are an excellent alternative to the Visitor pattern, as well
as being a great substitute for most common uses of adaptation.

This library tries to be the simplest possible implementation of generic
functions, and it therefore eschews the use of multiple or predicate
dispatch, as well as avoiding speedup techniques such as C dispatching
or code generation. But it has absolutely no dependencies, other than
Python 2.4, and the implementation is just a single Python module of
less than 100 lines.
