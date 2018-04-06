Comparison to Cython/Pyrex
==========================

It is certainly possible to write a solver in Cython or Pyrex. From
the
http://www.scipy.org/PerformancePython?highlight=\%28performance\%29
website you can find an example. One potential downside to Cython over
the previous solutions is it requires the user to understand how NumPy
arrays or Sage matrices are implemented so as to be able to access
their internal data. In contrast the scipy, and ctypes examples only 
require the user to know C or Fortran and from their perspective
the NumPy data magically gets passed to C or Fortran with no further
thought from them. In order for pyrex to be competitive as a way to
interactively write compiled code, the task of accessing the internal
structure of NumPy arrays or Sage matrices needs to be hidden.
