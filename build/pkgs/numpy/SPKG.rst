numpy: Package for scientific computing with Python
===================================================

Description
-----------

This package adds numerical linear algebra and other numerical computing
capabilities to python.


Upstream Contact
----------------

-  https://numpy.org/
-  Travis Oliphant
-  Fernando Perez
-  Brian Granger

Special Update/Build Instructions
---------------------------------

-  Scipy uses numpy's distutils to control its compilation of fortran
   code.

   Whenever numpy is updated it is necessary to make sure that scipy
   still builds ok.
