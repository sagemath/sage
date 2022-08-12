sage_bootstrap: System package database and build scripts for the SageMath distribution
=======================================================================================

Using the system package database
---------------------------------

https://doc.sagemath.org/html/en/developer/portability_testing.html#discovering-the-system-s-package-system

https://doc.sagemath.org/html/en/developer/portability_testing.html#using-sage-s-database-of-distribution-prerequisites

https://doc.sagemath.org/html/en/developer/portability_testing.html#using-sage-s-database-of-equivalent-distribution-packages



This directory contains the build system of Sage, the distribution.

Subdirectories:

 - bin: Various scripts.

    - sage-spkg-info: Format information about a Sage package

    - sage-get-system-packages

    - sage-guess-package-system

    - sage-print-system-package-command

    - sage-package




 - make: Makefiles and related scripts.

 - pkgs: New-style sage packages.

 - pkgs/distros:

 - sage_bootstrap: Python utility library for dealing with
   third-party tarballs and building Sage. See its README for
   more information.

 - test: Test suite for sage_bootstrap.


Developing the sage_bootstrap Python library
--------------------------------------------

This is a utility library for dealing with third-party tarballs and
building Sage. You should never import anything from the Sage
library here, nor should you import anything from sage_bootstrap into
Sage (because this would reconfigure logging). They must be kept
separate.

Everything here must support Python 2.6, 2.7, and 3.3+. Use tox
(https://testrun.org/tox/latest/) to automatically run the tests with
all relevant Python versions. Tests are written as unittest, not as
doctests, because the library is not meant to be used interactively.

Command-line utilities must be able to run as part of a pipe | filter
chain. So you have to be careful about what you send to stdout. You
should use:

* `print()` for anything that is meant to be sent to the output pipe.

* `log.info()` for human-readable messages about the normal program
  flow.
