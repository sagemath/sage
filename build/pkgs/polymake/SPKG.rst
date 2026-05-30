polymake: Computations with polyhedra, fans, simplicial complexes, matroids, graphs, tropical hypersurfaces
===========================================================================================================

Description
-----------

polymake is open source software for research in polyhedral geometry. It
deals with polytopes, polyhedra and fans as well as simplicial
complexes, matroids, graphs, tropical hypersurfaces, and other objects.
Supported platforms include various flavors of Linux, Free BSD and Mac
OS.

License
-------

-  GPL v3


Upstream Contact
----------------

-  https://polymake.org/

Dependencies
------------

Polymake needs a working installation of Perl, including its shared
library and some modules (``XML::Writer XML::LibXML XML::LibXSLT
Term::ReadLine::Gnu JSON SVG``). The Polymake interface in Sage
additionally needs ``File::Slurp``. For full functionality including
polymake's polyDB, also the Perl module ``MongoDB`` is needed.

These are not provided by a Sage package. The dummy package
``perl_cpan_polymake_prereq`` will signal an error at build time if the
required prerequisites are not met.

The ``configure`` script will inform you about the equivalent system
packages that you should install. Otherwise, you can use CPAN (see
below).

Sage might install the ``Term::ReadLine::Gnu`` module, however, when you
install polymake, if it is not provided by the system, or if Sage
installs its own ``readline`` library.


A distribution-independent way to install Perl modules (into a user's
home directory or ``/usr/local``) is using CPAN. This is also the way to
install the modules on macOS. For this, if you don't have root access,
you will need the ``local::lib`` Perl module installed::

   cpan -i XML::Writer XML::LibXML XML::LibXSLT File::Slurp Term::ReadLine::Gnu JSON SVG MongoDB

Before installing the ``polymake`` package, refer to the SPKG pages for the following packages to ensure a more featureful Polymake installation:

- [4ti2](https://doc.sagemath.org/html/en/reference/spkg/4ti2.html)
- [latte_int](https://doc.sagemath.org/html/en/reference/spkg/latte_int.html)
- [topcom](https://doc.sagemath.org/html/en/reference/spkg/topcom.html)
- [qhull](https://doc.sagemath.org/html/en/reference/spkg/qhull.html)

For additional software that may enhance your Polymake installation (but for which no Sage package is available), you can manually install the following:

- ``azove``
- ``porta``
- ``vinci``
- ``SplitsTree4``

Information on missing Polymake prerequisites after installing polymake::

   $ sage -sh
   (sage-sh) $ polymake
   polytope> show_unconfigured;

In order to use Polymake from Sage, please refer to the [Jupymake SPKG page](https://doc.sagemath.org/html/en/reference/spkg/jupymake.html) for installation instructions.



Debugging polymake install problems
-----------------------------------

::

  # apt-get install libdevel-trace-perl
  $ cd src
  $ perl -d:Trace support/configure.pl
