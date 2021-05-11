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
library and some modules (XML::Writer XML::LibXML XML::LibXSLT
Term::ReadLine::Gnu JSON SVG). The Polymake interface in Sage
additionally needs File::Slurp. For full functionality including
polymake's polyDB, also the Perl module MongoDB is required.

These are not provided by a Sage package. The script package
perl_cpan_polymake_prereq will signal an error at build time if these
prerequisites are not met.

The configure script will inform you about the equivalent system
packages that you should install. Otherwise, you can use CPAN (see
below).

Sage might install the Term::ReadLine::Gnu module, however, when you
install polymake, if it is not provided by the system, or if Sage
installs its own readline library.


A distribution-independent way to install Perl modules (into a user's
home directory or /usr/local) is using CPAN. This is also the way to
install the modules on macOS. For this, if you don't have root access,
you will need the local::lib Perl module installed::

   cpan -i XML::Writer XML::LibXML XML::LibXSLT File::Slurp Term::ReadLine::Gnu JSON SVG MongoDB

Several Sage packages should be installed before installing the polymake
package to give a more featureful Polymake installation:

   sage -i 4ti2 latte_int topcom qhull

Software that would need to be installed manually (no Sage package
available) for a more featureful Polymake installation: azove, porta,
vinci, SplitsTree4.

Information on missing Polymake prerequisites after installing polymake::

   $ sage -sh
   (sage-sh) $ polymake
   polytope> show_unconfigured;


Debugging polymake install problems
-----------------------------------

::

  # apt-get install libdevel-trace-perl
  $ cd src
  $ perl -d:Trace support/configure.pl
