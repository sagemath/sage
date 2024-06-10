This directory contains Pynac, a C++ library for symbolic manipulation based on GiNaC. It can use arbitrary Python objects for numeric types.

# HISTORY #

Until https://github.com/sagemath/sage/issues/32386 (Sage 9.5 development series, 2021), Pynac was maintained in a separate repository, https://github.com/pynac/pynac, and releases were included in Sage as SPKG.

## Earlier history ##

Until Sage 4.0, released on May 29, 2009, most symbolics functionality in Sage relied on Maxima. At that time, the interface to Maxima only worked through pexpect, it has a large overhead. Even though Sage used to have rather complicated wrappers to delay calling Maxima as much as possible, this overhead still made symbolics in Sage hard to use.

Relying on Maxima for even the simplest symbol manipulation also meant that developing new symbolic computation facilities in Sage, such as summation and integration algorithms, was nearly impossible. All real functionality had to be in Maxima, written in Lisp. This made it hard for Sage developers to fix problems in Maxima or enhance its functionality.

In August 2008, Burcin Erocal wrote an experimental wrapper to the C++ library GiNaC. After discussions, it was decided that using GiNaC as the symbolics backend would provide a better framework, and make implementation of new functionality in Sage easier. (See sage-devel archives for details.)

Since Sage already included state of the art libraries for arbitrary precision numeric types (integers, rationals, floating point numbers) GiNaC’s dependence on CLN was redundant. Again in August 2008, William Stein spent 2 weeks working very hard to rewrite GiNaC’s numeric class to use arbitrary Python objects and enhance the initial wrapper. Thus, Pynac was born. Burcin Erocal took over maintenance of the library after that point and worked as time permits to add other features (ordered printing, pickling, arithmetic with infinity, etc.).

In May 2009, with a tremendous effort by Mike Hansen, Robert Bradshaw, William Stein, Carl Witty and Nick Alexander, Sage finally switched to using Pynac as the default backend for symbolic expressions.

https://github.com/pynac/pynac/wiki/Changelog

https://github.com/pynac/pynac/releases


# AUTHORS #

## The Pynac Group ##

see also https://github.com/pynac/pynac/graphs/contributors
and https://github.com/pynac/pynac/wiki/Pynac-Contributors

François Bissey
Robert Bradshaw
Volker Braun
Erik Bray
Burcin Erocal
Jean-Pierre Flori
Benjamin Hackl
Carlos R. Mafra
Paul Masson
R. Andrew Ohana
Clément Pernet
William Stein
Ralf Stephan
Aaditya Thakkar

## The GiNaC Group ##

Christian Bauer <Christian.Bauer@uni-mainz.de>
Chris Dams <Chris.Dams@mi.infn.it>
Alexander Frink <Alexander.Frink@uni-mainz.de>
Vladimir V. Kisil <kisilv@maths.leeds.ac.uk>
Richard Kreckel <Richard.Kreckel@uni-mainz.de>
Alexei Sheplyakov <varg@theor.jinr.ru>
Jens Vollinga <vollinga@thep.physik.uni-mainz.de>

## Authors of optional.hpp ##

In this package the file optional.hpp is distributed under the MIT license.
Author: Martin Moene, see https://github.com/martinmoene/optional-lite

## Contributors of patches ##

Roberto Bagnara, Do Hoang Son, Markus Nullmeier, Pearu Peterson, Benedikt
Pluemper, Ben Sapp, Stefan Weinzierl.


# DOCUMENTATION #

GiNaC (Pynac's ancestor) has great documentation. It's a good idea to read the [tutorial](http://www.ginac.de/tutorial/) to get familiar with the library. The [reference manual](http://www.ginac.de/reference/) for GiNaC can also be useful to find your way around.

https://github.com/pynac/pynac/wiki/Which-part-of-Pynac-gets-actually-used%3F


# CONTRIBUTING #

When making changes to Pynac, keep in mind that GiNaC (Pynac's ancestor) is under continued development, so it should be avoided to make big changes to Pynac that make it harder to apply patches coming from GiNaC.

In particular, respect the existing C++ style and do not apply automatic reformatting of any kind.
 * 8-space expanded tabs (no TAB characters)
 * Code blocks:
    * The braces opening a **function block** should be on a new line.
    * The braces opening a **condition block** should be on the end of the line with the condition. The closing braces should be indented such that they are on the same column as the start of the condition statement.
