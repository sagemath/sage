ncurses: Classic terminal output library
========================================

Description
-----------

Ncurses (new curses, pronounced "enn-curses") started as a freely
distributable "clone" of System V Release 4.0 (SVr4) curses. It has
outgrown the "clone" description, and now contains many features which
are not in SVr4 curses. Curses is a pun on the term "cursor
optimization". It is a library of functions that manage an application's
display on character-cell terminals (e.g., VT100).

The name "ncurses" was first used as the name of the curses library in
Pavel Curtis's pcurses, dated 1982. It was apparently developed on a BSD
4.4 system, at Cornell. Parts of pcurses are readily identifiable in
ncurses, including the basics for the terminfo compiler (named compile
in that package):

-  the Caps, used to define the terminfo capabilities
-  awk scripts MKcaptab.awk, MKnames.awk
-  the library modules used for the terminfo compiler.

Besides ncurses, parts of pcurses still survive in 2010, in recognizable
form in Solaris.

Website: http://invisible-island.net/ncurses

License
-------

-  MIT-style


Upstream Contact
----------------

-  bug-ncurses@gnu.org

Special Update/Build Instructions
---------------------------------

None
