pygments: Generic syntax highlighter
====================================

Description
-----------

Pygments is a syntax highlighting package written in Python.

It is a generic syntax highlighter suitable for use in code hosting,
forums, wikis or other applications that need to prettify source code.
Highlights are:

-  a wide range of over 300 languages and other text formats is
   supported

-  special attention is paid to details, increasing quality by a fair
   amount

-  support for new languages and formats are added easily
-  a number of output formats, presently HTML, LaTeX, RTF, SVG, all
   image
   formats that PIL supports and ANSI sequences

-  it is usable as a command-line tool and as a library

License
-------

Modified BSD


Upstream Contact
----------------

- Author: Georg Brandl
- Home Page: https://pygments.org

Dependencies
------------

Python


Special Update/Build Instructions
---------------------------------

Patches included:

-  sage_prompt.patch: patch pygments/lexers/agile.py to treat the
   "sage:" prompt like Python's ">>>" prompt. This allows a very
   kludgy patch to be removed from the Sphinx package (see #10118).
