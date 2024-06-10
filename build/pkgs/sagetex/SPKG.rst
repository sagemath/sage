sagetex: Embed code, results of computations, and plots from Sage into LaTeX documents
======================================================================================

Description
-----------

The SageTeX package allows you to embed code, results of computations,
and plots from Sage into LaTeX documents.

License
-------

The *source code* of the SageTeX package may be redistributed and/or
modified under the terms of the GNU General Public License as published
by the Free Software Foundation, either version 2 of the License, or (at
your option) any later version. To view a copy of this license, see
http://www.gnu.org/licenses/ or send a letter to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301, USA.

The *documentation* of the SageTeX package is licensed under the
Creative Commons Attribution-Share Alike 3.0 License. To view a copy of
this license, visit http://creativecommons.org/licenses/by-sa/3.0/ or
send a letter to Creative Commons, 171 Second Street, Suite 300, San
Francisco, California, 94105, USA.


SPKG Maintainers
----------------

Dan Drake (dr.dan.drake at gmail) and SageMath developers
(sage-devel@googlegroups.com)


Upstream Contact
----------------

Author: Dan Drake.

Web: https://github.com/sagemath/sagetex

Dependencies
------------

To install, nothing more than a standard Sage install. The
``spkg-check`` script will exit without actually testing anything if it
cannot find "latex" in your path.

Notes
-----

To use SageTeX, both Sage and LaTeX need to know about it. SageTeX comes
standard with Sage, so you only need to make sure LaTeX can find what it
needs. Full details are in the Sage installation guide at
http://doc.sagemath.org/html/en/installation/ and
http://doc.sagemath.org/html/en/tutorial/sagetex.html .

The directory ``$SAGE_ROOT/venv/share/doc/sagetex`` contains
documentation and an example file. See
``$SAGE_ROOT/venv/share/texmf/tex/latex/sagetex`` for the source code
and some possibly useful scripts. If you have problems or suggestions
see `the sage-support
group <http://groups.google.com/group/sage-support>`__.

If you want to help develop SageTeX, please clone the github repository
(see the "Upstream Contact" above) and send me patches based on that.
