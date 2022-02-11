r"""
Interface to Groebner Fan

AUTHOR:

- Anders Nedergaard Jensen: Write gfan C++ program, which implements
  algorithms many of which were invented by Jensen, Komei
  Fukuda, and Rekha Thomas.
- William Stein (2006-03-18): wrote gfan interface (first version)
- Marshall Hampton (2008-03-17): modified to use gfan-0.3, subprocess instead of os.popen2

TODO -- much functionality of gfan-0.3 is still not exposed::

   * at most 52 variables:

       - use gfan_substitute to make easier (?)
       MH: I think this is now irrelevant since gfan can accept the original ring variables

   * --symmetry is really useful
            - permutations are 0-based *not* cycle notation; a <---> 0
     output is broken up much more nicely.

   * -- can work in Z/pZ for p <= 32749

   * -- can compute individual GB's for lex and revlex (via buchberger)
"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from subprocess import Popen, PIPE


class Gfan(object):
    """
    Interface to Anders Jensen's Groebner Fan program.
    """
    def __call__(self, I, cmd='', verbose=False, format=True):
        if cmd != '' and cmd.lstrip()[0] != '-':
            cmd = 'gfan_%s'%cmd
        else:
            cmd = 'gfan'

        if len(cmd.split(' ')) > 1:
            cmd = cmd.split(' ')

        if verbose:
            print("gfan command:\n%s" % cmd)
            print("gfan input:\n%s" % I)

        gfan_processes = Popen(cmd, stdin=PIPE, stdout=PIPE, stderr=PIPE,
                               encoding='latin-1')
        ans, err = gfan_processes.communicate(input=I)

        # sometimes, gfan outputs stuff to stderr even though everything is fine
        # we avoid interpreting this as an error
        if (len(err) > 0) and not (err.startswith('_application PolyhedralCone')):
            raise RuntimeError(err)

        return ans

# The instance
gfan = Gfan()
