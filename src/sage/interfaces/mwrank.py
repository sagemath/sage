r"""
Interface to mwrank
"""

# ****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
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
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import os
import weakref
from .expect import Expect

instances = {}
def Mwrank(options="", server=None, server_tmpdir=None):
    """
    Create and return an mwrank interpreter, with given options.

    INPUT:

    -  ``options`` - string; passed when starting mwrank.
       The format is::

       -h       help            prints this info and quits
       -q       quiet           turns OFF banner display and prompt
       -v n     verbosity       sets verbosity to n (default=1)
       -o       PARI/GP output  turns ON extra PARI/GP short output (default is OFF)
       -p n     precision       sets precision to n decimals (default=15)
       -b n     quartic bound   bound on quartic point search (default=10)
       -x n     n aux           number of aux primes used for sieving (default=6)
       -l       list            turns ON listing of points (default ON unless v=0)
       -s       selmer_only     if set, computes Selmer rank only (default: not set)
       -d       skip_2nd_descent        if set, skips the second descent for curves with 2-torsion (default: not set)
       -S n     sat_bd          upper bound on saturation primes (default=100, -1 for automatic)

    .. WARNING::

       Do not use the option "-q" which turns off the prompt.

    EXAMPLES::

        sage: M = Mwrank('-v 0 -l')
        sage: print(M('0 0 1 -1 0'))
        Curve [0,0,1,-1,0] :    Rank = 1
        Generator 1 is [0:-1:1]; height 0.051...
        Regulator = 0.051...

    """
    global instances
    try:
        X = instances[options]()
        if X:
            return X
    except KeyError:
        pass
    X = Mwrank_class(options, server=server,server_tmpdir=server_tmpdir)
    instances[options] = weakref.ref(X)
    return X

import re
# regex matching '[a1,a2,a3,a4,a6]', no spaces, each ai a possibly signed integer
AINVS_LIST_RE = re.compile(r'\[[+-]?(\d+)(,[+-]?\d+){4}]')
# regex matching ' a1 a2 a3 a4 a6 ', any whitespace, each ai a possibly signed integer
AINVS_PLAIN_RE = re.compile(r'^(\s*)([+-]?(\d+)(\s+)){4}([+-]?(\d+))(\s*)$')


def validate_mwrank_input(s):
    r"""
    Returns a string suitable for mwrank input, or raises an error.

    INPUT:

    - `s` -- one of the following:

        - a list or tuple of 5 integers [a1,a2,a3,a4,a6] or (a1,a2,a3,a4,a6)
        - a string of the form '[a1,a2,a3,a4,a6]' or 'a1 a2 a3 a4 a6' where a1, a2, a3, a4, a6 are integers

    OUTPUT:

    For valid input, a string of the form '[a1,a2,a3,a4,a6]'.  For invalid input a ValueError is raised.

    EXAMPLES:

    A list or tuple of 5 integers::

        sage: from sage.interfaces.mwrank import validate_mwrank_input
        sage: validate_mwrank_input([1,2,3,4,5])
        '[1, 2, 3, 4, 5]'
        sage: validate_mwrank_input((-1,2,-3,4,-55))
        '[-1, 2, -3, 4, -55]'
        sage: validate_mwrank_input([1,2,3,4])
        Traceback (most recent call last):
        ...
        ValueError: [1, 2, 3, 4] is not valid input to mwrank (should have 5 entries)
        sage: validate_mwrank_input([1,2,3,4,i])
        Traceback (most recent call last):
        ...
        ValueError: [1, 2, 3, 4, I] is not valid input to mwrank (entries should be integers)


    A string of the form '[a1,a2,a3,a4,a6]' with any whitespace and integers ai::

        sage: validate_mwrank_input('0 -1 1 -7 6')
        '[0,-1,1,-7,6]'
        sage: validate_mwrank_input("[0,-1,1,0,0]\n")
        '[0,-1,1,0,0]'
        sage: validate_mwrank_input('0\t -1\t 1\t 0\t 0\n')
        '[0,-1,1,0,0]'
        sage: validate_mwrank_input('0 -1 1 -7 ')
        Traceback (most recent call last):
        ...
        ValueError: 0 -1 1 -7  is not valid input to mwrank

    """
    if isinstance(s,(list,tuple)):
        from sage.rings.integer_ring import ZZ
        if len(s)!=5:
            raise ValueError("%s is not valid input to mwrank (should have 5 entries)" % s)
        try:
            ai = [ZZ(a) for a in s]
            return str(ai)
        except (TypeError,ValueError):
            raise ValueError("%s is not valid input to mwrank (entries should be integers)" % s)

    if isinstance(s,str):
        if AINVS_PLAIN_RE.match(s):
            ai = s.split()
            return "["+",".join(ai)+"]"
        ss = s.replace(' ','').replace('\n','').replace('\t','')
        if AINVS_LIST_RE.match(ss):
            return ss
    raise ValueError("%s is not valid input to mwrank" % s)

class Mwrank_class(Expect):
    """
    Interface to the Mwrank interpreter.
    """
    def __init__(self, options="", server=None,server_tmpdir=None):
        """
        INPUT:


        -  ``options`` - string; passed when starting mwrank.
           The format is::

           -h       help            prints this info and quits
           -q       quiet           turns OFF banner display and prompt
           -v n     verbosity       sets verbosity to n (default=1)
           -o       PARI/GP output  turns ON extra PARI/GP short output (default is OFF)
           -p n     precision       sets precision to n decimals (default=15)
           -b n     quartic bound   bound on quartic point search (default=10)
           -x n     n aux           number of aux primes used for sieving (default=6)
           -l       list            turns ON listing of points (default ON unless v=0)
           -s       selmer_only     if set, computes Selmer rank only (default: not set)
           -d       skip_2nd_descent        if set, skips the second descent for curves with 2-torsion (default: not set)
           -S n     sat_bd          upper bound on saturation primes (default=100, -1 for automatic)

        .. WARNING::

            Do not use the option "-q" which turns off the prompt.

        .. NOTE::

           Normally instances of this class would be created by
           calling the global function :meth:`Mwrank`.

        EXAMPLES::

            sage: from sage.interfaces.mwrank import Mwrank_class
            sage: M = Mwrank_class('-v 0 -l')
            sage: M('0 -1 1 0 0')
            'Curve [0,-1,1,0,0] :...Rank = 0...Regulator = 1...'

            sage: from sage.interfaces.mwrank import Mwrank_class
            sage: TestSuite(Mwrank_class).run()
        """
        Expect.__init__(self,
                        name = 'mwrank',
                        prompt = 'Enter curve: ',
                        command = "mwrank %s" % options,
                        server = server,
                        server_tmpdir = server_tmpdir,
                        restart_on_ctrlc = True,
                        verbose_start = False)

    def __getattr__(self, attrname):
        """
        Standard function to return an attribute.

        EXAMPLES::

            sage: mwrank.zzz
            Traceback (most recent call last):
            ...
            AttributeError
        """
        raise AttributeError

    def __reduce__(self):
        """
        EXAMPLES::

            sage: mwrank.__reduce__()
            (<function _reduce_load_mwrank at 0x...>, ())
        """

        return _reduce_load_mwrank, tuple([])

    def __call__(self, cmd):
        """
        Interface to eval method.

        INPUT:

        - ``cmd`` A string, or Sage object which when converted to a
          string gives valid input to ``mwrank``.  The conversion is
          done by :meth:`validate_mwrank_input`.

        EXAMPLES:

        The input can be five integers separated by whitespace::

            sage: mwrank('0 -1 1 0 0')
            'Curve [0,-1,1,0,0] :...Basic pair: I=16, J=-304...'

        Or a list or tuple of exactly five integers::

            sage: s = mwrank([0,-1,1,0,0])
            sage: "Rank = 0" in s and "been determined unconditionally" in s
            True

        TESTS:

        Invalid input raises an ValueError (see :trac:`10108`); this includes
        syntactically valid input which defines a singular curve::

            sage: mwrank(10)
            Traceback (most recent call last):
            ...
            ValueError: Invalid input: 10 is not valid input to mwrank

            sage: mwrank('0 0 0 0 0')
            Traceback (most recent call last):
            ...
            ValueError: Invalid input ([0,0,0,0,0]) to mwrank (singular curve)

            sage: mwrank('0 0 0 -3 2')
            Traceback (most recent call last):
            ...
            ValueError: Invalid input ([0,0,0,-3,2]) to mwrank (singular curve)

        """
        try:
            s = validate_mwrank_input(cmd)
        except ValueError as err:
            raise ValueError("Invalid input: %s" % err)
        try:
            return self.eval(s)
        except ValueError as err:
            raise ValueError(err)
        except RuntimeError:
            raise ValueError(cmd)

    def eval(self, s, **kwds):
        """
        Return mwrank's output for the given input.

        INPUT:

        - ``s`` (str) - a Sage object which when converted to a string
          gives valid input to ``mwrank``.  The conversion is done by
          :meth:`validate_mwrank_input`.  Possible formats are:

          - a string representing exactly five integers separated by
            whitespace, for example '1 2 3 4 5'

          - a string representing exactly five integers separated by
            commas, preceded by '[' and followed by ']' (with
            arbitrary whitespace), for example '[1 2 3 4 5]'

          - a list or tuple of exactly 5 integers.

        .. NOTE::

           If a RuntimeError exception is raised, then the mwrank
           interface is restarted and the command is retried once.

        EXAMPLES::

            sage: mwrank.eval('12 3 4 5 6')
            'Curve [12,3,4,5,6] :...'
            sage: mwrank.eval('[12, 3, 4, 5, 6]')
            'Curve [12,3,4,5,6] :...'
            sage: mwrank.eval([12, 3, 4, 5, 6])
            'Curve [12,3,4,5,6] :...'
            sage: mwrank.eval((12, 3, 4, 5, 6))
            'Curve [12,3,4,5,6] :...'
        """
        if self._expect is not None and not self._expect.isalive():
            # if mwrank is interrupted twice in rapid succession,
            # then it doesn't restart correctly, and we're left with:
            #   "RuntimeError: [Errno 9] Bad file descriptor"
            # Doing _start again fixes that always. See trac #5157.
            self._start()
        try:
            ss = validate_mwrank_input(s)
            return Expect.eval(self, ss, **kwds)
        except ValueError as err:
            raise ValueError('Invalid input: %s' % err)
        except RuntimeError:
            raise ValueError('Invalid input (%s) to mwrank (singular curve)' % s)

    def console(self):
        """
        Start the mwrank console.

        EXAMPLES::

            sage: mwrank.console() # not tested: expects console input
            Program mwrank: ...

        """
        mwrank_console()


# An instance
mwrank = Mwrank()

def _reduce_load_mwrank():
    """
    Return the standard mwrank instance

    EXAMPLES::

        sage: from sage.interfaces.mwrank import _reduce_load_mwrank
        sage: _reduce_load_mwrank()
        Mwrank
    """
    return mwrank


def mwrank_console():
    """
    Start the mwrank console.

    EXAMPLES::

        sage: mwrank_console() # not tested: expects console input
        Program mwrank: ...
    """
    from sage.repl.rich_output.display_manager import get_display_manager
    if not get_display_manager().is_in_terminal():
        raise RuntimeError('Can use the console only in the terminal. Try %%mwrank magics instead.')
    os.system('mwrank')

