r"""
Interface to Macaulay2

.. NOTE::

    You must have ``Macaulay2`` installed on your computer
    for this interface to work. Macaulay2 is not included with Sage,
    but you can obtain it from <http://www.math.uiuc.edu/Macaulay2/>.
    Note additional optional Sage packages are required.

Sage provides an interface to the Macaulay2 computational algebra
system. This system provides extensive functionality for commutative
algebra. You do not have to install any optional packages.

The Macaulay2 interface offers three pieces of functionality:

- ``Macaulay2_console()`` -- A function that dumps you
  into an interactive command-line Macaulay2 session.

- ``Macaulay2(expr)`` -- Evaluation of arbitrary Macaulay2
  expressions, with the result returned as a string.

- ``Macaulay2.new(expr)`` -- Creation of a Sage object that wraps a
  Macaulay2 object.  This provides a Pythonic interface to Macaulay2.  For
  example, if ``f=Macaulay2.new(10)``, then ``f.gcd(25)`` returns the
  GCD of `10` and `25` computed using Macaulay2.

EXAMPLES::

    sage: print(macaulay2('3/5 + 7/11')) # optional - macaulay2
    68
    --
    55
    sage: f = macaulay2('f = i -> i^3') # optional - macaulay2
    sage: f                             # optional - macaulay2
    f
    sage: f(5)                          # optional - macaulay2
    125

    sage: R = macaulay2('ZZ/5[x,y,z]')  # optional - macaulay2
    sage: print(R)                       # optional - macaulay2
    ZZ
    --[x..z, Degrees => {3:1}, Heft => {1}, MonomialOrder => {MonomialSize => 32}, DegreeRank => 1]
     5                                                       {GRevLex => {3:1}  }
                                                             {Position => Up    }
    sage: x = macaulay2('x')            # optional - macaulay2
    sage: y = macaulay2('y')            # optional - macaulay2
    sage: print((x+y)^5)                # optional - macaulay2
     5    5
    x  + y
    sage: parent((x+y)^5)               # optional - macaulay2
    Macaulay2

    sage: R = macaulay2('QQ[x,y,z,w]')  # optional - macaulay2
    sage: f = macaulay2('x^4 + 2*x*y^3 + x*y^2*w + x*y*z*w + x*y*w^2 + 2*x*z*w^2 + y^4 + y^3*w + 2*y^2*z*w + z^4 + w^4') # optional - macaulay2
    sage: print(f)                      # optional - macaulay2
     4       3    4    4      2     3                2           2         2    4
    x  + 2x*y  + y  + z  + x*y w + y w + x*y*z*w + 2y z*w + x*y*w  + 2x*z*w  + w
    sage: g = f * macaulay2('x+y^5')    # optional - macaulay2
    sage: print(g.factor())             # optional - macaulay2
      4       3    4    4      2     3                2           2         2    4   5
    (x  + 2x*y  + y  + z  + x*y w + y w + x*y*z*w + 2y z*w + x*y*w  + 2x*z*w  + w )(y  + x)


AUTHORS:

- Kiran Kedlaya and David Roe (2006-02-05, during Sage coding sprint)
- William Stein (2006-02-09): inclusion in Sage; prompt uses regexp,
  calling of Macaulay2 functions via __call__.
- William Stein (2006-02-09): fixed bug in reading from file and
  improved output cleaning.
- Kiran Kedlaya (2006-02-12): added ring and ideal constructors,
  list delimiters, is_Macaulay2Element, sage_polystring,
  __floordiv__, __mod__, __iter__, __len__; stripped extra
  leading space and trailing newline from output.

TODO:

- get rid of all numbers in output, e.g., in ideal function below.
"""

#*****************************************************************************
#       Copyright (C) 2006 Kiran S. Kedlaya <kedlaya@mit.edu>
#                          David Roe <roed@mit.edu>
#                          William Stein <wstein@gmail.com>
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
from __future__ import print_function

import os

from sage.interfaces.expect import (Expect, ExpectElement, ExpectFunction,
                                    AsciiArtString)

from sage.misc.multireplace import multiple_replace
from sage.interfaces.tab_completion import ExtraTabCompletion

import re

def remove_output_labels(s):
    r"""
    Remove output labels of Macaulay2 from a string.

    - s: output of Macaulay2

    - s: string

    Returns: the input string with `n` symbols removed from the beginning of
    each line, where `n` is the minimal number of spaces or symbols of
    Macaulay2 output labels (looking like 'o39 = ') present on every non-empty
    line.

    Return type: string

    .. note::

        If ``s`` consists of several outputs and their lables have
        different width, it is possible that some strings will have leading
        spaces (or maybe even pieces of output labels). However, this
        function will try not cut any messages.

    EXAMPLES::

        sage: from sage.interfaces.macaulay2 import remove_output_labels
        sage: output = 'o1 = QQ [x, y]\n\no1 : PolynomialRing\n'
        sage: remove_output_labels(output)
        'QQ [x, y]\n\nPolynomialRing\n'
    """
    label = re.compile("^o[0-9]+ (=|:) |^ *")
    lines = s.split("\n")
    matches = [label.match(l) for l in lines if l != ""]
    if len(matches) == 0:
        return s
    else:
        n = min(m.end() - m.start() for m in matches)
        return "\n".join(l[n:] for l in lines)


PROMPT = "_EGAS_ : "


class Macaulay2(ExtraTabCompletion, Expect):
    """
    Interface to the Macaulay2 interpreter.
    """
    def __init__(self, maxread=None, script_subdirectory=None,
                 logfile=None, server=None,server_tmpdir=None):
        """
        Initialize a Macaulay2 interface instance.

        We replace the standard input prompt with a strange one, so that
        we do not catch input prompts inside the documentation.

        We replace the standard input continuation prompt, which is
        just a bunch of spaces and cannot be automatically detected in a
        reliable way. This is necessary for allowing commands that occupy
        several strings.

        We also change the starting line number to make all the output
        labels to be of the same length. This allows correct stripping of
        the output of several commands.

        TESTS::

            sage: macaulay2 == loads(dumps(macaulay2))
            True
        """
        init_str = (
            # Prompt changing commands
            """ZZ#{Standard,Core#"private dictionary"#"InputPrompt"} = lineno -> "%s";""" % PROMPT +
            """ZZ#{Standard,Core#"private dictionary"#"InputContinuationPrompt"} = lineno -> "%s";""" % PROMPT +
            # Also prevent line wrapping in Macaulay2
            "printWidth = 0;" +
            # And make all output labels to be of the same width
            "lineNumber = 10^9;")
        Expect.__init__(self,
                        name = 'macaulay2',
                        prompt = PROMPT,
                        command = "M2 --no-debug --no-readline --silent -e '%s'" % init_str,
                        server = server,
                        server_tmpdir = server_tmpdir,
                        script_subdirectory = script_subdirectory,
                        verbose_start = False,
                        logfile = logfile,
                        eval_using_file_cutoff=500)

    # Macaulay2 provides no "clear" function. However, Macaulay2 does provide
    # garbage collection; since expect automatically reuses variable names,
    # garbage collection in Sage properly sets up garbage collection in
    # Macaulay2.

    def __reduce__(self):
        """
        Used in serializing an Macaulay2 interface.

        EXAMPLES::

            sage: rlm2, t = macaulay2.__reduce__()
            sage: rlm2(*t)
            Macaulay2
        """
        return reduce_load_macaulay2, tuple([])

    def _read_in_file_command(self, filename):
        """
        Load and *execute* the content of ``filename`` in Macaulay2.

        INPUT:

        - filename: the name of the file to be loaded and executed
          (type: string)

        OUTPUT:

        Returns Macaulay2 command loading and executing commands in
        ``filename``, that is, ``'load "filename"'``.
        Return type: string

        TESTS::

            sage: filename = tmp_filename()
            sage: f = open(filename, "w")
            sage: f.write("sage_test = 7;")
            sage: f.close()
            sage: command = macaulay2._read_in_file_command(filename)
            sage: macaulay2.eval(command)  # optional - macaulay2
            sage: macaulay2.eval("sage_test")  # optional - macaulay2
            7
            sage: import os
            sage: os.unlink(filename)
            sage: macaulay2._read_in_file_command("test")
            'load "test"'
            sage: macaulay2(10^10000) == 10^10000  # optional - macaulay2
            True
        """
        return 'load "%s"' % filename

    def __getattr__(self, attrname):
        """
        EXAMPLES::

            sage: gb = macaulay2.gb  # optional - macaulay2
            sage: type(gb)           # optional - macaulay2
            <class 'sage.interfaces.macaulay2.Macaulay2Function'>
            sage: gb._name           # optional - macaulay2
            'gb'
        """
        if attrname[:1] == "_":
            raise AttributeError
        return Macaulay2Function(self, attrname)

    def eval(self, code, strip=True, **kwds):
        """
        Send the code x to the Macaulay2 interpreter and return the output
        as a string suitable for input back into Macaulay2, if possible.

        INPUT:

        - code -- str
        - strip -- ignored

        EXAMPLES::

            sage: macaulay2.eval("2+2") # optional - macaulay2
            4
        """
        code = code.strip()
        # TODO: in some cases change toExternalString to toString??
        ans = Expect.eval(self, code, strip=strip, **kwds).strip('\n')
        if strip:
            ans = remove_output_labels(ans)
        return AsciiArtString(ans)

    def restart(self):
        r"""
        Restart Macaulay2 interpreter.

        TEST::

            sage: macaulay2.restart()  # optional - macaulay2
        """
        # If we allow restart to be called as a function, there will be
        # parasitic output
        self.eval("restart")

    def get(self, var):
        """
        Get the value of the variable var.

        EXAMPLES::

            sage: macaulay2.set("a", "2") # optional - macaulay2
            sage: macaulay2.get("a")      # optional - macaulay2
            2
        """
        return self.eval("describe %s"%var, strip=True)

    def set(self, var, value):
        """
        Set the variable var to the given value.

        EXAMPLES::

            sage: macaulay2.set("a", "2")  # optional - macaulay2
            sage: macaulay2.get("a")       # optional - macaulay2
            2
        """
        cmd = '%s=%s;'%(var,value)
        ans = Expect.eval(self, cmd)
        if ans.find("stdio:") != -1:
            raise RuntimeError("Error evaluating Macaulay2 code.\nIN:%s\nOUT:%s"%(cmd, ans))

    def _object_class(self):
        """
        Returns the class of Macaulay2 elements.

        EXAMPLES::

            sage: macaulay2._object_class()
            <class 'sage.interfaces.macaulay2.Macaulay2Element'>

        """
        return Macaulay2Element

    def console(self):
        """
        Spawn a new M2 command-line session.

        EXAMPLES::

            sage: macaulay2.console()                    # not tested
            Macaulay 2, version 1.1
            with packages: Classic, Core, Elimination, IntegralClosure, LLLBases, Parsing, PrimaryDecomposition, SchurRings, TangentCone
            ...

        """
        macaulay2_console()

    def _left_list_delim(self):
        """
        Returns the Macaulay2 left delimiter for lists.

        EXAMPLES::

            sage: macaulay2._left_list_delim()
            '{'
        """
        return '{'

    def _right_list_delim(self):
        """
        Returns the Macaulay2 right delimiter for lists.

        EXAMPLES::

            sage: macaulay2._right_list_delim()
            '}'
        """
        return '}'

    def _true_symbol(self):
        """
        Returns the Macaulay2 symbol for True.

        EXAMPLES::

            sage: macaulay2._true_symbol()
            'true'
        """
        return 'true'

    def _false_symbol(self):
        """
        Returns the Macaulay2 symbol for False.

        EXAMPLES::

            sage: macaulay2._false_symbol()
            'false'
        """
        return 'false'

    def _equality_symbol(self):
        """
        Returns the Macaulay2 symbol for equality.

        EXAMPLES::

            sage: macaulay2._false_symbol()
            'false'
        """
        return '=='

    def cputime(self, t=None):
        """
        EXAMPLES::

            sage: R = macaulay2("QQ[x,y]")  # optional - macaulay2
            sage: x,y = R.gens()            # optional - macaulay2
            sage: a = (x+y+1)^20            # optional - macaulay2
            sage: macaulay2.cputime()       # optional - macaulay2; random
            0.48393700000000001
        """
        _t = float(self.cpuTime().to_sage())
        if t:
            return _t - t
        else:
            return _t

    def version(self):
        """
        Returns the version of Macaulay2.

        EXAMPLES::

            sage: macaulay2.version() # optional - macaulay2
            (1, 3, 1)
        """
        s = self.eval("version")
        r = re.compile("VERSION => (.*?)\n")
        s = r.search(s).groups()[0]
        return tuple(int(i) for i in s.split("."))

### Constructors

    def ideal(self, *gens):
        """
        Return the ideal generated by gens.

        INPUT:

        - gens -- list or tuple of Macaulay2 objects (or objects that can be
          made into Macaulay2 objects via evaluation)

        OUTPUT:

        the Macaulay2 ideal generated by the given list of gens

        EXAMPLES::

            sage: R2 = macaulay2.ring('QQ', '[x, y]'); R2            # optional - macaulay2
            QQ[x..y, Degrees => {2:1}, Heft => {1}, MonomialOrder => {MonomialSize => 16}, DegreeRank => 1]
                                                                     {Lex => 2          }
                                                                     {Position => Up    }
            sage: I = macaulay2.ideal( ('y^2 - x^3', 'x - y') ); I   # optional - macaulay2
                      3    2
            ideal (- x  + y , x - y)
            sage: J = I^3; J.gb().gens().transpose()                 # optional - macaulay2
            {-9} | y9-3y8+3y7-y6             |
            {-7} | xy6-2xy5+xy4-y7+2y6-y5    |
            {-5} | x2y3-x2y2-2xy4+2xy3+y5-y4 |
            {-3} | x3-3x2y+3xy2-y3           |

        """
        if len(gens) == 1 and isinstance(gens[0], (list, tuple)):
            gens = gens[0]
        gens2 = []
        for g in gens:
            if not isinstance(g, Macaulay2Element):
                gens2.append(self(g))
            else:
                gens2.append(g)
        return self('ideal {%s}'%(",".join([g.name() for g in gens2])))

    def ring(self, base_ring='ZZ', vars='[x]', order='Lex'):
        r"""
        Create a Macaulay2 ring.

        INPUT:

        - base_ring -- base ring (see examples below)
        - vars -- a tuple or string that defines the variable names
        - order -- string -- the monomial order (default: 'Lex')

        OUTPUT:

        - a Macaulay2 ring (with base ring ZZ)

        EXAMPLES:

        This is a ring in variables named a through d over the finite field
        of order 7, with graded reverse lex ordering::

            sage: R1 = macaulay2.ring('ZZ/7', '[a..d]', 'GRevLex');  R1  # optional - macaulay2
            ZZ
            --[a..d, Degrees => {4:1}, Heft => {1}, MonomialOrder => {MonomialSize => 16}, DegreeRank => 1]
             7                                                       {GRevLex => {4:1}  }
                                                                     {Position => Up    }
            sage: R1.char()                                             # optional - macaulay2
            7

        This is a polynomial ring over the rational numbers::

            sage: R2 = macaulay2.ring('QQ', '[x, y]'); R2               # optional - macaulay2
            QQ[x..y, Degrees => {2:1}, Heft => {1}, MonomialOrder => {MonomialSize => 16}, DegreeRank => 1]
                                                                     {Lex => 2          }
                                                                     {Position => Up    }
        """
        varstr = str(vars)[1:-1]
        if ".." in varstr:
            varstr = "symbol " + varstr[0] + ".." + "symbol " + varstr[-1]
        else:
            varstr = ", ".join(["symbol " + v for v in varstr.split(", ")])
        return self.new('%s[%s, MonomialSize=>16, MonomialOrder=>%s]'%(base_ring, varstr, order))

    def help(self, s):
        """
        EXAMPLES::

            sage: macaulay2.help("load")  # optional - macaulay2
            load -- read Macaulay2 commands
            *******************************
            ...
              * "input" -- read Macaulay2 commands and echo
              * "notify" -- whether to notify the user when a file is loaded
        """
        r = self.eval("help %s" % s)
        end = r.rfind("\n\nDIV")
        if end != -1:
            r = r[:end]
        return AsciiArtString(r)

    def _tab_completion(self):
        """
        Return a list of tab completions for Macaulay2.

        Returns dynamically built sorted list of commands obtained using
        Macaulay2 "apropos" command.
        Return type: list of strings

        TESTS::

            sage: names = macaulay2._tab_completion() # optional - macaulay2
            sage: 'ring' in names                 # optional - macaulay2
            True
            sage: macaulay2.eval("abcabc = 4")    # optional - macaulay2
            4
            sage: names = macaulay2._tab_completion() # optional - macaulay2
            sage: "abcabc" in names               # optional - macaulay2
            True
        """
        # Get all the names from Macaulay2 except numbered outputs like
        # o1, o2, etc. and automatic Sage variable names sage0, sage1, etc.
        # It is faster to get it back as a string.
        r = macaulay2.eval(r"""
            toString select(
                apply(apropos "^[[:alnum:]]+$", toString),
                s -> not match("^(o|sage)[0-9]+$", s))
            """)
        # Now split this string into separate names
        r = sorted(r[1:-1].split(", "))
        # Macaulay2 sorts things like A, a, B, b, ...
        return r

    def use(self, R):
        """
        Use the Macaulay2 ring R.

        EXAMPLES::

            sage: R = macaulay2("QQ[x,y]")                  # optional - macaulay2
            sage: P = macaulay2("ZZ/7[symbol x, symbol y]") # optional - macaulay2
            sage: macaulay2("x").cls()                      # optional - macaulay2
            ZZ
            --[x..y, Degrees => {2:1}, Heft => {1}, MonomialOrder => {MonomialSize => 32}, DegreeRank => 1]
             7                                                       {GRevLex => {2:1}  }
                                                                     {Position => Up    }
            sage: macaulay2.use(R)                          # optional - macaulay2
            sage: macaulay2("x").cls()                      # optional - macaulay2
            QQ[x..y, Degrees => {2:1}, Heft => {1}, MonomialOrder => {MonomialSize => 32}, DegreeRank => 1]
                                                                     {GRevLex => {2:1}  }
                                                                     {Position => Up    }
        """
        R = self(R)
        self.eval("use %s"%R.name())

    def new_from(self, type, value):
        """
        Returns a new Macaulay2Element of type type constructed from
        value.

        EXAMPLES::

            sage: l = macaulay2.new_from("MutableList", [1,2,3]) # optional - macaulay2
            sage: l                                              # optional - macaulay2
            MutableList{...3...}
            sage: list(l)                                        # optional - macaulay2
            [1, 2, 3]

        """
        type = self(type)
        value = self(value)
        return self.new("new %s from %s"%(type.name(), value.name()))


class Macaulay2Element(ExtraTabCompletion, ExpectElement):
    def _latex_(self):
        """
        EXAMPLES::

            sage: m = macaulay2('matrix {{1,2},{3,4}}') # optional - macaulay2
            sage: m                                     # optional - macaulay2
            | 1 2 |
            | 3 4 |
            sage: latex(m) # optional - macaulay2
            \begin{pmatrix}1& {2}\\ {3}& {4}\\ \end{pmatrix}
        """
        s = self.tex().external_string().strip('"').strip('$').replace('\\\\','\\')
        s = s.replace(r"\bgroup","").replace(r"\egroup","")
        return s

    def __iter__(self):
        """
        EXAMPLES::

            sage: l = macaulay2([1,2,3]) # optional - macaulay2
            sage: list(iter(l))          # optional - macaulay2
            [1, 2, 3]
        """
        for i in range(len(self)):  # zero-indexed!
            yield self[i]

    def __str__(self):
        """
        EXAMPLES::

            sage: R = macaulay2("QQ[x,y,z]/(x^3-y^3-z^3)") # optional - macaulay2
            sage: x = macaulay2('x')                       # optional - macaulay2
            sage: y = macaulay2('y')                       # optional - macaulay2
            sage: print(x+y)                               # optional - macaulay2
            x + y
            sage: print(macaulay2("QQ[x,y,z]"))            # optional - macaulay2
            QQ[x..z, Degrees => {3:1}, Heft => {1}, MonomialOrder => {MonomialSize => 32}, DegreeRank => 1]
                                                                     {GRevLex => {3:1}  }
                                                                     {Position => Up    }
            sage: print(macaulay2("QQ[x,y,z]/(x+y+z)"))    # optional - macaulay2
            QQ[x, y, z]
            -----------
             x + y + z
        """
        P = self._check_valid()
        return P.get(self._name)

    repr = __str__

    def external_string(self):
        """
        EXAMPLES::

           sage: R = macaulay2("QQ[symbol x, symbol y]")  # optional - macaulay2
           sage: R.external_string()                      # optional - macaulay2
           'QQ[x..y, Degrees => {2:1}, Heft => {1}, MonomialOrder => VerticalList{MonomialSize => 32, GRevLex => {2:1}, Position => Up}, DegreeRank => 1]'
        """
        P = self._check_valid()
        code = 'toExternalString(%s)'%self.name()
        X = P.eval(code, strip=True)

        if 'stdio:' in X:
            if 'to external string' in X:
                return P.eval('%s'%self.name())
            raise RuntimeError("Error evaluating Macaulay2 code.\nIN:%s\nOUT:%s"%(code, X))

        s = multiple_replace({'\r':'', '\n':' '}, X)
        return s

    def __len__(self):
        """
        EXAMPLES::

            sage: l = macaulay2([1,2,3])  # optional - macaulay2
            sage: len(l)                  # optional - macaulay2
            3
            sage: type(_)                 # optional - macaulay2
            <type 'int'>
        """
        self._check_valid()
        return int(self.parent()("#%s"%self.name()))

    def __getitem__(self, n):
        """
        EXAMPLES::

            sage: l = macaulay2([1,2,3])  # optional - macaulay2
            sage: l[0]                    # optional - macaulay2
            1
        """
        self._check_valid()
        n = self.parent()(n)
        return self.parent().new('%s # %s'%(self.name(), n.name()))

    def __setitem__(self, index, value):
        """
        EXAMPLES::

            sage: l = macaulay2.new_from("MutableList", [1,2,3]) # optional - macaulay2
            sage: l[0] = 4               # optional - macaulay2
            sage: list(l)                # optional - macaulay2
            [4, 2, 3]

        """
        P = self.parent()
        index = P(index)
        value = P(value)
        res = P.eval("%s # %s = %s"%(self.name(), index.name(), value.name()))
        if "assignment attempted to element of immutable list" in res:
            raise TypeError("item assignment not supported")

    def __call__(self, x):
        """
        EXAMPLES::

            sage: R = macaulay2("QQ[x, y]")     # optional - macaulay2
            sage: x,y = R.gens()                # optional - macaulay2
            sage: I = macaulay2.ideal(x*y, x+y) # optional - macaulay2
            sage: gb = macaulay2.gb             # optional - macaulay2
            sage: gb(I)                         # optional - macaulay2
            GroebnerBasis[status: done; S-pairs encountered up to degree 1]
        """
        self._check_valid()
        P = self.parent()
        r = P(x)
        return P('%s %s'%(self.name(), r.name()))

    def __floordiv__(self, x):
        """
        Quotient of division of self by other.  This is denoted //.

        EXAMPLE::

            sage: R.<x,y> = GF(7)[]

        Now make the M2 version of R, so we can coerce elements of R to M2::

            sage: macaulay2(R)                           # optional - macaulay2
            ZZ
            --[x..y, Degrees => {2:1}, Heft => {1}, MonomialOrder => {MonomialSize => 16}, DegreeRank => 1]
             7                                                       {GRevLex => {2:1}  }
                                                                     {Position => Up    }
            sage: f = (x^3 + 2*y^2*x)^7; f
            x^21 + 2*x^7*y^14
            sage: h = macaulay2(f); h                    # optional - macaulay2
             21     7 14
            x   + 2x y
            sage: f1 = (x^2 + 2*y*x)                     # optional - macaulay2
            sage: h1 = macaulay2(f1)                     # optional - macaulay2
            sage: f2 = (x^3 + 2*y*x)                     # optional - macaulay2
            sage: h2 = macaulay2(f2)                     # optional - macaulay2
            sage: u = h // [h1,h2]                       # optional - macaulay2
            sage: h == u[0]*h1 + u[1]*h2 + (h % [h1,h2]) # optional - macaulay2
            True
        """
        if isinstance(x, (list, tuple)):
            y = self.parent(x)
            z = self.parent().new('%s // matrix{%s}'%(self.name(), y.name()))
            return list(z.entries().flatten())
        else:
            return self.parent().new('%s // %s'%(self.name(), x.name()))

    def __mod__(self, x):
        """
        Remainder of division of self by other.  This is denoted %.

        EXAMPLE::

            sage: R.<x,y> = GF(7)[]

        Now make the M2 version of R, so we can coerce elements of R to M2::

            sage: macaulay2(R)                              # optional - macaulay2
            ZZ
            --[x..y, Degrees => {2:1}, Heft => {1}, MonomialOrder => {MonomialSize => 16}, DegreeRank => 1]
             7                                                       {GRevLex => {2:1}  }
                                                                     {Position => Up    }
            sage: f = (x^3 + 2*y^2*x)^7; f                  # optional - macaulay2
            x^21 + 2*x^7*y^14
            sage: h = macaulay2(f); print(h)                # optional - macaulay2
             21     7 14
            x   + 2x y
            sage: f1 = (x^2 + 2*y*x)                        # optional - macaulay2
            sage: h1 = macaulay2(f1)                        # optional - macaulay2
            sage: f2 = (x^3 + 2*y*x)                        # optional - macaulay2
            sage: h2 = macaulay2(f2)                        # optional - macaulay2
            sage: h % [h1,h2]                               # optional - macaulay2
            -3x*y
            sage: u = h // [h1,h2]                          # optional - macaulay2
            sage: h == u[0]*h1 + u[1]*h2 + (h % [h1,h2])    # optional - macaulay2
            True
        """
        if isinstance(x, (list, tuple)):
            y = self.parent(x)
            return self.parent().new('%s %% matrix{%s}'%(self.name(), y.name()))
        if not isinstance(x, Macaulay2Element):
            x = self.parent(x)
        return self.parent().new('%s %% %s'%(self.name(), x.name()))

    def __nonzero__(self):
        """
        EXAMPLES::

            sage: a = macaulay2(0)  # optional - macaulay2
            sage: a == 0            # optional - macaulay2
            True
            sage: a.__nonzero__()   # optional - macaulay2
            False
        """
        P = self.parent()
        return P.eval('%s == 0'%self.name()) == 'false'

    def sage_polystring(self):
        """
        If this Macaulay2 element is a polynomial, return a string
        representation of this polynomial that is suitable for
        evaluation in Python.  Thus ``*`` is used for multiplication
        and ``**`` for exponentiation.   This function is primarily
        used internally.

        EXAMPLES::

            sage: R = macaulay2.ring('QQ','(x,y)')               # optional - macaulay2
            sage: f = macaulay2('x^3 + 3*y^11 + 5')              # optional - macaulay2
            sage: print(f)                                       # optional - macaulay2
             3     11
            x  + 3y   + 5
            sage: f.sage_polystring()                            # optional - macaulay2
            'x**3+3*y**11+5'
        """
        return self.external_string().replace('^','**')

    def structure_sheaf(self):
        """
        EXAMPLES::

            sage: S = macaulay2('QQ[a..d]')                     # optional - macaulay2
            sage: R = S/macaulay2('a^3+b^3+c^3+d^3')            # optional - macaulay2
            sage: X = R.Proj()                                  # optional - macaulay2
            sage: print(X.structure_sheaf())                    # optional - macaulay2
            OO
              sage...
        """
        return self.parent()('OO_%s'%self.name())

    def substitute(self, *args, **kwds):
        """
        Note that we have to override the substitute method so that we get
        the default one from Macaulay2 instead of the one provided by Element.

        EXAMPLES::

            sage: R = macaulay2("QQ[x]")            # optional - macaulay2
            sage: P = macaulay2("ZZ/7[symbol x]")   # optional - macaulay2
            sage: x, = R.gens()                     # optional - macaulay2
            sage: a = x^2 + 1                       # optional - macaulay2
            sage: a = a.substitute(P)               # optional - macaulay2
            sage: a.to_sage().parent()              # optional - macaulay2
            Univariate Polynomial Ring in x over Finite Field of size 7

        """
        return self.__getattr__("substitute")(*args, **kwds)

    subs = substitute

    def _tab_completion(self):
        """
        Return a list of tab completions for ``self``.

        Returns dynamically built sorted list of commands obtained using
        Macaulay2 "methods" command. All returned functions can take ``self``
        as their first argument

        Return type: list of strings

        TEST::

            sage: a = macaulay2("QQ[x,y]")   # optional - macaulay2
            sage: traits = a._tab_completion()   # optional - macaulay2
            sage: "generators" in traits     # optional - macaulay2
            True
        """
        # It is possible, that these are not all possible methods, but
        # there are still plenty and at least there are no definitely
        # wrong ones...
        r = self.parent().eval(
            """currentClass = class %s;
            total = {};
            while true do (
                -- Select methods with first argument of the given class
                r = select(methods currentClass, s -> s_1 === currentClass);
                -- Get their names as strings
                r = apply(r, s -> toString s_0);
                -- Keep only alpha-numeric ones
                r = select(r, s -> match("^[[:alnum:]]+$", s));
                -- Add to existing ones
                total = total | select(r, s -> not any(total, e -> e == s));
                if parent currentClass === currentClass then break;
                currentClass = parent currentClass;
                )
            toString total""" % self.name())
        r = sorted(r[1:-1].split(", "))
        return r

    def cls(self):
        """
        Since class is a keyword in Python, we have to use cls to call
        Macaulay2's class.  In Macaulay2, class corresponds to Sage's
        notion of parent.

        EXAMPLES::

            sage: macaulay2(ZZ).cls()  # optional - macaulay2
            Ring

        """
        return self.parent()("class %s"%self.name())

    ##########################
    #Aliases for M2 operators#
    ##########################
    def dot(self, x):
        """
        EXAMPLES::

            sage: d = macaulay2.new("MutableHashTable") # optional - macaulay2
            sage: d["k"] = 4                            # optional - macaulay2
            sage: d.dot("k")                            # optional - macaulay2
            4
        """
        parent = self.parent()
        x = parent(x)
        return parent("%s.%s"%(self.name(), x))

    def _operator(self, opstr, x):
        """
        Returns the infix binary operation specified by opstr applied
        to self and x.

        EXAMPLES::

            sage: a = macaulay2("3")     # optional - macaulay2
            sage: a._operator("+", a)    # optional - macaulay2
            6
            sage: a._operator("*", a)    # optional - macaulay2
            9
        """
        parent = self.parent()
        x = parent(x)
        return parent("%s%s%s"%(self.name(), opstr, x.name()))

    def sharp(self, x):
        """
        EXAMPLES::

            sage: a = macaulay2([1,2,3]) # optional - macaulay2
            sage: a.sharp(0)             # optional - macaulay2
            1
        """
        return self._operator("#", x)

    def starstar(self, x):
        """
        The binary operator ``**`` in Macaulay2 is usually used for tensor
        or Cartesian power.

        EXAMPLES::

            sage: a = macaulay2([1,2]).set()  # optional - macaulay2
            sage: a.starstar(a)               # optional - macaulay2
            set {(1, 1), (1, 2), (2, 1), (2, 2)}

        """
        return self._operator("**", x)

    def underscore(self, x):
        """
        EXAMPLES::

            sage: a = macaulay2([1,2,3])  # optional - macaulay2
            sage: a.underscore(0)         # optional - macaulay2
            1
        """
        return self._operator("_", x)

    ####################
    #Conversion to Sage#
    ####################
    def to_sage(self):
        """
        EXAMPLES::

            sage: macaulay2(ZZ).to_sage()      # optional - macaulay2
            Integer Ring
            sage: macaulay2(QQ).to_sage()      # optional - macaulay2
            Rational Field

            sage: macaulay2(2).to_sage()       # optional - macaulay2
            2
            sage: macaulay2(1/2).to_sage()     # optional - macaulay2
            1/2
            sage: macaulay2(2/1).to_sage()     # optional - macaulay2
            2
            sage: _.parent()                   # optional - macaulay2
            Rational Field
            sage: macaulay2([1,2,3]).to_sage() # optional - macaulay2
            [1, 2, 3]

            sage: m = matrix([[1,2],[3,4]])
            sage: macaulay2(m).to_sage()       # optional - macaulay2
            [1 2]
            [3 4]

            sage: macaulay2(QQ['x,y']).to_sage()    # optional - macaulay2
            Multivariate Polynomial Ring in x, y over Rational Field
            sage: macaulay2(QQ['x']).to_sage()      # optional - macaulay2
            Univariate Polynomial Ring in x over Rational Field
            sage: macaulay2(GF(7)['x,y']).to_sage() # optional - macaulay2
            Multivariate Polynomial Ring in x, y over Finite Field of size 7

            sage: macaulay2(GF(7)).to_sage()       # optional - macaulay2
            Finite Field of size 7
            sage: macaulay2(GF(49, 'a')).to_sage() # optional - macaulay2
            Finite Field in a of size 7^2

            sage: R.<x,y> = QQ[]
            sage: macaulay2(x^2+y^2+1).to_sage()   # optional - macaulay2
            x^2 + y^2 + 1

            sage: R = macaulay2("QQ[x,y]")         # optional - macaulay2
            sage: I = macaulay2("ideal (x,y)")     # optional - macaulay2
            sage: I.to_sage()                      # optional - macaulay2
            Ideal (x, y) of Multivariate Polynomial Ring in x, y over Rational Field

            sage: X = R/I       # optional - macaulay2
            sage: X.to_sage()   # optional - macaulay2
            Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x, y)

            sage: R = macaulay2("QQ^2")  # optional - macaulay2
            sage: R.to_sage()            # optional - macaulay2
            Vector space of dimension 2 over Rational Field

            sage: m = macaulay2('"hello"')  # optional - macaulay2
            sage: m.to_sage()               # optional - macaulay2
            'hello'

        """
        repr_str = str(self)
        cls_str = str(self.cls())
        cls_cls_str = str(self.cls().cls())

        if repr_str == "ZZ":
            from sage.rings.all import ZZ
            return ZZ
        elif repr_str == "QQ":
            from sage.rings.all import QQ
            return QQ

        if cls_cls_str == "Type":
            if cls_str == "List":
                return [entry.to_sage() for entry in self]
            elif cls_str == "Matrix":
                from sage.matrix.all import matrix
                base_ring = self.ring().to_sage()
                entries = self.entries().to_sage()
                return matrix(base_ring, entries)
            elif cls_str == "Ideal":
                parent = self.ring().to_sage()
                gens = self.gens().entries().flatten().to_sage()
                return parent.ideal(*gens)
            elif cls_str == "QuotientRing":
                #Handle the ZZ/n case
                if "ZZ" in repr_str and "--" in repr_str:
                    from sage.rings.all import ZZ, GF
                    external_string = self.external_string()
                    zz, n = external_string.split("/")

                    #Note that n must be prime since it is
                    #coming from Macaulay 2
                    return GF(ZZ(n))

                ambient = self.ambient().to_sage()
                ideal = self.ideal().to_sage()
                return ambient.quotient(ideal)
            elif cls_str == "PolynomialRing":
                from sage.rings.all import PolynomialRing
                from sage.rings.polynomial.term_order import inv_macaulay2_name_mapping

                #Get the base ring
                base_ring = self.coefficientRing().to_sage()

                #Get a string list of generators
                gens = str(self.gens())[1:-1]

                # Check that we are dealing with default degrees, i.e. 1's.
                if self.degrees().any("x -> x != {1}").to_sage():
                    raise ValueError("cannot convert Macaulay2 polynomial ring with non-default degrees to Sage")
                #Handle the term order
                external_string = self.external_string()
                order = None
                if "MonomialOrder" not in external_string:
                    order = "degrevlex"
                else:
                    for order_name in inv_macaulay2_name_mapping:
                        if order_name in external_string:
                            order = inv_macaulay2_name_mapping[order_name]
                if len(gens) > 1 and order is None:
                    raise ValueError("cannot convert Macaulay2's term order to a Sage term order")

                return PolynomialRing(base_ring, order=order, names=gens)
            elif cls_str == "GaloisField":
                from sage.rings.all import ZZ, GF
                gf, n = repr_str.split(" ")
                n = ZZ(n)
                if n.is_prime():
                    return GF(n)
                else:
                    gen = str(self.gens())[1:-1]
                    return GF(n, gen)
            elif cls_str == "Boolean":
                if repr_str == "true":
                    return True
                elif repr_str == "false":
                    return False
            elif cls_str == "String":
                return str(repr_str)
            elif cls_str == "Module":
                from sage.modules.all import FreeModule
                if self.isFreeModule().to_sage():
                    ring = self.ring().to_sage()
                    rank = self.rank().to_sage()
                    return FreeModule(ring, rank)
        else:
            #Handle the integers and rationals separately
            if cls_str == "ZZ":
                from sage.rings.all import ZZ
                return ZZ(repr_str)
            elif cls_str == "QQ":
                from sage.rings.all import QQ
                repr_str = self.external_string()
                if "/" not in repr_str:
                    repr_str = repr_str + "/1"
                return QQ(repr_str)

            m2_parent = self.cls()
            parent = m2_parent.to_sage()

            if cls_cls_str == "PolynomialRing":
                from sage.misc.sage_eval import sage_eval
                gens_dict = parent.gens_dict()
                return sage_eval(self.external_string(), gens_dict)

        from sage.misc.sage_eval import sage_eval
        try:
            return sage_eval(repr_str)
        except Exception:
            raise NotImplementedError("cannot convert %s to a Sage object"%repr_str)


class Macaulay2Function(ExpectFunction):
    def _sage_doc_(self):
        """
        EXAMPLES::

            sage: print(macaulay2.load._sage_doc_())  # optional - macaulay2
            load -- read Macaulay2 commands
            *******************************
            ...
              * "input" -- read Macaulay2 commands and echo
              * "notify" -- whether to notify the user when a file is loaded
        """
        return self._parent.help(self._name)

    def _sage_src_(self):
        """
        EXAMPLES::

            sage: print(macaulay2.gb._sage_src_())  # optional - macaulay2
            code(methods gb)
            ...
        """
        if self._parent._expect is None:
            self._parent._start()
        E = self._parent._expect
        E.sendline("code(methods %s)"%self._name)
        E.expect(self._parent._prompt)
        s = E.before
        self._parent.eval("2+2")
        return s

def is_Macaulay2Element(x):
    """
    EXAMPLES::

        sage: from sage.interfaces.macaulay2 import is_Macaulay2Element
        sage: is_Macaulay2Element(2)              # optional - macaulay2
        False
        sage: is_Macaulay2Element(macaulay2(2))   # optional - macaulay2
        True
    """
    return isinstance(x, Macaulay2Element)

# An instance
macaulay2 = Macaulay2()


def macaulay2_console():
    """
    Spawn a new M2 command-line session.

    EXAMPLES::

        sage: macaulay2_console()                    # not tested
        Macaulay 2, version 1.1
        with packages: Classic, Core, Elimination, IntegralClosure, LLLBases, Parsing, PrimaryDecomposition, SchurRings, TangentCone
        ...

    """
    from sage.repl.rich_output.display_manager import get_display_manager
    if not get_display_manager().is_in_terminal():
        raise RuntimeError('Can use the console only in the terminal. Try %%macaulay2 magics instead.')
    os.system('M2')



def reduce_load_macaulay2():
    """
    Used for reconstructing a copy of the Macaulay2 interpreter from a pickle.

    EXAMPLES::

        sage: from sage.interfaces.macaulay2 import reduce_load_macaulay2
        sage: reduce_load_macaulay2()
        Macaulay2
    """
    return macaulay2

