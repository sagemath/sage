r"""
Interface to Risa / Asir
"""

# $OpenXM: OpenXM/src/sage/asir.py,v 1.2 2019/03/06 02:38:33 takayama Exp $
import os
import pexpect

from sage.cpython.string import bytes_to_str
from sage.interfaces.expect import Expect, ExpectElement
from sage.modules.free_module_element import vector
from sage.matrix.matrix_space import MatrixSpace
from sage.matrix.constructor import matrix
from sage.misc.verbose import verbose
from sage.rings.cc import CC
from sage.rings.qqbar import QQbar
from sage.rings.rational_field import QQ
from sage.rings.real_mpfr import RR

# Ref: @s/2018/09/20180907-sage-asir-proj,
# Using External Libraries and Interfaces

# conversion according to OpenXM doc
types = {-1: "void",
         0: "nothing",
         1: "number",
         2: "polynomial",
         3: "rational expression",
         4: "list",
         5: "vector",
         6: "matrix",
         7: "string",
         8: "structure",
         9: "distributed polynomial",
         10: "32bit unsigned integer",
         11: "error object",
         12: "matrix over GF(2)",
         13: "MATHCAP object",
         14: "first order formula",
         15: "matrix over GF(p)",
         16: "byte array",
         26: "distributed module polynomial",
         }

# http://www.math.sci.kobe-u.ac.jp/OpenXM/Current/doc/asir2000/html-en/man/man_35.html#ntype
number_types = {0: "rational",
                1: "real",
                2: "algebraic",
                3: "real",
                4: "complex",
                5: "finitefield",
                6: "finitefield",
                7: "finitefield"
                }


class Asir(Expect):
    r"""
    Interface to the Asir interpreter.

    EXAMPLES::

        sage: asir.evall("F=fctr(x^10-1)")    # optional - asir
        '[[1,1],[x-1,1],[x+1,1],[x^4-x^3+x^2-x+1,1],[x^4+x^3+x^2+x+1,1]]'
    """
    def __init__(self, maxread=None, script_subdirectory=None, logfile=None,
                 server=None, server_tmpdir=None, seed=None, command=None):
        """
        EXAMPLES::

            sage: asir == loads(dumps(asir))
            True
        """
        texmacs = 'openxm ox_texmacs --view sage --quiet --noCopyright'
        if command is None:
            command = os.getenv('SAGE_ASIR_COMMAND') or texmacs
        if server is None:
            server = os.getenv('SAGE_ASIR_SERVER') or None
        Expect.__init__(self,
                        name='asir',
                        # We want the prompt sequence to be unique to avoid
                        # confusion with syntax error messages containing >>>
                        prompt='asir>',
                        # We don't want any pagination of output
                        command=texmacs,
                        maxread=maxread,
                        server=server,
                        server_tmpdir=server_tmpdir,
                        script_subdirectory=script_subdirectory,
                        restart_on_ctrlc=False,
                        verbose_start=False,
                        logfile=logfile,
                        eval_using_file_cutoff=100)
        self._seed = seed

    def set_seed(self, seed=None):
        """
        Not implemented. Set the seed for the random number generator
        for this asir interpreter.
        """
        return 0

    def __reduce__(self):
        """
        EXAMPLES::

            sage: Asir().__reduce__()
            (<function reduce_load_Asir at 0x...>, ())
        """
        return reduce_load_Asir, tuple()

    def _read_in_file_command(self, filename):
        """
        EXAMPLES::

            sage: filename = tmp_filename()
            sage: asir._read_in_file_command(filename)
            'load("...");'
        """
        return 'load("%s");' % filename

    def _quit_string(self):
        """
        EXAMPLES::

            sage: asir._quit_string()
            'quit();'
        """
        return 'quit();'

    def _install_hints(self):
        """
        Return hints on how to install Asir.

        EXAMPLES::

            sage: print(asir._install_hints())
            You must get ...
        """
        return """
        You must get the program "asir" and "ox_texmacs" in order to use Asir
        from Sage. You can read all about Asir at
                http://www.openxm.org
        The command "openxm" must be in the search path.
        """

    def evall(self, cmd):
        """
        Evaluate the argument immediately. The argument of eval is buffered
        and ; should be added.

        EXAMPLES::

            sage: asir.eval('1+2;'); asir.evall('3+3')   #optional - asir
            '3'
            '6'
        """
        return self.eval(cmd + ';;')

    def _default_var_name(self):
        # because the usual "sage" fails
        return "Sage"

    def _eval_line(self, line, reformat=True, allow_use_file=False,
                   wait_for_prompt=True, restart_if_needed=False):
        """
        EXAMPLES::

            sage: print(asir._eval_line('2+2'))  #optional - asir
            4
        """
        from pexpect.exceptions import EOF
        if not wait_for_prompt:
            return Expect._eval_line(self, line)
        if line == '':
            return ''
        if self._expect is None:
            self._start()
        if allow_use_file and len(line) > 3000:
            return self._eval_line_using_file(line)
        try:
            E = self._expect
            # debug
            # self._synchronize(cmd='1+%s\n')
            verbose("in = '%s'" % line, level=3)
            E.sendline(line)
            E.expect(self._prompt)
            out = bytes_to_str(E.before)
            # debug
            verbose("out = '%s'" % out, level=3)
        except EOF:
            if self._quit_string() in line:
                return ''
        except KeyboardInterrupt:
            self._keyboard_interrupt()
        try:
            if reformat:
                if 'syntax error' in out:
                    raise SyntaxError(out)
            out = "\n".join(out.splitlines()[1:])
            return out
        except NameError:
            return ''

    def _keyboard_interrupt(self):
        print("Ctrl-C: Interrupting %s..." % self)
        if self._restart_on_ctrlc:
            try:
                self._expect.close(force=1)
            except pexpect.ExceptionPexpect as msg:
                raise RuntimeError("THIS IS A BUG -- PLEASE REPORT."
                                   " This should never happen.\n" + msg)
            self._start()
            raise KeyboardInterrupt("Restarting %s (WARNING: all variables defined in previous session are now invalid)" % self)
        else:
            self._expect.send('\003')  # control-c
            raise KeyboardInterrupt("Ctrl-c pressed while running %s" % self)

    def quit(self, verbose=False):
        """
        EXAMPLES::

            sage: # optional - asir
            sage: o = Asir()
            sage: o._start()
            sage: o.quit(True)
            Exiting spawned Asir process.
        """
        # Don't bother, since it just hangs in some cases, and it
        # isn't necessary, since asir behaves well with respect
        # to signals.
        if self._expect is not None:
            if verbose:
                print("Exiting spawned %s process." % self)

    def _start(self):
        """
        Start the Asir process.

        EXAMPLES::

            sage: # optional - asir
            sage: o = Asir()
            sage: o.is_running()
            False
            sage: o._start()
            sage: o.is_running()
            True
        """
        Expect._start(self)
#        self.eval("page_screen_output=0;")
#        self.eval("format none;")
        # set random seed
#        self.set_seed(self._seed)

    def _equality_symbol(self):
        """
        EXAMPLES::

            sage: # optional - asir
            sage: asir('0 == 1')
            0
            sage: asir('1 == 1')
            1
        """
        return '=='

    def _true_symbol(self):
        """
        EXAMPLES::

            sage: asir('1 == 1')  # optional - asir
            1
        """
        return '1'

    def _false_symbol(self):
        """
        EXAMPLES::

            sage: asir('0 == 1')  # optional - asir
            0
        """
        return '0'

    def set(self, var, value):
        """
        Set the variable ``var`` to the given ``value``.

        EXAMPLES::

            sage: # optional - asir
            sage: asir.set('X', '2')
            sage: asir.get('X')
            '2'
        """
        cmd = '%s=%s' % (var, value)
        out = self.evall(cmd)
        if out.find("error") != -1 or out.find("Error") != -1:
            raise TypeError("Error executing code in Asir\nCODE:\n\t%s\nAsir ERROR:\n\t%s" % (cmd, out))

    def get(self, var):
        """
        Get the value of the variable ``var``.

        EXAMPLES::

            sage: asir.set('X', '2') # optional - asir
            sage: asir.get('X') # optional - asir
            '2'
        """
        s = self.evall('%s;' % var)
        i = s.find('=')
        return s[i + 1:]

    def console(self):
        """
        Spawn a new Asir command-line session.

        This requires that the optional ``asir`` program be installed and in
        your ``PATH``, but no optional Sage packages need be installed.

        EXAMPLES::

            sage: asir_console()   # not tested
            This is Risa/Asir, ....
            ...
            [nnnn] 2+3;
            5
            [nnnn] quit();

        The command ``quit();`` exits the ``asir`` console and returns you to Sage.
        """
        asir_console()

    def version(self) -> str:
        """
        Return the version of Asir.

        OUTPUT: string

        EXAMPLES::

            sage: v = asir.version()   # optional - asir
            sage: v                      # optional - asir; random
            '20240806'
        """
        return str(self("version()")).strip()

    def _object_class(self):
        """
        EXAMPLES::

            sage: asir._object_class()
            <class 'sage.interfaces.asir.AsirElement'>
        """
        return AsirElement


asir_functions = set()


# TODO the following class has not yet been fully implemented.
class AsirElement(ExpectElement):
    def _get_sage_ring(self):
        r"""
        TESTS::

            sage: # optional - asir
            sage: asir('1')._get_sage_ring()
            Rational Field
            sage: asir('@i')._get_sage_ring()
            Complex Double Field
            sage: asir('1.2')._get_sage_ring()
            Real Double Field
        """
        if self.asir_type() == 'rational':
            import sage.rings.rational_field
            return sage.rings.rational_field.QQ
        if self.asir_type() == 'real':
            import sage.rings.real_double
            return sage.rings.real_double.RDF
        if self.asir_type() == 'complex':
            import sage.rings.complex_double
            return sage.rings.complex_double.CDF
        if self.asir_type() == 'algebraic':
            import sage.rings.qqbar
            return sage.rings.qqbar.QQbar
        if self.asir_type() == 'finitefield':
            raise TypeError('finite fields are not handled yet')
        if self.asir_type() == 'polynomial':
            raise TypeError('polynomials are not handled yet')
        if self.asir_type() == 'rational expression':
            raise TypeError('rational expression are not handled yet')
        raise TypeError("no Sage ring associated to this element.")

    def __bool__(self) -> bool:
        r"""
        Test whether this element is nonzero.

        EXAMPLES::

            sage: # optional - asir
            sage: bool(asir('0'))
            False
            sage: bool(asir('[]'))
            False
            sage: bool(asir('[0,0]'))
            True
            sage: bool(asir('newmat(3,2,[[0,0,0],[0,0,0]])'))
            False

            sage: # optional - asir
            sage: bool(asir('0.1'))
            True
            sage: bool(asir('[0,1,0]'))
            True
            sage: bool(asir('newmat(3,2,[[0,0,-0.1],[0,0,0]])'))
            True
        """
        # TODO better
        if self.asir_type() == 'list':
            return int(self.length()) != 0

        s = str(self)
        return s != ' [](0x0)' and any(x != '0' for x in s.split())

    def _matrix_(self, R=None):
        r"""
        Return Sage matrix from this ``asir`` element.

        EXAMPLES::

            sage: # optional - asir
            sage: A = asir('newmat(2,2,[[1,2],[3,9/2]])')
            sage: matrix(A)
            [  1   2]
            [  3 9/2]
            sage: _.base_ring()
            Rational Field

            sage: # optional - asir
            sage: A = asir('newmat(2,2,[[@i,1],[-1,0]])')
            sage: matrix(A)
            [1.0*I   1.0]
            [ -1.0   0.0]
            sage: _.base_ring()
            Complex Double Field

            sage: # optional - asir
            sage: A = asir('newmat(2,2,[[1,2],[3,4]])')
            sage: matrix(ZZ, A)
            [1 2]
            [3 4]
            sage: A = asir('newmat(2,2,[[1,2],[3,4.5]])')
            sage: matrix(RR, A)
            [1.00000000000000 2.00000000000000]
            [3.00000000000000 4.50000000000000]
        """
        nrows, ncols = self.size().sage()

        w = [row.sage() for row in self]

        if R is None:
            return matrix(nrows, ncols, w)

        return MatrixSpace(R, nrows, ncols)(w)

    def _vector_(self, R=None):
        r"""
        Return Sage vector from this asir element.

        EXAMPLES::

            sage: # optional - asir
            sage: A = asir('[1,2,3,4]')
            sage: vector(ZZ, A)
            (1, 2, 3, 4)
            sage: A = asir('[1,2.3,4.5]')
            sage: vector(A)
            (1.00000000000000, 2.30000000000000, 4.50000000000000)
            sage: A = asir('[1,@i]')
            sage: vector(A)
            (1.00000000000000, 1.00000000000000*I)
        """
        w = [x.sage() for x in self]
        if R is None:
            return vector(w)
        nrows = int(self.length())
        from sage.modules.free_module import FreeModule
        return FreeModule(R, nrows)(w)

    def asir_type(self):
        try:
            number = int(str(self.type()))
        except ValueError:
            number = int(str(self.type()).splitlines()[-1])
        typ = types[number]
        if typ != "number":
            return typ
        try:
            number = int(str(self.ntype()))
        except ValueError:
            number = int(str(self.ntype()).splitlines()[-1])
        return number_types[number]

    def __iter__(self):
        typ = self.asir_type()
        if typ == "list":
            L = int(self.length())
            for i in range(L):
                yield self[i]
        elif typ == "vector":
            L = int(self.length())
            for i in range(L):
                yield self[i]
        elif typ == "matrix":
            m, n = self.size()
            m = int(m)
            for i in range(m):
                yield self[i]
        else:
            raise TypeError(f'can iterate over lists, vectors, matrices, not {typ}')

    def _sage_(self):
        """
        Try to parse the asir object and return a sage object.

        EXAMPLES::

            sage: # optional - asir
            sage: A = asir('2833')
            sage: A.sage()
            2833
            sage: B = A**2
            sage: B.sage()
            8025889
            sage: C = asir('3*@i')
            sage: C.sage()
            3.00000000000000*I
            sage: A = asir('[1,2,3,4]')
            sage: A.sage()
            [1, 2, 3, 4]
            sage: A = asir('[1,2.3,4.5]')
            sage: A.sage()
            [1, 2.30000000000000, 4.50000000000000]
            sage: A = asir('[1,2.3*@i,4.5]')
            sage: A.sage()
            [1, 2.30000000000000*I, 4.50000000000000]

            sage: # optional - asir
            sage: A = asir('2833')
            sage: As = A.sage(); As
            2833
            sage: As.parent()
            Rational Field
        """
        if self.asir_type() == "rational":
            return QQ(str(self))
        if self.asir_type() == "real":
            return RR(str(self))
        if self.asir_type() == "complex":
            return CC((str(self.real()), str(self.imag())))
        if self.asir_type() == "algebraic":
            return QQbar(self)
        if self.asir_type() == "finitefield":
            raise NotImplementedError
        if self.asir_type() == "vector":
            return self._vector_()
        if self.asir_type() == "matrix":
            return self._matrix_()
        if self.asir_type() == "list":
            return [elt.sage() for elt in self]
        raise NotImplementedError(f'asir type {self.asir_type()} is not yet recognized')


# An instance
asir = Asir()


def reduce_load_Asir():
    """
    Return the Asir interface object defined in
    sage.interfaces.asir.

    EXAMPLES::

        sage: from sage.interfaces.asir import reduce_load_Asir
        sage: reduce_load_Asir()
        Asir
    """
    return asir


def asir_console():
    """
    Spawn a new Asir command-line session.

    This requires that the optional ``asir`` program be installed and in
    your ``PATH``, but no optional Sage packages need be installed.

    EXAMPLES::

        sage: asir_console()         # not tested
        This is Risa/Asir ....
        ...
        [nnnn] 2+3;
        5
        [nnnn] quit();

    quit(); exits the asir console and returns you to Sage.
    """
    from sage.repl.rich_output.display_manager import get_display_manager
    if not get_display_manager().is_in_terminal():
        raise RuntimeError('Can use the console only in the terminal. '
                           'Try %%asir magics instead.')
    os.system('openxm fep asir')    # with asir prompt
#    os.system('openxm asir -quiet')
