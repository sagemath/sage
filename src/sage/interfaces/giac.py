r"""
Interface to Giac

(adapted by F. Han from William Stein and Gregg Musiker maple's interface)

You must have the optional  Giac interpreter installed
and available as the command ``giac`` in your PATH in
order to use this interface. You need a giac version
supporting "giac --sage" ( roughly after 0.9.1 ). In this case you do not have
to install any  optional Sage packages. If giac is not already installed, you can
download binaries or sources or spkg (follow the sources link) from the homepage:

Homepage <http://www-fourier.ujf-grenoble.fr/~parisse/giac.html>

Type ``giac.[tab]`` for a list of all the functions
available from your Giac install. Type
``giac.[tab]?`` for Giac's help about a given
function. Type ``giac(...)`` to create a new Giac
object, and ``giac.eval(...)`` to run a string using
Giac (and get the result back as a string).

If the giac spkg is installed, you should find the full html documentation there::

    $SAGE_LOCAL/share/giac/doc/en/cascmd_local/index.html


EXAMPLES::

    sage: giac('3 * 5')                                 # optional - giac
    15
    sage: giac.eval('ifactor(2005)')                    # optional - giac
    '5*401'
    sage: giac.ifactor(2005)                            # optional - giac
    2005
    sage: l=giac.ifactors(2005) ; l; l[2]               # optional - giac
    [5,1,401,1]
    401
    sage: giac.fsolve('x^2=cos(x)+4', 'x','0..5')         # optional - giac
    1.9140206190...
    sage: giac.factor('x^5 - y^5')                      # optional - giac
    (x-y)*(x^4+x^3*y+x^2*y^2+x*y^3+y^4)
    sage: R.<x,y>=QQ[];f=(x+y)^5;f2=giac(f);(f-f2).normal() #optional - giac
    0
    sage: x,y=giac('x,y'); giac.int(y/(cos(2*x)+cos(x)),x)         #random; optional - giac
    y*2*((-(tan(x/2)))/6+(-2*1/6/sqrt(3))*ln(abs(6*tan(x/2)-2*sqrt(3))/abs(6*tan(x/2)+2*sqrt(3))))


If the string "error" (case insensitive) occurs in the output of
anything from Giac, a RuntimeError exception is raised.

Tutorial
--------

AUTHORS:

- Gregg Musiker (2006-02-02): initial version.

(adapted to giac by F.Han)

This tutorial is based on the Maple Tutorial for number theory from
http://www.math.mun.ca/~drideout/m3370/numtheory.html.

There are several ways to use the Giac Interface in Sage. We will
discuss two of those ways in this tutorial.


#. If you have a giac expression such as

   ::

       factor( (x^5-1));

   We can write that in sage as

   ::

       sage: giac('factor(x^5-1)')                 # optional - giac
       (x-1)*(x^4+x^3+x^2+x+1)

   Notice, there is no need to use a semicolon.

#. Since Sage is written in Python, we can also import giac
   commands and write our scripts in a pythonic way. For example,
   ``factor()`` is a giac command, so we can also factor
   in Sage using

   ::

       sage: giac('(x^5-1)').factor()              # optional - giac
       (x-1)*(x^4+x^3+x^2+x+1)

   where ``expression.command()`` means the same thing as
   ``command(expression)`` in Giac. We will use this
   second type of syntax whenever possible, resorting to the first
   when needed.

   ::

       sage: giac('(x^12-1)/(x-1)').normal()     # optional - giac
       x^11+x^10+x^9+x^8+x^7+x^6+x^5+x^4+x^3+x^2+x+1


The normal command will reduce a rational function to the
lowest terms. In giac, simplify is slower than normal because it
tries more sophisticated simplifications (ex algebraic extensions)
The factor command will factor a polynomial with
rational coefficients into irreducible factors over the ring of
integers (if your default configuration of giac (cf .xcasrc) has not
allowed square roots). So for example,


::

    sage: giac('(x^12-1)').factor( )           # optional - giac
    (x-1)*(x+1)*(x^2+1)*(x^2-x+1)*(x^2+x+1)*(x^4-x^2+1)

::

    sage: giac('(x^28-1)').factor( )           # optional - giac
    (x-1)*(x+1)*(x^2+1)*(x^6-x^5+x^4-x^3+x^2-x+1)*(x^6+x^5+x^4+x^3+x^2+x+1)*(x^12-x^10+x^8-x^6+x^4-x^2+1)

Another important feature of giac is its online help. We can
access this through sage as well. After reading the description of
the command, you can press q to immediately get back to your
original prompt.

Incidentally you can always get into a giac console by the
command

::

    sage: giac.console()                       # not tested
    sage: !giac                                # not tested

Note that the above two commands are slightly different, and the
first is preferred.

For example, for help on the giac command factors, we type ::

    sage: giac.help('factors')                     # not tested

::

    sage: alpha = giac((1+sqrt(5))/2)            # optional - giac
    sage: beta = giac(1-sqrt(5))/2               # optional - giac
    sage: f19  = alpha^19 - beta^19/sqrt(5)      # optional - giac
    sage: f19                                    # optional - giac
    (sqrt(5)/2+1/2)^19-((-sqrt(5)+1)/2)^19/sqrt(5)
    sage: (f19-(5778*sqrt(5)+33825)/5).normal()                           # optional - giac
    0

Let's say we want to write a giac program now that squares a
number if it is positive and cubes it if it is negative. In giac,
that would look like

::

    mysqcu := proc(x)
    if x > 0 then x^2;
    else x^3; fi;
    end;

In Sage, we write

::

    sage: mysqcu = giac('proc(x) if x > 0 then x^2 else x^3 fi end')    # optional - giac
    sage: mysqcu(5)                                                      # optional - giac
    25
    sage: mysqcu(-5)                                                     # optional - giac
    -125

More complicated programs should be put in a separate file and
loaded.
"""

#############################################################################
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#############################################################################

import os

from sage.interfaces.expect import Expect, ExpectElement, ExpectFunction, FunctionElement, gc_disabled

import pexpect

from sage.env import DOT_SAGE
from sage.misc.pager import pager

COMMANDS_CACHE = '%s/giac_commandlist_cache.sobj'%DOT_SAGE

class Giac(Expect):
    r"""
    Interface to the Giac interpreter.

    You must have the optional  Giac interpreter installed and available as the command ``giac`` in your PATH in order to use this interface. Try the command: print giac._install_hints() for more informations on giac installation.

    Type ``giac.[tab]`` for a list of all the functions available from your Giac install.
    Type ``giac.[tab]?`` for Giac's help about a given function.
    Type ``giac(...)`` to create a new Giac object.

    Full html documentation for giac is avaible from your giac installation at ``$PREFIX``/share/giac/doc/en/cascmd_en/index.html

    EXAMPLES:

    Any Giac instruction can be evaluated as a string by the giac command. You can access the giac functions by adding the ``giac.`` prefix to the usual Giac name.

    ::

      sage: l=giac('normal((y+sqrt(2))^4)'); l   # optional - giac
      y^4+4*sqrt(2)*y^3+12*y^2+8*sqrt(2)*y+4
      sage: f=giac('(u,v)->{ if (u<v){ [u,v] } else { [v,u] }}');f(1,2),f(3,1)   # optional - giac
      ([1,2], [1,3])

    The output of the giac command is a Giac object, and it can be used for another giac command.

    ::

      sage: l.factors()                          #optional  - giac
      [y+sqrt(2),4]
      sage: giac('(x^12-1)').factor( )           # optional - giac
      (x-1)*(x+1)*(x^2+1)*(x^2-x+1)*(x^2+x+1)*(x^4-x^2+1)
      sage: giac('assume(y>0)'); giac('y^2=3').solve('y')  #optional - giac
      y
      [sqrt(3)]

    You can create some Giac elements and avoid many quotes like this:

    ::

      sage: x,y,z=giac('x,y,z');type(y)   # optional - giac
      <class 'sage.interfaces.giac.GiacElement'>
      sage: I1=(1/(cos(2*y)+cos(y))).integral(y,0,pi/4).simplify()  #optional - giac
      sage: (I1-((-2*ln((sqrt(3)-3*tan(1/8*pi))/(sqrt(3)+3*tan(1/8*pi)))*sqrt(3)-3*tan(1/8*pi))/9)).normal()       # optional - giac
      0
      sage: ((y+z*sqrt(5))*(y-sqrt(5)*z)).normal() # optional - giac
      y^2-5*z^2

    Polynomials or elements of SR can be evaluated directly by the giac interface.

    ::

      sage: R.<a,b>=QQ[];f=(2+a+b);p=giac.gcd(f^3+5*f^5,f^2+f^5);p;R(p); #optional - giac
      a^2+2*a*b+4*a+b^2+4*b+4
      a^2 + 2*a*b + b^2 + 4*a + 4*b + 4

    Variable names in python and giac are independant.

    ::

      sage: a=sqrt(2);giac('Digits:=30;a:=5');a,giac('a'),giac(a),giac(a).evalf()  # optional - giac
      [...]
      (sqrt(2), 5, sqrt(2), 1.414213562373095048801688724209)


    """
    def __init__(self, maxread=10000, script_subdirectory=None, server=None, server_tmpdir=None, logfile=None):
        """
        Create an instance of the Giac interpreter.

        EXAMPLES::

            sage: giac == loads(dumps(giac))                # optional - giac
            True
        """
        Expect.__init__(self,
                        name = 'giac',
                        prompt = '[0-9]*>> ',
                        command = "giac --sage",
                        init_code= ['maple_mode(0);I:=i;'],      #  coercion could be broken in maple_mode
                        maxread = maxread,
                        script_subdirectory = script_subdirectory,
                        restart_on_ctrlc = False,                        server = server,
                        server_tmpdir = server_tmpdir,
                        verbose_start = False,
                        logfile = logfile,
                        eval_using_file_cutoff=1000)

    def _function_class(self):
        """
        EXAMPLES::

            sage: giac._function_class()                 # optional - giac
            <class 'sage.interfaces.giac.GiacFunction'>

        ::

            sage: type(giac.diff)                        # optional - giac
            <class 'sage.interfaces.giac.GiacFunction'>
        """
        return GiacFunction

    def _keyboard_interrupt(self):
        """
        TESTS:

        We check to make sure that the gap interface behaves correctly
        after a keyboard interrupt.

            sage: giac(2)                                # optional - giac
            2
            sage: try:                                   # optional - giac
            ...     giac._keyboard_interrupt()
            ... except KeyboardInterrupt:
            ...     pass
            ...
            Interrupting Giac...
            sage: giac(2)                                # optional - giac
            2
        """
        print "Interrupting %s..."%self
        self._expect.sendline(chr(3))  # send ctrl-c
        self._expect.expect(self._prompt)
#        self._expect.expect(self._prompt)
        raise RuntimeError("Ctrl-c pressed while running %s"%self)

    def __reduce__(self):
        """
        EXAMPLES::

            sage: giac.__reduce__()
            (<function reduce_load_Giac at 0x...>, ())
            sage: f, args = _
            sage: f(*args)
            Giac
        """
        return reduce_load_Giac, tuple([])

    def _read_in_file_command(self, filename):
        r"""
        Returns the string used to read filename into Giac.

        EXAMPLES::

            sage: giac._read_in_file_command('test')   # optional - giac
            'read "test"'

        ::

            sage: filename = tmp_filename()             # optional - giac
            sage: f = open(filename,'w')                # optional - giac
            sage: f.write('xx := 22;\n')                # optional - giac
            sage: f.close()              # optional - giac
            sage: giac.read(filename)    # optional - giac
            sage: giac.get('xx').strip() # optional - giac
            '22'
        """
        return 'read "%s"'%filename

    def _quit_string(self):
        """
        EXAMPLES::

            sage: giac._quit_string()     # optional - giac
            '@d'

        ::

            sage: m = Giac()
            sage: a = m(2)           # optional - giac
            sage: m.is_running()     # optional - giac
            True
            sage: m.quit()           # optional - giac
            sage: m.is_running()     # optional - giac
            False
        """
        return '@d'

    def _install_hints(self):
        """
        Hints for installing Giac on your computer.

        EXAMPLES::

            sage: print giac._install_hints()
            In order...
        """
        return r"""

In order to use the Giac interface you need to have Giac installed
and have a program called "giac" in your PATH. You need a giac version
supporting "giac --sage" ( roughly after 0.9.1 of march 2011). Some giac
instructions  and the help's langage depend of you LANG variable. To obtain
inline help for  giac commands, you also need to have the program "cas_help"
in your PATH.


If giac is not already installed, you can download binaries or sources
or a spkg ( for the spkg follow the sources link) from the homepage:

Homepage http://www-fourier.ujf-grenoble.fr/~parisse/giac.html


Full html documentation for giac is avaible from your giac installation at:

    ``$PREFIX``/share/giac/doc/en/cascmd_en/index.html

If you got giac from the spkg then ``$PREFIX`` is ``$SAGE_LOCAL``

"""

    def expect(self):
        """
        Returns the pexpect object for this Giac session.

        EXAMPLES::

            sage: m = Giac()
            sage: m.expect() is None
            True
            sage: m._start()           # optional - giac
            sage: m.expect()           # optional - giac
            <pexpect.spawn instance at 0x...>
            sage: m.quit()             # optional - giac
        """
        return self._expect

    def console(self):
        """
        Spawn a new Giac command-line session.

        EXAMPLES::

            sage: giac_console()                   # not tested - giac
            ...
            Homepage http://www-fourier.ujf-grenoble.fr/~parisse/giac.html
            Released under the GPL license 3.0 or above
            See http://www.gnu.org for license details
            -------------------------------------------------
            Press CTRL and D simultaneously to finish session
            Type ?commandname for help
            0>>

        """
        giac_console()


    def completions(self, s):
        """
        Return all commands that complete the command starting with the
        string s.

        EXAMPLES::

            sage: c = giac.completions('cas')  # optional - giac
            sage: 'cas_setup' in c             # optional - giac
            True
        """
        bs = chr(8)*len(s)
        if self._expect is None:
            self._start()
        E = self._expect
        E.sendline('%s%s%s'%(s,chr(63),chr(13)))
        t = E.timeout
        E.timeout=0.3  # since some things have no completion
        try:
            E.expect('----')
        except pexpect.TIMEOUT:
            E.timeout = t
            return []
        E.timeout = t
        v = E.before
        E.expect(self._prompt)
        E.expect(self._prompt)
        return v.split()[1:]

    def _commands(self):
        """
        Return list of all commands defined in Giac.

        EXAMPLES::

            sage: c = giac._commands() # optional - giac
            sage: len(c) > 100          # optional - giac
            True
            sage: 'Psi' in c          # optional - giac
            True
        """
        try:
            v = sum([self.completions(chr(65+n)) for n in range(26)], []) + \
                sum([self.completions(chr(97+n)) for n in range(26)], [])
        except RuntimeError:
            print "\n"*3
            print "*"*70
            print "WARNING: You do not have a working version of Giac installed!"
            print "*"*70
            v = []
        v.sort()
        return v

    def trait_names(self, verbose=True, use_disk_cache=True):
        """
        Returns a list of all the commands defined in Giac and optionally
        (per default) store them to disk.

        EXAMPLES::

            sage: c = giac.trait_names(use_disk_cache=False, verbose=False) # optional - giac
            sage: len(c) > 100  # optional - giac
            True
            sage: 'factors' in c  # optional - giac
            True
        """
        try:
            return self.__trait_names
        except AttributeError:
            import sage.misc.persist
            if use_disk_cache:
                try:
                    self.__trait_names = sage.misc.persist.load(COMMANDS_CACHE)
                    return self.__trait_names
                except IOError:
                    pass
            if verbose:
                print "\nBuilding Giac command completion list (this takes"
                print "a few seconds only the first time you do it)."
                print "To force rebuild later, delete %s."%COMMANDS_CACHE
            v = self._commands()
            self.__trait_names = v
            if len(v) > 200:
                # Giac is actually installed.
                sage.misc.persist.save(v, COMMANDS_CACHE)
            return v


    def cputime(self, t=None):
        r"""
        Returns the amount of CPU time that the Giac session has used. If
        ``t`` is not None, then it returns the difference
        between the current CPU time and ``t``.

        EXAMPLES::

            sage: t = giac.cputime() # optional - giac
            sage: t                   # random; optional - giac
            0.02
            sage: x = giac('x')      # optional - giac
            sage: giac.diff(x^2, x)  # optional - giac
            2*x
            sage: giac.cputime(t)    # random; optional - giac
            0.0
        """
        if t is None:
            return float(self('time()'))
        else:
            return float(self('time() - %s'%float(t)))


    def _eval_line(self, line, allow_use_file=True, wait_for_prompt=True, restart_if_needed=False):
        """
        EXAMPLES::

            sage: giac._eval_line('2+2')  # optional - giac
            '4'

            sage: A=matrix([range(280)])  # optional - giac
            sage: GA=giac(A)              # optional - giac
        """
        with gc_disabled():
            z = Expect._eval_line(self, line, allow_use_file=allow_use_file,
                    wait_for_prompt=wait_for_prompt)
            if z.lower().find("error") != -1:
                raise RuntimeError("An error occurred running a Giac command:\nINPUT:\n%s\nOUTPUT:\n%s"%(line, z))
        return z


    def eval(self, code, strip=True, **kwds):
        r"""
        Send the code x to the Giac interpreter.
        Remark: To enable multi-lines codes in the notebook magic mode: ``%giac``,
        the ``\n`` are removed before sending the code to giac.

        INPUT:

        - code -- str
        - strip -- Default is True and removes ``\n``

        EXAMPLES::

            sage: giac.eval("2+2;\n3") #optional - giac
            '4,3'
            sage: giac.eval("2+2;\n3",False) # optional - giac
            '4\n3'
            sage: s='g(x):={\nx+1;\nx+2;\n}'
            sage: giac(s)                    # optional - giac
            (x)->{
            x+1;
            x+2;
            }
            sage: giac.g(5)                   # optional - giac
            7
        """
        #we remove \n to enable multiline code in the notebook magic mode %giac
        if strip:
             code = code.replace("\n","").strip()
        ans = Expect.eval(self, code, strip=strip, **kwds).strip()
        return ans



    def set(self, var, value):
        """
        Set the variable var to the given value.

        EXAMPLES::

            sage: giac.set('xx', '2') # optional - giac
            sage: giac.get('xx')      # optional - giac
            '2'
        """
        cmd = '%s:=%s:;'%(var,value)   #if giac is not in maple mode ( maple_mode(0))
        out = self.eval(cmd)
        if out.find("error") != -1:
            raise TypeError("Error executing code in Giac\nCODE:\n\t%s\nGiac ERROR:\n\t%s"%(cmd, out))


    def get(self, var):
        """
        Get the value of the variable var.

        EXAMPLES::

            sage: giac.set('xx', '2') # optional - giac
            sage: giac.get('xx')      # optional - giac
            '2'
        """
        s = self.eval('%s'%var)
        return s

    def _object_class(self):
        """
        Returns the class of GiacElements.

        EXAMPLES::

            sage: giac._object_class()
            <class 'sage.interfaces.giac.GiacElement'>

        ::

            sage: m = giac(2)  # optional - giac
            sage: type(m)       # optional - giac
            <class 'sage.interfaces.giac.GiacElement'>
        """
        return GiacElement

    def _function_element_class(self):
        """
        Returns the GiacFunctionElement class.

        EXAMPLES::

            sage: giac._function_element_class()
            <class 'sage.interfaces.giac.GiacFunctionElement'>

        ::

            sage: two = giac(2)  # optional - giac
            sage: type(two.gcd)   # optional - giac
            <class 'sage.interfaces.giac.GiacFunctionElement'>
        """
        return GiacFunctionElement

    def _equality_symbol(self):
        """
        Returns the symbol used for equality testing in Giac.

        EXAMPLES::

            sage: giac._equality_symbol()               # optional - giac
            '=='

            sage: giac(2) == giac(2) # optional - giac
            True
        """
        return '=='

    def _true_symbol(self):
        """
        Returns the symbol used for truth in Giac.

        EXAMPLES::

            sage: giac._true_symbol()
            '1'

        ::

            sage: giac(2) == giac(2) # optional - giac
            True
        """
        return '1'

    def _assign_symbol(self):
        """
        Returns the symbol used for assignment in Giac.

        EXAMPLES::

            sage: giac._assign_symbol()
            ':='
        """
        return ":="

    def _help(self, str):
        r"""
        Returns the Giac help on ``str``.

        EXAMPLES::

            sage: giac._help('gcd')  # not tested ; output may vary (LANG)
            "...gcd - greatest common divisor of polynomials...
        """
        return os.popen('cas_help %s'%str).read()
        # return os.popen('echo "?%s" | giac'%str).read()

    def help(self, str):
        """
        Display Giac help about str. This is the same as typing "?str" in
        the Giac console.

        INPUT:

        -  ``str`` - a string to search for in the giac help
           system


        EXAMPLES::

            sage: giac.help('Psi')         # not tested - depends of giac and $LANG
            Psi(a,n)=nth-derivative of the function DiGamma (=ln@Gamma) at point a (Psi(a,0)=Psi(a))...
        """
        pager()(self._help(str))

    def clear(self, var):
        """
        Clear the variable named var.

        EXAMPLES::

            sage: giac.set('xx', '2')  # optional - giac
            sage: giac.get('xx')       # optional - giac
            '2'
            sage: giac.clear('xx')     # optional - giac
            sage: giac.get('xx')       # optional - giac
            'xx'
        """
        self.eval('purge(%s)'%var)

    def version(self):
        """
        Wrapper for giac's version().

        EXAMPLES::

            sage: giac.version()  # optional - giac
            "giac...

        """
        return giac('version()')

class GiacFunction(ExpectFunction):
    def _sage_doc_(self):
        """
        Returns the Giac help for this function. This gets called when
        doing "?" on self.

        EXAMPLES::

            sage: giac.gcd._sage_doc_()  # not tested ; output may vary LANG
            "gcd - greatest common divisor of polynomials...
        """
        M = self._parent
        return M._help(self._name)

class GiacFunctionElement(FunctionElement):
    def _sage_doc_(self):
        """
        Returns the Giac help for this function. This gets called when
        doing "?" on self.

        EXAMPLES::

            sage: two = giac(2)  # optional - giac
            sage: two.gcd._sage_doc_() # not tested; output may vary LANG
            "...gcd - greatest common divisor of polynomials...
        """
        return self._obj.parent()._help(self._name)

class GiacElement(ExpectElement):
    def __float__(self):
        """
        Returns a floating point version of self.

        EXAMPLES::

            sage: float(giac(1/2))   # optional - giac
            0.5
            sage: type(_)            # optional - giac
            <type 'float'>
        """
        M = self.parent()
        return float(giac.eval('evalf(%s)'%self.name()))


    def unapply(self,var):
        """
        Creates a Giac function in the given arguments from a Giac symbol.

        EXAMPLES::

            sage: f=giac('y^3+1+t')        # optional - giac
            sage: g=(f.unapply('y,t'))     # optional - giac
            sage: g                        # optional - giac
            (y,t)->y^3+1+t
            sage: g(1,2)                   # optional - giac
            4
        """
        return giac('unapply(%s,%s)'%(self,var))


    def __hash__(self):
        """
        Returns a  integer representing the hash of self.

        These examples are optional, and require Giac to be installed. You
        don't need to install any Sage packages for this.

        EXAMPLES::

            sage: m = giac('x^2+y^2')                       # optional - giac
            sage: hash(m)                                   # random; optional - giac
            4614285348919569149
        """
        return hash(giac.eval('string(%s);'%self.name()))


    def __cmp__(self, other):
        """
        Compare equality between self and other, using giac.

        These examples are optional, and require Giac to be installed. You
        don't need to install any Sage packages for this.

        EXAMPLES::

            sage: a = giac(5)                              # optional - giac
            sage: b = giac(5)                              # optional - giac
            sage: a == b                                   # optional - giac
            True
            sage: a == 5                                   # optional - giac
            True

        ::

            sage: c = giac(3)                              # optional - giac
            sage: a == c                                   # optional - giac
            False
            sage: a < c                                    # optional - giac
            False
            sage: a < 6                                    # optional - giac
            True
            sage: c <= a                                   # optional - giac
            True

        ::

        TESTS::

            sage: x = var('x')                            # optional - giac
            sage: t = giac((x+1)^2)                       # optional - giac
            sage: u = giac(x^2+2*x+1)                     # optional - giac
            sage: u == t                                  # optional - giac
            False
        """
        P = self.parent()
        if P.eval("evalb(%s %s %s)"%(self.name(), P._equality_symbol(),
                                 other.name())) == P._true_symbol():
            return 0
        # (to be tested with giac). Maple  does not allow comparing objects of different types and
        # it raises an error in this case.
        # We catch the error, and return True for <
        try:
            if P.eval("evalb(%s %s %s)"%(self.name(), P._lessthan_symbol(), other.name())) == P._true_symbol():
                return -1
        except RuntimeError as e:
            msg = str(e)
            if 'is not valid' in msg and 'to < or <=' in msg:
                if (hash(str(self)) < hash(str(other))):
                    return -1
                else:
                    return 1
            else:
                raise RuntimeError(e)
        if P.eval("evalb(%s %s %s)"%(self.name(), P._greaterthan_symbol(), other.name())) == P._true_symbol():
            return 1
        # everything is supposed to be comparable in Python, so we define
        # the comparison thus when no comparable in interfaced system.
        if (hash(self) < hash(other)):
            return -1
        else:
            return 1

    def trait_names(self):
        """
        EXAMPLES::

            sage: a = giac(2) # optional - giac
            sage: 'sin' in a.trait_names() # optional - giac
            True
        """
        return self.parent().trait_names()


    def __len__(self):
        """
        EXAMPLES::

            sage: len(giac([1,2,3]))                # optional - giac
            3
        """
        return int(self.size())

    def __iter__(self):
        """
        EXAMPLES:
            sage: l = giac([1,2,3]) #optional
            sage: list(iter(l))          #optional
            [1, 2, 3]
        """
        for i in range(len(self)):  # zero-indexed if giac is maple_mode(0)
            yield self[i]

    def __del__(self):
        """
        Note that clearing object is pointless since it wastes time.
        (Ex: otherwise doing a=0 after a = (giac('x+y+z')^40).normal() is very slow )

        EXAMPLES::

            sage: a = giac(2)                        # optional - giac
            sage: a.__del__()                        # optional - giac
            sage: a                                  # optional - giac
            2
            sage: del a                              # optional - giac
            sage: a
            Traceback (most recent call last):
            ...
            NameError: name 'a' is not defined
        """
        return

    def __repr__(self):
        """
        Return a string representation of self.

        These examples are optional, and require Giac to be installed. You
        don't need to install any Sage packages for this.

        EXAMPLES::

            sage: x = var('x')
            sage: giac(x)                      # optional - giac
            x
            sage: giac(5)                      # optional - giac
            5
            sage: M = matrix(QQ,2,range(4))
            sage: giac(M)                      # optional - giac
            [[0,1],[2,3]]
        """
        self._check_valid()
        return self.parent().get(self._name)

    def _latex_(self):
        r"""
        You can output Giac expressions in latex.

        EXAMPLES::

            sage: print latex(giac('(x^4 - y)/(y^2-3*x)'))      # optional - giac
            "\frac{(x^{4}-y)}{(y^{2}-3 x)}"

        """
        return self.parent().eval('latex(%s)'%self.name())


    def _matrix_(self, R):
        r"""
        Return matrix over the (Sage) ring R determined by self, where self
        should be a  Giac matrix.
        Warning: It is slow, don't convert big matrices.

        EXAMPLES::

            sage: R.<x,y>=QQ[]
            sage: M=giac('matrix(4,4,(k,l)->(x^k-y^l))'); M      # optional - giac
            matrix[[0,1-y,1-y^2,1-y^3],[x-1,x-y,x-y^2,x-y^3],[x^2-1,x^2-y,x^2-y^2,x^2-y^3],[x^3-1,x^3-y,x^3-y^2,x^3-y^3]]
            sage: M.eigenvals()       # random; optional - giac
            0,0,(x^3+x^2+x-y^3-y^2-y+sqrt(x^6+2*x^5+3*x^4-14*x^3*y^3+2*x^3*y^2+2*x^3*y+6*x^3+2*x^2*y^3-14*x^2*y^2+2*x^2*y+5*x^2+2*x*y^3+2*x*y^2-14*x*y+4*x+y^6+2*y^5+3*y^4+6*y^3+5*y^2+4*y-12))/2,(x^3+x^2+x-y^3-y^2-y-sqrt(x^6+2*x^5+3*x^4-14*x^3*y^3+2*x^3*y^2+2*x^3*y+6*x^3+2*x^2*y^3-14*x^2*y^2+2*x^2*y+5*x^2+2*x*y^3+2*x*y^2-14*x*y+4*x+y^6+2*y^5+3*y^4+6*y^3+5*y^2+4*y-12))/2
            sage: Z=matrix(M,R);Z                                # optional - giac
            [         0     -y + 1   -y^2 + 1   -y^3 + 1]
            [     x - 1      x - y   -y^2 + x   -y^3 + x]
            [   x^2 - 1    x^2 - y  x^2 - y^2 -y^3 + x^2]
            [   x^3 - 1    x^3 - y  x^3 - y^2  x^3 - y^3]
            sage: parent(Z)                                      # optional - giac
            Full MatrixSpace of 4 by 4 dense matrices over Multivariate Polynomial Ring in x, y over Rational Field
        """
        P = self.parent()
        v = self.dim()
        n = int(v[0])
        m = int(v[1])

        from sage.matrix.matrix_space import MatrixSpace
        M = MatrixSpace(R, n, m)
        entries = [[R(self[r,c]) for c in range(m)] for r in range(n)]
        return M(entries)


    def _sage_(self):
        r"""
        Convert a giac expression back to a Sage expression.

        This currently does not implement a parser for the Giac output language,
        therefore only very simple expressions will convert successfully.
        Warning: List conversion is slow.

        EXAMPLE::

        sage: m = giac('x^2 + 5*y')                            # optional - giac
        sage: m.sage()                                          # optional - giac
        x^2 + 5*y

        ::

        sage: m = giac('sin(2*sqrt(1-x^2)) * (1 - cos(1/x))^2')  # optional - giac
        sage: m.trigexpand().sage()                              # optional - giac
        2*(cos(1/x) - 1)^2*sin(sqrt(-x^2 + 1))*cos(sqrt(-x^2 + 1))

        """
        result = repr(self)
        if str(self.type()) != 'DOM_LIST' :
            try:
                from sage.symbolic.all import SR
                return SR(result)
            except Exception:
                raise NotImplementedError("Unable to parse Giac output: %s" % result)
        else:
            return [entry.sage() for entry in self]

    def integral(self, var='x', min=None, max=None):
        r"""
        Return the integral of self with respect to the variable x.

        INPUT:


        -  ``var`` - variable

        -  ``min`` - default: None

        -  ``max`` - default: None


        Returns the definite integral if xmin is not None, otherwise
        returns an indefinite integral.

        EXAMPLES::

            sage: y=giac('y');f=(sin(2*y)/y).integral(y).simplify(); f        # optional - giac
            Si(2*y)
            sage: f.diff(y).simplify()                                        # optional - giac
            sin(2*y)/y

        ::

            sage: f = giac('exp(x^2)').integral('x',0,1) ; f                # optional - giac
            integra...
            sage: f.evalf(100)                                              # optional - giac
            1.4626517459071819025155758073473096674669301007326185820691973210905694258465619632003390265815626744
            sage: x,y=giac('x'),giac('y');integrate(cos(x+y),'x=0..pi').simplify()     # optional - giac
            -2*sin(y)
        """
        if min is None:
            return giac('int(%s,%s)'%(self.name(),var))
        else:
            if max is None:
                raise ValueError("neither or both of min/max must be specified.")
        return giac('int(%s,%s,%s,%s)'%(self.name(),var,giac(min),giac(max)))

    integrate=integral


    def sum(self, var, min=None, max=None):
        r"""
        Return the sum of self with respect to the variable x.

        INPUT:


        -  ``var`` - variable

        -  ``min`` - default: None

        -  ``max`` - default: None


        Returns the definite integral if xmin is not None, otherwise
        returns an indefinite integral.

        EXAMPLES::
        sage: giac('1/(1+k^2)').sum('k',-oo,+infinity).simplify()     # optional -  giac
        (pi*exp(pi)^2+pi)/(exp(pi)^2-1)

        """
        if min is None:
            return giac('sum(%s,%s)'%(self.name(),var))
        else:
            if max is None:
                raise ValueError("neither or both of min/max must be specified.")
            return giac('sum(%s,%s,%s,%s)'%(self.name(),var,giac(min),giac(max)))


# An instance
giac = Giac()

def reduce_load_Giac():
    """
    Returns the giac object created in sage.interfaces.giac.

    EXAMPLES::

        sage: from sage.interfaces.giac import reduce_load_Giac
        sage: reduce_load_Giac()
        Giac
    """
    return giac


import os
def giac_console():
    """
    Spawn a new Giac command-line session.

    EXAMPLES::

        sage: giac.console()  # not tested - giac
        ...
        Homepage http://www-fourier.ujf-grenoble.fr/~parisse/giac.html
        Released under the GPL license 3.0 or above
        See http://www.gnu.org for license details
        -------------------------------------------------
        Press CTRL and D simultaneously to finish session
        Type ?commandname for help
        """
    os.system('giac')


def __doctest_cleanup():
    """
    EXAMPLES::

        sage: from sage.interfaces.giac import __doctest_cleanup
        sage: m = giac(2)         # optional - giac
        sage: giac.is_running()   # optional - giac
        True
        sage: __doctest_cleanup()
        sage: giac.is_running()
        False
    """
    import sage.interfaces.quit
    sage.interfaces.quit.expect_quitall()



