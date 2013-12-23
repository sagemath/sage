r"""
Interface to R

The following examples try to follow "An Introduction to R" which can
be found at http://cran.r-project.org/doc/manuals/R-intro.html .

EXAMPLES:

Simple manipulations; numbers and vectors

The simplest data structure in R is the numeric vector which
consists of an ordered collection of numbers.  To create a
vector named $x$ using the R interface in Sage, you pass the
R interpreter object a list or tuple of numbers::

    sage: x = r([10.4,5.6,3.1,6.4,21.7]); x
    [1] 10.4  5.6  3.1  6.4 21.7

You can invert elements of a vector x in R by using the
invert operator or by doing 1/x::

    sage: ~x
    [1] 0.09615385 0.17857143 0.32258065 0.15625000 0.04608295
    sage: 1/x
    [1] 0.09615385 0.17857143 0.32258065 0.15625000 0.04608295

The following assignment creates a vector $y$ with 11 entries which
consists of two copies of $x$ with a 0 in between::

    sage: y = r([x,0,x]); y
    [1] 10.4  5.6  3.1  6.4 21.7  0.0 10.4  5.6  3.1  6.4 21.7

Vector Arithmetic

The following command generates a new vector $v$ of length 11 constructed
by adding together (element by element) $2x$ repeated 2.2 times, $y$
repeated just once, and 1 repeated 11 times::

    sage: v = 2*x+y+1; v
    [1] 32.2 17.8 10.3 20.2 66.1 21.8 22.6 12.8 16.9 50.8 43.5

One can compute the sum of the elements of an R vector in the following
two ways::

    sage: sum(x)
    [1] 47.2
    sage: x.sum()
    [1] 47.2

One can calculate the sample variance of a list of numbers::

    sage: ((x-x.mean())^2/(x.length()-1)).sum()
    [1] 53.853
    sage: x.var()
    [1] 53.853

    sage: x.sort()
    [1] 3.1  5.6  6.4 10.4 21.7
    sage: x.min()
    [1] 3.1
    sage: x.max()
    [1] 21.7
    sage: x
    [1] 10.4  5.6  3.1  6.4 21.7

    sage: r(-17).sqrt()
    [1] NaN
    sage: r('-17+0i').sqrt()
    [1] 0+4.123106i

Generating an arithmetic sequence::

    sage: r('1:10')
    [1] 1  2  3  4  5  6  7  8  9 10

Because ``from`` is a keyword in Python, it can't be used
as a keyword argument.  Instead, ``from_`` can be passed, and
R will recognize it as the correct thing::

    sage: r.seq(length=10, from_=-1, by=.2)
    [1] -1.0 -0.8 -0.6 -0.4 -0.2  0.0  0.2  0.4  0.6  0.8

    sage: x = r([10.4,5.6,3.1,6.4,21.7]);
    sage: x.rep(2)
    [1] 10.4  5.6  3.1  6.4 21.7 10.4  5.6  3.1  6.4 21.7
    sage: x.rep(times=2)
    [1] 10.4  5.6  3.1  6.4 21.7 10.4  5.6  3.1  6.4 21.7
    sage: x.rep(each=2)
    [1] 10.4 10.4  5.6  5.6  3.1  3.1  6.4  6.4 21.7 21.7

Missing Values::

    sage: na = r('NA')
    sage: z = r([1,2,3,na])
    sage: z
    [1]  1  2  3 NA
    sage: ind = r.is_na(z)
    sage: ind
    [1] FALSE FALSE FALSE  TRUE
    sage: zero = r(0)
    sage: zero / zero
    [1] NaN
    sage: inf = r('Inf')
    sage: inf-inf
    [1] NaN
    sage: r.is_na(inf)
    [1] FALSE
    sage: r.is_na(inf-inf)
    [1] TRUE
    sage: r.is_na(zero/zero)
    [1] TRUE
    sage: r.is_na(na)
    [1] TRUE
    sage: r.is_nan(inf-inf)
    [1] TRUE
    sage: r.is_nan(zero/zero)
    [1] TRUE
    sage: r.is_nan(na)
    [1] FALSE


Character Vectors::

    sage: labs = r.paste('c("X","Y")', '1:10', sep='""'); labs
    [1] "X1"  "Y2"  "X3"  "Y4"  "X5"  "Y6"  "X7"  "Y8"  "X9"  "Y10"


Index vectors; selecting and modifying subsets of a data set::

    sage: na = r('NA')
    sage: x = r([10.4,5.6,3.1,6.4,21.7,na]); x
    [1] 10.4  5.6  3.1  6.4 21.7   NA
    sage: x['!is.na(self)']
    [1] 10.4  5.6  3.1  6.4 21.7

    sage: x = r([10.4,5.6,3.1,6.4,21.7,na]); x
    [1] 10.4  5.6  3.1  6.4 21.7   NA
    sage: (x+1)['(!is.na(self)) & self>0']
    [1] 11.4  6.6  4.1  7.4 22.7
    sage: x = r([10.4,-2,3.1,-0.5,21.7,na]); x
    [1] 10.4 -2.0  3.1 -0.5 21.7   NA
    sage: (x+1)['(!is.na(self)) & self>0']
    [1] 11.4  4.1  0.5 22.7

Distributions::

    sage: r.options(width="60");
    $width
    [1] 100

    sage: rr = r.dnorm(r.seq(-3,3,0.1))
    sage: rr
     [1] 0.004431848 0.005952532 0.007915452 0.010420935
     [5] 0.013582969 0.017528300 0.022394530 0.028327038
     [9] 0.035474593 0.043983596 0.053990967 0.065615815
    [13] 0.078950158 0.094049077 0.110920835 0.129517596
    [17] 0.149727466 0.171368592 0.194186055 0.217852177
    [21] 0.241970725 0.266085250 0.289691553 0.312253933
    [25] 0.333224603 0.352065327 0.368270140 0.381387815
    [29] 0.391042694 0.396952547 0.398942280 0.396952547
    [33] 0.391042694 0.381387815 0.368270140 0.352065327
    [37] 0.333224603 0.312253933 0.289691553 0.266085250
    [41] 0.241970725 0.217852177 0.194186055 0.171368592
    [45] 0.149727466 0.129517596 0.110920835 0.094049077
    [49] 0.078950158 0.065615815 0.053990967 0.043983596
    [53] 0.035474593 0.028327038 0.022394530 0.017528300
    [57] 0.013582969 0.010420935 0.007915452 0.005952532
    [61] 0.004431848

Convert R Data Structures to Python/Sage::

    sage: rr = r.dnorm(r.seq(-3,3,0.1))
    sage: sum(rr._sage_())
    9.9772125168981...

Or you get a dictionary to be able to access all the information::

    sage: rs = r.summary(r.c(1,4,3,4,3,2,5,1))
    sage: rs
       Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
      1.000   1.750   3.000   2.875   4.000   5.000
      sage: d = rs._sage_()
      sage: d['DATA']
      [1, 1.75, 3, 2.875, 4, 5]
      sage: d['_Names']
      ['Min.', '1st Qu.', 'Median', 'Mean', '3rd Qu.', 'Max.']
      sage: d['_r_class']
      ['summaryDefault', 'table']

It is also possible to access the plotting capabilities of R
through Sage.  For more information see the documentation of
r.plot() or r.png().

AUTHORS:

- Mike Hansen (2007-11-01)
- William Stein (2008-04-19)
- Harald Schilly (2008-03-20)
- Mike Hansen (2008-04-19)
"""

##########################################################################
#
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#                     2007 Mike Hansen   <mhansen@gmail.com>
#                     2008 Harald Schilly <harald.schilly@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#
##########################################################################

from expect import Expect, ExpectElement, ExpectFunction, FunctionElement
from sage.misc.misc import DOT_SAGE
import re
import sage.rings.integer

COMMANDS_CACHE = '%s/r_commandlist.sobj'%DOT_SAGE
PROMPT = '__SAGE__R__PROMPT__> '
prompt_re = re.compile("^>", re.M)

#there is a mirror network, but lets take #1 for now
RRepositoryURL = "http://cran.r-project.org/"
RFilteredPackages = ['.GlobalEnv']

# crosscheck with https://svn.r-project.org/R/trunk/src/main/names.c
# but package:base should cover this. i think.
RBaseCommands = ['c', "NULL", "NA", "True", "False", "Inf", "NaN"]

class R(Expect):
    def __init__(self,
                 maxread=100000, script_subdirectory=None,
                 server_tmpdir = None,
                 logfile=None,
                 server=None,
                 init_list_length=1024):
        """
        An interface to the R interpreter.

        R is a comprehensive collection of methods for statistics,
        modelling, bioinformatics, data analysis and much more.
        For more details, see http://www.r-project.org/about.html

        Resources:

        * http://r-project.org/ provides more information about R.
        * http://rseek.org/ R's own search engine.

        EXAMPLES::

             sage: r.summary(r.c(1,2,3,111,2,3,2,3,2,5,4))
             Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
             1.00    2.00    3.00   12.55    3.50  111.00

        TESTS::

            sage: r == loads(dumps(r))
            True
        """
        Expect.__init__(self,

                  # The capitalized version of this is used for printing.
                  name = 'r',

                  # This is regexp of the input prompt.  If you can change
                  # it to be very obfuscated that would be better.   Even
                  # better is to use sequence numbers.
                  # options(prompt=\"<prompt> \")
                  prompt = '> ', #default, later comes the change

                  # This is the command that starts up your program
                  command = "R --vanilla --quiet",

                  maxread = maxread,

                  server=server,
                  server_tmpdir=server_tmpdir,

                  script_subdirectory = script_subdirectory,

                  # If this is true, then whenever the user presses Control-C to
                  # interrupt a calculation, the whole interface is restarted.
                  restart_on_ctrlc = False,

                  # If true, print out a message when starting
                  # up the command when you first send a command
                  # to this interface.
                  verbose_start = False,

                  logfile=logfile,

                  # If an input is longer than this number of characters, then
                  # try to switch to outputting to a file.
                  eval_using_file_cutoff=1024)

        self.__seq = 0
        self.__var_store_len = 0
        self.__init_list_length = init_list_length
        self._prompt_wait = [self._prompt]

    def _start(self):
        """
        Start up the R interpreter and sets the initial prompt and options.

        This is called the first time the R interface is actually used.

        EXAMPLES::

            sage: r = R()
            sage: r._start()
        """
        Expect._start(self)

        # width is line width, what's a good value? maximum is 10000!
        # pager needed to replace help view from less to printout
        # option device= is for plotting, is set to x11, NULL would be better?
        self._change_prompt(PROMPT)
        self.eval('options(prompt=\"%s\",continue=\"%s\", width=100,pager="cat",device="png")'%(PROMPT, PROMPT))
        self.expect().expect(PROMPT)
        self.eval('options(repos="%s")'%RRepositoryURL)
        self.eval('options(CRAN="%s")'%RRepositoryURL)

        # don't abort on errors, just raise them!
        # necessary for non-interactive execution
        self.eval('options(error = expression(NULL))')

    def png(self, *args, **kwds):
        """
        Creates an R PNG device.

        This should primarily be used to save an R graphic to a custom file.  Note
        that when using this in the notebook, one must plot in the same cell that
        one creates the device.  See r.plot() documentation for more information
        about plotting via R in Sage.

        These examples won't work on the many platforms where R still gets
        built without graphics support.

        EXAMPLES::

            sage: filename = tmp_filename() + '.png'
            sage: r.png(filename='"%s"'%filename)             # optional -- rgraphics
            NULL
            sage: x = r([1,2,3])
            sage: y = r([4,5,6])
            sage: r.plot(x,y) # This saves to filename, but is not viewable from command line; optional -- rgraphics
            null device
                      1
            sage: import os; os.unlink(filename) # We remove the file for doctesting; optional -- rgraphics

        We want to make sure that we actually can view R graphics, which happens
        differently on different platforms::

            sage: s = r.eval('capabilities("png")') # Should be on Linux and Solaris
            sage: t = r.eval('capabilities("aqua")') # Should be on all supported Mac versions
            sage: "TRUE" in s+t                      # optional -- rgraphics
            True
        """
        #Check to see if R has PNG support
        s = self.eval('capabilities("png")')
        t = r.eval('capabilities("aqua")')
        if "TRUE" not in s+t:
            raise RuntimeError, "R was not compiled with PNG support"

        from sage.server.support import EMBEDDED_MODE
        if EMBEDDED_MODE:
            self.setwd('"%s"'%os.path.abspath('.'))
        return RFunction(self, 'png')(*args, **kwds)

    def convert_r_list(self, l):
        r"""
        Converts an R list to a Python list.

        EXAMPLES::

            sage: s = 'c(".GlobalEnv", "package:stats", "package:graphics", "package:grDevices", \n"package:utils", "package:datasets", "package:methods", "Autoloads", \n"package:base")'
            sage: r.convert_r_list(s)
            ['.GlobalEnv',
             'package:stats',
             'package:graphics',
             'package:grDevices',
             'package:utils',
             'package:datasets',
             'package:methods',
             'Autoloads',
             'package:base']
        """
        pl = []
        l = "".join(l.split("\n"))
        l = l[2:len(l)-1]
        for l in l.split(","):
            pl.append(l.strip().strip('"'))
        return pl

    def install_packages(self, package_name):
        """
        Install an R package into Sage's R installation.

        EXAMPLES::

            sage: r.install_packages('aaMI')       # not tested
            ...
            R is free software and comes with ABSOLUTELY NO WARRANTY.
            You are welcome to redistribute it under certain conditions.
            Type 'license()' or 'licence()' for distribution details.
            ...
            Please restart Sage in order to use 'aaMI'.
        """
        cmd = """options(repos="%s"); install.packages("%s")"""%(RRepositoryURL, package_name)
        os.system("time echo '%s' | R --vanilla"%cmd)
        print "Please restart Sage in order to use '%s'."%package_name

    def __repr__(self):
        """
        Return string representation of this R interface.

        EXAMPLES::

            sage: r.__repr__()
            'R Interpreter'
        """
        return 'R Interpreter'

    def __reduce__(self):
        """
        Used in serializing an R interface.

        EXAMPLES::

            sage: rlr, t = r.__reduce__()
            sage: rlr(*t)
            R Interpreter
        """
        return reduce_load_R, tuple([])

    def __getattr__(self, attrname):
        """
        Called when you get an attribute of the R interface.  This
        manufactures an R function, which is a Python function that
        can then be called with various inputs.

        EXAMPLES::

            sage: c = r.c; c
            c
            sage: type(c)
            <class 'sage.interfaces.r.RFunction'>
        """
        if attrname[:1] == "_":
            raise AttributeError
        return RFunction(self, attrname)


    def _quit_string(self):
        r"""
        Return the string that when typed into R causes the R
        interpreter to exit.

        EXAMPLES::

            sage: r._quit_string()
            'quit(save="no")'
        """
        return 'quit(save="no")'

    def _read_in_file_command(self, filename):
        r"""
        Return the R command (as a string) to read in a file named
        filename into the R interpreter.

        EXAMPLES::

            sage: r._read_in_file_command('file.txt')
            'file=file("file.txt",open="r")\nsource(file)'
        """
        return 'file=file("%s",open="r")\nsource(file)'%filename

    def read(self, filename):
        r"""
        Read filename into the R interpreter by calling R's source function on a
        read-only file connection.

        EXAMPLES::

            sage: filename = tmp_filename()
            sage: f = open(filename, 'w')
            sage: f.write('a <- 2+2\n')
            sage: f.close()
            sage: r.read(filename)
            sage: r.get('a')
            '[1] 4'
        """
        self.eval( self._read_in_file_command(filename) )

    def _install_hints(self):
        """
        EXAMPLES::

            sage: print r._install_hints()
            R is currently installed with Sage.
        """
        return "R is currently installed with Sage.\n"

    def _source(self, s):
        """
        Returns the source code of an R function as a string.

        INPUT:

        - s -- the name of the function as a string

        EXAMPLES::

            sage: print r._source("print.anova")
            function (x, digits = max(getOption("digits") - 2L, 3L), signif.stars = getOption("show.signif.stars"),
            ...
        """
        if s[-2:] == "()":
            s = s[-2:]
        return self.eval('%s'%s)

    def source(self, s):
        """
        Display the R source (if possible) about the function named s.

        INPUT:

        - s -- a string representing the function whose source code you want to see

        OUTPUT: string -- source code

        EXAMPLES::

            sage: print r.source("print.anova")
            function (x, digits = max(getOption("digits") - 2L, 3L), signif.stars = getOption("show.signif.stars"),
            ...
        """
        return self._source(s)

    def version(self):
        """
        Return the version of R currently running.

        OUTPUT: tuple of ints; string

        EXAMPLES::

            sage: r.version() # not tested
            ((3, 0, 1), 'R version 3.0.1 (2013-05-16)')
            sage: rint, rstr = r.version()
            sage: rint[0] >= 3
            True
            sage: rstr.startswith('R version')
            True
        """
        major_re = re.compile('^major\s*(\d.*?)$', re.M)
        minor_re = re.compile('^minor\s*(\d.*?)$', re.M)
        version_string_re = re.compile('^version.string\s*(R.*?)$', re.M)

        s = self.eval('version')

        major = int( major_re.findall(s)[0].strip() )
        minor = tuple(int(i) for i in minor_re.findall(s)[0].strip().split(".") )
        version_string = version_string_re.findall(s)[0].strip()

        return ( (major,) + minor, version_string )

    def library(self, library_name):
        """
        Load the library library_name into the R interpreter.

        This function raises an ImportError if the given library
        is not known.

        INPUT:

        - library_name -- string

        EXAMPLES::

            sage: r.library('grid')
            sage: 'grid' in r.eval('(.packages())')
            True
            sage: r.library('foobar')
            Traceback (most recent call last):
            ...
            ImportError: ...
        """
        ret = self.eval('require("%s")'%library_name)
        # try hard to parse the message string in a locale-independent way
        if ' library(' in ret:       # locale-independent key-word
            raise ImportError, "%s"%ret
        else:
            try:
                # We need to rebuild keywords!
                del self.__trait_names
            except AttributeError:
                pass
            self.trait_names(verbose=False, use_disk_cache=False)

    require = library #overwrites require

    def available_packages(self):
        """
        Returns a list of all available R package names.

        This list is not necessarily sorted.

        OUTPUT: list of strings

        .. note::

            This requires an internet connection. The CRAN server is
            that is checked is defined at the top of sage/interfaces/r.py.

        EXAMPLES::

            sage: ap = r.available_packages()   # optional - internet
            sage: len(ap) > 20                  #optional
            True
        """
        p = self.new('available.packages("%s/src/contrib")'%RRepositoryURL)
        s = str(p).splitlines()[1:]
        v = [x.split()[0].strip("'") for x in s]
        return v
        #The following was more structural, but breaks on my machine.  (stein)
        #p = p._sage_()
        #s = p['_Dim'][0]
        #l = [[p['DATA'][i],p['DATA'][s+1+i]] for i in xrange(0,s)]
        #return l

    def _object_class(self):
        """
        Return the underlying class of elements of the R interface object.

        OUTPUT: a Python class

        EXAMPLES::

            sage: r._object_class()
            <class 'sage.interfaces.r.RElement'>
        """
        return RElement

    def _true_symbol(self):
        """
        Return the symbol that represents True in R.

        OUTPUT: string

        EXAMPLES::

            sage: r._true_symbol()
            '[1] TRUE'

        This is used behinds the scenes to implement comparison::

            sage: r('1') == r('1')
            [1] TRUE
            sage: r('1') == r('2')
            [1] FALSE
        """
        # return the string rep of truth, i.e., what the system outputs
        # when you type 1==1.
        return "[1] TRUE"

    def _false_symbol(self):
        """
        Return the symbol that represents True in R.

        OUTPUT: string

        EXAMPLES::

            sage: r._false_symbol()
            '[1] FALSE'
        """
        # return the string rep of false, i.e., what the system outputs
        # when you type 1==2.
        return "[1] FALSE"

    def _equality_symbol(self):
        """
        EXAMPLES::

            sage: r._equality_symbol()
            '=='
        """
        # return the symbol for checking equality, e.g., == or eq.
        return "=="

    def help(self, command):
        """
        Returns help string for a given command.

        INPUT:
        - command -- a string

        OUTPUT: HelpExpression -- a subclass of string whose __repr__ method is __str__, so it prints nicely

        EXAMPLES::

            sage: r.help('print.anova')
            anova                 package:stats                 R Documentation
            ...
                 Chambers, J. M. and Hastie, T. J. (1992) _Statistical Models in
                 S_, Wadsworth & Brooks/Cole.
            ...

        .. note::

            This is similar to typing r.command?.
        """
        s = self.eval('help("%s")'%command).strip()     # ?cmd is only an unsafe shortcut
        import sage.plot.plot
        if sage.plot.plot.EMBEDDED_MODE:
            s = s.replace('_\x08','')
        return HelpExpression(s)

    def _assign_symbol(self):
        """
        Return the symbol used in R for assignment, which is ' <- '.

        OUTPUT: string

        EXAMPLES::

            sage: r._assign_symbol()
            ' <- '
        """
        return " <- "

    def _left_list_delim(self):
        """
        Return the left delimiter for lists in R, which is 'c('

        OUTPUT: string

        EXAMPLES::

            sage: r._left_list_delim()
            'c('
        """
        return "c("

    def _right_list_delim(self):
        """
        Return the right delimiter for lists in R, which is 'c('

        OUTPUT: string

        EXAMPLES::

            sage: r._right_list_delim()
            ')'
        """
        return ")"

    def console(self):
        """
        Runs the R console as a separate new R process.

        EXAMPLES::

            sage: r.console()                    # not tested
                R version 2.6.1 (2007-11-26)
                Copyright (C) 2007 The R Foundation for Statistical Computing
                ISBN 3-900051-07-0
                ...
        """
        r_console()

    def function_call(self, function, args=None, kwds=None):
        """
        Return the result of calling an R function, with given args and keyword args.

        OUTPUT: RElement -- an object in R

        EXAMPLES::

            sage: r.function_call('length', args=[ [1,2,3] ])
            [1] 3
        """
        args, kwds = self._convert_args_kwds(args, kwds)
        self._check_valid_function_name(function)
        return self.new("%s(%s)"%(function, ",".join([s.name() for s in args] +
                                                     [self._sage_to_r_name(key)+'='+kwds[key].name() for key in kwds ] )))

    def call(self, function_name, *args, **kwds):
        r"""
        This is an alias for :meth:`function_call`.

        EXAMPLES::

            sage: r.call('length', [1,2,3])
            [1] 3
        """
        return self.function_call(function_name, args=args, kwds=kwds)

    def _an_element_impl(self):
        """
        Returns an element belonging to the R interpreter.  This is used
        behind the scenes when doing things like comparisons, etc.

        OUTPUT: RElement -- an R element.

        EXAMPLES::

            sage: r._an_element_impl()
            [1] 0
            sage: type(_)
            <class 'sage.interfaces.r.RElement'>
        """
        return self(0)

    def set(self, var, value):
        """
        Set the variable var in R to what the string value evaluates to in R.

        INPUT:

        - var -- a string
        - value -- a string

        EXAMPLES::

            sage: r.set('a', '2 + 3')
            sage: r.get('a')
            '[1] 5'
        """
        cmd = '%s <- %s'%(var,value)
        out = self.eval(cmd)
        if out.find("error") != -1:
            raise TypeError, "Error executing code in R\nCODE:\n\t%s\nR ERROR:\n\t%s"%(cmd, out)

    def get(self, var):
        """
        Returns the string representation of the variable var.

        INPUT:

        - var -- a string

        OUTPUT: string

        EXAMPLES::

            sage: r.set('a', 2)
            sage: r.get('a')
            '[1] 2'
        """
        s = self.eval('%s'%var)
        #return self._remove_indices_re.sub("", s).strip()
        return s

    def na(self):
        """
        Returns the NA in R.

        OUTPUT: RElement -- an element of R

        EXAMPLES::

            sage: r.na()
            [1] NA
        """
        return self('NA')

    def completions(self, s):
        """
        Return all commands names that complete the command starting with the
        string s.   This is like typing s[Ctrl-T] in the R interpreter.

        INPUT:

        - s -- string

        OUTPUT: list -- a list of strings

        EXAMPLES::

            sage: dummy = r.trait_names(use_disk_cache=False)    #clean doctest
            sage: r.completions('tes')
            ['testInheritedMethods', 'testPlatformEquivalence', 'testVirtual']
        """
        return [name for name in self.trait_names() if name[:len(s)] == s]

    def _commands(self):
        """
        Return list of all commands defined in R.

        OUTPUT: list -- a sorted list of strings

        EXAMPLES::

            sage: l = r._commands()
            sage: 'AIC' in l
            True
            sage: len(l) > 200
            True
        """
        v = RBaseCommands

        ll = self.eval('dput(search())') # loaded libs
        ll = self.convert_r_list(ll)

        for lib in ll:
            if lib in RFilteredPackages:
                continue

            if lib.find("package:") != 0:
                continue #only packages

            raw = self.eval('dput(objects("%s"))'%lib)
            raw = self.convert_r_list(raw)
            raw = [x.replace(".","_") for x in raw]

            #TODO are there others? many of them are shortcuts or
            #should be done on another level, like selections in lists
            #instead of calling obj.[[( fun-args) or other crazy stuff like that

            #TODO further filtering, check if strings are now
            #really functions, in R: exists(s, mode = "function"))
            # (apply to vector with sapply(vec,func))

            #filter only python compatible identifiers
            valid = re.compile('[^a-zA-Z0-9_]+')
            raw = [x for x in raw if valid.search(x) is None]
            v += raw

        v.sort()
        return v

    def trait_names(self, verbose=True, use_disk_cache=True):
        """
        Return list of all R functions.

        INPUT:

        - verbose -- bool (default: True); if True, display debugging information
        - use_disk_cache -- bool (default: True); if True, use the disk cache of
          trait names to save time.

        OUTPUT: list -- list of string

        EXAMPLES::

            sage: t = r.trait_names(verbose=False)
            sage: len(t) > 200
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
            if verbose and use_disk_cache:
                print "\nBuilding R command completion list (this takes"
                print "a few seconds only the first time you do it)."
                print "To force rebuild later, delete %s."%COMMANDS_CACHE
            v = self._commands()
            self.__trait_names = v
            if len(v) > 200 and use_disk_cache:
                sage.misc.persist.save(v, COMMANDS_CACHE)
            return v

    def plot(self, *args, **kwds):
        """
        The R plot function.  Type r.help('plot') for much more extensive
        documentation about this function.  See also below for a brief
        introduction to more plotting with R.

        If one simply wants to view an R graphic, using this function is
        is sufficient (because it calls dev.off() to turn off the device).

        However, if one wants to save the graphic to a specific file, it
        should be used as in the example below to write the output.

        EXAMPLES:

        This example saves a plot to the standard R output, usually
        a filename like ``Rplot001.png`` - from the command line, in
        the current directory, and in the cell directory in the notebook::

            sage: d=r.setwd('"%s"'%SAGE_TMP)    # for doctesting only; ignore if you are trying this;
            sage: r.plot("1:10")                # optional -- rgraphics
            null device
                      1

        To save to a specific file name, one should use :meth:`png` to set
        the output device to that file.  If this is done in the notebook, it
        must be done in the same cell as the plot itself::

            sage: filename = tmp_filename() + '.png'
            sage: r.png(filename='"%s"'%filename) # Note the double quotes in single quotes!; optional -- rgraphics
            NULL
            sage: x = r([1,2,3])
            sage: y = r([4,5,6])
            sage: r.plot(x,y)         # optional -- rgraphics
            null device
                      1
            sage: import os; os.unlink(filename) # For doctesting, we remove the file; optional -- rgraphics

        Please note that for more extensive use of R's plotting
        capabilities (such as the lattices package), it is advisable
        to either use an interactive plotting device or to use the
        notebook.  The following examples are not tested, because they
        differ depending on operating system::

            sage: r.X11() # not tested - opens interactive device on systems with X11 support
            sage: r.quartz() # not tested - opens interactive device on OSX
            sage: r.hist("rnorm(100)") # not tested - makes a plot
            sage: r.library("lattice") # not tested - loads R lattice plotting package
            sage: r.histogram(x = "~ wt | cyl", data="mtcars") # not tested - makes a lattice plot
            sage: r.dev_off() # not tested, turns off the interactive viewer

        In the notebook, one can use r.png() to open the device, but
        would need to use the following since R lattice graphics do
        not automatically print away from the command line::

            sage: filename = tmp_filename() + '.png' # Not needed in notebook, used for doctesting
            sage: r.png(filename='"%s"'%filename) # filename not needed in notebook, used for doctesting; optional -- rgraphics
            NULL
            sage: r.library("lattice")
            sage: r("print(histogram(~wt | cyl, data=mtcars))") # plot should appear; optional -- rgraphics
            sage: import os; os.unlink(filename) # We remove the file for doctesting, not needed in notebook; optional -- rgraphics
        """
        # We have to define this to override the plot function defined in the
        # superclass.
        from sage.server.support import EMBEDDED_MODE
        if EMBEDDED_MODE:
            self.setwd('"%s"'%os.path.abspath('.'))
        RFunction(self, 'plot')(*args, **kwds)
        return RFunction(self, 'dev.off')()

    def _strip_prompt(self, code):
        """
        Remove the standard R prompt from the beginning of lines in code.

        INPUT:

        - code -- a string

        OUTPUT: a string

        EXAMPLES::

            sage: s = '> a <- 2\n> b <- 3'
            sage: r._strip_prompt(s)
            ' a <- 2\n b <- 3'
        """
        return prompt_re.sub("", code)

    def eval(self, code, globals=None, locals=None, synchronize=True, *args, **kwds):
        """
        Evaluates a command inside the R interpreter and returns the output
        as a string.

        EXAMPLES::

            sage: r.eval('1+1')
            '[1] 2'
        """
        # TODO split code at ";" outside of quotes and send them as individual
        #      lines without ";".
        return Expect.eval(self, code, synchronize=synchronize, *args, **kwds)

    def _r_to_sage_name(self, s):
        """
        Returns a Sage/Python identifier from an R one.  This involves
        replacing periods with underscores, <- with __, and prepending
        _ in front of Python keywords.

        INPUT:

        - s -- a string

        OUTPUT: a string

        EXAMPLES::

            sage: f = r._r_to_sage_name
            sage: f('t.test')
            't_test'
            sage: f('attr<-')
            'attr__'
            sage: f('parent.env<-')
            'parent_env__'
            sage: f('class')
            'class_'
        """
        from keyword import iskeyword
        s = s.replace('.', '_')
        s = s.replace('<-', '__')
        if iskeyword(s):
            s += '_'
        return s

    def _sage_to_r_name(self, s):
        r"""
        The reverse of :meth:`_r_to_sage_name`.  See the docs for that method.

        EXAMPLES::

            sage: f = r._sage_to_r_name
            sage: f('t_test')
            't.test'
            sage: f('attr__')
            'attr<-'
            sage: f('parent_env__')
            'parent.env<-'
            sage: r._r_to_sage_name(f('parent_env__'))
            'parent_env__'
            sage: f('class_')
            'class'
        """
        if len(s) > 1 and s[-2:] == "__":
            s = s[:-2] + '<-'
        if len(s) > 0 and s[-1] == '_':
            s = s[:-1]
        s = s.replace('_', '.')
        return s

    def __getitem__(self, s):
        """
        Returns the RFunction with name s.

        INPUT:

        - s -- a string
        OUTPUT: RFunction -- the R function that in R has name s

        EXAMPLES::

            sage: r['as.data.frame']
            as.data.frame
            sage: r['print']
            print
        """
        return RFunction(self, s, r_name=True)

    def chdir(self, dir):
        """
        Changes the working directory to ``dir``

        INPUT:

        - ``dir`` -- the directory to change to.

        EXAMPLES::

            sage: import tempfile
            sage: tmpdir = tempfile.mkdtemp()
            sage: r.chdir(tmpdir)

        Check that ``tmpdir`` and ``r.getwd()`` refer to the same
        directory.  We need to use ``realpath()`` in case ``$TMPDIR``
        (by default ``/tmp``) is a symbolic link (see :trac:`10264`).

        ::

            sage: os.path.realpath(tmpdir) == sageobj(r.getwd())  # known bug (:trac:`9970`)
            True
        """
        self.execute('setwd(%r)' % dir)


# patterns for _sage_()
rel_re_param = re.compile('\s([\w\.]+)\s=')
rel_re_xrange = re.compile('([\d]+):([\d]+)')
rel_re_integer = re.compile('([^\d])([\d]+)L')
rel_re_terms = re.compile('terms\s*=\s*(.*?),')
rel_re_call = re.compile('call\s*=\s*(.*?)\),')

class RElement(ExpectElement):
    def __reduce__(self):
        """
        EXAMPLES::

            sage: a = r([1,2,3])
            sage: dumps(a)
            Traceback (most recent call last):
            ...
            NotImplementedError: pickling of R elements is not yet supported
        """
        raise NotImplementedError, "pickling of R elements is not yet supported"

    def trait_names(self):
        """
        Return a list of all methods of this object.

        .. note::

            Currently returns all R commands.

        EXAMPLES::

            sage: a = r([1,2,3])
            sage: t = a.trait_names()
            sage: len(t) > 200
            True
        """
        # TODO: rewrite it, just take methods(class=class(self))
        return self.parent().trait_names()

    def tilde(self, x):
        """
        The tilde regression operator in R.

        EXAMPLES::

            sage: x = r([1,2,3,4,5])
            sage: y = r([3,5,7,9,11])
            sage: a = r.lm( y.tilde(x) ) # lm( y ~ x )
            sage: d = a._sage_()
            sage: d['DATA']['coefficients']['DATA'][1]
            2
        """
        parent = self.parent()
        rx = parent(x)
        return parent.new("%s ~ %s"%(self.name(), rx.name()))

    stat_model = tilde

    def __len__(self):
        """
        Return the length of this object.

        OUTPUT: integer

        EXAMPLES::

            sage: x = r([10.4,5.6,3.1,6.4,21.7])
            sage: len(x)
            5
        """
        return int(self.parent().eval('dput(length(%s))'%self.name())[:-1] )

    def __getattr__(self, attrname):
        """
        Return attribute of this object, which is an R function with this object
        as the first input.

        INPUT:

        - attrname -- string

        OUTPUT: RFunctionElement

        EXAMPLES::

            sage: x = r([1,2,3])
            sage: length = x.length
            sage: type(length)
            <class 'sage.interfaces.r.RFunctionElement'>
            sage: length()
            [1] 3
        """
        self._check_valid()
        if attrname[:1] == "_":
            raise AttributeError
        return RFunctionElement(self, attrname)

    def __getitem__(self, n):
        """
        Return element(s) of self.

        INPUT:

        - n -- an integer, a tuple, a string that makes sense to R, or an RElement

        OUTPUT: RElement

        EXAMPLES::

            sage: x = r([10.4,5.6,3.1,6.4,21.7])
            sage: x[0]
            numeric(0)
            sage: x[1]
            [1] 10.4
            sage: x[-1]
            [1] 5.6  3.1  6.4 21.7
            sage: x[-2]
            [1] 10.4  3.1  6.4 21.7
            sage: x[-3]
            [1] 10.4  5.6  6.4 21.7
            sage: x['c(2,3)']
            [1]  5.6 3.1
            sage: key = r.c(2,3)
            sage: x[key]
            [1]  5.6 3.1
            sage: m = r.array('1:3',r.c(2,4,2))
            sage: m
            , , 1
                 [,1] [,2] [,3] [,4]
            [1,]    1    3    2    1
            [2,]    2    1    3    2
            , , 2
                 [,1] [,2] [,3] [,4]
            [1,]    3    2    1    3
            [2,]    1    3    2    1
            sage: m[1,2,2]
            [1] 2
            sage: m[1,r.c(1,2),1]
            [1] 1 3
        """
        P = self._check_valid()
        if isinstance(n, basestring):
            n = n.replace('self', self._name)
            return P.new('%s[%s]'%(self._name, n))
        elif (hasattr(n,'parent') and n.parent() is P): # the key is RElement itself
            return P.new('%s[%s]'%(self._name, n.name()))
        elif not isinstance(n,tuple):
            return P.new('%s[%s]'%(self._name, n))
        else:
            L = []
            for i in xrange(len(n)):
                if (hasattr(n[i],'parent') and n[i].parent() is P):
                    L.append(n[i].name())
                else:
                    L.append(str(n[i]))
            return P.new('%s[%s]'%(self._name, ','.join(L)))

    def __nonzero__(self):
        """
        Implements bool(self).

        .. note::

            bool(self) will only return True if self == 0 contains a FALSE in its representation.

        EXAMPLES::

            sage: x = r([10.4,5.6,3.1,6.4,21.7])
            sage: bool(x)
            True
            sage: y = r([0,0,0,0])
            sage: bool(y)
            False
            sage: bool(r(0))
            False
            sage: bool(r(1))
            True
        """
        return "FALSE" in repr(self==0)

    def _comparison(self, other, symbol):
        """
        Used to implement comparison of two objects.

        INPUT:

        - other -- RElement
        - symbol -- string

        OUTPUT: RElement -- output is an R element; not a bool!

        TESTS::

            sage: x = r([10.4,5.6,3.1,6.4,21.7])
            sage: x._comparison(10.4, "==")
            [1] TRUE FALSE FALSE FALSE FALSE
        """
        P = self.parent()
        other = P(other)
        return P('%s %s %s'%(self.name(), symbol, other.name()))

    def __eq__(self, other):
        """
        Equality testing term by term.

        INPUT:

        - other -- RElement

        OUTPUT: RElement -- an R element; not a bool!

        EXAMPLES:

        Notice that comparison is term by term and returns an R element. ::

            sage: x = r([10.4,5.6,3.1,6.4,21.7])
            sage: x == 10.4
            [1] TRUE FALSE FALSE FALSE FALSE
        """
        return self._comparison(other, "==")

    def __lt__(self, other):
        """
        Less than testing term by term.

        INPUT:

        - other -- RElement

        OUTPUT: RElement -- an R element; not a bool!

        EXAMPLES:

        Notice that comparison is term by term and returns an R element. ::

            sage: x = r([10.4,5.6,3.1,6.4,21.7])
            sage: x < 7
            [1] FALSE  TRUE  TRUE  TRUE FALSE
        """
        return self._comparison(other, "<")

    def __gt__(self, other):
        """
        Greater than testing term by term.

        INPUT:

        - other -- RElement

        OUTPUT: RElement -- an R element; not a bool!

        EXAMPLES:

        Notice that comparison is term by term and returns an R element. ::

            sage: x = r([10.4,5.6,3.1,6.4,21.7])
            sage: x > 8
            [1] TRUE FALSE FALSE FALSE  TRUE
        """
        return self._comparison(other, ">")

    def __le__(self, other):
        """
        Less than or equal testing term by term.

        INPUT:

        - other -- RElement

        OUTPUT: RElement -- an R element; not a bool!

        EXAMPLES::

            sage: x = r([10.4,5.6,3.1,6.4,21.7])
            sage: x <= 10.4
            [1] TRUE  TRUE  TRUE  TRUE FALSE
        """
        return self._comparison(other, "<=")

    def __ge__(self, other):
        """
        Greater than or equal testing term by term.

        INPUT:

        - other -- RElement

        OUTPUT: RElement -- an R element; not a bool!

        EXAMPLES::

            sage: x = r([10.4,5.6,3.1,6.4,21.7])
            sage: x >= 10.4
            [1] TRUE FALSE FALSE FALSE  TRUE
        """
        return self._comparison(other, ">=")

    def __ne__(self, other):
        """
        Not equal testing term by term.

        INPUT:

        - other -- RElement

        OUTPUT: RElement -- an R element; not a bool!

        EXAMPLES::

            sage: x = r([10.4,5.6,3.1,6.4,21.7])
            sage: x != 10.4
            [1] FALSE  TRUE  TRUE  TRUE  TRUE

        """
        return self._comparison(other, "!=")

    def __cmp__(self, other):
        r"""
        Return 0, 1, or -1 depending on how self and other compare.

        This is *not* called by the comparison operators, which
        do term-by-term comparison and return R elements.

        INPUT:

        - self, other -- R elements

        OUTPUT: 0, 1, or -1

        EXAMPLES::

            sage: one = r(1)
            sage: two = r(2)
            sage: cmp(one,one)
            0
            sage: cmp(one,two)
            -1
            sage: cmp(two,one)
            1
        """
        P = self.parent()
        if P.eval("%s %s %s"%(self.name(), P._equality_symbol(),
                                 other.name())) == P._true_symbol():
            return 0
        elif P.eval("%s %s %s"%(self.name(), P._lessthan_symbol(), other.name())) == P._true_symbol():
            return -1
        elif P.eval("%s %s %s"%(self.name(), P._greaterthan_symbol(), other.name())) == P._true_symbol():
            return 1
        else:
            return -1  # everything is supposed to be comparable in Python, so we define
                       # the comparison thus when no comparable in interfaced system.


    def dot_product(self, other):
        """
        Implements the notation self . other.

        INPUT:

        - self, other -- R elements

        OUTPUT: R element

        EXAMPLES::

            sage: c = r.c(1,2,3,4)
            sage: c.dot_product(c.t())
                 [,1] [,2] [,3] [,4]
            [1,]    1    2    3    4
            [2,]    2    4    6    8
            [3,]    3    6    9   12
            [4,]    4    8   12   16

            sage: v = r([3,-1,8])
            sage: v.dot_product(v)
                 [,1]
            [1,]   74
        """
        P = self._check_valid()
        Q = P(other)
        # the R operator is %*% for matrix multiplication
        return P('%s %%*%% %s'%(self.name(), Q.name()))

    def _subs_dots(self, x):
        """
        Replace dots by underscores; used internally to implement
        conversation from R expression to Sage objects.

        INPUT:

        - x -- regular expression match: ``<type '_sre.SRE_Match'>``

        OUTPUT: string

        EXAMPLES::

            sage: import re
            sage: a = r([1,2,3])
            sage: rel_re_param = re.compile('\s([\w\.]+)\s=')
            sage: rel_re_param.sub(a._subs_dots, ' test.test =')
             ' test_test ='
        """
        return x.group().replace('.','_')

    def _subs_xrange(self, x):
        """
        Change endpoints of xranges.  This is used internally in the
        code for converting R expressions to Sage objects.

        INPUT:

        - x -- regular expression match: ``<type '_sre.SRE_Match'>``

        OUTPUT: string

        EXAMPLES::

            sage: import re
            sage: a = r([1,2,3])
            sage: rel_re_xrange = re.compile('([\d]+):([\d]+)')
            sage: rel_re_xrange.sub(a._subs_xrange, ' 1:10')
            ' xrange(1,11)'
        """
        g = x.groups()
        g1 = int(g[1]) + 1
        return 'xrange(%s,%s)'%(g[0],g1)

    def _subs_integer(self, x):
        """
        Replaces strings like 'dL' with 'Integer(d)' where d is some
        integer.  This is used internally in the code for converting R
        expressions to Sage objects.

        EXAMPLES::

            sage: import re
            sage: a = r([1,2,3])
            sage: rel_re_integer = re.compile('([^\d])([\d]+)L')
            sage: rel_re_integer.sub(a._subs_integer, ' 1L 2L')
            ' Integer(1) Integer(2)'
            sage: rel_re_integer.sub(a._subs_integer, '1L 2L')
            '1L Integer(2)'

        """
        return '%sInteger(%s)'%x.groups()

    def _convert_nested_r_list(self, exp):
        """
        Converts a string representing a (possibly) nested list in R
        to a (possibly) nested Python list.  This is used internally
        in the code for converting R expressions to Sage objects.

        INPUT:

        - exp -- a string

        OUTPUT: a string

        EXAMPLES::

            sage: a = r([1,2,3])
            sage: s = 'c(1, 2, 3)'
            sage: a._convert_nested_r_list(s)
            '[1, 2, 3]'
        """
        from re import compile as re_compile
        from re import split   as re_split
        splt = re_compile('(c\(|\(|\))') # c( or ( or )
        lvl = 0
        ret = []
        for token in re_split(splt, exp):
            if token == 'c(':
                ret.append('[')
                lvl += 1
            elif token == '(':
                ret.append(token)
                if lvl > 0 : lvl += 1
            elif token == ')':
                if lvl == 1:
                    ret.append(']')
                    lvl -= 1
                else:
                    ret.append(token)
                    if lvl > 0:
                        lvl -= 1
            else:
                ret.append(token)

        return ''.join(ret)


    def _r_list(self, *args, **kwds):
        """
        This is used internally in the code for converting R
        expressions to Sage objects.

        EXAMPLES::

            sage: a = r([1,2,3])
            sage: list(sorted(a._r_list(1,2,3,k=5).items()))
            [('#0', 1), ('#1', 2), ('#2', 3), ('k', 5)]
        """
        ret = dict(kwds)
        i = 0
        for k in args:
            ret['#%s'%i] = k
            i += 1
        return ret

    def _r_structure(self, __DATA__, **kwds):
        """
        This is used internally in the code for converting R
        expressions to Sage objects.

        EXAMPLES::

            sage: a = r([1,2,3])
            sage: d = a._r_structure('data', a=1, b=2)
            sage: list(sorted(d.items()))
            [('DATA', 'data'), ('a', 1), ('b', 2)]
            sage: a._r_structure([1,2,3,4], _Dim=(2,2))
            [1 3]
            [2 4]

        """
        if '_Dim' in kwds: #we have a matrix
            # TODO what about more than 2 dimensions?
            #      additional checks!!
            try:
                from sage.matrix.constructor import matrix
                d = kwds.get('_Dim')
                # TODO: higher dimensions? happens often in statistics
                if len(d) != 2:
                    raise TypeError
                #since R does it the other way round, we assign
                #transposed and then transpose the matrix :)
                m = matrix(d[1], d[0], [i for i in __DATA__])
                return m.transpose()
            except TypeError:
                pass
        d = dict(DATA=__DATA__)
        d.update(kwds)
        return d

    def _sage_(self):
        r"""
        Returns Sage representation of the R object.

        R objects are basic C structures, of different kind, that can
        be stacked together.  This is similar to Python lists with
        variable objects, including lists of lists.  If R lists have
        names, they are translated to a Python dictionary, with anonymous
        list entries called ``#{number}``.

        OUTPUT: object -- Python object

        EXAMPLES::

            sage: rs = r.summary(r.c(1,4,3,4,3,2,5,1))
            sage: d = rs._sage_()
            sage: list(sorted(d.items()))
            [('DATA', [1, 1.75, 3, 2.875, 4, 5]),
             ('_Names', ['Min.', '1st Qu.', 'Median', 'Mean', '3rd Qu.', 'Max.']),
             ('_r_class', ['summaryDefault', 'table'])]
        """
        self._check_valid()
        P = self.parent()

        # This is the core of the trick: using dput
        # dput prints out the internal structure of R's data elements
        # options via .deparseOpts(control=...)
        # TODO: dput also works with a file, if things get huge!
        # [[However, using a file for output often isn't necessary
        # since pipe buffering works pretty well for output.
        # That said, benchmark this.  -- William Stein]]
        exp = P.eval('dput(%s)'%self.name())

        # Preprocess expression
        # example what this could be:
        # structure(list(statistic = structure(0.233549683248457, .Names = "t"),
        # parameter = structure(5.58461538461538, .Names = "df"), p.value = 0.823657802106985,
        # conf.int = structure(c(-2.41722062247400, 2.91722062247400
        # ), conf.level = 0.95), estimate = structure(c(2.75, 2.5), .Names = c("mean of x",
        # "mean of y")), null.value = structure(0, .Names = "difference in means"),
        # alternative = "two.sided", method = "Welch Two Sample t-test",
        # data.name = "c(1, 2, 3, 5) and c(1, 2, 3, 4)"), .Names = c("statistic",
        # "parameter", "p.value", "conf.int", "estimate", "null.value",
        # "alternative", "method", "data.name"), class = "htest")

        # R's structure (from help):
        # structure(.Data, ...)
        #    .Data: an object which will have various attributes attached to it.
        #    ...: attributes, specified in 'tag=value' form, which will be
        #         attached to data.
        #For historical reasons (these names are used when deparsing),
        # attributes '".Dim"', '".Dimnames"', '".Names"', '".Tsp"' and
        # '".Label"' are renamed to '"dim"', '"dimnames"', '"names"',
        # '"tsp"' and '"levels"'.



        # we want this in a single line
        exp.replace('\n','')
        exp = "".join(exp.split("\n"))

        # python compatible parameters
        exp = rel_re_param.sub(self._subs_dots, exp)

        # Rename class since it is a Python keyword
        exp = re.sub(' class = ', ' _r_class = ',exp)

        # Change 'structure' to '_r_structure'
        # TODO: check that we are outside of quotes ""
        exp = re.sub(' structure\(', ' _r_structure(', exp)
        exp = re.sub('^structure\(', '_r_structure(', exp) #special case

        # Change 'list' to '_r_list'
        exp = re.sub(' list\(', ' _r_list(', exp)
        exp = re.sub('\(list\(', '(_r_list(', exp)

        # Change 'a:b' to 'xrange(a,b+1)'
        exp = rel_re_xrange.sub(self._subs_xrange, exp)

        # Change 'dL' to 'Integer(d)'
        exp = rel_re_integer.sub(self._subs_integer, exp)

        # Wrap the right hand side of terms = ... in quotes since it
        # has a ~ in it.
        exp = rel_re_terms.sub(r'terms = "\1",', exp)


        # Change call = ..., to call = "...",
        exp = rel_re_call.sub(r'call = "\1",', exp)

        # seems to work for
        # rr = r.summary(r.princomp(r.matrix(r.c(1,2,3,4,3,4,1,2,2),4)))
        # rr._sage_()
        # but the call expression get's evaluated. why?!? TODO


        # translation:
        # c is an ordered list
        # list is a dictionary (where _Names give the entries names.
        #    map entries in names to (value, name) in each entry?
        # structure is .. see above .. strucuture(DATA,**kw)
        # TODO: thinking of just replacing c( with ( to get a long tuple?


        exp = self._convert_nested_r_list(exp)

        # Set up the globals
        globs = {'NA':None, 'NULL':None, 'FALSE':False, 'TRUE':True,
                 '_r_list':self._r_list, '_r_structure':self._r_structure,
                 'Integer':sage.rings.integer.Integer,
                 'character':str}

        return eval(exp, globs, globs)


    def _latex_(self):
        r"""
        Return LaTeX representation of this R object.

        This calls the ``latex`` command in R.

        OUTPUT: a latex expression (basically a string)

        EXAMPLES::

            sage: latex(r(2))  # optional - Hmisc R package
            2
        """
        from sage.misc.latex import LatexExpr
        self._check_valid()
        P = self.parent()
        # latex is in Hmisc, this is currently not part of Sage's R!!!
        try:
            P.library('Hmisc')
        except ImportError:
            raise RuntimeError, "The R package 'Hmisc' is required for R to LaTeX conversion, but it is not available."
        return LatexExpr(P.eval('latex(%s, file="");'%self.name()))



class RFunctionElement(FunctionElement):
    def __reduce__(self):
        """
        EXAMPLES::

            sage: a = r([1,2,3])
            sage: a.mean
            mean
            sage: dumps(a.mean)
            Traceback (most recent call last):
            ...
            NotImplementedError: pickling of R element methods is not yet supported
        """
        raise NotImplementedError, "pickling of R element methods is not yet supported"

    def _sage_doc_(self):
        """
        Returns the help for self as a string.

        EXAMPLES::

            sage: a = r([1,2,3])
            sage: length = a.length
            sage: print length._sage_doc_()
            length                 package:base                 R Documentation
            ...
            <BLANKLINE>
        """
        M = self._obj.parent()
        return M.help(self._name)

    def _sage_src_(self):
        """
        Returns the source code of self.

        EXAMPLES::

            sage: a = r([1,2,3])
            sage: length = a.length
            sage: print length._sage_src_()
            function (x)  .Primitive("length")
        """
        M = self._obj.parent()
        return M.source(self._name)

    def __call__(self, *args, **kwds):
        """
        EXAMPLES::

            sage: a = r([1,2,3])
            sage: length = a.length
            sage: length()
            [1] 3
        """
        return self._obj.parent().function_call(self._name, args=[self._obj] + list(args), kwds=kwds)


class RFunction(ExpectFunction):
    def __init__(self, parent, name, r_name=None):
        """
        A Function in the R interface.

        INPUT:

        - parent -- the R interface
        - name -- the name of the function for Python
        - r_name -- the name of the function in R itself (which can have dots in it)

        EXAMPLES::

            sage: length = r.length
            sage: type(length)
            <class 'sage.interfaces.r.RFunction'>
            sage: loads(dumps(length))
            length
        """
        self._parent = parent
        if r_name:
            self._name = name
        else:
            self._name = parent._sage_to_r_name(name)

    def __cmp__(self, other):
        """
        EXAMPLES::

            sage: r.mean == loads(dumps(r.mean))
            True
            sage: r.mean == r.lr
            False
        """
        if not isinstance(other, RFunction):
            return cmp(type(self), type(other))
        return cmp(self._name, other._name)

    def _sage_doc_(self):
        """
        Returns the help for self.

        EXAMPLES::

            sage: length = r.length
            sage: print length._sage_doc_()
            length                 package:base                 R Documentation
            ...
            <BLANKLINE>
        """
        M = self._parent
        return M.help(self._name)

    def _sage_src_(self):
        """
        Returns the source of self.

        EXAMPLES::

            sage: length = r.length
            sage: print length._sage_src_()
            function (x)  .Primitive("length")

        """
        M = self._parent
        return M.source(self._name)

    def __call__(self, *args, **kwds):
        """
        EXAMPLES::

            sage: length = r.length
            sage: length([1,2,3])
            [1] 3
        """
        return self._parent.function_call(self._name, args=list(args), kwds=kwds)

def is_RElement(x):
    """
    Return True if x is an element in an R interface.

    INPUT:

    - x -- object

    OUTPUT: bool

    EXAMPLES::

        sage: from sage.interfaces.r import is_RElement
        sage: is_RElement(2)
        False
        sage: is_RElement(r(2))
        True
    """
    return isinstance(x, RElement)

# An instance of R
r = R()

def reduce_load_R():
    """
    Used for reconstructing a copy of the R interpreter from a pickle.

    EXAMPLES::

        sage: from sage.interfaces.r import reduce_load_R
        sage: reduce_load_R()
        R Interpreter
    """
    return r

import os
def r_console():
    """
    Spawn a new R command-line session.

    EXAMPLES::

        sage: r.console()                    # not tested
            R version 2.6.1 (2007-11-26)
            Copyright (C) 2007 The R Foundation for Statistical Computing
            ISBN 3-900051-07-0
            ...
    """
    # This will only spawn local processes
    os.system('R --vanilla')

def r_version():
    """
    Return the R version.

    EXAMPLES::

        sage: r_version() # not tested
        ((3, 0, 1), 'R version 3.0.1 (2013-05-16)')
        sage: rint, rstr = r_version()
        sage: rint[0] >= 3
        True
        sage: rstr.startswith('R version')
        True
    """
    return r.version()

class HelpExpression(str):
    """
    Used to improve printing of output of r.help.
    """
    def __repr__(self):
        """
        Return string representation of self.

        OUTPUT: string

        EXAMPLES::

            sage: a = sage.interfaces.r.HelpExpression("This\nis\nR!")
            sage: type(a)
            <class 'sage.interfaces.r.HelpExpression'>
            sage: a
            This
            is
            R!
        """
        return self.__str__()

