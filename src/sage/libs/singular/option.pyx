"""
libSingular: Options

Singular uses a set of global options to determine verbosity and the
behavior of certain algorithms. We provide an interface to these
options in the most 'natural' python-ic way. Users who do not wish to
deal with Singular functions directly usually do not have to worry
about this interface or Singular options in general since this is
taken care of by higher level functions.

We compute a Groebner basis for Cyclic-5 in two different contexts::

    sage: P.<a,b,c,d,e> = PolynomialRing(GF(127))
    sage: I = sage.rings.ideal.Cyclic(P)
    sage: import sage.libs.singular.function_factory
    sage: std = sage.libs.singular.function_factory.ff.std

By default, tail reductions are performed::

    sage: from sage.libs.singular.option import opt, opt_ctx
    sage: opt['red_tail']
    True
    sage: red = std(I)[-1]; red
    d^2*e^6 + 28*b*c*d + ...

If we don't want this, we can create an option context, which disables
this::

    sage: with opt_ctx(red_tail=False, red_sb=False):
    ....:    notred = std(I)[-1]; notred
    d^2*e^6 + ...
    sage: len(list(red)) < len(list(notred))
    True

However, this does not affect the global state::

    sage: opt['red_tail']
    True

On the other hand, any assignment to an option object will immediately
change the global state::

    sage: opt['red_tail'] = False
    sage: opt['red_tail']
    False
    sage: opt['red_tail'] = True
    sage: opt['red_tail']
    True

Assigning values within an option context, only affects this context::

    sage: with opt_ctx:
    ....:    opt['red_tail'] = False

    sage: opt['red_tail']
    True

Option contexts can also be safely stacked::

    sage: with opt_ctx:
    ....:     opt['red_tail'] = False
    ....:     print(opt)
    ....:     with opt_ctx:
    ....:         opt['red_through'] = False
    ....:         print(opt)
    general options for libSingular (current value 0x00000082)
    general options for libSingular (current value 0x00000002)

    sage: print(opt)
    general options for libSingular (current value 0x02000082)

Furthermore, the integer valued options ``deg_bound`` and
``mult_bound`` can be used::

    sage: R.<x,y> = QQ[]
    sage: I = R*[x^3+y^2,x^2*y+1]
    sage: opt['deg_bound'] = 2
    sage: std(I)
    [x^2*y + 1, x^3 + y^2]
    sage: opt['deg_bound'] = 0
    sage: std(I)
    [y^3 - x, x^2*y + 1, x^3 + y^2]

The same interface is available for verbosity options::

    sage: from sage.libs.singular.option import opt_verb
    sage: opt_verb['not_warn_sb']
    False
    sage: opt.reset_default()  # needed to avoid side effects
    sage: opt_verb.reset_default()  # needed to avoid side effects

AUTHOR:

- Martin Albrecht (2009-08): initial implementation
- Martin Albrecht (2010-01): better interface, verbosity options
- Simon King (2010-07): Python-ic option names; deg_bound and mult_bound
"""
# ****************************************************************************
#       Copyright (C) 2010 Martin Albrecht <M.R.Albrecht@rhul.ac.uk>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.libs.singular.decl cimport singular_options, singular_verbose_options, Kstd1_deg, Kstd1_mu

from sage.libs.singular.decl cimport OPT_PROT, OPT_REDSB, OPT_NOT_BUCKETS, OPT_NOT_SUGAR, OPT_SUGARCRIT, OPT_REDTHROUGH, OPT_DEGBOUND, OPT_MULTBOUND
from sage.libs.singular.decl cimport OPT_RETURN_SB, OPT_FASTHC, OPT_OLDSTD, OPT_REDTAIL, OPT_INTSTRATEGY, OPT_NOTREGULARITY
from sage.libs.singular.decl cimport OPT_WEIGHTM, Sy_bit

from sage.libs.singular.decl cimport V_SHOW_MEM, V_YACC, V_REDEFINE, V_READING, V_LOAD_LIB, V_DEBUG_LIB
from sage.libs.singular.decl cimport V_LOAD_PROC, V_DEF_RES, V_SHOW_USE, V_IMAP, V_PROMPT
from sage.libs.singular.decl cimport V_NSB, V_CONTENTSB, V_CANCELUNIT

_options_py_to_singular={'return_sb':'returnSB',
                         'fast_hc':'fastHC',
                         'inf_red_tail':'infRedTail',
                         'int_strategy':'intStrategy',
                         'not_regularity':'notRegularity',
                         'not_sugar':'notSugar',
                         'not_buckets':'notBuckets',
                         'qring_nf':'qringNF',
                         'redsb':'redSB',
                         'red_sb':'redSB',
                         'red_tail':'redTail',
                         'red_through':'redThrough',
                         'sugar_crit':'sugarCrit',
                         'weight_m':'weightM',
                         'content_sb':'contentSB',
                         'mult_bound':'multBound',
                         'deg_bound':'degBound',
                         'imap':'Imap',
                         'debug_lib':'debugLib',
                         'def_res':'defRes',
                         'load_lib':'loadLib',
                         'load_proc':'loadProc',
                         'not_warn_sb':'notWarnSB'}

cdef class LibSingularOptions_abstract:
    """
    Abstract Base Class for libSingular options.
    """
    cdef unsigned int *global_options
    cdef object name
    cdef object name_map

    def __init__(self, **kwds):
        """
        INPUT:

        - ``**kwds`` -- all keyword parameters are immediately applied

        EXAMPLES::

            sage: from sage.libs.singular.option import LibSingularOptions
            sage: opt = LibSingularOptions(redTail=False)
            sage: opt['redTail']
            False
            sage: opt['redTail'] = True
            sage: opt['redTail']
            True
            sage: opt['deg_bound'] = 2

        The options can be named in Python or Singular style::

            sage: opt['degBound']
            2
        """
        for k,v in kwds.iteritems():
            self[k] = v

    def __getitem__(self, name):
        """
        EXAMPLES::

            sage: from sage.libs.singular.option import opt
            sage: opt['red_tail']
            True
            sage: opt['deg_bound'] = 2

        The options can be named in Python or Singular style::

            sage: opt['degBound']
            2
            sage: opt.reset_default()  # needed to avoid side effects
        """
        name = _options_py_to_singular.get(name,name)
        if name == "degBound":
            if bool(self.global_options[0] & self.name_map[name]):
                return Kstd1_deg
            return int(0)
        if name == "multBound":
            if bool(self.global_options[0] & self.name_map[name]):
                return Kstd1_mu
            return int(0)
        try:
            return bool(self.global_options[0] & self.name_map[name])
        except KeyError:
            raise NameError("Option '%s' unknown."%(name,))

    def __setitem__(self, name, value):
        """
        EXAMPLES::

            sage: from sage.libs.singular.option import opt, opt_ctx
            sage: opt['redTail']
            True
            sage: with opt_ctx:
            ....:    opt['redTail'] = False
            ....:    opt['redTail']
            False
            sage: opt['red_tail']
            True
            sage: opt.reset_default()  # needed to avoid side effects
        """
        name = _options_py_to_singular.get(name,name)
        try:
            if value:
                self.global_options[0] = self.global_options[0] | self.name_map[name]
            else:
                self.global_options[0] = self.global_options[0] & ~self.name_map[name]
            if name == 'degBound':
                global Kstd1_deg
                Kstd1_deg = value
            elif name == 'multBound':
                global Kstd1_mu
                Kstd1_mu = value
        except KeyError:
            raise NameError("Option '%s' unknown."%(name,))

    def __int__(self):
        """
        EXAMPLES::

            sage: from sage.libs.singular.option import opt
            sage: hex(int(opt))
            '0x6000082'
        """
        return int(self.global_options[0])

    def save(self):
        """
        Return a triple of integers that allow reconstruction of the options.

        EXAMPLES::

            sage: from sage.libs.singular.option import opt
            sage: opt['deg_bound']
            0
            sage: opt['red_tail']
            True
            sage: s = opt.save()
            sage: opt['deg_bound'] = 2
            sage: opt['red_tail'] = False
            sage: opt['deg_bound']
            2
            sage: opt['red_tail']
            False
            sage: opt.load(s)
            sage: opt['deg_bound']
            0
            sage: opt['red_tail']
            True
            sage: opt.reset_default()  # needed to avoid side effects
        """
        return (int(self.global_options[0]), self['deg_bound'], self['mult_bound'])

    def load(self, value=None):
        """
        EXAMPLES::

            sage: from sage.libs.singular.option import opt as sopt
            sage: bck = sopt.save(); hex(bck[0]), bck[1], bck[2]
            ('0x6000082', 0, 0)
            sage: sopt['redTail'] = False
            sage: hex(int(sopt))
            '0x4000082'
            sage: sopt.load(bck)
            sage: sopt['redTail']
            True
        """
        if value is None:
            value = (None,0,0)
        self.global_options[0] = int(value[0])
        global Kstd1_deg
        global Kstd1_mu
        Kstd1_deg = value[1]
        Kstd1_mu  = value[2]

    def __repr__(self):
        """
        EXAMPLES::

            sage: from sage.libs.singular.option import opt as sopt
            sage: sopt
            general options for libSingular (current value 0x06000082)
        """
        return "%s options for libSingular (current value 0x%08x)"%(self.name, self.global_options[0])


cdef class LibSingularOptions(LibSingularOptions_abstract):
    """
    Pythonic Interface to libSingular's options.

    Supported options are:

    - ``return_sb`` or ``returnSB`` -- the functions ``syz``,
      ``intersect``, ``quotient``, ``modulo`` return a standard
      base instead of a generating set if ``return_sb``
      is set. This option should not be used for ``lift``.

    - ``fast_hc`` or ``fastHC`` -- tries to find the highest corner
      of the staircase (HC) as fast as possible during a standard
      basis computation (only used for local orderings).

    - ``int_strategy`` or ``intStrategy`` -- avoids division of
      coefficients during standard basis computations. This option
      is ring dependent. By default, it is set for rings with
      characteristic 0 and not set for all other rings.

    - ``lazy`` -- uses a more lazy approach in std computations, which
      was used in SINGULAR version before 2-0 (and which may lead to
      faster or slower computations, depending on the example).

    - ``length`` -- select shorter reducers in std computations

    - ``not_regularity`` or ``notRegularity`` -- disables the
      regularity bound for ``res`` and ``mres``

    - ``not_sugar`` or ``notSugar`` -- disables the sugar strategy
      during standard basis computation

    - ``not_buckets`` or ``notBuckets`` -- disables the bucket
      representation of polynomials during standard basis
      computations. This option usually decreases the memory
      usage but increases the computation time. It should only
      be set for memory-critical standard basis computations.

    - ``old_std`` or ``oldStd`` -- uses a more lazy approach in std
      computations, which was used in SINGULAR version before 2-0
      (and which may lead to faster or slower computations, depending
      on the example).

    - ``prot`` -- shows protocol information indicating the progress
      during the following computations: ``facstd``, ``fglm``,
      ``groebner``, ``lres``, ``mres``, ``minres``, ``mstd``,
      ``res``, ``slimgb``, ``sres``, ``std``, ``stdfglm``,
      ``stdhilb``, ``syz``.

    - ``red_sb`` or ``redSB`` -- computes a reduced standard basis in
      any standard basis computation

    - ``red_tail`` or ``redTail`` -- reduction of the tails of
      polynomials during standard basis computations. This option
      is ring dependent. By default, it is set for rings with global
      degree orderings and not set for all other rings.

    - ``red_through`` or ``redThrough`` -- for inhomogeneous input,
      polynomial reductions during standard basis computations are
      never postponed, but always finished through. This option is
      ring dependent. By default, it is set for rings with global
      degree orderings and not set for all other rings.

    - ``sugar_crit`` or ``sugarCrit`` -- uses criteria similar to the
      homogeneous case to keep more useless pairs

    - ``weight_m`` or ``weightM`` -- automatically computes suitable
      weights for the weighted ecart and the weighted sugar method

    In addition, two integer valued parameters are supported, namely:

    - ``deg_bound`` or ``degBound`` -- the standard basis computation
      is stopped if the total (weighted) degree exceeds ``deg_bound``.
      ``deg_bound`` should not be used for a global ordering with
      inhomogeneous input. Reset this bound by setting ``deg_bound``
      to 0. The exact meaning of "degree" depends on the ring ordering
      and the command: ``slimgb`` uses always the total degree with
      weights 1, ``std`` does so for block orderings, only.

    - ``mult_bound`` or ``multBound`` -- the standard basis computation
      is stopped if the ideal is zero-dimensional in a ring with local
      ordering and its multiplicity is lower than ``mult_bound``.
      Reset this bound by setting ``mult_bound`` to 0.

    EXAMPLES::

        sage: from sage.libs.singular.option import LibSingularOptions
        sage: libsingular_options = LibSingularOptions()
        sage: libsingular_options
        general options for libSingular (current value 0x06000082)

    Here we demonstrate the intended way of using libSingular options::

        sage: R.<x,y> = QQ[]
        sage: I = R*[x^3+y^2,x^2*y+1]
        sage: I.groebner_basis(deg_bound=2)
        [x^3 + y^2, x^2*y + 1]
        sage: I.groebner_basis()
        [x^3 + y^2, x^2*y + 1, y^3 - x]

    The option ``mult_bound`` is only relevant in the local case::

        sage: from sage.libs.singular.option import opt
        sage: Rlocal.<x,y,z> = PolynomialRing(QQ, order='ds')
        sage: x^2<x
        True
        sage: J = [x^7+y^7+z^6,x^6+y^8+z^7,x^7+y^5+z^8, x^2*y^3+y^2*z^3+x^3*z^2,x^3*y^2+y^3*z^2+x^2*z^3]*Rlocal
        sage: J.groebner_basis(mult_bound=100)
        [x^3*y^2 + y^3*z^2 + x^2*z^3, x^2*y^3 + x^3*z^2 + y^2*z^3, y^5, x^6 + x*y^4*z^5, x^4*z^2 - y^4*z^2 - x^2*y*z^3 + x*y^2*z^3, z^6 - x*y^4*z^4 - x^3*y*z^5]
        sage: opt['red_tail'] = True # the previous commands reset opt['red_tail'] to False
        sage: J.groebner_basis()
        [x^3*y^2 + y^3*z^2 + x^2*z^3, x^2*y^3 + x^3*z^2 + y^2*z^3, y^5, x^6, x^4*z^2 - y^4*z^2 - x^2*y*z^3 + x*y^2*z^3, z^6, y^4*z^3 - y^3*z^4 - x^2*z^5, x^3*y*z^4 - x^2*y^2*z^4 + x*y^3*z^4, x^3*z^5, x^2*y*z^5 + y^3*z^5, x*y^3*z^5]
    """
    def __init__(self, **kwds):
        """
        Create a new option interface.

        EXAMPLES::

            sage: from sage.libs.singular.option import LibSingularOptions
            sage: libsingular_options = LibSingularOptions()
            sage: libsingular_options
            general options for libSingular (current value 0x...)
        """
        self.global_options = &singular_options
        self.name = "general"
        self.name_map = {"prot":          Sy_bit(OPT_PROT),
                         "redSB":         Sy_bit(OPT_REDSB),
                         "notBuckets":    Sy_bit(OPT_NOT_BUCKETS),
                         "notSugar":      Sy_bit(OPT_NOT_SUGAR),
                         "sugarCrit":     Sy_bit(OPT_SUGARCRIT),
                         "redThrough":    Sy_bit(OPT_REDTHROUGH),
                         "returnSB":      Sy_bit(OPT_RETURN_SB),
                         "fastHC":        Sy_bit(OPT_FASTHC),
                         "oldStd":        Sy_bit(OPT_OLDSTD),
                         "redTail":       Sy_bit(OPT_REDTAIL),
                         "intStrategy":   Sy_bit(OPT_INTSTRATEGY),
                         "notRegularity": Sy_bit(OPT_NOTREGULARITY),
                         "weightM":       Sy_bit(OPT_WEIGHTM),
                         "degBound":      Sy_bit(OPT_DEGBOUND),
                         "multBound":     Sy_bit(OPT_MULTBOUND)}
        LibSingularOptions_abstract.__init__(self, **kwds)

    def reset_default(self):
        """
        Reset libSingular's default options.

        EXAMPLES::

            sage: from sage.libs.singular.option import opt
            sage: opt['red_tail']
            True
            sage: opt['red_tail'] = False
            sage: opt['red_tail']
            False
            sage: opt['deg_bound']
            0
            sage: opt['deg_bound'] = 2
            sage: opt['deg_bound']
            2
            sage: opt.reset_default()
            sage: opt['red_tail']
            True
            sage: opt['deg_bound']
            0
        """
        from sage.libs.singular.singular import _saved_options
        self.load(_saved_options)


#############

cdef class LibSingularVerboseOptions(LibSingularOptions_abstract):
    """
    Pythonic Interface to libSingular's verbosity options.

    Supported options are:

    - ``mem`` -- shows memory usage in square brackets
    - ``yacc`` -- only available in debug version
    - ``redefine`` -- warns about variable redefinitions
    - ``reading`` -- shows the number of characters read from a file
    - ``loadLib`` or ``load_lib`` -- shows loading of libraries
    - ``debugLib`` or ``debug_lib`` -- warns about syntax errors
      when loading a library
    - ``loadProc`` or ``load_proc`` -- shows loading of procedures
      from libraries
    - ``defRes`` or ``def_res`` -- shows the names of the syzygy
      modules while converting ``resolution`` to ``list``
    - ``usage`` -- shows correct usage in error messages
    - ``Imap`` or ``imap`` -- shows the mapping of variables with
      the ``fetch`` and ``imap`` commands
    - ``notWarnSB`` or ``not_warn_sb`` -- do not warn if
      a basis is not a standard basis
    - ``contentSB`` or ``content_sb`` -- avoids to divide by the
      content of a polynomial in ``std`` and related algorithms.
      Should usually not be used.
    - ``cancelunit`` -- avoids to divide polynomials by non-constant
      units in ``std`` in the local case. Should usually not be used

    EXAMPLES::

        sage: from sage.libs.singular.option import LibSingularVerboseOptions
        sage: libsingular_verbose = LibSingularVerboseOptions()
        sage: libsingular_verbose
        verbosity options for libSingular (current value 0x00002851)
    """
    def __init__(self, **kwds):
        """
        Create a new option interface.

        EXAMPLES::

            sage: from sage.libs.singular.option import LibSingularVerboseOptions
            sage: libsingular_verb_options = LibSingularVerboseOptions()
            sage: libsingular_verb_options
            verbosity options for libSingular (current value 0x00002851)
        """
        self.global_options = &singular_verbose_options
        self.name = "verbosity"
        self.name_map = {"mem":      Sy_bit(V_SHOW_MEM),
                         "yacc":     Sy_bit(V_YACC),
                         "redefine": Sy_bit(V_REDEFINE),
                         "reading":  Sy_bit(V_READING),
                         "loadLib":  Sy_bit(V_LOAD_LIB),
                         "debugLib": Sy_bit(V_DEBUG_LIB),
                         "loadProc": Sy_bit(V_LOAD_PROC),
                         "defRes":   Sy_bit(V_DEF_RES),
                         "usage":    Sy_bit(V_SHOW_USE),
                         "Imap":     Sy_bit(V_IMAP),
                         "prompt":   Sy_bit(V_PROMPT),
                         "notWarnSB":Sy_bit(V_NSB),
                         "contentSB":Sy_bit(V_CONTENTSB),
                         "cancelunit":Sy_bit(V_CANCELUNIT),
                         }
        LibSingularOptions_abstract.__init__(self, **kwds)

    def reset_default(self):
        """
        Return to libSingular's default verbosity options.

        EXAMPLES::

            sage: from sage.libs.singular.option import opt_verb
            sage: opt_verb['not_warn_sb']
            False
            sage: opt_verb['not_warn_sb'] = True
            sage: opt_verb['not_warn_sb']
            True
            sage: opt_verb.reset_default()
            sage: opt_verb['not_warn_sb']
            False
        """
        from sage.libs.singular.singular import _saved_verbose_options
        self.global_options[0] = int(_saved_verbose_options)

cdef class LibSingularOptionsContext:
    """
    Option context

    This object localizes changes to options.

    EXAMPLES::

        sage: from sage.libs.singular.option import opt, opt_ctx
        sage: opt
        general options for libSingular (current value 0x06000082)

    ::

        sage: with opt_ctx(redTail=False):
        ....:     print(opt)
        ....:     with opt_ctx(redThrough=False):
        ....:         print(opt)
        general options for libSingular (current value 0x04000082)
        general options for libSingular (current value 0x04000002)

        sage: print(opt)
        general options for libSingular (current value 0x06000082)
    """
    cdef list bck
    cdef list bck_degBound
    cdef list bck_multBound
    cdef public LibSingularOptions_abstract opt
    cdef object options

    def __init__(self, LibSingularOptions_abstract opt, **kwds):
        """
        Create a new context.

        EXAMPLES::

            sage: from sage.libs.singular.option import LibSingularOptionsContext, opt
            sage: LibSingularOptionsContext(opt)
            general options context for libSingular
        """
        self.bck = []
        self.bck_degBound = []
        self.bck_multBound = []
        self.options = kwds
        self.opt = opt

    def __enter__(self):
        """
        EXAMPLES::

            sage: from sage.libs.singular.option import opt, opt_ctx
            sage: opt['redTail']
            True
            sage: with opt_ctx(redTail=False):
            ....:   opt['redTail']
            False
        """
        self.bck.append(self.opt.global_options[0])
        self.bck_degBound.append(Kstd1_deg)
        self.bck_multBound.append(Kstd1_mu)
        opt = self.opt.__class__()
        for k,v in self.options.iteritems():
            opt[k] = v

    def __call__(self, **kwds):
        """
        Return a new option context where ``**kwds`` are applied.

        EXAMPLES::

            sage: from sage.libs.singular.option import opt, opt_ctx
            sage: opt['redTail']
            True
            sage: with opt_ctx(redTail=False):
            ....:   opt['redTail']
            False
        """
        new = self.__class__(self.opt, **kwds)
        return new

    def __exit__(self, typ, value, tb):
        """
        EXAMPLES::

            sage: from sage.libs.singular.option import opt, opt_ctx
            sage: opt['redTail']
            True
            sage: with opt_ctx(redTail=False):
            ....:   opt['redTail']
            False
        """
        self.opt.global_options[0] = self.bck.pop()
        global Kstd1_deg
        global Kstd1_mu
        Kstd1_deg = self.bck_degBound.pop()
        Kstd1_mu  = self.bck_multBound.pop()

    def __repr__(self):
        """
        EXAMPLES::

            sage: from sage.libs.singular.option import opt_ctx
            sage: opt_ctx
            general options context for libSingular
        """
        return "%s options context for libSingular"%(self.opt.name)


opt = LibSingularOptions()
opt.reset_default()
opt_verb = LibSingularVerboseOptions()
opt_verb.reset_default()
opt_ctx = LibSingularOptionsContext(opt)
opt_verb_ctx = LibSingularOptionsContext(opt_verb)
