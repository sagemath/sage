# sage.doctest: optional - fricas
r"""
To add handling of a FriCAS type constructor, proceed as follows:

* in ``LazyParent``, add the translation of the constructor to a
  SageMath parent.

* in ``fricas.spad``, add an appropriate package, if necessary.

* in ``FRICAS_DOMAIN_DISPATCH``, map the type constructor to a function that
  constructs the package containing ``sexport`` and a function that
  constructs the element given the parsed string produced by
  ``sexport``.
"""
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import lazy_import
from sage.rings.integer import Integer
lazy_import('sage.calculus.var', ['var', 'function'])
lazy_import('sage.symbolic.expression', ['symbol_table', 'register_symbol'])
lazy_import('sage.symbolic.constants', ['I', 'e', 'pi'])

FRICAS_CONSTANTS = {'%i': I,
                    '%e': e,
                    '%pi': pi}

# the dispatch dictionary for SEXPorter and Evaluator
# if the first element of the pair is a function, it is used by
# SEXPorter to rewrite the domain
FRICAS_DOMAIN_DISPATCH = {
    "Record": ("_record", "_eval_record"),
    "Union": ("_union", "_eval_union"),
    "String": ("_inputform", "_eval_call"),
    "Integer": ("_inputform", "_eval_call"),
    "PositiveInteger": ("_inputform", "_eval_call"),
    "NonNegativeInteger": ("_inputform", "_eval_call"),
    "Float": ("_inputform", "_eval_float"),
    "Boolean": ("_inputform", "_eval_bool"),
    "AlgebraicNumber": ("_inputform", "_eval_qqbar"),
    "Expression": ("_inputform", "_eval_sr"),
    "OrderedCompletion": ("_inputform", "_eval_sr"),
    "PiDomain": ("_inputform", "_eval_sr"),
    "PrimeField": ("_finite", "_eval_gf"),
    "IntegerMod": ("_finite", "_eval_call"),
    "FiniteField": ("_finite", "_eval_gf"),
    "Fraction": ("_unary", "_eval_fraction"),
    "List": ("_unary", "_eval_list"),
    "Matrix": ("_unary", "_eval_matrix"),
    "Polynomial": ("_unary", "_eval_polynomialring"),
    "Factored": ("_unary", "_eval_factorization"),
    "Vector": ("_unary", "_eval_list"),
    "DirectProduct": ("_aggregate", "_eval_list"),
    "UnivariatePolynomial": (lambda x: ("Polynomial", x[2]),
                             "_eval_polynomialring"),  # coerce to Polynomial
    "DistributedMultivariatePolynomial": (lambda x: ("Polynomial", x[2]),
                                          "_eval_polynomialring"),
    "MultivariatePolynomial": (lambda x: ("Polynomial", x[2]),
                               "_eval_polynomialring")}


class SEXPorter:
    """
    A class constructing the FriCAS package call to export
    objects in the given (parsed) domain.
    """
    def __init__(self, domain):
        """
        INPUT:

        - ``domain`` -- a nested list of strings and integers,
          describing a FriCAS domain.

        EXAMPLES::

            sage: from sage.interfaces.fricas_translator import SEXPorter
            sage: SEXPorter(("Integer",))._domain
            ('Integer',)
        """
        if (isinstance(domain, tuple)
            and (head := domain[0]) in FRICAS_DOMAIN_DISPATCH
            and callable(fun := FRICAS_DOMAIN_DISPATCH[head][0])):
            self._domain = fun(domain)
        else:
            self._domain = domain

    def _unparse(self):
        """
        Return the FriCAS domain as a string.

        .. WARNING::

            This method does not work with arguments that are lists,
            such as `MultivariatePolynomial([x,y], Integer)`.

        EXAMPLES::

            sage: from sage.interfaces.fricas_translator import SEXPorter
            sage: SEXPorter(("Fraction", ("FiniteField", 2, 3)))._unparse()
            'Fraction(FiniteField(2,3))'

            sage: SEXPorter(('Union', ('Fraction', ('Integer',)), '"failed"'))._unparse()
            'Union(Fraction(Integer),"failed")'

            sage: dom = ('Union', (':', 'f1', ('Expression', ('Integer',))), (':', 'f2', ('List', ('Expression', ('Integer',)))), (':', 'fail', '"failed"'), (':', 'pole', '"potentialPole"'))
            sage: SEXPorter(dom)._unparse()
            'Union(f1: Expression(Integer),f2: List(Expression(Integer)),fail: "failed",pole: "potentialPole")'
        """
        if not isinstance(self._domain, tuple):
            return str(self._domain)

        if len(self._domain) == 1:
            return str(self._domain[0])

        if self._domain[0] == ':' and len(self._domain) == 3:
            tag = self._domain[1]
            dom = SEXPorter(self._domain[2])._unparse()
            return f"{tag}: {dom}"

        return (self._domain[0] + "("
                + ",".join(SEXPorter(e)._unparse() for e in self._domain[1:])
                + ")")

    def _inputform(self):
        """
        Return an ``InputFormExport`` package call for this domain.

        EXAMPLES::

            sage: from sage.interfaces.fricas_translator import SEXPorter
            sage: SEXPorter(("Integer",))._inputform()
            '(sexport$InputFormExport(Integer))'
        """
        return f"(sexport$InputFormExport({self._unparse()}))"

    def _finite(self):
        """
        Return a ``FiniteExport`` package call for this domain.

        EXAMPLES::

            sage: from sage.interfaces.fricas_translator import SEXPorter
            sage: SEXPorter(("FiniteField", 2, 3))._finite()
            '(sexport$FiniteExport(FiniteField(2,3)))'
        """
        return f"(sexport$FiniteExport({self._unparse()}))"

    def _aggregate(self):
        """
        Return an ``AggregateExport`` package call for this domain.

        EXAMPLES::

            sage: from sage.interfaces.fricas_translator import SEXPorter
            sage: SEXPorter(('DirectProduct', 2, ('Integer',)))._aggregate()
            '(sexport$AggregateExport(Integer, DirectProduct(2,Integer), (sexport$InputFormExport(Integer))))'
        """
        inner = SEXPorter(self._domain[2])
        base_str = inner._unparse()
        export_str = inner.export_call()
        return f"(sexport$AggregateExport({base_str}, {self._unparse()}, {export_str}))"

    def _unary(self):
        """
        Return a ``UnaryExport`` package call for this domain.

        EXAMPLES::

            sage: from sage.interfaces.fricas_translator import SEXPorter
            sage: SEXPorter(("List", ("Integer",)))._unary()
            '(sexport$UnaryExport(Integer, (sexport$InputFormExport(Integer))))'
        """
        inner = SEXPorter(self._domain[1])
        base_str = inner._unparse()
        export_str = inner.export_call()
        return f"(sexport$UnaryExport({base_str}, {export_str}))"

    def _record(self):
        """
        Return a function call for a record.

        This constructs a function of the form::

            sexportrecordname1name2(obj) ==
                convert([sexport(obj.name1)$Domain1Export(...),
                         sexport(obj.name2)$Domain1Export(...)])@SExpression

        EXAMPLES::

            sage: from sage.interfaces.fricas_translator import SEXPorter
            sage: dom = ('Record', (':', 'particular', ('Integer',)), (':', 'basis', ('List', ('Integer',))))
            sage: SEXPorter(dom)._record()
            'sexportRecordparticularbasis'
        """
        from sage.interfaces.fricas import fricas

        def make_call(field_name, field_domain):
            return SEXPorter(field_domain).export_call() + f"(obj.{field_name})"

        items = ",".join(make_call(field[1], field[2]) for field in self._domain[1:])
        name = "sexportRecord" + "".join(field[1] for field in self._domain[1:])
        fricas.eval(f"{name}(obj) == convert([{items}])@SExpression")
        return name

    def _union(self):
        """
        Return a function call for a union.

        This constructs a function of the form::

            sexportunionname1name2(obj) == (
                obj case Domain1 => (Domain1 sexport(obj)$Domain1Export(...));
                obj case Domain2 => (Domain2 sexport(obj)$Domain2Export(...));
                (Domain3 sexport(obj)$Domain3Export(...)))

        EXAMPLES::

            sage: from sage.interfaces.fricas_translator import SEXPorter
            sage: dom = ('Union', ('Expression', ('Integer',)), ('PrimeField', 3))
            sage: SEXPorter(dom)._union()
            'sexportUnionExpressionIntegerPrimeField3'

            sage: dom = ('Union', ('Integer',), "failed")
            sage: SEXPorter(dom)._union()
            'sexportUnionIntegerfailed'
        """
        from sage.interfaces.fricas import fricas
        from sage.misc.flatten import flatten
        from re import compile

        def is_tagged(dom):
            return isinstance(dom, tuple) and len(dom) == 3 and dom[0] == ':'

        # detect if this is a tagged union (all fields tagged)
        tagged = all(is_tagged(dom) for dom in self._domain[1:])

        def make_call(field_domain, tag_index):
            if isinstance(field_domain, str):  # e.g., '"failed"'
                s = f'convert({field_domain})@SExpression'
            else:
                domain = SEXPorter(field_domain)
                s = f'{domain.export_call()}(obj)'

            return f'convert([convert({tag_index})@SExpression, {s}])@SExpression'

        def make_case(field, tag_index):
            if tagged:
                # field = (':', tag_name, domain)
                tag_name = field[1]
                field_domain = field[2]
            else:
                tag_name = None
                field_domain = field

            # for the last case we want no "case" prefix
            if tag_index == len(self._domain) - 2:
                return make_call(field_domain, tag_index)

            if tagged:
                return f'obj case {tag_name} => {make_call(field_domain, tag_index)}'

            domain = SEXPorter(field_domain)
            return f'obj case {domain._unparse()} => {make_call(field_domain, tag_index)}'

        items = "; ".join(make_case(dom, tag)
                          for tag, dom in enumerate(self._domain[1:]))

        name = "sexport" + "".join(map(str, flatten(self._domain)))
        pattern = compile(r'[\W_]+')
        name = pattern.sub('', name)

        fricas.eval(f"{name}(obj:{self._unparse()}): SExpression == ({items});")
        return name


    def export_call(self):
        """
        Return the function call doing the export, as a string.

        EXAMPLES::

            sage: from sage.interfaces.fricas_translator import SEXPorter
            sage: SEXPorter(("List", ("FiniteField", 2, 3))).export_call()
            '(sexport$UnaryExport(FiniteField(2,3), (sexport$FiniteExport(FiniteField(2,3)))))'

            sage: SEXPorter(("List", ("UnivariatePolynomial", "x", ("Integer",)))).export_call()
            '(sexport$UnaryExport(Polynomial(Integer), (sexport$UnaryExport(Integer, (sexport$InputFormExport(Integer))))))'

            sage: SEXPorter(('DirectProduct', 2, ('Integer',))).export_call()
            '(sexport$AggregateExport(Integer, DirectProduct(2,Integer), (sexport$InputFormExport(Integer))))'
        """
        if not isinstance(self._domain, tuple):
            # it should be a FriCAS string
            return self._domain

        head = self._domain[0]
        if head not in FRICAS_DOMAIN_DISPATCH:
            raise NotImplementedError(f"{head} cannot be translated from FriCAS to SageMath yet")

        return getattr(self, FRICAS_DOMAIN_DISPATCH[head][0])()


class SEXEvaluator:
    def __init__(self, ast, dom):
        """
        INPUT:

        - ``ast`` -- a nested list of strings and integers describing a FriCAS object
        - ``dom`` -- a :class:`LazyParent` describing the domain of ``ast``

        EXAMPLES::

            sage: from sage.interfaces.fricas_translator import SEXEvaluator, LazyParent
            sage: SEXEvaluator(1, LazyParent(['Integer']))._ast
            1
        """
        self._ast = ast
        self._dom = dom

    def eval(self):
        r"""
        Return the evaluation of ``self`` in the given parent.

        EXAMPLES::

            sage: from sage.interfaces.fricas_translator import SEXEvaluator, LazyParent
            sage: SEXEvaluator("5", LazyParent(['Integer'])).eval()
            5
        """
        head = self._dom.head()
        return getattr(self, FRICAS_DOMAIN_DISPATCH[head][1])()

    # specific evaluators
    # the naming convention is the parent in sage, in lower case

    def _eval_call(self):
        r"""
        Return the evaluation by passing it to the parent.

        EXAMPLES::

            sage: from sage.interfaces.fricas_translator import LazyParent, SEXEvaluator
            sage: SEXEvaluator("-5", LazyParent(['Integer']))._eval_call()
            -5
        """
        return self._dom.parent()(self._ast)

    def _eval_record(self):
        r"""
        Return the evaluation as a record.

        EXAMPLES::

            sage: from sage.interfaces.fricas_translator import LazyParent, SEXEvaluator
            sage: fricas._start()  # to fill the symbol table
            sage: obj = (1, (('exp', ('*', -1, 'x')),))
            sage: EXPR = ('Expression', ('Integer',))
            sage: dom = ('Record', (':', 'particular', EXPR), (':', 'basis', ('List', EXPR)))
            sage: result = SEXEvaluator(obj, LazyParent(dom)).eval()
            sage: result["particular"]
            1
            sage: result["basis"]
            [e^(-x)]
        """
        assert len(self._dom._domain) == len(self._ast) + 1
        result = dict()
        for field, value in zip(self._dom._domain[1:], self._ast):
            key = field[1]
            parent = LazyParent(field[2])
            val = SEXEvaluator(value, parent).eval()
            result[key] = val
        return result

    def _eval_union(self):
        r"""
        Return the evaluation of a union.

        This has to dispatch to the correct type

        EXAMPLES::

            sage: from sage.interfaces.fricas_translator import LazyParent, SEXEvaluator
            sage: obj = (1, (('exp', ('*', -1, 'x')),))
            sage: EXPR = ('Expression', ('Integer',))
            sage: dom = ('Record', (':', 'particular', EXPR), (':', 'basis', ('List', EXPR)))
            sage: result = SEXEvaluator(obj, LazyParent(dom)).eval()
            sage: result["particular"]
            1
            sage: result["basis"]
            [e^(-x)]
        """
        tag, val = self._ast
        dom = self._dom._domain[1 + tag]
        if isinstance(dom, tuple) and dom[0] == ':':
            dom = dom[2]
        if not isinstance(dom, tuple):
            assert dom == val
            return val[1:-1]
        parent = LazyParent(dom)
        return SEXEvaluator(val, parent).eval()

    def _eval_bool(self):
        r"""
        Return the evaluation as a boolean.

        EXAMPLES::

            sage: from sage.interfaces.fricas_translator import LazyParent, SEXEvaluator
            sage: SEXEvaluator("true", LazyParent(['Boolean']))._eval_bool()
            True
        """
        return self._ast == "true"

    def _eval_float(self):
        r"""
        Return the evaluation as an arbitrary precision float.

        EXAMPLES::

            sage: from sage.interfaces.fricas_translator import SEXParser, SEXPorter, SEXEvaluator, LazyParent
            sage: f = fricas("3.1415::Float")
            sage: P = f._check_valid()
            sage: dom_str = P.get_string(f"sageprint(dom({f._name})::Any)")
            sage: dom = SEXParser(dom_str).parse()
            sage: fun = SEXPorter(dom).export_call()
            sage: obj_str = P.get_string(f"sageprint({fun}({f._name}))")
            sage: obj = SEXParser(obj_str).parse(); obj
            ('float', 231801786030234225607, -66, 2)
            sage: SEXEvaluator(obj, LazyParent(dom))._eval_float()
            3.1415000000000000000
        """
        _, x, e, b = self._ast
        # Warning: precision$Float gives the current precision,
        # whereas the bitlength of two's complement gives the
        # precision of self.
        x = int(x)
        prec = max(x.bit_length() if x >= 0 else (~x).bit_length(), 53)
        P = self._dom.parent(prec=prec)
        return P(Integer(x) * Integer(b) ** Integer(e))

    def _eval_gf(self):
        r"""
        Return the evaluation as an element of a finite field.

        EXAMPLES::

            sage: from sage.interfaces.fricas_translator import SEXEvaluator, LazyParent
            sage: SEXEvaluator(5, LazyParent(['FiniteField', 2, 3]))._eval_gf()
            z3^2 + 1
        """
        P = self._dom.parent()
        return P.from_integer(self._ast % P.order())

    def _eval_list(self):
        r"""
        Return the evaluation by coercing a list to the parent.

        EXAMPLES::

            sage: from sage.interfaces.fricas_translator import SEXParser, LazyParent, SEXEvaluator
            sage: ast = SEXParser("(1 2 3)").parse()
            sage: dom = LazyParent(('List', ('Integer',)))
            sage: SEXEvaluator(ast, dom)._eval_list()
            [1, 2, 3]
        """
        base = self._dom.base()
        P = self._dom.parent()
        return P([SEXEvaluator(a, base).eval() for a in self._ast])

    def _eval_fraction(self):
        r"""
        Return the evaluation as a fraction.

        EXAMPLES::

            sage: from sage.interfaces.fricas_translator import SEXParser, LazyParent, SEXEvaluator
            sage: ast = SEXParser("(3 2)").parse()
            sage: SEXEvaluator(ast, LazyParent(('Fraction', ('Integer',))))._eval_fraction()
            3/2
        """
        base = self._dom.base()
        if not isinstance(self._ast, tuple):
            num = SEXEvaluator(self._ast, base).eval()
            return dom.parent()(num)

        a, b = self._ast
        num = SEXEvaluator(a, base).eval()
        den = SEXEvaluator(b, base).eval()
        return num / den

    def _eval_matrix(self):
        r"""
        Return the evaluation of ``ast`` as a matrix.

        EXAMPLES::

            sage: from sage.interfaces.fricas_translator import SEXParser, LazyParent, SEXEvaluator
            sage: ast = SEXParser("((1 2) (3 4))").parse()
            sage: SEXEvaluator(ast, LazyParent(('Matrix', ('Integer',))))._eval_matrix()
            [1 2]
            [3 4]

            sage: SEXEvaluator(ast, LazyParent(('Matrix', ('PrimeField', 3))))._eval_matrix()
            [1 2]
            [0 1]

            sage: ast = SEXParser("((1 2) (3 (sin x)))").parse()
            sage: SEXEvaluator(ast, LazyParent(('Matrix', ('Expression', ('Integer',)))))._eval_matrix()
            [     1      2]
            [     3 sin(x)]

            sage: ast = SEXParser("(((float 191846138366579336806 -67 2)))").parse()
            sage: SEXEvaluator(ast, LazyParent(('Matrix', ('Float',))))._eval_matrix()
            [1.3000000000000000000]
        """
        base = self._dom.base()
        m = [[SEXEvaluator(e, base).eval() for e in row]
             for row in self._ast]
        P = self._dom.parent()
        return P(m)

    def _eval_polynomialring(self):
        r"""
        Return the evaluation as a polynomial.

        EXAMPLES::

            sage: from sage.interfaces.fricas_translator import SEXParser, SEXPorter, SEXEvaluator, LazyParent
            sage: f = fricas("x^2*y - 3*z + 1")
            sage: P = f._check_valid()
            sage: dom_str = P.get_string(f"sageprint(dom({f._name})::Any)")
            sage: dom = SEXParser(dom_str).parse()
            sage: fun = SEXPorter(dom).export_call()
            sage: obj_str = P.get_string(f"sageprint({fun}({f._name}))")
            sage: obj = SEXParser(obj_str).parse(); obj
            (((('z', 1),), -3), ((('y', 1), ('x', 2)), 1), ((), 1))
            sage: SEXEvaluator(obj, LazyParent(dom))._eval_polynomialring()
            x^2*y - 3*z + 1
        """
        base = self._dom.base()
        names = tuple(sorted(set(v for mon, _ in self._ast for v, _ in mon)))
        P = self._dom.parent(names=names)

        from sage.rings.polynomial.polynomial_ring import PolynomialRing_generic
        if isinstance(P, PolynomialRing_generic):

            def to_exponent(mon):
                if len(mon):
                    return mon[0][1]
                return 0

            return P._from_dict({to_exponent(mon): SEXEvaluator(c, base).eval()
                                 for mon, c in self._ast})

        def to_tuple(mon):
            t = [0]*len(P._names)
            for v, e in mon:
                t[P._names.index(v)] = e
            return tuple(t)

        return P._from_dict({to_tuple(mon): SEXEvaluator(c, base).eval()
                             for mon, c in self._ast})

    def _eval_factorization(self):
        r"""
        Return the evaluation as a factored object.

        EXAMPLES::

            sage: from sage.interfaces.fricas_translator import SEXParser, SEXPorter, SEXEvaluator, LazyParent
            sage: f = fricas("-48").factor()
            sage: P = f._check_valid()
            sage: dom_str = P.get_string(f"sageprint(dom({f._name})::Any)")
            sage: dom = SEXParser(dom_str).parse()
            sage: fun = SEXPorter(dom).export_call()
            sage: obj_str = P.get_string(f"sageprint({fun}({f._name}))")
            sage: obj = SEXParser(obj_str).parse(); obj
            (-1, ((2, 4), (3, 1)))
            sage: SEXEvaluator(obj, LazyParent(dom))._eval_factorization()
            -1 * 2^4 * 3
        """
        from sage.structure.factorization import Factorization
        base = self._dom.base()
        unit, factors = self._ast
        return Factorization([(SEXEvaluator(f, base).eval(), e)
                              for f, e in factors],
                             unit=SEXEvaluator(unit, base).eval(),
                             sort=False,
                             simplify=False)

    def _eval_sr_aux(self):
        r"""
        Return the evaluation as an element of the
        symbolic ring.

        Post processing by :meth:`_eval_sr` is necessary.

        EXAMPLES:

        This function cannot use the symbol table to translate
        symbols which are not function calls, as :issue:`31849` shows
        ``D`` would erroneously be interpreted as differential
        then::

            sage: var("D")
            D
            sage: integrate(D/x, x, algorithm='fricas')
            D*log(x)

        However, it does have to check for constants, for example
        ``%pi``::

            sage: from sage.interfaces.fricas_translator import SEXEvaluator
            sage: SEXEvaluator("%pi", None)._eval_sr_aux()
            pi
        """
        from sage.symbolic.ring import SR

        if isinstance(self._ast, (Integer, float)):
            return SR(self._ast)

        if isinstance(self._ast, str):
            try:
                return FRICAS_CONSTANTS[self._ast]
            except KeyError:
                return var(self._ast.replace("%", "_"))

        f, *args = self._ast
        try:
            fun = symbol_table["fricas"][(f, len(args))]
        except KeyError:
            try:
                fun = symbol_table["fricas"][(f, -1)]
            except KeyError:
                fun = function(f)

        args_eval = [SEXEvaluator(a, self._dom)._eval_sr_aux() for a in args]
        return fun(*args_eval)

    def _eval_qqbar(self):
        r"""
        Evaluate as element of the algebraic field.

        EXAMPLES::

            sage: from sage.interfaces.fricas_translator import SEXParser, SEXPorter, SEXEvaluator, LazyParent
            sage: f = fricas("rootOf(x^5 - x - 1)$AN")
            sage: P = f._check_valid()
            sage: dom_str = P.get_string(f"sageprint(dom({f._name})::Any)")
            sage: dom = SEXParser(dom_str).parse()
            sage: fun = SEXPorter(dom).export_call()
            sage: obj_str = P.get_string(f"sageprint({fun}({f._name}))")
            sage: obj = SEXParser(obj_str).parse(); obj
            ('::',
             ('rootOf', ('+', ('+', ('^', 'x', 5), ('*', -1, 'x')), -1), 'x'),
             ('AlgebraicNumber',))
            sage: SEXEvaluator(obj, LazyParent(dom))._eval_qqbar()
            1.167303978261419?
        """
        ex = SEXEvaluator(self._ast[1], None)._eval_sr()
        return self._dom.parent()(ex)

    def _eval_sr(self):
        r"""
        Evaluate as an element of the symbolic ring.


        EXAMPLES::

            sage: from sage.interfaces.fricas_translator import SEXEvaluator
            sage: SEXEvaluator("x", None)._eval_sr()
            x
        """
        # a FriCAS expressions may contain implicit references to a
        # rootOf expression within itself, as for example in the
        # result of integrate(1/(1+x^5), x).  Each algebraic number
        # appearing in the expression is only introduced once and
        # assigned a variable (usually of the form %%...).
        rootOf = dict()  # (variable, polynomial)
        rootOf_ev = dict()  # variable -> (complex) algebraic number

        def convert_rootOf(x, y):
            if y in rootOf:
                assert rootOf[y] == x
            else:
                rootOf[y] = x
            return y

        register_symbol(convert_rootOf, {'fricas': 'rootOf'}, 2)
        ex = self._eval_sr_aux()

        # postprocessing of rootOf
        from sage.rings.qqbar import QQbar
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        while rootOf:
            for var, poly in rootOf.items():
                pvars = poly.variables()
                rvars = [v for v in pvars if v not in rootOf_ev]  # remaining variables
                uvars = [v for v in rvars if v in rootOf]  # variables to evaluate
                if len(uvars) == 1:
                    assert uvars[0] == var, "the only variable in uvars should be %s but is %s" % (var, uvars[0])
                    break
            else:
                assert False, "circular dependency in rootOf expression"
            # substitute known roots
            poly = poly.subs(rootOf_ev)
            evars = set(poly.variables()).difference([var])
            del rootOf[var]
            if evars:
                # we just need any root per FriCAS specification -
                # however, if there are extra variables, we cannot
                # use QQbar.any_root
                rootOf_ev[var] = poly.roots(var, multiplicities=False)[0]
            else:
                R = PolynomialRing(QQbar, "x")
                # PolynomialRing does not accept variable names with
                # leading underscores
                poly = R(poly.subs({var: R.gen()}))
                # we just need any root per FriCAS specification
                rootOf_ev[var] = poly.any_root()

        return ex.subs({var: (val.radical_expression()
                              if val.parent() is QQbar else val)
                        for var, val in rootOf_ev.items()})


class LazyParent:
    def __init__(self, domain):
        """
        Initialize this lazy parent from a FriCAS domain description.

        EXAMPLES::

            sage: from sage.interfaces.fricas_translator import LazyParent
            sage: LazyParent(['Integer']).head()
            'Integer'
        """
        self._domain = domain

    def head(self):
        """
        Return the constructor.

        EXAMPLES::

            sage: from sage.interfaces.fricas_translator import LazyParent
            sage: LazyParent(('List', ('Integer',))).head()
            'List'
        """
        return self._domain[0]

    def args(self):
        """
        Return the arguments of the constructor.

        EXAMPLES::

            sage: from sage.interfaces.fricas_translator import LazyParent
            sage: LazyParent(('List', ('Integer',))).args()
            (('Integer',),)
        """
        return self._domain[1:]

    @cached_method
    def base(self):
        """
        Return the base domain as a :class:`LazyParent`.

        EXAMPLES::

            sage: from sage.interfaces.fricas_translator import LazyParent
            sage: LazyParent(('List', ('Integer',))).base().head()
            'Integer'
        """
        head = self.head()
        args = self.args()
        if head in ["List",
                    "Fraction",
                    "Matrix",
                    "Polynomial",
                    "Factored",
                    "Vector"]:
            return LazyParent(args[0])

        if head in ["UnivariatePolynomial",
                    "MultivariatePolynomial",
                    "DistributedMultivariatePolynomial",
                    "DirectProduct"]:
            return LazyParent(args[1])

        raise NotImplementedError(f"Cannot build base for {head}")

    @cached_method
    def parent(self, **kwargs):
        """
        Return the corresponding SageMath parent.

        EXAMPLES::

            sage: from sage.interfaces.fricas_translator import LazyParent
            sage: LazyParent(['Integer']).parent()
            Integer Ring
        """
        head = self.head()
        args = self.args()

        if head == "String":
            return lambda x: x[1:-1]

        if head == "List":
            return list

        if head in ["Vector", "DirectProduct"]:
            from sage.modules.free_module_element import vector
            return vector

        if head in ["Matrix", "RectangularMatrix", "SquareMatrix"]:
            from sage.matrix.constructor import matrix
            return matrix

        if head in ["FiniteField", "PrimeField"]:
            from sage.rings.finite_rings.finite_field_constructor import GF
            return GF(*args)

        if head == "IntegerMod":
            from sage.rings.finite_rings.integer_mod_ring import IntegerModRing
            return IntegerModRing(*args)

        if head in ["Integer", "PositiveInteger", "NonNegativeInteger"]:
            from sage.rings.integer_ring import ZZ
            return ZZ

        if head == "AlgebraicNumber":
            from sage.rings.qqbar import QQbar
            return QQbar

        if head == "Expression":
            from sage.symbolic.ring import SR
            return SR

        if head == "Float":
            from sage.rings.real_mpfr import RealField
            return RealField(kwargs["prec"])

        if head == "Fraction":
            from sage.rings.fraction_field import FractionField
            base = self.base().parent(**kwargs)
            return FractionField(base)

        if head in ["UnivariatePolynomial",
                    "MultivariatePolynomial",
                    "DistributedMultivariatePolynomial"]:
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            base = self.base().parent(**kwargs)
            return PolynomialRing(base, args[0])

        if head == "Polynomial":
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            base = self.base().parent(**kwargs)
            if not hasattr(self, "_polynomial_symbols"):
                self._polynomial_symbols = []
            self._polynomial_symbols = sorted(set(list(kwargs.get("names", []))
                                                  + self._polynomial_symbols))
            # we always want a multivariate polynomial ring here
            return PolynomialRing(base,
                                  len(self._polynomial_symbols),
                                  names=self._polynomial_symbols)

        raise NotImplementedError(f"Cannot build parent for {head}")


class SEXParser:
    _WHITESPACE = " "
    _LEFTBRACKET = "("
    _RIGHTBRACKET = ")"
    _STRINGMARKER = '"'
    _ESCAPEMARKER = '_'  # STRINGMARKER must be escaped in strings

    def __init__(self, s):
        """
        INPUT:

        - ``s`` -- string

        EXAMPLES::

            sage: from sage.interfaces.fricas_translator import SEXParser
            sage: SEXParser("abc")._s
            'abc'
        """
        self._s = s
        self._start = 0  # specifies where to start parsing

    def parse(self):
        r"""
        Parse the string.

        TESTS::

            sage: from sage.interfaces.fricas_translator import SEXParser, SEXPorter
            sage: SEXParser("abc").parse()
            'abc'

            sage: SEXParser("(asin c)").parse()
            ('asin', 'c')

            sage: SEXParser("(pi)").parse()
            ('pi',)

            sage: f = fricas('(x + y/2)::DMP(["x", "y", "z"], FRAC INT)')
            sage: P = f._check_valid()
            sage: dom_str = P.get_string(f"sageprint(dom({f._name})::Any)")
            sage: dom = SEXParser(dom_str).parse(); dom
            ('DistributedMultivariatePolynomial',
             ('x', 'y', 'z'),
             ('Fraction', ('Integer',)))
        """
        a = self._start
        while self._s[a] in self._WHITESPACE:
            a += 1

        self._start = a
        if self._s[a] == self._LEFTBRACKET:
            return self._parse_list()
        if self._s[a] == self._STRINGMARKER:
            return self._parse_string()
        return self._parse_atom()

    def _parse_list(self):
        """
        Parse the initial part of a string, assuming that it is a
        whitespace separated list.

        TESTS::

            sage: from sage.interfaces.fricas_translator import SEXParser
            sage: SEXParser("()").parse()  # indirect doctest
            ()

            sage: SEXParser("(a b c)").parse()
            ('a', 'b', 'c')

            sage: SEXParser('(bcd)').parse()
            ('bcd',)
        """
        # self._start specifies the position of the left bracket
        assert self._s[self._start] == self._LEFTBRACKET
        self._start += 1

        args = []
        while self._s[self._start] != self._RIGHTBRACKET:
            e = self.parse()
            args.append(e)
            self._start += 1

        return tuple(args)

    def _parse_atom(self):
        """
        Parse the initial part of a string, assuming that it is an
        atom, but not a string.

        Symbols and numbers must not contain ``_WHITESPACE`` and
        ``_RIGHTBRACKET``.

        TESTS::

            sage: from sage.interfaces.fricas_translator import SEXParser
            sage: SEXParser("abc").parse()  # indirect doctest
            'abc'
            sage: SEXParser("123 xyz").parse()
            123
        """
        a = self._start
        b = len(self._s)

        while (a < b
               and self._s[a] not in self._WHITESPACE
               and self._s[a] != self._RIGHTBRACKET):
            a += 1

        token = self._s[self._start:a]
        self._start = a - 1
        try:
            return Integer(token)
        except TypeError:
            try:
                return float(token)
            except ValueError:
                return token

    def _parse_string(self):
        r"""
        Parse the initial part of a string, assuming that it represents a
        string.

        TESTS::

            sage: from sage.interfaces.fricas_translator import SEXParser
            sage: SEXParser('"abc" 123')._parse_string()
            '"abc"'

            sage: SEXParser('"" 123')._parse_string()
            '""'

            sage: SEXParser('"____" 123')._parse_string()
            '"__"'

            sage: SEXParser('"_a" 123')._parse_string()
            '"a"'

            sage: SEXParser('"_" _"" 123')._parse_string()
            '"" ""'

            sage: SEXParser('"(b c)"')._parse_string()
            '"(b c)"'
        """
        a = self._start  # the position of the left quote
        assert self._s[a] == self._STRINGMARKER
        b = a + 1
        a = b
        S = []

        while self._s[b] != self._STRINGMARKER:
            if self._s[b] == self._ESCAPEMARKER:
                S.append(self._s[a:b])
                b += 1
                a = b
            b += 1

        self._start = b
        S.append(self._s[a:b])
        return '"' + "".join(S) + '"'
