# sage_setup: distribution = sagemath-brial
# distutils: language = c++
# distutils: extra_compile_args = -std=c++11

from libcpp.string cimport string as std_string
from libcpp.vector cimport vector


cdef extern from "polybori/pbori_defs.h" namespace "COrderEnums":
    cdef enum ordercodes:
        pblp            "CTypes::lp"
        pbdlex          "CTypes::dlex"
        pbdp_asc        "CTypes::dp_asc"
        pbblock_dlex    "CTypes::block_dlex"
        pbblock_dp_asc  "CTypes::block_dp_asc"
        pbdp            "static_cast<COrderEnums::ordercodes>(17)"
        pbblock_dp      "static_cast<COrderEnums::ordercodes>(19)"

    cdef enum comparecodes:
        less_than               "CTypes::less_than"
        equality                "CTypes::equality"
        greater_than            "CTypes::greater_than"
        less_or_equal_max       "CTypes::less_or_equal_max"
        greater_or_equal_min    "CTypes::greater_or_equal_min"


cdef extern from "sage/libs/polybori/pb_wrap.h":
    cdef cppclass PBNavigator "CCuddNavigator":
        PBNavigator thenBranch()
        PBNavigator elseBranch()
        int value "operator*"()
        bint isConstant()
        bint isTerminated()
        bint operator==(const PBNavigator&)
        size_t hash()

    cdef cppclass PBBlockIter "COrderingBase::block_iterator":
        long dereference "operator*"()
        void increment "operator++"()
        bint operator!=(const PBBlockIter&)

    cdef cppclass PBOrdering "COrderingBase":
        int getOrderCode()
        int getBaseOrderCode()
        void appendBlock(int nextBlockStartIndex)
        PBBlockIter blockBegin()
        PBBlockIter blockEnd()
        int lastBlockStart()
        bint isDegreeOrder()

    cdef cppclass PBSet "DefaultRinged<BooleSet> "
    cdef cppclass PBPoly "DefaultRinged<BoolePolynomial> "

    cdef cppclass PBBooleVariable "BooleVariable":
        PBBooleVariable(int idx, const PBRing& r)

    cdef cppclass PBVar "DefaultRinged<BooleVariable> ":
        void operator=(const PBBooleVariable&)
        int index()
        bint operator==(const PBVar&)

    cdef cppclass PBRing "BoolePolyRing":
        PBRing()
        PBRing(int nvars, int order)
        size_t id()
        size_t nVariables()
        PBVar variable(int n)
        PBOrdering ordering()
        void changeOrdering(int)
        size_t hash()
        void setVariableName(int idx, char *varname)
        char* getVariableName(int idx)
        PBRing clone()
        PBPoly coerce(PBPoly rhs)
        PBPoly zero()
        PBPoly one()

    cdef cppclass PBVarBlock "VariableBlock":
        PBVarBlock(int size, int start_index, int offset, bint reverse,
                   PBRing r)
        PBVar operator()(int index)

    cdef cppclass PBMonomIter "BooleMonomial::const_iterator":
        int dereference()
        void increment()
        bint operator==(const PBMonomIter&)

    cdef cppclass PBMonomVarIter "BooleMonomial::variable_iterator":
        PBVar dereference()
        void increment()
        bint operator==(const PBMonomVarIter&)

    cdef cppclass PBMonom "DefaultRinged<BooleMonomial> ":
        PBMonom()
        PBMonom(PBRing r)
        PBMonom(PBVar m)
        void operator=(const PBBooleSet&)
        bint reducibleBy(PBMonom rhs)
        int deg()
        size_t hash()
        size_t stableHash()
        int firstIndex()
        comparecodes compare(PBMonom rhs)
        PBSet set()
        PBSet divisors()
        PBSet multiples(PBMonom rhs)
        PBMonomIter begin()
        PBMonomIter end()
        PBMonomVarIter variableBegin()
        PBMonomVarIter variableEnd()
        void imul "operator*="(PBMonom right)
        void idiv "operator/="(PBMonom right)
        PBNavigator navigation()
        PBMonom GCD(PBMonom rhs)
        PBRing ring()

    cdef cppclass PBSetIter "BooleSet::const_iterator":
        PBSetIter(const PBSetIter&)
        PBMonom dereference()
        void increment()
        bint operator==(const PBSetIter&)

    cdef cppclass PBBooleSet "BooleSet":
        PBBooleSet()
        PBBooleSet(PBRing r)
        PBBooleSet(PBPoly p)
        PBBooleSet(PBNavigator d, PBRing r)
        PBBooleSet(int ind, PBNavigator a, PBNavigator b, PBRing r)

    cdef cppclass PBSet "DefaultRinged<BooleSet> ":
        void operator=(const PBBooleSet&)
        bint owns(PBMonom val)
        int nNodes()
        int nSupport()
        int size()
        size_t hash()
        size_t stableHash()
        PBNavigator navigation()
        PBSet cartesianProduct(PBSet rhs)
        PBSet diff(PBSet rhs)
        PBSet divide(PBMonom rhs)
        PBSet change(int idx)
        PBSet subset0(int i)
        PBSet subset1(int i)
        PBSet unite(PBSet rhs)
        PBSet intersect(PBSet rhs)
        PBMonom usedVariables()
        PBSet divisorsOf(PBMonom rhs)
        PBSet multiplesOf(PBMonom rhs)
        double sizeDouble()
        PBSetIter begin()
        PBSetIter end()
        bint isZero()
        bint isOne()

    PBSet pb_include_divisors "include_divisors"(PBSet p)
    PBSet pb_minimal_elements "minimal_elements"(PBSet p)

    cdef cppclass PBPolyIter "BoolePolynomial::ordered_iterator":
        PBPolyIter(const PBPolyIter&)
        PBMonom dereference()
        void increment()
        bint operator==(const PBPolyIter&)

    cdef cppclass PBBoolePolynomial "BoolePolynomial":
        void operator=(const PBPoly&)
        PBBoolePolynomial(PBRing r)
        PBBoolePolynomial(PBSet d)
        PBBoolePolynomial(PBMonom d)
        PBBoolePolynomial(PBVar d)
        PBBoolePolynomial(int d, PBRing r)

    cdef cppclass PBPoly "DefaultRinged<BoolePolynomial> ":
        PBPoly()
        void operator=(const PBBoolePolynomial&)
        int deg()
        int leadDeg()
        int lexLeadDeg()
        int totalDeg()
        int length()
        int eliminationLength()
        int nUsedVariables()
        int nNodes()
        size_t hash()
        size_t stableHash()
        bint isZero()
        bint isOne()
        bint isConstant()
        bint isSingleton()
        bint isPair()
        bint isSingletonOrPair()
        bint hasConstantPart()
        bint firstReducibleBy(PBPoly rhs)
        PBMonom lead()
        PBMonom lexLead()
        PBMonom firstTerm()
        PBMonom usedVariables()
        PBPoly gradedPart(int d)
        PBSet diagram()
        PBSet set()
        PBSet leadDivisors()
        PBNavigator navigation()
        PBPolyIter orderedBegin()
        PBPolyIter orderedEnd()
        void iadd "operator+="(PBPoly right)
        void iadd_PBMonom "operator+="(PBMonom right)
        void imul "operator*="(PBPoly right)
        void imul_monom "operator*="(PBMonom right)
        void idiv "operator/="(PBPoly right)
        void idiv_monom "operator/="(PBMonom right)
        bint operator==(const PBPoly&)

    PBSet pb_zeros "zeros"(PBPoly p, PBSet s)
    PBPoly pb_spoly "spoly"(PBPoly p, PBPoly r)
    PBPoly pb_map_every_x_to_x_plus_one "map_every_x_to_x_plus_one"(PBPoly)


ctypedef vector[PBBoolePolynomial] PBPolyVector
ctypedef vector[PBBoolePolynomial].iterator PBPolyVectorIter


cdef extern from *:
    cdef cppclass PBPolyEntry "PolyEntry":
        PBPoly p

    cdef cppclass PBRedStrategy "ReductionStrategy":
        PBRedStrategy()
        PBRedStrategy(const PBRing&)
        PBSet leadingTerms
        PBSet minimalLeadingTerms
        PBSet leadingTerms11
        PBSet leadingTerms00
        PBSet llReductor
        PBSet monomials
        PBSet monomials_plus_one
        bint optBrutalReductions
        bint optLL
        PBPoly nf(PBPoly p)
        bint optRedTailDegGrowth
        bint optRedTail
        void setupSetsForLastElement()
        bint canRewrite(PBPoly p)
        void addGenerator(PBPoly p)
        int select1(PBMonom p)

        int select_short(PBPoly p)
        PBPoly headNormalForm(PBPoly p)
        PBPoly reducedNormalForm(PBPoly p)

        size_t size()
        PBPolyEntry operator[](int)

    cdef cppclass PBFGLMStrategy "FGLMStrategy":
        PBFGLMStrategy(PBRing from_ring, PBRing to_ring, PBPolyVector vec)
        PBPolyVector main()

    cdef cppclass PBGBStrategy "GroebnerStrategy":
        PBGBStrategy(PBRing)
        bint reduceByTailReduced
        PBRedStrategy generators

        bint enabledLog
        unsigned int reductionSteps
        int normalForms
        int currentDegree
        int chainCriterions
        int variableChainCriterions
        int easyProductCriterions
        int extendedProductCriterions
        int averageLength

        bint optHFE
        bint optLazy

        bint optDelayNonMinimals

        bint optExchange
        bint optAllowRecursion

        bint optStepBounded
        bint optLinearAlgebraInLastBlock
        bint optModifiedLinearAlgebra
        bint optRedTailInLastBlock
        bint optDrawMatrices
        std_string matrixPrefix

        PBSet llReductor
        PBSet minimalLeadingTerms
        PBSet leadingTerms
        bint containsOne()
        bint variableHasValue(int idx)
        int addGenerator(PBPoly, bint is_impl)
        void addGeneratorDelayed(PBPoly)
        void addAsYouWish(PBPoly)
        void addNonTrivialImplicationsDelayed(const PBPolyEntry&)
        void symmGB_F2()
        void llReduceAll()
        void cleanTopByChainCriterion "pairs.cleanTopByChainCriterion"()
        size_t npairs "pairs.queue.size"()
        int suggestPluginVariable()
        PBPoly nextSpoly()
        PBPoly nf(PBPoly p)
        PBPolyVector allGenerators()
        PBPolyVector minimalize()
        PBPolyVector minimalizeAndTailReduce()
        PBPolyVector faugereStepDense(PBPolyVector v)

    int pairs_top_sugar(PBGBStrategy strat)
    PBPolyVector someNextDegreeSpolys(PBGBStrategy strat, size_t n)
    PBPolyVector nextDegreeSpolys(PBGBStrategy strat)
    PBPolyVector small_next_degree_spolys(PBGBStrategy strat, double f, size_t n)
    void implications(PBGBStrategy strat, int i)

    PBPoly cheap_reductions(PBRedStrategy strat, PBPoly p)

    PBPoly pb_add_up_polynomials "add_up_polynomials" \
        (PBPolyVector v, PBPoly p)

    PBPoly pb_nf3 "nf3" (PBRedStrategy strat, PBPoly p, PBMonom m)

    PBPoly pb_red_tail "red_tail" (PBRedStrategy strat, PBPoly p)

    PBPoly pb_interpolate "interpolate" (PBSet z, PBSet o)

    PBPoly pb_interpolate_smallest_lex "interpolate_smallest_lex" \
        (PBSet z, PBSet o)

    PBSet pb_contained_variables_cudd_style "contained_variables_cudd_style" \
        (PBSet m)

    PBSet pb_mod_var_set "mod_var_set" (PBSet a, PBSet v)

    PBPoly pb_mult_fast_sim "mult_fast_sim" (PBPolyVector v, PBRing r)

    PBSet pb_recursively_insert "recursively_insert"(PBNavigator p,
                                                int idx, PBSet m)
    PBPoly pb_ll_red_nf_noredsb "ll_red_nf_noredsb"(PBPoly p,
                                                PBSet reductors)
    PBPoly pb_ll_red_nf_noredsb_single_recursive_call "ll_red_nf_noredsb_single_recursive_call" \
        (PBPoly p, PBSet reductors)

    PBPoly pb_ll_red_nf "ll_red_nf"(PBPoly p, PBSet reductors)

    PBSet pb_mod_mon_set "mod_mon_set"(PBSet as, PBSet vs)

    PBPolyVector pb_parallel_reduce "parallel_reduce" \
        (PBPolyVector inp, PBGBStrategy strat, int average_steps, double delay_f)

    PBPolyVector pb_gauss_on_polys "gauss_on_polys" (PBPolyVector inp)

    PBPolyVector pb_easy_linear_factors "easy_linear_factors"(PBPoly p)

    PBPoly pb_substitute_variables "substitute_variables<BoolePolyRing, std::vector<BoolePolynomial>, BoolePolynomial>" \
        (PBRing ring, PBPolyVector vec, PBPoly poly) except +

    void pb_set_variable_name "BooleEnv::setVariableName" \
        (int idx, char *varname)

    char * pb_get_variable_name "BooleEnv::getVariableName" \
        (int idx)

    PBSet pb_random_set "random_set" (PBMonom variables, int length)
    void pb_set_random_seed "set_random_seed" (unsigned int seed)

    void PBPolyVector_set(PBPolyVector v, int i, PBPoly p)

    cdef cppclass PBConstant "BooleConstant":
        PBConstant()
        PBConstant(int val)
        bint isZero()
        bint isOne()
        bint isConstant()
        bint hasConstantPart()
        int deg()

    cdef cppclass PBVariableFactory "VariableFactory":
        PBVariableFactory(const PBRing&)

    cdef cppclass PBVarFactory "DefaultRinged<VariableFactory> ":
        void operator=(const PBVariableFactory&)
        PBVar operator()(long index)

    cdef cppclass PBMonomialFactory "MonomialFactory":
        PBMonomialFactory(const PBRing&)

    cdef cppclass PBMonomFactory "DefaultRinged<MonomialFactory> ":
        void operator=(const PBMonomialFactory&)
        PBMonom operator()()

    cdef cppclass PBPolynomialFactory "PolynomialFactory":
        PBPolynomialFactory(const PBRing &)

    cdef cppclass PBPolyFactory "DefaultRinged<PolynomialFactory> ":
        void operator=(const PBPolynomialFactory&)
        PBPoly operator()(long)
        PBPoly operator()(PBPoly)
        PBPoly operator()(PBMonom)

    cdef cppclass PBRefCounter:
        PBRefCounter()
        bint released()
