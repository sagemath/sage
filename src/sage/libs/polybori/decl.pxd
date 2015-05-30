

cdef extern from "pb_wrap.h":
    ctypedef struct std_string "std::string":
        char *(* c_str)()

    std_string new_stdstring "std::string"(char *str)

    cdef enum ordercodes "COrderEnums::ordercodes":
        pblp            "CTypes::lp"
        pbdlex          "CTypes::dlex"
        pbdp_asc        "CTypes::dp_asc"
        pbblock_dlex    "CTypes::block_dlex"
        pbblock_dp_asc  "CTypes::block_dp_asc"
        pbdp            "17"
        pbblock_dp      "19"

    cdef enum comparecodes "CCompareEnums::comparecodes":
        less_than               "CTypes::less_than"
        equality                "CTypes::equality"
        greater_than            "CTypes::greater_than"
        less_or_equal_max       "CTypes::less_or_equal_max"
        greater_or_equal_min    "CTypes::greater_or_equal_min"

    ctypedef struct PBNavigator "struct CCuddNavigator":
        PBNavigator (* thenBranch) ()
        PBNavigator (* elseBranch) ()
        int (* value "operator*")()
        bint (* isConstant)()
        bint (* isTerminated)()
        bint (* is_equal "operator==")(PBNavigator right)
        size_t (* hash)()

    # non-allocating versions
    PBNavigator* PBNavigator_construct \
            "Construct_p<CCuddNavigator, CCuddNavigator>" \
            (void* mem, PBNavigator N)
    void PBNavigator_destruct "Destruct<CCuddNavigator>"(PBNavigator *mem)

    ctypedef struct PBRing "BoolePolyRing"

    # non-allocating versions

    ctypedef struct PBBlockIter "COrderingBase::block_iterator":
        int (* value "operator*")()
        PBBlockIter (* next "operator++")()
        size_t (* hash)()

    bint PBBlockIter_equals "operator=="(PBBlockIter lhs, PBBlockIter rhs)

    ctypedef struct PBOrdering "COrderingBase":
        int (*getOrderCode)()
        int (*getBaseOrderCode)()
        void (*appendBlock)(int nextBlockStartIndex)
        PBBlockIter (*blockBegin)()
        PBBlockIter (*blockEnd)()
        int (*lastBlockStart)()
        bint isDegreeOrder()

    ctypedef struct PBSet "DefaultRinged<BooleSet>"
    ctypedef struct PBPoly "DefaultRinged<BoolePolynomial>"

    ctypedef struct PBVar "DefaultRinged<BooleVariable>":
        int (* index "index")()
        bint (* is_equal "operator==")(PBVar right)

    ctypedef struct PBRing "BoolePolyRing":
        size_t (*id) ()
        size_t (* nVariables )()
        PBVar (* variable )(int n)
        PBOrdering (*ordering )()
        void (*appendBlock "ordering().appendBlock")(int nextBlockStartIndex)
        void (*changeOrdering  )(int)
        size_t (*hash )()
        void (*setVariableName ) (int idx, char *varname)
        char* (*getVariableName )(int idx)
        PBRing (*clone )()
        PBPoly (*coerce )(PBPoly rhs)
        PBPoly (*zero )()
        PBPoly (*one )()

    cdef cppclass PBVarBlock "VariableBlock":
        PBVarBlock(int size, int start_index, int offset, bint reverse,\
                       PBRing r)
        PBVar operator()(int index)

    ctypedef struct PBMonomIter "BooleMonomial::const_iterator":
        int (* value "operator*")()
        int (* next "operator++")()
        bint (* equal "equal")(PBMonomIter rhs)
        size_t (* hash)()

    void PBMonomIter_destruct "Destruct<BooleMonomial::const_iterator>" \
            (PBMonomIter *mem)


    ctypedef struct PBMonomVarIter "BooleMonomial::variable_iterator":
        PBVar (* value "operator*")()
        int (* next "operator++")()
        bint (* equal "equal")(PBMonomVarIter rhs)

    void PBMonomVarIter_destruct "Destruct<BooleMonomial::variable_iterator>" \
            (PBMonomVarIter *mem)

    ctypedef struct PBMonom "DefaultRinged<BooleMonomial>":
        bint (* reducibleBy)(PBMonom rhs)
        int (* deg)()
        size_t (* hash)()
        size_t (* stableHash)()
        int (* firstIndex)()
        comparecodes (* compare)(PBMonom rhs)
        PBSet (* set)()
        PBSet (* divisors)()
        PBSet (* multiples)(PBMonom rhs)
        PBMonomIter (* begin)()
        PBMonomIter (* end)()
        PBMonomVarIter (*variableBegin)()
        PBMonomVarIter (*variableEnd)()
        void (* imul  "operator*=")(PBMonom right)
        void (* idiv  "operator/=")(PBMonom right)
        PBNavigator (* navigation)()
        PBMonom (* GCD)(PBMonom rhs)
        PBRing (* ring)()

    object PBMonom_to_str "_to_PyString<BooleMonomial>"(PBMonom *p)

    # Wrapping constructors
    PBMonom PBMonom_Constructor "BooleMonomial" (PBRing r)
    PBMonom PBMonom_Constructor_var "BooleMonomial" (PBVar m)

    ctypedef struct PBSetIter "BooleSet::const_iterator":
        PBMonom (* value "operator*")()
        int (* next "operator++")()
        bint (* equal)(PBSetIter rhs)

    PBSetIter* PBSetIter_construct_begin \
            "construct_bset_begin" (PBSetIter *mem, PBSet)

    PBSetIter* PBSetIter_construct_end \
            "construct_bset_end" (PBSetIter *mem, PBSet)

    void PBSetIter_destruct "Destruct<BooleSet::const_iterator>"(PBSetIter *mem)

    ctypedef struct PBSet "DefaultRinged<BooleSet>":
        bint (* owns)(PBMonom val)
        int (* nNodes)()
        int (* nSupport)()
        int (* size)()
        size_t (* hash)()
        size_t (* stableHash)()
        PBNavigator (* navigation)()
        PBSet (* cartesianProduct)(PBSet rhs)
        PBSet (* diff)(PBSet rhs)
        PBSet (* divide)(PBMonom rhs)
        PBSet (* change)(int idx)
        PBSet (* subset0)(int i)
        PBSet (* subset1)(int i)
        PBSet (* unite)(PBSet rhs)
        PBSet (* intersect)(PBSet rhs)
        PBMonom (* usedVariables)()
        PBSet (* divisorsOf)(PBMonom rhs)
        PBSet (* multiplesOf)(PBMonom rhs)
        double (* sizeDouble)()
        PBSetIter (* begin)()
        PBSetIter (* end)()
        bint (* isZero)()
        bint (* isOne)()

    PBSet pb_include_divisors "include_divisors" (PBSet p)
    PBSet pb_minimal_elements "minimal_elements" (PBSet p)

    object PBSet_to_str "_to_PyString<BooleSet>"(PBSet *p)

    # non-allocating versions
    PBSet PBSet_Constructor_ring "BooleSet"(PBRing r)
    PBSet PBSet_Constructor_poly "BooleSet"(PBPoly p)
    PBSet PBSet_Constructor_nav "BooleSet" (PBNavigator d, PBRing r)
    PBSet PBSet_Constructor_indsetset "BooleSet" \
            (int ind, PBNavigator a, PBNavigator b, PBRing r)


    ctypedef struct PBPolyIter "BoolePolynomial::ordered_iterator":
        PBMonom (* value "operator*")()
        int (* next "operator++")()
        bint (* equal)(PBPolyIter rhs)

    void PBPolyIter_destruct "Destruct<BoolePolynomial::ordered_iterator>"(PBPolyIter *mem)

    ctypedef struct PBPoly "DefaultRinged<BoolePolynomial>":
        int (* deg)()
        int (* leadDeg)()
        int (* lexLeadDeg)()
        int (* totalDeg)()
        int (* length)()
        int (* eliminationLength)()
        int (* nUsedVariables)()
        int (* nNodes)()
        size_t (* hash)()
        size_t (* stableHash)()
        bint (* isZero)()
        bint (* isOne)()
        bint (* isConstant)()
        bint (* isSingleton)()
        bint (* isPair)()
        bint (* isSingletonOrPair)()
        bint (* hasConstantPart)()
        bint (* firstReducibleBy)(PBPoly rhs)
        PBMonom (* lead)()
        PBMonom (* lexLead)()
        PBMonom (* firstTerm)()
        PBMonom (* usedVariables)()
        PBPoly (* gradedPart)(int d)
        PBSet (* diagram)()
        PBSet (* set)()
        PBSet (* leadDivisors)()
        PBNavigator (* navigation)()
        PBPolyIter (* orderedBegin)()
        PBPolyIter (* orderedEnd)()
        void (* iadd  "operator+=")(PBPoly right)
        void (* iadd_PBMonom  "operator+=")(PBMonom right)
        void (* imul  "operator*=")(PBPoly right)
        void (* imul_monom  "operator*=")(PBMonom right)
        void (* idiv  "operator/=")(PBPoly right)
        void (* idiv_monom  "operator/=")(PBMonom right)
        bint (* is_equal  "operator==")(PBPoly right)

    PBSet pb_zeros "zeros" (PBPoly p, PBSet s)
    PBPoly pb_spoly "spoly" (PBPoly p, PBPoly r)

    PBPoly pb_map_every_x_to_x_plus_one "map_every_x_to_x_plus_one" (PBPoly)

    # construction
    PBVar PBVar_Constructor "BooleVariable" (int idx, PBRing r)
    PBRing PBRing_Constructor "BoolePolyRing" (int nvars, int order)

    # construction
    PBPoly PBPoly_Constructor_ring "BoolePolynomial"(PBRing r)
    PBPoly PBPoly_Constructor_set "BoolePolynomial" (PBSet d)
    PBPoly PBPoly_Constructor_monom "BoolePolynomial" (PBMonom d)
    PBPoly PBPoly_Constructor_var "BoolePolynomial" (PBVar d)
    PBPoly PBPoly_Constructor_int_ring "BoolePolynomial" (int d, PBRing r)

    object PBPoly_to_str "_to_PyString<DefaultRinged<BoolePolynomial> >"(PBPoly *p)

    ctypedef struct PBPolyVectorIter \
            "std::vector<BoolePolynomial>::iterator ":
        PBPoly (* value "operator*")()
        int (* next "operator++")()

    bint PBPolyVectorIter_equal "operator=="(PBPolyVectorIter lhs, PBPolyVectorIter rhs)

    void PBPolyVectorIter_destruct "Destruct<std::vector<BoolePolynomial>::iterator>"(PBPolyVectorIter *mem)

    ctypedef struct PBPolyVector "std::vector<BoolePolynomial>":
        int (* size)()
        PBPoly (* get "operator[]")(int)
        PBPolyVectorIter (* begin)()
        PBPolyVectorIter (* end)()
        void (* push_back)(PBPoly val)

    PBPolyVector* PBPolyVector_construct \
            "Construct< std::vector<BoolePolynomial> >"(void *mem)
    void PBPolyVector_destruct "Destruct< std::vector<BoolePolynomial> >"\
            (PBPolyVector *mem)

    ctypedef struct PBPolyEntry "PolyEntry":
        PBPoly p

    ctypedef struct PBRedStrategy "ReductionStrategy":
        PBSet leadingTerms
        PBSet minimalLeadingTerms
        PBSet leadingTerms11
        PBSet leadingTerms00
        PBSet llReductor
        PBSet monomials
        PBSet monomials_plus_one
        bint optBrutalReductions
        bint optLL
        PBPoly (* nf)(PBPoly p)
        bint optRedTailDegGrowth
        bint optRedTail
        void (* setupSetsForLastElement)()
        bint (* canRewrite)(PBPoly p)
        void (* addGenerator)(PBPoly p)
        int select1(PBMonom p)

        int (* select_short)(PBPoly p)
        PBPoly (* headNormalForm)(PBPoly p)
        PBPoly (* reducedNormalForm)(PBPoly p)

        int (* size)()
        PBPolyEntry (* get "operator[]")(int)

    ctypedef struct PBFglmStrategy "WrappedPtr<FGLMStrategy>":
        PBPolyVector (* main "operator->()->main")()

    PBFglmStrategy PBFglmStrategy_Constructor "WrappedPtr<FGLMStrategy>" \
        (PBRing from_ring, PBRing to_ring, PBPolyVector vec)

    cdef cppclass PBGBStrategy "GroebnerStrategy":
        PBGBStrategy(PBRing)
        bint reduceByTailReduced
        PBRedStrategy generators

        bint enabledLog "enabledLog"
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
        bint (* containsOne)()
        bint (* variableHasValue)(int idx)
        int (* addGenerator)(PBPoly, bint is_impl)
        void (* addGeneratorDelayed)(PBPoly)
        void (* addAsYouWish)(PBPoly)
        void (* symmGB_F2)()
        void (* llReduceAll)()
        void (* cleanTopByChainCriterion  "pairs.cleanTopByChainCriterion")()
        int (* nGenerators  "generators.size")()
        int (* npairs  "pairs.queue.size")()
        int (* suggestPluginVariable)()
        PBPoly (* nextSpoly)()
        PBPoly (* nf)(PBPoly p)
        PBPolyVector (* allGenerators)()
        PBPolyVector (* minimalize)()
        PBPolyVector (* minimalizeAndTailReduce)()
        PBPolyVector (* faugereStepDense)(PBPolyVector v)
        bint (* generators_leadingTerms_owns  "generators.leadingTerms.owns")(PBMonom term)

    PBGBStrategy PBGBStrategy_Constructor "WrappedPtr<GroebnerStrategy>" \
        (PBRing r)


    int (* pairs_top_sugar)(PBGBStrategy strat)
    PBPolyVector (* someNextDegreeSpolys)(PBGBStrategy strat, int n)
    PBPolyVector (* nextDegreeSpolys)(PBGBStrategy strat)
    PBPolyVector (* small_next_degree_spolys)(PBGBStrategy strat, double f, int n)
    void (* implications)(PBGBStrategy strat, int i)

    PBPoly (* cheap_reductions)(PBRedStrategy strat, PBPoly p)

    PBPoly GB_get_ith_gen "get_ith_gen" (PBGBStrategy strat, int i)

    # non-allocating versions
    PBGBStrategy* PBGBStrategy_construct_ring \
         "Construct_p<GroebnerStrategy, BoolePolyRing>"(void *mem, PBRing r)

    PBGBStrategy* PBGBStrategy_construct_gbstrategy \
            "Construct_p<GroebnerStrategy, GroebnerStrategy>" \
            (void *mem, PBGBStrategy strat)
    void PBGBStrategy_destruct "Destruct<GroebnerStrategy>"(PBGBStrategy *mem)

    PBFglmStrategy* PBFglmStrategy_construct "Construct_ppp<FGLMStrategy, BoolePolyRing, BoolePolyRing, PolynomialVector>"(void *mem, PBRing from_ring, PBRing to_ring, PBPolyVector vec)
    void PBFglmStrategy_destruct "Destruct<FGLMStrategy>"(PBFglmStrategy *mem)


    # allocating versions

    PBRedStrategy* PBRedStrategy_new "New_p<ReductionStrategy, BoolePolyRing>"(PBRing)
    void PBRedStrategy_delete "Delete<ReductionStrategy>"(PBRedStrategy *mem)

    # Additional routines

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

    #int pb_select1 "select1" (GBStrategy s, PBMonom m)

    PBSet pb_recursively_insert "recursively_insert"(PBNavigator p,
                                                int idx, PBSet m)
    PBPoly pb_ll_red_nf_noredsb "ll_red_nf_noredsb"(PBPoly p,
                                                PBSet reductors)
    PBPoly pb_ll_red_nf_noredsb_single_recursive_call "ll_red_nf_noredsb_single_recursive_call"(PBPoly p,
                                                                                                PBSet reductors)


    PBPoly pb_ll_red_nf "ll_red_nf"(PBPoly p, PBSet reductors)

    PBSet pb_mod_mon_set "mod_mon_set"(PBSet as, PBSet vs)

    PBPolyVector pb_parallel_reduce "parallel_reduce" \
        (PBPolyVector inp, PBGBStrategy strat, int average_steps, double delay_f)

    PBPolyVector pb_gauss_on_polys "gauss_on_polys" (PBPolyVector inp)

    PBPolyVector pb_easy_linear_factors "easy_linear_factors"(PBPoly p)

    PBPoly pb_substitute_variables "substitute_variables<BoolePolyRing, std::vector<BoolePolynomial>, BoolePolynomial>" (PBRing ring, PBPolyVector vec, PBPoly poly) except +

    void pb_set_variable_name "BooleEnv::setVariableName" \
        (int idx, char *varname)

    char * pb_get_variable_name "BooleEnv::getVariableName" \
        (int idx)

    PBSet pb_random_set "random_set" (PBMonom variables, int length)
    void pb_set_random_seed "set_random_seed" (unsigned int seed)

    #M4RI initialization
    void m4ri_build_all_codes()
    void m4ri_destroy_all_codes()

    void PBPolyVector_set(PBPolyVector v, int i, PBPoly p)


    ctypedef struct PBConstant "struct BooleConstant":
        bint isZero()
        bint isOne()
        bint isConstant()
        bint hasConstantPart()
        int deg()

    PBConstant* PBConstant_construct "Construct_p<BooleConstant, int>" \
        (void* mem, int val)

    ctypedef struct PBVarFactory "DefaultRinged<VariableFactory>":
        PBVar (* call "operator()")(int index)

    ctypedef struct PBMonomFactory "DefaultRinged<MonomialFactory>":
       PBMonom (* call "operator()")()

    ctypedef struct PBPolyFactory "DefaultRinged<PolynomialFactory>":
        PBPoly (* call_int "operator()")(int)
        PBPoly (* call_poly "operator()")(PBPoly)
        PBPoly (* call_monom "operator()")(PBMonom)

    PBVarFactory PBVarFactory_Constructor "VariableFactory" (PBRing r)
    PBMonomFactory PBMonomFactory_Constructor "MonomialFactory" (PBRing r)
    PBPolyFactory PBPolyFactory_Constructor "PolynomialFactory" (PBRing r)

    long* init_counter "new long" (int)
    void kill_counter "Destruct<long>" (long*)

    ctypedef struct PBRefCounter:
        bint (* released) ()

    PBRefCounter PBRefCounter_Constructor "PBRefCounter" ()
