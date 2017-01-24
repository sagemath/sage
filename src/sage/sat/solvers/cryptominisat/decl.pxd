from libc.stdint cimport uint32_t, uint64_t
from libcpp.vector cimport vector

cdef extern from "Solver.h" namespace "CMSat":
    cdef cppclass lbool:
        bint getBool()
        bint operator==(lbool b)

    lbool l_True
    lbool l_False
    lbool l_Undef

    ctypedef uint32_t Var

    cdef cppclass Lit:
        Lit()
        Lit(Var var, bint sign)
        Var var()
        bint sign()

    cdef cppclass vec[T]:
        vec()
        vec(Py_ssize_t size)
        uint32_t size()
        void     pop()
        void push(T&)
        T& operator[](int)
        T& at(int)

    cdef cppclass Clause:
        uint32_t size()
        void pop()
        bint isXor()
        bint learnt()
        Lit& operator[](uint32_t i)

    cdef cppclass SolverConf:
        SolverConf()
        SolverConf(SolverConf &)
        float    random_var_freq    #The frequency with which the decision heuristic tries to choose
                                    #a random variable.  (default 0.02) NOTE: This is really
                                    #strange. If the number of variables set is large, then the
                                    #random chance is in fact _far_ lower than this value. This is
                                    #because the algorithm tries to set one variable randomly, but
                                    #if that variable is already set, then it _silently_ fails, and
                                    #moves on (doing non-random flip)!
        float    clause_decay       #Inverse of the clause activity decay factor. Only applies if
                                    #using MiniSat-style clause activities (default: 1 / 0.999)
        int      restart_first      #The initial restart limit.  (default 100)
        float    restart_inc        #The factor with which the restart limit is multiplied in each
                                    #restart.  (default 1.5)
        float    learntsize_factor  #The initial limit for learnt clauses is a factor of the
                                    #original clauses.  (default 1 / 3)
        float    learntsize_inc     #The limit for learnt clauses is multiplied with this factor
                                    #each restart.  (default 1.1)
        bint     expensive_ccmin    #Should clause minimisation qby Sorensson&Biere be used?
                                    #(default TRUE)
        int      polarity_mode      #Controls which polarity the decision heuristic chooses. Auto
                                    #means Jeroslow-Wang (default: polarity_auto)
        int      verbosity          #Verbosity level. 0=silent, 1=some progress report, 2=lots of
                                    #report, 3 = all report (default 2)
        Var      restrictPickBranch #Pick variables to branch on preferentally from the highest [0,
                                    #restrictedPickBranch]. If set to 0, preferentiality is turned
                                    #off (i.e. picked randomly between [0, all])
        uint32_t simpBurstSConf
        float simpStartMult
        float simpStartMMult
        bint doPerformPreSimp
        float failedLitMultiplier

        # Optimisations to do
        bint      doFindXors         #Automatically find non-binary xor clauses and convert them to
                                     #xor clauses
        bint      doFindEqLits       #Automatically find binary xor clauses (i.e. variable equi- and
                                     #antivalences)
        bint      doRegFindEqLits    #Regularly find binary xor clauses (i.e. variable equi- and
                                     #antivalences)
        bint      doReplace          #Should var-replacing be performed? If set to FALSE, equi- and
                                     #antivalent variables will not be replaced with one
                                     #another. NOTE: This precludes using a lot of the algorithms!
        bint      doConglXors        #Do variable elimination at the XOR-level (xor-ing 2 xor
                                     #clauses thereby removing a variable)
        bint      doHeuleProcess     #Perform local subsitutuion as per Heule's theis
        bint      doSchedSimp        #Should simplifyProblem() be scheduled regularly? (if set to
                                     #FALSE, a lot of optimisations are disabled)
        bint      doSatELite         #Should try to subsume & self-subsuming resolve &
                                     #variable-eliminate & block-clause eliminate?
        bint      doXorSubsumption   #Should try to subsume & local-subsitute xor clauses
        bint      doHyperBinRes      #Should try carry out hyper-binary resolution
        bint      doBlockedClause    #Should try to remove blocked clauses
        bint      doVarElim          #Perform variable elimination
        bint      doSubsume1         #Perform self-subsuming resolution
        bint      doClausVivif       #Perform asymmetric branching at the beginning of the solving
        bint      doSortWatched      #Sort watchlists according to size&type: binary, tertiary,
                                     #normal (>3-long), xor clauses
        bint      doMinimLearntMore  #Perform learnt-clause minimisation using watchists' binary and
                                     #tertiary clauses?  ("strong minimization" in PrecoSat)
        bint      doMinimLMoreRecur  #Always perform recursive/transitive on-the-fly self
                                     #self-subsuming resolution --> an enhancement of "strong
                                     #minimization" of PrecoSat
        bint      doFailedLit        #Carry out Failed literal probing + doubly propagated literal
                                     #detection + 2-long xor clause detection during failed literal
                                     #probing + hyper-binary resolution
        bint      doRemUselessBins   #Should try to remove useless binary clauses at the beginning
                                     #of solving?
        bint      doSubsWBins
        bint      doSubsWNonExistBins #Try to do subsumption and self-subsuming resolution with
                                      #non-existent binary clauses (i.e. binary clauses that don't
                                      #exist but COULD exists)
        bint      doRemUselessLBins   #Try to remove useless learnt binary clauses
        bint      doPrintAvgBranch
        bint      doCacheOTFSSR
        bint      doCacheOTFSSRSet
        bint      doExtendedSCC
        bint      doCalcReach #Calculate reachability, and influence variable decisions with that
        bint      doBXor
        uint64_t  maxConfl

        #  interrupting & dumping
        uint32_t  maxRestarts
        bint      needToDumpLearnts  #If set to TRUE, learnt clauses will be dumped to the file
                                     #speified by "learntsFilename"
        bint      needToDumpOrig     #If set to TRUE, a simplified version of the original
                                     #clause-set will be dumped to the file speified by
                                     #"origFilename". The solution to this file should perfectly
                                     #satisfy the problem
        uint32_t  maxDumpLearntsSize #When dumping the learnt clauses, this is the maximum clause
                                     #size that should be dumped
        bint      libraryUsage       #Set to true if not used as a library. In fact, this is TRUE by
                                     #default, and Main.cpp sets it to "FALSE". Disables some
                                     #simplifications at the beginning of solving (mostly
                                     #performStepsBeforeSolve() )
        bint      greedyUnbound      #If set, then variables will be greedily unbounded (set to
                                     #l_Undef). This is EXPERIMENTAL

        uint32_t origSeed

    cdef cppclass GaussConf:
        pass


    cdef struct RetClause:
        bint is_xor
        bint right_hand_side
        vector[Lit] lits

    cdef cppclass Solver:
        Solver()
        Solver(SolverConf& conf, GaussConf& _gaussconfig)

        # Problem specification:
        Var     newVar() except +
        Var     newVar(bint dvar) except +
        bint    addClause(vec[Lit] ps)
        bint    addLearntClause(vec[Lit] ps)
        bint    addXorClause(vec[Lit] ps, bint xorEqualFalse) except +

        # Solving
        lbool    solve       (vec[Lit]& assumps) # Search for a model that respects a given set of
                                                 # assumptions.
        lbool    solve       ()                  # Search without assumptions.
        void     handleSATSolution()             # Extends model, if needed, and fills "model"
        void     handleUNSATSolution()           # If conflict really was zero-length, sets OK to
                                                 # false
        bint     okay         ()                 # FALSE means solver is in a conflicting state

        # Variable mode:

        void    setDecisionVar (Var v, bint b) # Declare if a variable should be eligible for
                                               # selection in the decision heuristic.
        void    addBranchingVariable (Var v)

        # Read state:
        lbool   value      (Var x) #The current value of a variable.
        lbool   value      (Lit p) #The current value of a literal.
        lbool   modelValue (Lit p) #The value of a literal in the last model. The last call to solve
                                   #must have been satisfiable.
        uint32_t     nAssigns   () #The current number of assigned literals.
        uint32_t     nClauses   () #The current number of original clauses.
        uint32_t     nLiterals  () #The current number of total literals.
        uint32_t     nLearnts   () #The current number of learnt clauses.
        uint32_t     nVars      () #The current number of variables.

       # Extra results: (read-only member variable)
        vec[lbool] model            # If problem is satisfiable, this vector contains the model (if
                                    # any).
        vec[Lit]   conflict         # If problem is unsatisfiable (possibly under assumptions), this
                                    # vector represent the final conflict clause expressed in the
                                    # assumptions.

        # Logging
        void needStats()            # Prepares the solver to output statistics
        void needProofGraph()       # Prepares the solver to output proof graphs during solving

        vec[Lit] get_unitary_learnts() # return the set of unitary learnt clauses

        vector[RetClause] dumpOrigClauses() # return original simplified clauses

