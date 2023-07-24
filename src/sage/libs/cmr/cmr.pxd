# -*- python -*-
# distutils: libraries = cmr

# (progn (replace-regexp "/[*]\\(.\\|\n\\)*?[*]/" "" nil (point) (point-max)) (replace-regexp "[;{}]" ""  nil (point) (point-max)) (replace-regexp "CMR_EXPORT *" ""  nil (point) (point-max)) (replace-regexp "bool" "bint" nil (point) (point-max)))

cdef extern from "stdbool.h":

    ctypedef int bool

cdef extern from "cmr/env.h":
    const int CMR_OKAY
    const int CMR_ERROR_INPUT
    const int CMR_ERROR_OUTPUT
    const int CMR_ERROR_MEMORY
    const int CMR_ERROR_INVALID
    const int CMR_ERROR_OVERFLOW
    const int CMR_ERROR_TIMEOUT

    ctypedef int CMR_ERROR

    ctypedef struct CMR

    int CMRcreateEnvironment(CMR** pcmr)
    int CMRfreeEnvironment(CMR** pcmr)

    char* CMRgetErrorMessage(CMR* cmr)
    void CMRclearErrorMessage(CMR* cmr)

cdef extern from "cmr/matrix.h":

    ctypedef struct CMR_SUBMAT:
        size_t numRows
        size_t* rows
        size_t numColumns
        size_t* columns

    CMR_ERROR CMRsubmatCreate(CMR* cmr, size_t numRows, size_t numColumns, CMR_SUBMAT** psubmatrix)
    CMR_ERROR CMRsubmatCreate1x1(CMR* cmr, size_t row, size_t column, CMR_SUBMAT** psubmatrix)
    CMR_ERROR CMRsubmatFree(CMR* cmr, CMR_SUBMAT** psubmatrix)

    ctypedef struct CMR_CHRMAT:
        size_t numRows
        size_t numColumns
        size_t numNonzeros
        size_t* rowSlice
        size_t* entryColumns
        char* entryValues

    CMR_ERROR CMRchrmatCreate(CMR* cmr, CMR_CHRMAT** presult, int numRows, int numColumns, int numNonzeros)
    CMR_ERROR CMRchrmatSortNonzeros(CMR* cmr, CMR_CHRMAT* matrix)
    # CMR_ERROR CMRchrmatPrintDense(CMR* cmr, CMR_CHRMAT* matrix, FILE* stream, char zeroChar, bint header)
    CMR_ERROR CMRchrmatFindEntry(CMR_CHRMAT* matrix, size_t row, size_t column, size_t* pentry)
    CMR_ERROR CMRchrmatZoomSubmat(CMR* cmr, CMR_CHRMAT* matrix, CMR_SUBMAT* submatrix, CMR_CHRMAT** presult)
    CMR_ERROR CMRchrmatFree(CMR* cmr, CMR_CHRMAT** pmatrix)

cdef extern from "cmr/k_modular.h":

    CMR_ERROR CMRtestUnimodularity(CMR* cmr, CMR_CHRMAT* matrix, int* pisUnimodular)
    CMR_ERROR CMRtestStrongUnimodularity(CMR* cmr, CMR_CHRMAT* matrix, bool* pisStronglyUnimodular)
    CMR_ERROR CMRtestKmodularity(CMR* cmr, CMR_CHRMAT* matrix, bool* pisKmodular, size_t* pk)
    CMR_ERROR CMRtestStrongKmodularity(CMR* cmr, CMR_CHRMAT* matrix, bool* pisStronglyKmodular, size_t* pk)

cdef extern from "cmr/camion.h":

    ctypedef struct CMR_CAMION_STATISTICS:
        size_t totalCount
        double totalTime

    CMR_ERROR CMRstatsCamionInit(CMR_CAMION_STATISTICS* stats)
    # CMR_ERROR CMRstatsCamionPrint(FILE* stream, CMR_CAMION_STATISTICS* stats, const char* prefix)
    CMR_ERROR CMRtestCamionSigned(CMR* cmr, CMR_CHRMAT* matrix, bool* pisCamionSigned, CMR_SUBMAT** psubmatrix, CMR_CAMION_STATISTICS* stats, double timeLimit)
    CMR_ERROR CMRcomputeCamionSigned(CMR* cmr, CMR_CHRMAT* matrix, bool* pwasCamionSigned, CMR_SUBMAT** psubmatrix, CMR_CAMION_STATISTICS* stats, double timeLimit)

cdef extern from "cmr/matroid.h":

    ctypedef struct CMR_MINOR:
        size_t numPivots
        size_t* pivotRows
        size_t* pivotColumns
        CMR_SUBMAT* remainingSubmatrix

    CMR_ERROR CMRminorCreate(CMR* cmr, CMR_MINOR** pminor, size_t numPivots, CMR_SUBMAT* submatrix)
    CMR_ERROR CMRminorFree(CMR* cmr, CMR_MINOR** pminor)

cdef extern from "cmr/element.h":

    ctypedef int CMR_ELEMENT

    const char* CMRelementString(CMR_ELEMENT element, char* buffer)
    bint CMRelementIsValid(CMR_ELEMENT element)
    CMR_ELEMENT CMRrowToElement(size_t row)
    CMR_ELEMENT CMRcolumnToElement(size_t column)
    bint CMRelementIsRow(CMR_ELEMENT element)
    size_t CMRelementToRowIndex(CMR_ELEMENT element)
    bint CMRelementIsColumn(CMR_ELEMENT element)
    size_t CMRelementToColumnIndex(CMR_ELEMENT element)
    CMR_ELEMENT CMRelementTranspose(CMR_ELEMENT element)

cdef extern from "cmr/separation.h":

    ctypedef struct CMR_SEPA:
        unsigned char* rowsToPart
        unsigned char* columnsToPart
        size_t numRows[2]
        size_t numColumns[2]
        size_t* rows[2]
        size_t* columns[2]
        size_t extraRows[2][2]
        size_t extraColumns[2][2]
        unsigned char* indicatorMemory
        size_t* elementMemory

    CMR_ERROR CMRsepaCreate(CMR* cmr, size_t numRows, size_t numColumns, CMR_SEPA** psepa)
    CMR_ERROR CMRsepaInitialize(CMR* cmr, CMR_SEPA* separation, size_t firstExtraRow0, size_t firstExtraColumn1, size_t firstExtraRow1, size_t firstExtraColumn0, size_t secondExtraRow0, size_t secondExtraColumn1, size_t secondExtraRow1, size_t secondExtraColumn0)
    CMR_ERROR CMRsepaInitializeMatrix(CMR* cmr, CMR_SEPA* separation, CMR_CHRMAT* matrix, unsigned char totalRank)
    CMR_ERROR CMRsepaFree(CMR* cmr, CMR_SEPA** psepa)
    unsigned char CMRsepaRankBottomLeft(CMR_SEPA* sepa)
    unsigned char CMRsepaRankTopRight(CMR_SEPA* sepa)
    unsigned char CMRsepaRank(CMR_SEPA* sepa)
    CMR_ERROR CMRsepaCheckTernary(CMR* cmr, CMR_SEPA* sepa, CMR_CHRMAT* matrix, CMR_SUBMAT* submatrix, bool* pisTernary, CMR_SUBMAT** psubmatrix)
    CMR_ERROR CMRoneSum(CMR* cmr, CMR_CHRMAT* first, CMR_CHRMAT* second, CMR_CHRMAT** presult)
    CMR_ERROR CMRtwoSum(CMR* cmr, CMR_CHRMAT* first, CMR_CHRMAT* second, CMR_ELEMENT firstMarker, CMR_ELEMENT secondMarker, CMR_CHRMAT** presult)
    CMR_ERROR CMRthreeSum(CMR* cmr, CMR_CHRMAT* first, CMR_CHRMAT* second, CMR_ELEMENT firstMarker1, CMR_ELEMENT secondMarker1, CMR_ELEMENT firstMarker2, CMR_ELEMENT secondMarker2, CMR_CHRMAT** presult)

cdef extern from "cmr/graph.h":

    ctypedef int CMR_GRAPH_NODE
    ctypedef int CMR_GRAPH_EDGE
    ctypedef int CMR_GRAPH_ITER

    ctypedef struct CMR_GRAPH_NODE_DATA:
        int prev
        int next
        int firstOut

    ctypedef struct CMR_GRAPH_ARC_DATA:
        int target
        int prev
        int next

    ctypedef struct CMR_GRAPH:
        size_t numNodes
        size_t memNodes
        CMR_GRAPH_NODE_DATA* nodes
        int firstNode
        int freeNode
        size_t numEdges
        size_t memEdges
        CMR_GRAPH_ARC_DATA* arcs
        int freeEdge

    size_t CMRgraphMemNodes(CMR_GRAPH* graph)
    size_t CMRgraphNumNodes(CMR_GRAPH* graph)
    size_t CMRgraphMemEdges(CMR_GRAPH* graph)
    size_t CMRgraphNumEdges(CMR_GRAPH* graph)
    CMR_GRAPH_NODE CMRgraphEdgeU(CMR_GRAPH* graph, CMR_GRAPH_EDGE e)
    CMR_GRAPH_NODE CMRgraphEdgeV(CMR_GRAPH* graph, CMR_GRAPH_EDGE e)
    CMR_ERROR CMRgraphCreateEmpty(CMR* cmr, CMR_GRAPH** pgraph, int memNodes, int memEdges)
    CMR_ERROR CMRgraphFree(CMR* cmr, CMR_GRAPH** pgraph)
    CMR_ERROR CMRgraphClear(CMR* cmr, CMR_GRAPH* graph)
    CMR_ERROR CMRgraphAddNode(CMR* cmr, CMR_GRAPH* graph, CMR_GRAPH_NODE* pnode)
    CMR_ERROR CMRgraphAddEdge(CMR* cmr, CMR_GRAPH* graph, CMR_GRAPH_NODE u, CMR_GRAPH_NODE v, CMR_GRAPH_EDGE* pedge)
    CMR_ERROR CMRgraphDeleteNode(CMR* cmr, CMR_GRAPH* graph, CMR_GRAPH_NODE v)
    CMR_ERROR CMRgraphDeleteEdge(CMR* cmr, CMR_GRAPH* graph, CMR_GRAPH_EDGE e)
    CMR_GRAPH_NODE CMRgraphNodesFirst(CMR_GRAPH* graph)
    bint CMRgraphNodesValid(CMR_GRAPH* graph, CMR_GRAPH_NODE v)
    CMR_GRAPH_NODE CMRgraphNodesNext(CMR_GRAPH* graph, CMR_GRAPH_NODE v)
    CMR_GRAPH_ITER CMRgraphIncFirst(CMR_GRAPH* graph, CMR_GRAPH_NODE v)
    bint CMRgraphIncValid(CMR_GRAPH* graph, CMR_GRAPH_ITER i)
    CMR_GRAPH_ITER CMRgraphIncNext(CMR_GRAPH* graph, CMR_GRAPH_ITER i)
    CMR_GRAPH_EDGE CMRgraphIncEdge(CMR_GRAPH* graph, CMR_GRAPH_ITER i)
    CMR_GRAPH_NODE CMRgraphIncSource(CMR_GRAPH* graph, CMR_GRAPH_ITER i)
    CMR_GRAPH_NODE CMRgraphIncTarget(CMR_GRAPH* graph, CMR_GRAPH_ITER i)
    CMR_GRAPH_ITER CMRgraphEdgesFirst(CMR_GRAPH* graph)
    CMR_GRAPH_ITER CMRgraphEdgesNext(CMR_GRAPH* graph, CMR_GRAPH_ITER i)
    bint CMRgraphEdgesValid(CMR_GRAPH* graph, CMR_GRAPH_ITER i)
    CMR_GRAPH_EDGE CMRgraphEdgesEdge(CMR_GRAPH* graph, CMR_GRAPH_ITER i)
    # CMR_ERROR CMRgraphPrint(FILE* stream, CMR_GRAPH* graph)
    CMR_ERROR CMRgraphMergeNodes(CMR* cmr, CMR_GRAPH* graph, CMR_GRAPH_NODE u, CMR_GRAPH_NODE v)
    # CMR_ERROR CMRgraphCreateFromEdgeList(CMR* cmr, CMR_GRAPH** pgraph, CMR_ELEMENT** pedgeElements, char*** pnodeLabels, FILE* stream)

cdef extern from "cmr/graphic.h":

    ctypedef struct CMR_GRAPHIC_STATISTICS:
        size_t totalCount
        double totalTime
        size_t checkCount
        double checkTime
        size_t applyCount
        double applyTime
        size_t transposeCount
        double transposeTime

    CMR_ERROR CMRstatsGraphicInit(CMR_GRAPHIC_STATISTICS* stats)
    # CMR_ERROR CMRstatsGraphicPrint(FILE* stream, CMR_GRAPHIC_STATISTICS* stats, const char* prefix)
    CMR_ERROR CMRcomputeGraphicMatrix(CMR* cmr, CMR_GRAPH* graph, CMR_CHRMAT** pmatrix, CMR_CHRMAT** ptranspose, int numForestEdges, CMR_GRAPH_EDGE* forestEdges, int numCoforestEdges, CMR_GRAPH_EDGE* coforestEdges, bool* pisCorrectForest)
    CMR_ERROR CMRtestGraphicMatrix(CMR* cmr, CMR_CHRMAT* matrix, bool* pisGraphic, CMR_GRAPH** pgraph, CMR_GRAPH_EDGE** pforestEdges, CMR_GRAPH_EDGE** pcoforestEdges, CMR_SUBMAT** psubmatrix, CMR_GRAPHIC_STATISTICS* stats, double timeLimit)
    CMR_ERROR CMRtestCographicMatrix(CMR* cmr, CMR_CHRMAT* matrix, bool* pisCographic, CMR_GRAPH** pgraph, CMR_GRAPH_EDGE** pforestEdges, CMR_GRAPH_EDGE** pcoforestEdges, CMR_SUBMAT** psubmatrix, CMR_GRAPHIC_STATISTICS* stats, double timeLimit)
    CMR_ERROR CMRtestGraphicColumnSubmatrixGreedy(CMR* cmr, CMR_CHRMAT* transpose, size_t* orderedColumns, CMR_SUBMAT** psubmatrix)

cdef extern from "cmr/series_parallel.h":

    ctypedef struct CMR_SP_STATISTICS:
        size_t totalCount
        double totalTime
        size_t reduceCount
        double reduceTime
        size_t wheelCount
        double wheelTime
        size_t nonbinaryCount
        double nonbinaryTime

    CMR_ERROR CMRstatsSeriesParallelInit(CMR_SP_STATISTICS* stats)

    # CMR_ERROR CMRstatsSeriesParallelPrint(FILE* stream, CMR_SP_STATISTICS* stats, const char* prefix)

    ctypedef struct CMR_SP_REDUCTION:
        CMR_ELEMENT element
        CMR_ELEMENT mate

    char* CMRspReductionString(CMR_SP_REDUCTION reduction, char* buffer)

    bint CMRspIsRow(CMR_SP_REDUCTION reduction)
    bint CMRspIsColumn(CMR_SP_REDUCTION reduction)
    bint CMRspIsZero(CMR_SP_REDUCTION reduction)
    bint CMRspIsUnit(CMR_SP_REDUCTION reduction)
    bint CMRspIsCopy(CMR_SP_REDUCTION reduction)
    bint CMRspIsValid(CMR_SP_REDUCTION reduction)

    CMR_ERROR CMRtestBinarySeriesParallel(CMR* cmr, CMR_CHRMAT* matrix, bool* pisSeriesParallel, CMR_SP_REDUCTION* reductions, size_t* pnumReductions, CMR_SUBMAT** preducedSubmatrix, CMR_SUBMAT** pviolatorSubmatrix, CMR_SP_STATISTICS* stats, double timeLimit)
    CMR_ERROR CMRtestTernarySeriesParallel(CMR* cmr, CMR_CHRMAT* matrix, bool* pisSeriesParallel, CMR_SP_REDUCTION* reductions, size_t* pnumReductions, CMR_SUBMAT** preducedSubmatrix, CMR_SUBMAT** pviolatorSubmatrix, CMR_SP_STATISTICS* stats, double timeLimit)
    CMR_ERROR CMRdecomposeBinarySeriesParallel(CMR* cmr, CMR_CHRMAT* matrix, bool* pisSeriesParallel, CMR_SP_REDUCTION* reductions, size_t maxNumReductions, size_t* pnumReductions, CMR_SUBMAT** preducedSubmatrix, CMR_SUBMAT** pviolatorSubmatrix, CMR_SEPA** pseparation, CMR_SP_STATISTICS* stats, double timeLimit)
    CMR_ERROR CMRdecomposeTernarySeriesParallel(CMR* cmr, CMR_CHRMAT* matrix, bool* pisSeriesParallel, CMR_SP_REDUCTION* reductions, size_t maxNumReductions, size_t* pnumReductions, CMR_SUBMAT** preducedSubmatrix, CMR_SUBMAT** pviolatorSubmatrix, CMR_SEPA** pseparation, CMR_SP_STATISTICS* stats, double timeLimit)

cdef extern from "cmr/network.h":

    ctypedef struct CMR_NETWORK_STATISTICS:
        size_t totalCount
        double totalTime
        CMR_CAMION_STATISTICS camion
        CMR_GRAPHIC_STATISTICS graphic

    CMR_ERROR CMRstatsNetworkInit(CMR_NETWORK_STATISTICS* stats)
    # CMR_ERROR CMRstatsNetworkPrint(FILE* stream, CMR_NETWORK_STATISTICS* stats, const char* prefix)
    CMR_ERROR CMRcomputeNetworkMatrix(CMR* cmr, CMR_GRAPH* digraph, CMR_CHRMAT** pmatrix, CMR_CHRMAT** ptranspose, bool* arcsReversed, int numForestArcs, CMR_GRAPH_EDGE* forestArcs, int numCoforestArcs, CMR_GRAPH_EDGE* coforestArcs, bool* pisCorrectForest)
    CMR_ERROR CMRtestNetworkMatrix(CMR* cmr, CMR_CHRMAT* matrix, bool* pisNetwork, CMR_GRAPH** pdigraph, CMR_GRAPH_EDGE** pforestArcs, CMR_GRAPH_EDGE** pcoforestArcs, bool** parcsReversed, CMR_SUBMAT** psubmatrix, CMR_NETWORK_STATISTICS* stats, double timeLimit)
    CMR_ERROR CMRtestConetworkMatrix(CMR* cmr, CMR_CHRMAT* matrix, bool* pisConetwork, CMR_GRAPH** pdigraph, CMR_GRAPH_EDGE** pforestArcs, CMR_GRAPH_EDGE** pcoforestArcs, bool** parcsReversed, CMR_SUBMAT** psubmatrix, CMR_NETWORK_STATISTICS* stats, double timeLimit)

cdef extern from "cmr/dec.h":

    ctypedef struct CMR_DEC

    ctypedef int CMR_DEC_TYPE

    const int CMR_DEC_IRREGULAR
    const int CMR_DEC_UNKNOWN
    const int CMR_DEC_ONE_SUM
    const int CMR_DEC_TWO_SUM
    const int CMR_DEC_THREE_SUM
    const int CMR_DEC_GRAPHIC
    const int CMR_DEC_COGRAPHIC
    const int CMR_DEC_PLANAR
    const int CMR_DEC_SERIES_PARALLEL
    const int CMR_DEC_SPECIAL_R10
    const int CMR_DEC_SPECIAL_FANO
    const int CMR_DEC_SPECIAL_FANO_DUAL
    const int CMR_DEC_SPECIAL_K_5
    const int CMR_DEC_SPECIAL_K_5_DUAL
    const int CMR_DEC_SPECIAL_K_3_3
    const int CMR_DEC_SPECIAL_K_3_3_DUAL

    ctypedef int CMR_DEC_FLAGS

    const int CMR_DEC_MASK_REPRESENTATION
    const int CMR_DEC_IS_GRAPHIC
    const int CMR_DEC_IS_COGRAPHIC
    const int CMR_DEC_IS_REGULAR
    const int CMR_DEC_HAS_LOWER_LEFT_NONZEROS
    const int CMR_DEC_HAS_UPPER_RIGHT_NONZEROS

    CMR_ERROR CMRdecFree(CMR* cmr, CMR_DEC** pdec)
    bint CMRdecHasMatrix(CMR_DEC* dec)
    bint CMRdecHasTranspose(CMR_DEC* dec)
    CMR_CHRMAT* CMRdecGetMatrix(CMR_DEC* dec)
    CMR_CHRMAT* CMRdecGetTranspose(CMR_DEC* dec)
    size_t CMRdecNumChildren(CMR_DEC* dec)
    CMR_DEC* CMRdecChild(CMR_DEC* dec, size_t childIndex)
    int CMRdecIsSum(CMR_DEC* dec, bool* plowerLeftNonzeros, bool* pupperRightNonzeros)
    CMR_DEC_TYPE CMRdecIsSpecialLeaf(CMR_DEC* dec, int* prepresentationMatrix)
    bint CMRdecIsGraphicLeaf(CMR_DEC* dec)
    bint CMRdecIsCographicLeaf(CMR_DEC* dec)
    bint CMRdecIsGraphic(CMR_DEC* dec)
    bint CMRdecIsCographic(CMR_DEC* dec)
    bint CMRdecIsRegular(CMR_DEC* dec)
    bint CMRdecIsSeriesParallelReduction(CMR_DEC* dec)
    bint CMRdecIsUnknown(CMR_DEC* dec)
    bint CMRdecNumRows(CMR_DEC* dec)
    CMR_ELEMENT* CMRdecRowElements(CMR_DEC* dec)
    size_t* CMRdecRowsParent(CMR_DEC* dec)
    bint CMRdecNumColumns(CMR_DEC* dec)
    CMR_ELEMENT* CMRdecColumnElements(CMR_DEC* dec)
    size_t* CMRdecColumnsParent(CMR_DEC* dec)
    # CMR_ERROR CMRdecPrint(CMR* cmr, CMR_DEC* dec, FILE* stream, size_t indent, bint printMatrices, bint printGraphs, bint printReductions)
    char* CMRdecConsistency(CMR_DEC* dec, bint recurse)
    CMR_GRAPH* CMRdecGraph(CMR_DEC* dec)
    CMR_GRAPH_EDGE* CMRdecGraphForest(CMR_DEC* dec)
    size_t CMRdecGraphSizeForest(CMR_DEC* dec)
    CMR_GRAPH_EDGE* CMRdecGraphCoforest(CMR_DEC* dec)
    size_t CMRdecGraphSizeCoforest(CMR_DEC* dec)
    bool* CMRdecGraphArcsReversed(CMR_DEC* dec)
    CMR_GRAPH* CMRdecCograph(CMR_DEC* dec)
    CMR_GRAPH_EDGE* CMRdecCographForest(CMR_DEC* dec)
    CMR_GRAPH_EDGE* CMRdecCographCoforest(CMR_DEC* dec)
    bool* CMRdecCographArcsReversed(CMR_DEC* dec)


cdef extern from "cmr/regular.h":

    ctypedef int CMR_DEC_CONSTRUCT

    const int CMR_DEC_CONSTRUCT_NONE
    const int CMR_DEC_CONSTRUCT_LEAVES
    const int CMR_DEC_CONSTRUCT_ALL

    ctypedef struct CMR_REGULAR_PARAMETERS:
        bint directGraphicness
        bint seriesParallel
        bint planarityCheck
        bint completeTree
        CMR_DEC_CONSTRUCT matrices
        CMR_DEC_CONSTRUCT transposes
        CMR_DEC_CONSTRUCT graphs

    CMR_ERROR CMRparamsRegularInit(CMR_REGULAR_PARAMETERS* params)

    ctypedef struct CMR_REGULAR_STATISTICS:
        size_t totalCount
        double totalTime
        CMR_SP_STATISTICS seriesParallel
        CMR_GRAPHIC_STATISTICS graphic
        CMR_NETWORK_STATISTICS network
        size_t sequenceExtensionCount
        double sequenceExtensionTime
        size_t sequenceGraphicCount
        double sequenceGraphicTime
        size_t enumerationCount
        double enumerationTime
        size_t enumerationCandidatesCount

    CMR_ERROR CMRstatsRegularInit(CMR_REGULAR_STATISTICS* stats)
    # CMR_ERROR CMRstatsRegularPrint(FILE* stream, CMR_REGULAR_STATISTICS* stats, const char* prefix)
    CMR_ERROR CMRtestBinaryRegular(CMR* cmr, CMR_CHRMAT* matrix, bint *pisRegular, CMR_DEC** pdec, CMR_MINOR** pminor, CMR_REGULAR_PARAMETERS* params, CMR_REGULAR_STATISTICS* stats, double timeLimit)


cdef extern from "cmr/tu.h":

    ctypedef struct CMR_TU_PARAMETERS:
        CMR_REGULAR_PARAMETERS regular

    CMR_ERROR CMRparamsTotalUnimodularityInit(CMR_TU_PARAMETERS* params)

    ctypedef struct CMR_TU_STATISTICS:
        size_t totalCount
        double totalTime
        CMR_CAMION_STATISTICS camion
        CMR_REGULAR_STATISTICS regular

    CMR_ERROR CMRstatsTotalUnimodularityInit(CMR_TU_STATISTICS* stats)
    # CMR_ERROR CMRstatsTotalUnimodularityPrint(FILE* stream, CMR_TU_STATISTICS* stats, const char* prefix)
    CMR_ERROR CMRtestTotalUnimodularity(CMR* cmr, CMR_CHRMAT* matrix, bool* pisTotallyUnimodular, CMR_DEC** pdec, CMR_SUBMAT** psubmatrix, CMR_TU_PARAMETERS* params, CMR_TU_STATISTICS* stats, double timeLimit)


cdef extern from "cmr/ctu.h":

    ctypedef struct CMR_CTU_STATISTICS:
        size_t totalCount
        double totalTime
        CMR_TU_STATISTICS tu

    CMR_ERROR CMRstatsComplementTotalUnimodularityInit(CMR_CTU_STATISTICS* stats)
    # CMR_ERROR CMRstatsComplementTotalUnimodularityPrint(FILE* stream, CMR_CTU_STATISTICS* stats, const char* prefix)
    CMR_ERROR CMRcomplementRowColumn(CMR* cmr, CMR_CHRMAT* matrix, size_t complementRow, size_t complementColumn, CMR_CHRMAT** presult)
    CMR_ERROR CMRtestComplementTotalUnimodularity(CMR* cmr, CMR_CHRMAT* matrix, bool* pisComplementTotallyUnimodular, size_t* pcomplementRow, size_t* pcomplementColumn, CMR_CTU_STATISTICS* stats)


# Our global CMR environment
cdef CMR *cmr

cdef CMR_CALL(CMR_ERROR _cmr_error)
