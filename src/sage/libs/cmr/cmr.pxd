# -*- python -*-
# distutils: libraries = cmr

# (progn (replace-regexp "/[*]\\(.\\|\n\\)*?[*]/" "" nil (point) (point-max)) (replace-regexp "[;{}]" ""  nil (point) (point-max)) (replace-regexp "CMR_EXPORT *" ""  nil (point) (point-max)) (replace-regexp "bool" "bint" nil (point) (point-max)))

cdef extern from "cmr/env.h":
    int CMR_OKAY
    int CMR_ERROR_INPUT
    int CMR_ERROR_OUTPUT
    int CMR_ERROR_MEMORY
    int CMR_ERROR_INVALID
    int CMR_ERROR_OVERFLOW
    int CMR_ERROR_TIMEOUT

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

    CMR_ERROR CMRtestUnimodularity(CMR* cmr, CMR_CHRMAT* matrix, bint* pisUnimodular)
    CMR_ERROR CMRtestStrongUnimodularity(CMR* cmr, CMR_CHRMAT* matrix, bint* pisStronglyUnimodular)
    CMR_ERROR CMRtestKmodularity(CMR* cmr, CMR_CHRMAT* matrix, bint* pisKmodular, size_t* pk)
    CMR_ERROR CMRtestStrongKmodularity(CMR* cmr, CMR_CHRMAT* matrix, bint* pisStronglyKmodular, size_t* pk)

cdef extern from "cmr/camion.h":

    ctypedef struct CMR_CAMION_STATISTICS:
        size_t totalCount
        double totalTime

    CMR_ERROR CMRstatsCamionInit(CMR_CAMION_STATISTICS* stats)
    # CMR_ERROR CMRstatsCamionPrint(FILE* stream, CMR_CAMION_STATISTICS* stats, const char* prefix)
    CMR_ERROR CMRtestCamionSigned(CMR* cmr, CMR_CHRMAT* matrix, bint* pisCamionSigned, CMR_SUBMAT** psubmatrix, CMR_CAMION_STATISTICS* stats, double timeLimit)
    CMR_ERROR CMRcomputeCamionSigned(CMR* cmr, CMR_CHRMAT* matrix, bint* pwasCamionSigned, CMR_SUBMAT** psubmatrix, CMR_CAMION_STATISTICS* stats, double timeLimit)

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
    CMR_ERROR CMRsepaCheckTernary(CMR* cmr, CMR_SEPA* sepa, CMR_CHRMAT* matrix, CMR_SUBMAT* submatrix, bint* pisTernary, CMR_SUBMAT** psubmatrix)
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
    CMR_ERROR CMRcomputeGraphicMatrix(CMR* cmr, CMR_GRAPH* graph, CMR_CHRMAT** pmatrix, CMR_CHRMAT** ptranspose, int numForestEdges, CMR_GRAPH_EDGE* forestEdges, int numCoforestEdges, CMR_GRAPH_EDGE* coforestEdges, bint* pisCorrectForest)
    CMR_ERROR CMRtestGraphicMatrix(CMR* cmr, CMR_CHRMAT* matrix, bint* pisGraphic, CMR_GRAPH** pgraph, CMR_GRAPH_EDGE** pforestEdges, CMR_GRAPH_EDGE** pcoforestEdges, CMR_SUBMAT** psubmatrix, CMR_GRAPHIC_STATISTICS* stats, double timeLimit)
    CMR_ERROR CMRtestCographicMatrix(CMR* cmr, CMR_CHRMAT* matrix, bint* pisCographic, CMR_GRAPH** pgraph, CMR_GRAPH_EDGE** pforestEdges, CMR_GRAPH_EDGE** pcoforestEdges, CMR_SUBMAT** psubmatrix, CMR_GRAPHIC_STATISTICS* stats, double timeLimit)
    CMR_ERROR CMRtestGraphicColumnSubmatrixGreedy(CMR* cmr, CMR_CHRMAT* transpose, size_t* orderedColumns, CMR_SUBMAT** psubmatrix)


# Our global CMR environment
cdef CMR *cmr

cdef CMR_CALL(CMR_ERROR _cmr_error)
