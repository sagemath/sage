# Cython declarations for the HiGHS C API.
#
# This wraps the public header highs/interfaces/highs_c_api.h.
# The HiGHS headers rely on relative includes within the "highs" include tree,
# so consumers must compile with -I$includedir/highs (typically provided by
# pkg-config via `highs.get_compile_args()` in Meson).

cdef extern from "highs/interfaces/highs_c_api.h":
    ctypedef int HighsInt

    # Status constants
    const HighsInt kHighsStatusError
    const HighsInt kHighsStatusOk
    const HighsInt kHighsStatusWarning

    # Variable type constants
    const HighsInt kHighsVarTypeContinuous
    const HighsInt kHighsVarTypeInteger
    const HighsInt kHighsVarTypeSemiContinuous
    const HighsInt kHighsVarTypeSemiInteger

    # Objective sense constants
    const HighsInt kHighsObjSenseMinimize
    const HighsInt kHighsObjSenseMaximize

    # Model status constants
    const HighsInt kHighsModelStatusNotset
    const HighsInt kHighsModelStatusLoadError
    const HighsInt kHighsModelStatusModelError
    const HighsInt kHighsModelStatusOptimal
    const HighsInt kHighsModelStatusInfeasible
    const HighsInt kHighsModelStatusUnboundedOrInfeasible
    const HighsInt kHighsModelStatusUnbounded
    const HighsInt kHighsModelStatusTimeLimit
    const HighsInt kHighsModelStatusIterationLimit
    const HighsInt kHighsModelStatusInterrupt
    const HighsInt kHighsModelStatusMemoryLimit

    # Basis status constants
    const HighsInt kHighsBasisStatusLower
    const HighsInt kHighsBasisStatusBasic
    const HighsInt kHighsBasisStatusUpper
    const HighsInt kHighsBasisStatusZero
    const HighsInt kHighsBasisStatusNonbasic

    # Core functions
    void* Highs_create() nogil
    void Highs_destroy(void* highs) nogil
    HighsInt Highs_run(void* highs) nogil

    # Model building
    HighsInt Highs_addCol(
        void* highs,
        double cost,
        double lower,
        double upper,
        HighsInt num_nz,
        const HighsInt* index,
        const double* value,
    ) nogil
    HighsInt Highs_addRow(
        void* highs,
        double lower,
        double upper,
        HighsInt num_nz,
        const HighsInt* index,
        const double* value,
    ) nogil
    HighsInt Highs_addCols(
        void* highs,
        HighsInt num_new_col,
        const double* costs,
        const double* lower,
        const double* upper,
        HighsInt num_new_nz,
        const HighsInt* starts,
        const HighsInt* indices,
        const double* values,
    ) nogil
    HighsInt Highs_addRows(
        void* highs,
        HighsInt num_new_row,
        const double* lower,
        const double* upper,
        HighsInt num_new_nz,
        const HighsInt* starts,
        const HighsInt* indices,
        const double* values,
    ) nogil

    # Objective
    HighsInt Highs_changeObjectiveSense(void* highs, HighsInt sense) nogil
    HighsInt Highs_changeColCost(void* highs, HighsInt col, double cost) nogil
    HighsInt Highs_changeObjectiveOffset(void* highs, double offset) nogil

    # Bounds
    HighsInt Highs_changeColBounds(void* highs, HighsInt col, double lower, double upper) nogil
    HighsInt Highs_changeRowBounds(void* highs, HighsInt row, double lower, double upper) nogil

    # Integrality
    HighsInt Highs_changeColIntegrality(void* highs, HighsInt col, HighsInt integrality) nogil
    HighsInt Highs_getColIntegrality(const void* highs, HighsInt col, HighsInt* integrality) nogil

    # Solution queries
    HighsInt Highs_getSolution(
        void* highs,
        double* col_value,
        double* col_dual,
        double* row_value,
        double* row_dual,
    ) nogil
    HighsInt Highs_getBasis(void* highs, HighsInt* col_status, HighsInt* row_status) nogil

    # Info queries
    HighsInt Highs_getModelStatus(const void* highs) nogil
    HighsInt Highs_getNumCol(const void* highs) nogil
    HighsInt Highs_getNumRow(const void* highs) nogil
    HighsInt Highs_getNumNz(const void* highs) nogil
    double Highs_getInfinity(const void* highs) nogil
    double Highs_getObjectiveValue(const void* highs) nogil
    HighsInt Highs_getDoubleInfoValue(const void* highs, const char* info, double* value) nogil
    HighsInt Highs_getIntInfoValue(const void* highs, const char* info, HighsInt* value) nogil
    HighsInt Highs_getObjectiveSense(const void* highs, HighsInt* sense) nogil
    HighsInt Highs_getObjectiveOffset(const void* highs, double* offset) nogil

    # Get columns/rows data
    HighsInt Highs_getColsByRange(
        const void* highs,
        HighsInt from_col,
        HighsInt to_col,
        HighsInt* num_col,
        double* costs,
        double* lower,
        double* upper,
        HighsInt* num_nz,
        HighsInt* matrix_start,
        HighsInt* matrix_index,
        double* matrix_value,
    ) nogil
    HighsInt Highs_getRowsByRange(
        const void* highs,
        HighsInt from_row,
        HighsInt to_row,
        HighsInt* num_row,
        double* lower,
        double* upper,
        HighsInt* num_nz,
        HighsInt* matrix_start,
        HighsInt* matrix_index,
        double* matrix_value,
    ) nogil

    # Name functions
    HighsInt Highs_passColName(void* highs, HighsInt col, const char* name) nogil
    HighsInt Highs_passRowName(void* highs, HighsInt row, const char* name) nogil

    # Options
    HighsInt Highs_setBoolOptionValue(void* highs, const char* option, HighsInt value) nogil
    HighsInt Highs_setIntOptionValue(void* highs, const char* option, HighsInt value) nogil
    HighsInt Highs_setDoubleOptionValue(void* highs, const char* option, double value) nogil
    HighsInt Highs_setStringOptionValue(void* highs, const char* option, const char* value) nogil
    HighsInt Highs_getBoolOptionValue(const void* highs, const char* option, HighsInt* value) nogil
    HighsInt Highs_getIntOptionValue(const void* highs, const char* option, HighsInt* value) nogil
    HighsInt Highs_getDoubleOptionValue(const void* highs, const char* option, double* value) nogil
    HighsInt Highs_getStringOptionValue(const void* highs, const char* option, char* value) nogil

    # Write/read model
    HighsInt Highs_writeModel(void* highs, const char* filename) nogil
    HighsInt Highs_readModel(void* highs, const char* filename) nogil

    # Basis functions
    HighsInt Highs_setBasic(void* highs, HighsInt* col_status, HighsInt* row_status) nogil
    HighsInt Highs_setBasis(void* highs, const HighsInt* col_status, const HighsInt* row_status) nogil

    # Delete functions
    HighsInt Highs_deleteRowsByRange(void* highs, HighsInt from_row, HighsInt to_row) nogil
    HighsInt Highs_deleteColsByRange(void* highs, HighsInt from_col, HighsInt to_col) nogil
    HighsInt Highs_deleteRowsBySet(void* highs, HighsInt num_set_entries, const HighsInt* set) nogil
    HighsInt Highs_deleteColsBySet(void* highs, HighsInt num_set_entries, const HighsInt* set) nogil

    # Change coefficient
    HighsInt Highs_changeCoeff(void* highs, HighsInt row, HighsInt col, double value) nogil

    # Names (not declared nogil since they work with strings)
    HighsInt Highs_getColName(const void* highs, HighsInt col, char* name)
    HighsInt Highs_getRowName(const void* highs, HighsInt row, char* name)
