# distutils: libraries = lrcalc

cdef extern from "lrcalc/hashtab.h":
    ctypedef struct hashtab:
        pass

    ctypedef struct hash_itr:
        pass

    ctypedef unsigned long hashkey_t
    ctypedef int (*cmp_t) (void* a, void* b)
    ctypedef hashkey_t (*hash_t) (void* a)

    hashtab* hash_new(cmp_t cm, hash_t hsh)
    void hash_free(hashtab *ht)

    void* hash_lookup(hashtab *ht, void *key)
    void* hash_insert(hashtab *ht, void *key, void *value)

    bint hash_good(hash_itr)
    void hash_first(hashtab* s, hash_itr itr)
    void hash_next(hash_itr itr)
    void* hash_key(hash_itr itr)
    void* hash_value(hash_itr itr)
    int hash_intvalue(hash_itr itr)

cdef extern from "lrcalc/vector.h":
    ctypedef struct vector:
        size_t length
        int* array

    vector* v_new(int length)
    void v_free(vector* v)
    void v_print(vector *v)
    int v_length(vector* v)
    int v_elem(vector* v, int i)

    ctypedef struct vecpair:
        vector *first
        vector *second

    vector* vp_first(vecpair* vp)
    vector* vp_second(vecpair* vp)

cdef extern from "lrcalc/list.h":
    cdef struct _list:
        void **array
        size_t allocated
        size_t length
    void l_free(_list *lst)

cdef extern from "lrcalc/symfcn.h":
    long long lrcoef_c "lrcoef"(vector* outer, vector* inner1, vector* inner2)
    hashtab* mult_c "mult"(vector *sh1, vector *sh2, int maxrows)
    hashtab* skew_c "skew"(vector *outer, vector *inner, int maxrows)
    hashtab* coprod_c "coprod"(vector *part, int all)
    void fusion_reduce_c "fusion_reduce"(hashtab* ht, int rows, int cols, int opt_zero)
    _list *quantum_reduce_c "quantum_reduce"(hashtab* ht, int rows, int col)

    ctypedef struct skewtab:
        vector *outer
        vector *inner
        vector *conts
        int maxrows
        vector *conjugate
        int rows
        int cols
        int matrix[1]

    skewtab *st_new(vector *outer, vector *inner, vector *conts, int maxrows)
    int st_next(skewtab *st)
    void st_print(skewtab *st)
    void st_free(skewtab *st)


cdef extern from "lrcalc/schublib.h":
    hashtab* mult_schubert_c "mult_schubert"(vector *sh1, vector *sh2, int rank)
