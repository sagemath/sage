/*
 * Copy comparison methods from type "src" to type "dst" if "dst" has
 * no comparison defined.
 *
 * This is implemented in pure C instead of Cython because it needs to
 * be different in Python 2 and Python 3.
 */
static CYTHON_INLINE void inherit_comparison(PyTypeObject* dst, const PyTypeObject* src)
{
    /* Do nothing if "dst" already has comparison defined */
#if (PY_MAJOR_VERSION < 3)
    if (dst->tp_compare) return;
#endif
    if (dst->tp_richcompare) return;

    /* Copy comparison method(s) */
#if (PY_MAJOR_VERSION < 3)
    dst->tp_compare = src->tp_compare;
#endif
    dst->tp_richcompare = src->tp_richcompare;
}
