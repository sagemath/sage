#include <stdio.h>

#if PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION >= 9

static void _type_debug(PyTypeObject* tp)
{
    printf("Not implemented for CPython >= 3.9\n");
}

#else

#define HAVE_WEAKREFS(tp) (1)
#define HAVE_CLASS(tp) (1)
#define HAVE_ITER(tp) (1)
#define HAVE_RICHCOMPARE(tp) (1)
#define HAVE_INPLACEOPS(tp) (1)
#define HAVE_SEQUENCE_IN(tp) (1)
#define HAVE_GETCHARBUFFER(tp) (1)
#define HAVE_NEW_DIVISION(tp) (1)
#define HAVE_INDEX(tp) (1)
#define HAVE_NEWBUFFER(tp) (1)
#define HAVE_FINALIZE(tp) (tp->tp_flags & Py_TPFLAGS_HAVE_FINALIZE)

static void print_object(void* pyobj)
{
    if (pyobj == NULL)
    {
        printf("NULL\n");
        return;
    }

    PyObject* obj = (PyObject*)pyobj;

    if (PyTuple_Check(obj))
    {
        printf("%s:\n", Py_TYPE(obj)->tp_name);

        Py_ssize_t i, n = PyTuple_GET_SIZE(obj);
        for (i = 0; i < n; i++)
        {
            printf("    ");
            PyObject_Print(PyTuple_GET_ITEM(obj, i), stdout, 0);
            printf("\n");
        }
    }
    else if (PyDict_Check(obj))
    {
        printf("%s:\n", Py_TYPE(obj)->tp_name);

        Py_ssize_t pos = 0;
        PyObject* key;
        PyObject* value;
        while (PyDict_Next(obj, &pos, &key, &value))
        {
            printf("    ");
            PyObject_Print(key, stdout, 0);
            printf(": ");
            PyObject_Print(value, stdout, 0);
            printf("\n");
        }
    }
    else
    {
        PyObject_Print(obj, stdout, 0);
        printf("\n");
    }
}


#define pointer_check_constant(ptr, val) \
    if (ptr == (void*)(val)) \
    { \
        special = 0; \
        printf(#val "\n"); \
    }

#define subattr_pointer_value(sub, attr) \
    special = 1; \
    pointer_check_constant(tp->attr, NULL) \
    else pointer_check_constant(tp->attr, PyType_GenericAlloc) \
    else pointer_check_constant(tp->attr, PyType_GenericNew) \
    else pointer_check_constant(tp->attr, PyObject_Del) \
    else pointer_check_constant(tp->attr, PyObject_GC_Del) \
    else pointer_check_constant(tp->attr, PyObject_GenericGetAttr) \
    else pointer_check_constant(tp->attr, PyObject_GenericSetAttr) \
    else pointer_check_constant(tp->attr, _Py_HashPointer) \
    else pointer_check_constant(tp->attr, PyObject_HashNotImplemented) \
    else pointer_check_constant(tp->attr, subtype_traverse) \
    else pointer_check_constant(tp->attr, subtype_clear) \
    else pointer_check_constant(tp->attr, subtype_dealloc) \
    else \
    { \
        PyObject* mro = tp->tp_mro; \
        PyTypeObject* subtp; \
        Py_ssize_t i, n = PyTuple_GET_SIZE(mro); \
        for (i = n-1; i >= 0; i--) \
        { \
            subtp = (PyTypeObject*)PyTuple_GET_ITEM(mro, i); \
            if (subtp != tp && PyType_Check(subtp)) \
            { \
                if (subtp->sub && tp->attr == subtp->attr) \
                { \
                    special = 0; \
                    printf("== %s\n", subtp->tp_name); \
                    break; \
                } \
            } \
        } \
    } \
    if (special) printf("%p\n", tp->attr);

#define attr_pointer_value(attr) subattr_pointer_value(tp_name, attr);

#define attr_pointer(attr) \
    printf("  " #attr ": "); \
    attr_pointer_value(attr);
#define attr_pointer_meth(attr, method) \
    printf("  " #attr " (" method "): "); \
    attr_pointer_value(attr);
#define subattr_pointer(sub, attr) \
    printf("    " #attr ": "); \
    subattr_pointer_value(sub, sub->attr);
#define subattr_pointer_meth(sub, attr, method) \
    printf("    " #attr " (" method "): "); \
    subattr_pointer_value(sub, sub->attr);

#define attr_object(attr) \
    printf("  " #attr ": "); \
    print_object(tp->attr);
#define attr_object_meth(attr, method) \
    printf("  " #attr " (" method "): "); \
    print_object(tp->attr);

#define attr_flag(flag) \
    if (tp->tp_flags & Py_TPFLAGS_ ## flag) \
        printf("    " #flag "\n");


static void _type_debug(PyTypeObject* tp)
{
    int special;

    PyObject_Print((PyObject*)tp, stdout, 0);
    printf(" (%p)\n", tp);
    printf("  ob_refcnt: %ld\n", (long)Py_REFCNT(tp));
    printf("  ob_type: "); print_object(Py_TYPE(tp));
    printf("  tp_name: %s\n", tp->tp_name);
    printf("  tp_basicsize: %ld\n", (long)tp->tp_basicsize);
    printf("  tp_itemsize: %ld\n", (long)tp->tp_itemsize);
    printf("  tp_dictoffset: %ld\n", (long)tp->tp_dictoffset);
    if HAVE_WEAKREFS(tp)
    {
        printf("  tp_weaklistoffset: %ld\n", (long)tp->tp_weaklistoffset);
    }

    if HAVE_CLASS(tp)
    {
        attr_object_meth(tp_base, "__base__");
        attr_object_meth(tp_bases, "__bases__");
        attr_object_meth(tp_mro, "__mro__");
        attr_object_meth(tp_dict, "__dict__");
    }

    attr_pointer(tp_alloc);
    attr_pointer_meth(tp_new, "__new__");
    attr_pointer_meth(tp_init, "__init__");
    attr_pointer_meth(tp_dealloc, "__dealloc__");
    if (tp->tp_flags & Py_TPFLAGS_HEAPTYPE)
    {
        attr_pointer_meth(tp_del, "__del__");
    }
    #if PY_MAJOR_VERSION >= 3
    if HAVE_FINALIZE(tp)
    {
        attr_pointer_meth(tp_finalize, "__del__");
    }
    #endif
    attr_pointer(tp_free);

    attr_pointer_meth(tp_repr, "__repr__");
    attr_pointer(tp_print);
    attr_pointer_meth(tp_hash, "__hash__");
    attr_pointer_meth(tp_call, "__call__");
    attr_pointer_meth(tp_str, "__str__");
    #if PY_MAJOR_VERSION <= 2
        attr_pointer_meth(tp_compare, "cmp");
    #endif
    attr_pointer_meth(tp_richcompare, "__richcmp__");
    attr_pointer_meth(tp_getattr, "__getattribute__");
    attr_pointer_meth(tp_setattr, "__setattribute__");
    attr_pointer_meth(tp_getattro, "__getattribute__");
    attr_pointer_meth(tp_setattro, "__setattribute__");
    if HAVE_ITER(tp)
    {
        attr_pointer_meth(tp_iter, "__iter__");
        attr_pointer_meth(tp_iternext, "__next__");
    }
    if HAVE_CLASS(tp)
    {
        attr_pointer_meth(tp_descr_get, "__get__");
        attr_pointer_meth(tp_descr_set, "__set__");

        attr_object(tp_cache);
        attr_object(tp_weaklist);
    }
    if HAVE_RICHCOMPARE(tp)
    {
        attr_pointer(tp_traverse);
        attr_pointer(tp_clear);
    }
    if HAVE_CLASS(tp)
    {
        attr_pointer(tp_is_gc);
    }

    attr_pointer(tp_as_number);
    if (special)
    {
        subattr_pointer_meth(tp_as_number, nb_add, "__add__");
        subattr_pointer_meth(tp_as_number, nb_subtract, "__sub__");
        subattr_pointer_meth(tp_as_number, nb_multiply, "__mul__");
        #if PY_MAJOR_VERSION <= 2
            subattr_pointer_meth(tp_as_number, nb_divide, "__div__");
        #endif
        if HAVE_NEW_DIVISION(tp)
        {
            subattr_pointer_meth(tp_as_number, nb_floor_divide, "__floordiv__");
            subattr_pointer_meth(tp_as_number, nb_true_divide, "__truediv__");
        }
        subattr_pointer_meth(tp_as_number, nb_remainder, "__mod__");
        subattr_pointer_meth(tp_as_number, nb_divmod, "__divmod__");
        subattr_pointer_meth(tp_as_number, nb_power, "__pow__");
        subattr_pointer_meth(tp_as_number, nb_negative, "__neg__");
        subattr_pointer_meth(tp_as_number, nb_positive, "__pos__");
        subattr_pointer_meth(tp_as_number, nb_absolute, "__abs__");
        #if PY_MAJOR_VERSION <= 2
            subattr_pointer_meth(tp_as_number, nb_nonzero, "__nonzero__");
        #else
            subattr_pointer_meth(tp_as_number, nb_bool, "__bool__");
        #endif
        subattr_pointer_meth(tp_as_number, nb_invert, "__invert__");
        subattr_pointer_meth(tp_as_number, nb_lshift, "__lshift__");
        subattr_pointer_meth(tp_as_number, nb_rshift, "__rshift__");
        subattr_pointer_meth(tp_as_number, nb_and, "__and__");
        subattr_pointer_meth(tp_as_number, nb_or, "__or__");
        subattr_pointer_meth(tp_as_number, nb_xor, "__xor__");
        subattr_pointer_meth(tp_as_number, nb_int, "__int__");
        #if PY_MAJOR_VERSION <= 2
            subattr_pointer_meth(tp_as_number, nb_long, "__long__");
        #endif
        if HAVE_INDEX(tp)
        {
            subattr_pointer_meth(tp_as_number, nb_index, "__index__");
        }
        subattr_pointer_meth(tp_as_number, nb_float, "__float__");
        #if PY_MAJOR_VERSION <= 2
            subattr_pointer_meth(tp_as_number, nb_oct, "__oct__");
            subattr_pointer_meth(tp_as_number, nb_hex, "__hex__");
            subattr_pointer(tp_as_number, nb_coerce);
        #endif

        if HAVE_INPLACEOPS(tp)
        {
            subattr_pointer_meth(tp_as_number, nb_inplace_add, "__iadd__");
            subattr_pointer_meth(tp_as_number, nb_inplace_subtract, "__isub__");
            subattr_pointer_meth(tp_as_number, nb_inplace_multiply, "__imul__");
            #if PY_MAJOR_VERSION <= 2
                subattr_pointer_meth(tp_as_number, nb_inplace_divide, "__idiv__");
            #endif
            if HAVE_NEW_DIVISION(tp)
            {
                subattr_pointer_meth(tp_as_number, nb_inplace_floor_divide, "__ifloordiv__");
                subattr_pointer_meth(tp_as_number, nb_inplace_true_divide, "__itruediv__");
            }
            subattr_pointer_meth(tp_as_number, nb_inplace_remainder, "__imod__");
            subattr_pointer_meth(tp_as_number, nb_inplace_power, "__ipow__");
            subattr_pointer_meth(tp_as_number, nb_inplace_lshift, "__ilshift__");
            subattr_pointer_meth(tp_as_number, nb_inplace_rshift, "__irshift__");
            subattr_pointer_meth(tp_as_number, nb_inplace_and, "__iand__");
            subattr_pointer_meth(tp_as_number, nb_inplace_or, "__ior__");
            subattr_pointer_meth(tp_as_number, nb_inplace_xor, "__ixor__");
        }
    }

    attr_pointer(tp_as_sequence);
    if (special)
    {
        subattr_pointer_meth(tp_as_sequence, sq_length, "__len__");
        subattr_pointer_meth(tp_as_sequence, sq_concat, "__add__");
        if HAVE_INPLACEOPS(tp)
        {
            subattr_pointer_meth(tp_as_sequence, sq_inplace_concat, "__iadd__");
        }
        subattr_pointer_meth(tp_as_sequence, sq_repeat, "__mul__");
        if HAVE_INPLACEOPS(tp)
        {
            subattr_pointer_meth(tp_as_sequence, sq_inplace_repeat, "__imul__");
        }
        subattr_pointer_meth(tp_as_sequence, sq_item, "__getitem__");
        subattr_pointer_meth(tp_as_sequence, sq_ass_item, "__setitem__");
        if HAVE_SEQUENCE_IN(tp)
        {
            subattr_pointer_meth(tp_as_sequence, sq_contains, "__contains__");
        }
    }

    attr_pointer(tp_as_mapping);
    if (special)
    {
        subattr_pointer_meth(tp_as_mapping, mp_length, "__len__");
        subattr_pointer_meth(tp_as_mapping, mp_subscript, "__getitem__");
        subattr_pointer_meth(tp_as_mapping, mp_ass_subscript, "__setitem__");
    }

    attr_pointer(tp_as_buffer);
    if (special)
    {
        #if PY_MAJOR_VERSION <= 2
            subattr_pointer(tp_as_buffer, bf_getreadbuffer);
            subattr_pointer(tp_as_buffer, bf_getwritebuffer);
            subattr_pointer(tp_as_buffer, bf_getsegcount);
            if HAVE_GETCHARBUFFER(tp)
            {
                subattr_pointer(tp_as_buffer, bf_getcharbuffer);
            }
        #endif
        if HAVE_NEWBUFFER(tp)
        {
            subattr_pointer_meth(tp_as_buffer, bf_getbuffer, "__getbuffer__");
            subattr_pointer_meth(tp_as_buffer, bf_releasebuffer, "__releasebuffer__");
        }
    }

    printf("  tp_flags:\n");
    attr_flag(HEAPTYPE);
    attr_flag(BASETYPE);
    attr_flag(READY);
    attr_flag(READYING);
    attr_flag(HAVE_GC);
    #if PY_MAJOR_VERSION <= 2
        attr_flag(CHECKTYPES);
    #endif
    attr_flag(HAVE_VERSION_TAG);
    attr_flag(VALID_VERSION_TAG);
    attr_flag(IS_ABSTRACT);
    if (tp->tp_flags & Py_TPFLAGS_HAVE_VERSION_TAG)
    {
        printf("  tp_version_tag: %lu\n", (unsigned long)tp->tp_version_tag);
    }
}

#endif
