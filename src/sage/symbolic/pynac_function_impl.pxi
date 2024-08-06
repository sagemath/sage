cpdef call_registered_function(unsigned serial,
                               int nargs,
                               list args,
                               bint hold,
                               bint allow_numeric_result,
                               result_parent):
    r"""
    Call a function registered with Pynac (GiNaC).

    INPUT:

    - ``serial`` -- serial number of the function

    - ``nargs`` -- declared number of args (0 is variadic)

    - ``args`` -- list of arguments to pass to the function;
      each must be an :class:`Expression`

    - ``hold`` -- whether to leave the call unevaluated

    - ``allow_numeric_result`` -- if ``True``, keep numeric results numeric;
      if ``False``, make all results symbolic expressions

    - ``result_parent`` -- an instance of :class:`SymbolicRing`

    EXAMPLES::

        sage: from sage.symbolic.expression import find_registered_function, call_registered_function
        sage: s_arctan = find_registered_function('arctan', 1)
        sage: call_registered_function(s_arctan, 1, [SR(1)], False, True, SR)
        1/4*pi
        sage: call_registered_function(s_arctan, 1, [SR(1)], True, True, SR)
        arctan(1)
        sage: call_registered_function(s_arctan, 1, [SR(0)], False, True, SR)
        0
        sage: call_registered_function(s_arctan, 1, [SR(0)], False, True, SR).parent()
        Integer Ring
        sage: call_registered_function(s_arctan, 1, [SR(0)], False, False, SR).parent()
        Symbolic Ring
    """
    cdef Py_ssize_t i
    cdef GEx res
    cdef GExVector vec
    if nargs == 0 or nargs > 3:
        for i in range(len(args)):
            vec.push_back((<Expression>args[i])._gobj)
        res = g_function_evalv(serial, vec, hold)
    elif nargs == 1:
        res = g_function_eval1(serial,
                (<Expression>args[0])._gobj, hold)
    elif nargs == 2:
        res = g_function_eval2(serial, (<Expression>args[0])._gobj,
                (<Expression>args[1])._gobj, hold)
    elif nargs == 3:
        res = g_function_eval3(serial,
                (<Expression>args[0])._gobj, (<Expression>args[1])._gobj,
                (<Expression>args[2])._gobj, hold)

    if allow_numeric_result and is_a_numeric(res):
        return py_object_from_numeric(res)

    return new_Expression_from_GEx(result_parent, res)


cpdef unsigned find_registered_function(name, int nargs) except -1:
    r"""
    Look up a function registered with Pynac (GiNaC).

    Raise a :exc:`ValueError` if the function is not registered.

    OUTPUT: serial number of the function, for use in :func:`call_registered_function`

    EXAMPLES::

        sage: from sage.symbolic.expression import find_registered_function
        sage: find_registered_function('arctan', 1)  # random
        19
        sage: find_registered_function('archenemy', 1)
        Traceback (most recent call last):
        ...
        ValueError: cannot find GiNaC function with name archenemy and 1 arguments
    """
    try:
        return find_function(str_to_bytes(name), nargs)
    except RuntimeError as err:
        raise ValueError("cannot find GiNaC function with name %s and %s arguments" % (name, nargs))


cpdef unsigned register_or_update_function(self, name, latex_name, int nargs,
                                           evalf_params_first, bint update) except -1:
    r"""
    Register the function ``self`` with Pynac (GiNaC).

    OUTPUT: serial number of the function, for use in :func:`call_registered_function`

    EXAMPLES::

        sage: from sage.symbolic.function import BuiltinFunction
        sage: class Archosaurian(BuiltinFunction):
        ....:     def __init__(self):
        ....:         BuiltinFunction.__init__(self, 'archsaur', nargs=1)
        ....:     def _eval_(self, x):
        ....:         return x * exp(x)
        sage: archsaur = Archosaurian()  # indirect doctest
        sage: archsaur(2)
        2*e^2
    """
    cdef GFunctionOpt opt
    cdef unsigned serial

    if update:
        serial = find_registered_function(name, nargs)
        opt = g_registered_functions().index(serial)
    else:
        opt = g_function_options_args(str_to_bytes(name), nargs)

    if hasattr(self, '_eval_'):
        opt.eval_func(self)

    if not evalf_params_first:
        opt.do_not_evalf_params()

    if hasattr(self, '_subs_'):
        opt.subs_func(self)

    if hasattr(self, '_evalf_'):
        opt.evalf_func(self)

    if hasattr(self, '_conjugate_'):
        opt.conjugate_func(self)

    if hasattr(self, '_real_part_'):
        opt.real_part_func(self)

    if hasattr(self, '_imag_part_'):
        opt.imag_part_func(self)

    if hasattr(self, '_derivative_'):
        opt.derivative_func(self)

    if hasattr(self, '_tderivative_'):
        opt.do_not_apply_chain_rule()
        opt.derivative_func(self)

    if hasattr(self, '_power_'):
        opt.power_func(self)

    if hasattr(self, '_series_'):
        opt.series_func(self)

    # custom print functions are called from python
    # so we don't register them with the ginac function_options object

    if latex_name:
        opt.latex_name(str_to_bytes(latex_name))

    if not update:
        serial = g_register_new(opt)

    g_foptions_assign(g_registered_functions().index(serial), opt)

    sfunction_serial_dict[serial] = self

    return serial


# we keep a database of symbolic functions initialized in a session
# this also makes the .operator() method of symbolic expressions work
cdef dict sfunction_serial_dict = {}


cpdef get_sfunction_from_serial(unsigned int serial):
    """
    Return an already created :class:`SymbolicFunction` given the serial.

    These are stored in the dictionary ``sfunction_serial_dict``.

    EXAMPLES::

        sage: from sage.symbolic.function import get_sfunction_from_serial
        sage: get_sfunction_from_serial(65) #random
        f
    """
    global sfunction_serial_dict
    return sfunction_serial_dict.get(serial)


cpdef get_sfunction_from_hash(long myhash):
    """
    Return an already created :class:`SymbolicFunction` given the hash.

    EXAMPLES::

        sage: from sage.symbolic.expression import get_sfunction_from_hash
        sage: get_sfunction_from_hash(1)  # random
    """
    for sfunc in sfunction_serial_dict.itervalues():
        if isinstance(sfunc, SymbolicFunction) and \
                myhash == (<SymbolicFunction>sfunc)._hash_():
            # found one
            return sfunc
    return None
