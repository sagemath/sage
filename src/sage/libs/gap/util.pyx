"""
Utility functions for libGAP
"""

###############################################################################
#       Copyright (C) 2012, Volker Braun <vbraun.name@gmail.com>
#
#   Distributed under the terms of the GNU General Public License (GPL)
#   as published by the Free Software Foundation; either version 2 of
#   the License, or (at your option) any later version.
#                   http://www.gnu.org/licenses/
###############################################################################

from sage.env import SAGE_LOCAL
from sage.misc.misc import is_64_bit
from libc.stdint cimport uintptr_t
from element cimport *



############################################################################
### Hooking into the GAP memory management #################################
############################################################################

cdef class ObjWrapper(object):
    """
    Wrapper for GAP master pointers

    EXAMPLES::

        sage: from sage.libs.gap.util import ObjWrapper
        sage: x = ObjWrapper()
        sage: y = ObjWrapper()
        sage: x == y
        True
    """

    def __richcmp__(ObjWrapper self, ObjWrapper other, int op):
        r"""
        Comparison wrapped libGAP_Obj.

        INPUT:

        - ``lhs``, ``rhs`` -- :class:`ObjWrapper`.

        - ``op`` -- integer. The comparison operation to be performed.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: from sage.libs.gap.util import ObjWrapper
            sage: x = ObjWrapper()
            sage: y = ObjWrapper()
            sage: x == y
            True
        """
        cdef result
        cdef libGAP_Obj self_value = self.value
        cdef libGAP_Obj other_value = other.value
        if op==0:      # <   0
            return self_value < other_value
        elif op==1:    # <=  1
            return self_value <= other_value
        elif op==2:    # ==  2
            return self_value == other_value
        elif op==4:    # >   4
            return self_value > other_value
        elif op==5:    # >=  5
            return self_value >= other_value
        elif op==3:    # !=  3
            return self_value != other_value
        else:
            assert False  # unreachable

    def __hash__(self):
        """
        Return a hash value

        EXAMPLES::

            sage: from sage.libs.gap.util import ObjWrapper
            sage: x = ObjWrapper()
            sage: hash(x)
            0
        """
        return <int>(self.value)


cdef ObjWrapper wrap_obj(libGAP_Obj obj):
    """
    Constructor function for :class:`ObjWrapper`
    """
    cdef ObjWrapper result = ObjWrapper.__new__(ObjWrapper)
    result.value = obj
    return result


owned_objects_refcount = dict()

cpdef get_owned_objects():
    """
    Helper to access the refcount dictionary from Python code
    """
    return owned_objects_refcount


cdef void reference_obj(libGAP_Obj obj):
    """
    Reference ``obj``
    """
    cdef ObjWrapper wrapped = wrap_obj(obj)
    global owned_objects_refcount
    if wrapped in owned_objects_refcount:
        owned_objects_refcount[wrapped] += 1
    else:
        owned_objects_refcount[wrapped] = 1


cdef void dereference_obj(libGAP_Obj obj):
    """
    Reference ``obj``
    """
    cdef ObjWrapper wrapped = wrap_obj(obj)
    global owned_objects_refcount
    refcount = owned_objects_refcount.pop(wrapped)
    if refcount > 1:
        owned_objects_refcount[wrapped] = refcount - 1


cdef void gasman_callback():
    """
    Callback before each GAP garbage collection
    """
    global owned_objects_refcount
    for obj in owned_objects_refcount.iterkeys():
        libGAP_MARK_BAG((<ObjWrapper>obj).value)





############################################################################
### Initialization of libGAP ###############################################
############################################################################

def gap_root():
    """
    Find the location of the GAP root install which is stored in the gap
    startup script.

    EXAMPLES::

        sage: from sage.libs.gap.util import gap_root
        sage: gap_root()   # random output
        '/home/vbraun/opt/sage-5.3.rc0/local/gap/latest'
    """
    import os.path
    gapdir = os.path.join(SAGE_LOCAL, 'gap', 'latest')
    if os.path.exists(gapdir):
        return gapdir
    print 'The gap-4.5.5.spkg (or later) seems to be not installed!'
    gap_sh = open(os.path.join(SAGE_LOCAL, 'bin', 'gap')).read().splitlines()
    gapdir = filter(lambda dir:dir.strip().startswith('GAP_DIR'), gap_sh)[0]
    gapdir = gapdir.split('"')[1]
    gapdir = gapdir.replace('$SAGE_LOCAL', SAGE_LOCAL)
    return gapdir


cdef initialize():
    """
    Initialize the GAP library, if it hasn't already been
    initialized.  It is safe to call this multiple times.

    TESTS::

        sage: libgap(123)   # indirect doctest
        123
    """
    global _gap_is_initialized
    if _gap_is_initialized: return

    # Define argv and environ variables, which we will pass in to
    # initialize GAP. Note that we must pass define the memory pool
    # size!
    cdef char* argv[12]
    argv[0] = "sage"
    argv[1] = "-l"
    s = gap_root()
    argv[2] = s

    from sage.interfaces.gap import _get_gap_memory_pool_size_MB
    memory_pool = _get_gap_memory_pool_size_MB()
    argv[3] = "-o"
    argv[4] = memory_pool
    argv[5] = "-s"
    argv[6] = memory_pool

    argv[7] = "-m"
    argv[8] = "64m"

    argv[9] = "-q"    # no prompt!
    argv[10] = "-T"    # no debug loop
    argv[11] = NULL
    cdef int argc = 11   # argv[argc] must be NULL

    # Initialize GAP and capture any error messages
    # The initialization just prints error and does not use the error handler
    libgap_start_interaction('')
    libgap_initialize(argc, argv)
    gap_error_msg = str(libgap_get_output())
    libgap_finish_interaction()
    if gap_error_msg:
        raise RuntimeError('libGAP initialization failed\n' + gap_error_msg)

    # The error handler is called if a GAP evaluation fails, e.g. 1/0
    libgap_set_error_handler(&error_handler)

    # Prepare global GAP variable to hold temporary GAP objects
    global reference_holder
    libgap_enter()
    reference_holder = libGAP_GVarName("$SAGE_libgap_reference_holder")
    libgap_exit()

    # Finished!
    _gap_is_initialized = True



############################################################################
### Evaluate string in GAP #################################################
############################################################################

cdef libGAP_Obj gap_eval(str gap_string) except? NULL:
    r"""
    Evaluate a string in GAP.

    INPUT:

    - ``gap_string`` -- string. A valid statement in GAP.

    OUTPUT:

    The resulting GAP object or NULL+Python Exception in case of error.

    EXAMPLES::

        sage: libgap.eval('if 4>3 then\nPrint("hi");\nfi')
        NULL
        sage: libgap.eval('1+1')   # testing that we have sucessfully recovered
        2

        sage: libgap.eval('if 4>3 thenPrint("hi");\nfi')
        Traceback (most recent call last):
        ...
        ValueError: libGAP: Syntax error: then expected
        if 4>3 thenPrint("hi");
        fi;
                       ^
        sage: libgap.eval('1+1')   # testing that we have sucessfully recovered
        2
    """
    initialize()
    cdef libGAP_ExecStatus status

    cmd = gap_string + ';\n'
    try:
        libgap_enter()
        libgap_start_interaction(cmd)
        try:
            sig_on()
            status = libGAP_ReadEvalCommand(libGAP_BottomLVars)
            if status != libGAP_STATUS_END:
                libgap_call_error_handler()
            sig_off()
        except RuntimeError, msg:
            raise ValueError('libGAP: '+str(msg).strip())

        if libGAP_Symbol != libGAP_S_SEMICOLON:
            raise ValueError('did not end with semicolon')
        libGAP_GetSymbol()
        if libGAP_Symbol != libGAP_S_EOF:
            raise ValueError('can only evaluate a single statement')

    finally:
        libgap_finish_interaction()
        libgap_exit()

    if libGAP_ReadEvalResult != NULL:
        libgap_enter()
        libGAP_AssGVar(libGAP_Last3, libGAP_VAL_GVAR(libGAP_Last2))
        libGAP_AssGVar(libGAP_Last2, libGAP_VAL_GVAR(libGAP_Last))
        libGAP_AssGVar(libGAP_Last, libGAP_ReadEvalResult)
        libgap_exit()

    return libGAP_ReadEvalResult   # may be NULL, thats ok


############################################################################
### Helper to protect temporary objects from deletion ######################
############################################################################

cdef void hold_reference(libGAP_Obj obj):
    """
    Hold a reference (inside the GAP kernel) to obj

    This ensures that the GAP garbage collector does not delete
    ``obj``. This works by assigning it to a global variable. This is
    very simple, but you can't use it to keep two objects alive. Be
    careful.
    """
    libgap_enter()
    global reference_holder
    libGAP_AssGVar(reference_holder, obj)
    libgap_exit()


############################################################################
### Error handler ##########################################################
############################################################################

include 'sage/ext/interrupt.pxi'
from cpython.exc cimport PyErr_SetObject

cdef void error_handler(char* msg):
    """
    The libgap error handler

    We call ``sig_error()`` which causes us to jump back to the Sage
    signal handler. Since we wrap libGAP C calls in ``sig_on`` /
    ``sig_off`` blocks, this then jumps back to the ``sig_on`` where
    the ``RuntimeError`` we raise here will be seen.
    """
    msg_py = msg
    msg_py = msg_py.replace('For debugging hints type ?Recovery from NoMethodFound\n', '')
    PyErr_SetObject(RuntimeError, msg_py)
    sig_error()


############################################################################
### Debug functions ########################################################
############################################################################

cdef inline void DEBUG_CHECK(libGAP_Obj obj):
    """
    Check that ``obj`` is valid.

    This function is only useful for debugging.
    """
    libgap_enter()
    libGAP_CheckMasterPointers()
    libgap_exit()
    if obj == NULL:
        print 'DEBUG_CHECK: Null pointer!'




cpdef memory_usage():
    """
    Return information about the memory useage.

    See :meth:`~sage.libs.gap.libgap.Gap.mem` for details.
    """
    cdef size_t SizeMptrsArea = libGAP_OldBags - libGAP_MptrBags
    cdef size_t SizeOldBagsArea = libGAP_YoungBags - libGAP_OldBags
    cdef size_t SizeYoungBagsArea = libGAP_AllocBags - libGAP_YoungBags
    cdef size_t SizeAllocationArea = libGAP_StopBags - libGAP_AllocBags
    cdef size_t SizeUnavailableArea = libGAP_EndBags - libGAP_StopBags
    return (SizeMptrsArea, SizeOldBagsArea, SizeYoungBagsArea, SizeAllocationArea, SizeUnavailableArea)


cpdef error_enter_libgap_block_twice():
    """
    Demonstrate that we catch errors from entering a block twice.

    EXAMPLES::

        sage: from sage.libs.gap.util import error_enter_libgap_block_twice
        sage: error_enter_libgap_block_twice()
        Traceback (most recent call last):
        ...
        RuntimeError: Entered a critical block twice
    """
    from sage.libs.gap.libgap import libgap
    try:
        # The exception will be seen by this sig_on() after being
        # raised by the second libgap_enter().
        sig_on()
        libgap_enter()
        libgap_enter()
        sig_off()
    finally:
        libgap_exit()


cpdef error_exit_libgap_block_without_enter():
    """
    Demonstrate that we catch errors from omitting libgap_enter.

    EXAMPLES::

        sage: from sage.libs.gap.util import error_exit_libgap_block_without_enter
        sage: error_exit_libgap_block_without_enter()
        Traceback (most recent call last):
        ...
        RuntimeError: Called libgap_exit without previous libgap_enter
    """
    from sage.libs.gap.libgap import libgap
    sig_on()
    libgap_exit()
    sig_off()

############################################################################
### Auxilliary functions ###################################################
############################################################################


def command(command_string):
    """
    Playground for accessing Gap via libGap.

    You should not use this function in your own programs. This is
    just here for convenience if you want to play with the libgap
    libray code.

    EXAMPLES::

        sage: from sage.libs.gap.util import command
        sage: command('1')
        Output follows...
        1

        sage: command('1/0')
        Traceback (most recent call last):
        ...
        ValueError: libGAP: Error, Rational operations: <divisor> must not be zero

        sage: command('NormalSubgroups')
        Output follows...
        <Attribute "NormalSubgroups">

        sage: command('rec(a:=1, b:=2)')
        Output follows...
        rec( a := 1, b := 2 )
    """
    initialize()
    cdef libGAP_ExecStatus status

    cmd = command_string + ';\n'
    try:
        libgap_enter()
        libgap_start_interaction(cmd)
        try:
            sig_on()
            status = libGAP_ReadEvalCommand(libGAP_BottomLVars)
            if status != libGAP_STATUS_END:
                libgap_call_error_handler()
            sig_off()
        except RuntimeError, msg:
            raise ValueError('libGAP: '+str(msg).strip())

        assert libGAP_Symbol == libGAP_S_SEMICOLON, 'Did not end with semicolon?'
        libGAP_GetSymbol()
        if libGAP_Symbol != libGAP_S_EOF:
            raise ValueError('command() expects a single statement.')

        if libGAP_ReadEvalResult:
            libGAP_ViewObjHandler(libGAP_ReadEvalResult)
            s = libgap_get_output()
            print 'Output follows...'
            print s.strip()
        else:
            print 'No output.'

    finally:
        libgap_exit()
        libgap_finish_interaction()

    DEBUG_CHECK(libGAP_ReadEvalResult)

    if libGAP_ReadEvalResult != NULL:
        libgap_enter()
        libGAP_AssGVar(libGAP_Last3, libGAP_VAL_GVAR(libGAP_Last2))
        libGAP_AssGVar(libGAP_Last2, libGAP_VAL_GVAR(libGAP_Last))
        libGAP_AssGVar(libGAP_Last, libGAP_ReadEvalResult)
        libgap_exit()

