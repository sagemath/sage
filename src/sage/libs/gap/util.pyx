"""
Utility functions for GAP
"""

#*****************************************************************************
#       Copyright (C) 2012 Volker Braun <vbraun.name@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from libc.signal cimport signal, SIGCHLD, SIG_DFL
from posix.dlfcn cimport dlopen, dlclose, RTLD_NOW, RTLD_GLOBAL

from cpython.exc cimport PyErr_Fetch, PyErr_Restore
from cpython.object cimport Py_LT, Py_LE, Py_EQ, Py_NE, Py_GT, Py_GE
from cpython.ref cimport PyObject, Py_XINCREF, Py_XDECREF
from cysignals.signals cimport sig_on, sig_off

import os
import warnings
import sage.env

from .gap_includes cimport *
from .element cimport *
from sage.cpython.string import FS_ENCODING
from sage.cpython.string cimport str_to_bytes, char_to_str
from sage.interfaces.gap_workspace import prepare_workspace_dir


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
        Comparison wrapped Obj.

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
        cdef Obj self_value = self.value
        cdef Obj other_value = other.value
        if op == Py_LT:
            return self_value < other_value
        elif op == Py_LE:
            return self_value <= other_value
        elif op == Py_EQ:
            return self_value == other_value
        elif op == Py_GT:
            return self_value > other_value
        elif op == Py_GE:
            return self_value >= other_value
        elif op == Py_NE:
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
        return <Py_hash_t>(self.value)


cdef ObjWrapper wrap_obj(Obj obj):
    """
    Constructor function for :class:`ObjWrapper`
    """
    cdef ObjWrapper result = ObjWrapper.__new__(ObjWrapper)
    result.value = obj
    return result


# a dictionary to keep all GAP elements
# needed for GASMAN callbacks
#
cdef dict owned_objects_refcount = dict()

#
# used in Sage's libgap.Gap.count_GAP_objects
#
cpdef get_owned_objects():
    """
    Helper to access the refcount dictionary from Python code
    """
    return owned_objects_refcount


cdef void reference_obj(Obj obj):
    """
    Reference ``obj``
    """
    cdef ObjWrapper wrapped = wrap_obj(obj)
    global owned_objects_refcount
#    print("reference_obj called "+ crepr(obj) +"\n")
    if wrapped in owned_objects_refcount:
        owned_objects_refcount[wrapped] += 1
    else:
        owned_objects_refcount[wrapped] = 1


cdef void dereference_obj(Obj obj):
    """
    Reference ``obj``
    """
    cdef ObjWrapper wrapped = wrap_obj(obj)
    global owned_objects_refcount
    refcount = owned_objects_refcount.pop(wrapped)
    if refcount > 1:
        owned_objects_refcount[wrapped] = refcount - 1


cdef void gasman_callback() with gil:
    """
    Callback before each GAP garbage collection
    """
    global owned_objects_refcount
    for obj in owned_objects_refcount:
        MarkBag((<ObjWrapper>obj).value)


############################################################################
### Initialization of GAP ##################################################
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
    if os.path.exists(sage.env.GAP_ROOT_DIR):
        return sage.env.GAP_ROOT_DIR

    # Attempt to figure out the appropriate GAP_ROOT by reading the
    # local/bin/gap shell script; this is an ugly hack that exists for
    # historical reasons; the best approach to setting where Sage looks for
    # the appropriate GAP_ROOT is to set the GAP_ROOT_DIR variable
    SAGE_LOCAL = sage.env.SAGE_LOCAL
    with open(os.path.join(SAGE_LOCAL, 'bin', 'gap')) as f:
        gap_sh = f.read().splitlines()
    gapdir = next(x for x in gap_sh if x.strip().startswith('GAP_ROOT'))
    gapdir = gapdir.split('"')[1]
    gapdir = gapdir.replace('$SAGE_LOCAL', SAGE_LOCAL)
    return gapdir


# To ensure that we call initialize_libgap only once.
cdef bint _gap_is_initialized = False


cdef char* _reset_error_output_cmd = """\
libgap_errout := "";
MakeReadWriteGlobal("ERROR_OUTPUT");
ERROR_OUTPUT := OutputTextString(libgap_errout, false);
MakeReadOnlyGlobal("ERROR_OUTPUT");
"""

cdef char* _close_error_output_cmd = """\
CloseStream(ERROR_OUTPUT);
MakeReadWriteGlobal("ERROR_OUTPUT");
ERROR_OUTPUT := "*errout*";
MakeReadOnlyGlobal("ERROR_OUTPUT");
MakeImmutable(libgap_errout);
"""


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
    # Hack to ensure that all symbols provided by libgap are loaded into the
    # global symbol table
    # Note: we could use RTLD_NOLOAD and avoid the subsequent dlclose() but
    # this isn't portable

    cdef void* handle
    libgapname = str_to_bytes(sage.env.GAP_SO)
    handle = dlopen(libgapname, RTLD_NOW | RTLD_GLOBAL)
    if handle is NULL:
        raise RuntimeError(
                "Could not dlopen() libgap even though it should already "
                "be loaded!")
    dlclose(handle)

    # Define argv variable, which we will pass in to
    # initialize GAP. Note that we must pass define the memory pool
    # size!
    cdef char* argv[18]
    argv[0] = "sage"
    argv[1] = "-l"
    s = str_to_bytes(gap_root(), FS_ENCODING, "surrogateescape")
    argv[2] = s

    from sage.interfaces.gap import _get_gap_memory_pool_size_MB
    memory_pool = str_to_bytes(_get_gap_memory_pool_size_MB())
    argv[3] = "-o"
    argv[4] = memory_pool
    argv[5] = "-s"
    argv[6] = memory_pool

    argv[7] = "-m"
    argv[8] = "64m"

    argv[9] = "-q"    # no prompt!
    argv[10] = "-E"   # don't use readline as this will interfere with Python
    argv[11] = "--nointeract"  # Implies -T
    argv[12] = "-x"    # set the "screen" width so that GAP is less likely to
    argv[13] = "4096"  # insert newlines when printing objects
                       # 4096 unfortunately is the hard-coded max, but should
                       # be long enough for most cases

    cdef int argc = 14   # argv[argc] must be NULL

    from .saved_workspace import workspace
    workspace, workspace_is_up_to_date = workspace()
    ws = str_to_bytes(workspace, FS_ENCODING, "surrogateescape")
    if workspace_is_up_to_date:
        argv[argc] = "-L"
        argv[argc + 1] = ws
        argc += 2

    # Get the path to the sage.gaprc file and check that it exists
    sage_gaprc = os.path.join(os.path.dirname(__file__), 'sage.gaprc')
    if not os.path.exists(sage_gaprc):
        warnings.warn(f"Sage's GAP initialization file {sage_gaprc} is "
                       "is missing; some functionality may be limited")
    else:
        sage_gaprc = str_to_bytes(sage_gaprc, FS_ENCODING, "surrogateescape")
        argv[argc] = sage_gaprc
        argc += 1

    argv[argc] = NULL

    sig_on()
    # Initialize GAP but disable their SIGINT handler
    GAP_Initialize(argc, argv, gasman_callback, error_handler,
                   handleSignals=False)
    sig_off()

    # Disable GAP's SIGCHLD handler ChildStatusChanged(), which calls
    # waitpid() on random child processes.
    signal(SIGCHLD, SIG_DFL)

    # Set the ERROR_OUTPUT global in GAP to an output stream in which to
    # receive error output
    GAP_EvalString(_reset_error_output_cmd)

    # Prepare global GAP variable to hold temporary GAP objects
    global reference_holder
    reference_holder = GVarName("$SAGE_libgap_reference_holder")

    # Finished!
    _gap_is_initialized = True

    # Save a new workspace if necessary
    if not workspace_is_up_to_date:
        prepare_workspace_dir()
        from sage.misc.temporary_file import atomic_write
        with atomic_write(workspace) as f:
            f.close()
            gap_eval('SaveWorkspace("{0}")'.format(f.name))


############################################################################
### Evaluate string in GAP #################################################
############################################################################

cdef Obj gap_eval(str gap_string) except? NULL:
    r"""
    Evaluate a string in GAP.

    INPUT:

    - ``gap_string`` -- string. A valid statement in GAP.

    OUTPUT:

    The resulting GAP object or NULL+Python Exception in case of error.
    The result object may also be NULL without a Python exception set for
    statements that do not return a value.

    EXAMPLES::

        sage: libgap.eval('if 4>3 then\nPrint("hi");\nfi')
        sage: libgap.eval('1+1')   # testing that we have successfully recovered
        2

        sage: libgap.eval('if 4>3 thenPrint("hi");\nfi')
        Traceback (most recent call last):
        ...
        GAPError: Syntax error: then expected in stream:1
        if 4>3 thenPrint("hi");
               ^^^^^^^^^
        sage: libgap.eval('1+1')   # testing that we have successfully recovered
        2

    TESTS:

    A bad eval string that results in multiple statement evaluations by GAP
    and hence multiple errors should still result in a single exception
    with a message capturing all errors that occurrer::

        sage: libgap.eval('Complex Field with 53 bits of precision;')
        Traceback (most recent call last):
        ...
        GAPError: Error, Variable: 'Complex' must have a value
        Syntax error: ; expected in stream:1
        Complex Field with 53 bits of precision;;
         ^^^^^^^^^^^^
        Error, Variable: 'with' must have a value
        Syntax error: ; expected in stream:1
        Complex Field with 53 bits of precision;;
         ^^^^^^^^^^^^^^^^^^^^
        Error, Variable: 'bits' must have a value
        Syntax error: ; expected in stream:1
        Complex Field with 53 bits of precision;;
         ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        Error, Variable: 'precision' must have a value

    Test that on a subsequent attempt we get the same message (no garbage was
    left in the error stream)::

        sage: libgap.eval('Complex Field with 53 bits of precision;')
        Traceback (most recent call last):
        ...
        GAPError: Error, Variable: 'Complex' must have a value
        ...
        Error, Variable: 'precision' must have a value

        sage: libgap.eval('1+1')  # test that we successfully recover
        2
    """
    initialize()
    cdef Obj result
    cdef int i, j, nresults

    # Careful: We need to keep a reference to the bytes object here
    # so that Cython doesn't deallocate it before GAP is done with
    # its contents.
    cmd = str_to_bytes(gap_string + ';\n')
    sig_on()
    try:
        GAP_Enter()
        result = GAP_EvalString(cmd)
        # We can assume that the result object is a GAP PList (plain list)
        # and we should use functions for PLists directly for now; see
        # https://github.com/gap-system/gap/pull/2988/files#r233021437

        # If an error occurred in GAP_EvalString we won't even get
        # here if the error handler was set; but in case it wasn't
        # let's still check the result...
        nresults = LEN_LIST(result)
        if nresults > 1:  # to mimick the old libGAP
            # TODO: Get rid of this restriction eventually?
            raise GAPError("can only evaluate a single statement")

        # Get the result of the first statement
        result = ELM0_LIST(result, 1) # 1-indexed!

        if ELM0_LIST(result, 1) != GAP_True:
            # An otherwise unhandled error occurred in GAP (such as a
            # syntax error).  Try running the error handler manually
            # to capture the error output, if any.
            # This should result in a GAPError being set.
            error_handler_check_exception()

        # The actual resultant object, if any, is in the second entry
        # (which may be unassigned--see previous github comment; in this case
        # 0 is returned without setting a Python exception, so we should treat
        # this like returning None)

        return ELM0_LIST(result, 2)
    finally:
        GAP_Leave()
        sig_off()


###########################################################################
### Helper to protect temporary objects from deletion ######################
############################################################################

# Hold a reference (inside the GAP kernel) to obj so that it doesn't
# get deleted this works by assigning it to a global variable. This is
# very simple, but you can't use it to keep two objects alive. Be
# careful.
cdef UInt reference_holder

cdef void hold_reference(Obj obj):
    """
    Hold a reference (inside the GAP kernel) to obj

    This ensures that the GAP garbage collector does not delete
    ``obj``. This works by assigning it to a global variable. This is
    very simple, but you can't use it to keep two objects alive. Be
    careful.
    """
    global reference_holder
    AssGVar(reference_holder, obj)


############################################################################
### Error handler ##########################################################
############################################################################

class GAPError(ValueError):  # ValueError for historical reasons
    """
    Exceptions raised by the GAP library
    """


cdef str extract_libgap_errout():
    """
    Reads the global variable libgap_errout and returns a Python string
    containing the error message (with some boilerplate removed).
    """
    cdef Obj r
    cdef char *msg

    r = GAP_ValueGlobalVariable("libgap_errout")

    # Grab a pointer to the C string underlying the GAP string libgap_errout
    # then copy it to a Python str (char_to_str contains an implicit strcpy)
    msg = CSTR_STRING(r)
    if msg != NULL:
        msg_py = char_to_str(msg)
        msg_py = msg_py.replace('For debugging hints type ?Recovery from '
                                'NoMethodFound\n', '').strip()
    else:
        # Shouldn't happen but just in case...
        msg_py = ""

    return msg_py


cdef void error_handler() with gil:
    """
    The libgap error handler.

    If an error occurred, we raise a ``GAPError``; when the original
    ``GAP_EvalString`` returns, this exception will be seen.

    TODO: We should probably prevent re-entering this function if we
    are already handling an error; if there is an error in our stream
    handling code below it could result in a stack overflow.
    """
    cdef PyObject* exc_type = NULL
    cdef PyObject* exc_val = NULL
    cdef PyObject* exc_tb = NULL

    try:
        GAP_EnterStack()

        # Close the error stream: this flushes any remaining output and
        # closes the stream for further writing; reset ERROR_OUTPUT to
        # something sane just in case (trying to print to a closed
        # stream segfaults GAP)
        GAP_EvalStringNoExcept(_close_error_output_cmd)

        # Fetch any existing exception before calling
        # extract_libgap_errout() so that the exception indicator is
        # cleared
        PyErr_Fetch(&exc_type, &exc_val, &exc_tb)

        msg = extract_libgap_errout()
        # Sometimes error_handler() can be called multiple times
        # from a single GAP_EvalString call before it returns.
        # In this case, we just update the exception by appending
        # to the existing exception message
        if exc_type is <PyObject*>GAPError and exc_val is not NULL:
            msg = str(<object>exc_val) + '\n' + msg
        elif not msg:
            msg = "an unknown error occurred in GAP"

        # Raise an exception using PyErr_Restore().
        # This way, we can keep any existing traceback object.
        # Note that we manually need to deal with refcounts here.
        Py_XDECREF(exc_type)
        Py_XDECREF(exc_val)
        exc_type = <PyObject*>GAPError
        exc_val = <PyObject*>msg
        Py_XINCREF(exc_type)
        Py_XINCREF(exc_val)
        PyErr_Restore(exc_type, exc_val, exc_tb)
    finally:
        # Reset ERROR_OUTPUT with a new text string stream
        GAP_EvalStringNoExcept(_reset_error_output_cmd)
        GAP_LeaveStack()


cdef void error_handler_check_exception() except *:
    error_handler()
