r"""
Library interface to Embeddable Common Lisp (ECL)
"""
#*****************************************************************************
#       Copyright (C) 2009 Nils Bruin <nbruin@sfu.ca>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function

#This version of the library interface prefers to convert ECL integers and
#rationals to SAGE types Integer and Rational. These parts could easily be
#adapted to work with pure Python types.

include "cysignals/signals.pxi"
include "sage/ext/cdefs.pxi"

from libc.stdlib cimport abort
from libc.signal cimport SIGINT, SIGBUS, SIGSEGV, SIGCHLD
from libc.signal cimport raise_ as signal_raise
from posix.signal cimport sigaction, sigaction_t

from sage.rings.integer cimport Integer
from sage.rings.rational cimport Rational

#it would be preferrable to let bint_symbolp wrap an efficient macro
#but the macro provided in object.h doesn't seem to work
cdef bint bint_symbolp(cl_object obj):
    return not(cl_symbolp(obj) == Cnil)

#these type predicates are only provided in "cl_*" form, so we wrap them
#with the proper type cast.

cdef bint bint_numberp(cl_object obj):
    return not(cl_numberp(obj) == Cnil)
cdef bint bint_integerp(cl_object obj):
    return not(cl_integerp(obj) == Cnil)
cdef bint bint_rationalp(cl_object obj):
    return not(cl_rationalp(obj) == Cnil)

cdef extern from "eclsig.h":
    int ecl_sig_on() except 0
    void ecl_sig_off()
    cdef sigaction_t ecl_sigint_handler
    cdef sigaction_t ecl_sigbus_handler
    cdef sigaction_t ecl_sigsegv_handler
    cdef mpz_t ecl_mpz_from_bignum(cl_object obj)
    cdef cl_object ecl_bignum_from_mpz(mpz_t num)

cdef cl_object string_to_object(char * s):
    return ecl_read_from_cstring(s)

# We need to keep a list of objects bound to python, to protect them from being
# garbage collected. We want a list in which we can quickly add and remove
# elements. Lookup is not necessary. A doubly linked list seems
# most appropriate. A node looks like
#   N = ( value next . prev)
# so that car(N)=value, cadr(N)=next, cddr(N)=prev.
# we write routines to insert a node after a given node
# and to delete a given node. This can all be done with modifying pointers.
# note that circular structures are unpleasant for most lisp routines.
# perhaps this even puts a strain on the garbage collector?
# an alternative data structure would be an array where the free nodes get
# chained in a "free list" for quick allocation (and if the free list is empty
# upon allocating a node, the array needs to be extended)

cdef cl_object insert_node_after(cl_object node,cl_object value):
    cdef cl_object next,newnode

    next=cl_cadr(node)
    newnode=cl_cons(value,cl_cons(next,node))
    cl_rplaca(cl_cdr(node),newnode)
    if next != Cnil:
        cl_rplacd(cl_cdr(next),newnode)
    return newnode

cdef void remove_node(cl_object node):
    cdef cl_object next, prev
    next=cl_cadr(node)
    prev=cl_cddr(node)
    if next != Cnil:
        cl_rplacd(cl_cdr(next),prev)
    if prev != Cnil:
        cl_rplaca(cl_cdr(prev),next)

# our global list of pointers. This will be a pointer to a sentinel node,
# after which all new nodes can be inserted. list_of_object gets initialised
# by init_ecl() and bound to the global ECL variable *SAGE-LIST-OF-OBJECTS*

cdef cl_object list_of_objects

cdef cl_object safe_eval_clobj         #our own error catching eval
cdef cl_object safe_apply_clobj        #our own error catching apply
cdef cl_object safe_funcall_clobj      #our own error catching funcall
cdef cl_object read_from_string_clobj  #our own error catching reader

cdef bint ecl_has_booted = 0

# ECL signal handling

def test_sigint_before_ecl_sig_on():
    """
    TESTS:

    If an interrupt arrives *before* ecl_sig_on(), we should get an
    ordinary KeyboardInterrupt::

        sage: from sage.libs.ecl import test_sigint_before_ecl_sig_on
        sage: test_sigint_before_ecl_sig_on()
        Traceback (most recent call last):
        ...
        KeyboardInterrupt
    """
    # Raise a SIGINT *now*.  Since we are outside of sig_on() at this
    # point, this SIGINT will not be seen yet.
    signal_raise(SIGINT)
    # An ordinary KeyboardInterrupt should be raised by ecl_sig_on()
    # since ecl_sig_on() calls sig_on() before anything else.  This
    # will catch the pending SIGINT.
    ecl_sig_on()
    # We should never get here.
    abort()

def test_ecl_options():
    """
    Print an overview of the ECL options

    TESTS::

        sage: from sage.libs.ecl import test_ecl_options
        sage: test_ecl_options()
        ECL_OPT_INCREMENTAL_GC = 0
        ECL_OPT_TRAP_SIGSEGV = 1
        ECL_OPT_TRAP_SIGFPE = 1
        ECL_OPT_TRAP_SIGINT = 1
        ECL_OPT_TRAP_SIGILL = 1
        ECL_OPT_TRAP_SIGBUS = 1
        ECL_OPT_TRAP_SIGCHLD = 0
        ECL_OPT_TRAP_SIGPIPE = 1
        ECL_OPT_TRAP_INTERRUPT_SIGNAL = 1
        ECL_OPT_SIGNAL_HANDLING_THREAD = 0
        ECL_OPT_SIGNAL_QUEUE_SIZE = 16
        ECL_OPT_BOOTED = 1
        ECL_OPT_BIND_STACK_SIZE = ...
        ECL_OPT_BIND_STACK_SAFETY_AREA = ...
        ECL_OPT_FRAME_STACK_SIZE = ...
        ECL_OPT_FRAME_STACK_SAFETY_AREA = ...
        ECL_OPT_LISP_STACK_SIZE = ...
        ECL_OPT_LISP_STACK_SAFETY_AREA = ...
        ECL_OPT_C_STACK_SIZE = ...
        ECL_OPT_C_STACK_SAFETY_AREA = ...
        ECL_OPT_SIGALTSTACK_SIZE = 1
        ECL_OPT_HEAP_SIZE = ...
        ECL_OPT_HEAP_SAFETY_AREA = ...
        ECL_OPT_THREAD_INTERRUPT_SIGNAL = 0
        ECL_OPT_SET_GMP_MEMORY_FUNCTIONS = 0
    """
    print('ECL_OPT_INCREMENTAL_GC = {0}'.format(
        ecl_get_option(ECL_OPT_INCREMENTAL_GC)))
    print('ECL_OPT_TRAP_SIGSEGV = {0}'.format(
        ecl_get_option(ECL_OPT_TRAP_SIGSEGV)))
    print('ECL_OPT_TRAP_SIGFPE = {0}'.format(
        ecl_get_option(ECL_OPT_TRAP_SIGFPE)))
    print('ECL_OPT_TRAP_SIGINT = {0}'.format(
        ecl_get_option(ECL_OPT_TRAP_SIGINT)))
    print('ECL_OPT_TRAP_SIGILL = {0}'.format(
        ecl_get_option(ECL_OPT_TRAP_SIGILL)))
    print('ECL_OPT_TRAP_SIGBUS = {0}'.format(
        ecl_get_option(ECL_OPT_TRAP_SIGBUS)))
    print('ECL_OPT_TRAP_SIGCHLD = {0}'.format(
        ecl_get_option(ECL_OPT_TRAP_SIGCHLD)))
    print('ECL_OPT_TRAP_SIGPIPE = {0}'.format(
        ecl_get_option(ECL_OPT_TRAP_SIGPIPE)))
    print('ECL_OPT_TRAP_INTERRUPT_SIGNAL = {0}'.format(
        ecl_get_option(ECL_OPT_TRAP_INTERRUPT_SIGNAL)))
    print('ECL_OPT_SIGNAL_HANDLING_THREAD = {0}'.format(
        ecl_get_option(ECL_OPT_SIGNAL_HANDLING_THREAD)))
    print('ECL_OPT_SIGNAL_QUEUE_SIZE = {0}'.format(
        ecl_get_option(ECL_OPT_SIGNAL_QUEUE_SIZE)))
    print('ECL_OPT_BOOTED = {0}'.format(
        ecl_get_option(ECL_OPT_BOOTED)))
    print('ECL_OPT_BIND_STACK_SIZE = {0}'.format(
        ecl_get_option(ECL_OPT_BIND_STACK_SIZE)))
    print('ECL_OPT_BIND_STACK_SAFETY_AREA = {0}'.format(
        ecl_get_option(ECL_OPT_BIND_STACK_SAFETY_AREA)))
    print('ECL_OPT_FRAME_STACK_SIZE = {0}'.format(
        ecl_get_option(ECL_OPT_FRAME_STACK_SIZE)))
    print('ECL_OPT_FRAME_STACK_SAFETY_AREA = {0}'.format(
        ecl_get_option(ECL_OPT_FRAME_STACK_SAFETY_AREA)))
    print('ECL_OPT_LISP_STACK_SIZE = {0}'.format(
        ecl_get_option(ECL_OPT_LISP_STACK_SIZE)))
    print('ECL_OPT_LISP_STACK_SAFETY_AREA = {0}'.format(
        ecl_get_option(ECL_OPT_LISP_STACK_SAFETY_AREA)))
    print('ECL_OPT_C_STACK_SIZE = {0}'.format(
        ecl_get_option(ECL_OPT_C_STACK_SIZE)))
    print('ECL_OPT_C_STACK_SAFETY_AREA = {0}'.format(
        ecl_get_option(ECL_OPT_C_STACK_SAFETY_AREA)))
    print('ECL_OPT_SIGALTSTACK_SIZE = {0}'.format(
        ecl_get_option(ECL_OPT_SIGALTSTACK_SIZE)))
    print('ECL_OPT_HEAP_SIZE = {0}'.format(
        ecl_get_option(ECL_OPT_HEAP_SIZE)))
    print('ECL_OPT_HEAP_SAFETY_AREA = {0}'.format(
        ecl_get_option(ECL_OPT_HEAP_SAFETY_AREA)))
    print('ECL_OPT_THREAD_INTERRUPT_SIGNAL = {0}'.format(
        ecl_get_option(ECL_OPT_THREAD_INTERRUPT_SIGNAL)))
    print('ECL_OPT_SET_GMP_MEMORY_FUNCTIONS = {0}'.format(
        ecl_get_option(ECL_OPT_SET_GMP_MEMORY_FUNCTIONS)))

def init_ecl():
    r"""
    Internal function to initialize ecl. Do not call.

    This function initializes the ECL library for use within Python.
    This routine should only be called once and importing the ecl library
    interface already does that, so do not call this yourself.

    EXAMPLES::

        sage: from sage.libs.ecl import *

    At this point, init_ecl() has run. Explicitly executing it
    gives an error::

        sage: init_ecl()
        Traceback (most recent call last):
        ...
        RuntimeError: ECL is already initialized
    """
    global list_of_objects
    global safe_eval_clobj
    global safe_apply_clobj
    global safe_funcall_clobj
    global read_from_string_clobj
    global ecl_has_booted
    cdef char *argv[1]
    cdef sigaction_t sage_action[32]
    cdef int i

    if ecl_has_booted:
        raise RuntimeError("ECL is already initialized")

    # we need it to stop handling SIGCHLD
    ecl_set_option(ECL_OPT_TRAP_SIGCHLD, 0);

    #we keep our own GMP memory functions. ECL should not claim them
    ecl_set_option(ECL_OPT_SET_GMP_MEMORY_FUNCTIONS,0);

    #we need a dummy argv for cl_boot (we just don't give any parameters)
    argv[0]="sage"

    #get all the signal handlers before initializing Sage so we can
    #put them back afterwards.
    for i in range(1,32):
        sigaction(i, NULL, &sage_action[i])

    #initialize ECL
    ecl_set_option(ECL_OPT_SIGNAL_HANDLING_THREAD, 0)
    cl_boot(1, argv)

    #save signal handler from ECL
    sigaction(SIGINT, NULL, &ecl_sigint_handler)
    sigaction(SIGBUS, NULL, &ecl_sigbus_handler)
    sigaction(SIGSEGV, NULL, &ecl_sigsegv_handler)

    #verify that no SIGCHLD handler was installed
    cdef sigaction_t sig_test
    sigaction(SIGCHLD, NULL, &sig_test)
    assert sage_action[SIGCHLD].sa_handler == NULL  # Sage does not set SIGCHLD handler
    assert sig_test.sa_handler == NULL              # And ECL bootup did not set one 

    #and put the Sage signal handlers back
    for i in range(1,32):
        sigaction(i, &sage_action[i], NULL)

    #initialise list of objects and bind to global variable
    # *SAGE-LIST-OF-OBJECTS* to make it rooted in the reachable tree for the GC
    list_of_objects=cl_cons(Cnil,cl_cons(Cnil,Cnil))
    cl_set(string_to_object("*SAGE-LIST-OF-OBJECTS*"),list_of_objects)

    # We define our own error catching eval, apply and funcall/
    # Presently these routines are only converted to byte-code. If they
    # ever turn out to be a bottle neck, it should be easy to properly
    # compile them.

    read_from_string_clobj=cl_eval(string_to_object("(symbol-function 'read-from-string)"))

    cl_eval(string_to_object("""
        (defun sage-safe-eval (form)
            (handler-case
                (values (eval form))
                (serious-condition (cnd)
                    (values nil (princ-to-string cnd)))))
        """))
    safe_eval_clobj=cl_eval(string_to_object("(symbol-function 'sage-safe-eval)"))

    cl_eval(string_to_object("""
        (defun sage-safe-apply (func args)
            (handler-case
                (values (apply func args))
                (serious-condition (cnd)
                    (values nil (princ-to-string cnd)))))
        """))

    safe_apply_clobj=cl_eval(string_to_object("(symbol-function 'sage-safe-apply)"))
    cl_eval(string_to_object("""
        (defun sage-safe-funcall (func arg)
            (handler-case
                (values (funcall func arg))
                (serious-condition (cnd)
                    (values nil (princ-to-string cnd)))))
        """))
    safe_funcall_clobj=cl_eval(string_to_object("(symbol-function 'sage-safe-funcall)"))

    ecl_has_booted = 1

cdef cl_object ecl_safe_eval(cl_object form) except NULL:
    """
    TESTS:

    Test interrupts::

        sage: from sage.libs.ecl import *
        sage: from cysignals.tests import interrupt_after_delay
        sage: ecl_eval("(setf i 0)")
        <ECL: 0>
        sage: inf_loop = ecl_eval("(defun infinite() (loop (incf i)))")
        sage: interrupt_after_delay(1000)
        sage: inf_loop()
        Traceback (most recent call last):
        ...
        RuntimeError: ECL says: Console interrupt.
    """
    cdef cl_object s
    ecl_sig_on()
    cl_funcall(2,safe_eval_clobj,form)
    ecl_sig_off()

    if ecl_nvalues > 1:
        s = si_coerce_to_base_string(ecl_values(1))
        raise RuntimeError("ECL says: "+ecl_base_string_pointer_safe(s))
    else:
        return ecl_values(0)

cdef cl_object ecl_safe_funcall(cl_object func, cl_object arg) except NULL:
    cdef cl_object l, s
    l = cl_cons(func,cl_cons(arg,Cnil));

    ecl_sig_on()
    cl_apply(2,safe_funcall_clobj,cl_cons(func,cl_cons(arg,Cnil)))
    ecl_sig_off()

    if ecl_nvalues > 1:
        s = si_coerce_to_base_string(ecl_values(1))
        raise RuntimeError("ECL says: "+ecl_base_string_pointer_safe(s))
    else:
        return ecl_values(0)

cdef cl_object ecl_safe_apply(cl_object func, cl_object args) except NULL:
    cdef cl_object s
    ecl_sig_on()
    cl_funcall(3,safe_apply_clobj,func,args)
    ecl_sig_off()

    if ecl_nvalues > 1:
        s = si_coerce_to_base_string(ecl_values(1))
        raise RuntimeError("ECL says: "+ecl_base_string_pointer_safe(s))
    else:
        return ecl_values(0)

cdef cl_object ecl_safe_read_string(char * s) except NULL:
    cdef cl_object o
    o = ecl_cstring_to_base_string_or_nil(s)
    o = ecl_safe_funcall(read_from_string_clobj,o)
    return o

def shutdown_ecl():
    r"""
    Shut down ecl. Do not call.

    Given the way that ECL is used from python, it is very difficult to ensure
    that no ECL objects exist at a particular time. Hence, destroying ECL is a
    risky proposition.

    EXAMPLE::

        sage: from sage.libs.ecl import *
        sage: shutdown_ecl()
    """
    cl_shutdown()

#this prints the objects that sage wants the GC to keep track of.
#these should be all non-immediate EclObject wrapped objects
def print_objects():
    r"""
    Print GC-protection list

    Diagnostic function. ECL objects that are bound to Python objects need to
    be protected from being garbage collected. We do this by including them
    in a doubly linked list bound to the global ECL symbol
    *SAGE-LIST-OF-OBJECTS*. Only non-immediate values get included, so
    small integers do not get linked in. This routine prints the values
    currently stored.

    EXAMPLE::

        sage: from sage.libs.ecl import *
        sage: a=EclObject("hello")
        sage: b=EclObject(10)
        sage: c=EclObject("world")
        sage: print_objects() #random because previous test runs can have left objects
        NIL
        WORLD
        HELLO
    """

    cdef cl_object c, s
    c = list_of_objects
    while True:
        s = si_coerce_to_base_string(cl_write_to_string(1,cl_car(c)))
        print(ecl_base_string_pointer_safe(s))
        c = cl_cadr(c)
        if c == Cnil:
            break

cdef cl_object python_to_ecl(pyobj) except NULL:
    # conversion of a python object into an ecl object
    # most conversions are straightforward. Noteworthy are:
    # python lists -> lisp (NIL terminated) lists
    # tuples -> dotted lists
    # strings ->parsed by lisp reader

    cdef bytes s
    cdef cl_object L, ptr

    if isinstance(pyobj,bool):
        if pyobj:
            return Ct
        else:
            return Cnil
    elif pyobj is None:
        return Cnil
    elif isinstance(pyobj,int):
        return ecl_make_integer(pyobj)
    elif isinstance(pyobj,long):
        if pyobj >= MOST_NEGATIVE_FIXNUM and pyobj <= MOST_POSITIVE_FIXNUM:
            return ecl_make_integer(pyobj)
        else:
            return python_to_ecl(Integer(pyobj))
    elif isinstance(pyobj,float):
        return ecl_make_doublefloat(pyobj)
    elif isinstance(pyobj,unicode):
        s=<bytes>(str(pyobj))
        return ecl_safe_read_string(s)
    elif isinstance(pyobj,bytes):
        s=<bytes>pyobj
        return ecl_safe_read_string(s)
    elif isinstance(pyobj,Integer):
        if pyobj >= MOST_NEGATIVE_FIXNUM and pyobj <= MOST_POSITIVE_FIXNUM:
            return ecl_make_integer(pyobj)
        else:
            return ecl_bignum_from_mpz( (<Integer>pyobj).value )
    elif isinstance(pyobj,Rational):
        return ecl_make_ratio(
                python_to_ecl( (<Rational>pyobj).numerator()  ),
                python_to_ecl( (<Rational>pyobj).denominator()))
    elif isinstance(pyobj,EclObject):
        return (<EclObject>pyobj).obj
    elif isinstance(pyobj,list):
        if len(pyobj) == 0:
            return Cnil
        else:
            L=cl_cons(python_to_ecl(pyobj[0]),Cnil)
            ptr=L
            for a in pyobj[1:]:
                cl_rplacd(ptr,cl_cons(python_to_ecl(a),Cnil))
                ptr=cl_cdr(ptr)
            return L
    elif isinstance(pyobj,tuple):
        if len(pyobj) == 0:
            return Cnil
        elif len(pyobj) == 1:
            return python_to_ecl(pyobj[0])
        else:
            L=cl_cons(python_to_ecl(pyobj[0]),Cnil)
            ptr=L
            for a in pyobj[1:-1]:
                cl_rplacd(ptr,cl_cons(python_to_ecl(a),Cnil))
                ptr=cl_cdr(ptr)
            cl_rplacd(ptr,python_to_ecl(pyobj[-1]))
            return L
    else:
        raise TypeError("Unimplemented type for python_to_ecl")

cdef ecl_to_python(cl_object o):
    cdef cl_object s
    cdef Integer N
    # conversions from an ecl object to a python object.

    if o == Cnil:
        return None
    elif bint_fixnump(o):
        #SAGE specific conversion
        #return ecl_fixint(o)
        return Integer(ecl_fixint(o))
    elif bint_integerp(o):
        #SAGE specific conversion
        N = Integer.__new__(Integer)
        N.set_from_mpz(ecl_mpz_from_bignum(o))
        return N
    elif bint_rationalp(o):
        #SAGE specific conversion
        #vanilla python does not have a class to represent rational numbers
        return Rational((ecl_to_python(cl_numerator(o)),ecl_to_python(cl_denominator(o))))
    elif bint_floatp(o):
        #Python conversion
        #Since SAGE mainly uses mpfr, perhaps "double is not an appropriate return type
        return ecl_to_double(o)
    elif o == Ct:
        return True
    elif bint_consp(o):
        L=[]
        while o != Cnil:
            L.append(ecl_to_python(cl_car(o)))
            o = cl_cdr(o)
            if not(bint_listp(o)):
                L.append(ecl_to_python(o))
                return tuple(L)
        return L
    else:
        s = si_coerce_to_base_string(cl_write_to_string(1,o))
        return ecl_base_string_pointer_safe(s)

#Maxima's BFLOAT multiprecision float type can be read with:
#def bfloat_to_python(e):
#  prec=Integer(str(e.car().cddr().car()))
#  mant=Integer(str(e.cdr().car()))
#  exp=Integer(str(e.cddr().car()))
#  return 2^(exp-prec)*mant

cdef class EclObject:
    r"""
    Python wrapper of ECL objects

    The ``EclObject`` forms a wrapper around ECL objects. The wrapper ensures
    that the data structure pointed to is protected from garbage collection in
    ECL by installing a pointer to it from a global data structure within the
    scope of the ECL garbage collector. This pointer is destroyed upon
    destruction of the EclObject.

    EclObject() takes a Python object and tries to find a representation of it
    in Lisp.

    EXAMPLES:

    Python lists get mapped to LISP lists. None and Boolean values to
    appropriate values in LISP::

        sage: from sage.libs.ecl import *
        sage: EclObject([None,true,false])
        <ECL: (NIL T NIL)>

    Numerical values are translated to the appropriate type in LISP::

        sage: EclObject(1)
        <ECL: 1>
        sage: EclObject(10**40)
        <ECL: 10000000000000000000000000000000000000000>

    Floats in Python are IEEE double, which LISP has as well. However,
    the printing of floating point types in LISP depends on settings::

        sage: a = EclObject(float(10^40))
        sage: ecl_eval("(setf *read-default-float-format* 'single-float)")
        <ECL: SINGLE-FLOAT>
        sage: a
        <ECL: 1.d40>
        sage: ecl_eval("(setf *read-default-float-format* 'double-float)")
        <ECL: DOUBLE-FLOAT>
        sage: a
        <ECL: 1.e40>

    Tuples are translated to dotted lists::

        sage: EclObject( (false, true))
        <ECL: (NIL . T)>

    Strings are fed to the reader, so a string normally results in a symbol::

        sage: EclObject("Symbol")
        <ECL: SYMBOL>

    But with proper quotation one can construct a lisp string object too::

        sage: EclObject('"Symbol"')
        <ECL: "Symbol">

    EclObjects translate to themselves, so one can mix::

        sage: EclObject([1,2,EclObject([3])])
        <ECL: (1 2 (3))>

    Calling an EclObject translates into the appropriate LISP ``apply``,
    where the argument is transformed into an EclObject itself, so one can
    flexibly apply LISP functions::

        sage: car=EclObject("car")
        sage: cdr=EclObject("cdr")
        sage: car(cdr([1,2,3]))
        <ECL: 2>

    and even construct and evaluate arbitrary S-expressions::

        sage: eval=EclObject("eval")
        sage: quote=EclObject("quote")
        sage: eval([car, [cdr, [quote,[1,2,3]]]])
        <ECL: 2>

    TESTS:

    We check that multiprecision integers are converted correctly::

        sage: i = 10 ^ (10 ^ 5)
        sage: EclObject(i) == EclObject(str(i))
        True
        sage: EclObject(-i) == EclObject(str(-i))
        True
        sage: EclObject(i).python() == i
        True
        sage: EclObject(-i).python() == -i
        True
    """
    cdef cl_object obj   #the wrapped object
    cdef cl_object node  #linked list pointer: car(node) == obj

    cdef void set_obj(EclObject self, cl_object o):
        if self.node:
            remove_node(self.node)
            self.node=NULL
        self.obj=o
        if not(bint_fixnump(o) or bint_characterp(o) or bint_nullp(o)):
            self.node=insert_node_after(list_of_objects,o)

    def __init__(self,*args):
        r"""
        Create an EclObject

        See EclObject for full documentation.

        EXAMPLES::

            sage: from sage.libs.ecl import *
            sage: EclObject([None,true,false])
            <ECL: (NIL T NIL)>

        """
        if len(args) != 0:
            self.set_obj(python_to_ecl(args[0]))

    def __reduce__(self):
        r"""
        This is used for pickling. Not implemented

        Ecl does not natively support serialization of its objects, so the
        python wrapper class EclObject does not support pickling. There are
        independent efforts for developing serialization for Common Lisp, such as
        CL-STORE. Look at those if you need serialization of ECL objects.

        EXAMPLES::

            sage: from sage.libs.ecl import *
            sage: s=EclObject([1,2,3])
            sage: s.__reduce__()
            Traceback (most recent call last):
            ...
            NotImplementedError: EclObjects do not have a pickling method
            sage: s==loads(dumps(s))
            Traceback (most recent call last):
            ...
            NotImplementedError: EclObjects do not have a pickling method
        """
        raise NotImplementedError("EclObjects do not have a pickling method")

    def python(self):
        r"""
        Convert an EclObject to a python object.

        EXAMPLES::

            sage: from sage.libs.ecl import *
            sage: L=EclObject([1,2,("three",'"four"')])
            sage: L.python()
            [1, 2, ('THREE', '"four"')]

        """
        return ecl_to_python(self.obj)

    def __dealloc__(self):
        r"""
        Deallocate EclObject

        It is important to remove the GC preventing reference to the object upon
        deletion of the wrapper.

        EXAMPLES::

            sage: from sage.libs.ecl import *
            sage: L=EclObject("symbol")
            sage: del L

        """
        if self.node:
            remove_node(self.node)

    def __repr__(self):
        r"""
        Produce a string representation suitable for interactive printing.

        Converts the wrapped LISP object to a string, decorated in such a way that
        it can be recognised as a LISP object.

        EXAMPLES::

            sage: from sage.libs.ecl import *
            sage: L=EclObject("symbol")
            sage: repr(L)
            '<ECL: SYMBOL>'

        """
        return "<ECL: "+str(self)+">"

    def __str__(self):
        r"""
        Produce a string representation.

        Converts the wrapped LISP object to a string and returns that as a Python
        string.

        EXAMPLES::

            sage: from sage.libs.ecl import *
            sage: L=EclObject("symbol")
            sage: str(L)
            'SYMBOL'

        """
        cdef cl_object s
        s = si_coerce_to_base_string(cl_write_to_string(1,self.obj))
        return ecl_base_string_pointer_safe(s)

    def __hash__(self):
        r"""
        Return a hash value of the object

        Returns the hash value returned by SXHASH, which is a routine that is
        specified in Common Lisp. According to the specification, lisp objects that
        are EQUAL have the same SXHASH value. Since two EclObjects are equal if
        their wrapped objects are EQUAL according to lisp, this is compatible with
        Python's concept of hash values.

        It is not possible to enforce immutability of lisp objects, so care should
        be taken in using EclObjects as dictionary keys.

        EXAMPLES::

            sage: from sage.libs.ecl import *
            sage: L=EclObject([1,2])
            sage: L
            <ECL: (1 2)>
            sage: hash(L) #random
            463816586
            sage: L.rplacd(EclObject(3))
            sage: L
            <ECL: (1 . 3)>
            sage: hash(L) #random
            140404060

        """
        return ecl_fixint(cl_sxhash(self.obj))

    def __call__(self, *args):
        r"""
        Apply self to arguments.

        EXAMPLES::

            sage: from sage.libs.ecl import *
            sage: sqr=EclObject("(lambda (x) (* x x))").eval()
            sage: sqr(10)
            <ECL: 100>

        """
        lispargs = EclObject(list(args))
        return ecl_wrap(ecl_safe_apply(self.obj,(<EclObject>lispargs).obj))

    def __richcmp__(left, right, int op):
        r"""
        Comparison test.

        An EclObject is not equal to any non-EclObject. Two EclObjects are equal
        if their wrapped lisp objects are EQUAL. Since LISP has no univeral ordering,
        less than and greater than tests are not implemented for EclObjects.

        EXAMPLES::

            sage: from sage.libs.ecl import *
            sage: a=EclObject(1)
            sage: b=EclObject(2)
            sage: a==b
            False
            sage: a<b
            Traceback (most recent call last):
            ...
            NotImplementedError: EclObjects can only be compared for equality
            sage: EclObject("<")(a,b)
            <ECL: T>
        """
        if   op == 2: # "=="
            if not(isinstance(left,EclObject)) or not(isinstance(right,EclObject)):
                return False
            else:
                return bint_equal((<EclObject>left).obj,(<EclObject>right).obj)
        elif op == 3: # "!="
            if not(isinstance(left,EclObject)) or not(isinstance(right,EclObject)):
                return True
            else:
                return not(bint_equal((<EclObject>left).obj,(<EclObject>right).obj))

        #Common lisp only seems to be able to compare numeric and string types
        #and does not have generic routines for doing that.
        #we could dispatch based on type here, but that seems
        #inappropriate for an *interface*.
        raise NotImplementedError("EclObjects can only be compared for equality")

    def __iter__(self):
        r"""
        Implements the iterator protocol for EclObject.

        EclObject implements the iterator protocol for lists. This means
        one can use an EclObject in the context where an iterator is
        expected (for instance, in a list comprehension or in a for loop).
        The iterator produces EclObjects wrapping the members of the list that
        the original EclObject wraps.

        The wrappers returned are all newly constructed but refer to the
        original members of the list iterated over. This is usually what is
        intended but, just as in Python, can cause surprises if the original
        object is changed between calls to the iterator.

        Since EclObject translates Python Lists into LISP lists and Python
        tuples into LISP "dotted" lists (lists for which the final CDR is not
        necessarily NIL), and both these python structures are iterable, the
        corresponding EclObjects are iterable as well.

        EclObjects that are not lists are not iterable.

        EXAMPLES::

            sage: from sage.libs.ecl import *
            sage: [i for i in EclObject("(1 2 3)")]
            [<ECL: 1>, <ECL: 2>, <ECL: 3>]
            sage: [i for i in EclObject("(1 2 . 3)")]
            [<ECL: 1>, <ECL: 2>, <ECL: 3>]
            sage: [i for i in EclObject("NIL")]
            []

        TESTS:

        These show that Python lists and tuples behave as
        described above::

            sage: [i for i in EclObject([1,2,3])]
            [<ECL: 1>, <ECL: 2>, <ECL: 3>]
            sage: [i for i in EclObject((1,2,3))]
            [<ECL: 1>, <ECL: 2>, <ECL: 3>]

        This tests that we cannot iterate EclObjects we shouldn't,
        as described above::

            sage: [i for i in EclObject("T")]
            Traceback (most recent call last):
            ...
            TypeError: ECL object is not iterable

        """
        return EclListIterator(self)

    def eval(self):
        r"""
        Evaluate object as an S-Expression

        EXAMPLES::

            sage: from sage.libs.ecl import *
            sage: S=EclObject("(+ 1 2)")
            sage: S
            <ECL: (+ 1 2)>
            sage: S.eval()
            <ECL: 3>

        """
        cdef cl_object o
        o=ecl_safe_eval(self.obj)
        if o == NULL:
            raise RuntimeError("ECL runtime error")
        return ecl_wrap(o)

    def cons(self,EclObject d):
        r"""
        apply cons to self and argument and return the result.

        EXAMPLES::

            sage: from sage.libs.ecl import *
            sage: a=EclObject(1)
            sage: b=EclObject(2)
            sage: a.cons(b)
            <ECL: (1 . 2)>

        """
        return ecl_wrap(cl_cons(self.obj,d.obj))

    def rplaca(self,EclObject d):
        r"""
        Destructively replace car(self) with d.

        EXAMPLES::

            sage: from sage.libs.ecl import *
            sage: L=EclObject((1,2))
            sage: L
            <ECL: (1 . 2)>
            sage: a=EclObject(3)
            sage: L.rplaca(a)
            sage: L
            <ECL: (3 . 2)>

        """
        if not(bint_consp(self.obj)):
            raise TypeError("rplaca can only be applied to a cons")
        cl_rplaca(self.obj, d.obj)


    def rplacd(self,EclObject d):
        r"""
        Destructively replace cdr(self) with d.

        EXAMPLES::

            sage: from sage.libs.ecl import *
            sage: L=EclObject((1,2))
            sage: L
            <ECL: (1 . 2)>
            sage: a=EclObject(3)
            sage: L.rplacd(a)
            sage: L
            <ECL: (1 . 3)>

        """
        if not(bint_consp(self.obj)):
            raise TypeError("rplacd can only be applied to a cons")
        cl_rplacd(self.obj, d.obj)

    def car(self):
        r"""
        Return the car of self

        EXAMPLES::

            sage: from sage.libs.ecl import *
            sage: L=EclObject([[1,2],[3,4]])
            sage: L.car()
            <ECL: (1 2)>
            sage: L.cdr()
            <ECL: ((3 4))>
            sage: L.caar()
            <ECL: 1>
            sage: L.cadr()
            <ECL: (3 4)>
            sage: L.cdar()
            <ECL: (2)>
            sage: L.cddr()
            <ECL: NIL>
        """
        if not(bint_consp(self.obj)):
            raise TypeError("car can only be applied to a cons")
        return ecl_wrap(cl_car(self.obj))

    def cdr(self):
        r"""
        Return the cdr of self

        EXAMPLES::

            sage: from sage.libs.ecl import *
            sage: L=EclObject([[1,2],[3,4]])
            sage: L.car()
            <ECL: (1 2)>
            sage: L.cdr()
            <ECL: ((3 4))>
            sage: L.caar()
            <ECL: 1>
            sage: L.cadr()
            <ECL: (3 4)>
            sage: L.cdar()
            <ECL: (2)>
            sage: L.cddr()
            <ECL: NIL>
        """
        if not(bint_consp(self.obj)):
            raise TypeError("cdr can only be applied to a cons")
        return ecl_wrap(cl_cdr(self.obj))

    def caar(self):
        r"""
        Return the caar of self

        EXAMPLES::

            sage: from sage.libs.ecl import *
            sage: L=EclObject([[1,2],[3,4]])
            sage: L.car()
            <ECL: (1 2)>
            sage: L.cdr()
            <ECL: ((3 4))>
            sage: L.caar()
            <ECL: 1>
            sage: L.cadr()
            <ECL: (3 4)>
            sage: L.cdar()
            <ECL: (2)>
            sage: L.cddr()
            <ECL: NIL>
        """
        if not(bint_consp(self.obj) and bint_consp(cl_car(self.obj))):
            raise TypeError("caar can only be applied to a cons")
        return ecl_wrap(cl_caar(self.obj))

    def cadr(self):
        r"""
        Return the cadr of self

        EXAMPLES::

            sage: from sage.libs.ecl import *
            sage: L=EclObject([[1,2],[3,4]])
            sage: L.car()
            <ECL: (1 2)>
            sage: L.cdr()
            <ECL: ((3 4))>
            sage: L.caar()
            <ECL: 1>
            sage: L.cadr()
            <ECL: (3 4)>
            sage: L.cdar()
            <ECL: (2)>
            sage: L.cddr()
            <ECL: NIL>
        """
        if not(bint_consp(self.obj) and bint_consp(cl_cdr(self.obj))):
            raise TypeError("cadr can only be applied to a cons")
        return ecl_wrap(cl_cadr(self.obj))

    def cdar(self):
        r"""
        Return the cdar of self

        EXAMPLES::

            sage: from sage.libs.ecl import *
            sage: L=EclObject([[1,2],[3,4]])
            sage: L.car()
            <ECL: (1 2)>
            sage: L.cdr()
            <ECL: ((3 4))>
            sage: L.caar()
            <ECL: 1>
            sage: L.cadr()
            <ECL: (3 4)>
            sage: L.cdar()
            <ECL: (2)>
            sage: L.cddr()
            <ECL: NIL>
        """
        if not(bint_consp(self.obj) and bint_consp(cl_car(self.obj))):
            raise TypeError("cdar can only be applied to a cons")
        return ecl_wrap(cl_cdar(self.obj))

    def cddr(self):
        r"""
        Return the cddr of self

        EXAMPLES::

            sage: from sage.libs.ecl import *
            sage: L=EclObject([[1,2],[3,4]])
            sage: L.car()
            <ECL: (1 2)>
            sage: L.cdr()
            <ECL: ((3 4))>
            sage: L.caar()
            <ECL: 1>
            sage: L.cadr()
            <ECL: (3 4)>
            sage: L.cdar()
            <ECL: (2)>
            sage: L.cddr()
            <ECL: NIL>
        """
        if not(bint_consp(self.obj) and bint_consp(cl_cdr(self.obj))):
            raise TypeError("cddr can only be applied to a cons")
        return ecl_wrap(cl_cddr(self.obj))

    def fixnump(self):
        r"""
        Return True if self is a fixnum, False otherwise

        EXAMPLES::

            sage: from sage.libs.ecl import *
            sage: EclObject(2**3).fixnump()
            True
            sage: EclObject(2**200).fixnump()
            False

        """
        return bint_fixnump(self.obj)

    def characterp(self):
        r"""
        Return True if self is a character, False otherwise

        Strings are not characters

        EXAMPLES:

            sage: from sage.libs.ecl import *
            sage: EclObject('"a"').characterp()
            False

        """
        return bint_characterp(self.obj)

    def nullp(self):
        r"""
        Return True if self is NIL, False otherwise

        EXAMPLES::

            sage: from sage.libs.ecl import *
            sage: EclObject([]).nullp()
            True
            sage: EclObject([[]]).nullp()
            False
        """
        return bint_nullp(self.obj)

    def listp(self):
        r"""
        Return True if self is a list, False otherwise. NIL is a list.

        EXAMPLES::

            sage: from sage.libs.ecl import *
            sage: EclObject([]).listp()
            True
            sage: EclObject([[]]).listp()
            True
        """
        return bint_listp(self.obj)

    def consp(self):
        r"""
        Return True if self is a cons, False otherwise. NIL is not a cons.

        EXAMPLES::

            sage: from sage.libs.ecl import *
            sage: EclObject([]).consp()
            False
            sage: EclObject([[]]).consp()
            True
        """
        return bint_consp(self.obj)

    def atomp(self):
        r"""
        Return True if self is atomic, False otherwise.

        EXAMPLES::

            sage: from sage.libs.ecl import *
            sage: EclObject([]).atomp()
            True
            sage: EclObject([[]]).atomp()
            False

        """
        return bint_atomp(self.obj)

    def symbolp(self):
        r"""
        Return True if self is a symbol, False otherwise.

        EXAMPLES::

            sage: from sage.libs.ecl import *
            sage: EclObject([]).symbolp()
            True
            sage: EclObject([[]]).symbolp()
            False

        """
        return bint_symbolp(self.obj)

cdef class EclListIterator:
    r"""
    Iterator object for an ECL list

    This class is used to implement the iterator protocol for EclObject.
    Do not instantiate this class directly but use the iterator method
    on an EclObject instead. It is an error if the EclObject is not a list.

    EXAMPLES::

        sage: from sage.libs.ecl import *
        sage: I=EclListIterator(EclObject("(1 2 3)"))
        sage: type(I)
        <type 'sage.libs.ecl.EclListIterator'>
        sage: [i for i in I]
        [<ECL: 1>, <ECL: 2>, <ECL: 3>]
        sage: [i for i in EclObject("(1 2 3)")]
        [<ECL: 1>, <ECL: 2>, <ECL: 3>]
        sage: EclListIterator(EclObject("1"))
        Traceback (most recent call last):
        ...
        TypeError: ECL object is not iterable

    """
    cdef EclObject current

    def __init__(EclListIterator self, EclObject o):
        r"""
        Initialize EclListIterator

        EXAMPLES::

            sage: from sage.libs.ecl import *
            sage: I=EclListIterator(EclObject("(1 2 3)"))
            sage: type(I)
            <type 'sage.libs.ecl.EclListIterator'>

        """
        if not o.listp():
            raise TypeError("ECL object is not iterable")
        self.current = ecl_wrap(o.obj)

    def __iter__(EclListIterator self):
        r"""
        Return self

        It seems standard that iterators return themselves if asked to produce
        an iterator.

        EXAMPLES::

            sage: from sage.libs.ecl import *
            sage: I=EclListIterator(EclObject("(1 2 3)"))
            sage: id(I) == id(I.__iter__())
            True

        """
        return self

    def __next__(EclListIterator self):
        r"""
        Get next element from iterator

        EXAMPLES::

            sage: from sage.libs.ecl import *
            sage: I=EclListIterator(EclObject("(1 2 3)"))
            sage: next(I)
            <ECL: 1>
            sage: next(I)
            <ECL: 2>
            sage: next(I)
            <ECL: 3>
            sage: next(I)
            Traceback (most recent call last):
            ...
            StopIteration

        """

        if self.current.nullp():
            raise StopIteration
        elif self.current.consp():
            r = self.current.car()
            self.current = self.current.cdr()
        else:
            r = self.current
            self.current = ecl_wrap(Cnil)
        return r

#input: a cl-object. Output: EclObject wrapping that.
cdef EclObject ecl_wrap(cl_object o):
    cdef EclObject obj = EclObject.__new__(EclObject)
    obj.set_obj(o)
    return obj

#convenience routine to more easily evaluate strings
cpdef EclObject ecl_eval(bytes s):
    """
    Read and evaluate string in Lisp and return the result

    EXAMPLES::

        sage: from sage.libs.ecl import *
        sage: ecl_eval("(defun fibo (n)(cond((= n 0) 0)((= n 1) 1)(T (+ (fibo (- n 1)) (fibo (- n 2))))))")
        <ECL: FIBO>
        sage: ecl_eval("(mapcar 'fibo '(1 2 3 4 5 6 7))")
        <ECL: (1 1 2 3 5 8 13)>

    """
    cdef cl_object o
    o=ecl_safe_read_string(s)
    o=ecl_safe_eval(o)
    return ecl_wrap(o)

init_ecl()
