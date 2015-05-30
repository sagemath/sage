# Run by ``sage-CSI`` in gdb's Python interpreter.

import os
import sys
import glob

import gdb
from Cython.Debugger import libpython, libcython
from Cython.Debugger.libcython import cy, CythonCommand, CyGDBError

try:
    if not color:
        libcython.pygments = None  # disable escape-sequence coloring
except (NameError, AttributeError):
    pass


def cython_debug_files():
    """
    Cython extra debug information files
    """
    pattern = os.path.join(os.environ['SAGE_SRC'], 'cython_debug',
                           'cython_debug_info_*')
    return glob.glob(pattern)



print('\n\n')
print('Cython backtrace')
print('----------------')

# The Python interpreter in GDB does not do automatic backtraces for you
try:

    for f in cython_debug_files():
        cy.import_.invoke(f, None)


    class SageBacktrace(CythonCommand):
        name = 'cy fullbt'
        alias = 'cy full_backtrace'
        command_class = gdb.COMMAND_STACK
        completer_class = gdb.COMPLETE_NONE
        cy = cy

        def print_stackframe(self, frame, index, is_c=False):
            if not is_c and self.is_python_function(frame):
                pyframe = libpython.Frame(frame).get_pyop()
                if pyframe is None or pyframe.is_optimized_out():
                    # print this python function as a C function
                    return self.print_stackframe(frame, index, is_c=True)
                func_name = pyframe.co_name
                func_cname = 'PyEval_EvalFrameEx'
                func_args = []
            elif self.is_cython_function(frame):
                cyfunc = self.get_cython_function(frame)
                f = lambda arg: self.cy.cy_cvalue.invoke(arg, frame=frame)

                func_name = cyfunc.name
                func_cname = cyfunc.cname
                func_args = [] # [(arg, f(arg)) for arg in cyfunc.arguments]
            else:
                func_name = frame.name()
                func_cname = func_name
                func_args = []

            try:
                gdb_value = gdb.parse_and_eval(func_cname)
            except (RuntimeError, TypeError):
                func_address = 0
            else:
                func_address = int(str(gdb_value.address).split()[0], 0)

            source_desc, lineno = self.get_source_desc(frame)
            a = ', '.join('%s=%s' % (name, val) for name, val in func_args)
            out = '#%-2d 0x%016x in %s(%s)' % (index, func_address, func_name, a)
            if source_desc.filename is not None:
                out += 'at %s:%s' % (source_desc.filename, lineno)
            print(out)
            try:
                source = source_desc.get_source(lineno - 5, lineno + 5,
                                                mark_line=lineno, lex_entire=True)
                print(source)
            except gdb.GdbError:
                pass


        def invoke(self, args, from_tty):
            self.newest_first_order(args, from_tty)


        def newest_first_order(self, args, from_tty):
            frame = gdb.newest_frame()
            index = 0
            while frame:
                frame.select()
                self.print_stackframe(frame, index)
                index += 1
                frame = frame.older()


        def newest_last_order(self, args, from_tty):
            frame = gdb.newest_frame()
            n_frames = 0;
            while frame.older():
                frame = frame.older()
                n_frames += 1
            index = 0
            while frame:
                frame.select()
                self.print_stackframe(frame, index)
                index += 1
                frame = frame.newer()


    trace = SageBacktrace.register()
    trace.invoke(None, None)


except Exception as e:
    import traceback
    traceback.print_exc()

sys.stdout.flush()
