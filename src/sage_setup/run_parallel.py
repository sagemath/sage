########################################################################
##
## Parallel gcc execution
##
## This code is responsible for making distutils dispatch the calls to
## build_ext in parallel. Since distutils doesn't seem to do this by
## default, we create our own extension builder and override the
## appropriate methods.  Unfortunately, in distutils, the logic of
## deciding whether an extension needs to be recompiled and actually
## making the call to gcc to recompile the extension are in the same
## function. As a result, we can't just override one function and have
## everything magically work. Instead, we split this work between two
## functions. This works fine for our application, but it means that
## we can't use this modification to make the other parts of Sage that
## build with distutils call gcc in parallel.
##
########################################################################

import os
import sys
import time

keep_going = False

def run_command(cmd):
    """
    INPUT:

    - ``cmd`` -- a string; a command to run

    OUTPUT: prints ``cmd`` to the console and then runs
    ``os.system(cmd)``.
    """
    print(cmd)
    sys.stdout.flush()
    return os.system(cmd)

def apply_func_progress(p):
    """
    Given a triple p consisting of a function, value and a string,
    output the string and apply the function to the value.

    The string could for example be some progress indicator.

    This exists solely because we can't pickle an anonymous function
    in execute_list_of_commands_in_parallel below.
    """
    sys.stdout.write(p[2])
    sys.stdout.flush()
    return p[0](p[1])

def execute_list_of_commands_in_parallel(command_list, nthreads):
    """
    Execute the given list of commands, possibly in parallel, using
    ``nthreads`` threads.  Terminates ``setup.py`` with an exit code
    of 1 if an error occurs in any subcommand.

    INPUT:

    - ``command_list`` -- a list of commands, each given as a pair of
       the form ``[function, argument]`` of a function to call and its
       argument

    - ``nthreads`` -- integer; number of threads to use

    WARNING: commands are run roughly in order, but of course successive
    commands may be run at the same time.
    """
    # Add progress indicator strings to the command_list
    N = len(command_list)
    progress_fmt = "[{:%i}/{}] " % len(str(N))
    for i in range(N):
        progress = progress_fmt.format(i+1, N)
        command_list[i] = command_list[i] + (progress,)

    from multiprocessing import Pool
    # map_async handles KeyboardInterrupt correctly if an argument is
    # given to get().  Plain map() and apply_async() do not work
    # correctly, see Trac #16113.
    pool = Pool(nthreads)
    result = pool.map_async(apply_func_progress, command_list, 1).get(99999)
    pool.close()
    pool.join()
    process_command_results(result)

def process_command_results(result_values):
    error = None
    for r in result_values:
        if r:
            print("Error running command, failed with status %s."%r)
            if not keep_going:
                sys.exit(1)
            error = r
    if error:
        sys.exit(1)

def execute_list_of_commands(command_list):
    """
    INPUT:

    - ``command_list`` -- a list of strings or pairs

    OUTPUT:

    For each entry in command_list, we attempt to run the command.
    If it is a string, we call ``os.system()``. If it is a pair [f, v],
    we call f(v).

    If the environment variable :envvar:`SAGE_NUM_THREADS` is set, use
    that many threads.
    """
    t = time.time()
    # Determine the number of threads from the environment variable
    # SAGE_NUM_THREADS, which is set automatically by sage-env
    try:
        nthreads = int(os.environ['SAGE_NUM_THREADS'])
    except KeyError:
        nthreads = 1

    # normalize the command_list to handle strings correctly
    command_list = [ [run_command, x] if isinstance(x, str) else x for x in command_list ]

    # No need for more threads than there are commands, but at least one
    nthreads = min(len(command_list), nthreads)
    nthreads = max(1, nthreads)

    def plural(n,noun):
        if n == 1:
            return "1 %s"%noun
        return "%i %ss"%(n,noun)

    print("Executing %s (using %s)"%(plural(len(command_list),"command"), plural(nthreads,"thread")))
    execute_list_of_commands_in_parallel(command_list, nthreads)
    print("Time to execute %s: %.2f seconds."%(plural(len(command_list),"command"), time.time() - t))
