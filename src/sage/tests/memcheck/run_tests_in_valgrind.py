"""
Launch valgrind and run the memory leak tests


From the commandline, run

    sage -python -m sage.tests.memcheck.run_tests_in_valgrind

to launch valgrind and execute the memory leak tests. Requires valgrind
to be installed. Alternatively, run as a unit test:

    sage: from sage.tests.memcheck.run_tests_in_valgrind import run_tests_in_valgrind
    sage: run_tests_in_valgrind()    # optional - valgrind
"""

import subprocess


def run_tests_in_valgrind() -> None:
    """
    Run the sage.tests.memcheck.run_tests module inside valgrind
    """
    subprocess.check_call([
        'valgrind',
        '--suppressions=src/sage/ext_data/valgrind/valgrind-python.supp',
        '--show-possibly-lost=no',
        '--show-reachable=no',
        './venv/bin/python',
        '-m',
        'sage.tests.memcheck.run_tests'
    ])


if __name__ == '__main__':
    run_tests_in_valgrind()
