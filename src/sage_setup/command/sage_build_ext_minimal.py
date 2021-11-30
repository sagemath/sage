import os
import multiprocessing
from setuptools.command.build_ext import build_ext


class sage_build_ext_minimal(build_ext):
    """
    In contrast to :func:`~sage_setup.sage_build_ext.sage_build_ext`, this build extension is designed
    to be used in combination with Cython's cythonize function.
    Thus, we only take care of some options and letting Cython do the main work.
    """

    def initialize_options(self):
        build_ext.initialize_options(self)
        self.parallel = self.get_default_number_build_jobs()

    @staticmethod
    def get_default_number_build_jobs() -> int:
        """
        Get number of parallel build jobs used by default, i.e. unless explicitly
        set by the --parallel command line argument of setup.py.

        First, the environment variable `SAGE_NUM_THREADS` is checked.
        If that is unset, return the number of processors on the system,
        with a maximum of 10 (to prevent overloading the system if there a lot of CPUs).

        OUTPUT:
            number of parallel jobs that should be run
        """
        try:
            cpu_count = len(os.sched_getaffinity(0))
        except AttributeError:
            cpu_count = multiprocessing.cpu_count()
        cpu_count = min(cpu_count, 10)
        return int(os.environ.get("SAGE_NUM_THREADS", cpu_count))
