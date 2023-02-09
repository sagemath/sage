from sage.all import *

from sage.doctest.control import DocTestDefaults, DocTestController
from sage.doctest.forker import SageDocTestRunner, DocTestTask
from sage.doctest.parsing import parse_optional_tags

import timeit
import doctest

DEFAULTS = DocTestDefaults()
DEFAULTS.serial = True
DEFAULTS.long = True

PREFIX = 'track__'

def myglob(path, pattern):
    # python 2 does not have support for ** in glob patterns
    import fnmatch
    import os
    
    matches = []
    for root, dirnames, filenames in os.walk(path):
        for filename in fnmatch.filter(filenames, pattern):
            matches.append(os.path.join(root, filename))
    return matches

class BenchmarkMetaclass(type):
    _dir = None

    def _run(cls, fname=[SAGE_SRC + '/sage/']):
        old_runner = DocTestTask.runner
        DocTestTask.runner = BenchmarkRunner
        try:
            DocTestController(DEFAULTS, fname).run()
        finally:
            DocTestTask.runner = old_runner

    def __dir__(cls):
        if cls._dir is None:
            cls._dir = set()
            BenchmarkRunner._dir = cls._dir
            cls._run()
        return list(cls._dir)

    def __getattr__(cls, name):
        if not name.startswith(PREFIX):
            raise AttributeError
        cls.create_timer(name)
        return getattr(cls, name)

    def create_timer(cls, name):
        try:
            type.__getattr__(cls, name)
        except AttributeError:
            def time_doctest(self):
                BenchmarkRunner._selected = name
                BenchmarkRunner._time = 0
                cls._run(myglob(os.path.join(SAGE_SRC,'sage'), BenchmarkRunner.decode(name)+"*.*"))
                return BenchmarkRunner._time

            time_doctest.__name__ = name
            setattr(cls, name, time_doctest)

class Benchmarks(object):
    __metaclass__ = BenchmarkMetaclass

    def __getattr__(self, name):
        if not name.startswith(PREFIX):
            raise AttributeError
        type(self).create_timer(name)
        return getattr(self, name)

class BenchmarkRunner(SageDocTestRunner):
    _selected = None
    _dir = set()
    _time = None

    @classmethod
    def encode(cls, prefix, filename, name, digest):
        module = os.path.splitext(os.path.basename(filename))[0]
        method = name.split('.')[-1]
        return PREFIX+module+"__"+method+"__"+digest

    @classmethod
    def decode(cls, name):
        # This is not the full file name but only the bit up to the first _.
        # But that's good enough as we later seach for this with a glob pattern
        return name.split("__")[1]

    def run(self, test, clear_globs=True, *args, **kwargs):
        self._do_not_run_tests = True
        super(BenchmarkRunner, self).run(test, *args, clear_globs=False, **kwargs)

        name = BenchmarkRunner.encode(PREFIX, test.filename, test.name, self.running_doctest_digest.hexdigest())
        if name not in BenchmarkRunner._dir:
            for example in test.examples:
                if isinstance(example, doctest.Example):
                    if "long time" in parse_optional_tags(example.source):
                        BenchmarkRunner._dir.add(name)
                        break

        self._do_not_run_tests = False
        if type(self)._selected is not None:
            if name == type(self)._selected:
                pass
            else:
                self._do_not_run_tests = True
            if not self._do_not_run_tests:
                super(BenchmarkRunner, self).run(test, *args, clear_globs=clear_globs, **kwargs)
                return 

    def compile_and_execute(self, example, compiler, globs, *args, **kwargs):
        if self._do_not_run_tests:
            compiler = lambda example: compile("print(%s.encode('utf8'))"%(repr(example.want),), "want.py", "single")
            super(BenchmarkRunner, self).compile_and_execute(example, compiler, globs, *args, **kwargs)
        else:
            compiled = compiler(example)
            # TODO: According to https://asv.readthedocs.io/en/latest/writing_benchmarks.html, time.process_time would be better here but it's not available in Python 2
            start = timeit.default_timer()
            exec(compiled, globs)
            end = timeit.default_timer()
            type(self)._time += end-start

class BenchmarkRunningRunner(SageDocTestRunner):
    def run(self, test, *args, **kwargs):
        print(test,BenchmarkRunningRunner._test)
        if test == BenchmarkRunningRunner._test:
            super(BenchmarkRunningRunner).run(test, *args, **kwargs)
