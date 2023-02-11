r"""
Benchmarks for SageMath

This module searches for doctest timings created with ``sage -t
--asv_stats_path=$SAGE_ASV_STATS`` if the ``$SAGE_ASV_STATS`` environment
variable is set. For each timing, it dynamically creates a method in a class in
this module to report that timing to airspeed velocity as configured in
``asv.conf.json``.

EXAMPLES:

Since this variable is usually not set, this module does nothing::

    sage: import sage.benchmark
    sage: dir(sage.benchmark)
    ['__builtins__',
     '__cached__',
     '__doc__',
     '__file__',
     '__loader__',
     '__name__',
     '__package__',
     '__path__',
     '__spec__',
     'create_benchmark_class',
     'create_track_method',
     'create_trackers_from_doctest_stats',
     'json',
     'os']

"""
# ****************************************************************************
#       Copyright (C) 2023 Julian Rueth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import json
import os

# TODO: Explain here what this is for and put the replacement logic in a central space
# TODO: Can we load a plugin that undoes this replacement on the client?
DOT = "․"


def create_trackers_from_doctest_stats():
    r"""
    Create asv-compatible benchmark classes that contain a `track_…` method
    for each doctest block reported in the log produced by ``sage -t
    --asv_stats_path=$SAGE_ASV_STATS``

    We export these wrappers in this module. These get picked up by asv which
    will then run these as "benchmarks". Instead of running the doctests to
    benchmark them, we are just reporting the timings that we collected in a
    previous doctest run instead.

    EXAMPLES:

    When ``SAGE_ASV_STATS`` is not set, this creates an empty dict::

        sage: from sage.benchmark import create_trackers_from_doctest_stats
        sage: create_trackers_from_doctest_stats()
        {}

    """
    stats_file = os.environ.get("SAGE_ASV_STATS", None)
    if not stats_file:
        return {}

    runs_by_module = json.load(open(stats_file))

    return {
        module.replace(".", DOT): create_benchmark_class(module, runs)
        for (module, runs) in runs_by_module.items()
    }


def create_benchmark_class(name, runs):
    r"""
    Return an ASV compatible benchmark class called `name` that pretends to
    perform the ``runs`` representing all the doctests of a module.

    EXAMPLES::

        sage: from sage.benchmark import create_benchmark_class
        sage: benchmark = create_benchmark_class("sage.rings.padics.padic_generic", [
        ....:     ["", [], "pass", 1.337],
        ....:     ["pAdicGeneric.some_elements", ["long", "meataxe"], "pass\npass", 13.37]
        ....: ])
        sage: dir(benchmark)
        [...
         'track_[sage.rings.padics.padic_generic]',
         'track_pAdicGeneric.some_elements[long,meataxe]']

    """
    class Benchmark:
        pass

    # asv groups entries by top-level package name, so everything goes into "sage".
    # To break this logic, we replace "." with a One Dot Leader character that
    # looks similar.
    Benchmark.__name__ = name.replace(".", DOT)

    for identifier, flags, source, time in runs:
        method = "track_" + (identifier or f"[{name}]").replace(".", DOT)
        if flags:
            method += "[" + ",".join(flags) + "]"

        if hasattr(Benchmark, method):
            # raise NotImplementedError(f"cannot merge duplicate results for {identifier} with {flags}")
            print("duplicate results")
            continue

        setattr(Benchmark, method, create_track_method(name, identifier, flags, source, time))

    return Benchmark


def create_track_method(module, identifier, flags, source, time):
    r"""
    Return a function that can be added as a method of a benchmark class.

    The method returns tracking results that will be picked up by asv as timing
    benchmark results.

    EXAMPLES::

        sage: from sage.benchmark import create_track_method
        sage: method = create_track_method("sage.rings.padics.padic_generic", "", [], "pass", 1.337)
        sage: method.pretty_name
        '[sage.rings.padics.padic_generic]'
        sage: method.pretty_source
        'pass'
        sage: method(None)
        1337.00000000000

    ::

        sage: method = create_track_method("sage.rings.padics.padic_generic", "pAdicGeneric.some_elements", ["long", "meataxe"], "pass\npass", 13.37)
        sage: method.pretty_name
        'pAdicGeneric.some_elements'
        sage: print(method.pretty_source)
        # with long, meataxe
        pass
        pass
        sage: method(None)
        13370.0000000000

    """
    def method(self):
        return time * 1e3

    method.unit = "ms"

    # We do not need to run this "benchmark" more than once; we are just
    # reporting the same static data every time.
    method.repeat = 1
    method.number = 1
    method.min_run_count = 1

    # The doctest of a module itself has no identifier, so we write the module name instead.
    method.pretty_name = identifier or f"[{module}]"

    method.pretty_source=source
    if flags:
        method.pretty_source = f"# with {', '.join(flags)}\n" + method.pretty_source

    return method


locals().update(create_trackers_from_doctest_stats())
