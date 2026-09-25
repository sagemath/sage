import subprocess
import sys
import textwrap

import pytest

from sage.numerical.backends.generic_backend import GenericBackend
from sage.numerical.backends.generic_backend_test import GenericBackendTests
from sage.numerical.backends.glpk_backend import _glpk_uses_thread_local_env
from sage.numerical.mip import MixedIntegerLinearProgram


def _run_glpk_subprocess(source):
    """Run a GLPK regression in a process that may safely crash."""
    proc = subprocess.run(
        [sys.executable, "-c", textwrap.dedent(source)],
        capture_output=True,
        text=True,
        check=False,
        timeout=30,
    )
    assert proc.returncode == 0, (
        f"GLPK subprocess exited with status {proc.returncode}\n"
        f"stdout (tail):\n{proc.stdout[-4000:]}\n"
        f"stderr (tail):\n{proc.stderr[-4000:]}"
    )


class TestGLPKBackend(GenericBackendTests):

    @pytest.fixture
    def backend(self) -> GenericBackend:
        return MixedIntegerLinearProgram(solver='GLPK').get_backend()


requires_thread_local_glpk = pytest.mark.skipif(
    not _glpk_uses_thread_local_env(),
    reason="the loaded GLPK uses one process-global environment",
)


@pytest.mark.skipif(
    _glpk_uses_thread_local_env(),
    reason="the loaded GLPK uses a thread-local environment",
)
def test_process_global_glpk_allows_serialized_thread_handoff():
    """Non-reentrant GLPK objects remain usable after their creator exits."""
    import threading

    from sage.numerical.backends.glpk_backend import GLPKBackend
    from sage.numerical.backends.glpk_graph_backend import GLPKGraphBackend
    from sage.numerical.backends.glpk_thread_resources import glpk_thread_resources

    state = {}

    def create_backends():
        state["lp"] = GLPKBackend()
        state["graph"] = GLPKGraphBackend()

    thread = threading.Thread(target=create_backends)
    thread.start()
    thread.join()

    assert state["lp"].add_variable() == 0
    state["graph"].add_vertex("A")
    assert state["graph"].vertices() == ["A"]
    assert glpk_thread_resources().resources == []


@requires_thread_local_glpk
@pytest.mark.parametrize(
    ("backend_import", "constructor", "access"),
    [
        (
            "from sage.numerical.backends.glpk_backend import GLPKBackend",
            "GLPKBackend()",
            "backend.add_variable()",
        ),
        (
            "from sage.numerical.backends.glpk_graph_backend import GLPKGraphBackend",
            "GLPKGraphBackend()",
            'backend.add_vertex("A")',
        ),
    ],
    ids=("lp", "graph"),
)
def test_registration_failure_does_not_leave_a_dangling_environment_cleaner(
        backend_import, constructor, access):
    """A failed registry append must not retain a live native pointer."""
    _run_glpk_subprocess(f"""
        {backend_import}
        from sage.numerical.backends.glpk_thread_resources import glpk_thread_resources

        class FailingList(list):
            def append(self, resource):
                raise MemoryError("injected registry append failure")

        resources = glpk_thread_resources()
        original_resources = resources.resources
        assert original_resources == []
        assert resources.environment_cleaner is None
        resources.resources = FailingList()
        try:
            try:
                {constructor}
            except MemoryError:
                pass
            else:
                raise AssertionError("backend creation unexpectedly succeeded")
        finally:
            resources.resources = original_resources

        assert resources.environment_cleaner is None
        backend = {constructor}
        {access}
    """)


@requires_thread_local_glpk
@pytest.mark.parametrize(
    ("backend_import", "constructor", "foreign_access", "owner_access"),
    [
        (
            "from sage.numerical.backends.glpk_backend import GLPKBackend",
            "GLPKBackend()",
            "backend.ncols()",
            "survivor.add_variable()",
        ),
        (
            "from sage.numerical.backends.glpk_graph_backend import GLPKGraphBackend",
            "GLPKGraphBackend()",
            "backend.vertices()",
            'survivor.add_vertex("owner")',
        ),
    ],
    ids=("lp", "graph"),
)
def test_cross_thread_release_is_swept_by_the_live_owner(
        backend_import, constructor, foreign_access, owner_access):
    """A live owner must collect a backend released by another thread."""
    _run_glpk_subprocess(f"""
        import gc
        import queue
        import threading

        {backend_import}
        from sage.numerical.backends.glpk_thread_resources import glpk_thread_resources

        transferred = queue.Queue()
        continue_owner = threading.Event()
        state = {{}}

        def run_owner():
            registry = glpk_thread_resources()
            backend = {constructor}
            survivor = {constructor}
            state["registry"] = registry
            transferred.put(backend)
            del backend
            assert continue_owner.wait(10)
            {owner_access}
            state["remaining_after_sweep"] = len(registry.resources)
            del survivor
            gc.collect()
            state["remaining_after_release"] = len(registry.resources)

        owner = threading.Thread(target=run_owner)
        owner.start()
        backend = transferred.get(timeout=10)

        try:
            {foreign_access}
        except RuntimeError as error:
            assert "different thread" in str(error)
        else:
            raise AssertionError("foreign thread accessed a GLPK backend")

        del backend
        gc.collect()
        assert state["registry"].has_released
        continue_owner.set()
        owner.join(10)
        assert not owner.is_alive()
        assert state["remaining_after_sweep"] == 1
        assert state["remaining_after_release"] == 0
    """)


@requires_thread_local_glpk
@pytest.mark.parametrize(
    ("backend_import", "constructor", "access"),
    [
        (
            "from sage.numerical.backends.glpk_backend import GLPKBackend",
            "GLPKBackend()",
            "backend.add_variable()",
        ),
        (
            "from sage.numerical.backends.glpk_graph_backend import GLPKGraphBackend",
            "GLPKGraphBackend()",
            "backend.vertices()",
        ),
    ],
    ids=("lp", "graph"),
)
def test_recycled_thread_ident_does_not_reauthorize_resource(
        backend_import, constructor, access):
    """
    A recycled numeric thread ID must not authorize a stale GLPK resource.

    Retaining the private cleanup guard deliberately keeps the old resource
    live after thread exit.  On platforms that recycle thread IDs, a
    replacement thread must still be distinguished from the original owner.
    In particular, its destructor must not pass the stale pointer to the
    replacement thread's GLPK allocator.
    """
    _run_glpk_subprocess(f"""
        import gc
        import threading

        {backend_import}
        from sage.numerical.backends import glpk_thread_resources as resources_module

        state = {{}}

        def create_backend():
            state["owner_ident"] = threading.get_ident()
            state["backend"] = {constructor}
            state["registry"] = resources_module.glpk_thread_resources()
            # Retain the private guard so its destructor cannot free and NULL
            # the resource before the native thread ID can be recycled.
            state["cleanup"] = resources_module._glpk_thread_data.cleanup

        thread = threading.Thread(target=create_backend)
        thread.start()
        thread.join()
        assert "backend" in state
        assert len(state["registry"].resources) == 1

        # The main thread is not the owner and must reject the live stale
        # resource before passing its pointer to GLPK.
        backend = state["backend"]
        try:
            {access}
        except RuntimeError as error:
            assert "different thread" in str(error)
        else:
            raise AssertionError("main thread accessed the stale backend")
        del backend

        def access_from_recycled_ident():
            if threading.get_ident() != state["owner_ident"]:
                return
            state["ident_was_recycled"] = True
            backend = state.pop("backend")
            try:
                {access}
            except RuntimeError as error:
                state["outcome"] = str(error)
            except BaseException as error:
                state["outcome"] = type(error).__name__
            else:
                state["outcome"] = "accepted"
            del backend
            gc.collect()

        for unused in range(100):
            thread = threading.Thread(target=access_from_recycled_ident)
            thread.start()
            thread.join()
            if state.get("ident_was_recycled"):
                break

        if state.get("ident_was_recycled"):
            assert "different thread" in state.get("outcome", ""), state.get("outcome")
    """)


def test_mincost_okalg_accepts_negative_one_optimum():
    """A valid objective of ``-1`` must not be mistaken for an exception."""
    _run_glpk_subprocess("""
        from sage.numerical.backends.glpk_graph_backend import GLPKGraphBackend

        graph = GLPKGraphBackend()
        graph.add_vertices(["source", "sink"])
        graph.set_vertex_demand("source", 1)
        graph.set_vertex_demand("sink", -1)
        graph.add_edge("source", "sink", {"cap": 1, "cost": -1})
        assert graph.mincost_okalg() == -1.0
    """)


def test_get_row_dual_rejects_invalid_indices():
    """Invalid row indices must raise in Python instead of aborting in GLPK."""
    _run_glpk_subprocess("""
        from sage.numerical.backends.glpk_backend import (
            GLPKBackend,
            glp_simplex_only,
            glp_simplex_or_intopt,
        )

        backend = GLPKBackend()
        backend.add_variable()
        backend.add_linear_constraint([(0, 1)], None, 1)
        backend.solver_parameter(glp_simplex_or_intopt, glp_simplex_only)

        for index in (-1, backend.nrows()):
            try:
                backend.get_row_dual(index)
            except ValueError:
                pass
            else:
                raise AssertionError(f"get_row_dual accepted index {index}")
    """)


@requires_thread_local_glpk
def test_unavailable_lp_accessors_raise_consistently():
    """LP accessors must not turn an unavailable backend into valid values."""
    _run_glpk_subprocess("""
        import threading

        from sage.numerical.mip import MixedIntegerLinearProgram

        state = {}

        def create_backend():
            state["program"] = MixedIntegerLinearProgram(
                maximization=True, solver="GLPK"
            )

        thread = threading.Thread(target=create_backend)
        thread.start()
        thread.join()
        program = state["program"]
        backend = program.get_backend()

        operations = (
            ("ncols", backend.ncols),
            ("nrows", backend.nrows),
            ("is_maximization", backend.is_maximization),
            ("is_variable_binary", lambda: backend.is_variable_binary(0)),
            ("is_variable_integer", lambda: backend.is_variable_integer(0)),
            ("is_variable_continuous", lambda: backend.is_variable_continuous(0)),
            ("is_variable_basic", lambda: backend.is_variable_basic(0)),
            (
                "is_variable_nonbasic_at_lower_bound",
                lambda: backend.is_variable_nonbasic_at_lower_bound(0),
            ),
            ("is_slack_variable_basic", lambda: backend.is_slack_variable_basic(0)),
            (
                "is_slack_variable_nonbasic_at_lower_bound",
                lambda: backend.is_slack_variable_nonbasic_at_lower_bound(0),
            ),
            ("number_of_variables", program.number_of_variables),
            ("number_of_constraints", program.number_of_constraints),
            ("get_row_dual", lambda: backend.get_row_dual(0)),
            ("warm_up", backend.warm_up),
        )
        for name, operation in operations:
            try:
                operation()
            except RuntimeError:
                pass
            except Exception as error:
                raise AssertionError(
                    f"{name} raised {type(error).__name__}, not RuntimeError"
                ) from error
            else:
                raise AssertionError(f"{name} accepted an unavailable backend")
    """)


@requires_thread_local_glpk
def test_unavailable_graph_accessors_raise_consistently():
    """Graph lookup helpers must preserve backend-unavailable exceptions."""
    _run_glpk_subprocess("""
        import threading

        from sage.numerical.backends.glpk_graph_backend import GLPKGraphBackend

        state = {}

        def create_backend():
            state["backend"] = GLPKGraphBackend()

        thread = threading.Thread(target=create_backend)
        thread.start()
        thread.join()
        backend = state["backend"]

        operations = (
            ("get_vertex", lambda: backend.get_vertex("A")),
            ("get_edge", lambda: backend.get_edge("A", "B")),
            ("delete_edge", lambda: backend.delete_edge("A", "B")),
            ("set_vertex_demand", lambda: backend.set_vertex_demand("A", 1)),
            ("maxflow_ffalg", lambda: backend.maxflow_ffalg("A", "B")),
        )
        for name, operation in operations:
            try:
                operation()
            except RuntimeError:
                pass
            except Exception as error:
                raise AssertionError(
                    f"{name} raised {type(error).__name__}, not RuntimeError"
                ) from error
            else:
                raise AssertionError(f"{name} accepted an unavailable backend")
    """)
