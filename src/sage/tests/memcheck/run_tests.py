# sage_setup: distribution = sagemath-repl
import types


def run_tests() -> None:
    """
    Run all memcheck tests
    """
    from sage.tests.memcheck import symbolic_expression
    run_tests_in_module(symbolic_expression)


def run_tests_in_module(mod: types.ModuleType) -> None:
    """
    Run all memcheck tests in the given module
    """
    for entry in dir(mod):
        if not entry.startswith('test_'):
            continue
        test_func = getattr(mod, entry)
        test_func()


if __name__ == '__main__':
    run_tests()
