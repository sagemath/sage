# sage_setup: distribution = sagemath-repl
from typing import Callable, Any
import valgrind


def verify_no_leak(callback: Callable[[], Any],
                   repeat: int = 10000,
                   fuzzy: int = 10,
                   ) -> None:
    """
    Verify that the callback does not generate new definitely lost blocks

    Raises an assertion if the callback leaks memory
    """
    callback()   # warm_up
    initial_blocks = (0, 0, 0, 0)
    valgrind.memcheck_do_leak_check()
    initial_blocks = valgrind.memcheck_count_leak_blocks()
    for _ in range(repeat):
        callback()
    valgrind.memcheck_do_leak_check()
    leak_blocks = valgrind.memcheck_count_leak_blocks()
    leak = leak_blocks[0] - initial_blocks[0]
    if leak < repeat - fuzzy:
        return  # callback did not leak at least once per call
    blocks = round(leak / repeat, 2)
    message = f'{callback} leaked {blocks} block on average ({repeat} iterations)'
    raise AssertionError(message)
