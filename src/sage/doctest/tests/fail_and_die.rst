The :exc:`NameError` raised on the second line should be displayed, even
if we crash immediately afterwards (also test printing of line continuation)::

    sage: import time, signal
    sage: print(1,
    ....:       2)
    1 2
    sage: this_gives_a_NameError
    sage: os.kill(os.getpid(), signal.SIGKILL)
