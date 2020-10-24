Save the current PID to the file given by :envvar:DOCTEST_TEST_PID_FILE::

    sage: with open(os.environ['DOCTEST_TEST_PID_FILE'], "w") as file:
    ....:     file.write(str(os.getpid()))

Interrupt the doctester (the parent process) while blocking the quit
signal (used to kill this process)::

    sage: import signal
    sage: import time
    sage: from cysignals.pselect import PSelecter
    sage: with PSelecter([signal.SIGQUIT]):
    ....:     os.kill(os.getppid(), signal.SIGINT)
    ....:     time.sleep(30)
