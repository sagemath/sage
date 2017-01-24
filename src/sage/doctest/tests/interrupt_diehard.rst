Save the current PID to the file given by :envvar:DOCTEST_TEST_PID_FILE::

    sage: open(os.environ['DOCTEST_TEST_PID_FILE'], "w").write(str(os.getpid()))

Interrupt the doctester (the parent process) while blocking the hangup
signal (used to kill this process)::

    sage: import signal
    sage: import time
    sage: from cysignals.pselect import PSelecter
    sage: with PSelecter([signal.SIGHUP]):
    ....:     os.kill(os.getppid(), signal.SIGINT)
    ....:     time.sleep(30)
