ptyprocess: Python interaction with subprocesses in a pseudoterminal
====================================================================

Description
-----------

Launch a subprocess in a pseudo terminal (pty), and interact with both
the process and its pty.

Sometimes, piping stdin and stdout is not enough. There might be a
password prompt that doesn't read from stdin, output that changes when
it's going to a pipe rather than a terminal, or curses-style interfaces
that rely on a terminal. If you need to automate these things, running
the process in a pseudo terminal (pty) is the answer.

License
-------

Ptyprocess is under the ISC license, as code derived from Pexpect.

   http://opensource.org/licenses/ISC


Upstream Contact
----------------

https://github.com/pexpect/ptyprocess

