Sage startup scripts
====================

There are two kinds of startup scripts that Sage reads when starting:

The sagerc shell script
-----------------------

The *bash shell script* :file:`$DOT_SAGE/sagerc` (with the default
value of :envvar:`DOT_SAGE`, this is :file:`~/.sage/sagerc`) is read
by :sage_root:`src/bin/sage-env` after Sage has set its
environment variables.
It can be used to override some of the environment variables determined
by Sage, or it can contain other shell commands like creating
directories.
This script is sourced not only when running Sage itself, but also when
running any of the subcommands (like ``sage --python``, ``sage -b`` or
``sage -i <package>``).
In particular, setting ``PS1`` here overrides the default prompt for
the Sage shells ``sage --buildsh`` and ``sage --sh``.

.. note::

  This script is run with the Sage directories in its :envvar:`PATH`,
  so executing ``git`` for example will run the Git inside Sage.

The default location of this file can be changed using the
environment variable :envvar:`SAGE_RC_FILE`.

The init.sage script
--------------------

The *Sage script* :file:`$DOT_SAGE/init.sage` (with the default
value of :envvar:`DOT_SAGE`, this is :file:`~/.sage/init.sage`)
contains Sage commands to be executed every time Sage starts.
If you want symbolic variables ``y`` and ``z`` in every Sage session,
you could put ::

    var('y, z')

in this file.

The default location of this file can be changed using the
environment variable :envvar:`SAGE_STARTUP_FILE`.
