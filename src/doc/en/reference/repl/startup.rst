Sage startup scripts
====================

Sage can be configured through environment variables and the
:file:`init.sage` script for interactive sessions.

Environment variables
---------------------

Set environment variables in the environment used to launch Sage.
For persistent settings, add ``export`` commands to the appropriate startup
file for your shell, such as :file:`~/.bashrc` for interactive, non-login Bash
shells. If another application launches Sage, use that application's
environment configuration. These settings are inherited by Sage when it starts.

Variables needed during initialization must be set before starting Sage
or the Jupyter process that will launch its kernel. Setting them in
:file:`init.sage` is too late for settings already read during initialization.

The init.sage script
--------------------

The *Sage script* :file:`$DOT_SAGE/init.sage` (with the default
value of :envvar:`DOT_SAGE`, this is :file:`~/.sage/init.sage`)
contains Sage commands to be executed when an interactive Sage session
starts, including a Sage Jupyter kernel.
If you want symbolic variables ``y`` and ``z`` in every interactive Sage session,
you could put ::

    var('y, z')

in this file.

The default location of this file can be changed using the
environment variable :envvar:`SAGE_STARTUP_FILE`.

This file is suitable for defining variables and functions and configuring
the IPython session. It runs after the Sage library has been imported.
It is not read when the Python command-line launcher executes ``sage -c``
or runs a script.

The sagerc shell script
-----------------------

The Python command-line launcher does not read :file:`sagerc` and ignores
:envvar:`SAGE_RC_FILE`. Move environment settings from this file to the
environment used to launch Sage, as described above.

Only the legacy Bash launcher :sage_root:`src/bin/sage` reads this Bash
script, through :sage_root:`src/bin/sage-env`, after setting its environment
variables. Its default location is :file:`$DOT_SAGE/sagerc`
(:file:`~/.sage/sagerc` by default), and can be changed using
:envvar:`SAGE_RC_FILE`.

.. _sage_subcommands:

Sage Subcommands
----------------

The legacy Bash launcher provides several subcommands, including:

* ``sage -b``: Rebuilds the Sage library. This is intended for developer-only
  editable Meson builds. It runs ``ninja -C build`` to recompile. If a
  non-Meson build is detected, it will error out with a message.

* ``sage -br``: **Deprecated.** This command is no longer supported and
  will be removed in a future release. Use ``sage -b`` to build and
  ``./sage`` to run.

* ``sage --python``: Runs the Python interpreter included with Sage.

* ``sage -i <package>``: **Deprecated.** This command is obsolete. Use ``./configure --with-<package>`` followed by ``make`` instead.
