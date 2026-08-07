Environment variables used by Sage
==================================

Sage uses several environment variables when running.  These all have
sensible default values, so many users won't need to set any of these.
(There are also variables used to compile Sage; see the Sage
Installation Guide for more about those.)

.. envvar:: DOT_SAGE

  This is the directory, to which the user has
  read and write access, where Sage stores a number of files.  The
  default location is ``~/.sage/``, but you can change that by setting
  this variable.

.. envvar:: SAGE_RC_FILE

  This is a shell script which is sourced after
  Sage has determined its environment variables.  This script is
  executed before starting Sage or any of its subcommands (like
  ``sage -i <package>``).  The default value is
  :file:`$DOT_SAGE/sagerc`.

.. envvar:: SAGE_STARTUP_FILE

  This is a file including commands to be
  executed every time Sage starts.  The default value is
  :file:`$DOT_SAGE/init.sage`.

.. envvar:: SAGE_SERVER

  This is only used for installing
  packages. Alternative mirror from which to download sources, see the
  Installation Guide for details.

.. envvar:: BROWSER

  On most platforms, Sage will detect the
  command to run a web browser, but if this doesn't seem to work on
  your machine, set this variable to the appropriate command.

.. envvar:: TMPDIR

  This variable is used by Python, and hence by Sage, as the directory
  in which temporary files should be stored.

.. envvar:: PATH

  This is the shell search path used to find executable commands. Sage
  adjusts it for the Sage shell and subprocesses.

.. envvar:: EDITOR

  This selects the editor used by commands such as ``%edit``.

.. envvar:: TERM

  This describes the terminal type for command-line programs.

.. envvar:: SAGE_GAP_COMMAND

  This specifies how the GAP executable is called by the GAP interface.

.. envvar:: SAGE_GAP_MEMORY

  This specifies the amount of memory allocated to GAP by the GAP interface.

Relevant environment variables for other packages
=================================================

This is a non-exhaustive list of environment variables which influence
some package contained within the SageMath distribution.

In many cases, SageMath uses a custom default value if the variable is
not set, which is not the same default that the system-wide package
would use. So, if you would like to use your system-wide configuration,
you need to explicitly set the environment variable to the system-wide
default.

.. envvar:: IPYTHONDIR

  directory where the configuration of IPython is stored. By default,
  this is some directory inside :envvar:`DOT_SAGE`.
  See http://ipython.readthedocs.io/en/stable/development/config.html
  for more information.

.. envvar:: JUPYTER_CONFIG_DIR

  directory where the configuration of Jupyter is stored. By default,
  this is some directory inside :envvar:`DOT_SAGE`.
  See https://docs.jupyter.org/en/latest/use/jupyter-directories.html
  for more information.

.. envvar:: MPLCONFIGDIR

  directory where the configuration of Matplotlib is stored.
  See https://matplotlib.org/faq/environment_variables_faq.html#envvar-MPLCONFIGDIR
  By default, this is some directory inside :envvar:`DOT_SAGE`.
