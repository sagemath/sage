.. _sec-post-installation:

Post installation
=================

Now we assume that you installed Sage properly on your system. This
section quickly explains how to start the Sage console and the Jupyter
Notebook from the command line.

If you did install the Windows version or the macOS application you
should have icons available on your desktops or launching menus. Otherwise
you are strongly advised to create shortcuts for Sage as indicated in the part
6 of the "Installation steps" Section in :ref:`build-from-source-step-by-step`.
Assuming that you have this shortcut, running

.. CODE-BLOCK:: bash

    sage

in a console starts a Sage session.  To quit the session enter ``quit`` and
then press :kbd:`Enter`.

To start a Jupyter Notebook instead of a Sage console, run the command

.. CODE-BLOCK:: bash

    sage -n jupyter

instead of just ``sage``. To quit the Jupyter Notebook press :kbd:`Ctrl` + :kbd:`c`
twice in the console where you launched the command.

You can pass extra parameters to this command. For example,

.. CODE-BLOCK:: bash

    sage -n jupyter --port 8899

will run the Jupyter server on a port different from the default (8888).
In particular on WSL, this is very useful because Jupyter may not be able to
detect whether the default port is already taken by another instance of
Jupyter running in Windows.


Run-time environment variables
------------------------------

Sage uses the following environment variables when it runs:

- :envvar:`DOT_SAGE` - this is the directory, to which the user has read and
  write access, where Sage stores a number of files.
  The default location is :file:`$HOME/.sage/`.

- :envvar:`SAGE_STARTUP_FILE` - a file including commands to be executed every
  time Sage starts.
  The default value is :file:`$DOT_SAGE/init.sage`.

- :envvar:`BROWSER` - on most platforms, Sage will detect the command to
  run a web browser, but if this doesn't seem to work on your machine, set this
  variable to the appropriate command.

- :envvar:`TMPDIR` - this variable is used by Python, and hence by
  Sage; it gives the directory in which temporary files should be
  stored. This includes files used by the notebook. Some browsers have
  security settings which restrict the locations of files that they
  will access, and users may need to set this variable to handle this
  situation.

- See
  https://docs.python.org/3/using/cmdline.html#environment-variables
  for more variables used by Python (not an exhaustive list). With
  Python 3.11 or later, a brief summary can also be obtained by
  running ``python3 --help-env``.

Using a Jupyter Notebook remotely
---------------------------------

If Sage is installed on a remote machine to which you have ``ssh`` access, you
can launch a Jupyter Notebook using a command such as

.. CODE-BLOCK:: bash

    ssh -L localhost:8888:localhost:8888 -t USER@REMOTE sage -n jupyter --no-browser --port=8888

where ``USER@REMOTE`` needs to be replaced by the login details to the remote
machine. This uses local port forwarding to connect your local machine to the
remote one. The command will print a URL to the console which you can copy and
paste in a web browser.

Note that this assumes that a firewall which might be present between server
and client allows connections on port 8888. See details on port forwarding on
the internet, e.g. https://www.ssh.com/ssh/tunneling/example.


.. _sec-post-installation-wsl:

WSL post-installation steps
---------------------------

If you've installed Sage from source on WSL, there are a couple of extra steps you can do to make your life easier:


Create a notebook launch script
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you plan to use JupyterLab, install that first.

Now create a script called ``~/sage_nb.sh`` containing the following lines, and fill in the correct paths for your desired starting directory and ``SAGE_ROOT``


.. CODE-BLOCK:: bash

    #!/bin/bash
    # Switch to desired windows directory
    cd /mnt/c/path/to/desired/starting/directory
    # Start the Jupyter notebook
    SAGE_ROOT/sage --notebook
    # Alternatively you can run JupyterLab - delete the line above, and uncomment the line below
    #SAGE_ROOT/sage --notebook jupyterlab

Make it executable:

.. CODE-BLOCK:: bash

    chmod ug+x ~/sage_nb.sh

Run it to test:

.. CODE-BLOCK:: bash

    cd ~
    ./sage_nb.sh

The Jupyter(Lab) server should start in the terminal window, and you windows browser should open a page showing the Jupyter or JupyterLab starting page, at the directory you specified.

Create a shortcut
^^^^^^^^^^^^^^^^^

This is a final nicety that lets you start the Jupyter or JupyterLab server in one click:

* Open Windows explorer, and type ``%APPDATA%\Microsoft\Windows\Start Menu\Programs`` in the address bar and press enter. This is the folder that contains you start menu shortcuts. If you want the sage shortcut somewhere else (like your desktop), open that folder instead.
* Open a separate window and go to ``%LOCALAPPDATA%\Microsoft\WindowsApps\``
* Right-click-drag the ``ubuntu.exe`` icon from the second window into the first, then choose ``Create shortcuts here`` from the context menu when you drop it.
* To customize this shortcut, right-click on it and choose properties.

  * On the General tab:

    * Change the name to whatever you want, e.g. "Sage 9.2 JupyterLab"

  * On the Shortcut tab:

    * Change Target to: ``ubuntu.exe run ~/sage_nb.sh``
    * Change Start in to: ``%USERPROFILE%``
    * Change Run to: Minimised
    * Change the icon if you want

Now hit the start button or key and type the name you gave it. it should appear in the list, and should load the server and fire up your browser when you click on it.


.. _sec-post-installation-system-jupyter:

Setting up Sage as a Jupyter kernel in an existing Jupyter installation
-----------------------------------------------------------------------

You may already have a global installation of Jupyter.  For added
convenience, it is possible to link your installation of Sage into

your Jupyter installation, adding it to the list of available kernels
that can be selected in the notebook or JupyterLab interface.

Assuming that Sage can be invoked by typing ``sage``, you can use

.. CODE-BLOCK:: bash

    sage -sh -c 'ls -d $SAGE_VENV/share/jupyter/kernels/sagemath'

to find the location of the Sage kernel description.

Now pick a name for the kernel that identifies it clearly and uniquely.

For example, if you install Sage from source tarballs, you could decide
to include the version number in the name, such as ``sagemath-9.6``.
If you build Sage from a clone of the git repository, it is better to
choose a name that identifies the directory, perhaps ``sagemath-dev``
or ``sagemath-teaching`` because the version will change.

Now assuming that the Jupyter notebook can be started by typing
``jupyter notebook``, the following command will install Sage as a
new kernel named ``sagemath-dev``.

.. CODE-BLOCK:: bash

    jupyter kernelspec install --user $(sage -sh -c 'ls -d $SAGE_VENV/share/jupyter/kernels/sagemath') --name sagemath-dev

The ``jupyter kernelspec`` approach by default does lead to about 2Gb of
Sage documentation being copied into your personal jupyter configuration
directory. You can avoid that by instead putting a symlink in the relevant spot
and

.. CODE-BLOCK:: bash

    jupyter --paths

to find valid data directories for your Jupyter installation.
A command along the lines of

.. CODE-BLOCK:: bash

    ln -s $(sage -sh -c 'ls -d $SAGE_VENV/share/jupyter/kernels/sagemath') $HOME/.local/share/jupyter/kernels/sagemath-dev

can then be used to create a symlink to the Sage kernel description
in a location where your own ``jupyter`` can find it.

If you have installed Sage from source, the alternative command

.. CODE-BLOCK:: bash

    ln -s $(sage -sh -c 'ls -d $SAGE_ROOT/venv/share/jupyter/kernels/sagemath') $HOME/.local/share/jupyter/kernels/sagemath-dev

creates a symlink that will stay current even if you switch to a different Python version
later.

To get the full functionality of the Sage kernel in your global
Jupyter installation, the following Notebook Extension packages also
need to be installed (or linked) in the environment from which the
Jupyter installation runs.

You can check the presence of some of these packages using the command
``jupyter nbextension list``.

- For the Sage interacts, you will need the package
  ``widgetsnbextension`` installed in the Python environment of the
  Jupyter installation.  If your Jupyter installation is coming from
  the system package manager, it is best to install
  ``widgetsnbextension`` in the same way.  Otherwise, install it
  using ``pip``.

  To verify that interacts work correctly, you can evaluate the following code
  in the notebook::

    @interact
    def _(k=slider(vmin=-1.0, vmax= 3.0, step_size=0.1, default=0), auto_update=True):
    plot([lambda u:u^2-1, lambda u:u+k], (-2,2),
         ymin=-1, ymax=3, fill={1:[0]}, fillalpha=0.5).show()

- For 3D graphics using Three.js, by default, internet connectivity
  is needed, as Sage's custom build of the Javascript package
  Three.js is retrieved from a content delivery network.

  To verify that online 3D graphics with Three.js works correctly,
  you can evaluate the following code in the notebook::

    plot3d(lambda u,v:(u^2+v^2)/4-2,(-2,2),(-2,2)).show()

  However, it is possible to configure graphics with Three.js for
  offline use.  In this case, the Three.js installation from the Sage
  distribution needs to be made available in the environment of the
  Jupyter installation.  This can be done by copying or symlinking.
  The Three.js installation in the environment of the Jupyter
  installation must exactly match the version that comes from the
  Sage distribution.  It is not supported to use several Jupyter
  kernels corresponding to different versions of the Sage distribution.

  To verify that offline 3D graphics with Three.js works correctly,
  you can evaluate the following code in the notebook::

    plot3d(lambda u,v:(u^2+v^2)/4-2,(-2,2),(-2,2), online=False).show()

- For 3D graphics using jsmol, you will need the package
  ``jupyter-jsmol`` installed in the Python environment of the
  Jupyter installation. You can install it using ``pip``.
  (Alternatively, you can copy or symlink it.)

  To verify that jsmol graphics work correctly, you can evaluate the
  following code in the notebook::

    plot3d(lambda u,v:(u^2+v^2)/4-2,(-2,2),(-2,2)).show(viewer="jmol")


Using Jupyter notebook through Visual Studio Code in WSL
--------------------------------------------------------

If you have installed Sage on Windows using Windows Subsystem for
Linux (WSL), it is convenient to use Visual Studio Code (VS Code)
to interact with Sage.

Here are steps to use Sage in a Jupyter notebook in VS Code:

* Install and run `VS Code <https://code.visualstudio.com/download>`_
  in Windows.

* Click the "Extension" icon on the left (or press :kbd:`Ctrl` +
  :kbd:`Shift` + :kbd:`X`) to open a list of extensions. Install the
  "Remote - WSL" and "Jupyter" extensions.

* In the command palette (:kbd:`Ctrl` + :kbd:`Shift` + :kbd:`P`),
  enter "Remote-WSL: New Window", and hit :kbd:`Enter`.

* In the command palette, enter "Create: New Jupyter Notebook", and
  hit :kbd:`Enter`.

* Click "Select Kernel" on the right (or press :kbd:`Ctrl` +
  :kbd:`Alt` + :kbd:`Enter`), select Sage, and hit :kbd:`Enter`.
