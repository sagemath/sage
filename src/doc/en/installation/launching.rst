.. _sec-launching:

Launching SageMath
==================

Now we assume that you installed SageMath properly on your system. This
section quickly explains how to start the Sage console and the Jupyter
Notebook from the command line.

If you did install the Windows version or the macOS application you
should have icons available on your desktops or launching menus. Otherwise
you are strongly advised to create shortcuts for Sage as indicated in the part
6 of the "Installation steps" Section in :ref:`build-from-source-step-by-step`.
Assuming that you have this shortcut, running

.. code-block:: console

    $ sage

in a console starts a Sage session.  To quit the session enter ``quit`` and
then press ``<Enter>``.

To start a Jupyter Notebook instead of a Sage console, run the command

.. code-block:: console

    $ sage -n jupyter

instead of just ``sage``. To quit the Jupyter Notebook press ``<Ctrl> + <c>``
twice in the console where you launched the command.

You can pass extra parameters to this command. For example,

.. code-block:: console

    $ sage -n jupyter --port 8899

will run the Jupyter server on a port different from the default (8888).
In particular on WSL, this is very useful because Jupyter may not be able to
detect whether the default port is already taken by another instance of
Jupyter running in Windows.


Environment variables
---------------------

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
  for more variables used by Python (not an exhaustive list).
  A brief summary can also be obtained by running ``python3 --help-env``.

Using a Jupyter Notebook remotely
---------------------------------

If Sage is installed on a remote machine to which you have ``ssh`` access, you
can launch a Jupyter Notebook using a command such as

.. code-block:: console

    $ ssh -L localhost:8888:localhost:8888 -t USER@REMOTE sage -n jupyter --no-browser --port=8888

where ``USER@REMOTE`` needs to be replaced by the login details to the remote
machine. This uses local port forwarding to connect your local machine to the
remote one. The command will print a URL to the console which you can copy and
paste in a web browser.

Note that this assumes that a firewall which might be present between server
and client allows connections on port 8888. See details on port forwarding on
the internet, e.g. https://www.ssh.com/ssh/tunneling/example.

.. _sec-launching-wsl-post-installation:

WSL Post-installation steps
---------------------------

If you've installed SageMath from source on WSL, there are a couple of extra steps you can do to make your life easier:


Create a notebook launch script
"""""""""""""""""""""""""""""""

If you plan to use JupyterLab, install that first.

Now create a script called ``~/sage_nb.sh`` containing the following lines, and fill in the correct paths for your desired starting directory and ``SAGE_ROOT``


.. code-block:: bash

    #!/bin/bash
    # Switch to desired windows directory
    cd /mnt/c/path/to/desired/starting/directory
    # Start the Jupyter notebook
    SAGE_ROOT/sage --notebook
    # Alternatively you can run JupyterLab - delete the line above, and uncomment the line below
    #SAGE_ROOT/sage --notebook jupyterlab

Make it executable:

.. code-block:: console

    $ chmod ug+x ~/sage_nb.sh

Run it to test:

.. code-block:: console

    $ cd ~
    $ ./sage_nb.sh

The Jupyter(Lab) server should start in the terminal window, and you windows browser should open a page showing the Jupyter or JupyterLab starting page, at the directory you specified.

Create a shortcut
"""""""""""""""""

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

------------------------------------------------------------------------

For further reading you can have a look at the other documents in the
SageMath documentation at http://doc.sagemath.org/.


.. _sec-launching-system-jupyter:

Setting up SageMath as a Jupyter kernel in an existing Jupyter notebook or JupyterLab installation
--------------------------------------------------------------------------------------------------

By default, SageMath installs itself as a Jupyter kernel in the same
environment as the SageMath installation. This is the most convenient way to
use SageMath in a Jupyter notebook. If your notebook or JupyterLab server also
runs from that environment, no extra setup is needed.

To check if the Sage kernel is available, start the Jupyter notebook or
JupyterLab server that you plan to use and look for the kernel named
``SageMath <x.y.z>`` in the list of available kernels. Alternatively, run the
``jupyter`` command from that same Jupyter installation:

.. code-block:: console

    $ jupyter kernelspec list
    Available kernels:
      python3     <path to env>/share/jupyter/kernels/python3
      sagemath    <path to env>/share/jupyter/kernels/sagemath

In case the Sage kernel is not listed in that Jupyter installation, you can
install it with the ``sage --jupyter-kernel`` command described below.

You may already have a global installation of Jupyter. For added
convenience, it is possible to make your installation of SageMath available
in that Jupyter installation, adding it to the list of available kernels
that can be selected in the notebook or JupyterLab interface.

If that Jupyter installation lives in the *same* Python environment as
SageMath, the kernel is already available. Otherwise, assuming that the command
``sage`` starts the SageMath installation that you want to use as a kernel, run

.. code-block:: console

    $ sage --jupyter-kernel install --user

This installs a kernel spec into your per-user Jupyter directory. By default the
spec is *portable*: it records the absolute path to Sage's Python interpreter
and adds Sage's ``bin`` directories to ``PATH``, so it works even when the
Jupyter server runs in a different environment, for example a system-wide
JupyterLab or a JupyterHub deployment. (Pass ``--no-portable`` to instead write
a spec that uses a bare ``python3``, which only works when the Jupyter server
itself runs inside Sage's environment.)

The portable kernel spec is portable between Jupyter environments, but it is
not relocatable: it contains absolute paths to this SageMath installation. If
you move SageMath, switch to a different SageMath installation, or upgrade in a
way that changes Sage's Python environment, rerun the command above.

Only ``kernel.json`` and the kernel icons are installed, and the SageMath
documentation is symlinked rather than copied. This avoids copying gigabytes of
documentation into your Jupyter configuration directory.

Run

.. code-block:: console

    $ sage --jupyter-kernel install --help

to see the available options, such as ``--sys-prefix`` or ``--prefix`` to
install into a different location instead of the per-user directory. For a
shared JupyterHub or another managed Jupyter installation, an administrator may
prefer ``--prefix`` with the Jupyter data prefix used by that installation.

The command above registers the SageMath kernel. Some optional frontend
features may also need packages or assets in the environment from which the
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
   is needed, as SageMath's custom build of the Javascript package
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

Using Jupyter notebook through Visual Studio Code (VS Code) in WSL
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

If you have installed Sage on Windows using Windows Subsystem for
Linux (WSL), it is convenient to use Visual Studio Code (VS Code)
to interact with Sage.

Here are steps to use SageMath in a Jupyter notebook in VS Code:

* Install and run `VS Code <https://code.visualstudio.com/download>`_
  in Windows.

* Click the "Extension" icon on the left (or press :kbd:`Ctrl` +
  :kbd:`Shift` + :kbd:`X`) to open a list of extensions. Install the
  "WSL" and "Jupyter" extensions.

* In the command palette (:kbd:`Ctrl` + :kbd:`Shift` + :kbd:`P`),
  enter "Remote-WSL: New Window", and hit :kbd:`Enter`.

* In the command palette, enter "Create: New Jupyter Notebook", and
  hit :kbd:`Enter`.

* Click "Select Kernel" on the right (or press :kbd:`Ctrl` +
  :kbd:`Alt` + :kbd:`Enter`), select SageMath, and hit :kbd:`Enter`.
