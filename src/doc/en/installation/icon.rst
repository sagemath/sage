KDE Desktop icon
================

These instructions will help you make a KDE desktop icon which
starts the Sage notebook. Instructions for a Gnome desktop icon are
probably similar.


#. Create a ``notebook.sage`` file containing only the line

   ::

      notebook(openviewer=True)

#. In your Desktop subdirectory, create a file
   ``Sage-notebook.desktop`` containing the lines

   ::

       [Desktop Entry]
       Comment=
       Comment[de]=
       Encoding=UTF-8
       Exec=/usr/local/bin/sage /home/martin/notebook.sage
       GenericName=
       GenericName[de]=
       Icon=
       MimeType=
       Name=Sage
       Name[de]=Sage
       Path=$HOME
       StartupNotify=true
       Terminal=false
       TerminalOptions=
       Type=Application
       X-DCOP-ServiceType=
       X-KDE-SubstituteUID=false
       X-KDE-Username=

   You will have to edit the ``Exec=`` line to point to your ``sage``
   script and your ``notebook.sage`` file.

#. Right click on the Sage notebook desktop icon and click on
   ``Properties``, then ``Application``, then ``Advanced Options``,
   then ``Run in Terminal``. If you want to title the xwindow
   terminal, add in the terminal option box ``-T "sage notebook"``.


To quit the Sage notebook, first enter ``Ctrl-c`` in the xwindow
terminal running Sage, then enter ``Ctrl-d`` to quit Sage in the
terminal, and finally close the browser (or browser tab) which was
displaying the Sage notebook server.

For a picture for your icon, check out the Sage art at http://wiki.sagemath.org/art.
