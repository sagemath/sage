This is Sage -- open source mathematical software.

         http://www.sagemath.org

**
  These binaries are only for the OS X version that is indicated by the package .dmg name.
  They generally will not work on any other OS X version (unless explicitly stated otherwise).
**

There are two ways in which Sage can be distributed.  One is as a
"regular" OS X application named something like Sage-VERSION.  If you
see such an application, skip to the section about Sage.app.  If
instead you see a folder called "sage", proceed as follows.

1) Download the dmg somewhere and double click on it.

2) Drag the sage folder somewhere, e.g., /Applications

 ** WARNING ** If you get an error copying the folder do the following:
    Do not drag the folder out of the dmg image. Use the
    shell (via Terminal) and do a
     "cp -R -P /Volumes/sage-2.9.2-OSX10.4-intel-i386-Darwin/sage ."
    from the location where you want to install to. Adjust the name of the
    Volume as needed.

3) Use finder to visit the sage folder you just copied it and double click on the "sage" icon.

4) Select to run it with "Terminal":
     Choose Applications, then select "All Applications" in the
     "Enable:" drop down.  Change the "Applications" drop down
     to "Utilities".  On the left, scroll and select "Terminal".
     Click "Open", then in the next dialog select "Update".

5) Sage should pop up in a window.

6) For the graphical notebook, type
     notebook()
   You might have to open Firefox or Safari (your choice)
   to the URL
       http://localhost:8000
   to use Sage on your computer.

7) Email
      http://groups.google.com/group/sage-support
   with any questions.


========
Sage.app
========

Sage.app is a fairly simple application which starts the Sage server
when launched and stops it when quitting.  It also contains a
primitive browser which can be used to view the Sage notebook if
desired.  Sage.app can also be used to launch terminal sessions for
many of your favorite Sage commands.  Note that holding the option key
changes some of the commands that can be run, and holding the command
key will allow you to edit the command before running them.

Any feedback is welcome: http://groups.google.com/group/sage-devel

------------
Installation
------------

Simply copy the application to your hard drive as you would any other
application.

You may also wish to add the sage executable to your PATH for use in a
normal terminal session.  This can be done from the Preferences window
and the Add To PATH button.  It will ask whether to add to
/usr/local/bin or ~/bin.  To remove this, simply delete the symlink in
whichever location you chose or overwrite it with a new version.

WARNING: When Sage.app is run, it may take a _long_ time for the
server to start up.  Please don't give up.  If you think it has taken
too long, use View Log under the Development menu to see if there are
any errors.

-----
Setup
-----

If you have a copy of Sage.app which did not come bundled with a Sage
distribution (say you built a new version of Sage.app yourself), then
the first time it is run, it will prompt you for a Sage executable.
Please choose the executable that is called sage and is what you would
run from the command line.  You can change this later in the
preferences.  You can also drag a sage folder (SAGE_ROOT) onto
Sage.app to change what sage executable it will use.

-----------
Preferences
-----------

Sage.app can be run as a normal application, or as a menu extra (icon
on the right side of the menu), or both.  This can be set in
Preferences where you can set whether to use Sage.app as a browser,
whether to always start the server, and whether you wish to always
edit terminal session arguments before they are run.

You may also set default parameters for terminal sessions.  For
example to send -b to gap (to inhibit the initial banner), type gap as
the session (lowercase as you see it in the menus--use "sage" (without
quotes) for sage sessions, and "notebook" for starting the notebook
server), -b as the argument and hit Apply.  Unfortunately, due to a
bug in the Cocoa framework, choosing the session with the mouse does
not correctly update the Default Arguments field.  Instead please
select a previous session with the keyboard.  The same applies to the
Terminal Emulator selection.

You may select a method of running terminal sessions.  If you do not
see your favorite terminal emulator here, please let know what it is,
so we can add it.  You can also create/edit the applescript directly
for the ultimate in flexibility.  Place %@ where you wish the command
to be inserted (but only do so once).

------------
SAGE_BROWSER
------------

WARNING: setting SAGE_BROWSER in .bashrc will overwrite the value
which allows Sage.app to work.  Therefore please use something like

SAGE_BROWSER=${SAGE_BROWSER-/path/to/script/to/open/your/browser}

(or nothing at all) in your .bashrc if you wish to use Sage.app as the
Sage notebook browser.

-----------------
OS X 10.4 (Tiger)
-----------------

On Tiger (10.4) the method of making Sage.app appear in the Dock does
not work, hence the preference to show in the Dock is disabled.  If
you wish the application to appear in the Dock, you may run the
following in a shell (changing the path to Sage.app appropriately)

defaults write /Applications/Sage.app/Contents/Info LSUIElement 0

To go back to Dockless mode:

defaults write /Applications/Sage.app/Contents/Info LSUIElement 1

You will of course have to restart Sage.app

--------
Building
--------

To build Sage.app yourself, you just have to run "make" in the
SAGE_ROOT/src/mac-app directory. This app will only work locally
because the contained Sage is not relocatable. To build a binary
distribution that can be installed to a different directory you must
use https://github.com/sagemath/binary-pkg

If you wish to make changes, or create a version which
does not contain a distribution, then open
src/mac-app/Sage.xcodeproj in Xcode.  If you are
building on OS X 10.4, then you will need to change the SDK by opening
Project Settings, and in Cross-Develop in General set it to 10.4
(universal).
