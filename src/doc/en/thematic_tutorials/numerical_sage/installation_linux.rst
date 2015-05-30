Installing Visualization Tools on Linux
=======================================

This section assumes you are running linux. You may need
administrator rights to complete this section. First try

.. skip

::

    sage: import Tkinter

If this works great if not it means your sage was not compiled with
tcl/tk bindings. There are two possible reasons. If you used a
binary then this will be the case. If you built from source but
don't have the tcl/tk libraries on your system you will have the
same problem. To fix this install the tcl and tk libraries and
development (source and headers) packages for your linux
distribution. Now you need to rebuild Sage's python

.. skip

::

    sage: install_package('python-2.5<version>.spkg')

where you should replace :math:`<\text{version}>` by the version
of python. Type

.. skip

::

    sage: !ls spkg/standard | grep python-2.5

This will give you the name of the python package to use above. In
the case that it gives multiple names, choose the one with the
highest version number. Now again try

.. skip

::

    sage: import Tkinter

If this works we can install vtk, but first you need cmake. Test if
your system has it by typing cmake at the shell, (or !cmake in
sage). If cmake is on your system then you are ready. If not either
install cmake on your system using your distributions tools or do

.. skip

::

    sage: install_package('cmake-2.4.7')

Now we want to compile VTK which does all the hard work. To do this
make sure you have the opengl libraries for you system installed.
This will be something like libgl1-mesa-glx, libgl1-mesa-dev.

.. skip

::

    sage: install_package('vtk-5.0.3.p1')

This will take quite a while to compile, probably 20 min to an hour
or so. Once this is done we install the python wrappers, the next
part takes about 10 seconds,

.. skip

::

    sage: install_package('MayaVi-1.5')
    sage: install_package('scitools++')
    sage: install_package('PyVTK-0.4.74')

Now you're done.

