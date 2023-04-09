from configparser import ConfigParser
import os

config = ConfigParser()

config.add_section('directories')
config.set('directories', 'basedirlist', os.environ['SAGE_LOCAL'])

config.add_section('libs')
config.set('libs', 'system_freetype', 'True')
config.set('libs', 'system_qhull', 'True')
# lto is problematic if we mix libraries from the OS with our own libraries,
# which are not necessarily compiled with the same gcc version
# https://github.com/sagemath/sage/issues/27754
config.set('libs', 'enable_lto', 'False')

#####################################################################
# Sage code -- all this code just sets the graphical_backend variable.
# If True, that means we try to build GUI's; otherwise, we definitely
# will not even try, even if we could.  See trac #5301.
#####################################################################

print("NOTE: Set SAGE_MATPLOTLIB_GUI to anything but 'no' to try to build the Matplotlib GUI.")

graphical_backend = 'False'
if os.environ.get('SAGE_MATPLOTLIB_GUI', 'no').lower() != 'no':
    graphical_backend = 'auto'

if graphical_backend == 'auto':
    print("Building graphical backends.  WARNING: This may causing some annoying and confusing behavior")
    print("when using Sage + pylab, at least on OS X.")
else:
    print("Not building any matplotlib graphical backends.")

config.add_section('gui_support')
for backend in ('gtk', 'gtkagg', 'tkagg', 'wxagg', 'macosx', 'windowing'):
    config.set('gui_support', backend, graphical_backend)

with open('src/mplsetup.cfg', 'w') as configfile:
    config.write(configfile)
