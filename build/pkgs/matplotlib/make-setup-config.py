from six.moves import configparser as ConfigParser
import os

config = ConfigParser.SafeConfigParser()

config.add_section('directories')
config.set('directories', 'basedirlist', os.environ['SAGE_LOCAL'])



#####################################################################
# Sage code -- all this code just sets the graphical_backend variable.
# If True, that means we try to build GUI's; otherwise, we definitely
# will not even try, even if we could.  See trac #5301.
#####################################################################

print("NOTE: Set SAGE_MATPLOTLIB_GUI to anything but 'no' to try to build the Matplotlib GUI.")

graphical_backend='False'
if os.environ.get('SAGE_MATPLOTLIB_GUI', 'no').lower() != 'no':
    graphical_backend = 'auto'

if graphical_backend == 'auto':
    print("Building graphical backends.  WARNING: This may causing some annoying and confusing behavior")
    print("when using Sage + pylab, at least on OS X.")
else:
    print("Not building any matplotlib graphical backends.")

config.add_section('gui_support')
for backend in ('gtk', 'gtkagg', 'tkagg', 'wxagg', 'macosx', 'windowing'):
    config.set('gui_support', backend,  graphical_backend)

with open('src/setup.cfg', 'wb') as configfile:
    config.write(configfile)
