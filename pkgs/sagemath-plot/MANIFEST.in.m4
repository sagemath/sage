dnl MANIFEST.in is generated from this file by SAGE_ROOT/bootstrap via m4.
prune sage

include VERSION.txt

exclude *.m4
include requirements.txt

global-include all__sagemath_plot.py

graft sage/plot

include sage/interfaces/jmoldata.p*
include sage/interfaces/povray.p*
include sage/interfaces/tachyon.p*
include sage/interfaces/gnuplot.p*

include sage/ext_data/graphs/graph_plot_js.html
graft sage/ext_data/threejs

global-exclude __pycache__
global-exclude *.py[co]
global-exclude *.bak
global-exclude *.so
global-exclude *~
prune .tox
prune build
prune dist
