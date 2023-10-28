import importlib

from .builddir.build_info import CONDA_PREFIX

conf = importlib.import_module('.builddir.pkgs.sage-conf_meson._conf', 'sage_conf')
conf_options = [x for x in conf.__dict__ if not x.startswith("_")]
globals().update({k: getattr(conf, k) for k in conf_options})

SAGE_LOCAL = CONDA_PREFIX
