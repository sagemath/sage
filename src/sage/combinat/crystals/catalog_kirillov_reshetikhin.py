# sage_setup: distribution = sagemath-combinat
"""
Catalog Of Crystal Models For Kirillov-Reshetikhin Crystals

We currently have the following models:

* :func:`KashiwaraNakashimaTableaux
  <sage.combinat.crystals.kirillov_reshetikhin.KashiwaraNakashimaTableaux>`
* :class:`~sage.combinat.rigged_configurations.kr_tableaux.KirillovReshetikhinTableaux`
* :func:`LSPaths <sage.combinat.crystals.kirillov_reshetikhin.KirillovReshetikhinCrystalFromLSPaths>`
* :func:`RiggedConfigurations
  <sage.combinat.rigged_configurations.rigged_configurations.KirillovReshetikhinCrystal>`
"""

from .kirillov_reshetikhin import KashiwaraNakashimaTableaux
from .kirillov_reshetikhin import KirillovReshetikhinCrystalFromLSPaths as LSPaths
from sage.combinat.rigged_configurations.kr_tableaux import KirillovReshetikhinTableaux
from sage.combinat.rigged_configurations.rigged_configurations import KirillovReshetikhinCrystal as RiggedConfigurations
