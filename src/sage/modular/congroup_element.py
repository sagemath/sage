# sage_setup: distribution = sagemath-schemes
###########################################################
# Re-bindings for unpickling
#
# We want to ensure class
# sage.modular.congroup_element.CongruenceSubgroupElement still exists, so we
# can unpickle safely.
#
###########################################################

from sage.modular.arithgroup.arithgroup_element import ArithmeticSubgroupElement

CongruenceSubgroupElement = ArithmeticSubgroupElement
