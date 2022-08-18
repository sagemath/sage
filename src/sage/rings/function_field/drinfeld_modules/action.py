from sage.categories.action import Action

from sage.rings.function_field.drinfeld_modules.drinfeld_module import DrinfeldModule


class DrinfeldModuleAction(Action):
    def __init__(self, finite_drinfeld_module):
        # Verifications
        if not isinstance(finite_drinfeld_module, DrinfeldModule):
            raise TypeError('input must be a DrinfeldModule')
        # Work
        self.__finite_drinfeld_module = finite_drinfeld_module
        super().__init__(finite_drinfeld_module.polring(),
                finite_drinfeld_module.ore_polring().base_ring())

    ###########
    # Methods #
    ###########

    def finite_drinfeld_module(self):
        return self.__finite_drinfeld_module

    ##########################
    # Special Sage functions #
    ##########################

    def _act_(self, g, x):
        return self.finite_drinfeld_module()(g)(x)

    def _latex_(self):
        return f'\\text{{Action{{ }}on{{ }}}}' \
                f'{latex(self.extension())}\\text{{{{ }}' \
                f'induced{{ }}by{{ }}}}{self.finite_drinfeld_module()}'

    def _repr_(self):
        return f'Action on {self.domain()} induced by ' \
                f'{self.finite_drinfeld_module()}'
