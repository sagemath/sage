from sage.categories.action import Action

from sage.rings.function_field.drinfeld_modules.drinfeld_module import DrinfeldModule


class DrinfeldModuleAction(Action):

    def __init__(self, drinfeld_module):
        if not isinstance(drinfeld_module, DrinfeldModule):
            raise TypeError('input must be a DrinfeldModule')
        self._drinfeld_module = drinfeld_module
        super().__init__(drinfeld_module.function_ring(),
                drinfeld_module.base_ring())

    def _act_(self, pol, x):
        if pol not in self._drinfeld_module.function_ring():
            raise TypeError('first input must be in the function ring')
        if x not in self._drinfeld_module.base_ring():
            raise TypeError('second input must be in the base ring')
        return self._drinfeld_module(pol)(x)

    def _latex_(self):
        return f'\\text{{Action{{ }}on{{ }}}}' \
                f'{latex(self._drinfeld_module.base_ring())}\\text{{{{ }}' \
                f'induced{{ }}by{{ }}}}{self._drinfeld_module}'

    def _repr_(self):
        return f'Action on {self._drinfeld_module.base_ring()} induced by ' \
                f'{self._drinfeld_module}'

    def drinfeld_module(self):
        return self._drinfeld_module
