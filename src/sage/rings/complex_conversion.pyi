from typing import Any

from sage.structure.element import Element
from sage.categories.map import Map

class CCtoCDF(Map):
    def _call_(self, x: Element) -> Element:
        ...
