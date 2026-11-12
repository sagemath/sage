"""Regression tests for generic overloads in documentation cross-references."""

import sys
from importlib.util import find_spec
from types import ModuleType

import pytest

if find_spec('sage_docbuild') is None or find_spec('sphinx') is None:
    pytest.skip(
        'sage_docbuild and Sphinx are optional build dependencies',
        allow_module_level=True,
    )

from sage_docbuild.ext.crossrefs import (
    _is_type_parameter,
    _type_parameter_names,
    _type_parameters_of_imported_module,
)


@pytest.mark.parametrize('descriptor', ['staticmethod', 'classmethod'])
@pytest.mark.parametrize('overload_outermost', [False, True])
def test_overloaded_descriptor_parameters(monkeypatch, descriptor, overload_outermost):
    """Both decorator orders retain the overload's local type parameter."""
    name = f'_sage_docbuild_{descriptor}_{overload_outermost}'
    module = ModuleType(name)
    monkeypatch.setitem(sys.modules, name, module)
    outer, inner = ('overload', descriptor) if overload_outermost else (descriptor, 'overload')
    receiver = 'cls, ' if descriptor == 'classmethod' else ''
    # MissingBound must remain unevaluated when collecting parameter names.
    exec(f'''
from typing import overload

property_reads = 0

class Example:
    @{outer}
    @{inner}
    def method[DescriptorT: MissingBound]({receiver}value: DescriptorT) -> DescriptorT: ...

    @{descriptor}
    def method({receiver}value):
        return value

    @property
    def unread(self):
        global property_reads
        property_reads += 1
        raise AssertionError("A documentation lookup evaluated a property")
''', vars(module))
    try:
        assert _type_parameter_names(name) == {'DescriptorT'}
        assert module.property_reads == 0
        node = {'py:module': name, 'py:class': 'Example'}
        assert _is_type_parameter(node, 'DescriptorT')
        assert _is_type_parameter(node, 'typing.DescriptorT')
        for target in ('UnknownType', 'unrelated.DescriptorT', 'typing.Any'):
            assert not _is_type_parameter(node, target)
    finally:
        _type_parameters_of_imported_module.cache_clear()
