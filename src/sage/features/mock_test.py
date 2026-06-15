import pytest

from sage.features import FeatureNotPresentError
from sage.features.mock import MockInstance, mock


def test_mock_creation():
    m = mock('module', 'foo')
    assert isinstance(m, MockInstance)
    assert m._name == 'foo'
    assert m._module == 'module'

def test_mock_repr_returns_name():
    m = mock('module', 'foo')
    assert repr(m) == "<Mock name='foo' module='module'>"

def test_mock_call_raises_exception():
    m = mock('module', 'foo')
    with pytest.raises(FeatureNotPresentError):
        m()

def test_mock_call_function_raises_exception():
    m = mock('module', 'foo')
    with pytest.raises(FeatureNotPresentError):
        m.function() # type: ignore

def test_mock_get_attribute_returns_another_mock():
    m = mock('module', 'foo')
    attribute = m.bar # type: ignore
    assert isinstance(attribute, MockInstance)
    assert attribute._name == 'foo.bar'
    assert attribute._module == 'module'

def test_can_inherit_from_mock_without_raising_exception():
    m = mock('module', 'foo')
    class MyMock(m): # type: ignore
        pass
    with pytest.raises(FeatureNotPresentError):
        MyMock().bar()

def test_can_inherit_from_mock_with_other_bases():
    m = mock('module', 'foo')
    class Bar:
        pass
    class MyMock(m, Bar): # type: ignore
        pass
    with pytest.raises(FeatureNotPresentError):
        MyMock().bar()

def test_can_inherit_from_mock_with_other_metaclass():
    m = mock('module', 'foo')
    class MyMeta(type):
        def __new__(cls, name, bases, dct):
            return super().__new__(cls, name, bases, dct)
    class Bar(metaclass=MyMeta):
        pass
    class MyMock(m, Bar): # type: ignore
        pass
    with pytest.raises(FeatureNotPresentError):
        MyMock().bar()

# def test_can_inherit_from_mock_and_algebra_element():
#     m = mock('module', 'foo')
#     from sage.structure.element import AlgebraElement
#     class MyMock(m, AlgebraElement):
#         pass
#     with pytest.raises(FeatureNotPresentError):
#         MyMock()
