import pytest
from sage.rings.rational_field import ZZ
from sage.tensor.modules.comp import Components
from sage.tensor.modules.finite_rank_free_module import FiniteRankFreeModule
from sage.tensor.modules.free_module_tensor import FreeModuleTensor


class TestFreeModuleTensor:
    @pytest.fixture
    def module(self):
        module = FiniteRankFreeModule(ZZ, 3, name="M")
        module.basis("e")
        return module

    def test_init_throws_exception_when_tensor_configuration_has_wrong_length(
        self, module: FiniteRankFreeModule
    ):
        with pytest.raises(ValueError):
            FreeModuleTensor(module, (2, 1), config=("^", "^", "^", "^"))

    def test_display_comp_takes_tensor_configuration_into_account(
        self, module: FiniteRankFreeModule
    ):
        tensor = FreeModuleTensor(module, (2, 1), config=("_", "^", "^"))
        tensor[1, 2, 1], tensor[1, 2, 2] = 2, -1
        assert tensor.display_comp() == "X_1^21 = 2 \nX_1^22 = -1 "

    def test_components_transform_according_to_tensor_configuration(
        self, module: FiniteRankFreeModule
    ):
        tensor = FreeModuleTensor(module, (1, 1), config=("_", "^"))
        tensor[0, 1] = -3
        tensor[2, 2] = 2
        assert tensor.display() == "-3 e^0⊗e_1 + 2 e^2⊗e_2"

        a = module.automorphism()
        a[:] = [[0, 0, 1], [1, 0, 0], [0, -1, 0]]
        f = module.basis("e").new_basis(a, "f")
        # f = a*e, f_0=e_1, f_1=-e_2, f_2=e_0
        result = Components(ZZ, f, 2)
        result[:] = [[0, 0, 0], [0, 2, 0], [-3, 0, 0]]
        assert tensor.components(f) == result
        assert tensor.display(f) == "2 f^1⊗f_1 - 3 f^2⊗f_0"
