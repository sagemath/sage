# pylint: disable=C0115:missing-class-docstring,C0116:missing-function-docstring,W0212:protected-access
# pyright: reportPrivateUsage=false

from typing import List, Tuple
import pytest
from sage.tensor.modules.free_module_tensor import IndexCharacterNormalized
from sage.tensor.modules.tensor_with_indices import TensorWithIndices

ParseResult = List[Tuple[str, IndexCharacterNormalized]]


class TestTensorWithIndices:
    @pytest.mark.parametrize(
        "input_string, expected_output",
        [
            ("^a", [("a", "UP")]),
            ("^a^b", [("a", "UP"), ("b", "UP")]),
            ("^a^b^c", [("a", "UP"), ("b", "UP"), ("c", "UP")]),
        ],
    )
    def test_parse_indices_handles_up_indices(
        self, input_string: str, expected_output: ParseResult
    ):
        assert TensorWithIndices._parse_indices(input_string) == (
            expected_output,
            [],
            [],
            [],
        )

    @pytest.mark.parametrize(
        "input_string, expected_output",
        [
            ("_a", [("a", "DOWN")]),
            ("_a_b", [("a", "DOWN"), ("b", "DOWN")]),
            ("_a_b_c", [("a", "DOWN"), ("b", "DOWN"), ("c", "DOWN")]),
        ],
    )
    def test_parse_indices_handles_down_indices(
        self, input_string: str, expected_output: ParseResult
    ):
        assert TensorWithIndices._parse_indices(input_string) == (
            expected_output,
            [],
            [],
            [],
        )

    @pytest.mark.parametrize(
        "input_string, expected_output",
        [
            ("a", [("a", "UP")]),
            ("a_b", [("a", "UP"), ("b", "DOWN")]),
            ("a_b_c", [("a", "UP"), ("b", "DOWN"), ("c", "DOWN")]),
        ],
    )
    def test_parse_indicies_assumes_up_by_default(
        self, input_string: str, expected_output: ParseResult
    ):
        assert TensorWithIndices._parse_indices(input_string) == (
            expected_output,
            [],
            [],
            [],
        )

    @pytest.mark.parametrize(
        "input_string, expected_output",
        [
            ("^a", [("a", "UP")]),
            ("_a^b", [("a", "DOWN"), ("b", "UP")]),
            ("^a_b^c", [("a", "UP"), ("b", "DOWN"), ("c", "UP")]),
            ("^a_b_c", [("a", "UP"), ("b", "DOWN"), ("c", "DOWN")]),
            ("_a^b_c", [("a", "DOWN"), ("b", "UP"), ("c", "DOWN")]),
            ("_a_b^c", [("a", "DOWN"), ("b", "DOWN"), ("c", "UP")]),
        ],
    )
    def test_parse_indices_handles_mixed_indices(
        self, input_string: str, expected_output: ParseResult
    ):
        assert TensorWithIndices._parse_indices(input_string) == (
            expected_output,
            [],
            [],
            [],
        )

    @pytest.mark.parametrize(
        "input_string, expected_output",
        [
            ("^ab", [("a", "UP"), ("b", "UP")]),
            ("_ab", [("a", "DOWN"), ("b", "DOWN")]),
            ("^ab_c", [("a", "UP"), ("b", "UP"), ("c", "DOWN")]),
            ("^ab_cd", [("a", "UP"), ("b", "UP"), ("c", "DOWN"), ("d", "DOWN")]),
            ("_ab^cd", [("a", "DOWN"), ("b", "DOWN"), ("c", "UP"), ("d", "UP")]),
        ],
    )
    def test_parse_indices_assumes_following_indicies_have_same_character(
        self, input_string: str, expected_output: ParseResult
    ):
        assert TensorWithIndices._parse_indices(input_string) == (
            expected_output,
            [],
            [],
            [],
        )

    @pytest.mark.parametrize(
        "input_string, expected_output",
        [
            ("{a}", [("a", "UP")]),
            ("{a}{b}", [("a", "UP"), ("b", "UP")]),
            ("^{ab}_c", [("a", "UP"), ("b", "UP"), ("c", "DOWN")]),
            ("^a_{bc}", [("a", "UP"), ("b", "DOWN"), ("c", "DOWN")]),
        ],
    )
    def test_parse_indices_handles_latex_braces(
        self, input_string: str, expected_output: ParseResult
    ):
        assert TensorWithIndices._parse_indices(input_string) == (
            expected_output,
            [],
            [],
            [],
        )

    @pytest.mark.parametrize(
        "input_string", ["(..", "[..", "(..)[", "(..)]", ")..", "]..", "(..))", "()"]
    )
    def test_parse_indices_raises_error_on_unbalanced_parenthesis(
        self, input_string: str
    ):
        with pytest.raises(ValueError):
            TensorWithIndices._parse_indices(input_string)

    @pytest.mark.parametrize(
        "input_string",
        [
            "[(..)]",
            "([..])",
        ],
    )
    def test_parse_indices_raises_error_on_nested_parenthesis(
        self, input_string: str
    ) -> None:
        with pytest.raises(ValueError):
            TensorWithIndices._parse_indices(input_string)

    @pytest.mark.parametrize(
        "input_string",
        [
            "^17",  # numbers are not allowed
            "^;",  # special characters are not allowed
            "^\u00ae",  # special characters are not allowed
            "^\u25e2",  # special characters are not allowed
        ],
    )
    def test_parse_raises_error_on_invalid_input(self, input_string: str):
        with pytest.raises(ValueError):
            TensorWithIndices._parse_indices(input_string)

    @pytest.mark.parametrize("input_string", ["^a^a", "_a_a", "^a_b^a", "aa"])
    def test_parse_indices_raises_error_on_invalid_contraction(self, input_string: str):
        with pytest.raises(ValueError):
            TensorWithIndices._parse_indices(input_string)

    @pytest.mark.parametrize(
        "input_string, expected_output, expected_symmetries",
        [
            ("^(ab)", [("a", "UP"), ("b", "UP")], ([(0, 1)], [])),
            ("^(ab)_c", [("a", "UP"), ("b", "UP"), ("c", "DOWN")], ([(0, 1)], [])),
            (
                "^(ab)_(cd)",
                [("a", "UP"), ("b", "UP"), ("c", "DOWN"), ("d", "DOWN")],
                ([(0, 1), (2, 3)], []),
            ),
            (
                "^[ab]_[cd]",
                [("a", "UP"), ("b", "UP"), ("c", "DOWN"), ("d", "DOWN")],
                ([], [(0, 1), (2, 3)]),
            ),
        ],
    )
    def test_parse_indices_handles_symmetries(
        self,
        input_string: str,
        expected_output: ParseResult,
        expected_symmetries: Tuple[List[Tuple[int, int]], List[Tuple[int, int]]],
    ):
        assert TensorWithIndices._parse_indices(input_string) == (
            expected_output,
            expected_symmetries[0],
            expected_symmetries[1],
            [],
        )

    def test_parse_indices_raises_error_if_contraction_not_allowed_but_present(self):
        with pytest.raises(ValueError):
            TensorWithIndices._parse_indices("^a_a", allow_contraction=False)

    def test_parse_indices_raises_error_if_symmetries_not_allowed_but_present(self):
        with pytest.raises(ValueError):
            TensorWithIndices._parse_indices("^[ab]_[cd]", allow_symmetries=False)

    def test_parse_indices_raises_error_if_index_configuration_not_as_specified(self):
        with pytest.raises(IndexError):
            TensorWithIndices._parse_indices(
                "^(ab)_c", index_configuration=("UP", "DOWN", "DOWN")
            )
