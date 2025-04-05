from collections.abc import Sequence
from typing import Any

def call_registered_function(serial: int, nargs: int, args: list, hold: bool, allow_numeric_result: bool, result_parent: Any) -> Any: ...
def find_registered_function(name: str, nargs: int) -> int: ...
def register_or_update_function(self: Any, name: str, latex_name: str, nargs: int, evalf_params_first: bool, update: bool) -> int: ...
def get_sfunction_from_serial(serial: int) -> Any: ...
def get_sfunction_from_hash(myhash: int) -> Any: ...
