from cpython.object import PyTypeObject, newfunc, destructor

def hook_tp_functions_type(tp: PyTypeObject, tp_new: newfunc, tp_dealloc: destructor, useGC: bool) -> None:
    ...

def hook_tp_functions(global_dummy: object, tp_new: newfunc, tp_dealloc: destructor, useGC: bool) -> None:
    ...
