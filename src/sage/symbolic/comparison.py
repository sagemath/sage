from sage.misc.lazy_import import lazy_import
lazy_import('sage.symbolic.expression',
            ['print_order', '_print_key', 'print_sorted', '_math_key',
             'math_sorted', 'mixed_order', '_mixed_key', 'mixed_sorted'],
            deprecation=32386)
