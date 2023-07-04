from sage.misc.lazy_import import lazy_import
lazy_import('sage.tensor.modules.finite_rank_free_module',
            'FiniteRankFreeModule')
# NB: in Sage 8.8.beta2, the lazy import of FiniteRankFreeModule is necessary
# to avoid some import order issue when Chart is imported in
# free_module_tensor, see comments 12 to 18 in :trac:`27655`.
