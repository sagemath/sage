"""
Crystal features that are imported by default in the interpreter namespace
"""
from __future__ import absolute_import

from . import catalog as crystals

from sage.misc.lazy_import import lazy_import

lazy_import('sage.combinat.crystals.letters',
            'CrystalOfLetters',
            deprecation=(15882, "this is being removed from the global namespace. Use crystals.Letters instead"))

lazy_import('sage.combinat.crystals.fast_crystals',
            'FastCrystal',
            deprecation=(15882, "this is being removed from the global namespace. Use crystals.FastRankTwo instead"))

lazy_import('sage.combinat.crystals.highest_weight_crystals',
            'HighestWeightCrystal',
            deprecation=(15882, "this is being removed from the global namespace. Use crystals.HighestWeight instead"))

lazy_import('sage.combinat.crystals.kyoto_path_model',
            'KyotoPathModel',
            deprecation=(15882, "this is being removed from the global namespace. Use crystals.KyotoPathModel instead"))

lazy_import('sage.combinat.crystals.direct_sum',
            'DirectSumOfCrystals',
            deprecation=(15882, "this is being removed from the global namespace. Use crystals.DirectSum instead"))

lazy_import('sage.combinat.crystals.tensor_product',
            ['CrystalOfTableaux', 'TensorProductOfCrystals'],
            deprecation=(15882, "this is being removed from the global namespace. Use crystals.<tab> instead"))

lazy_import('sage.combinat.crystals.spins',
            ['CrystalOfSpins', 'CrystalOfSpinsPlus', 'CrystalOfSpinsMinus'],
            deprecation=(15882, "this is being removed from the global namespace. Use crystals.<tab> instead"))

lazy_import('sage.combinat.crystals.affine',
            ['AffineCrystalFromClassical', 'AffineCrystalFromClassicalAndPromotion'],
            deprecation=(15882, "this is being removed from the global namespace. Use crystals.<tab> instead"))

lazy_import('sage.combinat.crystals.elementary_crystals',
            ['TCrystal', 'RCrystal', 'ElementaryCrystal', 'ComponentCrystal'],
            deprecation=(15882, "this is being removed from the global namespace. Use crystals.elementary.<tab> instead"))

lazy_import('sage.combinat.crystals.kirillov_reshetikhin',
            'KirillovReshetikhinCrystal',
            deprecation=(15882, "this is being removed from the global namespace. Use crystals.KirillovResetikhin instead"))

lazy_import('sage.combinat.crystals.littelmann_path',
            ['CrystalOfLSPaths', 'CrystalOfProjectedLevelZeroLSPaths'],
            deprecation=(15882, "this is being removed from the global namespace. Use crystals.<tab> instead"))

lazy_import('sage.combinat.crystals.generalized_young_walls',
            ['InfinityCrystalOfGeneralizedYoungWalls', 'CrystalOfGeneralizedYoungWalls'],
            deprecation=(15882, "this is being removed from the global namespace. Use crystals.<tab> instead"))

lazy_import('sage.combinat.crystals.monomial_crystals',
            ['InfinityCrystalOfNakajimaMonomials', 'CrystalOfNakajimaMonomials'],
            deprecation=(15882, "this is being removed from the global namespace. Use crystals.<tab> instead"))

lazy_import('sage.combinat.crystals.infinity_crystals',
            'InfinityCrystalOfTableaux',
            deprecation=(15882, "this is being removed from the global namespace. Use crystals.infinity.Tableaux instead"))

