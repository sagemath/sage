# All the smooth reflexive polytopes in dimension 3

polytopes_3d = {
    'F.3D.0000': [[1, 1, 0], [1, 0, 0], [0, 1, -1], [1, -2, 1], [1, 1, 1], [0, 1, 1], [0, -3, 1], [-4, 1, -3]],
    'F.3D.0001': [[1, 1, -1], [1, 1, 1], [-4, 1, 1], [1, -4, 1], [0, 1, -1], [1, 0, -1]],
    'F.3D.0002': [[0, 1, -1], [0, 0, -1], [1, 1, 1], [1, -1, 1], [1, 0, 0], [1, 1, 0], [-1, 1, 1], [-1, -3, 1], [-1, -1, -1], [-1, 1, -1]],
    'F.3D.0003': [[0, -2, 1], [0, 1, 1], [1, 1, 0], [1, 0, 0], [1, -1, 1], [1, 1, 1], [0, 0, -1], [0, 1, -1], [-1, 1, -1], [-1, -1, -1], [-1, -2, 0], [-1, 1, 0]],
    'F.3D.0004': [[0, 1, -1], [0, 0, -1], [1, 1, 1], [1, -1, 1], [1, 0, 0], [1, 1, 0], [0, 1, 1], [0, -2, 1], [-2, -2, -1], [-2, 1, -1]],
    'F.3D.0005': [[1, 0, -1], [1, -2, 1], [1, 1, 1], [1, 1, -1], [-2, 1, 1], [0, 1, -1], [-2, -2, 1], [0, 0, -1]],
    'F.3D.0006': [[1, 1, 0], [1, 0, -1], [1, 0, 1], [1, 1, 1], [-1, 1, 0], [-1, 1, 1], [-1, -2, 1], [-1, -2, -3]],
    'F.3D.0007': [[0, -2, 1], [0, 1, 1], [1, 1, 1], [1, -1, 1], [1, 0, 0], [1, 1, 0], [-1, 0, -2], [-1, 1, -2], [-1, 1, 0], [-1, -2, 0]],
    'F.3D.0008': [[1, -1, 0], [1, 0, -1], [1, 1, -1], [1, 1, 1], [1, -1, 1], [-1, 1, 1], [-1, 1, -1], [-1, -1, 1], [-1, 0, -1], [-1, -1, 0]],
    'F.3D.0009': [[1, -1, 0], [1, -1, -1], [1, 0, -1], [1, 1, 1], [1, 1, 0], [1, 0, 1], [-1, 1, 0], [-1, 1, 1], [-1, 0, -1], [-1, 0, 1], [-1, -1, -1], [-1, -1, 0]],
    'F.3D.0010': [[1, 1, -2], [-2, 1, 1], [1, 1, 1], [1, -1, 1], [1, -1, -1], [-1, -1, 1], [-2, 0, 1], [1, 0, -2]],
    'F.3D.0011': [[1, 0, -1], [1, -2, 1], [1, 1, 1], [1, 1, -1], [-1, 1, 1], [-1, 1, -1], [-1, -2, 1], [-1, 0, -1]],
    'F.3D.0012': [[1, 0, -1], [1, -2, 1], [1, 1, 1], [1, 1, -1], [0, 1, 1], [-2, 1, -1], [0, -2, 1], [-2, 0, -1]],
    'F.3D.0013': [[1, 0, 1], [1, 1, 1], [1, 0, -2], [1, 1, -2], [-2, 1, 1], [-2, -3, 1]],
    'F.3D.0014': [[1, 1, -1], [1, 1, 1], [-3, 1, 1], [1, -3, 1], [-1, 1, -1], [1, -1, -1]],
    'F.3D.0015': [[1, -1, -1], [1, -1, 1], [1, 1, 1], [1, 1, -1], [-1, 1, 1], [-1, 1, -1], [-1, -1, 1], [-1, -1, -1]],
    'F.3D.0016': [[1, -1, 1], [1, 1, 1], [-2, 1, 1], [1, 1, -2], [-2, -1, 1], [1, -1, -2]],
    'F.3D.0017': [[1, 1, 1], [-3, 1, 1], [1, -3, 1], [1, 1, -3]],
}

# The table in Batyrev's paper, the values are (c_1^3, c_1 c_2, b_2, h^0, rank of the automorphism group), which can provide a good identification of the variety (only for smooth toric Fano 3-folds)
BATYREV_3FOLD_LOOKUP = {
    # \rho=1
    "P3": (3, 64, 1, 35, 15),

    # \rho=2
    "B1": (3, 62, 2, 34, 15),  # P_{P^2}(O \otimes O(2))
    "B2": (3, 56, 2, 31, 12),  # P_{P^2}(O \otimes O(1))
    "B3": (3, 54, 2, 30, 11),  # P_{P^1}(O \otimes O \otimes O(1))
    "B4": (3, 54, 2, 30, 11),  # P^2 × P^1  (same invariants as B3)

    # \rho=3
    "C1": (3, 52, 3, 29, 11),  # P_{P^1×P^1}(O \otimes O(1,1))
    "C2": (3, 50, 3, 28, 10),  # P_{S1}(O \otimes O(l)), l^2=1 on S1
    "C3": (3, 48, 3, 27, 9),   # P^1 × P^1 × P^1
    "C4": (3, 48, 3, 27, 9),   # S1 × P^1 (same invariants as C3)
    "C5": (3, 44, 3, 25, 7),   # P_{P^1×P^1}(O \otimes O(1,−1))

    # \rho=3 (blow-ups of P^3)
    "D1": (3, 50, 3, 28, 10),  # blow-up of P^3 (same invariants as C2)
    "D2": (3, 46, 3, 26, 8),   # blow-up of P^3

    # \rho=4
    "E1": (3, 46, 4, 26, 9),   # S2-bundle over P^1
    "E2": (3, 44, 4, 25, 8),   # S2-bundle over P^1
    "E3": (3, 42, 4, 24, 7),   # S2 × P^1
    "E4": (3, 40, 4, 23, 6),   # S2-bundle over P^1

    # \rho=5
    "F1": (3, 36, 5, 21, 5),   # S3 × P^1
    "F2": (3, 36, 5, 21, 5),   # S3-bundle over P^1 (same invariants as F1)
}

# AL: the lookup tables for classification of type of Fano 3-folds
# (i.e. projective spaces, products, projective bundles, blow-ups of projective spaces)

# --- 3-folds ---
CLASS_BY_LABEL_3FOLD = {
    # \rho = 1
    "P3": "product",  # treat P^3 as the trivial 'product' case

    # \rho = 2
    "B1": "projective_bundle",  # P_{P^2}(O \otimes O(2))
    "B2": "projective_bundle",  # P_{P^2}(O \otimes O(1))
    "B3": "projective_bundle",  # P_{P^1}(O \otimes O \otimes O(1))
    "B4": "product",            # P^2 × P^1

    # \rho = 3
    "C1": "projective_bundle",  # P_{P^1×P^1}(O \otimes O(1,1))
    "C2": "projective_bundle",  # P_{S1}(O \otimes O(ℓ))
    "C3": "product",            # (P^1)^3
    "C4": "product",            # S1 × P^1
    "C5": "projective_bundle",  # P_{P^1×P^1}(O \otimes O(1,−1))

    # \rho = 3 blow-ups
    "D1": "blow_up",            # blow-up of P^3
    "D2": "blow_up",            # blow-up of P^3

    # \rho = 4
    "E1": "projective_bundle",  # S2-bundle over P^1
    "E2": "projective_bundle",  # S2-bundle over P^1
    "E3": "product",            # S2 × P^1
    "E4": "projective_bundle",  # S2-bundle over P^1

    # \rho = 5
    "F1": "product",            # S3 × P^1
    "F2": "projective_bundle",  # S3-bundle over P^1
}
