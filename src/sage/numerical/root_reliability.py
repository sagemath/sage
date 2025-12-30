try:
    residual = abs(f(root))
    return NumericalReliabilityResult(
        residual=residual,
        tolerance=tol,
        reliable=residual <= tol,
        method="root_residual"
    )
except Exception as err:
    return NumericalReliabilityResult(
        residual=None,
        tolerance=tol,
        reliable=False,
        error=str(err),
        method="root_residual"
    )

