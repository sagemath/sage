from .reliability_certificate import NumericalReliabilityCertificate

def certify_residual(f, root, tol=1e-10):
    try:
        res = abs(f(root))
        passed = res <= tol
        return NumericalReliabilityCertificate(
            method="residual",
            metrics={"residual": res},
            tolerance=tol,
            passed=passed,
        )
    except Exception as e:
        return NumericalReliabilityCertificate(
            method="residual",
            metrics={"error": str(e)},
            tolerance=tol,
            passed=False,
        )

