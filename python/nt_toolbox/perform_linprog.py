"""Linear programming with explicit equality constraints and nonnegative variables."""

import numpy as np
from scipy.optimize import linprog


def perform_linprog(A, b, c, maxit=-1, tol=1e-10):
    """Minimize c @ x subject to A @ x = b and x >= 0.

    This preserves the Numerical Tours return convention (the primal vector),
    using the maintained HiGHS solver instead of the old dense simplex port.
    Infeasible, unbounded, and iteration-limited problems raise an error.
    """
    options = {
        "dual_feasibility_tolerance": max(tol, 1e-10),
        "primal_feasibility_tolerance": max(tol, 1e-10),
    }
    if maxit >= 0:
        options["maxiter"] = int(maxit)
    result = linprog(
        np.asarray(c).ravel(order="F"),
        A_eq=A,
        b_eq=np.asarray(b).ravel(order="F"),
        bounds=(0, None),
        method="highs",
        options=options,
    )
    if not result.success:
        raise RuntimeError(f"Linear programming failed: {result.message}")
    return result.x
