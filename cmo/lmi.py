"""Common-P D-region certificate over the nominal operating polytope (Eqs. 23-25).

The augmented closed-loop matrix is affine in (rho0(d0), Ed(IL0)); the four
vertices of d0 in [d(1 V), d(46 V, +0.5 A)] x IL0 in [0.1, 5.1] A are tested with
one common P. Bisection gives the largest decay alpha_b and the largest
damping zeta_b; the certified values carry the 1 % safety factor of Eq. (25).
"""
from __future__ import annotations

import numpy as np
import cvxpy as cp
from scipy.linalg import matrix_balance

from .params import PLANT
from .plant import equilibrium, gamma


def vertices(K, P=PLANT):
    g = gamma(P)
    _, dlo = equilibrium(1.0, 0.0, P)
    _, dhi = equilibrium(46.0, 0.5, P)
    ILs = (1.0 / P.R, 46.0 / P.R + 0.5)
    out = []
    for d0 in (dlo, dhi):
        for IL in ILs:
            rho0 = P.rL + P.rd + d0 * (P.Ron - P.rd)
            Ed = P.Vin + P.Vf - (P.Ron - P.rd) * IL
            A = np.array([[-(rho0 + g * P.rC) / P.L, -g / P.L, 0.0],
                          [g / P.C, -g / (P.R * P.C), 0.0],
                          [-g * P.rC, -g, 0.0]])
            B = np.array([[Ed / P.L], [0.0], [0.0]])
            out.append(A - B @ np.asarray(K)[None, :])
    return out, (dlo, dhi), ILs


def _feasible(Acls, kind, val, scale):
    n = Acls[0].shape[0]
    X = cp.Variable((n, n), symmetric=True)
    cons = [X >> 1e-3 * np.eye(n)]
    eps = 1e-4
    for A in Acls:
        A = A / scale
        M = A @ X
        if kind == 'alpha':
            cons.append(M + M.T + 2 * (val / scale) * X << -eps * np.eye(n))
        else:
            th = np.arccos(val)
            s, c = np.sin(th), np.cos(th)
            blk = cp.bmat([[s * (M + M.T), c * (M - M.T)],
                           [c * (M.T - M), s * (M + M.T)]])
            cons.append(0.5 * (blk + blk.T) << -eps * np.eye(2 * n))
    prob = cp.Problem(cp.Minimize(0), cons)
    try:
        prob.solve(solver='CLARABEL')
    except Exception:
        return False, None
    ok = prob.status in ('optimal', 'optimal_inaccurate') and prob.status == 'optimal'
    return ok, (X.value if ok else None)


def certificate(K, P=PLANT, iters=40):
    Acls, drange, ILs = vertices(K, P)
    # balance with the nominal vertex to condition the SDP
    _, (sc, perm) = matrix_balance(Acls[0], permute=False, separate=True)
    T = np.diag(sc)
    Ti = np.diag(1 / sc)
    Ab = [Ti @ A @ T for A in Acls]
    scale = max(np.max(np.abs(np.linalg.eigvals(A))) for A in Ab)
    lo, hi = 0.0, min(-np.max(np.linalg.eigvals(A).real) for A in Ab)
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        ok, _ = _feasible(Ab, 'alpha', mid, scale)
        lo, hi = (mid, hi) if ok else (lo, mid)
    alpha_b = lo
    zlo, zhi = 0.0, 1.0
    for _ in range(iters):
        mid = 0.5 * (zlo + zhi)
        ok, _ = _feasible(Ab, 'zeta', mid, scale)
        zlo, zhi = (mid, zhi) if ok else (zlo, mid)
    zeta_b = zlo
    a_c, z_c = 0.99 * alpha_b, 0.99 * zeta_b
    okA, XA = _feasible(Ab, 'alpha', a_c, scale)
    okZ, XZ = _feasible(Ab, 'zeta', z_c, scale)
    # independent strictness check on the recovered matrices
    def worst(X, kind, val):
        w = -np.inf
        for A in Ab:
            A = A / scale
            M = A @ X
            if kind == 'alpha':
                L = M + M.T + 2 * (val / scale) * X
            else:
                th = np.arccos(val); s, c = np.sin(th), np.cos(th)
                L = np.block([[s * (M + M.T), c * (M - M.T)], [c * (M.T - M), s * (M + M.T)]])
            w = max(w, np.max(np.linalg.eigvalsh(0.5 * (L + L.T))))
        return w
    return dict(alpha_b=alpha_b, zeta_b=zeta_b, alpha_cert=a_c, zeta_cert=z_c,
                lamminP_alpha=float(np.min(np.linalg.eigvalsh(XA))) if okA else np.nan,
                lmax_alpha=float(worst(XA, 'alpha', a_c)) if okA else np.nan,
                lamminP_zeta=float(np.min(np.linalg.eigvalsh(XZ))) if okZ else np.nan,
                lmax_zeta=float(worst(XZ, 'zeta', z_c)) if okZ else np.nan,
                d_range=drange, IL_range=ILs,
                vertex_poles=[np.linalg.eigvals(A) for A in Acls])
