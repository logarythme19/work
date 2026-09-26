"""Algorithms 1-7: six equal-budget global operators and refineBO.

Common rules: one counted call = one six-scenario averaged evaluation plus the
37-inequality audit; identical normalized box, matched seed and initial
population of 12 (xi0 inserted first, 60 % of the individuals within a
normalized radius 0.05 of xi0, the rest uniform); Deb feasibility-first
ordering (feasible by J, infeasible by aggregate violation). Global budget 240,
local budget 60, total 300. Labels pAEABC, pIGWO, pIGWO-DLH denote project
variants, not operator-identical reproductions of the cited algorithms.
"""
from __future__ import annotations

import numpy as np
from scipy.stats import norm
from scipy.linalg import cho_factor, cho_solve

NPOP = 12
B_GLOBAL = 240
B_LOCAL = 60
METHODS = ('PSO', 'GA', 'ABC', 'pAEABC', 'pIGWO', 'pIGWO-DLH')


class Budget(Exception):
    pass


def better(a, b):
    if b is None:
        return True
    if a.feasible != b.feasible:
        return a.feasible
    if a.feasible:
        return a.J < b.J
    return a.viol < b.viol


def best_of(recs):
    b = None
    for r in recs:
        if better(r, b):
            b = r
    return b


def order(recs):
    return sorted(range(len(recs)), key=lambda i: (0, recs[i].J) if recs[i].feasible else (1, recs[i].viol))


class Counter:
    """Normalized-space wrapper: u in [0,1]^d <-> x = lb + u (ub - lb)."""

    def __init__(self, fn, lb, ub, budget, stage):
        self.fn, self.lb, self.ub = fn, np.asarray(lb, float), np.asarray(ub, float)
        self.budget, self.stage = budget, stage
        self.n = 0
        self.hist = []

    def __call__(self, u):
        if self.n >= self.budget:
            raise Budget()
        u = np.clip(np.asarray(u, float), 0.0, 1.0)
        rec = self.fn(self.lb + u * (self.ub - self.lb))
        rec.u = u
        self.n += 1
        rec.n_eval = self.n
        rec.stage = self.stage
        self.hist.append(rec)
        return rec

    def tau(self):
        return min(1.0, self.n / self.budget)


def reflect(u):
    u = np.where(u < 0, -u, u)
    u = np.where(u > 1, 2 - u, u)
    return np.clip(u, 0.0, 1.0)


def init_pop(rng, u0, d):
    X = rng.uniform(0, 1, (NPOP, d))
    nloc = int(round(0.6 * NPOP))
    for i in range(1, nloc):
        z = rng.normal(size=d)
        X[i] = u0 + 0.05 * rng.random() * z / np.linalg.norm(z)
    X[0] = u0
    return reflect(X)


# ---------------------------------------------------------------- Alg. 1 PSO
def pso(ev, rng, u0):
    d = len(u0)
    X = init_pop(rng, u0, d)
    V = np.zeros_like(X)
    P = [ev(x) for x in X]
    Px = X.copy()
    G = best_of(P)
    while True:
        w = 0.9 - 0.5 * ev.tau()
        for i in range(NPOP):
            V[i] = w * V[i] + 1.7 * rng.random(d) * (Px[i] - X[i]) + 1.7 * rng.random(d) * (G.u - X[i])
            V[i] = np.clip(V[i], -0.2, 0.2)
            X[i] = reflect(X[i] + V[i])
            r = ev(X[i])
            if better(r, P[i]):
                P[i] = r; Px[i] = X[i].copy()
            if better(r, G):
                G = r


# ----------------------------------------------------------------- Alg. 2 GA
def ga(ev, rng, u0):
    d = len(u0)
    pop = [ev(x) for x in init_pop(rng, u0, d)]
    while True:
        def tour():
            i, j = rng.integers(0, NPOP, 2)
            return pop[i] if better(pop[i], pop[j]) else pop[j]
        kids = []
        for _ in range(NPOP // 2):
            a, b = tour().u, tour().u
            uu = rng.random(d)
            beta = np.where(uu <= 0.5, (2 * uu) ** (1 / 16), (1 / (2 * (1 - uu))) ** (1 / 16))
            for c in (0.5 * ((1 + beta) * a + (1 - beta) * b), 0.5 * ((1 - beta) * a + (1 + beta) * b)):
                c = c.copy()
                for j in range(d):
                    if rng.random() < 1.0 / d:
                        q = rng.random()
                        dl = (2 * q) ** (1 / 21) - 1 if q < 0.5 else 1 - (2 * (1 - q)) ** (1 / 21)
                        c[j] += dl
                kids.append(ev(reflect(c)))
        merged = pop + kids
        pop = [merged[i] for i in order(merged)[:NPOP]]


# ------------------------------------------------ Alg. 3 ABC / Alg. 4 pAEABC
def _abc(ev, rng, u0, adaptive):
    d = len(u0)
    pop = [ev(x) for x in init_pop(rng, u0, d)]
    trials = np.zeros(NPOP, int)
    limit = max(5, int(round(0.5 * NPOP * d)))

    def move(i, sc):
        k = rng.integers(0, NPOP - 1)
        k = k + 1 if k >= i else k
        j = rng.integers(0, d)
        y = pop[i].u.copy()
        y[j] += sc * rng.uniform(-1, 1) * (y[j] - pop[k].u[j])
        return y

    while True:
        tau = ev.tau()
        sc = 1.0
        if adaptive:
            Dk = float(np.mean(np.std([r.u for r in pop], axis=0)))
            sc = float(np.clip(0.5 + 0.15 / max(Dk, 1e-3), 0.5, 2.0))
        for i in range(NPOP):
            y = move(i, sc)
            if adaptive and rng.random() < 0.25 * (1 - tau):
                j2 = rng.integers(0, d)
                y[j2] += rng.uniform(-0.1, 0.1)
            r = ev(reflect(y))
            if better(r, pop[i]):
                pop[i] = r; trials[i] = 0
            else:
                trials[i] += 1
        o = order(pop)
        rk = np.empty(NPOP); rk[o] = np.arange(NPOP)
        pr = (NPOP - rk) / np.sum(NPOP - rk)
        for _ in range(NPOP):
            i = int(rng.choice(NPOP, p=pr))
            r = ev(reflect(move(i, sc)))
            if better(r, pop[i]):
                pop[i] = r; trials[i] = 0
            else:
                trials[i] += 1
        i = int(np.argmax(trials))
        if trials[i] >= limit:
            pop[i] = ev(rng.random(d)); trials[i] = 0


def abc(ev, rng, u0):
    _abc(ev, rng, u0, False)


def aeabc(ev, rng, u0):
    _abc(ev, rng, u0, True)


# -------------------------------------------- Alg. 5 pIGWO / Alg. 6 pIGWO-DLH
def _leaders(pop):
    o = order(pop)
    return pop[o[0]].u, pop[o[1]].u, pop[o[2]].u


def _gwo(rng, x, L, a, d):
    acc = np.zeros(d)
    for xm in L:
        A = 2 * a * rng.random(d) - a
        Cm = 2 * rng.random(d)
        acc += xm - A * np.abs(Cm * xm - x)
    return acc / 3


def igwo(ev, rng, u0):
    d = len(u0)
    pop = [ev(x) for x in init_pop(rng, u0, d)]
    while True:
        tau = ev.tau()
        a = 2 * (1 - tau ** 2)
        L = _leaders(pop)
        cand = []
        for i in range(NPOP):
            y = _gwo(rng, pop[i].u, L, a, d)
            if rng.random() < 0.15 * (1 - tau) + 0.02:
                j = rng.integers(0, d)
                y[j] += rng.normal(0, 0.1)
            cand.append(ev(reflect(y)))
        merged = pop + cand
        pop = [merged[i] for i in order(merged)[:NPOP]]


def igwo_dlh(ev, rng, u0):
    d = len(u0)
    pop = [ev(x) for x in init_pop(rng, u0, d)]
    K = int(np.ceil((ev.budget - NPOP) / (2 * NPOP)))
    k = 0
    while True:
        a = max(0.0, 2 - 2 * k / K)
        L = _leaders(pop)
        for i in range(NPOP):
            r1 = ev(np.clip(_gwo(rng, pop[i].u, L, a, d), 0, 1))
            Ri = 0.5 * np.linalg.norm(pop[i].u - L[0])
            z = rng.normal(size=d)
            try:
                r2 = ev(np.clip(pop[i].u + Ri * rng.random() * z / np.linalg.norm(z), 0, 1))
            except Budget:
                pop[i] = r1
                raise
            pop[i] = r1 if better(r1, r2) else r2      # non-elitist replacement
        k += 1


GLOBAL = {'PSO': pso, 'GA': ga, 'ABC': abc, 'pAEABC': aeabc, 'pIGWO': igwo, 'pIGWO-DLH': igwo_dlh}


# -------------------------------------------------------- Alg. 7 refineBO
def _matern52(A, B, ell):
    r = np.sqrt(np.maximum(((A[:, None, :] - B[None, :, :]) ** 2).sum(-1), 0.0)) / ell
    s5 = np.sqrt(5) * r
    return (1 + s5 + 5 * r * r / 3) * np.exp(-s5)


class GP:
    """Shared-kernel GP for J and every constraint (one Cholesky per refit)."""

    def __init__(self, U, Y):
        self.U = U
        mu = Y.mean(0); sd = Y.std(0); sd[sd < 1e-12] = 1.0
        self.mu, self.sd = mu, sd
        Z = (Y - mu) / sd
        best, bestll = None, -np.inf
        n = len(U)
        for ell in (0.05, 0.1, 0.2, 0.4, 0.8):
            Kxx = _matern52(U, U, ell) + 1e-6 * np.eye(n)
            try:
                cf = cho_factor(Kxx)
            except np.linalg.LinAlgError:
                continue
            a = cho_solve(cf, Z[:, 0])
            ll = -0.5 * Z[:, 0] @ a - np.sum(np.log(np.diag(cf[0])))
            if ll > bestll:
                best, bestll = (ell, cf), ll
        self.ell, self.cf = best
        self.alpha = cho_solve(self.cf, Z)

    def predict(self, V):
        Ks = _matern52(V, self.U, self.ell)
        m = Ks @ self.alpha
        v = cho_solve(self.cf, Ks.T)
        var = np.maximum(1.0 - np.sum(Ks * v.T, 1), 1e-12)
        return m * self.sd + self.mu, np.sqrt(var)[:, None] * self.sd


def refine_bo(ev, rng, hist):
    d = len(hist[0].u)
    allr = list(hist)
    for k in range(ev.budget):
        U = np.array([r.u for r in allr])
        Jv = np.array([min(r.J, 5.0) for r in allr])
        G = np.array([np.clip(r.g, -5, 5) for r in allr])
        gp = GP(U, np.column_stack([Jv, G]))
        rad = 0.04 * (0.003 / 0.04) ** (k / max(ev.budget - 1, 1))
        o = order(allr)[:10]
        C = np.array([allr[i].u for i in o])
        z = rng.normal(size=(1000, d))
        z *= (rng.random((1000, 1)) ** (1 / d)) / np.linalg.norm(z, axis=1, keepdims=True)
        pool = np.clip(C[rng.integers(0, len(C), 1000)] + rad * z, 0, 1)
        m, s = gp.predict(pool)
        mJ, sJ = m[:, 0], s[:, 0]
        mg, sg = m[:, 1:], s[:, 1:]
        pj = norm.cdf(-mg / sg)
        Pf = np.prod(pj, axis=1)
        inc = best_of(allr)
        if inc.feasible:
            zz = (inc.J - mJ) / sJ
            EI = (inc.J - mJ) * norm.cdf(zz) + sJ * norm.pdf(zz)
            acq = Pf * (EI + 0.02 * sJ) + 1e-4 * Pf
            i = int(np.argmax(acq))
        else:          # predeclared restoration mode
            acq = np.sum(np.maximum(mg, 0), 1) - 0.1 * np.mean(sg, 1)
            i = int(np.argmin(acq))
        r = ev(pool[i])
        allr.append(r)
    return allr


def run_method(method, fn, lb, ub, x0, seed):
    """One matched-seed run: 240 global + 60 refineBO calls. Returns history."""
    lb, ub = np.asarray(lb, float), np.asarray(ub, float)
    u0 = (np.asarray(x0, float) - lb) / (ub - lb)
    rng = np.random.default_rng(seed)
    evg = Counter(fn, lb, ub, B_GLOBAL, 'global')
    try:
        GLOBAL[method](evg, rng, u0)
    except Budget:
        pass
    evl = Counter(fn, lb, ub, B_LOCAL, 'local')
    try:
        refine_bo(evl, rng, evg.hist)
    except Budget:
        pass
    for r in evl.hist:
        r.n_eval += B_GLOBAL
    return evg.hist + evl.hist
