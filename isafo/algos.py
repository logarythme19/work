"""Les six metaheuristiques, a budget identique et graines appariees.

Operateurs repris de l'article anterieur (table 6), transposes a la boite de
decision v5 de dimension 7. Les prefixes "p" (pAEABC, pIGWO, pIGWO-DLH) y
marquaient l'absence d'identite prouvee avec les publications citees; la meme
reserve vaut ici et doit etre conservee dans la redaction.

Regles communes, non negociables (par.3.2):
  - un appel d'evaluateur = une mission complete, jamais une mise a jour
    interne de vitesse, de croisement, d'abeille ou de loup;
  - meme boite normalisee, meme population initiale pour une graine donnee,
    meme comparateur de Deb;
  - le budget est compte en appels, pas en iterations: pIGWO-DLH consomme deux
    appels par comparaison, ce qui lui donne moins d'iterations a budget egal.
    C'est voulu: c'est le budget qui est apparie, pas le nombre de tours.
"""

from __future__ import annotations

from typing import Callable, List, Optional

import numpy as np

from .deb import Record, deb_better, argsort_deb, best_of

METHODS = ("PSO", "GA", "ABC", "AEABC", "IGWO", "IGWO_DLH")
DEFAULT_POP = 10          # par.3.2: population 10 sauf justification


class Budget(Exception):
    """Levee des qu'un algorithme epuise son budget d'appels."""


class Evaluator:
    """Compte les appels, archive l'historique, coupe net a l'epuisement."""

    def __init__(self, fn: Callable[[np.ndarray], Record], budget: int,
                 stage: str = ""):
        self.fn = fn
        self.budget = int(budget)
        self.n = 0
        self.history: List[Record] = []
        self.stage = stage

    def __call__(self, x: np.ndarray) -> Record:
        if self.n >= self.budget:
            raise Budget()
        self.n += 1
        rec = self.fn(np.asarray(x, dtype=float))
        rec.n_eval = self.n
        rec.stage = self.stage
        self.history.append(rec)
        return rec

    @property
    def remaining(self) -> int:
        return self.budget - self.n

    def tau(self) -> float:
        """Temps de budget normalise, utilise par les schedules."""
        return min(1.0, self.n / max(self.budget, 1))


def _reflect_clip(x, lb, ub):
    """Reflexion unique puis ecretage (regle de l'article anterieur)."""
    x = np.where(x < lb, 2 * lb - x, x)
    x = np.where(x > ub, 2 * ub - x, x)
    return np.clip(x, lb, ub)


def _init_pop(rng, lb, ub, n_pop, seeds: Optional[np.ndarray] = None):
    X = rng.uniform(lb, ub, size=(n_pop, len(lb)))
    if seeds is not None and len(seeds) > 0:
        s = np.atleast_2d(np.asarray(seeds, dtype=float))
        k = min(len(s), n_pop)
        X[:k] = np.clip(s[:k], lb, ub)
    return X


def _eval_pop(ev, X):
    return [ev(x) for x in X]


# ----------------------------------------------------------------- PSO
def run_pso(ev, rng, lb, ub, n_pop=DEFAULT_POP, seeds=None, mode="margin"):
    d = len(lb)
    X = _init_pop(rng, lb, ub, n_pop, seeds)
    V = np.zeros_like(X)
    vmax = 0.2 * (ub - lb)
    c1 = c2 = 1.7
    recs = _eval_pop(ev, X)
    P = [r for r in recs]
    Px = X.copy()
    G = best_of(recs, mode)
    Gx = G.x.copy()
    while True:
        w = 0.9 - 0.5 * ev.tau()
        for i in range(n_pop):
            r1 = rng.random(d); r2 = rng.random(d)
            V[i] = (w * V[i] + c1 * r1 * (Px[i] - X[i]) + c2 * r2 * (Gx - X[i]))
            V[i] = np.clip(V[i], -vmax, vmax)
            X[i] = _reflect_clip(X[i] + V[i], lb, ub)
            rec = ev(X[i])
            if deb_better(rec, P[i], mode):
                P[i] = rec; Px[i] = X[i].copy()
            if deb_better(rec, G, mode):
                G = rec; Gx = X[i].copy()


# ------------------------------------------------------------------ GA
def _sbx(rng, p1, p2, lb, ub, eta=15.0):
    u = rng.random(len(p1))
    beta = np.where(u <= 0.5, (2 * u) ** (1 / (eta + 1)),
                    (1 / (2 * (1 - u))) ** (1 / (eta + 1)))
    c1 = 0.5 * ((1 + beta) * p1 + (1 - beta) * p2)
    c2 = 0.5 * ((1 - beta) * p1 + (1 + beta) * p2)
    return np.clip(c1, lb, ub), np.clip(c2, lb, ub)


def _poly_mut(rng, x, lb, ub, eta=20.0, pm=None):
    d = len(x)
    pm = (1.0 / d) if pm is None else pm
    y = x.copy()
    for j in range(d):
        if rng.random() < pm:
            u = rng.random()
            delta = ((2 * u) ** (1 / (eta + 1)) - 1 if u < 0.5
                     else 1 - (2 * (1 - u)) ** (1 / (eta + 1)))
            y[j] = x[j] + delta * (ub[j] - lb[j])
    return np.clip(y, lb, ub)


def run_ga(ev, rng, lb, ub, n_pop=DEFAULT_POP, seeds=None, mode="margin"):
    X = _init_pop(rng, lb, ub, n_pop, seeds)
    pop = _eval_pop(ev, X)
    while True:
        def tour():
            i, j = rng.integers(0, len(pop), size=2)
            return pop[i] if deb_better(pop[i], pop[j], mode) else pop[j]
        children = []
        for _ in range(n_pop // 2):
            a, b = tour(), tour()
            c1, c2 = _sbx(rng, a.x, b.x, lb, ub)
            for c in (c1, c2):
                children.append(ev(_poly_mut(rng, c, lb, ub)))
        merged = pop + children
        order = argsort_deb(merged, mode)[:n_pop]
        pop = [merged[i] for i in order]


# ------------------------------------------------------- ABC et pAEABC
def _abc_core(ev, rng, lb, ub, n_pop, seeds, mode, adaptive):
    d = len(lb)
    X = _init_pop(rng, lb, ub, n_pop, seeds)
    pop = _eval_pop(ev, X)
    trials = np.zeros(n_pop, dtype=int)
    limit = max(5, int(round(0.5 * n_pop * d)))

    def move(i, scale):
        k = rng.integers(0, n_pop)
        while k == i:
            k = rng.integers(0, n_pop)
        j = rng.integers(0, d)
        y = pop[i].x.copy()
        phi = rng.uniform(-1.0, 1.0) * scale
        y[j] = y[j] + phi * (y[j] - pop[k].x[j])
        return np.clip(y, lb, ub)

    while True:
        tau = ev.tau()
        if adaptive:
            Xa = np.array([r.x for r in pop])
            Dk = float(np.mean(np.std(Xa, axis=0)))
            scale = float(np.clip(0.5 + 0.15 / max(Dk, 1e-3), 0.5, 2.0))
        else:
            scale = 1.0

        for i in range(n_pop):                      # abeilles employees
            y = move(i, scale)
            if adaptive and rng.random() < 0.25 * (1.0 - tau):
                j2 = rng.integers(0, d)
                y[j2] = np.clip(y[j2] + rng.normal(0.0, 0.1) * (ub[j2] - lb[j2]),
                                lb[j2], ub[j2])
            rec = ev(y)
            if deb_better(rec, pop[i], mode):
                pop[i] = rec; trials[i] = 0
            else:
                trials[i] += 1

        order = argsort_deb(pop, mode)              # abeilles spectatrices
        ranks = np.empty(n_pop); ranks[order] = np.arange(n_pop)
        p = (n_pop - ranks); p = p / p.sum()
        for _ in range(n_pop):
            i = int(rng.choice(n_pop, p=p))
            rec = ev(move(i, scale))
            if deb_better(rec, pop[i], mode):
                pop[i] = rec; trials[i] = 0
            else:
                trials[i] += 1

        for i in range(n_pop):                      # eclaireuses
            if trials[i] >= limit:
                pop[i] = ev(rng.uniform(lb, ub)); trials[i] = 0


def run_abc(ev, rng, lb, ub, n_pop=DEFAULT_POP, seeds=None, mode="margin"):
    _abc_core(ev, rng, lb, ub, n_pop, seeds, mode, adaptive=False)


def run_aeabc(ev, rng, lb, ub, n_pop=DEFAULT_POP, seeds=None, mode="margin"):
    _abc_core(ev, rng, lb, ub, n_pop, seeds, mode, adaptive=True)


# ---------------------------------------------------- pIGWO et DLH
def _leaders(pop, mode):
    o = argsort_deb(pop, mode)
    return pop[o[0]].x, pop[o[min(1, len(o)-1)]].x, pop[o[min(2, len(o)-1)]].x


def _gwo_trial(rng, xi, xa, xb, xd, a, d):
    acc = np.zeros(d)
    for xm in (xa, xb, xd):
        A = 2 * a * rng.random(d) - a
        Cc = 2 * rng.random(d)
        D = np.abs(Cc * xm - xi)
        acc += xm - A * D
    return acc / 3.0


def run_igwo(ev, rng, lb, ub, n_pop=DEFAULT_POP, seeds=None, mode="margin"):
    d = len(lb)
    X = _init_pop(rng, lb, ub, n_pop, seeds)
    pop = _eval_pop(ev, X)
    while True:
        tau = ev.tau()
        a = 2.0 * (1.0 - tau ** 2)                  # schedule non lineaire
        xa, xb, xd = _leaders(pop, mode)
        cand = []
        for i in range(n_pop):
            y = _gwo_trial(rng, pop[i].x, xa, xb, xd, a, d)
            if rng.random() < 0.15 * (1.0 - tau) + 0.02:
                j = rng.integers(0, d)
                y[j] += rng.normal(0.0, 0.1) * (ub[j] - lb[j])
            cand.append(ev(np.clip(y, lb, ub)))
        merged = pop + cand
        order = argsort_deb(merged, mode)[:n_pop]
        pop = [merged[i] for i in order]


def run_igwo_dlh(ev, rng, lb, ub, n_pop=DEFAULT_POP, seeds=None, mode="margin"):
    """pIGWO-DLH: schedule lineaire + chasse locale radiale EVALUEE.

    Deux propositions comptees par loup et par tour, donc deux fois moins de
    tours a budget egal. Le remplacement est non elitiste: le loup courant cede
    la place a la meilleure des deux propositions.
    """
    d = len(lb)
    X = _init_pop(rng, lb, ub, n_pop, seeds)
    pop = _eval_pop(ev, X)
    K = max(1, (ev.budget - n_pop) // (2 * n_pop))
    k = 0
    while True:
        a = max(0.0, 2.0 - 2.0 * k / K)             # schedule lineaire
        xa, xb, xd = _leaders(pop, mode)
        for i in range(n_pop):
            y1 = np.clip(_gwo_trial(rng, pop[i].x, xa, xb, xd, a, d), lb, ub)
            r1 = ev(y1)

            Ri = 0.5 * float(np.linalg.norm(pop[i].x - xa))
            z = rng.normal(size=d)
            nz = np.linalg.norm(z)
            q = z / nz if nz > 0 else np.zeros(d)
            y2 = np.clip(pop[i].x + Ri * rng.random() * q, lb, ub)
            r2 = ev(y2)

            better = r1 if deb_better(r1, r2, mode) else r2
            if deb_better(better, pop[i], mode):
                pop[i] = better
        k += 1


RUNNERS = {
    "PSO": run_pso,
    "GA": run_ga,
    "ABC": run_abc,
    "AEABC": run_aeabc,
    "IGWO": run_igwo,
    "IGWO_DLH": run_igwo_dlh,
}


def run(method: str, fn, budget: int, seed: int, lb, ub,
        n_pop: int = DEFAULT_POP, seeds=None, mode: str = "margin",
        stage: str = "") -> List[Record]:
    """Execute une methode jusqu'a epuisement EXACT du budget d'appels."""
    if method not in RUNNERS:
        raise ValueError(f"methode inconnue: {method}")
    ev = Evaluator(fn, budget, stage=stage)
    rng = np.random.default_rng(seed)
    try:
        RUNNERS[method](ev, rng, np.asarray(lb), np.asarray(ub),
                        n_pop=n_pop, seeds=seeds, mode=mode)
    except Budget:
        pass
    return ev.history
