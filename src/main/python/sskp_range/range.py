#!/usr/bin/env python3
"""
Stochastic Knapsack Problem –  Shortest-Path Dynamic Programming Solver

Reference
----------
Range, T. M., Kozlowski, D., Petersen, N. C. (2018).
	“A shortest-path-based approach for the stochastic knapsack problem     
 	with non-decreasing expected overfilling costs.” Computers & OR 97.    
"""

###############################################################################
# 0. Preliminaries – imports and normal-distribution utilities                #
###############################################################################
"""
The solver needs only a handful of standard modules:

• math      – elementary mathematics ( π , exp , erf , sqrt , …)
• time      – wall-clock timing for performance measurement
• json      – loading instance data for the batch runner
• collections.defaultdict – light-weight adjacency lists and label sets
• numpy     – used only by the batch-runner to present solutions nicely

Two tiny helper functions ─  φ  and  Φ  ─ wrap the PDF and CDF of the unit
normal distribution; these are required by the  h_T(·,·)  approximation
(see Section 3.1 of the paper).
"""

import math
import time
import json
import collections
from typing import Dict, List, Tuple, Set, FrozenSet, Optional

import numpy as np   # optional – only for the tiny batch driver


# ---  ϕ  and  Φ  --------------------------------------------------------------
def phi(x: float) -> float:
    """Standard normal probability density function  ϕ(x) ."""
    return (1.0 / math.sqrt(2.0 * math.pi)) * math.exp(-x * x / 2.0)


def Phi(x: float) -> float:
    """Standard normal cumulative distribution function  Φ(x) ."""
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


###############################################################################
# 1. Expected overfill approximation  h_T(μ,σ)                                 #
###############################################################################
"""
Whenever we know the current *mean* μ and *standard deviation* σ of the
(knapsack) weight, we can approximate the *expected* overfill O = max(0, X−T)
via Proposition 1 of the paper:

    h_T(μ,σ) = σ · [ ϕ(k) − k · (1 − Φ(k)) ]      with k = (T − μ)/σ.

Edge case σ≈0 is interpreted as “deterministic weight”, in which the
expected overfill coincides with the positive part of μ − T.
"""


def h_T(mu: float, sigma: float, T: float) -> float:
    """Return the approximation  E[O]  for given (μ,σ,T)."""
    if sigma <= 1e-9:                      # treat almost-deterministic weight
        return max(0.0, mu - T)
    k = (T - mu) / sigma
    return sigma * (phi(k) - k * (1.0 - Phi(k)))


###############################################################################
# 2. Penalty functions  f(x)                                                   #
###############################################################################
"""
Overfilling the knapsack is not prohibited outright – it is *penalised*.
The current implementation supports two convex, non-decreasing variants:

    • linear    f(x) = slope · x
    • quadratic f(x) = ½ x²                        (default)

Both respect   f(0) = 0   and ignore negative arguments.
"""


def penalty_function(x: float,
                     kind: str = "quadratic",
                     slope: float = 10.0) -> float:
    """Evaluate the penalty function f(x)."""
    if x <= 0.0:
        return 0.0
    if kind == "linear":
        return slope * x
    if kind == "quadratic":
        return 0.5 * x * x
    raise ValueError("penalty_function: 'kind' must be 'linear' or 'quadratic'")


###############################################################################
# 3. Feasibility of a *partial* path                                           #
###############################################################################
"""
A label (state) is *feasible* iff the chance-constraint approximation
(Equation 7) is satisfied

        μ(P) + β·σ(P)  ≤  T + S .

The helper below is intentionally small because it is executed *very*
often during the DP search.
"""


def is_feasible(mu: float, sigma_sq: float,
                T: float, S: float, beta: float) -> bool:
    """Return True iff the chance constraint holds."""
    return mu + beta * math.sqrt(max(0.0, sigma_sq)) <= T + S


###############################################################################
# 4. Cost of a label                                                           #
###############################################################################
"""
The dynamic programme *minimises*   C(P) = f(h_T(…)) − R(P).
Hence the *profit* we are really after equals  −C(P).
"""


def label_cost(revenue: float, mu: float, sigma_sq: float,
               T: float, pen_kind: str, pen_slope: float) -> float:
    """Compute C(P) for a (partial) path."""
    exp_overfill = h_T(mu, math.sqrt(sigma_sq), T)
    return penalty_function(exp_overfill, pen_kind, pen_slope) - revenue


###############################################################################
# 5. Data structure for a *state / label*                                     #
###############################################################################
"""
Every node keeps a *set* of undominated states (§4, Proposition 4).  A state
remembers exactly the information needed for further extension and cost
computation.  For hashing and quick subset tests, the chosen items are stored
as a **frozenset**.
"""


class PathState:
    """Immutable container for all data associated with one DP label."""

    __slots__ = ("items", "r", "mu", "sigma_sq",
                 "cost", "last_item", "groups")

    def __init__(self,
                 items: FrozenSet[int],
                 r: float,
                 mu: float,
                 sigma_sq: float,
                 cost: float,
                 last_item: int,
                 groups: FrozenSet[int]):
        self.items: FrozenSet[int] = items
        self.r: float = r
        self.mu: float = mu
        self.sigma_sq: float = sigma_sq
        self.cost: float = cost
        self.last_item: int = last_item       # −1: origin
        self.groups: FrozenSet[int] = groups  # groups already used (if any)

    # ---------------------------------------------------------------------
    # Pretty printing and ordering (the latter is needed for “best-first”
    # selection inside the incremental extension loop)
    # ---------------------------------------------------------------------
    def __repr__(self) -> str:
        return (f"State(R={self.r:.2f}, μ={self.mu:.2f}, σ²={self.sigma_sq:.2f}, "
                f"C={self.cost:.2f}, last={self.last_item}, "
                f"items={sorted(self.items)})")

    def __lt__(self, other: 'PathState') -> bool:
        return self.cost < other.cost


###############################################################################
# 6. Dominance check  (Proposition 4)                                          #
###############################################################################
"""
State  s₁  *dominates*  s₂  if it is **no worse** in μ and σ²
   and **no smaller** in revenue.  In that case s₂ can be discarded.
"""


def dominates(s1: PathState, s2: PathState) -> bool:
    return (s1.r >= s2.r
            and s1.mu <= s2.mu
            and s1.sigma_sq <= s2.sigma_sq)


def eliminate_dominated(label_set: List[PathState],
                        new_state: PathState) -> bool:
    """
    Insert  new_state  into  label_set  *unless* it is dominated.
    Simultaneously remove any existing labels that are dominated by
    the newcomer.

    Returns  True  iff  new_state  has *already* been dominated.
    """
    dominated_by_new = []
    for old in label_set:
        if dominates(new_state, old):
            dominated_by_new.append(old)
        elif dominates(old, new_state):
            return True                     # new_state is useless

    for old in dominated_by_new:
        label_set.remove(old)              # destructive update
    return False


###############################################################################
# 7. Bounding the cost of *all* extensions  (Proposition 5)                    #
###############################################################################
"""
When the penalty function f is convex, we may obtain a *lower* bound on the
*best possible* additional cost for completing the current partial path.
If that bound added to the *current* cost is still ≥ incumbent best, the
label can be pruned *immediately* (without trying to extend it even once).

The following helper routines implement the mathematical construction
given in Section 4.1 of the paper.
"""
# -----------------------------------------------------------------------------


def build_remaining_groups(state: PathState,
                           items: Dict[int, Dict],
                           groups: Optional[Dict[int, List[int]]],
                           group_of_item: Dict[int, int]) -> List[Tuple[float, float]]:
    """
    For the fractional knapsack relaxation we need for every *unused* class k
        • Mk  =  max μ
        • Lk  =  min( −r/μ )

    We finally return a list [(Lk, Mk), …] already *sorted* by the product
    Mk·Lk (ascending) – cf. Equation 17, z(P,m).
    """
    selected_items = state.items
    last_gid = group_of_item.get(state.last_item, -1)

    tuples = []  # (Mk·Lk, Lk, Mk)

    if groups is None:                     # binary case – one item per group
        for j, dat in items.items():
            if j in selected_items or j <= state.last_item:
                continue
            Mk = dat['mu']
            Lk = -dat['r'] / Mk
            tuples.append((Mk * Lk, Lk, Mk))
    else:
        for gid, itms in groups.items():
            if gid <= last_gid or gid in state.groups:
                continue
            if any(i in selected_items for i in itms):
                continue
            Mk = max(items[i]['mu'] for i in itms)
            Lk = min(-items[i]['r'] / items[i]['mu'] for i in itms)
            tuples.append((Mk * Lk, Lk, Mk))

    tuples.sort()
    return [(Lk, Mk) for _, Lk, Mk in tuples]


def g_value(m: float, mu_P: float, T: float,
            pen_kind: str, pen_slope: float,
            segs: List[Tuple[float, float]]) -> float:
    """Evaluate  g_P(m)  = z(P,m) + f(max(0, μ(P)+m−T))."""
    z_val, rest = 0.0, m
    for Lk, Mk in segs:
        take = min(rest, Mk)
        z_val += Lk * take
        rest -= take
        if rest <= 0:
            break
    f_term = penalty_function(max(0.0, mu_P + m - T), pen_kind, pen_slope)
    return z_val + f_term


def min_g(mu_P: float, T: float, S: float,
          pen_kind: str, pen_slope: float,
          segs: List[Tuple[float, float]]) -> float:
    """
    Minimisation of g_P(m) on 0 ≤ m ≤ T+S−μ(P).
    Because g is convex and piecewise quadratic / linear, it suffices to
    test a finite set of candidate points (breakpoints plus at most one
    interior minimiser on each segment).
    """
    m_max = max(0.0, T + S - mu_P)
    if m_max == 0.0:
        return g_value(0.0, mu_P, T, pen_kind, pen_slope, segs)

    # candidate set initialised with both interval boundaries
    cand: Set[float] = {0.0, m_max}

    # breakpoints after every Mk
    cum = 0.0
    for Lk, Mk in segs:
        cum += Mk
        if cum < m_max:
            cand.add(cum)
        else:
            break

    # point where f becomes positive
    if 0.0 < T - mu_P < m_max:
        cand.add(T - mu_P)

    # analytic minima for the quadratic variant
    if pen_kind == "quadratic":
        cum = 0.0
        for Lk, Mk in segs:
            m_crit = T - mu_P - Lk        # derivative g'(m)=0
            seg_start = max(cum, T - mu_P)
            seg_end = min(cum + Mk, m_max)
            if seg_start <= m_crit <= seg_end:
                cand.add(m_crit)
            cum += Mk
            if cum >= m_max:
                break

    return min(g_value(m, mu_P, T, pen_kind, pen_slope, segs) for m in cand)


def extension_lower_bound(state: PathState, T: float, S: float,
                           pen_kind: str, pen_slope: float,
                           items: Dict[int, Dict],
                           groups: Optional[Dict[int, List[int]]],
                           group_of_item: Dict[int, int]) -> float:
    """Return the Proposition-5 bound LB(P)."""
    segs = build_remaining_groups(state, items, groups, group_of_item)
    if not segs:                               # nothing left to pack
        return 0.0

    min_g_val = min_g(state.mu, T, S, pen_kind, pen_slope, segs)
    f_hP = penalty_function(
        h_T(state.mu, math.sqrt(state.sigma_sq), T),
        pen_kind, pen_slope)
    return min_g_val - f_hP


###############################################################################
# 8. The main solver – Incremental label setting (Algorithm 2)                 #
###############################################################################
"""
The heart of the programme:  *solve_skp* creates the DAG, performs the
incremental label-setting DP, and finally returns the best solution found.
"""


def solve_skp(items: Dict[int, Dict[str, float]],
              T: float,
              S: float,
              beta: float,
              penalty_kind: str = "linear",
              penalty_slope: float = 10.0,
              item_groups: Optional[Dict[int, List[int]]] = None,
              L_init: int = 10,
              max_stages: int = 100,
              time_limit_sec: int = 600) -> Tuple[FrozenSet[int], float, float]:
    """
    Shortest-path dynamic programme for the (static) stochastic knapsack.

    Parameters
    ----------
    items : dict             item_id → {'r','mu','sigma_sq'}
    T,S,beta : float         capacity, allowed overfill, chance-constraint param
    penalty_kind : str       'linear' | 'quadratic'
    penalty_slope : float    slope for the linear variant
    item_groups : dict       group_id → list(item_id)     (None = binary case)
    L_init : int             initial per-node extension limit (Algorithm 2)
    max_stages : int         safety stopper, usually irrelevant
    time_limit_sec : int     wall-clock time limit

    Returns
    -------
    (items_selected,  profit (=R−f),  objective_cost (=−profit) )
    """

    # ------------------------------------------------------------------
    # 8.1  Build successor lists once & for all
    # ------------------------------------------------------------------
    sorted_ids = sorted(items)
    successors = collections.defaultdict(list)
    successors[-1] = sorted_ids[:]                   # from artificial origin
    for pos, itm in enumerate(sorted_ids):
        successors[itm] = sorted_ids[pos + 1:]

    if item_groups is None:
        group_of_item = {}
    else:
        group_of_item = {i: gid
                         for gid, lst in item_groups.items()
                         for i in lst}

    # ------------------------------------------------------------------
    # 8.2  Label containers and incumbents
    # ------------------------------------------------------------------
    Lambda: Dict[int, List[PathState]] = collections.defaultdict(list)
    Extended: Dict[int, Set[PathState]] = collections.defaultdict(set)

    origin = PathState(frozenset(), 0.0, 0.0, 0.0, 0.0,
                       last_item=-1, groups=frozenset())
    Lambda[-1].append(origin)

    best_profit = -math.inf
    best_cost = math.inf
    best_items: FrozenSet[int] = frozenset()

    # ------------------------------------------------------------------
    # 8.3  Outer “stage” loop   (L doubles every round)
    # ------------------------------------------------------------------
    L = L_init
    stage = 0
    t0 = time.time()

    while stage < max_stages:
        stage += 1
        some_extension = False

        # Θ_i – states not yet processed in *this* stage
        Theta = {v: [s for s in lst if s not in Extended[v]]
                 for v, lst in Lambda.items()}
        for lst in Theta.values():
            lst.sort()                             # ascending by cost

        for node in [-1] + sorted_ids:

            to_extend = Theta.get(node, [])[:L]
            if not to_extend:
                continue

            for state in to_extend:
                # --- time guard
                if time.time() - t0 > time_limit_sec:
                    return best_items, best_profit, best_cost

                Extended[node].add(state)
                some_extension = True

                # --- incumbent improvement
                if -state.cost > best_profit:
                    best_profit = -state.cost
                    best_cost = state.cost
                    best_items = state.items

                # --- bounding (convex f)
                if penalty_kind in ("linear", "quadratic"):
                    lb = extension_lower_bound(state, T, S,
                                               penalty_kind, penalty_slope,
                                               items, item_groups,
                                               group_of_item)
                    if state.cost + lb >= best_cost:
                        continue                    # safely prune

                # --------------------------------------------------
                # 8.4  *Extend* state to all admissible successors
                # --------------------------------------------------
                for nxt in successors[state.last_item]:

                    if nxt in state.items:
                        continue                    # no cycles

                    gid = group_of_item.get(nxt)
                    if gid is not None and gid in state.groups:
                        continue                    # at most one per class

                    new_r = state.r + items[nxt]['r']
                    new_mu = state.mu + items[nxt]['mu']
                    new_sigma_sq = state.sigma_sq + items[nxt]['sigma_sq']

                    if not is_feasible(new_mu, new_sigma_sq, T, S, beta):
                        continue

                    # incremental cost (Equation 14)
                    old_pen = penalty_function(
                        h_T(state.mu, math.sqrt(state.sigma_sq), T),
                        penalty_kind, penalty_slope)
                    new_pen = penalty_function(
                        h_T(new_mu, math.sqrt(new_sigma_sq), T),
                        penalty_kind, penalty_slope)
                    new_cost = state.cost - items[nxt]['r'] + (new_pen - old_pen)

                    new_groups = (state.groups | {gid}) if gid is not None else state.groups

                    new_state = PathState(state.items | {nxt},
                                          new_r, new_mu, new_sigma_sq,
                                          new_cost, nxt, new_groups)

                    if not eliminate_dominated(Lambda[nxt], new_state):
                        Lambda[nxt].append(new_state)

        # 8.5  Stage termination & housekeeping
        if not some_extension:
            break
        L *= 2
    # end outer loop

    return best_items, best_profit, best_cost


###############################################################################
# 9.  Simple built-in demo                                                     #
###############################################################################
"""
A small (10-item) instance taken from the paper’s computational section.
Run this module directly (“python skp_solver.py”) to see the DP in action.
"""


def demo_instance() -> None:
    items = {
        0: {'r': 111, 'mu': 44, 'sigma_sq': 77.44},
        1: {'r': 111, 'mu': 42, 'sigma_sq': 70.56},
        2: {'r': 21,  'mu': 73, 'sigma_sq': 213.16},
        3: {'r': 117, 'mu': 15, 'sigma_sq': 9.0},
        4: {'r': 123, 'mu': 71, 'sigma_sq': 201.64},
        5: {'r': 34,  'mu': 12, 'sigma_sq': 5.76},
        6: {'r': 3,   'mu': 13, 'sigma_sq': 6.76},
        7: {'r': 121, 'mu': 14, 'sigma_sq': 7.84},
        8: {'r': 112, 'mu': 23, 'sigma_sq': 21.16},
        9: {'r': 12,  'mu': 15, 'sigma_sq': 9.0},
    }

    T, S, β = 200.0, 100.0, 1.0
    pen_kind, pen_slope = "linear", 10.0

    sol_items, profit, obj_cost = solve_skp(items, T, S, β,
                                            pen_kind, pen_slope,
                                            item_groups=None,
                                            L_init=10)
    print("Selected items:", sorted(sol_items))
    print(f"Profit  R−f  : {profit:.2f}")
    print(f"Objective C : {obj_cost:.2f}")


###############################################################################
# 10.  Minimal batch driver for JSON instance files                            #
###############################################################################
"""
Utility function used during the authors’ computational study.  Reads the
 standard “normal_instances.json” format and calls the solver once.
"""


def run_batch(file_name: str, instance_no: int = 0):
    with open(file_name, 'r') as fh:
        data = json.load(fh)

    inst = data[instance_no]
    items = {i: {'r': r,
                 'mu': mu,
                 'sigma_sq': sigma ** 2}
             for i, (r, mu, sigma)
             in enumerate(zip(inst['expectedValues'],
                              inst['expectedWeights'],
                              inst['stdWeights']))}

    T = inst['capacity']
    S = 2 * T
    β = 1.0
    slope = inst.get('shortageCost', 10.0)

    start_runtime = time.time() * 1000
    sol, total_profit, _ = solve_skp(items, T, S, β,
                               penalty_kind="linear",
                               penalty_slope=slope,
                               item_groups=None,
                               L_init=10)
    end_runtime = time.time() * 1000
    runtime_millisec = end_runtime - start_runtime

    return np.sort(np.array(list(sol))), total_profit, runtime_millisec

# -----------------------------------------------------------------------------


if __name__ == "__main__":
    # demo_instance()
    results = run_batch('normal_instances_50.json', 5)  # Example for running the first instance in the batch file
    print(f"Optimal Expected Profit: {results[1]:.2f}, Runtime: {results[2]:.4f} milliseconds")
