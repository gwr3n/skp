#!/usr/bin/env python3
###############################################################################
#                                                                             #
#  ┌──────────────────────────────────────────────────────────────────────┐     #
#  │  Heuristic for the Static Stochastic Knapsack Problem (SSKP)        │     #
#  │  ───────────────────────────────────────────────────────────────────  │     #
#  │  Reference implementation of the ideas in                           │     #
#  │     Merzifonluoğlu, Geunes & Romeijn (2012)                         │     #
#  │     “The static stochastic knapsack problem with                    │     #
#  │      normally distributed item sizes”, Math. Prog. A, 134:459-489   │     #
#  └──────────────────────────────────────────────────────────────────────┘     #
#                                                                             #
#  Philosophy (literate-programming style)                                    #
#  ---------------------------------------                                     #
#  • The *explanation* precedes the code it motivates.  Each major block is    #
#    introduced by an extensive doc-string or comment so that the file may be #
#    read linearly as prose.                                                   #
#  • All equation numbers, symbols (μ, σ2, ρ, L, …) and section references     #
#    follow exactly the notation of the paper, so that a reader can keep the   #
#    article and this file side-by-side.                                       #
#                                                                             #
#  Scope restrictions implemented                                             #
#  ---------------------------------                                          #
#  The complete algorithm in the paper is quite elaborate; here we include    #
#  only the parts strictly needed for a **heuristic solver** under the        #
#  following assumptions:                                                     #
#                                                                             #
#     (1)  Every item is *truly* stochastic:  σ²ᵢ > 0  (I₀ ≡ ∅).              #
#          ⇒ At most one fractional item in any KKT candidate (§4.2).         #
#     (2)  Salvage value  s  is taken as 0, so  p = p′  (no back-ordering).   #
#     (3)  Item means satisfy  μᵢ > 0  and the capacity constant  C > 0.      #
#     (4)  Capacity itself is deterministic and not for sale/choice.          #
#     (5)  Thus, “multi-fractional” KKT patterns involving deterministic      #
#          items never arise and need not be generated.                       #
#                                                                             #
#  These assumptions match those in our study. Note that this implementation  #
#  has lower complexity than the original KKT heuristic, as it only needs to  #
#  consider a restricted number of cases.                                     #
#                                                                             #
#  High-level structure                                                       #
#  ---------------------                                                      #
#  ┌ Stage I ─ generate_candidate_solutions ┐                                 #
#  │  Exhaustively enumerates **all** continuous-relaxation KKT points        │
#  │  allowed by the above scope (§3).                                        │
#  └────────────────────────────────────────┘                                 #
#  ┌ Stage II ─ heuristic_binary_knapsack_solver ┐                            #
#  │  Rounds every fractional pattern up/down (f→0 / f→1)     (§4.2)          │
#  │  and returns the best discrete solution found.                           │
#  └─────────────────────────────────────────────┘                            #
#                                                                             #
#  External requirements                                                      #
#  ----------------------                                                     #
#      • numpy   – vectorised linear algebra                                  #
#      • scipy   – norm distribution, Brent root-finder                       #
#                                                                             #
#  Author : Roberto Rossi                                                     #
#  First  : 20-Jun-2025                                                       #
#                                                                             #
###############################################################################

import numpy as np
from scipy.stats import norm
from scipy.optimize import brentq
import json

# ════════════════════════════════════════════════════════════════════════════
# Utility mathematics                                                      ══
# ════════════════════════════════════════════════════════════════════════════

def L_func(z):
    """
    Standard normal *loss function*      L(z) = φ(z) − z·(1 − Φ(z))
    (sometimes called the “overshoot function” in inventory theory).

    It appears in Eq. (P) of the article once the overflow expectation
    E[(D−C)⁺] is expressed in closed form for  D ∼ N(μ, σ²).
    """
    phi = norm.pdf(z)
    Phi = norm.cdf(z)
    return phi - z * (1 - Phi)


def attractiveness(i, z, r, mu, v, p):
    """
    Compute ρᵢ(z) – the *attractiveness* of item *i* at service level *z*
    (Eq. 19).  The definition depends on whether the item is stochastic or
    deterministic; in this restricted implementation vᵢ>0 holds for every i,
    but we keep the general formula for completeness.
    """
    Phi = norm.cdf(z)
    if v[i] > 0:                                                # i ∈ I⁺
        return (r[i] - p * mu[i] * (1 - Phi)) / v[i]
    else:                                                       # i ∈ I⁰
        return -np.inf if r[i] < p * mu[i] * (1 - Phi) else np.inf


def solve_fractional_value(z, mu_k, v_k, S_mu, S_v, C):
    """
    Solve for the unique fractional level  f ∈ (0,1)  that satisfies the
    capacity-matching equations (23)–(24) for a pattern of type (21):

           f·μ_k + S_μ + z·√(f·v_k + S_v)  =  C.

    The function returns  (f*, V*)  with  V* = f*·v_k + S_v, or (None,None)
    if no feasible fraction exists (e.g. the equation has no root in [0,1]).
    """
    # Build the scalar equation  F(f) = 0.
    def F(f):
        return f * mu_k + S_mu + z * np.sqrt(f * v_k + S_v) - C

    # Check for sign change over [0,1] – otherwise Brent would fail.
    try:
        if F(0.0) * F(1.0) > 0:
            return None, None
        f_star = brentq(F, 0.0, 1.0)             # root inside the segment
    except Exception:                            # numerical mishap
        return None, None

    V_star = f_star * v_k + S_v
    return f_star, V_star


def compute_continuous_obj(r, x, p, z, V):
    """
    Continuous-relaxation objective  rᵀx − p·L(z)·√V   (Eq. (P) with s=0).
    """
    return np.dot(r, x) - p * L_func(z) * np.sqrt(V)


def binary_objective(x, r, mu, v, C, p):
    """
    Expected profit of *binary* vector x according to Sect. 2.2:

        obj(x) = rᵀx  −  p·E[(D(x) − C)⁺].

    The overflow expectation uses the closed form with L if V>0,
    else reverts to the deterministic expression.
    """
    m = np.dot(mu, x)
    V = np.dot(v, x)
    if V > 1e-10:
        penalty = p * np.sqrt(V) * L_func((C - m) / np.sqrt(V))
    else:
        penalty = p * max(0, m - C)
    return np.dot(r, x) - penalty


# ════════════════════════════════════════════════════════════════════════════
# Stage I – exhaustive enumeration of KKT-consistent continuous solutions  ══
# ════════════════════════════════════════════════════════════════════════════

def generate_candidate_solutions(r, mu, v, C, p):
    """
    Enumerate **all** KKT candidates described in Sect. 3.3.

    For each z ∈ 𝒵 (= 𝒵⁺ ∪ 𝒵⁰ plus mid-points) we order items by ρ(z),
    slide the threshold across ties (pattern 21) and, whenever necessary,
    solve for the single fractional variable via Brent.

    The resulting list of dictionaries is subsequently consumed by the
    rounding heuristic (Stage II).

    Returned fields
    --------------
    x                 full selection vector (may include one fractional)
    z, V              service-level parameter and total variance
    threshold_index   the (unique) fractional item index
    fractional_value  its value f ∈ (0,1)
    ordering          permutation of items by non-decreasing ρ
    rho               array of ρᵢ(z) values
    cont_obj          continuous objective value of the candidate
    """
    n       = len(r)
    idx_all = np.arange(n)
    I_plus  = idx_all[v > 0]                     # stochastic items
    I_zero  = idx_all[v == 0]                    # deterministic (unused here)

    # ──────────────────────────────────────────────────────────────
    # 1. Build the *service-level* set 𝒵 = 𝒵⁺ ∪ 𝒵⁰  (plus mid points)
    # ──────────────────────────────────────────────────────────────
    candZ = set()

    # 𝒵⁺ : solve z(i,j) for all stochastic pairs with non-degenerate denom.
    for i in I_plus:
        for j in I_plus:
            if i >= j:
                continue
            denom = p * (mu[i] * v[j] - mu[j] * v[i])
            if abs(denom) < 1e-12:
                continue
            ratio = (r[i] * v[j] - r[j] * v[i]) / denom
            if 0.0 < ratio < 1.0:
                candZ.add(norm.ppf(1.0 - ratio))

    # 𝒵⁰ : each deterministic item contributes a single z = zᵢ
    for i in I_zero:
        ratio = r[i] / (p * mu[i])
        if 0.0 < ratio < 1.0:
            candZ.add(norm.ppf(1.0 - ratio))

    candZ.add(0.0)                              # always safe / often optimal

    z_vals = np.array(sorted(candZ))
    if len(z_vals) >= 2:                         # add mid-points of intervals
        mids   = (z_vals[:-1] + z_vals[1:]) / 2.0
        z_vals = np.sort(np.concatenate((z_vals, mids)))

    # ──────────────────────────────────────────────────────────────
    # 2. For each z explore every threshold position (tied blocks, etc.)
    # ──────────────────────────────────────────────────────────────
    eps_tie  = 1e-6                              # tolerance for equal ρ
    candidates = []

    for z in z_vals:
        rho   = np.array([attractiveness(i, z, r, mu, v, p) for i in range(n)])
        order = np.argsort(rho)                  # σ(1)…σ(n) non-decreasing

        pos = 0
        while pos < n:
            # Identify contiguous block of equal-ρ items [start … end].
            start = pos
            while (pos + 1 < n and
                   abs(rho[order[pos+1]] - rho[order[start]]) <= eps_tie):
                pos += 1
            end   = pos
            block = order[start:end+1]
            t     = len(block)

            # Enumerate every *suffix* of that block (size s = 0…t)
            for s in range(t + 1):
                sel_suffix = block[t-s:] if s else np.empty(0, dtype=int)
                full_selected = np.concatenate((sel_suffix, order[end+1:]))

                # Determine fractional item k and pattern-specific sums
                if s == t:                       # threshold before the block
                    continue
                idx_frac = block[0] if s == 0 else sel_suffix[0]

                x = np.zeros(n)
                x[full_selected] = 1.0

                S_mu = mu[full_selected].sum() if full_selected.size else 0.0
                S_v  = v [full_selected].sum() if full_selected.size else 0.0
                mu_k, v_k = mu[idx_frac], v[idx_frac]

                # Solve capacity equation for fractional level f*
                f_star, V_star = solve_fractional_value(z, mu_k, v_k,
                                                         S_mu, S_v, C)
                if f_star is None:
                    continue                    # infeasible pattern

                x[idx_frac] = f_star
                cont_obj = compute_continuous_obj(r, x, p, z, V_star)

                # Collect everything for Stage II
                candidates.append(dict(
                    x=x, z=z, V=V_star,
                    threshold_index=idx_frac,
                    fractional_value=f_star,
                    ordering=order,
                    rho=rho,
                    cont_obj=cont_obj
                ))
            pos += 1                            # move to next block

    return candidates


# ════════════════════════════════════════════════════════════════════════════
# Stage II – integer heuristic (rounding of Stage I candidates)             ══
# ════════════════════════════════════════════════════════════════════════════

def heuristic_binary_knapsack_solver(r, mu, v, C, p):
    """
    Implement the rounding heuristic of Sect. 4.2.

    For every continuous candidate produced by Stage I we create up to two
    integer vectors – one with the fractional component set to 0, one with it
    set to 1 – then compute the exact expected profit (binary_objective).  The
    best of all rounded patterns is returned.
    """
    candidate_list = generate_candidate_solutions(r, mu, v, C, p)

    best_bin_obj   = -np.inf
    best_x_binary  = None
    best_details   = None

    for cand in candidate_list:
        x_cont  = cand['x']
        thr_idx = cand['threshold_index']

        # Build list of integer variants to test
        if np.all(np.isin(x_cont, (0.0, 1.0))):     # already integral
            x_variants = [x_cont.copy()]
        else:
            x0 = x_cont.copy(); x0[thr_idx] = 0.0
            x1 = x_cont.copy(); x1[thr_idx] = 1.0
            x_variants = [x0, x1]

        # Evaluate each integer variant
        for x_bin in x_variants:
            obj = binary_objective(x_bin, r, mu, v, C, p)
            if obj > best_bin_obj:
                best_bin_obj  = obj
                best_x_binary = x_bin.copy()
                best_details  = dict(base_candidate=cand,
                                     rounding=('none' if len(x_variants)==1
                                               else '0/1'),
                                     evaluated_obj=obj)

    if best_x_binary is None:
        print("Heuristic failed to produce any integer solution.")
        return None

    return dict(x_binary=best_x_binary,
                bin_obj=best_bin_obj,
                candidate_details=best_details)


# ════════════════════════════════════════════════════════════════════════════
# I/O helpers – load JSON benchmark files                                   ══
# ════════════════════════════════════════════════════════════════════════════

def load_problem_data(file_path):
    """
    Read a list of problem instances from *file_path* (JSON) and convert
    numeric lists to NumPy arrays for vectorised computation.
    """
    with open(file_path, 'r') as file:
        problem_data_list = json.load(file)

    for problem in problem_data_list:
        for key in ['expectedWeights', 'stdWeights', 'expectedValues']:
            if key in problem:
                problem[key] = np.array(problem[key])

    return problem_data_list


# ════════════════════════════════════════════════════════════════════════════
# Example driver and pretty printing                                       ══
# ════════════════════════════════════════════════════════════════════════════

def get_problem_data():
    """
    A tiny hard-coded demo instance (25 items) useful when no JSON file is
    available.  Values chosen arbitrarily.
    """
    return {
        'expectedWeights': np.array([...]),        # (values trimmed here)
        'stdWeights'     : np.array([...]),
        'expectedValues' : np.array([...]),
        'capacity'       : 116.10846413274393,
        'shortageCost'   : 10,
        'instanceID'     : 'demo-instance'
    }


def solve_instance(problem_data, logging=False):
    """
    Solve *one* instance and optionally print every intermediate artefact.
    """
    id   = problem_data['instanceID']
    mu   = problem_data['expectedWeights']
    sigma= problem_data['stdWeights']
    v    = sigma**2
    r    = problem_data['expectedValues']
    C    = problem_data['capacity']
    p    = problem_data['shortageCost']

    if logging:
        print(f"Parameters for instance {id}:")
        print(f" μ  = {mu}\n σ² = {v}\n r  = {r}\n C  = {C}\n p  = {p}")
        print("="*60)

    # Stage I – just for curiosity / debugging
    cand_list = generate_candidate_solutions(r, mu, v, C, p)
    if logging:
        print(f"Generated {len(cand_list)} continuous candidates.")
        for k, cand in enumerate(cand_list, 1):
            print(f"#{k:3d}: z={cand['z']:+6.3f}  k={cand['threshold_index']:2d} "
                  f"f={cand['fractional_value']:.3f}  obj={cand['cont_obj']:.3f}")
        print("="*60)

    # Stage II – rounding heuristic
    result = heuristic_binary_knapsack_solver(r, mu, v, C, p)
    if result is None:
        print(f"{id}: *no* integer solution found!"); return

    if logging:
        print(f"Best binary solution for {id}:")
        print(" x =", result['x_binary'])
        print(f" objective = {result['bin_obj']:.4f}")
        print(" Derived from continuous candidate:")
        print(result['candidate_details'])
    else:
        print(f"{id[:10]}: {result['bin_obj']:.3f}\t{result['x_binary']}")


def solve_batch():
    """
    Iterate over all instances in the JSON file *normal_instances.json*
    and solve each with the heuristic.
    """
    problem_set = load_problem_data('normal_instances_25.json')
    for problem in problem_set:
        solve_instance(problem)


# ════════════════════════════════════════════════════════════════════════════
# Script entry point                                                       ══
# ════════════════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    solve_batch()
