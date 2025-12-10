#!/usr/bin/env python3
"""
Branch-and-Bound implementation for the static stochastic knapsack problem
with normally distributed item sizes. This code uses the continuous-relaxation
algorithm (Section 3 of Merzifonluoğlu et al. 2012) as a bounding procedure and
implements an exact B&B algorithm for the binary problem (Section 4.1).

Two branching rules are implemented:
  (a) "fractional" - branch on the free variable with value closest to 0.5.
  (b) "attractiveness" - branch on the free variable having the highest attractiveness measure.

A switch (branch_rule) lets the user choose the preferred rule.

Open source routines from SciPy (brentq, norm) are used.

Author :  Roberto Rossi
Date   :  20-Jun-2025
"""

import json
import numpy as np
from scipy.stats import norm
from scipy.optimize import brentq
import time

# =============================================================================
# Helper functions
# =============================================================================

def L_func(z):
    """
    Standard normal loss function, as given in the paper:
      L(z) = φ(z) - z · (1 - Φ(z))
    """
    phi = norm.pdf(z)
    Phi = norm.cdf(z)
    return phi - z*(1 - Phi)

def attractiveness(i, z, r, mu, v, p):
    """
    Compute the measure of attractiveness ρ_i(z) for item i.
    For items with positive variance (i in I+): 
        ρ_i(z) = [ r_i - p*μ_i*(1 - Φ(z)) ] / v_i.
    For deterministic items (v[i]==0): returns -∞ (if item is unattractive) or +∞.
    """
    Phi = norm.cdf(z)
    if v[i] > 0:
        return (r[i] - p * mu[i] * (1 - Phi)) / v[i]
    else:
        # Deterministic item.
        if r[i] < p * mu[i] * (1 - Phi):
            return -np.inf
        else:
            return np.inf

def solve_fractional_value(z, mu_k, v_k, S_mu_free, S_v_free, S_mu_fixed, S_v_fixed, C):
    """
    For a candidate solution of the form for free items:
       x[i] = 0 for items before the threshold,
       x[k] = f (fractional, in (0,1)) for the threshold free item,
       x[i] = 1 for free items after the threshold,
    and given that fixed items already contribute S_mu_fixed to the sum and S_v_fixed to V,
    the binding constraint (derived from Eq. (1)) is:
    
         f * mu_k + S_mu_fixed + S_mu_free + z * sqrt( f*v_k + S_v_fixed + S_v_free ) = C.
    
    That is, we wish to find f in [0,1] such that:
         F(f) = f * mu_k + (S_mu_fixed + S_mu_free) + z * sqrt( f*v_k + (S_v_fixed + S_v_free) ) - C = 0.
         
    Returns:
       f, V_total   (where V_total = f*v_k + S_v_fixed + S_v_free)
       or (None, None) if no solution is found.
    """
    def F(f):
        denom = f * v_k + (S_v_fixed + S_v_free)
        return f * mu_k + (S_mu_fixed + S_mu_free) + z * np.sqrt(denom) - C

    # Check endpoints
    try:
        F0 = F(0.0)
        F1 = F(1.0)
    except Exception:
        return None, None

    if F0 * F1 > 0:
        # No sign change => no root in [0,1]
        return None, None

    if F0 == 0:
        # f=0 is a root
        f_star = 0.0
    elif F1 == 0:
        # f=1 is a root
        f_star = 1.0
    else:
        try:
            f_star = brentq(F, 0.0, 1.0)
        except Exception:
            return None, None

    V_total = f_star * v_k + (S_v_fixed + S_v_free)
    return f_star, V_total

def compute_objective(r, x, p, z, V):
    """
    Compute the objective value:
         Obj = r^T x - p * L(z)* sqrt(V)
    """
    return np.dot(r, x) - p * L_func(z) * np.sqrt(V)

# =============================================================================
# Main algorithm for continuous relaxation (with fixed variables)
# =============================================================================

def continuous_knapsack_solver(r, mu, v, C, p, fixed=None, tol=1e-5):
    """
    Solve the continuous relaxation (R) as in Section 3 of the paper.
    
    Input parameters:
      r      : np.array of expected revenues (shape (n,))
      mu     : np.array of expected sizes (shape (n,))
      v      : np.array of variances (shape (n,)); v=0 means deterministic.
      C      : capacity (scalar)
      p      : effective penalty cost (with salvage value already subtracted).
      fixed  : (optional) dict mapping indices (0-indexed) to fixed values {0,1}. 
               These variables will be fixed to the corresponding value.
    
    Returns a dict with:
      'x'  : best candidate continuous selection vector (values in [0,1]) of length n
      'z'  : corresponding dummy variable z
      'V'  : corresponding V value (variance sum)
      'obj': objective function value
      'details': additional details.
    
    Procedure:
      • Generate candidate z values from free items (those not fixed).
      • For each candidate z, rank the free items by attractiveness.
      • Then, for each candidate threshold k (in free order),
         solve for the fractional value f of the threshold free item (using the binding constraints).
      • Merge the solution with the fixed items and select the candidate with best objective.
    """
    n = len(r)
    all_indices = set(range(n))
    if fixed is None:
        fixed = {}
    fixed_indices = set(fixed.keys())
    free_indices = np.array(sorted(list(all_indices - fixed_indices)))
    
    # Precompute fixed contributions.
    S_mu_fixed = sum(mu[i] * fixed[i] for i in fixed)
    S_v_fixed  = sum(v[i] * fixed[i] for i in fixed)

    # If no free items exist, then x is completely fixed; we must simply check feasibility.
    if free_indices.size == 0:
        x_candidate = np.array([fixed[i] for i in range(n)])
        total_mu = np.dot(mu, x_candidate)
        total_v  = np.dot(v, x_candidate)
        # We can compute the implied z from Eq. (1): z = (C - total_mu) * sqrt(total_v)
        # (if total_v > 0; otherwise, one can set z arbitrarily, say zero)
        if total_v > 0:
            z_val = (C - total_mu) / np.sqrt(total_v)
        else:
            z_val = 0.0
        best_solution = x_candidate
        best_obj = compute_objective(r, x_candidate, p, z_val, total_v)
        best_details = {
            'z': z_val,
            'V': total_v,
            'message': 'All variables fixed'
        }
        result = {
            'x': best_solution,
            'z': best_details['z'],
            'V': best_details['V'],
            'obj': best_obj,
            'details': best_details
        }
        return result

    # In the remainder we work only with free indices.
    # For candidate z generation, define I+_free and I0_free:
    free_arr = free_indices  # free indices as a numpy array
    I_plus_free = free_arr[ v[free_arr] > 0 ]
    I_zero_free = free_arr[ v[free_arr] == 0 ]
    
    candidate_z = set()
    
    # (a) Candidate z from pairs of free items (I⁺). For each pair (i,j) with i < j:
    for i in I_plus_free:
        for j in I_plus_free:
            if i < j:
                denom = p * (mu[i]*v[j] - mu[j]*v[i])
                if np.abs(denom) < 1e-10:
                    continue
                ratio = (r[i]*v[j] - r[j]*v[i]) / denom
                if 0.0 < ratio < 1.0:
                    try:
                        z_candidate = norm.ppf(1 - ratio)
                        candidate_z.add(z_candidate)
                    except Exception:
                        pass
    # (b) Candidate z from deterministic free items (I0_free):
    for i in I_zero_free:
        if mu[i] <= 0:
            continue
        ratio = r[i] / (p * mu[i])
        if 0.0 < ratio < 1.0:
            z_candidate = norm.ppf(1 - ratio)
            candidate_z.add(z_candidate)
    
    # ----- NEW: Explicitly add the candidate corresponding to "all free items = 1" -----
    S_mu_free_all = np.sum(mu[free_indices])
    S_v_free_all = np.sum(v[free_indices])
    if S_v_fixed + S_v_free_all > 0:
        z_all1 = (C - (S_mu_fixed + S_mu_free_all)) / np.sqrt(S_v_fixed + S_v_free_all)
    else:
        # In the rare case where all free items are deterministic, set z_all arbitrarily.
        z_all1 = 0.0
    candidate_z.add(z_all1)
    # -----------------------------------------------------------------------------------

    # ----- NEW: Explicitly add the candidate corresponding to "all free items = 0" -----
    # When all free items are set to 0, then the free contributions are zero.
    S_mu_free_zero = 0.0
    S_v_free_zero  = 0.0
    if S_v_fixed + S_v_free_zero > 0:
        z_all0 = (C - (S_mu_fixed + S_mu_free_zero)) / np.sqrt(S_v_fixed + S_v_free_zero)
    else:
        z_all0 = 0.0
    candidate_z.add(z_all0)
    # -----------------------------------------------------------------------------------
    
    # (c) If no candidate z generated, add a default interval.
    if candidate_z:
        z_lower = min(candidate_z)
        z_upper = max(candidate_z)
    else:
        z_lower = -3.0
        z_upper = 3.0
        candidate_z.add(z_lower)
        candidate_z.add(z_upper)
    
    candidate_z = np.array(sorted(candidate_z))
    
    best_obj = -np.inf
    best_solution = None
    best_details = {}
    
    # Loop over each candidate z value.
    for z_val in candidate_z:
        # Compute attractiveness ρ_i(z) for each free item.
        rho = {}
        for i in free_indices:
            rho[i] = attractiveness(i, z_val, r, mu, v, p)
        # Order free indices in increasing order of rho.
        free_order = np.array(sorted(rho.keys(), key=lambda i: rho[i]))
        
        # Now consider candidate threshold indices (in free_order).
        for k in range(len(free_order)):
            # Construct a candidate x vector (of full dimension n).
            x_candidate = np.zeros(n)
            # For fixed indices, set as given.
            for i in fixed:
                x_candidate[i] = fixed[i]
            # Among free items, the candidate solution is assumed to have:
            #   free-order indices before the threshold: x = 0,
            #   threshold free item: fractional f (to be determined),
            #   free-order indices after threshold: x = 1.
            threshold_index = free_order[k]  # index of free item to be fractional
            free_full = free_order[k+1:]     # free items set fully to 1.
            
            # Set free items in free_full.
            for i in free_full:
                x_candidate[i] = 1.0
            # Compute sums from free items that are fully selected.
            S_mu_free = np.sum(mu[free_full]) if free_full.size > 0 else 0.0
            S_v_free  = np.sum(v[free_full]) if free_full.size > 0 else 0.0

            # Parameters for the threshold item.
            mu_k = mu[threshold_index]
            v_k  = v[threshold_index]
            
            # Solve for fractional value f using modified equation:
            #   f * mu_k + (S_mu_fixed + S_mu_free) + z_val * sqrt( f*v_k + (S_v_fixed + S_v_free) ) = C.
            f_star, V_total = solve_fractional_value(z_val, mu_k, v_k, S_mu_free, S_v_free, S_mu_fixed, S_v_fixed, C)
            if f_star is None:
                continue
            # Enforce candidate structure for the threshold free item.
            x_candidate[threshold_index] = f_star
            
            # (Optional) Check feasibility of f_star in [0,1]
            if f_star < 0 or f_star > 1:
                continue

            # Compute candidate objective value.
            obj_val = compute_objective(r, x_candidate, p, z_val, V_total)

            # Check boundary conditions, in case the fractional value is very close to 0 or 1; if so, we can set the item to 0 or 1.
            if f_star <= 0.5:
                y = x_candidate.copy()
                y[threshold_index] = 0.0
                total_mu_boundary = np.dot(mu, y)
                total_v_boundary = np.dot(v, y)
                if total_v_boundary > 0:
                    z_val_boundary = (C - total_mu_boundary) / np.sqrt(total_v_boundary)
                else:
                    z_val_boundary = 0.0

                obj_val_boundary = compute_objective(r, y, p, z_val_boundary, total_v_boundary)

                if obj_val_boundary > obj_val + tol:
                    obj_val = obj_val_boundary
                    x_candidate[threshold_index] = 0.0
                    V_total = total_v_boundary
                    z_val = z_val_boundary
            elif f_star >= 1-0.5:   
                y = x_candidate.copy()
                y[threshold_index] = 1.0
                total_mu_boundary = np.dot(mu, y)
                total_v_boundary = np.dot(v, y)
                if total_v_boundary > 0:
                    z_val_boundary = (C - total_mu_boundary) / np.sqrt(total_v_boundary)
                else:
                    z_val_boundary = 0.0

                obj_val_boundary = compute_objective(r, y, p, z_val_boundary, total_v_boundary)

                if obj_val_boundary > obj_val + tol:
                    obj_val = obj_val_boundary
                    x_candidate[threshold_index] = 1.0
                    V_total = total_v_boundary
                    z_val = z_val_boundary
            
            if obj_val > best_obj:
                best_obj = obj_val
                best_solution = x_candidate.copy()
                best_details = {
                    'z': z_val,
                    'V': V_total,
                    'fractional_index': threshold_index,
                    'fractional_value': f_star,
                    'free_order': free_order,
                    'rho': rho
                }
    
    if best_solution is None:
        # print("No feasible candidate solution found!")
        return None
    
    result = {
        'x': best_solution,
        'z': best_details['z'],
        'V': best_details['V'],
        'obj': best_obj,
        'details': best_details
    }
    return result

# =============================================================================
# Branch-and-Bound algorithm (exact) for the binary (0-1) SSKP
# =============================================================================

class BBNode:
    def __init__(self, fixed):
        # fixed is a dictionary mapping item indices to {0,1}
        self.fixed = fixed.copy()  # decisions that have been fixed in this node

def is_integer(x, tol=1e-5):
    """Return True if all entries of x are integer (0 or 1) within tolerance."""
    return np.all(np.abs(x - np.round(x)) <= tol)

def branch_and_bound_solver(r, mu, v, C, p, branch_rule="fractional", tol=1e-5, time_limit=600, logging=False):
    """
    Solve the full binary stochastic knapsack problem using branch-and-bound.
    
    r, mu, v, C, p: problem parameters.
    branch_rule: either "fractional" or "attractiveness" (see description above).
    tol: tolerance for integrality.
    time_limit: maximum allowed seconds.
    
    Returns a dictionary with the best binary solution found.
    """
    n = len(r)
    global_best = {'obj': -np.inf, 'x': None, 'node_fixed': None}
    start_time = time.time()
    
    # Node list (use a simple list as a stack for depth-first search)
    node_list = [BBNode(fixed={})]
    
    # Keep track of the total node count, etc.
    node_count = 0
    
    while node_list:
        if (time.time() - start_time) > time_limit:
            print("Time limit reached.")
            break
        node = node_list.pop()  # depth-first
        node_count += 1
        
        # Solve the continuous relaxation at this node.
        sol = continuous_knapsack_solver(r, mu, v, C, p, node.fixed, tol=tol)
        if sol is None:
            # Node is infeasible.
            continue
        # If the continuous relaxation objective is not better than current best, prune.
        if sol['obj'] <= global_best['obj'] + tol:
            continue
        
        x_relax = sol['x']
        # Check if solution is binary (within tolerance)
        if is_integer(x_relax, tol=tol):
            # Feasible integer solution; update global best if improved.
            obj_int = sol['obj']
            if obj_int > global_best['obj']:
                global_best = {'obj': obj_int, 'x': np.round(x_relax), 'node_fixed': node.fixed.copy()}
                if logging:
                    print(f"Found new incumbent: obj = {obj_int:.4f} at node with fixed {node.fixed}")
                # print("x =", x_relax)
                # #### CHECKS #####
                # fixed_mu_sum = sum(mu[i]*node.fixed[i] for i in node.fixed)
                # fixed_v_sum  = sum(v[i]*node.fixed[i] for i in node.fixed)
                # fixed_obj    = sum(r[i]*node.fixed[i] for i in node.fixed)
                # V_total = fixed_v_sum
                # obj = fixed_obj - p * np.sqrt(V_total)* L_func((C - fixed_mu_sum)/np.sqrt(V_total))  # choose z=0 arbitrarily.
                # print("obj =", obj)
                # #### CHECKS #####
            continue
        
        # Otherwise, we need to branch.
        # Determine free (not fixed) indices:
        free_indices = [i for i in range(n) if i not in node.fixed]
        
        # Select branching variable according to chosen rule:
        branch_var = None
        
        if branch_rule == "fractional":
            # Branch on the variable whose value is most fractional (closest to 0.5)
            best_score = np.inf
            for i in free_indices:
                xi = x_relax[i]
                # Only consider if not nearly integer.
                if abs(xi - round(xi)) > tol:
                    score = abs(xi - 0.5)
                    if score < best_score:
                        best_score = score
                        branch_var = i
        elif branch_rule == "attractiveness":
            # Use the attractiveness measure computed at the node with the obtained z value.
            z_val = sol['z']
            best_attr = -np.inf
            for i in free_indices:
                attr = attractiveness(i, z_val, r, mu, v, p)
                if attr > best_attr:
                    best_attr = attr
                    branch_var = i
        else:
            raise ValueError("branch_rule must be either 'fractional' or 'attractiveness'")
        
        if branch_var is None:
            # Should not happen
            continue
        
        # Branch: try fixing branch_var to 1 and 0.
        for fix_val in [1, 0]:
            new_fixed = node.fixed.copy()
            new_fixed[branch_var] = fix_val
            node_list.append(BBNode(new_fixed))
    if logging:    
        print(f"Total nodes processed: {node_count}")
    if global_best['x'] is None:
        print("No feasible integer solution found within the time limit.")
    else:
        global_best['node_count'] = node_count
    return global_best

# =============================================================================
# Batch processing function to load problem data and solve instances
# =============================================================================

def load_problem_data(file_path):
    """
    Load problem data from a JSON file and convert lists to numpy arrays.

    Args:
        file_path (str): Path to the JSON file containing problem data.

    Returns:
        list: A list of problem data dictionaries with numpy arrays for relevant keys.
    """
    with open(file_path, 'r') as file:
        problem_data_list = json.load(file)
    
    # Convert lists within problem_data_list to np.array
    for problem_data in problem_data_list:
        for key in ['expectedWeights', 'stdWeights', 'expectedValues']:
            if key in problem_data:
                problem_data[key] = np.array(problem_data[key])
    
    return problem_data_list

# =============================================================================
# Example usage
# =============================================================================

def get_problem_data():
    """
    Returns a dictionary containing the problem parameters for the knapsack problem.
    """
    return {
        'expectedWeights': np.array([48.17338524597328, 38.85191859076709, 53.96883030177949, 56.36756573697871, 9.287916575997754, 88.6193656711439, 79.81676010179477, 24.544027687087134, 96.5133509216311, 86.44159308328558, 73.98955375860147, 56.76435609836739, 48.78789772998606, 26.560437820519102, 41.21652156092145, 18.941536469813343, 69.4581366049341, 4.478076281128421, 29.0470924986981, 34.75687407665649, 76.10648914709449, 67.71076535499637, 79.81465808172918, 12.80047508643447, 54.175520973864096]),
        'stdWeights': np.array([4.817338524597328, 3.885191859076709, 5.396883030177949, 5.636756573697871, 0.9287916575997754, 8.86193656711439, 7.981676010179477, 2.4544027687087134, 9.65133509216311, 8.644159308328558, 7.398955375860147, 5.676435609836739, 4.878789772998606, 2.6560437820519104, 4.121652156092145, 1.8941536469813343, 6.94581366049341, 0.4478076281128421, 2.90470924986981, 3.475687407665649, 7.610648914709449, 6.771076535499637, 7.981465808172918, 1.280047508643447, 5.41755209738641]),
        'expectedValues': np.array([27.374614896234103, 78.4130955596743, 64.65813191442086, 23.266107454511886, 67.4093146834377, 33.59693259866023, 92.16485054308757, 37.79821214150361, 8.236141621162524, 42.6355367987863, 65.54323036758973, 36.508548229182615, 15.968166079227474, 38.83708767781832, 39.062391091831344, 75.03159125698056, 77.3931165367291, 82.83151230238252, 17.34237700822168, 40.044541913379156, 49.92889765445393, 85.56454225080667, 27.426648865161223, 17.446202319304014, 94.46271171962006]),
        'capacity': 116.10846413274393,
        'shortageCost': 10
    }

def solve_instance(problem_data, logging = False):
    # Extract individual parameters for use in the algorithm
    id = problem_data['instanceID']
    mu = problem_data['expectedWeights']
    sigma = problem_data['stdWeights']
    v = sigma**2
    r = problem_data['expectedValues']
    C = problem_data['capacity']
    p = problem_data['shortageCost']
    
    if logging:
        print("=== SSKP Problem Data ===")
        print("mu =", mu)
        print("v =", v)
        print("r =", r)
        print("Capacity C =", C)
        print("Penalty cost p =", p)
        print("="*50)

    # x = [2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,0]
    # fixed = {i: val for i, val in enumerate(x) if val in [0, 1]}
    # print(continuous_knapsack_solver_node(r, mu, v, C, p, fixed)["obj"])
    
    # First, solve the continuous relaxation (for debugging/demo)
    cont_sol = continuous_knapsack_solver(r, mu, v, C, p)
    if logging:
        if cont_sol is not None:
            print("Continuous relaxation solution (full problem):")
            print("x =", cont_sol['x'])
            print("z =", cont_sol['z'])
            print("V (free variance) =", cont_sol['V'])
            print("Objective =", cont_sol['obj'])
            print("Details:", cont_sol['details'])
        else:
            print("No continuous solution found!")
    
    # Now, call the branch-and-bound solver.
    # Choose branch_rule = "fractional" or "attractiveness"
    branch_rule = "fractional"  # change as desired
    if logging:
        print("\nRunning branch-and-bound (branching rule =", branch_rule, ")...")
    bb_sol = branch_and_bound_solver(r, mu, v, C, p, branch_rule=branch_rule, tol=1e-3, time_limit=600, logging=logging)
    if logging:
        if bb_sol['x'] is not None:
            print("\nBest integer solution found:")
            print("x =", bb_sol['x'])
            print("Objective =", bb_sol['obj'])
            print("Fixed decisions at node =", bb_sol['node_fixed'])
        else:
            print("No integer solution was found.")
    else:
        print(f"{id[:10]}: {bb_sol['obj']}\t{bb_sol['x']}")

def solve_batch():

    problem_data_list = load_problem_data('normal_instances_25.json')
    
    for problem_data in problem_data_list:
        solve_instance(problem_data)

if __name__ == '__main__':
    solve_batch()
    

        