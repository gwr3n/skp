import json

import numpy as np
from binary_heuristic import heuristic_binary_knapsack_solver
from branch_and_bound import branch_and_bound_solver
import time
import os

# =============================================================================
# Helper functions
# =============================================================================

def compute_running_time(func):
        def wrapper(*args, **kwargs):
            start_time = time.time()
            result = func(*args, **kwargs)
            end_time = time.time()
            print(f"Execution time: {end_time - start_time:.4f} seconds")
            return result
        return wrapper

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

def solve_bab_binary_knapsack_solver(problem_data):
     # Extract individual parameters for use in the algorithm
    id = problem_data['instanceID']
    mu = problem_data['expectedWeights']
    sigma = problem_data['stdWeights']
    v = sigma**2
    r = problem_data['expectedValues']
    C = problem_data['capacity']
    p = problem_data['shortageCost']

    branch_rule = "fractional"  # change as desired
    start_time = time.time() * 1000
    bb_sol = branch_and_bound_solver(r, mu, v, C, p, branch_rule=branch_rule, tol=1e-4, time_limit=600)
    end_time = time.time() * 1000
    bb_sol['solutionTime'] = end_time - start_time

    return bb_sol

def solve_heuristic_binary_knapsack_solver(problem_data):
    # Extract individual parameters for use in the algorithm
    id = problem_data['instanceID']
    mu = problem_data['expectedWeights']
    sigma = problem_data['stdWeights']
    v = sigma**2
    r = problem_data['expectedValues']
    C = problem_data['capacity']
    p = problem_data['shortageCost']

    start_time = time.time() * 1000
    heuristic_result = heuristic_binary_knapsack_solver(r, mu, v, C, p)
    end_time = time.time() * 1000
    heuristic_result['solutionTime'] = end_time - start_time

    return heuristic_result

def rcwt(input):
    """
    Replace all commas in the input string with tab characters.

    Args:
        input_string (str): The input string.

    Returns:
        str: The modified string with commas replaced by tabs.
    """
    if not isinstance(input, str):
        input_string = np.array2string(input, separator='\t', max_line_width=np.inf)
    return input_string

def solve_knapsack_instances_HEUR(problem_data_list, directory=None):
    # This function is used to solve knapsack instances using only the heuristic approach.
    # This is helpful for large instances, as the B&B solver is not efficient for large instances.
    header = (
        "instanceID,expectedValues,expectedWeights,stdWeights,capacity,shortageCost,"
        "HEURoptimalKnapsack,HEURsolutionValue,HEURsolutionTimeMs\n"
    )

    print("Solving knapsack instances...")

    file_path = 'solved_normal_instances_Merzifonluoglu_heuristic.csv' if directory is None else os.path.join(directory, 'solved_normal_instances_Merzifonluoglu_heuristic.csv')

    # Create CSV file and write header
    with open(file_path, 'w') as file:
        file.write(header)
        
    for problem_data in problem_data_list:

        # Solve using Heuristic
        heuristic_solution = solve_heuristic_binary_knapsack_solver(problem_data)

        result = (
             f"{problem_data['instanceID']},{rcwt(problem_data['expectedValues'])},{rcwt(problem_data['expectedWeights'])},"
                f"{rcwt(problem_data['stdWeights'])},{problem_data['capacity']},{problem_data['shortageCost']},"
                f"{rcwt(heuristic_solution['x_binary'])},"
                f"{heuristic_solution['bin_obj']},{heuristic_solution['solutionTime']}\n"
        )

        # Write results to CSV file
        with open(file_path, 'a') as file:
            file.write(''.join(result))
        
        print(".", end="", flush=True)

def solve_knapsack_instances(problem_data_list, directory=None):
    # This function is used to solve knapsack instances using both B&B and heuristic approaches.
    header = (
        "instanceID,expectedValues,expectedWeights,stdWeights,capacity,shortageCost,"
        "B&BoptimalKnapsack,B&BsolutionValue,B&BsolutionTimeMs,B&BexploredNodes,"
        "HEURoptimalKnapsack,HEURsolutionValue,HEURsolutionTimeMs,HEURoptGap\n"
    )

    print("Solving knapsack instances...")

    file_path = 'solved_normal_instances_Merzifonluoglu.csv' if directory is None else os.path.join(directory, 'solved_normal_instances_Merzifonluoglu.csv')

    # Create CSV file and write header
    with open(file_path, 'w') as file:
        file.write(header)
        
    for problem_data in problem_data_list:
        # Solve using Branch and Bound
        bb_solution = solve_bab_binary_knapsack_solver(problem_data)
        # Solve using Heuristic
        heuristic_solution = solve_heuristic_binary_knapsack_solver(problem_data)

        opt_gap = (bb_solution['obj'] - heuristic_solution['bin_obj']) / bb_solution['obj'] * 100 if bb_solution['obj'] != 0 else 0

        result = (
             f"{problem_data['instanceID']},{rcwt(problem_data['expectedValues'])},{rcwt(problem_data['expectedWeights'])},"
                f"{rcwt(problem_data['stdWeights'])},{problem_data['capacity']},{problem_data['shortageCost']},"
                f"{rcwt(bb_solution['x'])},{bb_solution['obj']},{bb_solution['solutionTime']},"
                f"{bb_solution['node_count']},{rcwt(heuristic_solution['x_binary'])},"
                f"{heuristic_solution['bin_obj']},{heuristic_solution['solutionTime']},"
                f"{opt_gap}\n"
        )

        # Write results to CSV file
        with open(file_path, 'a') as file:
            file.write(''.join(result))
        
        print(".", end="", flush=True)

def recursive_solve(directory):
    for root, dirs, files in os.walk(directory):
        # print(f"Processing directory: {root}")
        if 'normal_instances.json' in files:
            print(f"Found normal_instances.json in {root}")
            file_path = os.path.join(root, 'normal_instances.json')
            problem_data_list = load_problem_data(file_path)
            if len(problem_data_list[0]['expectedValues']) < 50:
                solve_knapsack_instances(problem_data_list, root)
            else:    
                solve_knapsack_instances_HEUR(problem_data_list, root)

def solve_normal_instances():
    problem_data_list = load_problem_data('normal_instances.json')
    solve_knapsack_instances_HEUR(problem_data_list)

    # problem_data_list = load_problem_data('normal_instances_500.json')
    # solve_knapsack_instances_HEUR(problem_data_list)

if __name__ == '__main__':
    recursive_solve("~Downloads/normal")