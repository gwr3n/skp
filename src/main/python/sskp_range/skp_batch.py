import os
import json
import numpy as np
from range import run_batch
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

def solve_knapsack_instances(problem_data_list, json_file_path, directory=None):
    # This function is used to solve knapsack instances using both B&B and heuristic approaches.
    header = (
        "instanceID,expectedValues,expectedWeights,stdWeights,capacity,shortageCost,"
        "optimalKnapsack,solutionValue,solutionTimeMs\n"
    )

    print("Solving knapsack instances...")

    file_path = 'solved_normal_instances_Range.csv' if directory is None else os.path.join(directory, 'solved_normal_instances_Range.csv')

    # Create CSV file and write header
    with open(file_path, 'w') as file:
        file.write(header)
        
    for idx, problem_data in enumerate(problem_data_list):
        # Solve using Shortest Path SSKP Solver
        bb_solution = run_batch(json_file_path, idx)

        result = (
             f"{problem_data['instanceID']},{rcwt(problem_data['expectedValues'])},{rcwt(problem_data['expectedWeights'])},"
                f"{rcwt(problem_data['stdWeights'])},{problem_data['capacity']},{problem_data['shortageCost']},"
                f"{rcwt(bb_solution[0])},{bb_solution[1]},{bb_solution[2]}\n"
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
                solve_knapsack_instances(problem_data_list, file_path, root)

if __name__ == '__main__':
    recursive_solve("/Users/gwren/Downloads/normal")