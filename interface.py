from generalized_dijkstra import dijkstra_solve
from ilp_solver import ilp_solve
from tree_width_solver import TWInterface


def localize(matrix, vector, weights, method):
    if method == "ILP":
        return ilp_localize(matrix, vector, weights)
    if method == "DIJKSTRA":
        return dijkstra_localize(matrix, vector, weights)
    if method == "TREEWIDTH":
        return treewidth_localize(matrix, vector, weights)
    if method == "COMPUTE_TREEWIDTH":
        return compute_treewidth(matrix, vector, weights)
    return None

def ilp_localize(matrix, vector, weights):
    matrix = [list(col) for col in matrix]
    vector = list(vector)
    weights = list(weights)
    solution_size, solution, setup, memory_peak = ilp_solve(matrix, vector, weights)
    return solution_size, solution

def dijkstra_localize(matrix, vector, weights):
    m = [set(col) for col in matrix]
    v = set(vector)
    solution_size, solution, setup, memory_peak = dijkstra_solve(m, v, weights)
    return solution_size, solution

def treewidth_localize(matrix, vector, weights, full_setup=True, memory_limit=2 ** 20):
    matrix = [frozenset(col) for col in matrix]
    vector = frozenset(vector)
    weights = list(weights)
    treewidth_interface = TWInterface(matrix, vector, weights, full_setup, memory_limit)
    solution_size, solution, memory_peak = treewidth_interface.homology_localization(print_status=False)
    return solution_size, solution


def compute_treewidth(matrix, vector, weights, full_setup=True, memory_limit=2 ** 20):
    matrix = [frozenset(col) for col in matrix]
    vector = frozenset(vector)
    weights = list(weights)
    treewidth_interface = TWInterface(matrix, vector, weights, full_setup, memory_limit)
    return treewidth_interface.treewidth_of_spine_graph()
