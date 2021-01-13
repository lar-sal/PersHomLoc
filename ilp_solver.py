import numpy as np
from gurobipy import Model, GRB, quicksum

import time

def to_row_matrix(matrix, vector):
    nr_of_rows = max([max([max(col) for col in matrix if col != set([])]), max(list(vector))])+1
    matrix_out = [[j for j in range(len(matrix)) if i in matrix[j]] for i in range(nr_of_rows)]
    vector_out = vector_to_zero_one(vector, nr_of_rows)
    return matrix_out, vector_out

def vector_to_zero_one(vector, max_value):
    vector_out = [0 for i in range(max_value+1)]
    for entry in vector:
        vector_out[entry] = 1
    return vector_out

def ilp_solve(matrix, vector, weights):
    if len(matrix) == 0:
        return 0.0, [], "NaN", "NaN"
    start_time = time.time()
    row_matrix, target = to_row_matrix(matrix, vector)
    b = list(range(len(matrix)))
    v = target
    c = weights
    mdl, xVar, yVar = generate_row_model(row_matrix, b, v, c)
    setup = time.time()-start_time
    mdl.optimize()
    solution = [int(xVar[x].getAttr(GRB.Attr.X)+0.2) for x in xVar]
    #print(solution)
    solution = [x for x in range(len(solution)) if solution[x] ==1]
    solution_size = sum([c[x] for x in solution])
    memory_peak = "NaN"
    return solution_size, solution, setup, memory_peak


# Code for generating a model
def generate_row_model(row_matrix, b, v, c):
    """
    rows_matrix: The matrix, represented as a list of rows
    b: List of variables in b
    v: The target vector
    c: The weights of each column
    """

    # Initializes a model
    mdl = Model('persistenthomologylocalization')

    # Add matrix variables
    x = mdl.addVars(list(b), vtype=GRB.BINARY)

    # Add the dummy variable needed to do arithmetics mod 2
    y = mdl.addVars(list(range(len(row_matrix))), vtype=GRB.INTEGER)

    # Set model to minimization
    mdl.modelSense = GRB.MINIMIZE

    # Set objective function to minimize
    mdl.setObjective(quicksum(x[j] * c[j] for j in b))

    # Set the constrains
    for i in list(range(len(row_matrix))):
        mdl.addConstr(quicksum(x[j] for j in row_matrix[i]) + v[i] == y[i] * 2)
    return mdl, x, y