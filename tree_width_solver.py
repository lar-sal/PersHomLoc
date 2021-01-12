# Simplices: A list of the d+1 simplices of the complex, where each simplex is represented as a frozenset containing its vertices.
# Weights: A dictionary with d-simplices (represented as frozensets) as keys and weights of simplex as value

from treewidth_algorithms.per_hom_loc_sg_rep import per_hom_loc_sg_rep
from treewidth_algorithms.tools.chain_complex import ChainComplex


class TWInterface(object):

    def __init__(self, columns, rows, column_weights, cycle, full_setup=True, memory_limit=2 ** 20):
        self.memory_limit = memory_limit
        self.column_weights = column_weights
        self.columns = columns
        self.rows = rows
        self.chain_complex = None
        self.cycle = cycle
        if full_setup:
            self.compute_treewidth()

    def set_cycle(self, cycle):
        self.cycle = cycle

    def compute_treewidth(self, connectivity_graph=True, spine_graph=True):
        self.chain_complex = ChainComplex(self.columns, self.rows, self.column_weights, connectivity_graph = connectivity_graph, spine_graph=spine_graph)

    def treewidth_of_spine_graph(self):
        return self.chain_complex.sgTD.tree_width

    def treewidth_of_connectivity_graph(self):
        return self.chain_complex.cgTD.tree_width

    def number_of_bags_in_decomposition_of_spine_graph(self):
        return self.chain_complex.sg_NTD_size()

    def number_of_bags_in_decomposition_of_connectivity_graph(self):
        return self.chain_complex.cg_NTD_size()

    def homology_localization(self, print_status=False):
        return per_hom_loc_sg_rep(self.chain_complex, self.cycle, print_status, self.memory_limit)

    def nice_tree_decomposition(self, spine_graph=True):
        if spine_graph:
            return self.chain_complex.sgTD.nice_TD
        else:
            return self.chain_complex.cgTD.nice_TD

