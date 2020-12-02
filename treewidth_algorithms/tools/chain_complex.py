# File containing the ChainComplex class and the TD
import networkx as nx
from networkx.algorithms import approximation


# TODO: change name to SimplicialComplex
class ChainComplex(object):
    """An interface for storing, accessing and generating information about a simplicial complex.

    Instances of this class organizes all relevant information about the input (boundary?)-matrix needed to run the
    homology localization treewidth_algorithms in this package. It builds upon the homloc code and serves as the interface between
    these treewidth_algorithms and the input (now represented as an interface object). As an example it can be used to get the
    image of a vector, to get information about the tree decomposition of the spine graph and the connectivity graph of
    the boundary matrix and to get the weight of simplices.

    Attributes:
        column_size: The number of columns in the matrix
        row_size: The number of rows in the matrix
        column_vector: A dictionary for converting the ID's of columns to its column vecotr
        column_weights: A dictionary for finding the weight of a column vector by its ID
        cg: The connectivity graph of the simplicial complex
        cgTD: A TreeDecomposition of the connectivity graph of the simplicial complex
        sg: The spine graph of the simplicial complex
        sgTD: A TreeDecomposition of the spine graph of the simplicial complex
    """

    def __init__(self, columns, rows, column_weights, connectivity_graph=True, spine_graph=True):

        # The number of d+1 simplices
        self.column_size = len(columns)

        # The number of d simplices
        self.row_size = len(rows)

        # The dict takes the ID of a column to its column vector represented as a frozen set
        self.column_vector = {index: frozenset(columns[index]) for index in range(0, len(columns))}

        # The dict takes a column to its column weight
        self.column_weights = column_weights

        if connectivity_graph:
            # The conectivity graph of the simplicial complex
            self.cg = self.con_graph()

            # The (nice and normal) treedecomposition of the connectivity graph
            self.cgTD = TD(self.cg)
            self.cgTD.get_niceTD()

        if spine_graph:
            # The spine graph of the simplicial complex
            self.sg = self.spine_graph()

            # The treedecomposition of the d-spine graph
            self.sgTD = TD(self.sg)
            self.sgTD.get_niceTD()

    # Gets the number of d + 1 simplices
    def number_of_d_pluss_1_simplices(self):
        """Returns the number of d+1 simplices."""
        return self.column_size

    # Gets the weight of a column
    def get_weight(self, column):
        """Returns the weight of a column."""
        return self.column_weights[column]

    def get_weight_set(self, columns):
        """Returns the weight of a set/frozenset/list of columns."""
        return sum([self.get_weight(column) for column in columns])

    def get_bnd(self, column_id):
        """Returns the vector corresponding to a column ID."""
        return self.column_vector[column_id]

    def boundary_map(self, columns):
        """Returns the boundary (represented as a frozenset of row ID's) of a set/frozenset/list of columns as ID's."""
        out = frozenset([])
        for column in columns:
            column = frozenset(self.get_bnd(column))
            out = frozenset.symmetric_difference(out, column)
        return out

    def d_skeleton(self, columns):
        """Finds the rows adjacent to a set of column ID's."""
        skel = frozenset([])
        for column in columns:
            skel = frozenset.union(skel, self.get_bnd(column))
        return skel

    def is_cg_edge(self, simp1, simp2):
        """Checks if a pair of d+1 columns are adjacent in the intersection graph."""
        intersection = frozenset.intersection(self.get_bnd(simp1), self.get_bnd(simp2))
        return intersection != frozenset([]) and simp1 != simp2

    def con_graph(self):
        """Makes the connectivity graph from the simplicies."""
        nodes = list(range(0, self.column_size))
        edges = [(a, b) for a in nodes for b in nodes if self.is_cg_edge(a, b)]
        graph = nx.Graph()
        graph.add_nodes_from(nodes)
        graph.add_edges_from(edges)
        return graph

    def spine_graph(self):
        """Makes the spine graph from the simplicies."""

        # Makes the nodes for the d+1 simplices
        nodes_up = list(zip(list(range(0, self.column_size)), [True for i in range(0, self.column_size)]))

        # Makes the nodes for the d simplices
        nodes_down = list(zip(list(range(0, self.row_size)), [False for i in range(0, self.row_size)]))

        # Initializes the graph and adds nodes
        graph = nx.Graph()
        graph.add_nodes_from(nodes_up)
        graph.add_nodes_from(nodes_down)

        # Adds the edges for every d+1 simplex
        for node in nodes_up:
            bnd = self.get_bnd(node[0])
            for down in bnd:
                d = tuple([down, False])
                graph.add_edge(node, d)
        return graph

    def get_NTD_cg(self):
        """Returns the chosen nice tree decomposition of the connectivity graph."""
        return self.cgTD.nice_TD

    def get_NTD_sg(self):
        """Returns the chosen nice tree decomposition of the spine graph."""
        return self.sgTD.nice_TD

    def cg_NTD_size(self):
        """Returns the number of bags in the nice tree decomposition of the connectivity graph."""
        return len(self.cgTD.nice_TD.nodes)

    def sg_NTD_size(self):
        """Returns the number of bags in the nice tree decomposition of the spine graph."""
        return len(self.sgTD.nice_TD.nodes)

    def cg_bag_type(self, bag_id):
        """Returns the bag type of a bag id in the nice tree decomposition of the connectivity graph."""
        return self.cgTD.bag_type(bag_id)

    def sg_bag_type(self, bag_id):
        """Returns the bag type of a bag id in the nice tree decomposition of the spine graph."""
        return self.sgTD.bag_type(bag_id)

    def cg_bag_content(self, bag_id):
        """Returns the ID's of vertices in the bag in the nice tree decomposition of the connectivity graph."""
        return self.cgTD.bag(bag_id)

    def sg_bag_content(self, bag_id):
        """Returns the ID's of vertices in the bag in the nice tree decomposition of the spine graph."""
        return self.sgTD.bag(bag_id)

    def cg_children(self, bag_id):
        """Returns the ID's of the children of a bag in the nice tree decomposition of the connectivity graph."""
        return self.cgTD.children(bag_id)

    def sg_children(self, bag_id):
        """Returns the ID's of the children of a bag in the nice tree decomposition of the spine graph."""
        return self.sgTD.children(bag_id)

    def cg_difference(self, bag_id1, bag_id2):
        """Returns the elements of the first bag that are not in the second bag of the connectivity graph."""
        return self.cgTD.get_difference(bag_id1, bag_id2)

    def sg_difference(self, bag_id1, bag_id2):
        """Returns the elements of the first bag that are not in the second bag of the spine graph."""
        return self.sgTD.get_difference(bag_id1, bag_id2)

    def sg_simplex_up(self, simplex):
        """Checks if a node in the spine graph corresponds to a d+1-dimensional simplex."""
        x, y = simplex
        return y

    def sg_cofaces(self, simplex):
        """Returns the d+1-cofaces of the input d-simplex in a frozen set."""
        return frozenset(self.sg.neighbors(simplex))


# TODO: Change name to TreeDecomposition
class TD(object):
    """An interface for storing, accessing and generating information about the tree decompositions of a graph.

     Objects of this class computes and stores (nice) tree decompositions and provides an interface for getting
     information about the bags of these tree decompositions.

    Attributes:
        G: The underlying graph that we get a tree decomposition of.
        tree_width: The tree width of the tree decomposition we are currently working with.
        tree_decomposition: The tree decomposition we are currently working with which is generally not nice.
        nice_TD: A nice tree decomposition made from the tree decomposition above.
        id_to_bag: A dictionary for associating a bag to each bag ID.
        id_to_type: A dictionary for associating a bag type to each bag ID.
        topord: A list giving a topological ordering of the bags in the nice tree decomposition rooted at an empty bag.
        toporddict: A dictionary giving the position of a bag in the topological ordering.
        """

    def __init__(self, graph):
        self.G = graph

        # Get treedecomposition using the first (min_degree) heuristic
        tw1 = approximation.treewidth.treewidth_min_degree(self.G)

        # Get treedecomposition using the second (fill_in) heuristic
        tw2 = approximation.treewidth.treewidth_min_fill_in(self.G)

        # We let the treedecomposition and treewidth be the smallest of the two heuristics
        if tw1[0] < tw2[0]:
            # The tree width of the input
            self.tree_width = tw1[0]

            # The tree decomposition of the input
            self.tree_decomposition = tw1[1]
        else:
            # The tree width of the input
            self.tree_width = tw2[0]

            # The tree decomposition of the input
            self.tree_decomposition = tw2[1]

        self.nice_TD = None

        self.id_to_bag = None

        self.id_to_type = None

        self.topord = None

        self.toporddict = None

    def get_niceTD(self):
        """Computes and stores a nice tree decomposition"""
        self.nice_TD = nx.Graph()

        # Topologically sorts the nodes of the "bad" treedecompositions
        topord = list(nx.dfs_postorder_nodes(self.tree_decomposition))

        # Makes a dictionary, saying which bag has which position in the "bad" treedecompositions for easy comparisons later on
        toporddict = dict(zip(topord, list(range(0, len(topord)))))

        # Initializes a new index dictionary to look up which index the child of a bag will have in the end.
        index_dict = {}
        # Initializes bag id to bag dictionary
        self.id_to_bag = {}
        # Initializes bag id to bag type (leaf, introduce, forget etc.)
        self.id_to_type = {}
        # Initializes the index of the bag we are currently working on

        index = 0
        for node in topord:
            children = [neigh for neigh in self.tree_decomposition.neighbors(node) if
                        toporddict[neigh] < toporddict[node]]
            childless = False
            if len(children) == 0:
                childless = True
                children.append(frozenset())
                self.nice_TD.add_node(index)
                self.id_to_bag[index] = frozenset()
                self.id_to_type[index] = "leaf"
                index = index + 1

            join_id = []
            for child in children:

                if childless:
                    child_index = index - 1
                else:
                    child_index = index_dict[child]

                # The elements that are in child and not in node
                remove = frozenset.difference(child, node)

                # The elements that are in node but not in child
                add = frozenset.difference(node, child)

                for element in remove:
                    child = frozenset.difference(child, frozenset([element]))
                    self.id_to_bag[index] = child
                    self.id_to_type[index] = "forget"
                    self.nice_TD.add_node(index)
                    self.nice_TD.add_edge(child_index, index)
                    child_index = index
                    index = index + 1

                for element in add:
                    child = frozenset.union(child, frozenset([element]))
                    self.id_to_bag[index] = child
                    self.id_to_type[index] = "introduce"
                    self.nice_TD.add_node(index)
                    self.nice_TD.add_edge(child_index, index)
                    child_index = index
                    index = index + 1

                if child != node:
                    print("Something is wrong.")

                join_id.append(index - 1)

            while len(join_id) != 1:
                a = join_id.pop()
                b = join_id.pop()
                join_id.append(index)
                self.id_to_bag[index] = node
                self.id_to_type[index] = "join"
                self.nice_TD.add_node(index)
                self.nice_TD.add_edge(a, index)
                self.nice_TD.add_edge(b, index)
                index = index + 1

            index_dict[node] = join_id[0]

        child = topord[-1]
        index = index_dict[child] + 1

        for element in child:
            current = frozenset.difference(child, frozenset([element]))
            self.id_to_bag[index] = current
            self.id_to_type[index] = "forget"
            self.nice_TD.add_node(index)
            self.nice_TD.add_edge(index - 1, index)
            child = current
            index = index + 1

        self.topord = list(nx.dfs_postorder_nodes(self.nice_TD, len(self.nice_TD.nodes()) - 1))
        self.toporddict = dict(zip(self.topord, list(range(0, len(self.topord)))))

    def get_children(self, node):
        """Returns a list of the children of a bag in the topological ordering of the (not nice) tree decomposition."""
        return [neigh for neigh in self.tree_decomposition.neighbors(node) if
                self.toporddict[neigh] < self.toporddict[node]]

    def children(self, node):
        """Returns a list of the children of a bag in the topological ordering of the nice tree decomposition."""
        return [neigh for neigh in self.nice_TD.neighbors(node) if neigh < node]

    def bag(self, bag_id):
        """Finds the  elements of a bag."""
        return self.id_to_bag[bag_id]

    def bag_type(self, bag_id):
        """Finds the type of a bag."""
        return self.id_to_type[bag_id]

    def get_difference(self, bag_id1, bag_id2):
        """Finds the difference between the elements of two bags."""
        return frozenset.difference(self.id_to_bag[bag_id1], self.id_to_bag[bag_id2])
