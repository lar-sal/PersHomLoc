class MyUniverse(object):
    """A translator for sending (small) frozensets to their bit representation.

    Instances of this class organizes all relevant information about the translation procedure as well as keeping track
    of which elements are in the universe. It is also a wrapper for typical set operations.

    Attributes:
        element_to_position: A dictionary mapping an element to the position of its digit.
        position_to_element: A dictionary mapping the position of a digit to its element.
        element_cost: A dictionary mapping an element to its cost.
        position_cost: A dictionary mapping the position of a digit to its cost
        free_positions: A list keeping track of the unoccupied digit positions.
        """

    def __init__(self):
        self.element_to_position = {}
        self.position_to_element = {}
        self.element_cost = {}
        self.position_cost = {}
        self.free_positions = [1]

    def position(self, elem):
        """Finds the digit position associated to an element."""
        return self.element_to_position[elem]

    def element(self, pos):
        """Finds the element associated to a digit position."""
        return self.position_to_element[pos]

    def add_element(self, elem, cost):
        """Adds an element to the universe with a given cost"""
        pos = self.free_positions.pop(0)
        if len(self.free_positions) == 0:
            self.free_positions.append(pos * 2)
        self.element_to_position[elem] = pos
        self.position_to_element[pos] = elem
        self.element_cost[elem] = cost
        self.position_cost[pos] = cost

    def add_elements(self, elems_costs):
        """Adds several elements from a set to the universe"""
        for (elem, cost) in elems_costs:
            self.add_element(elem, cost)

    def remove_element(self, elem):
        """Removes an element from the universe"""
        pos = self.position(elem)
        del self.element_to_position[elem]
        del self.position_to_element[pos]
        del self.element_cost[elem]
        del self.position_cost[pos]
        self.free_positions.append(pos)
        self.free_positions.sort()

    def remove_elements(self, elems):
        """Removes several element from a set from the universe"""
        for elem in elems:
            self.remove_element(elem)

    def set_to_int(self, elements):
        """Takes a frozenset and returns an integer representation of the set"""
        outint = 0
        for elem in elements:
            outint = outint | self.position(elem)
        return outint

    def int_to_set(self, in_int):
        """Takes an integer representation of a set and returns a frozenset containing its element."""
        out_set = []
        for pos in self.position_to_element:
            if (pos & in_int) != 0:
                out_set.append(self.position_to_element[pos])
        return frozenset(out_set)

    def sym_dif(self, int1, int2):
        """Takes the symmetric difference of two sets represented as integers."""
        return int1 ^ int2

    def dif(self, int1, int2):
        """Takes the difference of one set with another where both sets are represented as integers."""
        return int1 ^ (int1 & int2)

    def union(self, int1, int2):
        """Takes the union of two sets both represented as integers."""
        return int1 | int2

    def intersection(self, int1, int2):
        """Takes the intersection of two sets both represented as integers."""
        return int1 & int2

    def size(self, in_int):
        """Computes the size of a set represented as an integer."""
        size = 0
        for pos in self.position_to_element:
            if (pos & in_int) != 0:
                size = size + 1
        return size

    def cost(self, in_int):
        """Computes the cost of a set represented as integers."""
        price = 0
        for pos in self.position_to_element:
            if pos & in_int != 0:
                price = price + self.position_cost[pos]
        return price


class Translator(object):
    """A translator for transfering information of one universe giving bit representations of a set to another universe.

    Attributes:
        u1: The input universe.
        u2: The output universe.
        in_pos_to_out_pos: A dictionary sending the digit position of an element in the input universe to the position
        of the corresponding digit position of the element in the output universe.
        """
    def __init__(self, InputUniverse, OutputUniverse):
        self.u1 = InputUniverse
        self.u2 = OutputUniverse
        self.in_pos_to_out_pos = self.p_to_p()

    def p_to_p(self):
        """Method for generating the position to position translator dictionary"""
        p_to_p_dict = {}
        for pos_in in self.u1.position_to_element:
            elem = self.u1.position_to_element[pos_in]
            pos_out = self.u2.element_to_position[elem]
            p_to_p_dict[pos_in] = pos_out
        return p_to_p_dict

    def int_to_int(self, in_int):
        """Method for sending a set represented as an integer in the first universe to a set represented as an integer
        in the second universe."""
        out_int = 0
        for pos in self.u1.position_to_element:
            if (pos & in_int) != 0:
                out_int = out_int | self.in_pos_to_out_pos[pos]
        return out_int