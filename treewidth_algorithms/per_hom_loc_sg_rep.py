
import math

from treewidth_algorithms.error_file import MemoryLimitViolation
from treewidth_algorithms.tools.bit_magic import MyUniverse, Translator


def per_hom_loc_sg_rep(cc, v, give_status, memory_limit):
    """Finds the weight of the shortest cycle in the homology class of the input cycle in the given chain complex."""

    memory = 0
    memory_peak = 0


    # Initialises a dictionary of tables that will be filled in recursively
    tables = {}
    u = {}

    for bag_id in range(0, cc.sg_NTD_size()):
        u0 = MyUniverse()
        u1 = MyUniverse()
        u[bag_id] = {0: u0, 1: u1}

    for bag_id in range(0, cc.sg_NTD_size()):  # This ordering of bag ID's works because this ordering is topological.
        bag_type = cc.sg_bag_type(bag_id)
        children_id = cc.sg_children(bag_id)
        parsed = False

        # The base case where there is only one table entry to fill, the one containing one entry where Q,P=emptyset:
        if bag_type == "leaf" and len(children_id) == 0:
            # Looking up the value at bag id in the table we get, (Q = Ø -> (P = Ø -> solution value = 0,representative solution).
            tables[bag_id] = {0: {0: (0.0, [])}}
            memory = memory + 1
            parsed = True

        # Introducing a new simplex to the bag
        if bag_type == "introduce" and len(children_id) == 1:
            child_table = tables[children_id[0]]
            tables[bag_id], memory = introduce_sg_representant(child_table, bag_id, children_id[0], cc, v, u, memory,
                                                               memory_limit)
            parsed = True

        # Forgetting an old simplex in the bag
        if bag_type == "forget" and len(children_id) == 1:
            child_table = tables[children_id[0]]
            tables[bag_id], memory = forget_sg_representant(child_table, bag_id, children_id[0], cc, v, u, memory,
                                                            memory_limit)
            parsed = True

        # Joining two equal bags
        if bag_type == "join" and len(children_id) == 2:
            child_bag_id = cc.sg_children(bag_id)
            child_table1 = tables[children_id[0]]
            child_table2 = tables[children_id[1]]
            tables[bag_id], memory = join_sg_representant(child_table1, child_table2, bag_id, child_bag_id, cc, v, u,
                                                          memory, memory_limit)
            parsed = True

        # Make this an error message:
        if not parsed:
            print("Mismatch between bagtype and number of children at bag", bag_id, bag_type)

        if give_status:
            optimal = math.inf
            for q in tables[bag_id]:
                for p in tables[bag_id][q]:
                    (value, rep) = tables[bag_id][q][p]
                    optimal = min(optimal, value)

            print("The ", bag_id, "'th bag of type ", cc.sg_bag_type(bag_id), " and of size ",
                  len(cc.sg_bag_content(bag_id)), " has an optimal solution of size ", optimal, " and children ",
                  children_id)

        memory_peak = max(memory, memory_peak)


        for child_id in children_id:
            child_tab = tables.pop(child_id)
            for q in child_tab:
                memory = memory - len(child_tab[q])

    check1 = False

    check2 = False

    check3 = False

    # We get the table for the root bag which should contain precisely one entry for the empty set
    root_table = tables[cc.sg_NTD_size() - 1]
    if len(root_table) == 1:
        check1 = True

    # This table should contain a list containing a triple p = Ø, a rep. solution and the value of that solution
    entries = root_table[0]

    if len(entries) == 1:
        # print("The number of entries at Q = Ø is correct")
        check2 = True
    p = 0

    val, rep = entries[p]

    if p == 0:
        # print("p is empty at root")
        check3 = True

    if check1 and check2 and check3:
        allgood = True
        # print("All seemingly good:")

    return val, rep, memory_peak


def introduce_sg_representant(child_table, bag_id, child_id, cc, v, u, memory, memory_limit):
    """Fills out the table for the introduce-bag specified at input"""
    # The introduced simplex
    sigma_set = cc.sg_difference(bag_id, child_id)

    # A check on the length of sigma_set
    if len(sigma_set) != 1:
        print("The introduce bag introduces either more or fewer than one simplex")

    # "Unpack" the simplex from the set
    sigma, *_ = sigma_set

    # The dimension of the introduced simplex
    sigma_up = cc.sg_simplex_up(sigma)

    if sigma_up:
        return introduce_up_representant(child_table, sigma, sigma_set, bag_id, child_id, cc, v, u, memory,
                                         memory_limit)
    else:
        return introduce_down_representant(child_table, sigma, sigma_set, bag_id, child_id, cc, v, u, memory,
                                           memory_limit)


def introduce_up_representant(child_table, sigma_node, sigma_set_node, bag_id, child_id, cc, v, u, memory,
                              memory_limit):
    """Fills out the table for the introduce-bag specified at input when the introduced simplex is of dimension d+1."""
    sigma, x = sigma_node
    sigma_set = frozenset([a for (a, b) in sigma_set_node])

    u[bag_id][0] = u[child_id][0]
    u[bag_id][1] = u[child_id][1]

    # The boundary of the introduced simplex (as a frozenset)
    bnd = cc.get_bnd(sigma)

    # The simplices of the child bag
    simps = cc.sg_bag_content(child_id)

    # The d-simplices of the child bag
    up_simps_node = frozenset([simp for simp in simps if cc.sg_simplex_up(simp)])
    up_simps = frozenset([a for (a, b) in up_simps_node])

    # The d+1-simplices of the child bag
    down_simps_node = frozenset([simp for simp in simps if not cc.sg_simplex_up(simp)])
    down_simps = frozenset([a for (a, b) in down_simps_node])

    # The boundary of the introduced simplex in the bag (as a frozenset)
    bnd_sig = frozenset.intersection(down_simps, bnd)

    # The skeleton of the child bag

    # Initializing the table that will be filled out
    table = {}

    u[bag_id][0].add_element(sigma, 0.0)
    sigma_set_int = u[bag_id][0].set_to_int(sigma_set)

    bnd_sig_int = u[bag_id][1].set_to_int(bnd_sig)

    # For every subset of nodes in the child bag
    for q in child_table:

        # The set q union sigma
        q_sig_int = u[bag_id][0].union(q, sigma_set_int)

        # Make a dictionary for q and q_sigma

        table[q] = {}
        table[q_sig_int] = {}

        # For every entry in the dictionary of q in the child bag
        for p in child_table[q]:
            value, rep = child_table[q][p]
            rep = rep.copy()
            p_sig_int = u[bag_id][1].sym_dif(p, bnd_sig_int)
            value_sig = value
            table[q_sig_int][p_sig_int] = (value_sig, rep)
            table[q][p] = (value, rep)
            memory = memory + 2
            if memory > memory_limit:
                raise MemoryLimitViolation

    return table, memory


# Keep in mind that we are introducing a rho not a sigma here.
def introduce_down_representant(child_table, rho_node, rho_set_node, bag_id, child_id, cc, v, u, memory, memory_limit):
    """Fills out the table for the introduce-bag specified at input when the introduced simplex is of dimension d."""
    rho, x = rho_node
    rho_set = frozenset([a for (a, b) in rho_set_node])

    u[bag_id][0] = u[child_id][0]
    u[bag_id][1] = u[child_id][1]

    # The simplices of the child bag
    simps = cc.sg_bag_content(child_id)

    # The d-simplices of the child bag
    up_simps = frozenset([simp for simp in simps if cc.sg_simplex_up(simp)])

    # The parents of the introduced simplex in the bag (as a frozenset)
    cofaces_node = cc.sg_cofaces(rho_node)

    cofaces_rho_node = frozenset.intersection(up_simps, cofaces_node)
    cofaces_rho = frozenset([a for (a, b) in cofaces_rho_node])

    # Initializing the table that will be filled out
    table = {}

    u[bag_id][1].add_element(rho, 0)
    rho_set_int = u[bag_id][1].set_to_int(rho_set)

    cofaces_rho_int = u[bag_id][0].set_to_int(cofaces_rho)

    v_contribution_nr = 0
    if rho in v:
        v_contribution_nr = 1

    # For every subset of nodes in the child bag
    for q in child_table:

        q_contribution_int = u[bag_id][0].intersection(q, cofaces_rho_int)
        q_contribution_nr = u[bag_id][0].size(q_contribution_int)

        # Make a dictionary for q and q_sigma
        rho_in = (1 == ((q_contribution_nr + v_contribution_nr) % 2))

        table[q] = {}
        # For every entry in the dictionary of q in the child bag
        for p in child_table[q]:
            value, rep = child_table[q][p]
            rep = rep.copy()
            if rho_in:
                table[q][u[bag_id][1].union(p, rho_set_int)] = (value, rep)
            else:
                table[q][p] = child_table[q][p]
            memory = memory + 1
            if memory > memory_limit:
                raise MemoryLimitViolation

    return table, memory


def forget_sg_representant(child_table, bag_id, child_id, cc, v, u, memory, memory_limit):
    """Fills out the table for the forget-bag specified at input."""
    u[bag_id][0] = u[child_id][0]
    u[bag_id][1] = u[child_id][1]

    # The forgotten simplex
    sigma_set = cc.sg_difference(child_id, bag_id)

    # A check on the length of sigma_set
    if len(sigma_set) != 1:
        print("The introduce bag introduces either more or fewer than one simplex")

    # "Unpack" the simplex from the set
    sigma, *_ = sigma_set

    # The dimension of the introduced simplex
    sigma_up = cc.sg_simplex_up(sigma)

    if sigma_up:
        return forget_up_representant(child_table, sigma, sigma_set, bag_id, child_id, cc, v, u, memory, memory_limit)
    else:
        return forget_down_representant(child_table, sigma, sigma_set, bag_id, child_id, cc, v, u, memory, memory_limit)


def forget_up_representant(child_table, sigma_node, sigma_set_node, bag_id, child_id, cc, v, u, memory, memory_limit):
    """Fills out the table for the forget-bag specified at input when the introduced simplex is of dimension d+1."""
    sigma, x = sigma_node
    sigma_set = frozenset([a for (a, b) in sigma_set_node])

    # Initializing the table that will be filled out
    table = {}

    sigma_set_int = u[bag_id][0].set_to_int(sigma_set)

    # For every entry in the dictionary of q in the child bag
    for q in child_table:
        q_sig_int = u[bag_id][0].dif(q, sigma_set_int)

        # Checks if there is a list at table[q_sig] already
        if not q_sig_int in table:
            # If this is not the case, initializes this list
            table[q_sig_int] = {}

        # For every triple in the list at the child, q, we compute the new p_sig
        for p in child_table[q]:
            value, rep = child_table[q][p]


            # We check if the new p_sig is already in the list of tuples stored at q_sig. If it is we compare it to the old one and update it if neccessary
            if q_sig_int != q:
                new_value = value + cc.get_weight(sigma)
                new_rep = rep + list(u[bag_id][0].int_to_set(sigma_set_int))
                if p in table[q_sig_int]:
                    value_other, repx = table[q_sig_int][p]
                    if new_value < value_other:
                        table[q_sig_int][p] = (new_value, new_rep)

            # If no match for p_sig was found we add it to the list.
                else:
                    table[q_sig_int][p] = (new_value, new_rep)
                    memory = memory + 1
                    if memory > memory_limit:
                        raise MemoryLimitViolation
            else:
                if p in table[q_sig_int]:
                    value_other, repx = table[q_sig_int][p]
                    if value < value_other:
                        table[q_sig_int][p] = (value, rep)
                # If no match for p_sig was found we add it to the list.
                else:
                    table[q_sig_int][p] = (value, rep)
                    memory = memory + 1
                    if memory > memory_limit:
                        raise MemoryLimitViolation

    u[bag_id][0].remove_element(sigma)
    return table, memory


def forget_down_representant(child_table, rho_node, rho_set_node, bag_id, child_id, cc, v, u, memory, memory_limit):
    """Fills out the table for the forget-bag specified at input when the introduced simplex is of dimension d."""
    rho, x = rho_node
    rho_set = frozenset([a for (a, b) in rho_set_node])

    # Initializing the table that will be filled out
    table = {}

    rho_set_int = u[bag_id][1].set_to_int(rho_set)

    # For every entry in the dictionary of q in the child bag
    for q in child_table:

        table[q] = {}

        # For every triple in the list at the child, q, we compute the new p_sig
        for p in child_table[q]:

            p_rho_old_int = u[bag_id][1].dif(p, rho_set_int)

            value, rep = child_table[q][p]
            rep = rep.copy()

            if u[bag_id][1].intersection(p, rho_set_int) == 0:
                table[q][p_rho_old_int] = (value, rep)
                memory = memory + 1
                if memory > memory_limit:
                    raise MemoryLimitViolation

    u[bag_id][1].remove_element(rho)
    return table, memory


def join_sg_representant(child_table1, child_table2, bag_id, child_bag_id, cc, v, u, memory, memory_limit):
    """Fills the table for the join-bag specified at input."""
    u[bag_id][0] = u[child_bag_id[0]][0]
    u[bag_id][1] = u[child_bag_id[0]][1]

    u2q = u[child_bag_id[1]][0]
    u2p = u[child_bag_id[1]][1]

    # Translates opposite ways!
    tq = Translator(u[bag_id][0], u2q)
    tp = Translator(u2p, u[bag_id][1])

    # The simplices of the child bag
    simps_node = cc.sg_bag_content(bag_id)

    # The d-simplices of the child bag
    up_simps_node = frozenset([simp for simp in simps_node if cc.sg_simplex_up(simp)])

    # The d+1-simplices of the child bag
    down_simps_node = frozenset([simp for simp in simps_node if not cc.sg_simplex_up(simp)])
    down_simps = frozenset([a for (a, b) in down_simps_node])


    v_bag = frozenset.intersection(v, down_simps)

    # Initializes the table that will be filled out during this call
    table = {}

    # Get integer representation of the simplices of v in the bag
    v_bag_int = u[bag_id][1].set_to_int(v_bag)

    down_simps_int = u[bag_id][1].set_to_int(down_simps)

    for q in child_table1:

        # Initializes the list at q
        table[q] = {}

        q_elems = u[bag_id][0].int_to_set(q)

        # Gets the boundary of q
        q_elems_node = frozenset(list(zip(q_elems, [True for a in q_elems])))
        bound_q = cc.boundary_map(q_elems)
        bound_q_in_bag = frozenset.intersection(bound_q, down_simps)

        bound_q_in_bag_int = u[bag_id][1].set_to_int(bound_q_in_bag)

        q2 = tq.int_to_int(q)

        # For every entry at Q in left bag
        for p1 in child_table1[q]:
            val1, rep1 = child_table1[q][p1]
            rep1 = rep1.copy()

            # For every entry at Q in the right bag
            for p2 in child_table2[q2]:

                p2_transformed = tp.int_to_int(p2)
                val2, rep2 = child_table2[q2][p2]
                rep2 = rep2.copy()

                # Compute the new representative
                p_new_int = u[bag_id][1].sym_dif(p1, u[bag_id][1].sym_dif(p2_transformed,
                                                                          u[bag_id][1].sym_dif(bound_q_in_bag_int,
                                                                                               v_bag_int)))

                # Compute the value of that solution
                new_value = val1 + val2

                # Keep track of if this entry is already filled
                if p_new_int in table[q]:
                    value_other, repx = table[q][p_new_int]
                    if new_value < value_other:
                        table[q][p_new_int] = (new_value, rep1 + rep2)

                # If no match for p_sig was found we add it to the list.
                else:
                    table[q][p_new_int] = (new_value, rep1 + rep2)
                    memory = memory + 1
                    if memory > memory_limit:
                        raise MemoryLimitViolation

    return table, memory
