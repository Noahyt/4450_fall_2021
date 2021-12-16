"""Code for computing FDM."""

import numpy as np
import copy


def compute_fdm(node_coord, edge_node_ids, nodal_loads, force_densities, fixed_node_ids):
    """Using the Force Density Method (FDM) to compute a surface structure
    that exhibit membrane behavior only, with no bending moments

    References: Linkwitz (2014) article, pp. 64-65, on the class website.

    Parameters
    ----------
    node_coord : a list of 3-float list
        a list of nodes' coordinates
        [[node_0_x, node_0_y, node_0_z],
         [node_1_x, node_1_y, node_1_z],
         ...
         [node_{nn-1}_x, node_{nn-1}_y, node_{nn-1}_z],
         ]
        Here nn means the number of nodes.
    edge_node_ids : a list of 2-int list
        a list indicating the edge connectivity in the base grid
        [[e 0 start node id, e 0 end node id],
         [e 1 start node id, e 1 end node id],
         ...,
         [e {ne-1} start node id, e {ne-1} end node id],
        ]
        Here ne means the number of egdes.
    nodal_loads : a list of 3-float list
        a list of nodal point loads
        If given [[load_x, load_y, load_z]], a list of length 1,
        we will apply uniform nodal load to all the nodes.

        Or, if the given loads has length > 1, it must has length nn (number of nodes):
        [[load_0_x, load_0_y, load_0_z],
         [load_1_x, load_1_y, load_1_z],
         ...
         [load_{nn-1}_x, load_{nn-1}_y, load_{nn-1}_z],
         ]
    force_densities : a list of float
        a list of force densities for each edge
        If given [q], a list of length 1,
        we will apply uniform force densities to all the edges.

        Or, if the given force densities has length > 1, it must has length ne (number of edges):
        [q_0,
         q_1,
         ...
         q_{ne-1}]
    fixed_node_ids : a list of int
        Fixed nodes' indices
        [fixed node 1 id, fixed node 2 id, ...]

    Returns
    -------
    fdm_node_coord : a list of 3-float lists
        a list of optimized nodes' coordinates
        [[opt_node_0_x, opt_node_0_y, opt_node_0_z],
         [opt_node_1_x, opt_node_1_y, opt_node_1_z],
         ...
         [opt_node_{nn-1}_x, opt_node_{nn-1}_y, opt_node_{nn-1}_z],
         ]

    """
    # convert list into a numpy array
    node_coord = np.array(node_coord)
    original_coord = copy.copy(node_coord)

    # number of nodes
    nn = len(node_coord)
    assert node_coord.shape[1] == 3, 'node_coord must have 3 columns (x, y, z)!'

    # number of edges
    ne = len(edge_node_ids)

    # Free & Fix nodal indices
    free = [i for i in range(nn) if i not in fixed_node_ids]

    # Parse input point loads
    p = np.zeros([nn, 3])
    if len(nodal_loads) == 1:
        # if input loads only have one vector,
        # assume uniform load on every node
        for i in range(nn):
            p[i, 0] = nodal_loads[0][0]
            p[i, 1] = nodal_loads[0][1]
            p[i, 2] = nodal_loads[0][2]
    else:
        if len(nodal_loads) != nn:
            raise ValueError(
                "Error: loads vector length not equal to num of nodes")
        else:
            for i in range(nn):
                p[i, 0] = nodal_loads[i][0]
                p[i, 1] = nodal_loads[i][1]
                p[i, 2] = nodal_loads[i][2]

    # Parse force Densities
    q = force_densities
    if len(q) == 1:
        Q = np.identity(ne) * q[0]
    else:
        if len(q) != ne:
            raise ValueError(
                "Error: force densities vector length not equal num of edges")
        else:
            Q = np.identity(ne)
            for i, qq in enumerate(q):
                Q[i, i] = qq

    ##################################
    # YOUR CODE STARTS FROM HERE
    # Note: I put some comments to outline some steps you have to fill in,
    # but you don't have to strictly follow the order.
    ##################################
    # rename some variables.
    node_count = nn
    nodes = np.arange(node_count)
    edge_count = ne
    edges = np.array(edge_node_ids, dtype=int)
    edge_index = np.arange(edge_count)

    C = np.zeros([edge_count, node_count])

    # Set Start.
    C[edge_index, edges[:, 0]] = -1

    # Set End.
    C[edge_index, edges[:, 1]] = 1

    # Define permutation to take free nodes to first n positions.
    # Sort so that the index does not change for later elements.
    fixed = sorted(fixed_node_ids)
    free_ = list(set(nodes) - set(fixed))
    free_count = len(free_)
    permutation = free_ + fixed

    C = C[:, permutation]

    # Reorder C's rows into two matrices C_i (elements with free nodes),
    # and C_f (elements with fixed nodes)
    C_n = C[:, :free_count]
    C_f = C[:, free_count:]

    # reorder node_coordinate's rows into free and fixed in the same fashion
    node_coord = node_coord[permutation]
    x_f = node_coord[free_count:]

    # reorder load p's rows into free and fixed in the same fashion
    p = p[permutation]
    p_n = p[:free_count]

    # Using convention from Linkwitz 6.31
    D_n = np.dot(np.transpose(C_n), np.dot(Q, C_n))
    D_f = np.dot(np.transpose(C_n), np.dot(Q, C_f))

    a = D_n
    b = p_n - np.dot(D_f, x_f)

    coord_ = np.linalg.solve(a, b)

    # node_coord[free] = coord_
    # Replace in original index in coordinate array.
    for i, o_id in enumerate(free):
        original_coord[o_id] = coord_[i]

    return original_coord
