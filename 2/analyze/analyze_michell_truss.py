"""Analyzes Michell Truss."""
import numpy as np
from anastruct import SystemElements


class MichellTruss(object):

    def __init__(self, positions):

        self.ss = SystemElements()

        ur, lr = assign_node_reference(positions)
        bonds = build_connections(ur, lr)

        upper_positions = positions
        lower_positions = _generate_lower_positions(positions)

        pos, ref = flatten_and_remove_duplicate(
            upper_positions, ur, lower_positions, lr)

        pos = [p for _, p in sorted(zip(ref, pos))]
        ref = [r for r, _ in sorted(zip(ref, pos))]

        pos = np.array(pos)
        bonds = np.array(bonds)

        # Generate list of bonds w coordiantes.
        self.bond_in_coordinate = pos[bonds]

        # Add structure.
        for b in self.bond_in_coordinate:
            print("adding {}".format(b))
            self.ss.add_truss_element(b)

        self.tip_node = self.ss.find_node_id(upper_positions[-1][-1])
        self.upper_support = self.ss.find_node_id(upper_positions[0][0])
        self.lower_support = self.ss.find_node_id(lower_positions[0][0])

        # # Add support.
        self.ss.add_support_fixed(self.upper_support)
        self.ss.add_support_fixed(self.lower_support)

    def show(self):
        self.ss.show_structure()

    def load(self, fx=0, fy=-1):
        self.ss.point_load(self.tip_node, fx, fy)

    def solve(self):
        self.ss.solve()

    def metric(self):
        displacement = self.bond_in_coordinate[:, 1, :] - self.bond_in_coordinate[:, 0, :]
        lengths = np.linalg.norm(displacement, axis=-1)
        forces = np.array([a["N"] for a in self.ss.get_element_results(0)])

        print(lengths)
        print(forces)

        return np.sum(np.abs(lengths * forces))


def _generate_lower_positions(arr):
    return [[[n[0], n[1] * -1] for n in a] for a in arr]


def _flatten(l_):
    return [item for sl in l_ for item in sl]


def _flatten_and_remove_duplicate(upper, lower):
    top_fix = upper[0][0]
    upper = [a[1:] for a in upper]

    lower_fix = lower[0][0]
    lower = [a[1:] for a in lower]

    # Remove center from lower, since already included in upper.
    lower = [a[:-1] for a in lower]

    # flatten.
    upper_flat = [n for a in upper for n in a]
    lower_flat = [n for a in lower for n in a]
    return [top_fix] + upper_flat + [lower_fix] + lower_flat


def flatten_and_remove_duplicate(upper_positions, ur, lower_positions, lr):
    pos = _flatten_and_remove_duplicate(upper_positions, lower_positions)
    references = _flatten_and_remove_duplicate(ur, lr)
    return pos, references


def assign_node_reference(
        nodes, upper_fix_ref=0, lower_fix_ref=1, start_number=2):
    """Creates node number assignment."""
    a = 2
    upper_references = []
    lower_references = []
    for sheet in nodes:
        ur = [upper_fix_ref]
        lr = [lower_fix_ref]

        # Remove support node.
        sheet = sheet[1:]

        if len(sheet) > 1:
            for n in sheet[:-1]:
                # add to uper reference.
                ur.append(a)
                a = a + 1

                # add to lower_reference.
                lr.append(a)
                a = a + 1

        # Add center node to both.
        ur.append(a)
        lr.append(a)
        a = a + 1

        upper_references.append(ur)
        lower_references.append(lr)

    return upper_references, lower_references


def connect_radial(upper_references, lower_references):
    """Describes radial connections."""
    bonds = []
    if len(upper_references) > 1:
        for inner, outer in zip(upper_references[:-1], upper_references[1:]):
            for n1, n2 in zip(inner[1:], outer[1:-1]):
                bonds.append([n1, n2])
        for inner, outer in zip(lower_references[:-1], lower_references[1:]):
            for n1, n2 in zip(inner[1:], outer[1:-1]):
                bonds.append([n1, n2])
    return bonds


def _connect_linear(row):
    return [[a, b] for a, b in zip(row[1:], row[:-1])]


def connect_extension(upper_references, lower_references):
    """Describes extension connections."""
    bonds = []
    bonds += [_connect_linear(r) for r in upper_references]
    bonds += [_connect_linear(r) for r in lower_references]
    return bonds


def build_connections(upper_references, lower_references):
    """stuff"""
    rad = connect_radial(upper_references, lower_references)
    ext = _flatten(connect_extension(upper_references, lower_references))
    return rad + ext


def analyze(surface_node):
    """Analyzes a michell truss parameterized by the upper set of nodes.

    Args:
        nodes: List of nodes in subsequent layers of michell truss.

    Returns:
    """
