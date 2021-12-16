"""Creates Flask server for calling `fdm.py` from Grasshopper."""


from flask import Flask
import ghhops_server as hs

import fdm
import re

import rhino3dm


# register hops app as middleware
app = Flask(__name__)
hops = hs.Hops(app)


def _parse_tree_key(k):
    """Extracts coordinates from key"""
    # Extract coordinates.
    k = re.findall('{(.*?)}', k)[0]

    # Split by `;` to get individual indices.
    k = k.split(r';')

    return [int(k_) for k_ in k]


def _insert_at(mod_list, coor, data):
    """Inserts `data` into `mod_list` at `coor`.

    `coor` in general may be a list of coordinates. This function then
    recursively inserts lists until the appropriate level is reached.
    """
    if len(mod_list) < coor[0]:
        raise ValueError("List too short.")
    # TODO: Need to make clear whether this will ALWAYS add new items or
    # add to existing lists if they exis.
    if len(coor) > 1:
        # If list doesn't exist, insert one.
        if len(mod_list) <= coor[0]:
            mod_list.insert(coor[0], [])
        if not isinstance(mod_list[coor[0]], list):
            raise ValueError(
                "Cannot insert data at requested level. Not a list.")
        _insert_at(mod_list[coor[0]], coor[1:], data)
    else:
        mod_list.insert(coor[0], data)


def list_from_tree(tree):
    """Constructs list of list containing data in `tree`.

    It is useful to construct a list representation first, since different
    final indices may contain different.
    """
    address = [_parse_tree_key(k) for k in tree.keys()]
    data = tree.values()

    # Instantiate list which we will fill in whith data at `address`.
    lt = []

    [_insert_at(lt, a, d) for a, d in zip(address, data)]
    return lt


def _point_to_numpy(points):
    """Generates list of coordinates from list of `rhino3dm.point`."""
    return [[p.X, p.Y, p.Z] for p in points]


@hops.component(
    "/fdm",
    name="FDM",
    description="Computes final node positions using Force Density Method.",
    inputs=[
        hs.HopsBoolean("Run", "Bool", "Runs Computation."),
        hs.HopsPoint("Nodes", "N", access=hs.HopsParamAccess.LIST),
        hs.HopsNumber("Edges", "E", access=hs.HopsParamAccess.TREE),
        hs.HopsVector("Loads", "L", access=hs.HopsParamAccess.LIST),
        hs.HopsNumber("Force Densities", "FD", access=hs.HopsParamAccess.TREE),
        hs.HopsInteger("Fixed Node ID", "Fixed",
                       access=hs.HopsParamAccess.TREE)
        ],
    outputs=[
        hs.HopsPoint("Nodes", "N", access=hs.HopsParamAccess.TREE),
        # hs.HopsLine("Edges", "E", access=hs.HopsParamAccess.TREE)
        ]
    )
def hops_fdm(trigger, nodes, edges, loads, force_density, fixed):
    """Process hops/grasshopper objects and returns fdm solution."""
    if trigger:
        # Parse into numpy array.
        nodes = _point_to_numpy(nodes)
        print(f"nodes {nodes}")

        edges = list_from_tree(edges)
        edges = [list(map(int, e)) for e in edges]
        print(f"edges {edges}")

        loads = _point_to_numpy(loads)
        print(f"loads {loads}")

        force_densities = list_from_tree(force_density)[0]
        print(f"force_densities {force_densities}")

        fixed = list_from_tree(fixed)[0]
        print(f"fixed {fixed}")

        a = fdm.compute_fdm(nodes, edges, loads, force_densities, fixed)

        # Convert to points.
        points = [rhino3dm.Point3d(*p) for p in a]

        # Create new mesh edges.
        edges = [rhino3dm.Line(points[e[0]], points[e[1]]) for e in edges]

        print(edges)

        return points


def run_app():
    app.run(threaded=True)


if __name__ == "__main__":
    run_app()
