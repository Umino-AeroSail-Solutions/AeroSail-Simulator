import json
from railteam.frame_fea.assembly import Assembly
from railteam.frame_fea.element import Node, Element
import numpy as np

def load_from_json(filepath):
    with open(filepath) as f:
        data = json.load(f)
    
    # Build nodes
    nodes = {int(k): Node(int(k), np.array(v)) 
             for k, v in data["nodes"].items()}
    
    # Build elements
    props = data["properties"]
    elements = []
    for e in data["elements"]:
        ni = nodes[e["nodes"][0]]
        nj = nodes[e["nodes"][1]]
        p = props[e["prop"]]
        elements.append(Element(e["id"], ni, nj, p["E"], p["A"], p["I"]))
    
    bc_list = {int(k): v for k, v in data["boundary_conditions"].items()}
    force_dict = {int(k): v for k, v in data["loads"].items()}
    node_list = list(nodes.values())

    return elements, node_list, bc_list, force_dict

def run(filepath):
    elements, node_list, bc_list, force_dict = load_from_json(filepath)
    assembly = Assembly(elements, node_list)
    assembly.construct_big_K()
    assembly.build_force_vector(force_dict)
    assembly.apply_boundary_conditions(bc_list)
    assembly.solve()
    assembly.give_forces()
    return assembly

if __name__ == "__main__":
    assembly = run("railteam/frame_fea/cantilever.json")
    print(assembly.u_sol)
    for e in assembly.element_list:
        print(f"Element {e.id}: {e.forces}")