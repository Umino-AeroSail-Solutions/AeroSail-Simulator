from railteam.frame_fea.assembly import Assembly
from railteam.frame_fea.element import Node, Element
import numpy as np


# unit test case, cantilever beam
node1 = Node(0, np.array([0,0]))
node2 = Node(1, np.array([1,0]))

bc_list = {
    0   :   [1,1,1],
    1   :   [0,0,0]
}

force_dict = {
    0   :   [0,0,0],
    1   :   [0,-50,0]
}


element1 = Element(0, node1, node2, 70*10**9, 0.01**2, 1/12*0.01**4)
element_list = [element1]
node_list = [0, 1]

assembly = Assembly(element_list, node_list)
assembly.construct_big_K()
assembly.build_force_vector(force_dict)
assembly.apply_boundary_conditions(bc_list)
assembly.solve()
assembly.give_forces()

print(element1.forces)
print(assembly.u_sol)