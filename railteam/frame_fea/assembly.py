import numpy as np


class Assembly:

    def __init__(self, element_list, node_list):
        self.element_list = element_list
        self.node_list = node_list

        node_len = len(node_list)
        self.big_K = np.zeros((node_len*3, node_len*3))
        self.f      = np.zeros(node_len*3)
    

    def construct_big_K(self):
        for element in self.element_list:
            # split k_global into it's constituent parts, and allocate them into the big k based on node numbers
            matrix_spot_1 = element.node_i.id*3
            matrix_spot_2 = element.node_j.id*3
            self.big_K[matrix_spot_1:matrix_spot_1+3, matrix_spot_1:matrix_spot_1+3] += element.k_global[0:3, 0:3]
            self.big_K[matrix_spot_1:matrix_spot_1+3, matrix_spot_2:matrix_spot_2+3] += element.k_global[0:3, 3:6]
            self.big_K[matrix_spot_2:matrix_spot_2+3, matrix_spot_1:matrix_spot_1+3] += element.k_global[3:6, 0:3]
            self.big_K[matrix_spot_2:matrix_spot_2+3, matrix_spot_2:matrix_spot_2+3] += element.k_global[3:6, 3:6]

    def build_force_vector(self, force_dict_in):
        for node, forces in force_dict_in.items():
            for i in range(3):
                self.f[node*3+i] = forces[i]
    
    def apply_boundary_conditions(self, bc_list):
        bc_indices = []
        for node, bcs in bc_list.items():
            for j, bc in enumerate(bcs):
                if bc == 1:
                    bc_indices.append(int(node)*3+j)
        
        f_free = np.copy(self.f)
        K_free = np.copy(self.big_K)
        f_free = np.delete(f_free, bc_indices)
        K_free = np.delete(K_free, bc_indices, 0)
        K_free = np.delete(K_free, bc_indices, 1)
        
        self.bc_indices = np.sort(bc_indices)
        
        self.f_f = f_free
        self.big_K_f = K_free


        
    def solve(self):
        u_sol_f = np.linalg.solve(self.big_K_f, self.f_f)
        u_sol = u_sol_f
        for index in self.bc_indices:
            u_sol = np.insert(u_sol, index, 0)
        
        self.u_sol = u_sol

    
    def give_forces(self):
        for element in self.element_list:
            u_spot_1 = element.node_i.id*3
            u_spot_2 = element.node_j.id*3
            u_global_element = np.hstack((self.u_sol[u_spot_1:u_spot_1+3], self.u_sol[u_spot_2:u_spot_2+3]))

            element.forces = element.element_forces(u_global_element)
