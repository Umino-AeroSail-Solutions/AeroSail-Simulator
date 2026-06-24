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
        for node, forces in force_dict_in:
            for i in range(3):
                self.f[node+i] = forces[i]
    
    def apply_boundary_conditions(self, bc_list):
        bc_indexes = []
        for i, node in enumerate(bc_list):
            for j, bc in enumerate(node):
                if bc == 1:
                    bc_indexes.append(int(node)+j)
        
        f_free = np.copy(self.f)
        K_free = np.copy(self.big_K)
        for index in bc_index:
            f_free = np.delete(f_fixed, index)
            K_free = np.delete(K_free, index, 0)
            K_free = np.delete(K_free, index, 1)
        
        self.f_f = f_free
        self.big_K_f = big_K_f
        

        

    
    def solve(self):
        u_sol = np.linalg.solve(self.big_K_f, self.f_f)

            
