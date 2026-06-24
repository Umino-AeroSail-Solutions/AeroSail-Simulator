import numpy as np

class Node:
    def __init__(self, node_id, node_coords):
        self.id = node_id
        self.coords = node_coords

class Element:
    # x is to the right, y is upwards
    def __init__(self, id, node_i, node_j, E, A, I):
        self.id = id
        self.E = E
        self.A = A
        self.I = I
        self.node_i = node_i
        self.node_j = node_j
        self.L = np.linalg.norm(node_j.coords - node_i.coords)
        self.angle = np.arctan2(node_j.coords[1]-node_i.coords[1], node_j.coords[0]-node_i.coords[0])

        self.forces = np.zeros(6)
        # precompute and store
        self.T = self.transformation_matrix()
        self.k_local = self.local_stiffness()
        self.k_global = self.T.T @ self.k_local @ self.T


    def local_stiffness(self):
        E = self.E
        A = self.A
        I = self.I
        L = self.L
        return np.array([
            [E*A/L, 0, 0, -E*A/L, 0, 0],
            [0, 12*E*I/L**3, 6*E*I/L**2, 0, -12*E*I/L**3, 6*E*I/L**2],
            [0, 6*E*I/L**2, 4*E*I/L, 0, -6*E*I/L**2, 2*E*I/L],
            [-E*A/L, 0, 0, E*A/L, 0, 0],
            [0, -12*E*I/L**3, -6*E*I/L**2, 0, 12*E*I/L**3, -6*E*I/L**2],
            [0, 6*E*I/L**2, 2*E*I/L, 0, -6*E*I/L**2, 4*E*I/L]
        ])

    def transformation_matrix(self):
        angle = self.angle
        return np.array([
            [np.cos(angle), -np.sin(angle), 0, 0, 0, 0],
            [np.sin(angle), np.cos(angle), 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0],
            [0, 0, 0, np.cos(angle), -np.sin(angle), 0],
            [0, 0, 0, np.sin(angle), np.cos(angle), 0],
            [0, 0, 0, 0, 0, 1]
        ])


    def element_forces(self, u_global_element):
        return np.transpose(self.T) @ self.k_local @ self.T @ u_global_element