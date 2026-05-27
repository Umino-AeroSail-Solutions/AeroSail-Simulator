
def local_stiffness(E, A, I, L):
    np.array(
        [E*A/L, 0, 0, -E*A/L, 0, 0],
        [0, 12*E*I/L**3, 6*E*I/L**2, 0, -12*E*I/L**3, 6*E*I/L**2],
        [0, 6*E*I/L**2, 4*E*I/L, 0, -6*E*I/L**2, 2*E*I/L],
        [-E*A/L, 0, 0, E*A/L, 0, 0],
        [0, -12*E*I/L**3, -6*E*I/L**2, 0, 12*E*I/L**3, -6*E*I/L**2],
        [0, 6*E*I/L**2, 2*E*I/L, 0, -6*E*I/L**2, 4*E*I/L],
    )

def transformation_matrix(angle):
    np.array(
        [np.cos(angle), -np.sin(angle), 0, np.cos(angle), -np.sin(angle), 0],
        [np.sin(angle), np.cos(angle), 0, np.sin(angle), np.cos(angle), 0],
        [0, 0, 1, 0, 0, 1],
        [np.cos(angle), -np.sin(angle), 0, np.cos(angle), -np.sin(angle), 0],
        [np.sin(angle), np.cos(angle), 0, np.sin(angle), np.cos(angle), 0],
        [0, 0, 1, 0, 0, 1]
    )


def element_forces(u_global_element, E, A, I, node_i, node_j):
