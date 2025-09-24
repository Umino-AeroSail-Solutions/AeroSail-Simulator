import numpy as np
import matplotlib.pyplot as plt

class Material(object):
    def __init__(self, material_name, EModulus=None, YieldStrength=None, Poisson=None, Density=None, max_shear=None):
        self.material_name = material_name
        self.EMod = EModulus
        self.Yield = YieldStrength
        self.Poisson = Poisson
        self.Density = Density
        self.max_shear = max_shear

class Square_beam(object):
    def __init__(self, length, w=None, h=None, t=None, material=None):
        self.length = length
        self.w = w
        self.h = h
        self.t = t
        self.material = material
        self.EMod = material.EMod
        self.Yield = material.Yield
        self.max_shear_yield = material.max_shear
    def get_Ixx_Iyy_Ixy(self):
        self.b = self.w - 2*self.t
        self.d = self.h - 2*self.t
        self.Ixx = (self.w*(self.h**3)/12) - (self.b*(self.d**3)/12)
        self.Iyy = (self.h * (self.w ** 3) / 12) - (self.d * (self.b ** 3) / 12)
        self.Ixy = 0
    def get_bending_tension(self, moment1, moment2, x, y):
        tension = ((moment1*self.Iyy-moment2*self.Ixy)*y + (moment2*self.Ixx-moment1*self.Ixy)*x)/(self.Ixx*self.Iyy-(self.Ixy**2))
        return tension
    def gettensionSF(self, moment1, moment2, tension=0):
        corner_coordinates = [[-self.w/2, self.h/2],[self.w/2, self.h/2],[self.w/2, -self.h/2],[-self.w/2, -self.h/2]]
        tensions = []
        for coordinates in corner_coordinates:
            tensions.append(self.get_bending_tension(moment1, moment2, *coordinates)+tension)
        tensions = np.array(tensions)
        minSF = self.Yield / np.max(np.abs(tensions))
        return minSF

    def get_shear(self, shear1, shear2, loc):
        Vx = shear1
        Vy = shear2
        self.b = self.w - 2 * self.t
        self.d = self.h - 2 * self.t

        qb_1 = lambda s: -Vy / self.Ixx * self.t * (-self.h / 2 * s + s ** 2 / 2) - Vx / self.Iyy * self.t * (
                    self.b / 2) * s
        qb_2 = lambda s: -Vy / self.Ixx * self.t * (self.h / 2) * s - Vx / self.Iyy * self.t * (
                    self.b / 2 * s - s ** 2 / 2)
        qb_3 = lambda s: -Vy / self.Ixx * self.t * (self.h / 2 * s - s ** 2 / 2) - Vx / self.Iyy * self.t * -(
                    self.b / 2) * s
        qb_4 = lambda s: -Vy / self.Ixx * self.t * (-self.h / 2) * s - Vx / self.Iyy * self.t * (
                    -self.b / 2 * s + s ** 2 / 2)

        # # Corrected qs0
        # qs0 = torsion_shear / (2 * ((self.w + self.b) / 2) * ((self.h + self.d) / 2))
        # Wroooooong INTEGRATION CONSTANTS

        # Corrected perimeter wrap
        perimeter = 2 * self.h + 2 * self.b
        loc = loc % perimeter

        if 0 <= loc <= self.h:
            return qb_1(loc)
        elif self.h < loc <= self.h + self.b:
            loc -= self.h
            return qb_2(loc) + qb_1(self.h)
        elif self.h + self.b < loc <= 2 * self.h + self.b:
            loc -= self.h + self.b
            return qb_3(loc) + qb_2(self.b) + qb_1(self.h)
        elif 2 * self.h + self.b< loc <= perimeter:
            loc -= 2 * self.h + self.b
            return qb_4(loc) + qb_3(self.h) + qb_2(self.b) + qb_1(self.h)

    def create_shearflow_list(self, shear1, shear2, side_subdivisions=80):
        # Perimeter for shear sampling
        P = 2 * self.h + 2 * self.b

        # Compute initial shear safety factor
        shearflows = []
        for loc in np.linspace(0, P, side_subdivisions, endpoint=False):
            shearflows.append(self.get_shear(shear1, shear2, loc))
        self.shearflows = np.array([shearflows, np.linspace(0, P, side_subdivisions, endpoint=False)]).T
        # print(self.shearflows)
        return self.shearflows

    def adjust_for_torsion(self, shear1, shear2, torsion=0, side_subdivisions=80):
        self.b = self.w - 2 * self.t
        self.d = self.h - 2 * self.t
        P = 2 * self.h + 2 * self.b
        N = side_subdivisions
        ds = P / N

        self.create_shearflow_list(shear1, shear2, side_subdivisions)

        total_inner_qs0 = 0

        for q, loc in self.shearflows:
            if 0 <= loc <= self.h:
                lever_arm = self.b / 2 # left wall (vertical)
            elif self.h < loc <= self.h + self.b:
                lever_arm = self.h / 2  # bottom wall (horizontal)
            elif self.h + self.b < loc <= 2 * self.h + self.b:
                lever_arm = self.b / 2  # right wall
            elif 2 * self.h + self.b < loc <= P:
                lever_arm = self.h / 2  # top wall
            else:
                continue  # shouldn't happen, but safe

            total_inner_qs0 += q * lever_arm * ds

        print("Total torsional moment (from q):", total_inner_qs0)

        # If you want to correct for torsion (optional):
        A = self.h*self.b
        q0 = total_inner_qs0 / (2 * A) 
        self.shearflows[:, 0] -= q0

        # Plot
        # plt.plot(self.shearflows[:, 1], self.shearflows[:, 0])
        # corners = [0, self.h, self.h+self.b, 2*self.h+self.b, 2*self.h+2*self.b]
        # plt.vlines(corners, min(self.shearflows[:,0]), max(self.shearflows[:,0]), linestyles="dotted")
        # plt.xlabel("Location along perimeter, s")
        # plt.ylabel("Shear flow q(s)")
        # plt.title("Shear‐flow distribution around beam perimeter")
        # plt.tight_layout()
        # plt.show()


    # def get_shears(self, shear1, shear2, torsion_shear, loc):
    #     # Great shit anhong what the fuck is self.d and self.b it was not defined :(
    #     Vx = shear1
    #     Vy = shear2
    #     self.b = self.w - 2 * self.t
    #     self.d = self.h - 2 * self.t
    #     #     w
    #     # ||=======|      x
    #     # ||   b   |    <---+ (not the shear center, just for coordinate directions)
    #     # ||     d | h      |
    #     # ||       |        | y
    #     # ||-------|        v
    #     # side 1 || s: 0 --> h
    #     qb_1 = lambda s: -Vy/self.Ixx * self.t * (-self.h/2*s + (s**2) /2) - Vx/self.Iyy * self.t * ((self.w+self.b)/4)*s
    #     # side 2 -  s: 0 --> b
    #     qb_2 = lambda s: -Vy/self.Ixx * self.t * ((self.h+self.d)/4)*s - Vx/self.Iyy * self.t * (self.b/2*s - (s**2)/2)
    #     # side 3 |  s: 0 --> h
    #     qb_3 = lambda s: -Vy/self.Ixx * self.t * (self.h/2*s - (s**2) /2) - Vx/self.Iyy * self.t * -((self.w+self.b)/4)*s
    #     # side 4 _  s: 0 --> b
    #     qb_4 = lambda s: -Vy/self.Ixx * self.t * (-(self.h+self.d)/4)*s - Vx/self.Iyy * self.t * (-self.b/2*s + (s**2)/2)
    #     # qs0
    #     qs0 = torsion_shear/2/(((self.w+self.b)/2)*((self.h+self.d)/2))
    #
    #     # what's below could have been a nested-if but i made it like this for readability :3 (oh my god wtf is that face ew furry ew what)
    #     loc = loc % 2*self.h + 2*self.b
    #     if loc >= 0 and loc <= self.h:
    #         return qb_1(loc) + qs0
    #     elif loc > self.h and loc <= (self.h + self.b):
    #         loc = loc - self.h
    #         return qb_2(loc) + qb_1(self.h) + qs0
    #     elif loc > (self.h + self.b) and loc <= (2*self.h + self.b):
    #         loc = loc - self.h - self.b
    #         return qb_3(loc) + qb_2(self.b) + qb_1(self.h) + qs0
    #     elif loc > (2*self.h + self.b) and loc <= (2*self.h + 2*self.b):
    #         loc = loc - 2*self.h - self.b
    #         return qb_4(loc) + qb_3(self.h) + qb_2(self.b) + qb_1(self.h) + qs0

    def size_profile(self, shear1, shear2, tension, moment1, moment2, SF,
                     initial_thickness, side_subdivisions=50, plot_convergence=False):
        """
        Iteratively adjust thickness to meet a target safety factor based on shear and tension.
        If plot_convergence is True, displays how thickness and safety factors evolve each iteration.
        """
        # Check required parameters
        if self.w is None or self.h is None or self.material is None \
                or self.EMod is None or self.Yield is None:
            raise ValueError(
                "Not enough fixed parameters, ensure Materials, W, H and EMod are set and the material has a yield"
            )

        # Initialize thickness and inertia
        self.t = initial_thickness
        self.get_Ixx_Iyy_Ixy()

        # Track history
        if plot_convergence:
            last_ts = [self.t]
            tension_sfs = []
            shear_sfs = []
            actual_sfs = []

        # Compute initial safety factors
        self.tensionSF = self.gettensionSF(moment1, moment2, tension)
        self.adjust_for_torsion(shear1, shear2)
        max_shear = np.max(np.abs(self.shearflows[:, 0]/self.t))
        self.shearSF = self.max_shear_yield / max_shear
        self.actualSF = min(self.shearSF, self.tensionSF)
        design_sf = self.actualSF / SF

        if plot_convergence:
            tension_sfs.append(self.tensionSF)
            shear_sfs.append(self.shearSF)
            actual_sfs.append(self.actualSF)

        # Iteratively adjust thickness until within [1.0, 1.1] of target SF
        iteration = 0
        print(f"[Iteration {iteration}] Thickness: {self.t:.6f} m, Tension SF: {self.tensionSF:.3f}, "
              f"Shear SF: {self.shearSF:.3f}, Overall SF: {self.actualSF:.3f}")
        while design_sf < 1 or design_sf > 1.1:
            iteration += 1
            self.t = self.t / design_sf
            if plot_convergence:
                last_ts.append(self.t)

            self.get_Ixx_Iyy_Ixy()
            self.tensionSF = self.gettensionSF(moment1, moment2, tension)

            self.adjust_for_torsion(shear1, shear2)
            max_shear = np.max(np.abs(self.shearflows[:, 0]/self.t))
            self.shearSF = self.max_shear_yield / max_shear

            self.actualSF = min(self.shearSF, self.tensionSF)
            design_sf = self.actualSF / SF

            if plot_convergence:
                tension_sfs.append(self.tensionSF)
                shear_sfs.append(self.shearSF)
                actual_sfs.append(self.actualSF)

            print(f"[Iteration {iteration}] Thickness: {self.t:.6f} m, Tension SF: {self.tensionSF:.3f}, "
                  f"Shear SF: {self.shearSF:.3f}, Overall SF: {self.actualSF:.3f}")

            if iteration > 1000:
                raise RuntimeError("Thickness iteration did not converge after 1000 steps")

        # Plot convergence if requested
        if plot_convergence:
            fig, axs = plt.subplots(2, 1, figsize=(8, 6), sharex=True)

            axs[0].plot(range(len(last_ts)), last_ts, marker='o')
            axs[0].set_ylabel('Thickness t [m]')
            axs[0].grid(True)

            axs[1].plot(range(len(tension_sfs)), tension_sfs, label='Tension SF', linestyle='--')
            axs[1].plot(range(len(shear_sfs)), shear_sfs, label='Shear SF', linestyle='-.')
            axs[1].plot(range(len(actual_sfs)), actual_sfs, label='Overall SF', color='black')
            axs[1].axhline(SF, color='gray', linestyle=':', label='Target SF')
            axs[1].set_ylabel('Safety Factors')
            axs[1].set_xlabel('Iteration')
            axs[1].legend()
            axs[1].grid(True)

            plt.suptitle('Convergence of Thickness and Safety Factors')
            plt.tight_layout()
            plt.show()

        return self.t
    def compute_reaction_loads(self, loads, use_self_weight=np.array([0, 0, 0])):
        self.mass = (self.w*self.h-(self.w-2*self.t)*(self.h-2*self.t))*self.material.Density
        weight_load = [self.mass*use_self_weight[0], self.mass*use_self_weight[1], self.mass*use_self_weight[2], self.length/2]
        self.loads=loads
        simplified_loads = self.loads
        simplified_loads.append(weight_load)
        self.total_load = [0,0,0]
        unweighted_average_position = [0, 0]
        for load in simplified_loads:
            self.total_load += [load[0], load[1], load[2]]
            unweighted_average_position += [load[3]*load[0], load[3]*load[1]]
        self.total_load = np.array(self.total_load)
        self.average_position = np.array(unweighted_average_position)/np.sum(np.array(simplified_loads)[:,3])

        M_left = [self.total_load[0]*self.average_position[0], self.total_load[1]*self.average_position[1]]
        M_right = [self.total_load[0]*(self.length-self.average_position[0]), self.total_load[0]*(self.length-self.average_position[3])]

        self.rload_right = [M_left[0]/self.length, M_left[1]/self.length, 0, self.length]
        self.rload_left =[M_right[0]/self.length, M_right[1]/self.length, self.total_load[2], 0]
        self.loads.append(self.rload_left)
        self.loads.append(self.rload_right)
        print("Reaction loads: ", [self.rload_right, self.rload_left])

    def compute_internal_loads(self, loads, use_self_weight=np.array([0, 0, 0]), length_subdivisions=30):
        # first compute reactions and seed self.loads with your point‐loads + reactions
        self.compute_reaction_loads(loads, use_self_weight)
        # now append self‐weight as a series of small point‐loads,
        # each with [wx, wy, wz, x_position]
        for i in range(1, length_subdivisions + 1):
            frac = i / length_subdivisions
            x_i = frac * self.length
            wx, wy, wz = (self.mass * use_self_weight).tolist()  # now a 3‑tuple
            self.loads.append([wx, wy, wz, x_i])  # now 4 elements

        # convert to array and sort _rows_ by x (column 3)
        arr = np.array(self.loads)  # shape (N,4)
        idx = arr[:, 3].argsort()  # indices that sort by column 3
        self.loads_array = arr[idx, :]  # now a properly sorted (N×4) array

        # build cumulative internal loads…
        internal = [[0.0, 0.0, 0.0, 0.0]]
        for Vx, Vy, N, x in self.loads_array:
            pVx, pVy, pN, _ = internal[-1]
            internal.append([pVx + Vx, pVy + Vy, pN + N, x])
        self.internal_loads = np.array(internal)
        # …and so on for moments

        # build cumulative internal moments
        moments = [[0.0, 0.0, 0.0]]  # [Mx, My, Tz] at x=0

        # loop over each station in internal_loads
        for i, (Vx_i, Vy_i, N_i, xi) in enumerate(self.internal_loads):
            Mx = 0.0
            My = 0.0
            # sum contributions from all loads to the right of xi
            for Vx_j, Vy_j, N_j, xj in self.internal_loads[i:]:
                dx = xj - xi
                Mx += Vy_j * dx  # bending moment about x‐axis due to vertical shear
                My += Vx_j * dx  # bending moment about y‐axis due to horizontal shear
            moments.append([Mx, My, xi])

        # convert to numpy array (shape will match internal_loads)
        self.internal_moments = np.array(moments[1:])
        print("Moments: ", self.internal_moments)
        print("Shears: ", self.internal_loads)

    def thickness_simple_yield_size(self, loads, SF, initial_thickness, use_self_weight=np.array([0,0,0]),
                                    length_subdivisions=30):
        '''loads are defined as [shear1, shear2, normal tension, position across length], use self weight is an
        acceleration vector [x, y, z] where x is shear 1, y is shear 2 and z is normal tension'''
        self.t = initial_thickness
        self.best_t = self.t
        last_weight = self.mass
        tensionsfs = []
        shearsfs = []
        index = 0
        while abs(self.mass-last_weight) > 1e-3:
            for load in self.internal_loads:
                moment = self.internal_moments[index,:]
                self.compute_internal_loads(loads, use_self_weight, length_subdivisions=length_subdivisions)
                self.size_profile(load[0], load[1], load[2], moment[0], moment[1], moment[2], self.t)
                if self.t>self.best_t:
                    self.best_t = self.t
                tensionsfs.append(self.tensionSF)
                shearsfs.append(self.shearSF)


        print("Final thickness is: ", self.t, " meters")
        return self.t

    def plot_shear_heatmap(self, shear1, shear2, torsion_shear=0.0, subdivisions=300):
        """
        Plot a heatmap of shear flow around the beam perimeter and overlay the resultant shear vector.

        Parameters:
          shear1         – Vx, shear in x‐direction
          shear2         – Vy, shear in y‐direction
          torsion_shear  – total torsional shear flow (default 0)
          subdivisions   – number of points around perimeter to sample (default 300)
        """
        # Recompute geometry and inertias
        self.get_Ixx_Iyy_Ixy()
        self.b = self.w - 2 * self.t
        self.d = self.h - 2 * self.t

        # Total mid‐line perimeter
        P = 2 * self.h + 2 * self.b
        locs = np.linspace(0, P, subdivisions, endpoint=False)

        xs, ys, qs = [], [], []
        for loc in locs:
            # Wrap loc into [0, P)
            s = loc % P
            # Map s to (x, y) on mid‐line
            if s <= self.h:
                x = -self.w / 2
                y = self.h / 2 - s
            elif s <= self.h + self.b:
                x = -self.w / 2 + (s - self.h)
                y = -self.h / 2
            elif s <= 2 * self.h + self.b:
                x = self.w / 2
                y = -self.h / 2 + (s - self.h - self.b)
            else:
                x = self.w / 2 - (s - 2 * self.h - self.b)
                y = self.h / 2
            xs.append(x)
            ys.append(y)
            qs.append(self.get_shears(shear1, shear2, torsion_shear, loc))

        # Plot shear flow heatmap
        plt.scatter(xs, ys, c=qs)
        plt.colorbar(label='Shear flow q(s)')

        # Outline of the beam cross‐section
        beam_outline = plt.Rectangle(
            (-self.w / 2, -self.h / 2), self.w, self.h,
            fill=False, linewidth=1.5
        )
        plt.gca().add_patch(beam_outline)

        # Draw resultant shear vector at centroid
        # Normalize arrow to 80% of half‐width/height
        total = abs(shear1) + abs(shear2) if (abs(shear1) + abs(shear2)) != 0 else 1
        dx = (shear1 / total) * (self.w / 2) * 0.8
        dy = (shear2 / total) * (self.h / 2) * 0.8
        plt.arrow(0, 0, dx, dy, head_width=0.02 * max(self.w, self.h), length_includes_head=True)

        plt.gca().set_aspect('equal', 'box')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('Shear‐flow heatmap and resultant shear vector')
        plt.tight_layout()
        plt.show()

    def plot_shearflow(self, shear1, shear2, torsion_shear=0.0, subdivisions=200):
        """
        Compute and plot the shear‐flow q(s) around the full perimeter of the thin‐walled square beam.

        Parameters:
          shear1         – Vx, shear in x‑direction
          shear2         – Vy, shear in y‑direction
          torsion_shear  – total torsional shear flow (default 0)
          subdivisions   – number of points around perimeter to sample (default 200)
        """
        # Ensure inertia and geometry are up‑to‑date
        self.get_Ixx_Iyy_Ixy()
        # recompute inner dims
        self.b = self.w - 2 * self.t
        self.d = self.h - 2 * self.t

        # total mid‑line perimeter
        P = 2 * self.h + 2 * self.b
        # sample along [0, P)
        locs = np.linspace(0, P, subdivisions, endpoint=False)

        # compute q at each location
        q_vals = [self.get_shears(shear1, shear2, torsion_shear, loc) for loc in locs]

        # make the plot
        plt.plot(locs, q_vals)
        # plt.scatter(locs, q_vals)
        plt.xlabel("Location along perimeter, s")
        plt.ylabel("Shear flow q(s)")
        plt.title("Shear‐flow distribution around beam perimeter")
        plt.tight_layout()
        plt.show()

    def plot_NVM(self, loads, use_self_weight=np.array([0, 0, 0]), length_subdivisions=100):
        """
        Plot Shear‐Force (V) and Bending‐Moment (M) diagrams with reaction forces included.

        Parameters:
          loads               : list of [Vx, Vy, N, x_pos] applied loads
          use_self_weight     : 3‐vector [ax, ay, az] for distributed self‐weight
          length_subdivisions : number of slices to discretize self‐weight
        """
        # 1) Compute internal forces and moments
        self.compute_internal_loads(loads, use_self_weight, length_subdivisions)

        # 2) Extract arrays
        x = self.internal_loads[:, 3]
        V = self.internal_loads[:, 1]  # Shear force (Vy)
        M = self.internal_moments[:, 0]  # Moment (My)

        if len(x) != len(M):
            raise ValueError(f"Length mismatch: internal_loads has {len(x)} points, internal_moments has {len(M)}")

        # 3) Plot
        fig, (axV, axM) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)

        # --- Shear force diagram ---
        axV.step(x, V, where='post', label='Shear Force V')
        axV.set_ylabel('Shear V [N]')
        axV.set_title('Shear‐Force Diagram')
        axV.grid(True)

        # Plot reaction shear vectors
        for rload in [self.rload_left, self.rload_right]:
            x_pos = rload[3]
            Vy = rload[1]
            axV.annotate('', xy=(x_pos, Vy), xytext=(x_pos, 0),
                         arrowprops=dict(facecolor='red', shrink=0.05, width=1.5, headwidth=8),
                         fontsize=8)
            axV.text(x_pos, Vy * 1.05, f'R={Vy:.0f}', color='red', ha='center', fontsize=8)

        # --- Moment diagram ---
        axM.plot(x, M, label='Moment M', color='tab:blue')
        axM.set_ylabel('Moment M [N·m]')
        axM.set_xlabel('Position along beam [m]')
        axM.set_title('Bending‐Moment Diagram')
        axM.grid(True)

        plt.tight_layout()
        plt.show()


# Testing code

def test_inertia():
    mat = Material("Steel", EModulus=210e9, YieldStrength=250e6, max_shear=150e6)
    beam = Square_beam(length=1.0, w=0.2, h=0.15, t=0.005, material=mat)
    beam.get_Ixx_Iyy_Ixy()
    # For a thin wall, Ixx should be positive and on the order of 1e-6 to 1e-5 m^4
    assert beam.Ixx > 0
    assert beam.Iyy > 0
    print(f"Ixx = {beam.Ixx:.3e}, Iyy = {beam.Iyy:.3e}")

def test_tension_SF():
    mat = Material("Aluminum", EModulus=70e9, YieldStrength=300e6, max_shear=200e6)
    beam = Square_beam(length=2.0, w=0.1, h=0.1, t=0.01, material=mat)
    beam.get_Ixx_Iyy_Ixy()
    # Pure bending about x-axis: moment1=1000 Nm, moment2=0, no axial tension
    sf = beam.gettensionSF(moment1=1000, moment2=0, tension=0)
    assert sf > 0
    print(f"Tension safety factor = {sf:.2f}")

def test_shearflows():
    mat = Material("Composite", EModulus=50e9, YieldStrength=500e6, max_shear=300e6)
    beam = Square_beam(length=0.5, w=0.12, h=0.08, t=0.008, material=mat)
    beam.get_Ixx_Iyy_Ixy()
    # sample a few locs and ensure no NaNs
    locs = [0, beam.h/2, beam.h + beam.b/2, 2*beam.h + beam.b - 1e-6]
    for loc in locs:
        q = beam.get_shear(shear1=2000, shear2=100, loc=loc)
        beam.create_shearflow_list(shear1=200, shear2=10)
        assert not np.isnan(q)
        print(f"q(s={loc:.4f}) = {q:.2f}")

def test_heatmap_plot():
    mat = Material("TestMat", EModulus=100e9, YieldStrength=400e6, max_shear=250e6)
    beam = Square_beam(length=1.0, w=0.15, h=0.15, t=0.005, material=mat)
    # This will pop up a window (or inline) for you to visually inspect
    # beam.plot_shear_heatmap(shear1=500, shear2=30, torsion_shear=0)
    # beam.plot_shearflow(shear1=500, shear2=30, torsion_shear=0)
    beam.get_Ixx_Iyy_Ixy()
    beam.adjust_for_torsion(shear1=500, shear2=30)

def test_size_profile_convergence():
    mat = Material("Alu 6063t66", EModulus=68e9, YieldStrength=95e6, max_shear=150e6)
    beam = Square_beam(length=1.0, w=0.2, h=0.2, t=0.003, material=mat)
    # Try to size for SF=1.2 with some loads; will plot convergence
    final_t = beam.size_profile(
        shear1=4000, shear2=20000,
        tension=0, moment1=1000, moment2=500,
        SF=1.2, initial_thickness=0.05,
        side_subdivisions=80, plot_convergence=True
    )
    print(f"Final thickness = {final_t:.4f} m")


def test_reaction_and_internal_diagrams():
    # --- 1. Define material and beam geometry ---
    mat = Material("Steel", EModulus=210e9, YieldStrength=250e6, max_shear=150e6, Density=7850)
    beam = Square_beam(length=2.0, w=0.1, h=0.1, t=0.005, material=mat)

    # --- 2. Define applied loads: [Vx, Vy, N, x_pos] ---
    loads = [
        [0.0, -5000.0, 0.0, 0.5],   # downward point-load at x=0.5 m
        [0.0, 30000.0, 0.0, 1.2],   # downward point-load at x=1.2 m
    ]

    # --- 3. Plot using new plot_NVM method ---
    beam.plot_NVM(loads, use_self_weight=np.array([0,9.81,0]), length_subdivisions=500)


def test_thickness_simple_yield_size():
    # --- 1. Material and beam setup ---
    mat = Material("Aluminum", EModulus=70e9, YieldStrength=300e6, max_shear=200e6, Density=2700)
    beam = Square_beam(length=1.0, w=0.05, h=0.05, t=0.002, material=mat)

    # --- 2. Define point-load and compute internals ---
    loads = [[0.0, -1000.0, 0.0, 0.5]]
    beam.compute_internal_loads(loads, use_self_weight=np.array([0,0,0]), length_subdivisions=20)

    # --- 3. Run thickness sizing ---
    beam.thickness_simple_yield_size(loads,
                                     SF=1.5,
                                     initial_thickness=0.002,
                                     use_self_weight=np.array([0,0,0]),
                                     length_subdivisions=20)
    print(f"Computed minimum thickness: {beam.best_t:.4f} m")


def test_size_profile_convergence():
    mat = Material("Alu6063-T6", EModulus=68e9, YieldStrength=95e6, max_shear=150e6, Density=2700)
    beam = Square_beam(length=1.0, w=0.2, h=0.2, t=0.003, material=mat)
    final_t = beam.size_profile(
        shear1=4000, shear2=20000,
        tension=0, moment1=1000, moment2=500,
        SF=1.2, initial_thickness=0.05,
        side_subdivisions=80, plot_convergence=True
    )
    print(f"Final thickness = {final_t:.4f} m")

if __name__ == "__main__":
    test_size_profile_convergence()
    # test_reaction_and_internal_diagrams()
    # test_thickness_simple_yield_size()