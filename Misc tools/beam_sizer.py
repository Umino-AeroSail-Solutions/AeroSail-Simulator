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
        minSF = self.Yield / np.max(tensions)
        return minSF

    def get_shear(self, shear1, shear2, loc):
        Vx = shear1
        Vy = shear2
        self.b = self.w - 2 * self.t
        self.d = self.h - 2 * self.t

        qb_1 = lambda s: -Vy / self.Ixx * self.t * (-self.h / 2 * s + s ** 2 / 2) - Vx / self.Iyy * self.t * (
                    (self.w + self.d) / 4) * s
        qb_2 = lambda s: -Vy / self.Ixx * self.t * ((self.h + self.b) / 4) * s - Vx / self.Iyy * self.t * (
                    self.d / 2 * s - s ** 2 / 2)
        qb_3 = lambda s: -Vy / self.Ixx * self.t * (self.h / 2 * s - s ** 2 / 2) - Vx / self.Iyy * self.t * -(
                    (self.w + self.d) / 4) * s
        qb_4 = lambda s: -Vy / self.Ixx * self.t * (-(self.h + self.b) / 4) * s - Vx / self.Iyy * self.t * (
                    -self.d / 2 * s + s ** 2 / 2)

        # # Corrected qs0
        # qs0 = torsion_shear / (2 * ((self.w + self.d) / 2) * ((self.h + self.b) / 2))
        # Wroooooong INTEGRATION CONSTANTS

        # Corrected perimeter wrap
        perimeter = 2 * self.h + 2 * self.d
        loc = loc % perimeter

        if 0 <= loc <= self.h:
            return qb_1(loc)
        elif self.h < loc <= self.h + self.d:
            loc -= self.h
            return qb_2(loc) + qb_1(self.h)
        elif self.h + self.d < loc <= 2 * self.h + self.d:
            loc -= self.h + self.d
            return qb_3(loc) + qb_2(self.d) + qb_1(self.h)
        elif 2 * self.h + self.d < loc <= perimeter:
            loc -= 2 * self.h + self.d
            return qb_4(loc) + qb_3(self.h) + qb_2(self.d) + qb_1(self.h)

    def create_shearflow_list(self, shear1, shear2, side_subdivisions=80):
        # Perimeter for shear sampling
        P = 2 * self.h + 2 * self.d

        # Compute initial shear safety factor
        shearflows = []
        for loc in np.linspace(0, P, side_subdivisions, endpoint=False):
            shearflows.append(self.get_shear(shear1, shear2, loc))
        self.shearflows = np.array([shearflows, np.linspace(0, P, side_subdivisions, endpoint=False)]).T
        # print(self.shearflows)
        return self.shearflows

    def adjust_for_torsion(self, shear1, shear2, torsion=0, side_subdivisions=80):
        self.b = self.w - self.t
        self.d = self.h - self.t
        P = 2 * self.b + 2 * self.d
        N = side_subdivisions
        ds = P / N

        self.create_shearflow_list(shear1, shear2, side_subdivisions)

        total_inner_qs0 = 0

        for q, loc in self.shearflows:
            if 0 <= loc < self.h:
                lever_arm = self.d / 2  # left wall (vertical)
            elif self.h <= loc < self.h + self.b:
                lever_arm = self.b / 2  # bottom wall (horizontal)
            elif self.h + self.b <= loc < 2 * self.h + self.b:
                lever_arm = self.d / 2  # right wall
            elif 2 * self.h + self.b <= loc <= P:
                lever_arm = self.b / 2  # top wall
            else:
                continue  # shouldn't happen, but safe

            total_inner_qs0 += q * lever_arm * ds

        print("Total torsional moment (from q):", total_inner_qs0)

        # If you want to correct for torsion (optional):
        A = self.h*self.w
        q0 = total_inner_qs0 / (2 * A)
        self.shearflows[:, 0] -= q0

        # Plot
        plt.plot(self.shearflows[:, 1], self.shearflows[:, 0])
        plt.xlabel("Location along perimeter, s")
        plt.ylabel("Shear flow q(s)")
        plt.title("Shear‐flow distribution around beam perimeter")
        plt.tight_layout()
        plt.show()


    # def get_shears(self, shear1, shear2, torsion_shear, loc):
    #     # Great shit anhong what the fuck is self.d and self.b it was not defined :(
    #     Vx = shear1
    #     Vy = shear2
    #     self.b = self.w - 2 * self.t
    #     self.d = self.h - 2 * self.t
    #     #     w
    #     # ||=======|      x
    #     # ||   d   |    <---+ (not the shear center, just for coordinate directions)
    #     # ||     b | h      |
    #     # ||       |        | y
    #     # ||-------|        v
    #     # side 1 || s: 0 --> h
    #     qb_1 = lambda s: -Vy/self.Ixx * self.t * (-self.h/2*s + (s**2) /2) - Vx/self.Iyy * self.t * ((self.w+self.d)/4)*s
    #     # side 2 -  s: 0 --> d
    #     qb_2 = lambda s: -Vy/self.Ixx * self.t * ((self.h+self.b)/4)*s - Vx/self.Iyy * self.t * (self.d/2*s - (s**2)/2)
    #     # side 3 |  s: 0 --> h
    #     qb_3 = lambda s: -Vy/self.Ixx * self.t * (self.h/2*s - (s**2) /2) - Vx/self.Iyy * self.t * -((self.w+self.d)/4)*s
    #     # side 4 _  s: 0 --> d
    #     qb_4 = lambda s: -Vy/self.Ixx * self.t * (-(self.h+self.b)/4)*s - Vx/self.Iyy * self.t * (-self.d/2*s + (s**2)/2)
    #     # qs0
    #     qs0 = torsion_shear/2/(((self.w+self.d)/2)*((self.h+self.b)/2))
    #
    #     # what's below could have been a nested-if but i made it like this for readability :3 (oh my god wtf is that face ew furry ew what)
    #     loc = loc % 2*self.h + 2*self.d
    #     if loc >= 0 and loc <= self.h:
    #         return qb_1(loc) + qs0
    #     elif loc > self.h and loc <= (self.h + self.d):
    #         loc = loc - self.h
    #         return qb_2(loc) + qb_1(self.h) + qs0
    #     elif loc > (self.h + self.d) and loc <= (2*self.h + self.d):
    #         loc = loc - self.h - self.d
    #         return qb_3(loc) + qb_2(self.d) + qb_1(self.h) + qs0
    #     elif loc > (2*self.h + self.d) and loc <= (2*self.h + 2*self.d):
    #         loc = loc - 2*self.h - self.d
    #         return qb_4(loc) + qb_3(self.h) + qb_2(self.d) + qb_1(self.h) + qs0



    def size_profile(self, shear1, shear2, tension, moment1, moment2, SF,
                     initial_thickness, side_subdivisions=50, plot_convergence=False):
        """
        Iteratively adjust thickness to meet a target safety factor based on shear and tension.
        If plot_convergence is True, displays how thickness evolves each iteration.
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

        # Store thickness history if requested
        if plot_convergence:
            last_ts = [self.t]

        # Compute initial tension safety factor
        self.tensionSF = self.gettensionSF(moment1, moment2, tension)

        # Perimeter for shear sampling
        P = 2*self.h + 2*self.d

        # Compute initial shear safety factor
        shearflows = []
        for loc in np.linspace(0, P, side_subdivisions, endpoint=False):
            shearflows.append(self.get_shears(shear1, shear2, 0, loc))
        max_shear = np.max(shearflows)
        self.shearSF = max_shear / self.max_shear_yield

        # Determine actual safety factor and design ratio
        self.actualSF = min(self.shearSF, self.tensionSF)
        design_sf = self.actualSF / SF

        # Iteratively adjust thickness until within [1.0, 1.1] of target
        iteration = 0
        while design_sf < 1 or design_sf > 1.1:
            iteration += 1
            # Scale thickness
            self.t = self.t/design_sf
            # Record
            if plot_convergence:
                last_ts.append(self.t)
            # Update inertia and factors
            self.get_Ixx_Iyy_Ixy()
            self.tensionSF = self.gettensionSF(moment1, moment2, tension)

            shearflows = []
            for loc in np.linspace(0, P, side_subdivisions, endpoint=False):
                shearflows.append(self.get_shears(shear1, shear2, 0, loc))
            max_shear = np.max(shearflows)
            self.shearSF = max_shear / self.max_shear_yield

            self.actualSF = min(self.shearSF, self.tensionSF)
            print(design_sf)
            # Safety break
            if iteration > 1000:
                raise RuntimeError("Thickness iteration did not converge after 1000 steps")

        # Plot convergence if requested
        if plot_convergence:
            plt.plot(range(len(last_ts)), last_ts)
            plt.xlabel('Iteration')
            plt.ylabel('Thickness t')
            plt.title('Convergence of thickness sizing')
            plt.tight_layout()
            plt.show()

        return self.t
    # def thickness_simple_yield_size(self, loads, SF, initial_thickness, use_self_weight=np.array([0,0,0]),
    #                                 length_subdivisions=30):
    #     '''loads are defined as [shear1, shear2, normal tension, position across length], use self weight is an
    #     acceleration vector [x, y, z] where x is shear 1, y is shear 2 and z is normal tension'''


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
        P = 2 * self.h + 2 * self.d
        locs = np.linspace(0, P, subdivisions, endpoint=False)

        xs, ys, qs = [], [], []
        for loc in locs:
            # Wrap loc into [0, P)
            s = loc % P
            # Map s to (x, y) on mid‐line
            if s <= self.h:
                x = -self.w / 2
                y = self.h / 2 - s
            elif s <= self.h + self.d:
                x = -self.w / 2 + (s - self.h)
                y = -self.h / 2
            elif s <= 2 * self.h + self.d:
                x = self.w / 2
                y = -self.h / 2 + (s - self.h - self.d)
            else:
                x = self.w / 2 - (s - 2 * self.h - self.d)
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
        P = 2 * self.h + 2 * self.d
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
    locs = [0, beam.h/2, beam.h + beam.d/2, 2*beam.h + beam.d - 1e-6]
    for loc in locs:
        q = beam.get_shear(shear1=200, shear2=10, loc=loc)
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
    mat = Material("Titanium", EModulus=110e9, YieldStrength=900e6, max_shear=600e6)
    beam = Square_beam(length=1.0, w=0.5, h=0.5, t=0.003, material=mat)
    # Try to size for SF=1.2 with some loads; will plot convergence
    final_t = beam.size_profile(
        shear1=40, shear2=200,
        tension=0, moment1=800, moment2=200,
        SF=1.2, initial_thickness=0.01,
        side_subdivisions=80, plot_convergence=True
    )
    print(f"Final thickness = {final_t:.4f} m")

if __name__ == "__main__":
    test_inertia()
    test_tension_SF()
    test_shearflows()
    # Plot-based tests (uncomment if running interactively)
    test_heatmap_plot()
    test_size_profile_convergence()
    print("All numeric tests passed.")




