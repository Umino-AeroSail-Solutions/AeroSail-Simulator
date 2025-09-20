import numpy as np
import matplotlib.pyplot as plt

class Rivet():
    def __init__(self, Shear_Yield_Stress, diameter, Tension_Yield, cost, name="Unnamed Chinese Kid"):
        self.Shear_Yield = Shear_Yield_Stress
        self.diameter = diameter
        self.Tension_Yield = Tension_Yield
        self.cost = cost
        self.name = name

class Panel():
    def __init__(self, thickness, shear_yield, tension_yield, E_mod, max_bearing_stress):
        self.thickness = thickness
        self.min_rivet_d = 3*self.thickness
        self.max_bearing_stress = max_bearing_stress
        self.E_mod = E_mod
        self.shear_yield = shear_yield
        self.tension_yield = tension_yield

    def Get_Min_Diameter_Bearing(self, force):
        min_diameter_bearing = force / (self.max_bearing_stress*self.thickness) # Bearing min diameter comp https://roymech.org/Useful_Tables/Rivets.html
        return min(min_diameter_bearing, self.min_rivet_d)
    def Get_Max_Total_Hole_D(self, force, width):
        max_D = (force/(self.tension_yield*self.thickness)) + width
        return max_D
    def Get_Max_allowed_force(self, diameter, n, width):
        force_bearing = n*diameter*(self.max_bearing_stress*self.thickness)
        force_tension = (width-n*diameter)*(self.tension_yield*self.thickness)
        self.maxallowedforce = min(force_bearing, force_tension)
        return min(force_bearing, force_tension)
    def Get_Max_tearing_force(self, spacing, rivets_per_row):
        force_tearing = self.shear_yield*2*spacing*self.thickness*rivets_per_row
        return force_tearing
class Connection():
    def __init__(self, stress, Panel1, Panel2, width, SF=4.5): # Normal safety factor is huge due to 3x stress concentrations in holes
        self.SF = SF
        self.stress = stress*self.SF
        self.Panel1 = Panel1
        self.Panel2 = Panel2
        self.width = width
        self.shearflow= self.stress*self.width

    def get_required_spacing(self, Rivet, rivets_per_row=None, initialspacing=None, plot=True):
        print("COMPUTING RIVET: \n______________________________________\n", Rivet.name, "\n______________________________________\n")
        if initialspacing is None:
            initialspacing = Rivet.diameter * 3

        if rivets_per_row is None:
            minimum_spacing = 3*Rivet.diameter
            rivets_per_row = (self.width//minimum_spacing) -1

        if rivets_per_row <= 0:
            rivets_per_row = 1
        self.rivet = Rivet
        self.rivet_max_force = self.rivet.Shear_Yield
        self.connection_d = self.rivet.diameter

        spacing = initialspacing

        # Iteration 0
        iteration_n = 0
        force = self.stress * spacing * self.width
        self.rivets_per_row = rivets_per_row
        self.max_allowed_panel_force = min(
            self.Panel1.Get_Max_allowed_force(self.connection_d, self.rivets_per_row, self.width), self.Panel2.Get_Max_allowed_force(self.connection_d, self.rivets_per_row, self.width)
        )
        self.max_tearing_force = min(
            self.Panel1.Get_Max_tearing_force(spacing, self.rivets_per_row),
            self.Panel2.Get_Max_tearing_force(spacing, self.rivets_per_row),
        )

        rivetSF = self.rivet_max_force / force
        panelSF = self.max_allowed_panel_force / force
        tearingSF = self.max_tearing_force / force
        minSF = min(rivetSF, panelSF, tearingSF)
        # while tearingSF < 1 and rivets_per_row > 1:
        #     rivets_per_row -= 1
        #     self.max_tearing_force = min(
        #         self.Panel1.Get_Max_tearing_force(spacing, self.rivets_per_row),
        #         self.Panel2.Get_Max_tearing_force(spacing, self.rivets_per_row),
        #     )
        #     force = self.stress * spacing*self.width
        #     tearingSF = self.max_tearing_force / force

        # store convergence history
        iterations, rivet_history, panel_history, tearing_history, spacing_history = [], [], [], [], []

        while (((minSF - 1) > 0.1) or (minSF - 1 < 0) or (iteration_n<2)) and iteration_n < 100:

            force = self.stress * spacing * self.width
            self.max_tearing_force = min(
                self.Panel1.Get_Max_tearing_force(spacing, self.rivets_per_row),
                self.Panel2.Get_Max_tearing_force(spacing, self.rivets_per_row),
            )
            rivetSF = self.rivet_max_force / force
            panelSF = self.max_allowed_panel_force / force
            tearingSF = self.max_tearing_force / force
            minSF = min(rivetSF, panelSF)


            # record values
            iterations.append(iteration_n)
            rivet_history.append(rivetSF)
            panel_history.append(panelSF)
            tearing_history.append(tearingSF)
            spacing_history.append(spacing)

            spacing = spacing * minSF

            iteration_n += 1
        if plot:
            # Plot convergence
            plt.figure(figsize=(8, 5))
            plt.plot(iterations, rivet_history, label="Rivet SF")
            plt.plot(iterations, panel_history, label="Panel SF")
            plt.plot(iterations, tearing_history, label="Tearing SF")
            plt.axhline(1.0, color="red", linestyle="--", label="Target SF = 1")
            plt.xlabel("Iteration")
            plt.ylabel("Safety Factor")
            plt.title("Convergence of Safety Factors: " + Rivet.name)
            plt.legend()
            # Set integer ticks on x
            plt.xticks(np.arange(0, max(iterations) + 1, 1))

            # Vertical grid lines only at x-ticks (iterations)
            plt.grid(axis="x", which="major")

            # Horizontal grid on y as usual
            plt.grid(axis="y")
            plt.show()

            # --- Plot 2: Spacing history ---
            plt.figure(figsize=(8, 5))
            plt.plot(iterations, np.array(spacing_history) * 1000, marker="o")
            plt.axhline(3*self.rivet.diameter*1000, color="red", linestyle="--", label="Minumum recommended spacing")
            plt.legend()
            plt.xlabel("Iteration")
            plt.ylabel("Spacing [mm]")
            plt.title("Convergence of Rivet Spacing: " + Rivet.name)
            # Set integer ticks on x
            plt.xticks(np.arange(min(iterations), max(iterations) + 1, 1))

            # Vertical grid lines only at x-ticks (iterations)
            plt.grid(axis="x", which="major")

            # Horizontal grid on y as usual
            plt.grid(axis="y")
            plt.show()

        tearing = False

        if tearingSF < 1:
            print("WARNING! Tearing sf to small: ", tearingSF)
            tearing = True

        if spacing < 3*self.rivet.diameter:
            print("WARNING! Spacing to small: ", spacing*1000, "mm, it should be more than: ", 3*self.rivet.diameter*1000, "mm")

        return spacing, rivets_per_row, min(tearingSF, panelSF, rivetSF)*self.SF, tearing

    def find_The_Chosen_One(self, rivet_list, cost_multiplier=1, number_multiplier=0):
        best_optimization = 10000000000000000000000
        best_rivet = rivet_list[0]
        for rivet in rivet_list:
            spacing, rivets_per_row, minSF, tearing = self.get_required_spacing(rivet, plot=True)
            rivets_per_meter = rivets_per_row /spacing

            optimization = (cost_multiplier * rivets_per_meter * rivet.cost) + number_multiplier*rivets_per_meter
            print(f"Cost per meter = {optimization:.1f}")
            print(f"Rivets per row = {rivets_per_row:.1f}")
            print(f"Min SF = {minSF:.1f}")
            print(f"Converged rivet spacing ≈ {spacing * 1000:.1f} mm")
            if optimization < best_optimization and not tearing:
                best_optimization = optimization
                best_rivet = rivet

        print("\n\nBest rivet is: ", best_rivet.name, "---------------------------------")

        spacing, rivets_per_row, minSF, tearing = self.get_required_spacing(best_rivet, plot=True)
        print(f"Cost per meter = {best_optimization:.1f}")
        print(f"Rivets per row = {rivets_per_row:.1f}")
        print(f"Min SF = {minSF:.1f}")
        print(f"Converged rivet spacing ≈ {spacing * 1000:.1f} mm")
        return best_rivet


# EXAMPLE CODE ----------------------------------------------------------------------------
# Example rivet – 4 mm steel rivet
# steel_rivet = Rivet(
#     Shear_Yield_Stress = 3300.0,  # N (shear strength in single shear ≈ 3.3 kN)
#     diameter = 0.004,             # m
#     Tension_Yield = 4500.0,       # N (approx. tension capacity)
#     cost = 0.10
# )
#
# # Panels – mild steel
# thickness1 = 0.001  # 5 mm
# thickness2 = 0.001  # 3 mm
# shear_yield = 250e6       # Pa
# tension_yield = 250e6     # Pa
# E_mod = 200e9             # Pa
# max_bearing_stress = 300e6  # Pa
#
# panel5 = Panel(thickness1, shear_yield, tension_yield, E_mod, max_bearing_stress)
# panel3 = Panel(thickness2, shear_yield, tension_yield, E_mod, max_bearing_stress)
#
# # Connection setup
# # Applied stress = 100 MPa, width = 0.1 m
# conn = Connection(
#     stress = 10e5,
#     Panel1 = panel5,
#     Panel2 = panel3,
#     width = 0.1,
#     SF = 2
# )
#
# # Compute required spacing with 3 rivets per row
# spacing, rivets_per_row, minSF = conn.get_required_spacing(steel_rivet)
#
# print(f"Rivets per row = {rivets_per_row:.1f}")
# print(f"Min SF = {minSF:.1f}")
# print(f"Converged rivet spacing ≈ {spacing*1000:.1f} mm")