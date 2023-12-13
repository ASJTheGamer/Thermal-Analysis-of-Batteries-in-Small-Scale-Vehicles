from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import math
import matplotlib.pyplot as plt
import pandas as pd
import mplcursors

pdf = PdfPages('PDF_1.pdf')

vehicle_weight = 120  # kg
frontal_area = 0.6  # m^2
rolling_resistance_coefficient = 0.02
air_density = 1.2  # kg/m^3
air_drag_coefficient = 0.5
wheel_radius = 0.3  # m
# Define driving conditions

x = []
y = []

linear_acceleration = 3  # m/s^2
angular_acceleration = 1  # rad/s^2
# Define load range
load_range = np.arange(0, 101, 50)  # kg
velocity_range = np.arange(0, 81, 10) / 3.6  # m/s
grade_range = np.arange(0, 11, 1)  # %
# Define constants
gravity = 9.81  # m/s^2
# Define the constants for the calculation
air_density = 1.2  # kg/m^3
air_specific_heat = 1005  # J/kg-K
air_dynamic_viscosity = 1.84e-5  # Pa-s
air_thermal_conductivity = 0.026  # W/m-K
hot_plate_material_thermal_conductivity = 200  # W/m-K
hot_plate_thickness = 0.005  # m
ambient_temperature = 30  # °C
# Define the parameters for the cold plate
hot_plate_surface_area = 0.032  # m^2
fin_height = 0.003  # m
fin_thickness = 0.003  # m,
fin_density = 100  # fins per meter
# Define the ranges to search through
air_inflow_rates = [0.0212436]  # m^3/s
air_inlet_temperatures = range(30, 40, 3)  # °C
heating_loads = []
nums = 1;
# Loop through load range
for load in load_range:
    # Loop through velocity range
    for velocity in velocity_range:
        # Loop through grade range
        for grade in grade_range:
            # Calculate rolling resistance force
            rolling_resistance_force = rolling_resistance_coefficient * (vehicle_weight + load) * gravity
            # Calculate aerodynamic drag force
            aerodynamic_drag_force = 0.5 * air_density * frontal_area * air_drag_coefficient * velocity ** 2
            # Calculate grade resistance force
            grade_resistance_force = (vehicle_weight + load) * grade / 100 * gravity
            # Calculate linear acceleration force
            linear_acceleration_force = (vehicle_weight + load) * linear_acceleration
            # Calculate angular acceleration torque
            angular_acceleration_torque = (vehicle_weight + load) * wheel_radius * angular_acceleration
            # Calculate total force required
            total_force = rolling_resistance_force + aerodynamic_drag_force + grade_resistance_force + linear_acceleration_force
            # Calculate total torque required
            total_torque = total_force * wheel_radius
            # Calculate power required
            power = total_torque * angular_acceleration
            # Display results
            print(f"Number of calculation: {nums} ")
            nums += 1;
            print("Load:", load, "kg")
            print("Velocity:", velocity, "m/s")
            print("Grade:", grade, "%")
            print("Rolling resistance force:", rolling_resistance_force, "N")
            print("Aerodynamic drag force:", aerodynamic_drag_force, "N")
            print("Grade resistance force:", grade_resistance_force, "N")
            print("Linear acceleration force:", linear_acceleration_force, "N")
            print("Total force required:", total_force, "N")
            print("Total torque required:", total_torque, "Nm")
            print("Power required:", power, "W")

            TORQUE_TO_POWER_RATIO = 9.5488  # Conversion factor from torque (Nm) to power (W) at 1 rad/s
            # Input parameters
            torque_requirement = total_torque
            power_requirement = power
            power_requirement *= 1000.0
            # Calculation of motor specifications
            rated_torque = torque_requirement / 0.8  # Assuming 80% efficiency
            rated_power = power_requirement / TORQUE_TO_POWER_RATIO
            rated_speed = rated_power / rated_torque
            motor_size = math.ceil(rated_power / 1000.0)  # Assuming 1 kW per motor size
            # Display the motor specifications
            print(f"Rated torque: {rated_torque:.2f} Nm")
            print(f"Rated power: {rated_power:.2f} W")
            print(f"Rated speed: {rated_speed:.2f} rad/s")
            print(f"Motor size: {motor_size} kW")

            motor_power = rated_power  # Motor power in watts
            motor_voltage = 48  # Motor voltage in volts
            range_requirement = 100  # Range requirement in km
            cell_capacity = 20  # Cell capacity in Ah
            cell_voltage = 3.3  # Cell voltage in volts
            # Calculate total energy requirement
            energy_requirement = motor_power * (range_requirement / 40)  # Assuming an average speed of 40 km/h
            print("Total energy requirement:", energy_requirement, "Wh")
            # Calculate total voltage requirement
            total_voltage = motor_voltage
            print("Total voltage requirement:", total_voltage, "V")
            # Calculate number of cells required
            cells_in_series = round(total_voltage / cell_voltage)
            print("Cells in series:", cells_in_series)
            # Calculate capacity requirement
            capacity_requirement = energy_requirement / total_voltage
            print("Capacity requirement:", capacity_requirement, "Ah")
            # Calculate cells in parallel
            cells_in_parallel = round(capacity_requirement / cell_capacity)
            print("Cells in parallel:", cells_in_parallel)
            # Calculate total cells required
            total_cells = cells_in_series * cells_in_parallel
            print("Total cells required:", total_cells)
            # Select a suitable battery module
            battery_capacity = cells_in_parallel * cell_capacity
            print("Battery capacity:", battery_capacity, "Ah")
            # Display the battery module configuration
            print("Battery module configuration:", cells_in_series, "S", cells_in_parallel, "P")


            def heat_generation_rate(I, R, dV_dt, Q, R_temp, T, I_int, R_int, V_nominal, V_oc):
                V = V_nominal - I * R - dV_dt * (V_oc - V_nominal)  # Calculate the voltage across the cell
                R_temp_T0 = R * (1 + R_temp * (T - 25))  # Calculate the resistance at the given temperature
                return I * 2 * R + I * dV_dt * (V_oc - V) - Q + I * 2 * R_temp_T0 + I * I_int * R_int


            # Calculate the current flowing through the cell and the battery

            P = motor_power  # motor power in Watts
            V_cell = 3.3  # voltage of each cell in Volts
            Ah_cell = 20  # capacity of each cell in Amp-hours
            num_cells = cells_in_series  # number of cells connected in series in the battery to provide 48V
            V_battery = V_cell * num_cells  # voltage of the battery in Volts
            Ah_battery = Ah_cell * num_cells  # capacity of the battery in Amp-hours
            I = P / V_battery  # current in Amps
            # Calculate the heat generation rate in each individual cell
            R = 0.05  # internal resistance in Ohms
            dV_dt = -0.001  # variation in open circuit voltage with time in Volts per second
            Q = 2  # heat dissipated in Watts
            R_temp = 0.002  # temperature coefficient of resistance in Ohms per degree Celsius
            T = 35  # temperature in degrees Celsius
            I_int = 0.1  # internal current of the cell in Amps
            R_int = 0.02  # internal resistance of the cell in Ohms
            V_oc = 3.5  # open circuit voltage of the cell in Volts
            V_nominal = 3.3  # nominal voltage of the cell in Volts
            heat_gen_rate_cell = heat_generation_rate(I, R, dV_dt, Q, R_temp, T, I_int, R_int, V_nominal, V_oc)
            heating_loads.append(heat_gen_rate_cell)
            print("Heat generation in Cell: {:.2f} Watts".format(heat_gen_rate_cell))
            print("\n")
            # Loop through each air inflow rate
            for air_inflow_rate in air_inflow_rates:
                # Loop through each air inlet temperature
                for air_inlet_temperature in air_inlet_temperatures:
                    # Loop through each cooling load
                    for heating_load in heating_loads:
                        # Calculate the mass flow rate of air through the cold plate
                        air_mass_flow_rate = air_inflow_rate / 10 * air_density

            # Calculate the Reynolds number for the air flow over the fin
            reynolds_number = air_density * air_mass_flow_rate / 10 * fin_height / air_dynamic_viscosity
            # Calculate the Prandtl number for the air flow over the fin
            prandtl_number = air_specific_heat * air_dynamic_viscosity / air_thermal_conductivity
            # Calculate the Nusselt number for the air flow over the fin
            nusselt_number = 0.664 * reynolds_number ** 0.5 * prandtl_number ** 0.33 * (
                    fin_thickness / fin_height) ** 0.25
            # Calculate the heat transfer coefficient for the fin
            h_fin = nusselt_number * air_thermal_conductivity / fin_height
            # Calculate the overall heat transfer coefficient for the cold plate
            h_overall = h_fin + hot_plate_material_thermal_conductivity / hot_plate_thickness
            # Calculate the required cold plate temperature
            hot_plate_temperature = ambient_temperature + heating_load / (h_overall * hot_plate_surface_area)
            # Calculate the outlet temperature of the air through the cold plate
            air_outlet_temperature = hot_plate_temperature + heating_load / (air_mass_flow_rate * air_specific_heat)
            # Print the results
            print(f"Air inflow rate: {air_inflow_rate / 10} m^3/s")
            print(f"Air inlet temperature: {air_inlet_temperature} °C")
            print(f"heating load: {heating_load} W/m^2")
            x.append(heating_load)
            print(f"hot plate temperature: {hot_plate_temperature:.2f} °C")
            print(f"Air outlet temperature: {air_outlet_temperature:.2f} °C")
            y.append(air_outlet_temperature)

            print("\n")


def show_info(sel):
    ind = sel.target.index
    sel.annotation.set_text(f'{df.loc[ind, "Name"]}:\nHeating load:{df.loc[ind, "x"]}\nAOT:{df.loc[ind, "y"]}')


df = pd.DataFrame({'Name': ["".join('Pointing at') for _ in x],
                   'x': x,
                   'y': y})

plt.scatter(x, y, label="data point", color="green", marker=".", s=10)

# x-axis label
plt.xlabel('heating-load W/m^2')
# frequency label
plt.ylabel('air-outlet-temperature °C')
# plot title
plt.title('Heat load vs Air outlet')
# showing legend
plt.legend()
mplcursors.cursor(hover=True).connect('add', show_info)

plt.savefig('This.pdf')
# function to show the plot
plt.show()
