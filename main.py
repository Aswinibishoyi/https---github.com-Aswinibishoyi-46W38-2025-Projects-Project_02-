import os

# Change the current working directory to the project folder
os.chdir('D:/DTU 46W38 Scientific Programming in Wind Energy/turbie_project')

# Verify that the directory has been changed
print(f"Current working directory: {os.getcwd()}")

import turbie
from turbie import import_parameters
from turbie import *
import os
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
# import turbie_module  # Assuming turbie_module.py is in the same directory

# --- Project Setup ---
# Define folder paths
project_dir = os.getcwd()
turbie_inputs_dir = os.path.join(project_dir, 'turbie_inputs')
wind_files_dir = os.path.join(project_dir, 'wind_files')
results_dir = os.path.join(project_dir, 'results')

# Create results directory if it doesn't exist
if not os.path.exists(results_dir):
    os.makedirs(results_dir)

# Define parameter files
parameters_file = os.path.join(turbie_inputs_dir, 'turbie_parameters.txt')
ct_lookup_file = os.path.join(turbie_inputs_dir, 'C_T.txt')

# --- Simulation and Analysis ---
def run_simulation():
    """Main function to run the simulation and analyze results."""
    # Define wind speeds and turbulence intensities to simulate
    wind_speeds = [4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25] # Example wind speeds
    ti_categories = [0.05,0.1,0.15]          # Example TIs

    # Load turbine parameters once
    turbine_params = import_parameters(parameters_file)
    if turbine_params is None:
        return

    M = turbine_params['M']
    C = turbine_params['C']
    K = turbine_params['K']
    rho = turbine_params['rho']
    A = turbine_params['A']
    M_inv = np.linalg.inv(M)

    all_results = {}

    for ti in ti_categories:
        all_results[ti] = {'means': [], 'stds': []}
        ti_results_file = os.path.join(results_dir, f'ti_{ti}_stats.txt')
        with open(ti_results_file, 'w') as f:
            f.write("Wind Speed (m/s),Blade Mean (m),Blade Std Dev (m),Tower Mean (m),Tower Std Dev (m)\n")

        for ws in wind_speeds:
            print(f"Simulating wind speed: {ws} m/s, TI: {ti}")
            wind_file = os.path.join(wind_files_dir, f'wind_{ws}_ms_TI_{ti}.txt')
            wind_data_raw = np.loadtxt(wind_file, skiprows=1)
            time_series = wind_data_raw[:, 0]
            wind_speed_series = wind_data_raw[:, 1]

            # Create a wind speed interpolation function for the ODE solver
            wind_speed_interp = interp1d(time_series, wind_speed_series, kind='linear', fill_value='extrapolate')

            # Get the mean wind speed and corresponding thrust coefficient
            mean_wind_speed = np.mean(wind_speed_series)
            
            Ct = get_thrust_coefficient(mean_wind_speed, ct_lookup_file)
            if Ct is None:
                return

            # Initial conditions: y = [x1, x2, x1_dot, x2_dot]
            y0 = np.array([0.0, 0.0, 0.0, 0.0])

            # Time span for the simulation
            t_span = (time_series[0], time_series[-1])

            # Parameters dictionary for the ODE solver
            ode_params = {
                'wind_speed_interp': wind_speed_interp,
                'Ct': Ct,
                'rho': rho,
                'A': A
            }

            # Solve the ordinary differential equation
            sol = solve_ivp(
                lambda t, y: turbie_dynamics(t, y, M_inv, C, K, aerodynamic_forcing_vector, ode_params),
                t_span,
                y0,
                t_eval=time_series
            )

            # Check if the solution was successful
            if sol.success:
                time = sol.t
                blade_displacement = sol.y[0, :]
                tower_displacement = sol.y[1, :]

                # Save the time-series results to a file
                output_file = os.path.join(results_dir, f'results_ws_{ws}_ti_{ti}.txt')
                np.savetxt(
                    output_file,
                    np.vstack((time, wind_speed_series, blade_displacement, tower_displacement)).T,
                    header="Time (s), Wind Speed (m/s), Blade Displacement (m), Tower Displacement (m)",
                    fmt="%.6f",
                    delimiter="\t"
                )

                # Calculate means and standard deviations
                blade_mean = np.mean(blade_displacement)
                blade_std = np.std(blade_displacement)
                tower_mean = np.mean(tower_displacement)
                tower_std = np.std(tower_displacement)

                all_results[ti]['means'].append((ws, blade_mean, tower_mean))
                all_results[ti]['stds'].append((ws, blade_std, tower_std))

                # Append stats to the TI results file
                with open(ti_results_file, 'a') as f:
                    f.write(f"{ws},{blade_mean},{blade_std},{tower_mean},{tower_std}\n")

                # Plot the time-marching variation for one specific case
                if ws == 10 and ti == 0.1:
                    plot_time_series(time, wind_speed_series, blade_displacement, tower_displacement, ws, ti)
            else:
                print(f"Solver failed for ws={ws}, ti={ti}")

    # Plot means and standard deviations for all wind speeds and TIs
    plot_summary_results(all_results)

    # Discussion
    discuss_results(all_results)

def plot_time_series(time, wind_speed, blade_disp, tower_disp, ws, ti):
    """Plots the time series of wind speed and displacements for a single case."""
    fig, ax1 = plt.subplots(figsize=(12, 6))

    color = 'tab:blue'
    ax1.set_xlabel('Time (s)')
    ax1.set_ylabel('Displacement (m)', color=color)
    ax1.plot(time, blade_disp, label='Blade Displacement', color=color)
    ax1.plot(time, tower_disp, label='Tower Displacement', color='tab:green', linestyle='--')
    ax1.tick_params(axis='y', labelcolor=color)
    ax1.grid(True)

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    color = 'tab:red'
    ax2.set_ylabel('Wind Speed (m/s)', color=color)  # we already handled the x-label with ax1
    ax2.plot(time, wind_speed, label='Wind Speed', color=color, linestyle=':')
    ax2.tick_params(axis='y', labelcolor=color)

    # Add legends
    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines + lines2, labels + labels2, loc='upper left')

    plt.title(f'Time-Marching Variation for Wind Speed: {ws} m/s, TI: {ti}')
    plt.tight_layout()
    plt.savefig(os.path.join(results_dir, f'time_series_ws_{ws}_ti_{ti}.png'))
    plt.show()

def plot_summary_results(all_results):
    """Plots the mean and standard deviation of displacements vs. wind speed."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

    # Plot Means
    ax1.set_title('Mean Displacements vs. Wind Speed')
    ax1.set_xlabel('Wind Speed (m/s)')
    ax1.set_ylabel('Mean Displacement (m)')
    for ti, data in all_results.items():
        ws = [d[0] for d in data['means']]
        blade_means = [d[1] for d in data['means']]
        tower_means = [d[2] for d in data['means']]
        ax1.plot(ws, blade_means, 'o-', label=f'Blade (TI={ti})')
        ax1.plot(ws, tower_means, 'x-', label=f'Tower (TI={ti})')
    ax1.grid(True)
    ax1.legend()

    # Plot Standard Deviations
    ax2.set_title('Standard Deviations of Displacements vs. Wind Speed')
    ax2.set_xlabel('Wind Speed (m/s)')
    ax2.set_ylabel('Standard Deviation (m)')
    for ti, data in all_results.items():
        ws = [d[0] for d in data['stds']]
        blade_stds = [d[1] for d in data['stds']]
        tower_stds = [d[2] for d in data['stds']]
        ax2.plot(ws, blade_stds, 'o-', label=f'Blade (TI={ti})')
        ax2.plot(ws, tower_stds, 'x-', label=f'Tower (TI={ti})')
    ax2.grid(True)
    ax2.legend()

    plt.tight_layout()
    plt.savefig(os.path.join(results_dir, 'summary_plots.png'))
    plt.show()

def discuss_results(all_results):
    """
    Prints a discussion of the results based on the plots.
    This discussion would be a starting point and should be completed with
    observations from the generated plots.
    """
    print("\n\nDiscussion and Explanation of Results:")
    print("--------------------------------------")

    # Discussion on Means
    print("1. How means of displacements change with wind speed and TI:")
    print("   - As wind speed increases, the mean displacements for both the blade and tower are expected to increase.")
    print("   - This is because the aerodynamic force on the turbine blades, which is proportional to the square of the effective wind speed, increases with higher average wind speed. This results in a larger static deflection (mean displacement) of the blades and tower.")
    print("   - The mean displacements are not expected to change significantly with the turbulence intensity (TI), as TI primarily affects the fluctuations around the mean wind speed, not the mean itself. The plots should show the mean displacements for different TIs following a similar trend with increasing wind speed.")

    # Discussion on Standard Deviations
    print("\n2. How standard deviations of displacements change with wind speed and TI:")
    print("   - The standard deviation of the displacements, which represents the amplitude of fluctuations, is expected to increase with both increasing wind speed and increasing turbulence intensity.")
    print("   - With a higher average wind speed, the magnitude of the aerodynamic force is larger, so any fluctuations in wind speed will result in larger fluctuations in the applied force, leading to greater displacement variations.")
    print("   - A higher turbulence intensity (TI) means there are greater fluctuations in the wind speed around the mean. This directly translates to a more variable aerodynamic forcing, causing larger oscillations and a higher standard deviation in the blade and tower displacements.")
    print("   - The plots should illustrate a clear separation between the TI categories, with higher TI cases showing a larger standard deviation for the same wind speed.")

if __name__ == '__main__':
    # You would need to provide the actual files in the corresponding folders
    # before running this script.
    run_simulation()