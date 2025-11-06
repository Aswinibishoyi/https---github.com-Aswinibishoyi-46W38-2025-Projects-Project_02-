import os

# Change the current working directory to the project folder
os.chdir('D:/DTU 46W38 Scientific Programming in Wind Energy/turbie_project')

# Verify that the directory has been changed
print(f"Current working directory: {os.getcwd()}")



import numpy as np
from scipy.interpolate import interp1d

def import_wind_file(filepath):
    """
    Imports wind speed data from a text file.

    Args:
        filepath (str): The path to the wind data file.

    Returns:
        np.ndarray: A NumPy array of wind speeds over time.
    """
    try:
        wind_data = np.loadtxt(filepath, skiprows=1)
        return wind_data[:, 1]
    except FileNotFoundError:
        print(f"Error: Wind file not found at {filepath}")
        return None

def get_thrust_coefficient(mean_wind_speed, lookup_table_filepath):
    """
    Determines the thrust coefficient (CT) via interpolation.

    Args:
        mean_wind_speed (float): The average wind speed for the simulation.
        lookup_table_filepath (str): Path to the CT lookup table file.

    Returns:
        float: The interpolated thrust coefficient.
    """
    try:
        ct_data = np.loadtxt(lookup_table_filepath, skiprows=1)
        wind_speeds = ct_data[:, 0]
        ct_values = ct_data[:, 1]

        # Create an interpolation function
        ct_interpolator = interp1d(wind_speeds, ct_values, kind='linear', fill_value='extrapolate')

        # Return the interpolated CT value
        return ct_interpolator(mean_wind_speed).item()
    except FileNotFoundError:
        print(f"Error: C_T.txt lookup table not found at {lookup_table_filepath}")
        return None

def import_parameters(filepath):
    """
    Imports turbine parameters from a text file.

    Args:
        filepath (str): Path to the parameters file.

    Returns:
        dict: A dictionary containing the system matrices and other parameters.
    """
    params = {}
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
            for line in lines:
                parts = line.split('=')
                if len(parts) == 2:
                    key = parts[0].strip()
                    value = parts[1].strip()
                    try:
                        # Try to convert to float
                        params[key] = float(value)
                    except ValueError:
                        # Handle matrices
                        if '[' in value:
                            # Parse a 2x2 matrix
                            rows = value.replace('[', '').replace(']', '').split(';')
                            matrix = []
                            for row in rows:
                                # Split by whitespace and handle potential commas
                                row_values = [float(v.strip().replace(',', '')) for v in row.split()]
                                matrix.append(row_values)
                            params[key] = np.array(matrix)
                        else:
                            params[key] = value

        # Assume M, C, K are 2x2 matrices, and the rest are scalars
        M = params.get('M')
        C = params.get('C')
        K = params.get('K')

        if M is not None and C is not None and K is not None:
            return {
                'M': M,
                'C': C,
                'K': K,
                'rho': params.get('rho', 1.225),  # Default air density
                'A': params.get('A')
            }
        else:
            print("Error: Missing M, C, or K matrix in turbie_parameters.txt")
            return None
    except FileNotFoundError:
        print(f"Error: Parameters file not found at {filepath}")
        return None

def aerodynamic_forcing_vector(t, y, wind_speed_interp, Ct, rho, A):
    """
    Calculates the aerodynamic forcing vector F(t).

    Args:
        t (float): Current time.
        y (np.ndarray): Current state vector [x1, x2, x1_dot, x2_dot].
        wind_speed_interp (callable): Interpolation function for wind speed.
        Ct (float): Thrust coefficient.
        rho (float): Air density.
        A (float): Rotor area.

    Returns:
        np.ndarray: The forcing vector [faero, 0].
    """
    # Extract blade displacement and velocity from state vector
    x1_dot = y[2]

    # Get wind speed at time t
    u_t = wind_speed_interp(t).item()

    # Calculate aerodynamic forcing
    faero = 0.5 * rho * Ct * A * (u_t - x1_dot) * abs(u_t - x1_dot)

    return np.array([faero, 0.0])

def turbie_dynamics(t, y, M_inv, C, K, F_func, params):
    """
    Defines the time derivative of the state vector for the ODE solver.

    Args:
        t (float): Current time.
        y (np.ndarray): Current state vector [x1, x2, x1_dot, x2_dot].
        M_inv (np.ndarray): Inverse of the mass matrix.
        C (np.ndarray): Damping matrix.
        K (np.ndarray): Stiffness matrix.
        F_func (callable): Function to calculate the forcing vector.
        params (dict): Dictionary of parameters for the forcing function.

    Returns:
        np.ndarray: The time derivative of the state vector.
    """
    x = y[:2]
    x_dot = y[2:]

    # Calculate the forcing vector F(t)
    F = F_func(t, y, params['wind_speed_interp'], params['Ct'], params['rho'], params['A'])

    # Calculate x_ddot using the equation of motion
    # M*x_ddot + C*x_dot + K*x = F
    x_ddot = M_inv @ (F - C @ x_dot - K @ x)

    # Return the full derivative vector
    return np.concatenate((x_dot, x_ddot))