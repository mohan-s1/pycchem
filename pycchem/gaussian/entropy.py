#/* --------------------------------------------------------------------------------
#   Paolucci Group
#   University of Virginia
#   Mohan Shankar
#
#   entropy.py
#
#   This file takes in data from log files of frequency calculations from Gaussian
#   and applies statistical mechanics to calculate partition functions 
#-------------------------------------------------------------------------------- */
import numpy as np
import re
#-------------------------------------------------------------------------------- */
# HELPER FUNCTIONS TO READ IN DATA
def gauss_freq(file_path:str):
    """
    Reads in all frequences (units of cm^-1) from Gaussian log file and discards anything less than 100 cm^-1

    Args:
        file_path (str): Path to gaussian frequency log file

    Returns:
        _type_: NumPy array of frequencies greater than or equal to 100 cm^-1
    """
    pattern = r'Frequencies\s+--\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)' #regex to find desires lines
    frequencies = []

    with open(file_path, 'r') as file:
        for line in file:
            match = re.search(pattern, line)
            if match:
                freq_values = match.groups()
                frequencies.extend(freq_values)
    
    # turn strings to floats and cutoff at 100 cm-1
    
    no_low_freq = [float(x) for x in frequencies if float(x) >= 100]
    
    # Convert the list to a NumPy array with float values
    frequencies_array = np.array(no_low_freq, dtype=float)
    
    return frequencies_array

def gauss_low_freq(file_path: str):
    """
    Reads in all frequences (units of cm^-1) from Gaussian log file and set anything less than 100 cm^-1 to 100 cm^-1

    Args:
        file_path (str): Path to gaussian frequency log file

    Returns:
        _type_: NumPy array of frequencies greater than or equal to 100 cm^-1 where everything less than 100 cm^-1 was set to 100 cm^-1
    """
    pattern = r'Frequencies\s+--\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)'  # regex to find desired lines
    frequencies = []

    with open(file_path, 'r') as file:
        for line in file:
            match = re.search(pattern, line)
            if match:
                freq_values = match.groups()
                frequencies.extend(freq_values)
    
    # turn strings to floats and adjust values < 100 to 100
    adjusted_freq = [float(x) if float(x) >= 100 else 100 for x in frequencies]
    
    # Convert the list to a NumPy array with float values
    frequencies_array = np.array(adjusted_freq, dtype=float)
    
    return frequencies_array


def gauss_rt_data(file_path:str):
    """
    Define patterns to match the lines with rotational temperatures, symmetry number, and molecular mass

    Args:
        file_path (string): Path to Gaussian frequency log file 

    Returns:
        _type_: symmetry number, rotational temperatures, and molecular mass needed for rotational partition function
    """
    temp_pattern = r"Rotational temperatures \(Kelvin\)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)"
    symmetry_pattern = r"Rotational symmetry number\s+(\d+)"
    mass_pattern = r"Molecular mass:\s+(\d+\.\d+)\s+amu\."
    
    # Initialize variables to store the extracted data
    rotational_temperatures = []
    symmetry_number = None
    molecular_mass = None
    
    # Open the file and read it line-by-line
    with open(file_path, 'r') as file:
        for line in file:
            # Search for the temperature pattern in the current line
            temp_match = re.search(temp_pattern, line)
            if temp_match:
                # Extract the temperatures from the matched line and convert them to floats
                temperatures = [float(temp_match.group(1)), float(temp_match.group(2)), float(temp_match.group(3))]
                rotational_temperatures.extend(temperatures)
            
            # Search for the symmetry pattern in the current line
            symmetry_match = re.search(symmetry_pattern, line)
            if symmetry_match:
                # Extract the symmetry number from the matched line and convert it to an integer
                symmetry_number = int(symmetry_match.group(1))
            
            # Search for the molecular mass pattern in the current line
            mass_match = re.search(mass_pattern, line)
            if mass_match:
                # Extract the molecular mass from the matched line and convert it to a float
                molecular_mass = float(mass_match.group(1))
    
    return symmetry_number, rotational_temperatures, molecular_mass
#-------------------------------------------------------------------------------- */
def vib_entropy(frequencies_cm1:np.ndarray, temperatures:np.ndarray):
    """
    Calculate vibrational entropy from NumPy array of frequencies (units of cm^-1) and temperature (Kelvin)
    Entropy is calculated using the last equation for Sv from https://gaussian.com/wp-content/uploads/dl/thermo.pdf

    Args:
        frequencies_cm1 (np.ndarray): NumPy array with frequency values in cm^-1
        temperatures (np.ndarray): NumPy array with temperatures to calculate entropy over

    Returns:
        _type_: NumPy array with vibrational entropy with units of J/mol*K
    """
    # Constants
    R = 8.314462618  # Gas constant in J/(mol*K)
    h = 6.62607015e-34  # Planck constant in Js
    k_B = 1.380649e-23  # Boltzmann constant in J/K
    c = 2.99792458e10  # Speed of light in cm/s
    
    # Convert frequencies from cm^-1 to Hz
    frequencies_hz = frequencies_cm1 * c
    
    # Convert frequencies and temperatures to numpy arrays if they are not already
    frequencies_hz = np.asarray(frequencies_hz)
    temperatures = np.asarray(temperatures)
    
    # Reshape frequencies for broadcasting
    frequencies_hz = frequencies_hz[:, np.newaxis]
    
    print(f'Frequencies shape:{frequencies_hz.shape}')

    # Calculate x = (h * nu) / (k_B * T) for all frequencies and temperatures
    x = (h * frequencies_hz) / (k_B * temperatures)

    # print(f'x shape:{x.shape}')
    
    # Calculate the vibrational entropy using broadcasting
    entropy = x / (np.exp(x) - 1) - np.log(1 - np.exp(-x))
    S_vib = R * np.sum(entropy, axis=0)
    
    return S_vib

def rot_partition(temp: np.ndarray, sym:float,  theta_one:float, theta_two: float, theta_three:float):
    """
    Rotational partition function calculated using the equation from https://gaussian.com/wp-content/uploads/dl/thermo.pdf
    for a non-linear polyatomic molecule 

    Args:
        temp (np.ndarray): NumPy array for temperatures to calculate rotational entropy over
        sym (float): Symmetry number for the molecule
        theta_one (float): First rotational temperature (order of assignment is arbitrary)
        theta_two (float): Second rotational temperature (order of assignment is arbitrary)
        theta_three (float): Third rotational temperature (order of assignment is arbitrary)

    Returns:
        _type_: NumPy array with rotational entropy as a function of temperature with units of J/mol*K
    """
    pi = np.pi 
    return np.array(np.sqrt(pi)/sym * (temp**(3/2) / np.sqrt(theta_one * theta_two * theta_three)))

def trans_partition(mass:float, temp:np.ndarray, pres:float) -> np.ndarray:
    """
    _summary_

    Args:
        mass (float): mass in amu
        temp (np.ndarray): temperature in kelvin
        pres (float): pressure in atm

    Returns:
        np.ndarray: translational partition function at specified temperature(s) with units of J/mol*K
    """
    # Constants
    k_B = 1.380649e-23  # Boltzmann constant in J/K
    h = 6.62607015e-34  # Planck constant in Js
    N_A = 6.02214076e23  # Avogadro's number in mol^-1
    R = 8.314462618  # Gas constant in J/(mol*K)
    
    # Convert mass from amu to kg
    mass_kg = mass * 1.66053906660e-27  # 1 amu = 1.66053906660e-27 kg
    
    # Convert pressure from atm to Pa
    pressure_pa = pres * 101325  # 1 atm = 101325 Pa
    
    # Calculate volume using ideal gas law: V = (N_A * k_B * T) / P
    volume = (N_A * k_B * temp) / pressure_pa
    
    # Calculate translation partition function
    q_trans = ((2 * np.pi * mass_kg * k_B * temp) / (h**2))**(3/2) * (volume / N_A)
    
    return np.array(q_trans)
#-------------------------------------------------------------------------------- */
# MASTER FUNCTION TO CALCULATE ENTROPY

def calc_entropy(infile:str, temperature:np.ndarray, use_low_freq = False) -> np.ndarray:
    """
    _summary_

    Args:
        infile (str): Path to Gaussian frequency log file
        temperature (np.ndarray): NumPy array (that can be 1D) of temperature values to calculate entropy 
        use_low_freq (bool, optional): Set to True if you want frequencies < 100 cm^-1 to be included and 
        set to 100 cm^-1 when reading in log file 

    Returns:
        np.ndarray: Return s_trans, s_vib, s_rot over specified temperature in units of J/mol*K
    """

    R = 8.314462618  # Gas constant in J/(mol*K)
    h = 6.62607015e-34  # Planck constant in Js
    k_B = 1.380649e-23  # Boltzmann constant in J/K
    c = 2.99792458e10  # Speed of light in cm/s
    

    symmetry_number, rotational_temperatures, molecular_mass = gauss_rt_data(infile)
    
    if use_low_freq == False:
        found_frequencies = gauss_freq(infile)
    else:
        found_frequencies = gauss_low_freq(infile)
    
    q_trans = trans_partition(temp = temperature, mass = molecular_mass, pres = 1)
    q_rot = rot_partition(temp = temperature, sym = symmetry_number, theta_one = rotational_temperatures[0], theta_two = rotational_temperatures[1], theta_three = rotational_temperatures[2])

    s_trans = R * (np.log(q_trans) + 5/2)
    s_rot = R * (np.log(q_rot) + 3/2)
    s_vib = vib_entropy(frequencies_cm1 = found_frequencies, temperatures = temperature)

    return s_trans, s_vib, s_rot