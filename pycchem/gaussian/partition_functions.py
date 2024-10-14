#/* --------------------------------------------------------------------------------
#   Paolucci Group
#   University of Virginia
#   Mohan Shankar
#
#   partition_functions.py
#
#   This file takes in data from log files of frequency calculations from Gaussian
#   and applies statistical mechanics to calculate partition functions 
#-------------------------------------------------------------------------------- */
import numpy as np
import re
#-------------------------------------------------------------------------------- */
# HELPER FUNCTIONS TO READ IN DATA
def gaussian_parse_frequencies(file_path:str):
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

import re
import numpy as np

def gaussian_parse_low_frequencies(file_path: str):
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


def extract_rt_data(file_path:str):
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
def vibrational_entropy(frequencies_cm1:np.ndarray, temperatures:np.ndarray):
    """
    Calculate vibrational entropy from NumPy array of frequencies (units of cm^-1) and temperature (Kelvin)
    Entropy is calculated using the last equation for Sv from https://gaussian.com/wp-content/uploads/dl/thermo.pdf

    Args:
        frequencies_cm1 (np.ndarray): NumPy array with frequency values in cm^-1
        temperatures (np.ndarray): NumPy array with temperatures to calculate entropy over

    Returns:
        _type_: NumPy array with vibrational entropy 
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