#/* --------------------------------------------------------------------------------
#   Paolucci Group
#   University of Virginia
#   Mohan Shankar
#
#   vasp_utils.py
#
#   This file takes in data from OUTCAR files
#-------------------------------------------------------------------------------- */
# DEPENDENCIES
import re
import numpy as np
#-------------------------------------------------------------------------------- */
def vasp_parse_frequencies(file_path:str)  -> np.ndarray:
    """
    Parse an OUTCAR file for all frequencies >= 100 cm^-1 

    Args:
        file_path (str): Path to OUTCAR

    Returns:
        np.ndarray: NumPy array with vibrational frequencies (units of cm^-1)
    """
    # Regular expression to match numbers followed by 'cm-1'
    pattern = r"(\d+\.\d+)\s*cm-1"
    
    # Initialize an empty list to store the extracted values
    cm1_values = []
    
    # Open and read the file
    with open(file_path, 'r') as file:
        for line in file:
            # Find all matches in the current line
            matches = re.findall(pattern, line)
            # Convert matches to floats and append to the list if the frequency is >= 100 cm^-1
            cm1_values.extend([float(match) for match in matches if float(match) >= 100])
    
    return np.array(cm1_values)

def vasp_parse_low_frequencies(file_path:str)  -> np.ndarray:
    """
    Parse an OUTCAR file for all frequencies >= 100 cm^-1 and set all remaining frequencies to 100 cm^-1

    Args:
        file_path (str): Path to OUTCAR

    Returns:
        np.ndarray: NumPy array with vibrational frequencies (units of cm^-1)
    """
    # Regular expression to match numbers followed by 'cm-1'
    pattern = r"(\d+\.\d+)\s*cm-1"
    
    # Initialize an empty list to store the extracted values
    cm1_values = []
    
    # Open and read the file
    with open(file_path, 'r') as file:
        for line in file:
            # Find all matches in the current line
            matches = re.findall(pattern, line)
            # Convert matches to floats and append to the list if the frequency is >= 100 cm^-1 else set the frequency to 100 cm^-1
            # cm1_values.extend([max(float(match), 100) for match in matches]) # smarter, pythonic line 
            cm1_values.extend([float(match) if float(match) >= 100 else 100 for match in matches])
    return np.array(cm1_values)

def last_vasp_energy(file_path:str) -> float:
    """
    _summary_

    Args:
        file_path (str): Path to OUTCAR

    Returns:
        float: Last energy reported (eV)
    """
    pattern = r"energy\(sigma->0\)\s*=\s*(-?\d+\.\d+)"
    last_energy = None
    
    # Read the file line by line
    with open(file_path, 'r') as file:
        for line in file:
            match = re.search(pattern, line)
            if match:
                last_energy = float(match.group(1))
    
    return last_energy