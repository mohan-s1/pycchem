#/* --------------------------------------------------------------------------------
#   Paolucci Group
#   University of Virginia
#   Mohan Shankar
#
#   gauss_utils.py
#
#   This file has various utility functions for processing Gaussian data
#-------------------------------------------------------------------------------- */
# DEPENDENCIES
import re
import numpy as np
#-------------------------------------------------------------------------------- */
def last_gauss_energy(filename: str) -> float:
    """
    Finds the last energy (in Hartree) from a Gaussian Log file.

    Args:
        filename (str): Path to the log file.

    Returns:
        float: Last energy mentioned in the log file (in Hartree) or 0.0 if not found.
    """
    # Adjust the pattern to allow for varying spaces and match the format in the log file
    pattern = r"SCF Done:\s*.*=\s*(-?\d+\.\d+)\s+[aA]\.U\.\s+after"
    last_energy = None
    
    with open(filename, 'r') as file:
        for line in file:
            match = re.search(pattern, line)
            if match:
                last_energy = float(match.group(1))
    
    # Return 0.0 if no match was found
    return last_energy if last_energy is not None else 0.0