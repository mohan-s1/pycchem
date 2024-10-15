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
def last_gauss_energy(filename:str) -> float:
    """
    Finds last energy (units of Hartree) from Gaussian Log file 

    Args:
        filename (str): Path to log file

    Returns:
        float: Last energy mentioned in Log file (units of Hartree)
    """
    pattern = r"SCF Done:  E\(UPBE-PBE\) =  (-?\d+\.\d+)     A\.U\. after    \d+ cycles"
    last_energy = None
    
    with open(filename, 'r') as file:
        for line in file:
            match = re.search(pattern, line)
            if match:
                last_energy = float(match.group(1))
    
    return last_energy