#/* --------------------------------------------------------------------------------
#   Paolucci Group
#   University of Virginia
#   Mohan Shankar
#
#   misc_utils.py
#
#   This file houses vaarious, random functions that I've written and found somewhat
#   useful 
#-------------------------------------------------------------------------------- */
# DEPENDENCIES
import numpy as np
from scipy.interpolate import CubicSpline 
import requests
from bs4 import BeautifulSoup
#-------------------------------------------------------------------------------- */
def cubicspline_interp(x_orig: np.ndarray, y_orig: np.ndarray, x_new: np.ndarray) -> np.ndarray:
    """
    Args:
        x_orig (np.ndarray): original x data to interpolate over
        y_orig (np.ndarray): original y data to interpolate over
        x_new (np.ndarray): new range of x values to interpolate over

    Returns:
        np.ndarray: _description_
    """
    interp = CubicSpline(x_orig, y_orig) 
    return interp(x_new)

def nist_collector(url: str) -> np.ndarray:
    """
    Args:
        url (str): NIST Janaf table for compound of choice

    Returns:
        np.array: array with temperature (T/K), entropy (S0), enthalpy (H - H0)
    """
    
    # Send a GET request to the webpage
    response = requests.get(url)
    response.raise_for_status()  # Check if the request was successful

    # Parse the HTML content using BeautifulSoup
    soup = BeautifulSoup(response.content, 'html.parser')

    # Find the table in the HTML
    table = soup.find('table')

    # Extract table headers
    headers = [header.text.strip() for header in table.find_all('th')]

    # Extract table rows
    rows= []
    for row in table.find_all('tr')[2:-4]:
        cols = row.find_all('td')
        cols = [ele.text.strip() for ele in cols]
        if cols:
            rows.append(cols)

    temp = [] # temperature values for compound
    enthalpy = [] # enthalpy values for compound
    entropy = [] # entropy values for compound

    for row in rows:
        if len(row) == 15: # make sure the line isn't blank
            temp.append(row[0])
            entropy.append(row[4])
            enthalpy.append(row[8])
    

    std_enthalpy_index = temp.index('298.15') # find index where the standard temperature is used

    # std_enthalpy = enthalpy[std_enthalpy_index] # find change in enthalpy associated with standard temp

    temp.remove('298.15') # remove standard temp 

    temp.append('298.15')

    enthalpy.append(enthalpy[std_enthalpy_index]) # add standard enthalpy to end of array

    enthalpy.remove(enthalpy[std_enthalpy_index]) # remove standard enthalpy value from middle of array

    entropy.append(entropy[std_enthalpy_index]) # add standard entropy to end of array

    entropy.remove(entropy[std_enthalpy_index]) # remove standard entropy value from middle of array

    temp = np.array(temp, dtype = np.float32) # convert temperature list to np.array w/ units of Kelvin

    entropy = np.array(entropy, dtype = np.float32)/1000 # units of J/K*mol --> kJ/K*mol

    enthalpy = np.array(enthalpy, dtype = np.float32) # units of kJ/mol

    return np.transpose(np.vstack((temp, entropy, enthalpy)))