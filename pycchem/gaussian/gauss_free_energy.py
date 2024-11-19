import numpy as np
import re

# HELPER FUNCTIONS TO READ IN DATA
def gauss_freq_gauss(file_path:str):
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

def gauss_low_freq_gauss(file_path: str):
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


def gauss_rt_data_gauss(file_path:str):
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
def vib_entropy_gauss(frequencies_cm1:np.ndarray, temperatures:np.ndarray):
    """
    Calculate vibrational entropy from NumPy array of frequencies (units of cm^-1) and temperature (Kelvin)
    Entropy is calculated using the last equation for Sv from https://gaussian.com/wp-content/uploads/dl/thermo.pdf

    Args:
        frequencies_cm1 (np.ndarray): NumPy array with frequency values in cm^-1
        temperatures (np.ndarray): NumPy array with temperatures to calculate entropy over

    Returns:
        S_vib: NumPy array with vibrational entropy with units of J/mol*K
        ZPE: Float with zero-point energy with units of J/mol
        Ev: NumPy array of vibrational component of thermal correction to energy based on Gaussian approach in J/mol
    """
    # Constants
    R = 8.314462618  # Gas constant in J/(mol*K)
    h = 6.62607015e-34  # Planck constant in Js
    k_B = 1.380649e-23  # Boltzmann constant in J/K
    c = 2.99792458e10  # Speed of light in cm/s
    NA = 6.02214076e23  # Avogadro's number
    
    # Convert frequencies from cm^-1 to Hz
    frequencies_hz = frequencies_cm1 * c
    
    # Convert frequencies and temperatures to numpy arrays if they are not already
    frequencies_hz = np.asarray(frequencies_hz)
    temperatures = np.asarray(temperatures)
    
    ZPE = 1/2 * h * np.sum(frequencies_hz) * NA # ZPE in J/mol

    # Reshape frequencies for broadcasting
    frequencies_hz = frequencies_hz[:, np.newaxis]
    
    # print(f'Frequencies shape:{frequencies_hz.shape}')

    # Calculate x = (h * nu) / (k_B * T) for all frequencies and temperatures
    x = (h * frequencies_hz) / (k_B * temperatures)

    vib_temp = (h * frequencies_hz) / (k_B)

    Ev = R * np.sum(x * (0.5 + 1 / (np.exp(x) - 1)), axis = 0) # units of 

    # print(f"Thermal correction to vibrational energy: {Ev}")

    # print(f'x shape:{x.shape}')
    
    # Calculate the vibrational entropy using broadcasting
    entropy = x / (np.exp(x) - 1) - np.log(1 - np.exp(-x))
    S_vib = R * np.sum(entropy, axis=0)
    
    return S_vib, ZPE, Ev

def rot_partition_gauss(temp: np.ndarray, sym:float,  theta_one:float, theta_two: float, theta_three:float):
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

def trans_partition_gauss(mass:float, temperature:np.ndarray, pres:np.ndarray) -> np.ndarray:
    """
    _summary_

    Args:
        mass (float): mass in amu
        temp (np.ndarray): temperature in kelvin
        pres (float): pressure in atm

    Returns:
        np.ndarray: translational partition function at specified temperature(s) and pressure(s) 
        Constant temp. and increasing pressure along row; constant pressure and increasing temp. down column
        with units of J/mol*K
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

    if not isinstance(temperature, np.ndarray):
        raise TypeError(f"Expected temperature to be a NumPy array, but got {type(temperature).__name__}. Fix it by wrapping `np.array([temperature])`.")

    if temperature.ndim == 1:
        temperature = temperature[:, np.newaxis]
    
    # Calculate volume using ideal gas law: V = (N_A * k_B * T) / P
    volume = (N_A * k_B * temperature) / pressure_pa
    
    # Calculate translation partition function
    q_trans = ((2 * np.pi * mass_kg * k_B * temperature) / (h**2))**(3/2) * (volume / N_A)
    
    return np.array(q_trans)

def calc_entropy_gauss(infile:str, temperature:np.ndarray, pressure:np.ndarray, use_low_freq = False) -> np.ndarray:
    """
    _summary_

    Args:
        infile (str): Path to Gaussian frequency log file
        temperature (np.ndarray): NumPy array (that can be 1D) of temperature values to calculate entropy 
        pressure (np.ndarray): NumPy array of pressure values in units of atm
        use_low_freq (bool, optional): Set to True if you want frequencies < 100 cm^-1 to be included and 
        set to 100 cm^-1 when reading in log file 

    Returns:
        np.ndarray: Return s_trans, s_vib, s_rot, Ev (thermal correction to energy based on vibration) over specified temperature in units of [J/mol*K]
        If pressure is also varied (i.e. the array is 2D): 
            constant temp. and increasing pressure along row; constant pressure and increasing temp. down column
    """

    R = 8.314462618  # Gas constant in J/(mol*K)
    h = 6.62607015e-34  # Planck constant in Js
    k_B = 1.380649e-23  # Boltzmann constant in J/K
    c = 2.99792458e10  # Speed of light in cm/s
    

    symmetry_number, rotational_temperatures, molecular_mass = gauss_rt_data_gauss(infile)
    
    if use_low_freq == False:
        found_frequencies = gauss_freq_gauss(infile)
    else:
        found_frequencies = gauss_low_freq_gauss(infile)
    
    q_trans = trans_partition_gauss(temperature = temperature, mass = molecular_mass, pres = pressure)
    q_rot = rot_partition_gauss(temp = temperature, sym = symmetry_number, theta_one = rotational_temperatures[0], theta_two = rotational_temperatures[1], theta_three = rotational_temperatures[2])

    s_trans = R * (np.log(q_trans) + 5/2)
    s_rot = R * (np.log(q_rot) + 3/2)
    s_vib, ZPE, Ev = vib_entropy_gauss(frequencies_cm1 = found_frequencies, temperatures = temperature)

    return s_trans, s_vib, s_rot, ZPE, Ev

def last_gauss_energy_gauss(filename: str) -> float:
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

def calc_free_energy(file:str, temperature: np.ndarray, pressure: np.ndarray) -> np.ndarray:

    """
    Temperature must be in a NumPy array even if you have a single value!

    Returns:
        ∆H in kJ/mol, ∆S in J/mol*K, ∆G in kJ/mol
        These values should match NIST values if you perform 
        [{Hº(T) - Hº(Tr)} - {H(0) - Hº(Tr)}], Sº, and  [{Hº(T) - Hº(Tr)} - {H(0) - Hº(Tr)}] - TSº
    """

    R = 8.314 # J/mol*K

    har_to_kjmol = 2625.5
    
    s_trans, s_vib, s_rot, ZPE, Ev = calc_entropy_gauss(infile = file, temperature = temperature, pressure = pressure) # S [=] of J/mol*k, Ev [=] J/mol

    energy = last_gauss_energy_gauss(file)

    s_tot = s_trans + s_vib[:, np.newaxis] + s_rot[:, np.newaxis] # units of J/mol*K

    calculated_thermal_energy_corr = (3 * (R) + (Ev)) * temperature/(1000 *har_to_kjmol) # Er + Et + Ev in Hartree; should match Gaussian value for 
    # Thermal correction to Energy
    
    calculated_thermal_enthalpy_corr = calculated_thermal_energy_corr + (R * temperature)/(1000 * har_to_kjmol) # in Hartree; should match Gaussian value for
    # Thermal correction to Enthalpy

    delta_h = calculated_thermal_enthalpy_corr - (ZPE/(1000*har_to_kjmol)) # ∆H in Hartree; NIST H values DO NOT take zero-point energy into account 

    if s_tot.shape[1] == 1: # transpose from column to row vector to avoid weird broadcasting 
        delta_g = (delta_h * har_to_kjmol) - (temperature * (s_tot.T/1000)) + (8.314/1000) * temperature * np.log((pressure)) # ∆G in kJ/mol
        return (delta_h * har_to_kjmol), s_tot, delta_g
    
    elif s_tot.shape[1] > 1: # if more than one pressure value is passed
        delta_g = np.zeros((len(temperature), len(pressure)), dtype = np.float32) # constant temp along row; constant pres. down column )

        for t, temp in enumerate(temperature): # S_tot is 2D with constant temp along row; constant pres down column 
            for p, pres in enumerate(pressure): # S_tot has units of J/mol*K so div. by 1000 --> kJ
                delta_g[t][p] = (delta_h[t] * har_to_kjmol) - (temp * (s_tot[t][p])/1000) + (0.008314 * temp * np.log(pres))
        return (delta_h * har_to_kjmol), s_tot, delta_g