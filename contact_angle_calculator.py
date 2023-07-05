import sys
import numpy as np
from iocontrol.options import get_options as getopt
import MDAnalysis as md
from scipy.optimize import curve_fit
from scipy.interpolate import UnivariateSpline

def linear_func(x, m, b):
    return x * m + b

def get_empty_radial_density_histogram(rnum, znum):
	rdens = [np.zeros(rnum) for _ in range(znum)]
	return rdens

def parse_contact_angle_arguments():
    """
    Parse command-line arguments for contact angle calculation script.
    
    Command-line options
    --------------------
    -f : str, optional
        Input XTC trajectory file, default='run_hydrophilic_production_layer.xtc'
    -s : str, optional
        Input TPR topology file, default='run_hydrophilic_production_layer.tpr'
    -o : str, optional
        Output file prefix, default='zbin'
    -o2 : str, optional
        Output file prefix and type of surface, default='hydrophilic'
    -bin : float, optional
        Bin width for z-axis (in Angstrom), default=4.0
    -b : int, optional
        First frame to analyze, default=0
    -e : int, optional
        Last frame to analyze, default=-1 (meaning it will analyze until the last frame)
    -s : int, optional
        Number of frames to step_interval between analyses, default=1 (analyze every frame)

    Returns
    -------
    tuple
        A tuple containing the values of the parsed command-line arguments.
    """
    options = ['-f', '-s', '-o', '-o2', '-bin', '-b', '-e', '-s']
    types = ['str', 'str', 'str', 'str', 'float', 'int', 'int', 'int']
    defaults = ['run_hydrophilic_production_layer.xtc', 'run_hydrophilic_production_layer.tpr', 'zbin', 'hydrophilic', 4.0, 0, -1, 1]
    return getopt(sys.argv, options, types, defaults)

def calculate_num_density_water_atoms_spce_300K_per_cubic_angstrom():
	"""
	Calculate the number density of SPC/E water_atoms at 300 K in particles per cubic angstrom (Angstrom^-3).

	Returns
	-------
	float
	    Number density of SPC/E water_atoms at 300 K in particles per cubic angstrom (Angstrom^-3), scaled with the Lennard-Jones
	    parameter sigma.
	"""
	# Initial density of SPC/E water_atoms at 300 K in kg per liter (kg/L); this is a common reference value
	density_water_atoms_spce_300K_kg_per_L = 1.008

	# Convert density from kg/L to g/L; more convenient units for subsequent calculations
	density_water_atoms_spce_300K_g_per_L = density_water_atoms_spce_300K_kg_per_L * 1000

	# Molar mass of water_atoms in g/mol; needed to convert mass density to molar concentration
	molar_mass_water_atoms = 18.01528

	# Convert mass density to molar concentration in mol/L; enables conversion to number density
	concentration_water_atoms_spce_300K_mol_per_L = density_water_atoms_spce_300K_g_per_L / molar_mass_water_atoms

	# Avogadro's number: number of particles per mol; essential for converting moles to particles
	avogadro_number = 6.02214076 * 10**23

	# Conversion factor: 1 m^3 equals 10^9^3 nm^3; allows conversion from m^3 to nm^3
	cubic_meter_to_cubic_nanometer = (10**9)**3

	# Convert molar concentration to number density in particles per cubic nanometer (nm^-3)
	# This unit conversion allows for comparison with molecular dynamics simulation results
	num_density_water_atoms_spce_300K_per_nm3 = (
		concentration_water_atoms_spce_300K_mol_per_L
		* avogadro_number
		/ cubic_meter_to_cubic_nanometer
	)

	# Conversion factor: 1 nm equals 10 Angstrom; facilitates conversion from nm to Angstrom
	nanometer_to_angstrom = 10

	# Lennard-Jones parameter sigma for SPC/E water_atoms model in Angstrom; used for scaling the number density
	sigma_spce = 3.2

	# Convert number density from particles/nm^3 to particles/Angstrom^3 and scale with sigma^3
	# This step provides number density in units convenient for molecular dynamics simulations
	num_density_water_atoms_spce_300K_per_angstrom3 = (
		num_density_water_atoms_spce_300K_per_nm3
		/ (nanometer_to_angstrom**3)
		* (sigma_spce**3)
	)
    
	return num_density_water_atoms_spce_300K_per_angstrom3 * 1000

def index_water_atoms_to_z_layers_in_sessile_droplet(water_atoms, surface_z_position, zbin, znum):
	"""
	Assign each water_atoms atom to the corresponding z-layer in the sessile droplet.
		
	Parameters:
	----------
	water_atoms : MDAnalysis.core.groups.AtomGroup
		water_atoms atoms in the sessile droplet.
	surface_z_position : float
		The z-height of the first mono layer above the surface.
	zbin : float
		The width of the z-layer bins.
	znum : int
		The number of z-layers.

	Returns:
	-------
	atomid_in_layers : list of lists
		List of lists containing atom indices in each z-layer.
	"""
      
	atomid_in_layers = [[] for _ in range(znum)]
	for atom in water_atoms:
		if(atom.position[2] >= surface_z_position):
			index = int((atom.position[2] - surface_z_position) / zbin)
			atomid_in_layers[index].append(atom.index)
	return atomid_in_layers

def calculate_centers_of_mass_of_z_layer_in_sessile_droplet(univ, atomid_in_layers, znum):
	"""
	Calculate the centers of mass of each z-layer in the sessile droplet.

	Parameters:
	----------
	univ : MDAnalysis.core.universe.Universe
		MDAnalysis Universe containing the simulation data.
	atomid_in_layers : list of lists
		List of lists containing atom indices of water_atoms in each z-layer of the sessile droplet.
	znum : int
		The total number of z-layers.

	Returns:
	-------
	coms : list of numpy.ndarray
		List of lists containing the center of mass coordinates for each z-layer.
	"""

	coms = [[] for _ in range(znum)]
	for layer in range(znum):
		coms[layer].append(univ.atoms[atomid_in_layers[layer]].center_of_mass())
	return coms

def calculate_atom_layer_COM_distances_sessile_droplet(univ, atomid_in_layers, coms):
	"""
	Calculate the x distances of atoms to the center of mass (COM) in each z-layer of the sessile droplet.

	Parameters:
	----------
	univ : MDAnalysis.core.universe.Universe
		MDAnalysis Universe containing the simulation data.
	atomid_in_layers : list of lists
		List of lists containing atom indices in each z-layer.
	coms : list of numpy.ndarray
		List of lists containing the center of mass coordinates for each z-layer.

	Returns:
	-------
	atom_COM_distances_in_z_layer : list of lists
		A list containing tuples with atom indices and their distances to the COM for each z-layer.
	"""

	atom_COM_distances_in_z_layer = []
	for layer in range(len(atomid_in_layers)):
		layer_distances = []
		if len(atomid_in_layers[layer]) >= 1:
			for atom in atomid_in_layers[layer]:
				dist = np.abs(univ.atoms[atom].position[0] - coms[layer][0][0])
				layer_distances.append((atom, dist))
		atom_COM_distances_in_z_layer.append(layer_distances)
	return atom_COM_distances_in_z_layer

def update_central_z_layer_density(univ, atom_COM_distances_in_z_layer, surface_z_position, zbin, zcentdens, rbin):
	"""
	Update the central z-layer density based on atom distances to the center of mass in each z-layer.

	Parameters
	----------
	atom_COM_distances_in_z_layer : list of lists
		A list containing tuples with atom indices and their distances to the COM for each z-layer.
	surface_z_position : float
		The z-height of the first mono layer above the surface.
	zbin : float
		The width of the z-layer bins.
	zcentdens : numpy.ndarray
		A 1D array containing the central density of each z-layer.
	rbin : float
		The radial bin size used for density calculations.

	Returns
	-------
	zcentdens : numpy.ndarray
		The updated 1D array containing the central density of each z-layer.
	"""

	for layer in range(len(atom_COM_distances_in_z_layer)):
		for atom, dist in atom_COM_distances_in_z_layer[layer]:
			if dist <= 4 * rbin: # use x distance smaller than 0.8 nm from center of mass to compute the average central density of the droplet
				centindex = int((univ.atoms[atom].position[2] - surface_z_position) / zbin)
				zcentdens[centindex] += 1.0
	return zcentdens

def calculate_z_density_and_central_density(univ, water_atoms, surface_z_position, zbin, rbin, spce_numdens, end_step, step_interval, num_steps):
	"""
	Calculate the z-density and central density of water_atoms molecules in the sessile droplet.

	Parameters
	----------
	univ : MDAnalysis.core.universe.Universe
		MDAnalysis Universe containing the simulation data.
	water_atoms : AtomGroup
		Atom group containing the water_atoms molecules.
	surface_z_position : float
		The z-height of the first mono layer above the surface.
	zbin : float
		The width of the z-layer bins.
	rbin : float
		The radial bin size used for density calculations.
	spce_numdens : float
		The number density of SPC/E water_atoms at 300 K in particles per cubic angstrom.
	end_step : int
		The last step to consider in the trajectory.
	step_interval : int
		The number of steps to step_interval between analysis points.
	num_steps : int
		The total number of steps for analysis.

	Returns
	-------
	zdensity : numpy.ndarray
		A 1D array containing the z-density of water_atoms molecules in the sessile droplet.
	zcentdens : numpy.ndarray
		A 1D array containing the central density of water_atoms molecules in the sessile droplet.
	"""
    
	box = univ.dimensions
	znum = int((box[2]) / zbin) + 1
	zdensity = np.zeros(znum)
	zcentdens = np.zeros(znum)
	count = 0

	for i in univ.trajectory:
		if i.frame >= end_step / 2 and i.frame % step_interval == 0 and i.frame < end_step:
			count += 1

			atomid_in_layers = index_water_atoms_to_z_layers_in_sessile_droplet(water_atoms, surface_z_position, zbin, znum)
			coms = calculate_centers_of_mass_of_z_layer_in_sessile_droplet(univ, atomid_in_layers, znum)
			atom_COM_distances_in_z_layer = calculate_atom_layer_COM_distances_sessile_droplet(univ, atomid_in_layers, coms)
			zcentdens = update_central_z_layer_density(univ, atom_COM_distances_in_z_layer, surface_z_position, zbin, zcentdens, rbin)

			index = ((water_atoms.positions[:, 2] - surface_z_position) / zbin).astype(int)
			zdens = np.bincount(index[index >= 0], minlength=znum).astype(int)
			zdens = zdens.astype(float) / (box[0] * box[1] * zbin * spce_numdens)
			zdensity += zdens / num_steps

	zcentdens /= (4 * rbin) * zbin * box[1] * count  # number density per angstrom cubed

	return zdensity, zcentdens

def calculate_water_air_interface_density(zbin_preprocessing, zcentdens):
	"""
	Calculate the average central density for the region between 1 nm and 2 nm above the surface, 
	then divide it by 2 to determine the density at the water-air interface.

	Parameters
	----------
	zbin_preprocessing : float
		The size of the z-bin used for pre-processing.
	zcentdens : numpy.ndarray
		A 1D array that contains the central density data.

	Returns
	-------
	float
		The density at the water-air interface.

	"""
	# Determine the indices of z-layers at 1 nm and 2 nm above the surface
	indexpos_1nm_layer = int(10.0 / zbin_preprocessing) # index of the z-layer at 1 nm above the surface
	indexpos_2nm_layer = int(20.0 / zbin_preprocessing) # index of the z-layer at 2 nm above the surface, after the first atomic layer

	# Calculate the average central density for the region between 1 nm and 2 nm above the surface and divide it by 2
	water_air_interface_density = np.mean(zcentdens[indexpos_1nm_layer:indexpos_2nm_layer]) / 2.0

	return water_air_interface_density



def update_radial_density_histogram(atom_COM_distances_in_z_layer, rbin, zbin, box, averaging_window, rdens):
	"""
	Calculates the radial density histogram of atoms in each z-layer of a droplet, given their center of mass 
	(COM) distances to the corresponding z-layer COM. Updates the provided histogram with the new densities. 

	Parameters
	----------
	atom_COM_distances_in_z_layer : list
		List containing the COM distances of atoms in each z-layer of the droplet.
	rbin : float
		The radial bin width used for calculating the radial density histogram.
	zbin : float
		The z bin width of the z-layer
	box : list
		A list containing the dimensions of the simulation box.
	averaging_window : int
		The number of frames used for averaging the radial density histogram.
	rdens : list
		A list containing the radial density histograms for each z-layer of the droplet.

	Returns
	-------
	rdens : list
		A list containing the updated radial density histograms for each z-layer of the droplet.
	"""

	for layer in range(len(atom_COM_distances_in_z_layer)):
		for atom, dist in atom_COM_distances_in_z_layer[layer]:
				rdens[layer][int(dist / rbin)] += 1.0 / rbin / zbin / box[1] / averaging_window
	return rdens

def find_linear_interpolated_interface_location(func, rdens, zlayer, water_air_interface_density, rbin):
	"""
	Calculates the exact location of the water-air interface in a given z-layer of the droplet, 
	by fitting a linear function to the radial density histogram in that layer and finding the 
	intersection with the water-air interface density. 

    Parameters
    ----------
    func : function
        A function used for fitting the radial density histogram in the given z-layer.
    rdens : list
        A list containing the radial density histograms for each z-layer of the droplet.
    zlayer : int
        The index of the z-layer in which to find the water-air interface.
    water_air_interface_density : float
        The density value at the water-air interface computed with calculate_z_density_and_central_density().
    rbin : float
        The radial bin width used for calculating the radial density histogram.

    Returns
    -------
    interface_location : float
        The linearly interpolated (exact) location of the water-air interface in the given z-layer.

	"""

	# Find the index where the radial density changes from being greater than 
	# or equal to the water-air interface density to being less than it
	
	interface_index = np.argmax(rdens[zlayer] < water_air_interface_density)
	# Prepare two data points for curve fitting
	density_values_to_fit = rdens[zlayer][interface_index-1:interface_index+1]
	radial_positions_to_fit = np.array([(interface_index - 1) * rbin, interface_index * rbin])

	# Perform curve fitting using the provided function (linear_func) and the two data points
	fitted_params, fitted_covariance = curve_fit(func, radial_positions_to_fit, density_values_to_fit)
	
	# Extract the slope and intercept from the fitted parameters
	slope, intercept = fitted_params[0], fitted_params[1]
	
	# Calculate the exact location of the water-air interface by finding the intersection of 
	# the fitted line with the water_air_interface_density
	interface_location = (water_air_interface_density - intercept) / slope
	
	return interface_location

def select_coordinates_and_parameters(x_coordinates, z_coordinates, prefix):
	"""
	Selects the coordinates and parameters for the contact angle calculation based 
	on the curvature of the sessile droplet.

    Parameters
    ----------
    x_coordinates : list
        A list containing the x coordinates of the interface.
    z_coordinates : list
        A list containing the z coordinates of the interface.
    prefix : str
        Prefix to specify the droplet shape. For a convex droplet, prefix is "hydrophilic", 
		and for a concave droplet, prefix is "hydrophobic".

    Returns
    -------
    selected_interface_coordinates : numpy.ndarray
        A NumPy array containing the selected x and z coordinates of the interface.
    phasefactor : float
        A
	"""

	if prefix == "hydrophilic": # If the surface is hydrophilic and contact angle is expected to be < 90 degrees
		xb = x_coordinates[int(len(x_coordinates)/2):]
		yb = z_coordinates[int(len(z_coordinates)/2):]
		selected_interface_coordinates = np.array(list(zip(xb,yb)))
		selected_interface_coordinates = selected_interface_coordinates[np.where(selected_interface_coordinates[:,1]<20.0)]
		phasefactor = 0
		k = 3 # Convex droplet is best modeled cubic
	else:
		xb = x_coordinates[:int(len(x_coordinates)/2)-1]
		yb = z_coordinates[:int(len(z_coordinates)/2)-1]
		selected_interface_coordinates = np.array(list(zip(xb,yb)))
		selected_interface_coordinates = selected_interface_coordinates[:np.argmax(selected_interface_coordinates[:,0])-1]
		phasefactor = 90
		k = 2 # Concave droplet is best modeled quadratic

	return selected_interface_coordinates, phasefactor, k

def create_spline_weights(selected_interface_coordinates):
	"""
	Creates weights for the spline interpolation.

	Parameters:
	-----------
	selected_interface_coordinates : ndarray, shape=(n, 2)
		The selected coordinates of the droplet interface.

	Returns:
	--------
	spline_weights : ndarray, shape=(n,)
		The weights for the spline interpolation.
	"""

	spline_weights = np.copy(-selected_interface_coordinates[:,0])
	spline_weights = np.full(len(selected_interface_coordinates[:,0]), 1.0)
	spline_weights[0:2] *= 1.25
	spline_weights = spline_weights / np.sum(spline_weights)
	return spline_weights

def get_spline_interpolations(selected_interface_coordinates, spline_weights, spline_smoothing_factor, k):
	"""
	Calculates the first-order derivative of the spline interpolation and returns its value at the first point.

    Parameters:
    -----------
    selected_interface_coordinates : ndarray, shape=(n, 2)
        The selected coordinates of the droplet interface.
	spline_weights : ndarray, shape=(n,)
		A NumPy array containing the weights for each of the selected interface coordinates.
	spline_smoothing_factor : float
		A float representing the smoothing factor for the spline.
	k : int
		An integer representing the degree of the spline interpolation.

	Returns:
	--------
	x_splp1 : float
		 A float representing the difference between the x coordinates 
		 of the first two points in the selected interface coordinates.
	y_splp1 : float
		A float representing the difference between the y coordinates 
		of the first two points in the selected interface coordinates.
	"""

	y_spl = UnivariateSpline(selected_interface_coordinates[:,0], selected_interface_coordinates[:,1],
		s=spline_smoothing_factor, w=spline_weights,k=k)
	x_splp1 = selected_interface_coordinates[1,0] - selected_interface_coordinates[0,0]
	y_splp1 = y_spl.__call__(selected_interface_coordinates[1,0]) - y_spl.__call__(selected_interface_coordinates[0,0])
	return x_splp1, y_splp1


def calculate_contact_angle_spline(selected_interface_coordinates, phasefactor, k):
	"""
	Calculates the contact angle of the sessile droplet at the water-air interface
	using a cubic spline interpolation.

	Parameters:
	-----------
	selected_interface_coordinates : numpy.ndarray
		An array of x,y coordinates representing the water-air interface of the sessile droplet.
	phasefactor : float
		Phase factor to adjust the contact angle calculation (depends on curvature).
	k : int
		Order of the spline interpolation. For a convex droplet, k=3 (cubic), and for a concave 
		droplet, k=2 (quadratic).

	Returns:
	--------
	float
		The calculated contact angle in degrees.
	"""

	spline_smoothing_factor = len(selected_interface_coordinates[:, 0]) * 100
	spline_weights = create_spline_weights(selected_interface_coordinates)
	x_splp1, y_splp1 = get_spline_interpolations(selected_interface_coordinates, spline_weights, spline_smoothing_factor, k)
	contact_angle = phasefactor + np.arctan(y_splp1/x_splp1) * 180.0 / np.pi 
	return contact_angle

def calculate_contact_angle_linear(interface_radial_positions, zbin):
	"""
	Calculates the contact angle of the sessile droplet at the water-air interface
	using linear interpolation.

	Parameters:
	-----------
	interface_radial_positions : list
		A list of the radial position of the water-air interface at each z-layer.
	zbin : float
		The size of the z-bins

	Returns:
	--------
	float
		The calculated contact angle in degrees.
	"""

	if interface_radial_positions[0] != 0 and interface_radial_positions[1] != 0:
		contact_angle = np.arctan2(zbin, interface_radial_positions[0] - interface_radial_positions[1])*180.0/np.pi
	return contact_angle

def process_averaging_window(rdens, znum, rbin, zbin, water_air_interface_density, linear_func, contact_angles_spline, contact_angles_linear, prefix):
	"""
	Processes the averaging window by finding the interface location, selecting the relevant coordinates and parameters,
	and calculating the contact angles using both spline interpolation and linear interpolation.
    
	Parameters:
	-----------
	rdens : numpy.ndarray
		The radial density histogram of the sessile droplet.
	znum : int
		The number of z-layers in the histogram.
	rbin : float
		The size of the radial bins in the histogram.
	zbin : float
		The size of the z-bins in the histogram.
	water_air_interface_density : float
		The density at the water-air interface computed with calculate_z_density_and_central_density().
	linear_func : function
		A function to use for linear fitting and interpolation.
	contact_angles_spline : list
		A list to append the calculated contact angle using spline interpolation.
	contact_angles_linear : list
		A list to append the calculated contact angle using linear interpolation.
	prefix : str
		Prefix to specify the droplet shape. For a convex droplet, prefix is "hydrophilic", 
		and for a concave droplet, prefix is "hydrophobic".

	Returns:
	--------
	contact_angles_spline : list
		The list of calculated contact angles using spline interpolation.
	contact_angles_linear : list
		The list of calculated contact angles using linear interpolation.
	"""

	# Initialize an empty list to store the radial positions of the water-air interface
	interface_radial_positions = []

	# Iterate over each z-layer in the radial density histogram
	for zlayer in range(znum):
		# If the first radial density bin is inside the droplet, find the exact interface location (the "edge" of 
		# the droplet) using a linear interpolation function and store it in the interface_radial_positions list.
		if rdens[zlayer][0] >= water_air_interface_density:
			interface_location = find_linear_interpolated_interface_location(linear_func, rdens, zlayer, water_air_interface_density, rbin)
			interface_radial_positions.append(interface_location)

	# Create arrays for x and z coordinates based on the interface radial positions
	z_coordinates = np.arange(zbin / 2, len(interface_radial_positions) * zbin, zbin)
	z_coordinates = np.array((z_coordinates,z_coordinates)).reshape(1,-1)[0]
	x_coordinates = interface_radial_positions
	x_coordinates = np.array((np.array(x_coordinates), -np.array(x_coordinates))).reshape(1, -1)[0]

	# Select relevant coordinates and parameters based on the droplet shape prefix
	selected_interface_coordinates, phasefactor, k = select_coordinates_and_parameters(x_coordinates, z_coordinates, prefix)

	# Calculate the contact angle using spline interpolation
	contact_angle_spline = calculate_contact_angle_spline(selected_interface_coordinates, phasefactor, k)
	contact_angles_spline.append(contact_angle_spline)
	print("Spline contact angle: ", contact_angle_spline)
	
	# Calculate the contact angle using linear interpolation
	contact_angle_linear = calculate_contact_angle_linear(interface_radial_positions, zbin)
	contact_angles_linear.append(contact_angle_linear)
	print("Linear contact angle: ", contact_angle_linear)
	
	# Return the updated lists of contact angles
	return contact_angles_spline, contact_angles_linear


def process_frame(univ, water_atoms, surface_z_position, zbin, znum, rbin, rdens, averaging_window):
	"""
	Processes a single frame of the simulation trajectory by calculating the 
	COMs of the z-layers in the sessile droplet, calculating the x distances of 
	atoms to the center of mass (COM) in each z-layer of the sessile droplet, and
	updating the radial density histogram in each z-layer of a droplet with the new data.

	Parameters
	----------
	univ : MDAnalysis.Universe
		The MDAnalysis Universe object containing the simulation trajectory data.
	water_atoms : MDAnalysis.AtomGroup
		An AtomGroup object containing the water atoms in the simulation.
	surface_z_position : float
		The z-height of the first mono layer above the surface.
	zbin : float
		The z bin width
	znum : int
		The number of z-layers
	rbin : float
		The radial bin width used for calculating the radial density histogram.
	rdens : list
		A list containing the radial density histograms for each z-layer of the droplet.

	Returns
	-------
	rdens : list
		A list containing the updated radial density histograms for each z-layer of the droplet.

	"""

	box = univ.dimensions
	atomid_in_layers = index_water_atoms_to_z_layers_in_sessile_droplet(water_atoms, surface_z_position, zbin, znum)
	coms = calculate_centers_of_mass_of_z_layer_in_sessile_droplet(univ, atomid_in_layers, znum)
	atom_COM_distances_in_z_layer = calculate_atom_layer_COM_distances_sessile_droplet(univ, atomid_in_layers, coms)
	rdens = update_radial_density_histogram(atom_COM_distances_in_z_layer, rbin, zbin, box, averaging_window, rdens)
	return rdens


def calculate_contact_angles(univ, water_atoms, surface_z_position, water_air_interface_density, zbin, znum, rbin, contact_angles_spline, contact_angles_linear, prefix, start_step, step_interval, end_step, averaging_window, rnum):
    """
    Analyze simulation data by computing the radial density histogram and 
	calculating contact angles for each frame of the simulation trajectory.

    Parameters:
    -----------
    univ : Universe
        MDAnalysis Universe object containing the simulation trajectory.
    water_atoms : AtomGroup
        MDAnalysis AtomGroup object containing water atoms.
    surface_z_position : float
        The z-coordinate of the first monolayer above the surface of the substrate.
    zbin : float
        The size of the z-bins
    znum : int
        The number of z-layers
    rbin : float
        The size of the radial bins in the radial density histogram.
    rdens : numpy.ndarray
        The radial density histogram of the sessile droplet.
    contact_angles_spline : list
        A list to append the calculated contact angle using spline interpolation.
    contact_angles_linear : list
        A list to append the calculated contact angle using linear interpolation.
    prefix : str
        Prefix to specify the droplet shape. For a convex droplet, prefix is "hydrophilic", 
        and for a concave droplet, prefix is "hydrophobic".
    start_step : int, optional
        The first frame to analyze. Default is 0.
    step_interval : int, optional
        The interval between analyzed frames. Default is 1.
    end_step : int, optional
        The last frame to analyze.
    averaging_window : int, optional
        The number of frames to average over. Default is 50.
    rnum : int, optional
        The number of radial bins in the radial density histogram.

    Returns:
    --------
    contact_angles_spline : list
        The list of calculated contact angles using spline interpolation.
    contact_angles_linear : list
        The list of calculated contact angles using linear interpolation.
    """
    
    frame_count = 0
    rdens = get_empty_radial_density_histogram(rnum, znum)
    for frame in univ.trajectory:
        if frame.frame > start_step and frame.frame % step_interval == 0 and frame.frame < end_step:
            rdens = process_frame(univ, water_atoms, surface_z_position, zbin, znum, rbin, rdens, averaging_window)
            frame_count += 1

            if frame_count == averaging_window:
                print("time: ", frame.time, " ps")
                contact_angles_spline, contact_angles_linear = process_averaging_window(rdens, znum, rbin, zbin, water_air_interface_density, linear_func, contact_angles_spline, contact_angles_linear, prefix)
                frame_count = 0
                rdens = get_empty_radial_density_histogram(rnum, znum)
                
    return contact_angles_spline, contact_angles_linear

def get_contact_angles():
	# Parse command-line arguments and assign values to respective variables
	xtc, tpr, output, prefix, user_zbin, start_step, end_step, step_interval = parse_contact_angle_arguments()

	# Calculate the number density of SPC/E water_atoms at 300 K in particles per cubic angstrom
	spce_numdens = calculate_num_density_water_atoms_spce_300K_per_cubic_angstrom()

	# Create an MDAnalysis Universe object using the input topology and trajectory files and assign the box dimensions
	univ = md.Universe(tpr, xtc)
	box = univ.dimensions

	# Select water_atoms oxygen atoms (type OW) from the Universe and store them in an AtomGroup
	water_atoms = univ.select_atoms("type OW")

	# Set the end_step to the total number of frames in the trajectory if it is set to -1
	if end_step == -1:
		end_step = len(univ.trajectory) - 1

	# Calculate the total number of steps for analysis
	num_steps = (end_step - start_step) / step_interval + 1

	# Define averaging window, radial bin size, z-bin size, surface height, surface z-height, and related parameters
	averaging_window = 50
	averaging_window_in_ps = (univ.trajectory[start_step + 1].time - univ.trajectory[start_step].time) * averaging_window 
	rbin = 2.0
	zbin_preprocessing = 2.0
	surface_height = 20.0
	CH2_water_lj_sigma = 3.57 # Lennard-Jones sigma value between GROMOS 54A7 CH2 and OW water
	surface_z_position = surface_height + CH2_water_lj_sigma + zbin_preprocessing / 2.0

	# Calculate the z-density and central density of water_atoms molecules in the sessile droplet
	zdensity, zcentdens = calculate_z_density_and_central_density(
		univ, water_atoms, surface_z_position, zbin_preprocessing, rbin, spce_numdens, end_step, step_interval, num_steps
	)

	# Save the central density data to an output file with the specified output file name
	output_file_name = prefix + "zcentdens_" + output + ".xvg"
	np.savetxt(output_file_name, zcentdens, "%12.12e", comments='#', header="z-density profile in units of per Angstrom cubed through the droplet center") 

	# Load the saved central density data into a numpy array
	zcentdens = np.loadtxt(output_file_name, comments='#')

	# Calculate the interface density (average central density for the region between 1 nm and 2 nm above the surface divided by 2)
	water_air_interface_density = calculate_water_air_interface_density(zbin_preprocessing, zcentdens)


	zbin = user_zbin
	znum = int((box[2]) / zbin) + 1
	rmax = (box[0] * box[0] + box[1] * box[1]) ** 0.5
	rnum = int(rmax / rbin) + 1

	contact_angles_spline, contact_angles_linear = [], []
	contact_angles_spline, contact_angles_linear = calculate_contact_angles(univ, water_atoms, surface_z_position, water_air_interface_density, zbin, znum, rbin, contact_angles_spline, contact_angles_linear, prefix, start_step, step_interval, end_step, averaging_window, rnum)

	contact_angles_linear = np.array(list(zip(np.arange(univ.trajectory[start_step].time+averaging_window_in_ps/2.0,univ.trajectory[end_step].time,averaging_window_in_ps),contact_angles_linear)))
	contact_angles_spline = np.array(list(zip(np.arange(univ.trajectory[start_step].time+averaging_window_in_ps/2.0,univ.trajectory[end_step].time,averaging_window_in_ps),contact_angles_spline)))

	print("Mean linear contact angle: ", np.mean(contact_angles_linear[:,1]))
	print("Mean spline contact angle: ", np.mean(contact_angles_spline[:,1]))

	np.savetxt("contact_angles_linear_"+prefix+output+".xvg", contact_angles_linear, "%12.12e", comments='#', header="contact angle with linear extrapolation. timestep in ps") 
	np.savetxt("contangle_spline_"+prefix+output+".xvg", contact_angles_spline, "%12.12e", comments='#', header="contact angle with spline fitting. timestep in ps") 


if __name__ == "__main__":
    get_contact_angles()
