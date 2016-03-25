#######################################################################
#
#  UNIFIED ANALYSIS SCRIPT FOR DISK SIMULATION FOR THE AGORA PROJECT
#
#  FOR SCRIPT HISTORY SEE VERSION CONTROL CHANGELOG
#
#  Note: This script requires yt-3.2 or yt-dev. Older versions may 
#        yield incorrect results. 
#
#######################################################################
import matplotlib
matplotlib.use('Agg')
import sys
import math
import numpy as np
import matplotlib.colorbar as cb
import matplotlib.lines as ln
import yt.utilities.physical_constants as constants
import matplotlib.pyplot as plt
from yt.mods import *
from yt.units.yt_array import YTQuantity
from mpl_toolkits.axes_grid1 import AxesGrid
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredText
from yt.fields.particle_fields import add_volume_weighted_smoothed_field
from yt.fields.particle_fields import add_nearest_neighbor_field
mylog.setLevel(1)

#file_location = '../../AGORA-DISK-repository-for-use/Grackle+noSF/'
file_location = '../../AGORA-DISK-repository-for-use/Grackle+SF+ThermalFbck/'
#file_location = '/global/project/projectdirs/agora/AGORA-DISK-repository-for-use/Grackle+noSF/'
#file_location = '/global/project/projectdirs/agora/AGORA-DISK-repository-for-use/Grackle+SF+ThermalFbck/'

codes = ['ART-I', 'ART-II', 'CHANGA', 'ENZO', 'GADGET-3', 'GASOLINE', 'GEAR', 'GIZMO', 'RAMSES']
# filenames = [[file_location+'ART-I/IC/AGORA_Galaxy_LOW.d', file_location+'ART-I/t0.5Gyr/10MpcBox_csf512_02350.d'],
# 	     [file_location+'ART-II/noSF_def/OUT/AGORA_LOW_000000.art', file_location+'ART-II/noSF_def/OUT/AGORA_LOW_000087.art'],
# 	     [file_location+'CHANGA/disklow/disklow.000000', file_location+'CHANGA/disklow/disklow.000500'],
# 	     [file_location+'ENZO/DD0000/DD0000', file_location+'ENZO/DD0100/DD0100'],
# 	     [file_location+'GADGET-3/AGORA_ISO_LOW_ZSolar3/snap_iso_dry_000.hdf5', file_location+'GADGET-3/AGORA_ISO_LOW_ZSolar3/snap_iso_dry_010.hdf5'],
# 	     [file_location+'GASOLINE/LOW_nosf_nofb_gasoline_pfloor_jeanssoft_0myr.00001', file_location+'GASOLINE/LOW_nosf_nofb_gasoline_pfloor_jeanssoft.00335'],
# 	     [file_location+'GEAR/snapshot_0000', file_location+'GEAR/snapshot_0500'],
# 	     [file_location+'GIZMO/snapshot_temp_000', file_location+'GIZMO/snapshot_temp_100'],
# 	     [file_location+'RAMSES/output_00001/info_00001.txt', file_location+'RAMSES/output_00068/info_00068.txt']]
filenames = [[file_location+'ART-I/IC/AGORA_Galaxy_LOW.d', file_location+'ART-I/t0.5Gyr/10MpcBox_csf512_02350.d'],
	     [file_location+'ART-II/SF_FBth_def/OUT/AGORA_LOW_000000.art', file_location+'ART-II/SF_FBth_def/OUT/AGORA_LOW_000315.art'],
	     [file_location+'CHANGA/disklow/disklow.000000', file_location+'CHANGA/disklow/disklow.000500'],
	     [file_location+'ENZO/DD0000/DD0000', file_location+'ENZO/DD0050/DD0050'],
 	     [file_location+'GADGET-3/snap_iso_sf_000.hdf5', file_location+'GADGET-3/snap_iso_sf_010.hdf5'],
	     [file_location+'GASOLINE/LOW_dataset2_timezero.00001', file_location+'GASOLINE/LOW_dataset2.00335'],
  	     [file_location+'GEAR/snapshot_0000', file_location+'GEAR/snapshot_0500'],
	     [file_location+'GIZMO/snapshot_temp_000', file_location+'GIZMO/snapshot_temp_100'],
 	     [file_location+'RAMSES/output_00001/info_00001.txt', file_location+'RAMSES/output_00236/info_00236.txt']]

# codes = ['ART-I']
# filenames = [[file_location+'ART-I/IC/AGORA_Galaxy_LOW.d', file_location+'ART-I/t0.5Gyr/10MpcBox_csf512_02350.d']]
# codes = ['ART-II']
# filenames = [[file_location+'ART-II/SF_FBth_def/OUT/AGORA_LOW_000000.art', file_location+'ART-II/SF_FBth_def/OUT/AGORA_LOW_000315.art']]
# codes = ['CHANGA']
# filenames = [[file_location+'CHANGA/disklow/disklow.000000', file_location+'CHANGA/disklow/disklow.000500']]
# codes = ['ENZO']
# filenames = [[file_location+'ENZO/DD0000/DD0000', file_location+'ENZO/DD0050/DD0050']]
# codes = ['GADGET-3']
# filenames = [[file_location+'GADGET-3/snap_iso_sf_000.hdf5', file_location+'GADGET-3/snap_iso_sf_010.hdf5']]
# codes = ['GASOLINE']
# filenames = [[file_location+'GASOLINE/LOW_dataset2_timezero.00001', file_location+'GASOLINE/LOW_dataset2.00335']]
# codes = ['GEAR']
# filenames = [[file_location+'GEAR/snapshot_0000', file_location+'GEAR/snapshot_0500']]
# codes = ['GIZMO']
# filenames = [[file_location+'GIZMO/snapshot_temp_000', file_location+'GIZMO/snapshot_temp_100']]
# codes = ['RAMSES']
# filenames = [[file_location+'RAMSES/output_00001/info_00001.txt', file_location+'RAMSES/output_00236/info_00236.txt']] 
# codes = ['ENZO', 'GADGET-3']
# filenames = [[file_location+'ENZO/DD0000/DD0000', file_location+'ENZO/DD0050/DD0050'],
#  	     [file_location+'GADGET-3/snap_iso_sf_000.hdf5', file_location+'GADGET-3/snap_iso_sf_010.hdf5']]
gadget_default_unit_base = {'UnitLength_in_cm'         : 3.08568e+21,
			    'UnitMass_in_g'            :   1.989e+43,
			    'UnitVelocity_in_cm_per_s' :      100000}
color_names              = ['red', 'magenta', 'gold', 'lime', 'green', 'cyan', 'blue', 'blueviolet', 'black']
linestyle_names          = ['-', '--', '-.']

draw_density_map       = 1         # 0/1   = OFF/ON
draw_temperature_map   = 1         # 0/1   = OFF/ON
draw_cellsize_map      = 1         # 0/1   = OFF/ON
draw_metal_map         = 1         # 0/1   = OFF/ON
draw_PDF               = 1         # 0/1   = OFF/ON
draw_pos_vel_PDF       = 1         # 0/1/2 = OFF/ON/ON with 1D profile
draw_rad_height_PDF    = 0         # 0/1/2 = OFF/ON/ON with analytic ftn subtracted
draw_metal_PDF         = 1         # 0/1   = OFF/ON
draw_density_DF        = 0         # 0/1   = OFF/ON
draw_radius_DF         = 0         # 0/1   = OFF/ON
draw_height_DF         = 0         # 0/1   = OFF/ON
draw_cut_through       = 0         # 0/1   = OFF/ON
add_nametag            = 1         # 0/1   = OFF/ON
times                  = [0, 500]  # in Myr
figure_width           = 30        # in kpc
n_ref                  = 256       # for SPH codes
over_refine_factor     = 1         # for SPH codes
disk_normal_vector     = [0.0, 0.0, 1.0]

fig_density_map        = [] 
fig_temperature_map    = []
fig_cellsize_map       = []
fig_metal_map          = [] 
fig_PDF                = []
fig_pos_vel_PDF        = []
fig_rad_height_PDF     = []
fig_metal_PDF          = []
grid_density_map       = []
grid_temperature_map   = []
grid_cellsize_map      = []
grid_metal_map         = []
grid_PDF               = []
grid_pos_vel_PDF       = []
grid_rad_height_PDF    = []
grid_metal_PDF         = []
pos_vel_xs             = []
pos_vel_profiles       = []
density_DF_xs          = []
density_DF_profiles    = []
radius_DF_xs           = []
radius_DF_profiles     = []
height_DF_xs           = []
height_DF_profiles     = []
cut_through_zs         = []
cut_through_zvalues    = []
cut_through_xs         = []
cut_through_xvalues    = []

for time in range(len(times)):
	if draw_density_map == 1:
		fig_density_map      += [plt.figure(figsize=(100,20))]
		grid_density_map     += [AxesGrid(fig_density_map[time], (0.01,0.01,0.99,0.99), nrows_ncols = (2, len(codes)), axes_pad = 0.02, add_all = True, share_all = True,
						  label_mode = "1", cbar_mode = "single", cbar_location = "right", cbar_size = "2%", cbar_pad = 0.02)]
	if draw_temperature_map == 1:
		fig_temperature_map  += [plt.figure(figsize=(100,20))]
		grid_temperature_map += [AxesGrid(fig_temperature_map[time], (0.01,0.01,0.99,0.99), nrows_ncols = (2, len(codes)), axes_pad = 0.02, add_all = True, share_all = True,
						  label_mode = "1", cbar_mode = "single", cbar_location = "right", cbar_size = "2%", cbar_pad = 0.02)]
	if draw_cellsize_map == 1:
		fig_cellsize_map     += [plt.figure(figsize=(100,20))]
		grid_cellsize_map    += [AxesGrid(fig_cellsize_map[time], (0.01,0.01,0.99,0.99), nrows_ncols = (2, len(codes)), axes_pad = 0.02, add_all = True, share_all = True,
						  label_mode = "1", cbar_mode = "single", cbar_location = "right", cbar_size = "2%", cbar_pad = 0.02)]
	if draw_metal_map == 1:
		fig_metal_map        += [plt.figure(figsize=(100,20))]
		grid_metal_map       += [AxesGrid(fig_metal_map[time], (0.01,0.01,0.99,0.99), nrows_ncols = (2, len(codes)), axes_pad = 0.02, add_all = True, share_all = True,
						  label_mode = "1", cbar_mode = "single", cbar_location = "right", cbar_size = "2%", cbar_pad = 0.02)]
	if draw_PDF == 1:
		fig_PDF              += [plt.figure(figsize=(50, 80))]
		grid_PDF             += [AxesGrid(fig_PDF[time], (0.01,0.01,0.99,0.99), nrows_ncols = (3, int(math.ceil(len(codes)/3.0))), axes_pad = 0.05, add_all = True, share_all = True,
						  label_mode = "1", cbar_mode = "single", cbar_location = "right", cbar_size = "2%", cbar_pad = 0.05, aspect = False)]
	if draw_pos_vel_PDF >= 1:
		fig_pos_vel_PDF      += [plt.figure(figsize=(50, 80))]
		grid_pos_vel_PDF     += [AxesGrid(fig_pos_vel_PDF[time], (0.01,0.01,0.99,0.99), nrows_ncols = (3, int(math.ceil(len(codes)/3.0))), axes_pad = 0.05, add_all = True, share_all = True,
						  label_mode = "1", cbar_mode = "single", cbar_location = "right", cbar_size = "2%", cbar_pad = 0.05, aspect = False)]
		pos_vel_xs.append([])
		pos_vel_profiles.append([])
	if draw_rad_height_PDF >= 1:
		fig_rad_height_PDF   += [plt.figure(figsize=(50, 80))]
		grid_rad_height_PDF  += [AxesGrid(fig_rad_height_PDF[time], (0.01,0.01,0.99,0.99), nrows_ncols = (3, int(math.ceil(len(codes)/3.0))), axes_pad = 0.05, add_all = True, share_all = True,
						  label_mode = "1", cbar_mode = "single", cbar_location = "right", cbar_size = "2%", cbar_pad = 0.05, aspect = False)]
	if draw_metal_PDF == 1:
		fig_metal_PDF        += [plt.figure(figsize=(50, 80))]
		grid_metal_PDF       += [AxesGrid(fig_metal_PDF[time], (0.01,0.01,0.99,0.99), nrows_ncols = (3, int(math.ceil(len(codes)/3.0))), axes_pad = 0.05, add_all = True, share_all = True,
						  label_mode = "1", cbar_mode = "single", cbar_location = "right", cbar_size = "2%", cbar_pad = 0.05, aspect = False)]
	if draw_density_DF == 1:
		density_DF_xs.append([])
		density_DF_profiles.append([])
	if draw_radius_DF == 1:
		radius_DF_xs.append([])
		radius_DF_profiles.append([])
	if draw_height_DF == 1:
		height_DF_xs.append([])
		height_DF_profiles.append([])
	if draw_cut_through == 1:
		cut_through_zs.append([])
		cut_through_zvalues.append([])
		cut_through_xs.append([])
		cut_through_xvalues.append([])

for time in range(len(times)):
	for code in range(len(codes)):
		if filenames[code] == []:
			continue

		####################################
		#        PRE-ANALYSIS STEPS        #
		####################################

		# LOAD DATASETS
		if codes[code] == 'ART-I': # ART frontend doesn't find the accompanying files, so we specify them; see http://yt-project.org/docs/dev/examining/loading_data.html#art-data
			dirnames = filenames[code][time][:filenames[code][time].rfind('/')+1]
			if time == 0:
				timestamp = ''
			else:
				timestamp = filenames[code][time][filenames[code][time].rfind('_'):filenames[code][time].rfind('.')] 
			pf = load(filenames[code][time], file_particle_header=dirnames+'PMcrd'+timestamp+'.DAT', file_particle_data=dirnames+'PMcrs0'+timestamp+'.DAT', file_particle_stars=dirnames+'stars'+timestamp+'.dat')
	        elif codes[code] == 'CHANGA' or codes[code] == 'GASOLINE': # For TIPSY frontend, always make sure to place your parameter file in the same directory as your datasets
			pf = load(filenames[code][time], n_ref=n_ref, over_refine_factor=over_refine_factor)
	        elif codes[code] == 'GADGET-3': # For GADGET-3 2nd dataset, there somehow exist particles very far from the center; so we use [-2000, 2000] for a bounding_box
			pf = load(filenames[code][time], unit_base = gadget_default_unit_base, bounding_box=[[-2000.0, 2000.0], [-2000.0, 2000.0], [-2000.0, 2000.0]], n_ref=n_ref, over_refine_factor=over_refine_factor)
                elif codes[code] == 'GEAR': 
			from yt.frontends.gadget.definitions import gadget_header_specs
			from yt.frontends.gadget.definitions import gadget_ptype_specs
			from yt.frontends.gadget.definitions import gadget_field_specs
			gadget_header_specs["chemistry"] = (('h1',  4, 'c'),('h2',  4, 'c'),('empty',  256, 'c'),)
			gear_ptype_specs = ("Gas", "Stars", "Halo", "Disk", "Bulge", "Bndry")
			# For GEAR 1st dataset (Grackle+noSF)
			# pf = GadgetDataset(filenames[code][time], unit_base = gadget_default_unit_base, bounding_box=[[-1000.0, 1000.0], [-1000.0, 1000.0], [-1000.0, 1000.0]], header_spec="default+chemistry", 
			#                    ptype_spec=gear_ptype_specs, n_ref=n_ref, over_refine_factor=over_refine_factor)

			# For GEAR 2nd dataset (Grackle+SF+ThermalFbck); here "Metals" is acutally 10-species field; check how metallicity field is added below. 
			agora_gear  = ( "Coordinates",       
					"Velocities",
					"ParticleIDs",
					"Mass",
					("InternalEnergy", "Gas"),
					("Density", "Gas"),
					("SmoothingLength", "Gas"),
					("Metals", "Gas"),
					("StellarFormationTime", "Stars"),
					("StellarInitMass", "Stars"),					
					("StellarIDs", "Stars"),
					("StellarDensity", "Stars"),
					("StellarSmoothingLength", "Stars"),
					("StellarMetals", "Stars"),
					("Opt1", "Stars"),
					("Opt2", "Stars"),
					)
			gadget_field_specs["agora_gear"] = agora_gear
			pf = GadgetDataset(filenames[code][time], unit_base = gadget_default_unit_base, bounding_box=[[-1000.0, 1000.0], [-1000.0, 1000.0], [-1000.0, 1000.0]], header_spec="default+chemistry", 
					   ptype_spec=gear_ptype_specs, field_spec="agora_gear", n_ref=n_ref, over_refine_factor=over_refine_factor)
	        elif codes[code] == 'GIZMO':
			from yt.frontends.gadget.definitions import gadget_field_specs
			# For GIZMO 1st dataset (Grackle+noSF)
			# agora_gizmo = ( "Coordinates",      
			# 		"Velocities",
			# 		"ParticleIDs",
			# 		"Mass",
			# 		("Temperature", "Gas"),
			# 		("Density", "Gas"),
			# 		("Electron_Number_Density", "Gas"),
			# 		("HI_NumberDensity", "Gas"),
			# 		("SmoothingLength", "Gas"),
			# 		)

			# For GIZMO 2nd dataset (Grackle+SF+ThermalFbck)
			agora_gizmo = ( "Coordinates",       
					"Velocities",
					"ParticleIDs",
					"Mass",
					("Temperature", "Gas"),
					("Density", "Gas"),
					("ElectronAbundance", "Gas"),
					("NeutralHydrogenAbundance", "Gas"),
					("SmoothingLength", "Gas"),
					("StarFormationRate", "Gas"),
					("StellarFormationTime", "Stars"),
					("Metallicity", ("Gas", "Stars")),
					("DelayTime", "Stars"),
					("StellarInitMass", "Stars"),					
					)
			gadget_field_specs["agora_gizmo"] = agora_gizmo
			pf = load(filenames[code][time], unit_base = gadget_default_unit_base, bounding_box=[[-1000.0, 1000.0], [-1000.0, 1000.0], [-1000.0,1000.0]], field_spec="agora_gizmo", 
				  n_ref=n_ref, over_refine_factor=over_refine_factor)
		else:
			pf = load(filenames[code][time])

		# GAS PARTICLE FILEDS FOR SPH CODES
		if codes[code] == 'CHANGA' or codes[code] == 'GASOLINE':
			PartType_Gas_to_use = "Gas"
			MassType_to_use = "Mass"
			MetallicityType_to_use = "Metals"
		elif codes[code] == 'GIZMO' or codes[code] == 'GEAR':
			PartType_Gas_to_use = "Gas"
			MassType_to_use = "Mass"
			MetallicityType_to_use = "Metallicity"
		elif codes[code] == "GADGET-3":
			PartType_Gas_to_use = "PartType0"				
			MassType_to_use = "Masses"
			MetallicityType_to_use = "Metallicity"

		# AXIS SWAP FOR PLOT COLLECTION
		pf.coordinates.x_axis[1] = 0
		pf.coordinates.y_axis[1] = 2
		pf.coordinates.x_axis['y'] = 0
		pf.coordinates.y_axis['y'] = 2

		# ADDITIONAL FIELDS I
		def _density_squared(field, data):  
			return data[("gas", "density")]**2
		pf.add_field(("gas", "density_squared"), function=_density_squared, units="g**2/cm**6")
		def _CellSizepc(field,data): 
			return (data[("index", "cell_volume")].in_units('pc**3'))**(1/3.)
		pf.add_field(("index", "cell_size"), function=_CellSizepc, units='pc', display_name="$\Delta$ x", take_log=True )
		def _Inv2CellVolumeCode(field,data): 
			return data[("index", "cell_volume")]**-2
		pf.add_field(("index", "cell_volume_inv2"), function=_Inv2CellVolumeCode, units='code_length**(-6)', display_name="Inv2CellVolumeCode", take_log=True)	
		
		# ADDITIONAL FIELDS II: TEMPERATURE
                if codes[code] == 'GEAR' or codes[code] == 'GADGET-3' or codes[code] == 'RAMSES': 
			# From grackle/src/python/utilities/convenience.py: Calculate a tabulated approximation to mean molecular weight (valid for data that used Grackle 2.0 or below)
			def calc_mu_table_local(temperature):
				tt = np.array([1.0e+01, 1.0e+02, 1.0e+03, 1.0e+04, 1.3e+04, 2.1e+04, 3.4e+04, 6.3e+04, 1.0e+05, 1.0e+09])
				mt = np.array([1.18701555, 1.15484424, 1.09603514, 0.9981496, 0.96346395, 0.65175895, 0.6142901, 0.6056833, 0.5897776, 0.58822635])
				logttt= np.log(temperature)
				logmu = np.interp(logttt,np.log(tt),np.log(mt)) # linear interpolation in log-log space
				return np.exp(logmu)
			temperature_values = []
			mu_values = []
			T_over_mu_values = []
			current_temperature = 1e1
			final_temperature = 1e7
			dlogT = 0.1
			while current_temperature < final_temperature:
				temperature_values.append(current_temperature)
				current_mu = calc_mu_table_local(current_temperature)
				mu_values.append(current_mu)
				T_over_mu_values.append(current_temperature/current_mu)
				current_temperature = np.exp(np.log(current_temperature)+dlogT)
			def convert_T_over_mu_to_T(T_over_mu):
				logT_over_mu = np.log(T_over_mu)
				logT = np.interp(logT_over_mu, np.log(T_over_mu_values), np.log(temperature_values)) # linear interpolation in log-log space
				return np.exp(logT)
			if codes[code] == 'GEAR' or codes[code] == 'GADGET-3':
				def _Temperature_3(field, data):  
					gamma = 5.0/3.0
					T_over_mu = (data[PartType_Gas_to_use, "InternalEnergy"] * (gamma-1) * constants.mass_hydrogen_cgs / constants.boltzmann_constant_cgs).in_units('K').d # T/mu
					return YTArray(convert_T_over_mu_to_T(T_over_mu), 'K') # now T
				pf.add_field((PartType_Gas_to_use, "Temperature"), function=_Temperature_3, particle_type=True, force_override=True, units="K")
			elif codes[code] == 'RAMSES':
				# The pressure field includes the artificial pressure support term, so one needs to be careful (compare with the exsiting yt/frontends/ramses/fields.py)
				def _temperature_3(field, data):  
					T_J = 1800.0  # in K
					n_J = 8.0     # in H/cc
					gamma_0 = 2.0 
					x_H = 0.76    
					mH = 1.66e-24      # from pymses/utils/constants/__init__.py  (vs. in yt, mass_hydrogen_cgs = 1.007947*amu_cgs = 1.007947*1.660538921e-24 = 1.6737352e-24)
					kB = 1.3806504e-16 # from pymses/utils/constants/__init__.py  (vs. in yt, boltzmann_constant_cgs = 1.3806488e-16)
					if time != 0:
						T_over_mu = data["gas", "pressure"].d/data["gas", "density"].d * mH / kB - T_J * (data["gas", "density"].d * x_H / mH / n_J)**(gamma_0 - 1.0) # T/mu = T2 in Ramses
					else:
						T_over_mu = data["gas", "pressure"].d/data["gas", "density"].d * mH / kB # IC: no pressure support
					return YTArray(convert_T_over_mu_to_T(T_over_mu), 'K') # now T
				pf.add_field(("gas", "temperature"), function=_temperature_3, force_override=True, units="K")

		# ADDITIONAL FIELDS III: METALLICITY (IN MASS FRACTION, NOT IN ZSUN)
                if codes[code] == 'ART-I' or codes[code] == 'ART-II': # metallicity field in ART-I has a different meaning (see yt/frontends/art/fields.py), and metallicity field in ART-II is missing
			def _metallicity_2(field, data):  
				return data["gas", "metal_ii_density"] / data["gas", "density"]
			pf.add_field(("gas", "metallicity"), function=_metallicity_2, force_override=True, display_name="Metallicity", take_log=True, units="") 
                elif codes[code] == 'ENZO': # metallicity field in ENZO is in Zsun
			def _metallicity_2(field, data):  
				return data["gas", "metal_density"] / data["gas", "density"]
			pf.add_field(("gas", "metallicity"), function=_metallicity_2, force_override=True, display_name="Metallicity", take_log=True, units="") 
                elif codes[code] == 'GEAR': # "Metals" in GEAR is 10-species field (last one being the total metal fraction), so ("Metals", 10) needs to be added to _vector_fields in frontends/gadget/io.py 
 		 	def _metallicity_2(field, data):  
 		  		if len(data[PartType_Gas_to_use, "Metals"].shape) == 1:
 		  			return data[PartType_Gas_to_use, "Metals"]
 		  		else:
 		  			gear_total_metal_fraction_index = data[PartType_Gas_to_use, "Metals"].shape[-1] - 1 # normally = 9
 		  			return data[PartType_Gas_to_use, "Metals"][:,gear_total_metal_fraction_index].in_units("") # in_units("") turned out to be crucial!  
			# We are creating ("Gas", "Metallicity") here, different from ("Gas", "metallicity") which is auto-generated by yt but doesn't work properly
 		  	pf.add_field((PartType_Gas_to_use, MetallicityType_to_use), function=_metallicity_2, display_name="Metallicity", particle_type=True, take_log=True, units="")
 		  	# Also creating smoothed field following an example in yt-project.org/docs/dev/cookbook/calculating_information.html; use hardcoded num_neighbors as in frontend/gadget/fields.py
 		  	fn = add_volume_weighted_smoothed_field(PartType_Gas_to_use, "Coordinates", MassType_to_use, "SmoothingLength", "Density", MetallicityType_to_use, pf.field_info, nneighbors=64)
 		  	# Alias doesn't work -- e.g. pf.field_info.alias(("gas", "metallicity"), fn[0]) -- probably because pf=GadgetDataset(), not load(); so I add and replace existing ("gas", "metallicity")
			def _metallicity_3(field, data):  
				return data["deposit", PartType_Gas_to_use+"_smoothed_"+MetallicityType_to_use]
 		  	pf.add_field(("gas", "metallicity"), function=_metallicity_3, force_override=True, display_name="Metallicity", particle_type=False, take_log=True, units="")

		# ADDITIONAL FIELDS IV: FAKE PARTICLE FIELDS, CYLINDRICAL COORDINATES, etc.
		def rho_agora_disk(r, z):
			r_d = YTArray(3.432, 'kpc')
			z_d = 0.1*r_d
			M_d = YTArray(4.297e10, 'Msun')
			f_gas = 0.2
			r_0 = (f_gas*M_d / (4.0*np.pi*r_d**2*z_d)).in_units('g/cm**3')
			return r_0*numpy.exp(-r/r_d)*numpy.exp(-z/z_d)

                if codes[code] == "ART-I" or codes[code] == "ART-II" or codes[code] == "ENZO"  or codes[code] == "RAMSES":
			def _cylindrical_z_abs(field, data):
				return numpy.abs(data[("index", "cylindrical_z")])
			pf.add_field(("index", "cylindrical_z_abs"), function=_cylindrical_z_abs, take_log=False, particle_type=False, units="cm") 
			def _density_minus_analytic(field, data):	
				return data[("gas", "density")] - rho_agora_disk(data[("index", "cylindrical_r")], data[("index", "cylindrical_z_abs")])
			pf.add_field(("gas", "density_minus_analytic"), function=_density_minus_analytic, take_log=False, particle_type=False, display_name="density residual abs", units="g/cm**3") 
		else:
			def _particle_position_cylindrical_z_abs(field, data):
				return numpy.abs(data[(PartType_Gas_to_use, "particle_position_cylindrical_z")])
			pf.add_field((PartType_Gas_to_use, "particle_position_cylindrical_z_abs"), function=_particle_position_cylindrical_z_abs, take_log=False, particle_type=True, units="cm") 
			# particle_type=False doesn't make sense, but is critical for PhasePlot/ProfilePlot to work
			# requires a change in data_objects/data_container.py: remove raise YTFieldTypeNotFound(ftype)
			def _Density_2(field, data):
				return data[(PartType_Gas_to_use, "Density")].in_units('g/cm**3')
			pf.add_field((PartType_Gas_to_use, "Density_2"), function=_Density_2, take_log=True, particle_type=False, display_name="Density", units="g/cm**3") 
			def _Temperature_2(field, data):
				return data[(PartType_Gas_to_use, "Temperature")].in_units('K')
			pf.add_field((PartType_Gas_to_use, "Temperature_2"), function=_Temperature_2, take_log=True, particle_type=False, display_name="Temperature", units="K") 
			def _Mass_2(field, data):
				return data[(PartType_Gas_to_use, MassType_to_use)].in_units('Msun')
			pf.add_field((PartType_Gas_to_use, "Mass_2"), function=_Mass_2, take_log=True, particle_type=False, display_name="Mass", units="Msun")				
			def _Metallicity_2(field, data):
			 	return data[(PartType_Gas_to_use, MetallicityType_to_use)]
			pf.add_field((PartType_Gas_to_use, "Metallicity_2"), function=_Metallicity_2, take_log=True, particle_type=False, display_name="Metallicity", units="")				
			def _Density_2_minus_analytic(field, data):	
				return data[(PartType_Gas_to_use, "density")] - rho_agora_disk(data[(PartType_Gas_to_use, "particle_position_cylindrical_radius")], data[(PartType_Gas_to_use, "particle_position_cylindrical_z_abs")])
			pf.add_field((PartType_Gas_to_use, "Density_2_minus_analytic"), function=_Density_2_minus_analytic, take_log=False, particle_type=False, display_name="Density Residual", units="g/cm**3") 

		####################################
		#      MAIN ANALYSIS ROUTINES      #
		####################################

		# FIND CENTER AND PROJ_REGION
		v, cen = pf.h.find_max(("gas", "density")) # find the center to keep the galaxy at the center of all the images.
		sp = pf.sphere(cen, (figure_width, "kpc"))
		center = sp.quantities.center_of_mass(use_gas=True, use_particles=True).in_units("kpc")
                if codes[code] == "ART-I" or codes[code] == "ART-II" or codes[code] == "ENZO"  or codes[code] == "RAMSES":
			proj_region = pf.box(center - YTArray([figure_width, figure_width, figure_width], 'kpc'),
		 			     center + YTArray([figure_width, figure_width, figure_width], 'kpc')) # projected images made using a (2*figure_width)^3 box for AMR codes
		else:
		 	proj_region = pf.all_data()
			
		# DENSITY MAPS
		if draw_density_map == 1: 
			for ax in range(1, 3):  
				p = ProjectionPlot(pf, ax, ("gas", "density"), center = center, data_source=proj_region, width = (figure_width, 'kpc'), weight_field = None, fontsize=9)
				p.set_zlim(("gas", "density"), 1e-6, 1e-1)			
				plot = p.plots[("gas", "density")]

				plot.figure = fig_density_map[time]
				plot.axes = grid_density_map[time][(ax-1)*len(codes)+code].axes
				if code == 0: plot.cax = grid_density_map[time].cbar_axes[0]
				p._setup_plots()

			if add_nametag == 1:
				at = AnchoredText("%s" % codes[code], loc=2, prop=dict(size=6), frameon=True)
				grid_density_map[time][code].axes.add_artist(at)

		# TEMPERATURE MAPS
		if draw_temperature_map == 1:
			for ax in range(1, 3):  
				p2 = ProjectionPlot(pf, ax, ("gas", "temperature"), center = center, data_source=proj_region, width = (figure_width, 'kpc'), weight_field = ("gas", "density_squared"), fontsize=9)
				p2.set_zlim(("gas", "temperature"), 1e1, 1e6)
				plot2 = p2.plots[("gas", "temperature")]
				
				plot2.figure = fig_temperature_map[time]
				plot2.axes = grid_temperature_map[time][(ax-1)*len(codes)+code].axes
				if code == 0: plot2.cax = grid_temperature_map[time].cbar_axes[0]
				p2._setup_plots()

			if add_nametag == 1:
				at = AnchoredText("%s" % codes[code], loc=2, prop=dict(size=6), frameon=True)
				grid_temperature_map[time][code].axes.add_artist(at)

		# CELL-SIZE MAPS
		if draw_cellsize_map == 1:
			for ax in range(1, 3):  
				p25 = ProjectionPlot(pf, ax, ("index", "cell_size"), center = center, data_source=proj_region, width = (figure_width, 'kpc'), weight_field = ("index", "cell_volume_inv2"), fontsize=9)
				p25.set_zlim(("index", "cell_size"), 50, 500)
				plot25 = p25.plots[("index", "cell_size")]
				
				plot25.figure = fig_cellsize_map[time]
				plot25.axes = grid_cellsize_map[time][(ax-1)*len(codes)+code].axes
				if code == 0: plot25.cax = grid_cellsize_map[time].cbar_axes[0]
				p25._setup_plots()

			if add_nametag == 1:
				at = AnchoredText("%s" % codes[code], loc=2, prop=dict(size=6), frameon=True)
				grid_cellsize_map[time][code].axes.add_artist(at)

		# METAL MAPS
		if draw_metal_map == 1: 
			for ax in range(1, 3):  
				p = ProjectionPlot(pf, ax, ("gas", "metallicity"), center = center, data_source=proj_region, width = (figure_width, 'kpc'), weight_field = ("gas", "density_squared"), fontsize=9)
				p.set_zlim(("gas", "metallicity"), 0.005, 0.05)			
				plot = p.plots[("gas", "metallicity")]

				plot.figure = fig_metal_map[time]
				plot.axes = grid_metal_map[time][(ax-1)*len(codes)+code].axes
				if code == 0: plot.cax = grid_metal_map[time].cbar_axes[0]
				p._setup_plots()

			if add_nametag == 1:
				at = AnchoredText("%s" % codes[code], loc=2, prop=dict(size=6), frameon=True)
				grid_metal_map[time][code].axes.add_artist(at)

		# DENSITY-TEMPERATURE PDF
		if draw_PDF == 1:
			sp = pf.sphere(center, (0.5*figure_width, "kpc"))
			if codes[code] == "ART-I" or codes[code] == "ART-II" or codes[code] == "ENZO"  or codes[code] == "RAMSES":
				p3 = PhasePlot(sp, ("gas", "density"), ("gas", "temperature"), ("gas", "cell_mass"), weight_field=None, fontsize=12, x_bins=500, y_bins=500)
				p3.set_unit("cell_mass", 'Msun')
				p3.set_zlim(("gas", "cell_mass"), 1e3, 1e8)
				p3.set_colorbar_label(("gas", "cell_mass"), "Mass ($\mathrm{M}_{\odot}$)")
				plot3 = p3.plots[("gas", "cell_mass")]
			else:
				# Because ParticlePhasePlot doesn't yet work for a log-log PDF for some reason, I will do the following trick.  
				p3 = PhasePlot(sp, (PartType_Gas_to_use, "Density_2"), (PartType_Gas_to_use, "Temperature_2"), (PartType_Gas_to_use, "Mass_2"), weight_field=None, fontsize=12, x_bins=500, y_bins=500)
				p3.set_zlim((PartType_Gas_to_use, "Mass_2"), 1e3, 1e8)
				plot3 = p3.plots[(PartType_Gas_to_use, "Mass_2")]

			p3.set_xlim(1e-29, 1e-21)
			p3.set_ylim(10, 1e7)

			plot3.figure = fig_PDF[time]
			plot3.axes = grid_PDF[time][code].axes
			if code == 0: plot3.cax = grid_PDF[time].cbar_axes[0]
			p3._setup_plots()

			if add_nametag == 1:
				at = AnchoredText("%s" % codes[code], loc=3, prop=dict(size=10), frameon=True)
				grid_PDF[time][code].axes.add_artist(at)

		# POSITION-VELOCITY PDF
		if draw_pos_vel_PDF >= 1:
			sp = pf.sphere(center, (0.5*figure_width, "kpc"))
			sp.set_field_parameter("normal", disk_normal_vector) 
			if codes[code] == "ART-I" or codes[code] == "ART-II" or codes[code] == "ENZO"  or codes[code] == "RAMSES":
				sp_dense = sp.cut_region(["obj['gas', 'density'].in_units('g/cm**3') > 1.e-25"]) # For AMR codes, consider only the cells that are dense enough
				pf.field_info[("index", "cylindrical_r")].take_log = False
				pf.field_info[("gas", "cylindrical_tangential_velocity")].take_log = False
				p4 = PhasePlot(sp_dense, ("index", "cylindrical_r"), ("gas", "cylindrical_tangential_velocity"), ("gas", "cell_mass"), weight_field=None, fontsize=12, x_bins=300, y_bins=300)
				p4.set_unit("cylindrical_r", 'kpc')
				p4.set_unit("cylindrical_tangential_velocity", 'km/s')
				p4.set_unit("cell_mass", 'Msun')
				p4.set_zlim(("gas", "cell_mass"), 1e3, 1e8)
				p4.set_colorbar_label(("gas", "cell_mass"), "Mass ($\mathrm{M}_{\odot}$)")
				plot4 = p4.plots[("gas", "cell_mass")]
			else:
				pf.field_info[(PartType_Gas_to_use, "particle_position_cylindrical_radius")].take_log = False
				pf.field_info[(PartType_Gas_to_use, "particle_velocity_cylindrical_theta")].take_log = False
				p4 = ParticlePhasePlot(sp, (PartType_Gas_to_use, "particle_position_cylindrical_radius"), (PartType_Gas_to_use, "particle_velocity_cylindrical_theta"), \
							       (PartType_Gas_to_use, MassType_to_use), weight_field=None, fontsize=12, x_bins=300, y_bins=300)
				p4.set_unit("particle_position_cylindrical_radius", 'kpc')
				p4.set_unit("particle_velocity_cylindrical_theta", 'km/s')
				p4.set_unit(MassType_to_use, 'Msun')
				p4.set_zlim((PartType_Gas_to_use, MassType_to_use), 1e3, 1e8)
				p4.set_colorbar_label((PartType_Gas_to_use, MassType_to_use), "Mass ($\mathrm{M}_{\odot}$)")
				plot4 = p4.plots[(PartType_Gas_to_use, MassType_to_use)]

			p4.set_xlabel("Cylindrical Radius (kpc)")
			p4.set_ylabel("Rotational Velocity (km/s)")
 			p4.set_xlim(0, 14)
			p4.set_ylim(-100, 350)

			plot4.figure = fig_pos_vel_PDF[time]
			plot4.axes = grid_pos_vel_PDF[time][code].axes
			if code == 0: plot4.cax = grid_pos_vel_PDF[time].cbar_axes[0]
			p4._setup_plots()

			# Add 1D profile line if requested
			if draw_pos_vel_PDF == 2 and time != 0:
				if codes[code] == "ART-I" or codes[code] == "ART-II" or codes[code] == "ENZO"  or codes[code] == "RAMSES":
					p5 = ProfilePlot(sp_dense, ("index", "cylindrical_r"),  ("gas", "cylindrical_tangential_velocity"), \
								 weight_field=("gas", "cell_mass"), n_bins=50, x_log=False)
					p5.set_log("cylindrical_tangential_velocity", False)
					p5.set_log("cylindrical_r", False)
					p5.set_unit("cylindrical_r", 'kpc')
					p5.set_xlim(1e-3, 14)
					p5.set_ylim("cylindrical_tangential_velocity", -100, 350)
					line = ln.Line2D(p5.profiles[0].x.in_units('kpc'), p5.profiles[0]["cylindrical_tangential_velocity"].in_units('km/s'), linestyle="-", linewidth=2, color='k', alpha=0.7)
					pos_vel_xs[time].append(p5.profiles[0].x.in_units('kpc').d)
					pos_vel_profiles[time].append(p5.profiles[0]["cylindrical_tangential_velocity"].in_units('km/s').d)
				else:
					p5 = ProfilePlot(sp, (PartType_Gas_to_use, "particle_position_cylindrical_radius"), (PartType_Gas_to_use, "particle_velocity_cylindrical_theta"), \
								 weight_field=(PartType_Gas_to_use, MassType_to_use), n_bins=50, x_log=False)
					p5.set_log("particle_velocity_cylindrical_theta", False)
					p5.set_log("particle_position_cylindrical_radius", False)
					p5.set_unit("particle_position_cylindrical_radius", 'kpc')
					p5.set_xlim(1e-3, 14)
					p5.set_ylim("particle_velocity_cylindrical_theta", -100, 350)
					line = ln.Line2D(p5.profiles[0].x.in_units('kpc'), p5.profiles[0]["particle_velocity_cylindrical_theta"].in_units('km/s'), linestyle="-", linewidth=2, color='k', alpha=0.7)
					pos_vel_xs[time].append(p5.profiles[0].x.in_units('kpc').d)
					pos_vel_profiles[time].append(p5.profiles[0]["particle_velocity_cylindrical_theta"].in_units('km/s').d)
				grid_pos_vel_PDF[time][code].axes.add_line(line) 

			if add_nametag == 1:
				at = AnchoredText("%s" % codes[code], loc=4, prop=dict(size=10), frameon=True)
				grid_pos_vel_PDF[time][code].axes.add_artist(at)

		# RADIUS-HEIGHT PDF
		if draw_rad_height_PDF >= 1:
			sp = pf.sphere(center, (0.5*figure_width, "kpc"))
			if codes[code] == "ART-I" or codes[code] == "ART-II" or codes[code] == "ENZO"  or codes[code] == "RAMSES":
				if draw_rad_height_PDF == 1:
					p55 = PhasePlot(sp, ("index", "cylindrical_r"), ("index", "cylindrical_z_abs"), ("gas", "density"), weight_field=("gas", "cell_mass"), fontsize=12, x_bins=200, y_bins=200)
					p55.set_zlim(("gas", "density"), 1e-26, 1e-21)
					p55.set_log("cylindrical_r", False)
					p55.set_log("cylindrical_z_abs", False)
					p55.set_unit("cylindrical_r", 'kpc')
					p55.set_unit("cylindrical_z_abs", 'kpc')
					plot55 = p55.plots[("gas", "density")]
				elif draw_rad_height_PDF == 2:
					p55 = PhasePlot(sp, ("index", "cylindrical_r"), ("index", "cylindrical_z_abs"), ("gas", "density_minus_analytic"), weight_field=("gas", "cell_mass"), fontsize=12, x_bins=200, y_bins=200)
					p55.set_zlim(("gas", "density_minus_analytic"), -1e-24, 1e-24)
					p55.set_log("cylindrical_r", False)
					p55.set_log("cylindrical_z_abs", False)
					p55.set_unit("cylindrical_r", 'kpc')
					p55.set_unit("cylindrical_z_abs", 'kpc')
					plot55 = p55.plots[("gas", "density_minus_analytic")]
			else:
				# Because ParticlePhasePlot doesn't yet work for a log-log PDF for some reason, I will do the following trick.  
				if draw_rad_height_PDF == 1:
					p55 = PhasePlot(sp, (PartType_Gas_to_use, "particle_position_cylindrical_radius"), (PartType_Gas_to_use, "particle_position_cylindrical_z_abs"), (PartType_Gas_to_use, "Density_2"), weight_field=(PartType_Gas_to_use, "Mass_2"), fontsize=12, x_bins=200, y_bins=200)
					p55.set_zlim((PartType_Gas_to_use, "Density_2"), 1e-26, 1e-21)
					p55.set_log("particle_position_cylindrical_radius", False)
					p55.set_log("particle_position_cylindrical_z_abs", False)
					p55.set_unit("particle_position_cylindrical_radius", 'kpc')
					p55.set_unit("particle_position_cylindrical_z_abs", 'kpc')
					plot55 = p55.plots[(PartType_Gas_to_use, "Density_2")]
				elif draw_rad_height_PDF == 2:
					p55 = PhasePlot(sp, (PartType_Gas_to_use, "particle_position_cylindrical_radius"), (PartType_Gas_to_use, "particle_position_cylindrical_z_abs"), (PartType_Gas_to_use, "Density_2_minus_analytic"), weight_field=(PartType_Gas_to_use, "Mass_2"), fontsize=12, x_bins=200, y_bins=200)
					p55.set_zlim((PartType_Gas_to_use, "Density_2_minus_analytic"), -1e-24, 1e-24)
					p55.set_log("particle_position_cylindrical_radius", False)
					p55.set_log("particle_position_cylindrical_z_abs", False)
					p55.set_unit("particle_position_cylindrical_radius", 'kpc')
					p55.set_unit("particle_position_cylindrical_z_abs", 'kpc')
					plot55 = p55.plots[(PartType_Gas_to_use, "Density_2_minus_analytic")]

			p55.set_xlabel("Cylindrical Radius (kpc)")
			p55.set_ylabel("Vertical Height (kpc)")
			p55.set_xlim(0, 14)
			p55.set_ylim(0, 1.4)

			plot55.figure = fig_rad_height_PDF[time]
			plot55.axes = grid_rad_height_PDF[time][code].axes
			if code == 0: plot55.cax = grid_rad_height_PDF[time].cbar_axes[0]
			p55._setup_plots()

			if add_nametag == 1:
				at = AnchoredText("%s" % codes[code], loc=1, prop=dict(size=10), frameon=True)
				grid_rad_height_PDF[time][code].axes.add_artist(at)

		# DENSITY-TEMPERATURE-METALLICITY PDF
		if draw_metal_PDF == 1:
			sp = pf.sphere(center, (0.5*figure_width, "kpc"))
			if codes[code] == "ART-I" or codes[code] == "ART-II" or codes[code] == "ENZO"  or codes[code] == "RAMSES":
				p3 = PhasePlot(sp, ("gas", "density"), ("gas", "temperature"), ("gas", "metallicity"), weight_field=("gas", "cell_mass"), fontsize=12, x_bins=500, y_bins=500)
				p3.set_zlim(("gas", "metallicity"), 0.005, 0.05)
				p3.set_colorbar_label(("gas", "metallicity"), "Metallicity (Mass-weighted average of mass fraction)")
				plot3 = p3.plots[("gas", "metallicity")]
			else:
				# Because ParticlePhasePlot doesn't yet work for a log-log PDF for some reason, I will do the following trick.  
				p3 = PhasePlot(sp, (PartType_Gas_to_use, "Density_2"), (PartType_Gas_to_use, "Temperature_2"), (PartType_Gas_to_use, "Metallicity_2"), weight_field=(PartType_Gas_to_use, "Mass_2"), fontsize=12, x_bins=500, y_bins=500)
				p3.set_zlim((PartType_Gas_to_use, "Metallicity_2"), 0.005, 0.05)
				plot3 = p3.plots[(PartType_Gas_to_use, "Metallicity_2")]

			p3.set_xlim(1e-29, 1e-21)
			p3.set_ylim(10, 1e7)

			plot3.figure = fig_metal_PDF[time]
			plot3.axes = grid_metal_PDF[time][code].axes
			if code == 0: plot3.cax = grid_metal_PDF[time].cbar_axes[0]
			p3._setup_plots()

			if add_nametag == 1:
				at = AnchoredText("%s" % codes[code], loc=3, prop=dict(size=10), frameon=True)
				grid_metal_PDF[time][code].axes.add_artist(at)

		# DENSITY DF (DISTRIBUTION FUNCTION)
		if draw_density_DF == 1:
			sp = pf.sphere(center, (0.5*figure_width, "kpc"))
			if codes[code] == "ART-I" or codes[code] == "ART-II" or codes[code] == "ENZO"  or codes[code] == "RAMSES":
				p6 = ProfilePlot(sp, ("gas", "density"),  ("gas", "cell_mass"), weight_field=None, n_bins=50, x_log=True, accumulation=False)
#				p6 = ProfilePlot(sp, ("gas", "density"),  ("gas", "cell_mass"), weight_field=None, n_bins=50, x_log=True, accumulation=True)
				p6.set_log("cell_mass", True)
				p6.set_xlim(1e-29, 1e-21)
				density_DF_xs[time].append(p6.profiles[0].x.in_units('g/cm**3').d)
				density_DF_profiles[time].append(p6.profiles[0]["cell_mass"].in_units('Msun').d)
			else:
				# Because ParticleProfilePlot doesn't exist, I will do the following trick.  
				p6 = ProfilePlot(sp, (PartType_Gas_to_use, "Density_2"),  (PartType_Gas_to_use, "Mass_2"), weight_field=None, n_bins=50, x_log=True, accumulation=False)
#				p6 = ProfilePlot(sp, (PartType_Gas_to_use, "Density_2"),  (PartType_Gas_to_use, "Mass_2"), weight_field=None, n_bins=50, x_log=True, accumulation=True)
				p6.set_log("Mass_2", True)
				p6.set_xlim(1e-29, 1e-21)
				density_DF_xs[time].append(p6.profiles[0].x.in_units('g/cm**3').d)
				density_DF_profiles[time].append(p6.profiles[0]["Mass_2"].in_units('Msun').d)

		# CYLINDRICAL RADIUS DF + RADIALLY-BINNED GAS SURFACE DENSITY 
		if draw_radius_DF == 1:
			sp = pf.sphere(center, (0.5*figure_width, "kpc"))
			sp.set_field_parameter("normal", disk_normal_vector) 
			if codes[code] == "ART-I" or codes[code] == "ART-II" or codes[code] == "ENZO"  or codes[code] == "RAMSES":
				p7 = ProfilePlot(sp, ("index", "cylindrical_r"),  ("gas", "cell_mass"), weight_field=None, n_bins=50, x_log=False, accumulation=False)
				p7.set_log("cell_mass", True)
				p7.set_log("cylindrical_r", False)
				p7.set_unit("cylindrical_r", 'kpc')
				p7.set_xlim(1e-3, 15) 
				radius_DF_xs[time].append(p7.profiles[0].x.in_units('kpc').d)
				radius_DF_profiles[time].append(p7.profiles[0]["cell_mass"].in_units('Msun').d)
			else:
				p7 = ProfilePlot(sp, (PartType_Gas_to_use, "particle_position_cylindrical_radius"),  (PartType_Gas_to_use, "Mass_2"), weight_field=None, n_bins=50, x_log=False, accumulation=False)
				p7.set_log("Mass_2", True)
				p7.set_log("particle_position_cylindrical_radius", False)
				p7.set_unit("particle_position_cylindrical_radius", 'kpc')
				p7.set_xlim(1e-3, 15) 
				radius_DF_xs[time].append(p7.profiles[0].x.in_units('kpc').d)
				radius_DF_profiles[time].append(p7.profiles[0]["Mass_2"].in_units('Msun').d)

		# VERTICAL HEIGHT DF + VERTICALLY-BINNED GAS SURFACE DENSITY 
		if draw_height_DF == 1:
			sp = pf.sphere(center, (0.5*figure_width, "kpc"))
			sp.set_field_parameter("normal", disk_normal_vector) 
			if codes[code] == "ART-I" or codes[code] == "ART-II" or codes[code] == "ENZO"  or codes[code] == "RAMSES":
				p8 = ProfilePlot(sp, ("index", "cylindrical_z_abs"),  ("gas", "cell_mass"), weight_field=None, n_bins=10, x_log=False, accumulation=False)
				p8.set_log("cell_mass", True)
				p8.set_log("cylindrical_z_abs", False)
				p8.set_unit("cylindrical_z_abs", 'kpc')
				p8.set_xlim(1e-3, 1.4)
				height_DF_xs[time].append(p8.profiles[0].x.in_units('kpc').d)
				height_DF_profiles[time].append(p8.profiles[0]["cell_mass"].in_units('Msun').d)
			else:
				p8 = ProfilePlot(sp, (PartType_Gas_to_use, "particle_position_cylindrical_z_abs"),  (PartType_Gas_to_use, "Mass_2"), weight_field=None, n_bins=10, x_log=False, accumulation=False)
				p8.set_log("Mass_2", True)
				p8.set_log("particle_position_cylindrical_z_abs", False)
				p8.set_unit("particle_position_cylindrical_z_abs", 'kpc')
				p8.set_xlim(1e-3, 1.4)
				height_DF_xs[time].append(p8.profiles[0].x.in_units('kpc').d)
				height_DF_profiles[time].append(p8.profiles[0]["Mass_2"].in_units('Msun').d)

		# DENSITY ALONG THE ORTHO-RAY OBJECT CUTTING THROUGH THE CENTER
		if draw_cut_through == 1:
			ray = pf.ortho_ray(2, (center[0].in_units('code_length'), center[1].in_units('code_length'))) # see: http://yt-project.org/doc/visualizing/manual_plotting.html#line-plots
			srt = np.argsort(ray['z'])
			cut_through_zs[time].append(np.array(ray['z'][srt].in_units('kpc').d - center[2].in_units('kpc').d)) 
			cut_through_zvalues[time].append(np.array(ray[("gas", "density")][srt].in_units('g/cm**3').d))

			ray = pf.ortho_ray(0, (center[1].in_units('code_length'), center[2].in_units('code_length'))) 
			srt = np.argsort(ray['x'])
			cut_through_xs[time].append(np.array(ray['x'][srt].in_units('kpc').d - center[2].in_units('kpc').d)) 
			cut_through_xvalues[time].append(np.array(ray[("gas", "density")][srt].in_units('g/cm**3').d))

	####################################
	#        POST-ANALYSIS STEPS       #
	####################################

	# SAVE FIGURES
	if draw_density_map == 1:
		fig_density_map[time].savefig("Sigma_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)
	if draw_temperature_map == 1:
		fig_temperature_map[time].savefig("Temp_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)
	if draw_cellsize_map == 1:
		fig_cellsize_map[time].savefig("Cell_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)
	if draw_metal_map == 1:
		fig_metal_map[time].savefig("Metal2_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)
	if draw_PDF == 1:
		fig_PDF[time].savefig("PDF_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)
	if draw_pos_vel_PDF >= 1:
		fig_pos_vel_PDF[time].savefig("pos_vel_PDF_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)
		if draw_pos_vel_PDF == 2 and time != 0:
			plt.clf()
			plt.subplot(111, aspect=0.02)
			for code in range(len(codes)):
				lines = plt.plot(pos_vel_xs[time][code], pos_vel_profiles[time][code], color=color_names[code], linestyle=linestyle_names[np.mod(code, len(linestyle_names))])
			plt.xlim(0, 14)
			plt.ylim(-100, 350)
			plt.xlabel("$\mathrm{Cylindrical\ Radius\ (kpc)}$")
			plt.ylabel("$\mathrm{Rotational\ Velocity\ (km/s)}$")
			plt.legend(codes, loc=4, frameon=True)
			leg = plt.gca().get_legend()
			ltext = leg.get_texts()
			plt.setp(ltext, fontsize='small')
			plt.savefig("pos_vel_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)
			plt.clf()
	if draw_rad_height_PDF >= 1:
		fig_rad_height_PDF[time].savefig("rad_height_PDF_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)
	if draw_metal_PDF == 1:
		fig_metal_PDF[time].savefig("metal_PDF_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)
	if draw_density_DF == 1:
		plt.clf()
		plt.subplot(111, aspect=1)
		for code in range(len(codes)):
			lines = plt.plot(density_DF_xs[time][code], density_DF_profiles[time][code], color=color_names[code], linestyle=linestyle_names[np.mod(code, len(linestyle_names))])
		plt.semilogx()
		plt.semilogy()
		plt.xlim(1e-29, 1e-21) 
		plt.ylim(1e4, 1e9) # accumulation=False
		# plt.ylim(1e5, 2e10) # accumulation=True
		plt.xlabel("$\mathrm{Density\ (g/cm^3)}$")
		plt.ylabel("$\mathrm{Mass\ (M_{\odot})}$")
		plt.legend(codes, loc=2, frameon=True)
		leg = plt.gca().get_legend()
		ltext = leg.get_texts()
		plt.setp(ltext, fontsize='xx-small')
		plt.savefig("density_DF_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)
		plt.clf()
	if draw_radius_DF == 1:
		plt.clf()
		plt.subplot(111, aspect=1)
		for code in range(len(codes)):
			lines = plt.plot(radius_DF_xs[time][code], np.add.accumulate(radius_DF_profiles[time][code]), color=color_names[code], linestyle=linestyle_names[np.mod(code, len(linestyle_names))])
		plt.semilogy()
		plt.xlim(0, 14)
		plt.ylim(1e7, 2e10)
		plt.xlabel("$\mathrm{Cylindrical\ Radius\ (kpc)}$")
		plt.ylabel("$\mathrm{Mass\ (M_{\odot})}$")
		plt.legend(codes, loc=4, frameon=True)
		leg = plt.gca().get_legend()
		ltext = leg.get_texts()
		plt.setp(ltext, fontsize='small')
		plt.savefig("radius_DF_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)

		plt.clf()
		plt.subplot(111, aspect=1)
		for code in range(len(codes)):
			surface_density = []
			for radius in range(len(radius_DF_profiles[time][code])):
				if radius == 0:
					surface_area = 4*np.pi*(radius_DF_xs[time][code][radius]*1e3)**2 # in pc^2
				else:
					surface_area = 4*np.pi*((radius_DF_xs[time][code][radius]*1e3)**2 - (radius_DF_xs[time][code][radius-1]*1e3)**2)
				surface_density.append(radius_DF_profiles[time][code][radius] / surface_area)
			lines = plt.plot(radius_DF_xs[time][code], surface_density, color=color_names[code], linestyle=linestyle_names[np.mod(code, len(linestyle_names))])
		plt.semilogy()
		plt.xlim(0, 14)
		plt.ylim(1e-1, 5e3)
		plt.xlabel("$\mathrm{Cylindrical\ Radius\ (kpc)}$")
		plt.ylabel("$\mathrm{Surface\ Density\ (M_{\odot}/pc^2)}$")
		plt.legend(codes, loc=1, frameon=True)
		leg = plt.gca().get_legend()
		ltext = leg.get_texts()
		plt.setp(ltext, fontsize='small')
		plt.savefig("gas_surface_density_radial_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)
		plt.clf()
	if draw_height_DF == 1:
		plt.clf()
		plt.subplot(111, aspect=0.5)
		for code in range(len(codes)):
			lines = plt.plot(height_DF_xs[time][code], np.add.accumulate(height_DF_profiles[time][code]), color=color_names[code], linestyle=linestyle_names[np.mod(code, len(linestyle_names))])
		plt.semilogy()
		plt.xlim(0, 1.4)
		plt.ylim(1e9, 2e10)
		plt.xlabel("$\mathrm{Vertical\ Height\ (kpc)}$")
		plt.ylabel("$\mathrm{Mass\ (M_{\odot})}$")
		plt.legend(codes, loc=4, frameon=True)
		leg = plt.gca().get_legend()
		ltext = leg.get_texts()
		plt.setp(ltext, fontsize='small')
		plt.savefig("height_DF_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)

		plt.clf()
		plt.subplot(111, aspect=0.8)
		for code in range(len(codes)):
			surface_density = []
			for height in range(len(height_DF_profiles[time][code])):
				if height == 0:
					surface_area = height_DF_xs[time][code][height]*1e3 * figure_width*1e3 * 2 # surface_area = 2*d(height)*figure_width in pc^2
				else:
					surface_area = (height_DF_xs[time][code][height] - height_DF_xs[time][code][height-1])*1e3 * figure_width*1e3 * 2
				surface_density.append(height_DF_profiles[time][code][height] / surface_area)
			lines = plt.plot(height_DF_xs[time][code], surface_density, color=color_names[code], linestyle=linestyle_names[np.mod(code, len(linestyle_names))])
		plt.semilogy()
		plt.xlim(0, 1.4)
		plt.ylim(1e-1, 5e3)
		plt.xlabel("$\mathrm{Vertical\ Height\ (kpc)}$")
		plt.ylabel("$\mathrm{Surface\ Density\ (M_{\odot}/pc^2)}$")
		plt.legend(codes, loc=1, frameon=True)
		leg = plt.gca().get_legend()
		ltext = leg.get_texts()
		plt.setp(ltext, fontsize='small')
		plt.savefig("gas_surface_density_vertical_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)
		plt.clf()
	if draw_cut_through == 1:
		plt.clf()
		plt.subplot(111, aspect=0.5)
		for code in range(len(codes)):
			lines = plt.plot(cut_through_zs[time][code], cut_through_zvalues[time][code], color=color_names[code], linestyle=linestyle_names[np.mod(code, len(linestyle_names))])
		plt.semilogy()
		plt.xlim(-1.4, 1.4)
		plt.ylim(1e-26, 1e-22)
		plt.xlabel("$\mathrm{Vertical\ Height\ (kpc)}$")
		plt.ylabel("$\mathrm{Density\ (g/cm^3)}$")
		plt.legend(codes, loc=1, frameon=True)
		leg = plt.gca().get_legend()
		ltext = leg.get_texts()
		plt.setp(ltext, fontsize='x-small')
		z = np.arange(-1.4, 1.4, 0.05)
		plt.plot(z, rho_agora_disk(0, np.abs(z)), linestyle="--", linewidth=2, color='k', alpha=0.7)
		plt.savefig("cut_through_z_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)
		plt.clf()
		plt.subplot(111, aspect=0.5)
		for code in range(len(codes)):
			lines = plt.plot(cut_through_xs[time][code], cut_through_xvalues[time][code], color=color_names[code], linestyle=linestyle_names[np.mod(code, len(linestyle_names))])
		plt.semilogy()
		plt.xlim(-14, 14)
		plt.ylim(1e-26, 1e-22)
		plt.xlabel("$\mathrm{Cylindrical\ Radius\ (kpc)}$")
		plt.ylabel("$\mathrm{Density\ (g/cm^3)}$")
		plt.legend(codes, loc=1, frameon=True)
		leg = plt.gca().get_legend()
		ltext = leg.get_texts()
		plt.setp(ltext, fontsize='x-small')
		x = np.arange(-14, 14, 0.05)
		plt.plot(x, rho_agora_disk(np.abs(x), 0), linestyle="--", linewidth=2, color='k', alpha=0.7)
		plt.savefig("cut_through_x_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)

