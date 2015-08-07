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

codes = ['ART-I', 'ART-II', 'CHANGA', 'ENZO', 'GADGET-3', 'GASOLINE', 'GEAR', 'GIZMO', 'RAMSES']
filenames = [['./ART-I/IC/AGORA_Galaxy_LOW.d', './ART-I/t0.5Gyr/10MpcBox_csf512_02350.d'],
	     ['./ART-II/noSF_aggrRef/OUT/AGORA_LOW_000000.art', './ART-II/noSF_aggrRef/OUT/AGORA_LOW_000084.art'],
	     ['./CHANGA/disklow/disklow.000000', './CHANGA/disklow/disklow.000500'], 
	     ['./ENZO/DD0000/DD0000', './ENZO/DD0100/DD0100'],
 	     ['./GADGET-3/AGORA_ISO_LOW_ZSolar/snap_iso_dry_000.hdf5', './GADGET-3/AGORA_ISO_LOW_ZSolar/snap_iso_dry_010.hdf5'],
 	     ['./GASOLINE/LOW_nosf_nofb_gasoline_pfloor_jeanssoft_0myr.00001', './GASOLINE/LOW_nosf_nofb_gasoline_pfloor_jeanssoft.00335'],
  	     ['./GEAR/snapshot_0000', './GEAR/snapshot_0500'],
 	     ['./GIZMO/snapshot_000', './GIZMO/snapshot_100_hsml'],
 	     ['./RAMSES/output_00001/info_00001.txt', './RAMSES/output_00068/info_00068.txt']]
# codes = ['ART-I']
# filenames = [['./ART-I/IC/AGORA_Galaxy_LOW.d', './ART-I/t0.5Gyr/10MpcBox_csf512_02350.d']]
# codes = ['ART-II']
# filenames = [['./ART-II/noSF_aggrRef/OUT/AGORA_LOW_000000.art', './ART-II/noSF_aggrRef/OUT/AGORA_LOW_000084.art']]
# codes = ['CHANGA']
# filenames = [['./CHANGA/disklow/disklow.000000', './CHANGA/disklow/disklow.000500']]
# codes = ['ENZO']
# filenames = [['./ENZO/DD0000/DD0000', './ENZO/DD0100/DD0100']]
# codes = ['GADGET-3']
# filenames = [['./GADGET-3/AGORA_ISO_LOW_ZSolar/snap_iso_dry_000.hdf5', './GADGET-3/AGORA_ISO_LOW_ZSolar/snap_iso_dry_010.hdf5']]
# codes = ['GASOLINE']
# filenames = [['./GASOLINE/LOW_nosf_nofb_gasoline_pfloor_jeanssoft_0myr.00001', './GASOLINE/LOW_nosf_nofb_gasoline_pfloor_jeanssoft.00335']]
# codes = ['GEAR']
# filenames = [['./GEAR/snapshot_0000', './GEAR/snapshot_0500']]
# codes = ['GIZMO']
# filenames = [['./GIZMO/snapshot_000', './GIZMO/snapshot_100_hsml']]
# codes = ['RAMSES']
# filenames = [['./RAMSES/output_00001/info_00001.txt', './RAMSES/output_00068/info_00068.txt']] 
gadget_default_unit_base = {'UnitLength_in_cm'         : 3.08568e+21,
			    'UnitMass_in_g'            :   1.989e+43,
			    'UnitVelocity_in_cm_per_s' :      100000}

draw_density_map     = 1         # 1/0 = ON/OFF
draw_temperature_map = 1         # 1/0 = ON/OFF
draw_PDF             = 1         # 1/0 = ON/OFF
draw_pos_vel_PDF     = 1         # 1/0 = ON/OFF
add_nametag          = 1         # 1/0 = ON/OFF
times                = [0, 500]  # in Myr
figure_width         = 30        # in kpc
n_ref                = 256       # for SPH codes
over_refine_factor   = 1         # for SPH codes
disk_normal_vector   = [0.0, 0.0, 1.0]

fig_density_map      = [] 
fig_temperature_map  = []
fig_PDF              = []
fig_pos_vel_PDF      = []
grid_density_map     = []
grid_temperature_map = []
grid_PDF             = []
grid_pos_vel_PDF     = []

for time in range(len(times)):
	fig_density_map      += [plt.figure(figsize=(100,20))]
	fig_temperature_map  += [plt.figure(figsize=(100,20))]
	fig_PDF              += [plt.figure(figsize=(50, 80))]
	fig_pos_vel_PDF      += [plt.figure(figsize=(50, 80))]
	grid_density_map     += [AxesGrid(fig_density_map[time], (0.01,0.01,0.99,0.99), nrows_ncols = (2, len(codes)), axes_pad = 0.02, add_all = True, share_all = True,
					  label_mode = "1", cbar_mode = "single", cbar_location = "right", cbar_size = "2%", cbar_pad = 0.02)]
	grid_temperature_map += [AxesGrid(fig_temperature_map[time], (0.01,0.01,0.99,0.99), nrows_ncols = (2, len(codes)), axes_pad = 0.02, add_all = True, share_all = True,
					  label_mode = "1", cbar_mode = "single", cbar_location = "right", cbar_size = "2%", cbar_pad = 0.02)]
	grid_PDF             += [AxesGrid(fig_PDF[time], (0.01,0.01,0.99,0.99), nrows_ncols = (3, int(math.ceil(len(codes)/3.0))), axes_pad = 0.05, add_all = True, share_all = True,
					  label_mode = "1", cbar_mode = "single", cbar_location = "right", cbar_size = "2%", cbar_pad = 0.05, aspect = False)]
	grid_pos_vel_PDF     += [AxesGrid(fig_pos_vel_PDF[time], (0.01,0.01,0.99,0.99), nrows_ncols = (3, int(math.ceil(len(codes)/3.0))), axes_pad = 0.05, add_all = True, share_all = True,
					  label_mode = "1", cbar_mode = "single", cbar_location = "right", cbar_size = "2%", cbar_pad = 0.05, aspect = False)]

for time in range(len(times)):
	for code in range(len(codes)):
		# LOAD DATASETS
		if codes[code] == 'ART-I': # ART frontend doesn't find accompanying files, so we specify them; see http://yt-project.org/docs/dev/examining/loading_data.html#art-data
			dirnames = filenames[code][time][:filenames[code][time].rfind('/')+1]
			if time == 0:
				timestamp = ''
			else:
				timestamp = filenames[code][time][filenames[code][time].rfind('_'):filenames[code][time].rfind('.')] 
			pf = load(filenames[code][time], file_particle_header=dirnames+'PMcrd'+timestamp+'.DAT', file_particle_data=dirnames+'PMcrs0'+timestamp+'.DAT', file_particle_stars=dirnames+'stars'+timestamp+'.dat')
	        elif codes[code] == 'CHANGA' or codes[code] == 'GASOLINE': # For TIPSY frontend, always make sure to place your parameter file in the same directory as your datasets
			pf = load(filenames[code][time], n_ref=n_ref, over_refine_factor=over_refine_factor)
	        elif codes[code] == 'GADGET-3':
			pf = load(filenames[code][time], unit_base = gadget_default_unit_base, bounding_box=[[-1000.0, 1000.0], [-1000.0, 1000.0], [-1000.0,1000.0]], n_ref=n_ref, over_refine_factor=over_refine_factor)
                elif codes[code] == 'GEAR': # GEAR frontend currently not implemented; for this to work I removed the check for RuntimeError in read_record() of yt/utilities/fortran_utils.py
			pf = GadgetDataset(filenames[code][time], unit_base = gadget_default_unit_base, bounding_box=[[-1000.0, 1000.0], [-1000.0, 1000.0], [-1000.0,1000.0]], header_spec="default+pad256", n_ref=n_ref, over_refine_factor=over_refine_factor)
	        elif codes[code] == 'GIZMO':
			pf = load(filenames[code][time], unit_base = gadget_default_unit_base, bounding_box=[[-1000.0, 1000.0], [-1000.0, 1000.0], [-1000.0,1000.0]], field_spec="agora_unlv", n_ref=n_ref, over_refine_factor=over_refine_factor)
		else:
			pf = load(filenames[code][time])

		def _density_squared(field, data):  
			return data[("gas", "density")]**2
		pf.add_field(("gas", "density_squared"), function=_density_squared, units="g**2/cm**6")

		# ADDITIONAL FIELDS: TEMPERATURE FIELDS FOR CERTAIN CODES
                if codes[code] == 'RAMSES': 
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
				print "current_temperature = %e " % current_temperature
			def convert_T_over_mu_to_T(T_over_mu):
				logT_over_mu = np.log(T_over_mu)
				logT = np.interp(logT_over_mu, np.log(T_over_mu_values), np.log(temperature_values)) # linear interpolation in log-log space
				return np.exp(logT)
			# The pressure field includes the artificial pressure support term, so one needs to be careful (compare with the exsiting yt/frontends/ramses/fields.py)
			def _temperature_2(field, data):  
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
				return YTArray(convert_T_over_mu_to_T(T_over_mu), "K") # now T
		        pf.add_field(("gas", "temperature"), function=_temperature_2, force_override=True, units="K")

		# GAS PARTICLE FILEDS FOR SPH CODES
		if codes[code] == 'CHANGA' or codes[code] == 'GASOLINE' or codes[code] == 'GEAR' or codes[code] == 'GIZMO':
			PartType_to_use = "Gas"
			MassType_to_use = "Mass"
		elif codes[code] == "GADGET-3":
			PartType_to_use = "PartType0"				
			MassType_to_use = "Masses"

		# AXIS SWAP & NORMAL VECTOR
                if codes[code] == 'GEAR':
			# Temporary fix until GEAR group figures out the issue
			# pf.coordinates.x_axis = {0: 0, 1: 2, 2: 2, 'x': 0, 'y': 2, 'z': 2} 
			# pf.coordinates.y_axis = {0: 1, 1: 1, 2: 0, 'x': 1, 'y': 1, 'z': 0}
			pf.coordinates.x_axis[0] = 0 
			pf.coordinates.y_axis[0] = 1
			pf.coordinates.x_axis['x'] = 0
			pf.coordinates.y_axis['x'] = 1
			pf.coordinates.y_axis[1] = 1
			pf.coordinates.y_axis['y'] = 1
			pf.coordinates.x_axis[2] = 2 
			pf.coordinates.y_axis[2] = 0
			pf.coordinates.x_axis['z'] = 2
			pf.coordinates.y_axis['z'] = 0
		else:
			# From http://nbviewer.ipython.org/gist/ngoldbaum/a753d83a7f8e123b0a2c
			# pf.coordinates.x_axis = {0: 1, 1: 0, 2: 0, 'x': 1, 'y': 0, 'z': 0} 
			# pf.coordinates.y_axis = {0: 2, 1: 2, 2: 1, 'x': 2, 'y': 2, 'z': 1}
			pf.coordinates.x_axis[1] = 0
			pf.coordinates.y_axis[1] = 2
			pf.coordinates.x_axis['y'] = 0
			pf.coordinates.y_axis['y'] = 2

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
				# particle_type=False doesn't make sense, but is critical for PhasePlot to work
				# requires a change in data_objects/data_container.py: remove raise YTFieldTypeNotFound(ftype)
				def _Density_2(field, data):
					return data[(PartType_to_use, "Density")].in_units('g/cm**3')
				pf.add_field((PartType_to_use, "Density_2"), function=_Density_2, take_log=True, particle_type=False, display_name="Density", units="g/cm**3") 
				def _Temperature_2(field, data):
					return data[(PartType_to_use, "Temperature")].in_units('K')
				pf.add_field((PartType_to_use, "Temperature_2"), function=_Temperature_2, take_log=True, particle_type=False, display_name="Temperature", units="K") 
				def _Mass_2(field, data):
					return data[(PartType_to_use, MassType_to_use)].in_units('Msun')
				pf.add_field((PartType_to_use, "Mass_2"), function=_Mass_2, take_log=True, particle_type=False, display_name="Mass", units="Msun")				
				p3 = PhasePlot(sp, (PartType_to_use, "Density_2"), (PartType_to_use, "Temperature_2"), (PartType_to_use, "Mass_2"), weight_field=None, fontsize=12, x_bins=500, y_bins=500)
				p3.set_zlim((PartType_to_use, "Mass_2"), 1e3, 1e8)
				plot3 = p3.plots[(PartType_to_use, "Mass_2")]

				# p3 = ParticlePhasePlot(sp, (PartType_to_use, "Density"), (PartType_to_use, "Temperature"), (PartType_to_use, MassType_to_use), weight_field=None, fontsize=12, x_bins=500, y_bins=500)
				# p3.set_unit("Density", 'g/cm**3')
				# p3.set_unit("Temperature", 'K')
				# p3.set_unit(MassType_to_use, 'Msun') # this doesn't work?
				# p3.set_log("Density", True)
				# p3.set_log("Temperature", True)
				# p3.set_zlim((PartType_to_use, MassType_to_use), 1e3, 1e8)
				# p3.set_colorbar_label((PartType_to_use, MassType_to_use), "Mass ($\mathrm{M}_{\odot}$)")
				# plot3 = p3.plots[(PartType_to_use, MassType_to_use)]

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
		if draw_pos_vel_PDF == 1:
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
				pf.field_info[(PartType_to_use, "particle_position_cylindrical_radius")].take_log = False
				pf.field_info[(PartType_to_use, "particle_velocity_cylindrical_theta")].take_log = False
				p4 = ParticlePhasePlot(sp, (PartType_to_use, "particle_position_cylindrical_radius"), (PartType_to_use, "particle_velocity_cylindrical_theta"), \
							       (PartType_to_use, MassType_to_use), weight_field=None, fontsize=12, x_bins=300, y_bins=300)
				p4.set_unit("particle_position_cylindrical_radius", 'kpc')
				p4.set_unit("particle_velocity_cylindrical_theta", 'km/s')
				p4.set_unit(MassType_to_use, 'Msun')
				p4.set_zlim((PartType_to_use, MassType_to_use), 1e3, 1e8)
				p4.set_colorbar_label((PartType_to_use, MassType_to_use), "Mass ($\mathrm{M}_{\odot}$)")
				plot4 = p4.plots[(PartType_to_use, MassType_to_use)]

			p4.set_xlabel("Cylindrical Radius (kpc)")
			p4.set_ylabel("Rotational Velocity (km/s)")
			p4.set_xlim(0, 14)
			p4.set_ylim(-50, 350)

			plot4.figure = fig_pos_vel_PDF[time]
			plot4.axes = grid_pos_vel_PDF[time][code].axes
			if code == 0: plot4.cax = grid_pos_vel_PDF[time].cbar_axes[0]
			p4._setup_plots()

			if add_nametag == 1:
				at = AnchoredText("%s" % codes[code], loc=4, prop=dict(size=10), frameon=True)
				grid_pos_vel_PDF[time][code].axes.add_artist(at)


if draw_density_map == 1:
	for time in range(len(times)):
		fig_density_map[time].savefig("Sigma_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)
if draw_temperature_map == 1:
	for time in range(len(times)):
		fig_temperature_map[time].savefig("Temp_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)
if draw_PDF == 1:
	for time in range(len(times)):
		fig_PDF[time].savefig("PDF_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)
if draw_pos_vel_PDF == 1:
	for time in range(len(times)):
		fig_pos_vel_PDF[time].savefig("pos_vel_PDF_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)
