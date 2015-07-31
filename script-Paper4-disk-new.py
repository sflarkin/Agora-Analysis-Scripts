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
 	     ['./GASOLINE/LOW_nosf_nofb_gasoline_pfloor_Zsolar_0Myr.00001', './GASOLINE/LOW_nosf_nofb_gasoline_pfloor_Zsolar.00335'],
  	     ['./GEAR/snapshot_0000', './GEAR/snapshot_0500'],
 	     ['./GIZMO/snapshot_000', './GIZMO/snapshot_100_hsml'],
 	     ['./RAMSES/output_00001/info_00001.txt', './RAMSES/output_00068/info_00068.txt']]
gadget_default_unit_base = {'UnitLength_in_cm'         : 3.08568e+21,
			    'UnitMass_in_g'            :   1.989e+43,
			    'UnitVelocity_in_cm_per_s' :      100000}

draw_density_map     = 1         # 1/0 = ON/OFF
draw_temperature_map = 1         # 1/0 = ON/OFF
draw_PDF             = 1         # 1/0 = ON/OFF
add_nametag          = 1         # 1/0 = ON/OFF
times                = [0, 500]  # in Myr
figure_width         = 30        # in kpc
n_ref                = 256       # for SPH codes
over_refine_factor   = 1         # for SPH codes

fig_density_map      = [] 
fig_temperature_map  = []
fig_PDF              = []
grid_density_map     = []
grid_temperature_map = []
grid_PDF             = []

for time in range(len(times)):
	fig_density_map      += [plt.figure(figsize=(100,20))]
	fig_temperature_map  += [plt.figure(figsize=(100,20))]
	fig_PDF              += [plt.figure(figsize=(50, 80))]
	grid_density_map     += [AxesGrid(fig_density_map[time], (0.01,0.01,0.99,0.99), nrows_ncols = (2, len(codes)), axes_pad = 0.02, add_all = True, share_all = True,
					  label_mode = "1", cbar_mode = "single", cbar_location = "right", cbar_size = "2%", cbar_pad = 0.02)]
	grid_temperature_map += [AxesGrid(fig_temperature_map[time], (0.01,0.01,0.99,0.99), nrows_ncols = (2, len(codes)), axes_pad = 0.02, add_all = True, share_all = True,
					  label_mode = "1", cbar_mode = "single", cbar_location = "right", cbar_size = "2%", cbar_pad = 0.02)]
	grid_PDF             += [AxesGrid(fig_PDF[time], (0.01,0.01,0.99,0.99), nrows_ncols = (3, int(math.ceil(len(codes)/3.0))), axes_pad = 0.05, add_all = True, share_all = True,
					  label_mode = "1", cbar_mode = "single", cbar_location = "right", cbar_size = "2%", cbar_pad = 0.05, aspect = False)]

for time in range(len(times)):
	for code in range(len(codes)):
		# LOAD DATASETS
		if codes[code] == 'ART-I': 
			dirnames = filenames[code][time][:filenames[code][time].rfind('/')+1]
			if time == 0:
				timestamp = ''
			else:
				timestamp = filenames[code][time][filenames[code][time].rfind('_'):filenames[code][time].rfind('.')] 
			pf = load(filenames[code][time], file_particle_header=dirnames+'PMcrd'+timestamp+'.DAT', file_particle_data=dirnames+'PMcrs0'+timestamp+'.DAT', file_particle_stars=dirnames+'stars'+timestamp+'.dat')
	        elif codes[code] == 'CHANGA' or codes[code] == 'GASOLINE': 
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
		# As of July 2015 Ramses units have ongoing issues in yt-3.2; see https://bitbucket.org/yt_analysis/yt/issues/1055/ramses-units-error-with-boxlen-1 
                if codes[code] == 'RAMSES': 
			# from utilities:convenience.py:   Calculate a tabulated approximation to mean molecular weight (valid for data that used Grackle 2.0 or below)
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
				
			def _temperature_2(field, data):  # the pressure field includes the artificial pressure support term, so one needs to be careful (compare with the exsiting yt/frontends/ramses/fields.py)
				T_J = 1800.0  # in K
				n_J = 8.0     # in H/cc
				gamma_0 = 2.0 
				x_H = 0.76
				if time != 0:
					T_over_mu = data["gas", "pressure"].d/data["gas", "density"].d * constants.mass_hydrogen_cgs.d/constants.boltzmann_constant_cgs.d \
					    - T_J * (data["gas", "density"].d * x_H /constants.mass_hydrogen_cgs.d / n_J)**(gamma_0 - 1.0) # T/mu = T2 in Ramses
				else:
					T_over_mu = data["gas", "pressure"].d/data["gas", "density"].d * constants.mass_hydrogen_cgs.d/constants.boltzmann_constant_cgs.d  # no pressure support in IC
				return YTArray(convert_T_over_mu_to_T(T_over_mu), "K") # now T
		        pf.add_field(("gas", "temperature_2"), function=_temperature_2, units="K")

		# AXIS SWAP
                if codes[code] == 'GEAR':
			# temporary fix until GEAR group figures out the issue
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
			# from http://nbviewer.ipython.org/gist/ngoldbaum/a753d83a7f8e123b0a2c
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
					     center + YTArray([figure_width, figure_width, figure_width], 'kpc'))
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
				if codes[code] == 'RAMSES': 
					p2 = ProjectionPlot(pf, ax, ("gas", "temperature_2"), center = center, data_source=proj_region, width = (figure_width, 'kpc'), weight_field = ("gas", "density_squared"), fontsize=9)
					p2.set_zlim(("gas", "temperature_2"), 1e1, 1e6)
					plot2 = p2.plots[("gas", "temperature_2")]
				else:
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
			if codes[code] == 'RAMSES': 
				p3 = PhasePlot(sp, ("gas", "density"), ("gas", "temperature_2"), ("gas", "cell_mass"), weight_field=None, fontsize=12, x_bins=500, y_bins=500)
				p3.set_unit("cell_mass", 'Msun')
				p3.set_zlim(("gas", "cell_mass"), 1e3, 1e8)
				plot3 = p3.plots[("gas", "cell_mass")]
			else:
				p3 = PhasePlot(sp, ("gas", "density"), ("gas", "temperature"), ("gas", "cell_mass"), weight_field=None, fontsize=12, x_bins=500, y_bins=500)
				p3.set_unit("cell_mass", 'Msun')
				p3.set_zlim(("gas", "cell_mass"), 1e3, 1e8)
				plot3 = p3.plots[("gas", "cell_mass")]
			p3.set_xlim(1e-29, 1e-21)
			p3.set_ylim(10, 1e7)

			plot3.figure = fig_PDF[time]
			plot3.axes = grid_PDF[time][code].axes
			if code == 0: plot3.cax = grid_PDF[time].cbar_axes[0]

			p3._setup_plots()

			if add_nametag == 1:
				at = AnchoredText("%s" % codes[code], loc=3, prop=dict(size=10), frameon=True)
				grid_PDF[time][code].axes.add_artist(at)


if draw_density_map == 1:
	for time in range(len(times)):
		fig_density_map[time].savefig("Sigma_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)
if draw_temperature_map == 1:
	for time in range(len(times)):
		fig_temperature_map[time].savefig("Temp_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)
if draw_PDF == 1:
	for time in range(len(times)):
		fig_PDF[time].savefig("PDF_%dMyr" % times[time], bbox_inches='tight', pad_inches=0.03, dpi=300)
