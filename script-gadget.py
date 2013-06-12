#####################################################################
#
#  YT SCRIPT TO PRODUCE PROJECTED DARK MATTER DENSITY 
#  USING ONLY THE FINEST DM PARTICLES ON YT-3.0 alpha release 2
#  
#  PLEASE SEE:  https://hub.yt-project.org/nb/abu5nb
#
#  FUNCTIONALITY AVAILABLE GREATLY THANKS TO:  Matthew Turk et al.
#
#  SCRIPT WRITTEN BY:  Ji-hoon Kim on June 6, 2013
#
#####################################################################

import sys
sys.path.insert(0, "/home/mornkr/yt-3.0")
from yt.config import ytcfg; ytcfg["yt","loglevel"] = "20"
from yt.mods import *
import yt.utilities.lib as au
import numpy as np
import copy
import pylab
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from yt.analysis_modules.halo_finding.api import *
from yt.geometry.oct_container import ParticleOctreeContainer

def finest_DM_density(field, data):
    filter = np.where(data["ParticleMassMsun"] <= 340000)

    pos = data["all", "Coordinates"][filter]
    d = data.deposit(pos, [data["all", "Mass"][filter]], method = "cic")
    d /= data["CellVolume"]
    return d

GadgetFieldInfo.add_field(("deposit", "finest_DM_density"),
                          function = finest_DM_density,
                          validators = [ValidateSpatial()],
                          display_name = "\\mathrm{Finest DM Density}",
                          units = r"\mathrm{g}/\mathrm{cm}^{3}",
                          projected_units = r"\mathrm{g}/\mathrm{cm}^{2}",
                          projection_conversion = 'cm')

#center = np.array([29.11572,  31.61874,  29.3679])
center = np.array([29.754, 32.14, 28.29]) # Gadget unit system: [0, 60]
#center = np.array([ 29.76579666,  32.57726669,  28.76054764]) # all_density maximum location found with ds.h.find_max(('deposit', 'all_density')) in yt 
#center = np.array([ 29.81300354,  32.12242126,  28.33488464])  
ds = GadgetStaticOutput("snapshot_010", unit_base = {"mpchcm": 1.0})
print ds.h.oct_handler.n_ref
ds.h.oct_handler.n_ref = 2
print ds.h.oct_handler.n_ref
#ds = GadgetStaticOutput("snapshot_agora_adapt_noUNEQUAL_011", unit_base = {"mpchcm": 1.0})
#ds = load("snapshot_010", root_dimensions=[2, 2, 2])
print ds.parameters["Npart"], ds.parameters["Nall"]
print ds.units["mpchcm"]
print ds.units["mpch"]
print ds.units["cm"]

sp = ds.h.sphere(center, (1.0, 'mpc'))
total_particle_mass = sp.quantities["TotalQuantity"]( [("all","ParticleMassMsun")] )[0]
print "Total particle mass within a radius of 1 Mpc of the center: %0.3e Msun" % total_particle_mass
print ds.h.derived_field_list

pw = ProjectionPlot(ds, "z", ("deposit", "all_density"), weight_field=None, center=center, width=(1.0, 'mpch')).save()

w = (1.0, "mpch")
source = ds.h.region(center, center - (w[0]/ds[w[1]])/2.0, center + (w[0]/ds[w[1]])/2.0)
proj = ds.h.proj( ("deposit", "all_density"), 2, weight_field = ("deposit", "all_density"), data_source = source)
pw = proj.to_pw(fields = [("deposit", "all_density")], center = center, width = w)
pw.set_zlim(("deposit","all_density"), 1e-32, 1e-25)
pw.save("snapshot_010_Projection_z_all_density_subset.png")
#pw.save("snapshot_agora_adapt_noUNEQUAL_011_Projection_z_all_density_subset.png")

proj = ds.h.proj( ("deposit", "finest_DM_density"), 2, weight_field = ("deposit", "finest_DM_density"), data_source = source)
pw = proj.to_pw(fields = [("deposit", "finest_DM_density")], center = center, width = w)
pw.set_zlim(("deposit","finest_DM_density"), 1e-32, 1e-25)
pw.save("snapshot_010_Projection_z_finest_DM_density_subset.png")
#pw.save("snapshot_agora_adapt_noUNEQUAL_011_Projection_z_finest_DM_density_subset.png")



sphere_radius        = 300  # kpc
inner_radius         = 1    # kpc
total_bins           = 10
PI                   = 3.141592
dlogRadius           = (np.log10(sphere_radius) - np.log10(inner_radius)) / (total_bins-1)
prof_radius     = numpy.zeros([total_bins], float)
prof_DM         = numpy.zeros([total_bins], float)
prof_DM_shell   = numpy.zeros([total_bins], float)

for k in range(0, total_bins):
    prof_radius[k] = pow(10, np.log10(inner_radius) + float(k)*dlogRadius);
    prof_sphere = ds.h.sphere(center, (prof_radius[k], 'kpc'))
    prof_DM[k] = prof_sphere.quantities["TotalQuantity"]( [("all","ParticleMassMsun")] )[0]
    my_sphere = ds.h.sphere(center, (sphere_radius, 'kpc'))
    if k == 0:
        my_shell = numpy.where(my_sphere["Radiuskpc"] < prof_radius[k])
        prof_DM_shell[k] = prof_DM[k] / (4.0/3.0*PI*(prof_radius[k])**3) * 6.77e-32
    else:
        my_shell = numpy.where((my_sphere["Radiuskpc"] < prof_radius[k]) & (my_sphere["Radiuskpc"] > prof_radius[k-1]))
        prof_DM_shell[k] = (prof_DM[k] - prof_DM[k-1]) / (4.0/3.0*PI*(prof_radius[k])**3 - (4.0/3.0*PI*(prof_radius[k-1])**3)) * 6.77e-32

fout = open("profile-gadget.dat","a")
fout.write("# sphere_radius:"+str(sphere_radius)+"\n")
fout.write("# inner_radius:"+str(inner_radius)+"\n")
fout.write("# \n")
fout.write("# radius(kpc)   gas_DM_enclosed(Msun)   prof_DM_shell(g/cm^3)\n")
fout.write("# \n")
for k in range(0, total_bins):
    fout.write(str(prof_radius[k])+"   "+str(prof_DM[k])+"   "+str(prof_DM_shell[k]))
    fout.write("\n")
fout.close()



# halos = HaloFinder(ds, subvolume = source, threshold=80.)
# print halos[0].center_of_mass() 
# print halos[0].total_mass()
# halos.dump("./MergerHalos-gadget")
# pw = ProjectionPlot(ds, "z", ("deposit", "all_density"), weight_field=None, center=center, width=(1.0,'mpch'))
# pw.annotate_hop_circles(halos)
# pw.save()
