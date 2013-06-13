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

import sys, os
for spec in ["~/yt/yt-3.0", "~/yt-3.0"]:
    if os.path.isdir(os.path.expanduser(spec)):
        sys.path.insert(0, os.path.expanduser(spec))
        break
from yt.config import ytcfg; ytcfg["yt","loglevel"] = "20"
from yt.mods import *
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize

for name, method in [("CIC", "cic"), ("Density", "sum")]:
    def _func(_method):
        def finest_DM_func(field, data): # user-defined field
            filter = data["ParticleMassMsun"] <= 340000
            pos = data["all", "Coordinates"][filter, :]
            d = data.deposit(pos, [data["all", "Mass"][filter]],
                             method = _method)
            d /= data["CellVolume"]
            return d
        return finest_DM_func
    GadgetFieldInfo.add_field(("deposit", "finest_DM_%s" % name.lower()),
                              function = _func(method),
                              validators = [ValidateSpatial()],
                              display_name = "\\mathrm{Finest DM %s}" % name,
                              units = r"\mathrm{g}/\mathrm{cm}^{3}",
                              projected_units = r"\mathrm{g}/\mathrm{cm}^{2}",
                              projection_conversion = 'cm')

center = np.array([29.754, 32.14, 28.29]) # Gadget unit system: [0, 60]
ds = GadgetStaticOutput("snapshot_010", unit_base = {"mpchcm": 1.0})

#=======================
#  [1] TOTAL MASS
#=======================

sp = ds.h.sphere(center, (1.0, 'mpc'))
total_particle_mass = sp.quantities["TotalQuantity"]( ("all","ParticleMassMsun") )[0]
print "Total particle mass within a radius of 1 Mpc of the center: %0.3e Msun" % total_particle_mass


#=======================
#  [2] BASIC PLOTS
#=======================

w = (1.0, "mpch")
axis = 2

res = [1024] * 3
res[axis] = 16

LE = center - 0.5/ds['mpch']
RE = center + 0.5/ds['mpch']

source = ds.h.arbitrary_grid(LE, RE, res)

fields = [("deposit", "all_cic"), ("deposit", "finest_DM_cic")]

for field in fields:
    # Manually do this until we have a solution in place to do it from
    # arbitrary_grid objects
    num = (source[field] * source[field]).sum(axis=axis)
    num *= (RE[axis] - LE[axis])*ds['cm'] # dl
    den = (source[field]).sum(axis=axis)
    den *= (RE[axis] - LE[axis])*ds['cm'] # dl
    proj = (num/den)
    plt.clf()
    norm = LogNorm(1e-32, 1e-25)
    plt.imshow(proj, interpolation='nearest', origin='lower',
               norm = norm, extent = [-0.5, 0.5, -0.5, 0.5])
    plt.xlabel(r"$\mathrm{Mpc} / h (\mathrm{comoving})$")
    plt.ylabel(r"$\mathrm{Mpc} / h (\mathrm{comoving})$")
    cb = plt.colorbar()
    cb.set_label(r"$\mathrm{Density}\/\/[\mathrm{g}/\mathrm{cm}^3]$")
    plt.savefig("figures/%s_%s.png" % (ds, field[1]))

#=======================
#  [3] BASIC PROFILE
#=======================

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


#=======================
#  [4] HOP HALOFINDER
#=======================

# halos = HaloFinder(ds, subvolume = source, threshold=80.)
# print halos[0].center_of_mass() 
# print halos[0].total_mass()
# halos.dump("./MergerHalos-gadget")
# pw = ProjectionPlot(ds, "z", ("deposit", "all_density"), weight_field=None, center=center, width=(1.0,'mpch'))
# pw.annotate_hop_circles(halos)
# pw.save()
