#####################################################################
#
#  AGORA SCRIPT
#  
#  PLEASE SEE:  https://hub.yt-project.org/nb/abu5nb
#
#  FOR SCRIPT HISTORY SEE VERSION CONTROL CHANGELOG
#
#####################################################################

import sys, os
for spec in ["~/yt/yt-3.0", "~/yt-3.0"]:
    if os.path.isdir(os.path.expanduser(spec)):
        sys.path.insert(0, os.path.expanduser(spec))
        break
if not os.path.isdir("figures"):
    os.makedirs("figures")
from yt.config import ytcfg; ytcfg["yt","loglevel"] = "20"
from yt.mods import *
from yt.utilities.physical_constants import kpc_per_cm
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

def particle_count(field, data):
    return np.ones(data["all","ParticleMass"].shape, dtype="float64")
GadgetFieldInfo.add_field(("all", "particle_count"), function=particle_count,
                          particle_type = True)

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
    # arbitrary_grid objects.
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
total_bins           = 128
PI                   = 3.141592
dlogRadius           = (np.log10(sphere_radius) - np.log10(inner_radius)) / (total_bins-1)
prof_radius     = numpy.zeros([total_bins], float)
prof_DM         = numpy.zeros([total_bins], float)
prof_DM_shell   = numpy.zeros([total_bins], float)

sp = ds.h.sphere(center, (300.0, 'kpc'))
prof = BinnedProfile1D(sp, total_bins, "Radiuskpc",
                       inner_radius, sphere_radius,
                       end_collect = True)
prof.add_fields([("all","ParticleMassMsun"), ("all", "particle_count")],
                weight = None, accumulation=True)
prof["AverageDMDensity"] = (prof["all","ParticleMassMsun"] /
                           ((4.0/3.0) * np.pi * prof["Radiuskpc"]**3))

plt.clf()
plt.loglog(prof["Radiuskpc"], prof["AverageDMDensity"], '-k')
plt.xlabel(r"$\mathrm{Radius}\/\/[\mathrm{kpc}]$")
plt.ylabel(r"$\mathrm{Dark}\/\mathrm{Matter}\/\mathrm{Density}\/\/[M_\odot/\mathrm{kpc}^3]$")
plt.savefig("figures/%s_radprof.png" % ds)

plt.clf()
plt.loglog(prof["Radiuskpc"], prof["all", "particle_count"], '-k')
plt.xlabel(r"$\mathrm{Radius}\/\/[\mathrm{kpc}]$")
plt.ylabel(r"$\mathrm{N}$")
plt.savefig("figures/%s_pcount.png" % ds)

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
