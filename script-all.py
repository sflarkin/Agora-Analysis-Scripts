#####################################################################
#
#  UNIFIED ANALYSIS SCRIPT EXAMPLE FOR THE AGORA PROJECT
#
#  FOR SCRIPT HISTORY SEE VERSION CONTROL CHANGELOG
#
#####################################################################

import sys, os
for spec in ["~/yt/yt-3.0", "~/yt-3.0"]:
    if os.path.isdir(os.path.expanduser(spec)):
        sys.path.insert(0, os.path.expanduser(spec))
        break

if not os.path.isdir("images"): os.makedirs("images")
import h5py
from yt.config import ytcfg; ytcfg["yt","loglevel"] = "20"
from yt.mods import *
from yt.utilities.physical_constants import kpc_per_cm
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from yt.data_objects.particle_filters import \
        particle_filter, filter_registry
from yt.data_objects.particle_fields import \
        particle_deposition_functions
from yt.frontends.stream.data_structures import \
        load_particles
from yt.data_objects.particle_fields import \
    particle_deposition_functions, \
    particle_vector_functions

@particle_filter("finest", ["ParticleMassMsun"])
def finest(pfilter, data):
    return data["ParticleMassMsun"] < 340000

def main():
    # CHOOSE DATASETS YOU WANT TO PROCESS
    do_art2 = False
    do_enzo = False
    do_ramses = False
    do_gadget = False
    do_gasoline = False
    do_pkdgrav = False
    do_particles = False

    if do_art2:
        ds_art2 = load("A11QR1/s11Qzm1h2_a1.0000.art")
        ds_art2.add_particle_filter("finest")
        center = 128 * np.array([0.492470,  0.533444,  0.476942]) 
        process_dataset(ds_art2, center)

    if do_enzo:
        ds_enzo = load("./DD0040/data0040")
        ds_enzo.add_particle_filter("finest")
        center = np.array([ 0.49297869, 0.50791068, 0.50727271])
        process_dataset(ds_enzo, center)
        
    if do_ramses:
        ds_ramses = load("./output_00101/info_00101.txt")
        ds_ramses.add_particle_filter("finest")
        center = np.array([ 0.48598457, 0.52665735, 0.48984628])
        process_dataset(ds_ramses, center)

    if do_gadget:
        ds_gadget = GadgetStaticOutput("./snapshot_010", unit_base = {"mpchcm": 1.0})
        ds_gadget.add_particle_filter("finest")
        center = np.array([ 29.75540543, 32.12417221, 28.28912735])
        process_dataset(ds_gadget, center)

    if do_gasoline:
        cosmology_parameters = dict(current_redshift = 0.0, omega_lambda = 0.728,
                                    omega_matter = 0.272, hubble_constant = 0.702)
        ds_gasoline = TipsyStaticOutput("./agora_1e11.00400",
                                        cosmology_parameters = cosmology_parameters,
                                        unit_base = {'mpchcm': 1.0/60.0})
        ds_gasoline.add_particle_filter("finest")
        center = np.array([-0.014738, 0.026979, -0.010535])
        #center = np.array([-0.01477163, 0.02694199, -0.0105199])
        process_dataset(ds_gasoline, center)

    if do_pkdgrav:
        cosmology_parameters = dict(current_redshift = 0.0, omega_lambda = 0.728,
                                    omega_matter = 0.272, hubble_constant = 0.702)
        ds_pkdgrav = TipsyStaticOutput("./halo1e11_run1.00400", endian="<",
                                       field_dtypes = {"Coordinates": "d"},
                                       cosmology_parameters = cosmology_parameters,
                                       unit_base = {'mpchcm': 1.0/60.0})
        ds_pkdgrav.add_particle_filter("finest")
        center = np.array([-0.01434195, 0.027505, -0.01086525])
        process_dataset(ds_pkdgrav, center)

    if do_particles:
        f = h5py.File("./s11Qzm1h2_a1.0000.art.h5")
        data = dict((k, f[k][:].astype("float64")) for k in f)
        bbox = np.array([[0.0, 128.0], [0.0, 128.0], [0.0, 128.0]])
        data["particle_mass"] *= 2.2023338912587828e+43 # Convert to grams
        ds_particles = load_particles(data, 2.0604661199638546e+24, bbox=bbox)
        ds_particles.add_particle_filter("finest")
        particle_vector_functions("all", ["particle_position_%s" % ax for ax in 'xyz'],
                                  ["particle_velocity_%s" % ax for ax in 'xyz'], StreamFieldInfo)
        particle_vector_functions("finest", ["particle_position_%s" % ax for ax in 'xyz'],
                                  ["particle_velocity_%s" % ax for ax in 'xyz'], StreamFieldInfo)
        particle_deposition_functions("all", "Coordinates", "particle_mass", StreamFieldInfo)
        particle_deposition_functions("finest", "Coordinates", "particle_mass", StreamFieldInfo)
        center = np.array([0.492470,  0.533444,  0.476942]) 
        center *= 128
        process_dataset(ds_particles, center)

def process_dataset(ds, center):

    #=======================
    #  [1] TOTAL MASS
    #=======================
    sp = ds.h.sphere(center, (1.0, 'mpc'))
    total_particle_mass = sp.quantities["TotalQuantity"]( ("all","ParticleMassMsun") )[0]
    print "Total particle mass within a radius of 1 Mpc of the center: %0.3e Msun" \
          % total_particle_mass
    
    #=======================
    #  [2] PLOTS
    #=======================
    for ptype in ["finest", "all"]:
        p = ProjectionPlot(ds, "z", ("deposit", "%s_density" % ptype), center = center)
        p.save("./images/%s_z1_%s.png" % (ds, ptype))
        p.zoom(60)
        p.save("./images/%s_z2_%s.png" % (ds, ptype))

    #=======================
    #  [3] PLOTS-2
    #=======================
    w = (1.0/0.702, "mpc")
    #w = (1.0, "mpch")
    axis = 2
    colorbounds = (1e-32, 1e-25)
    res = [1024] * 3
    res[axis] = 16

    LE = center - 0.5*(w[0]/ds[w[1]])
    RE = center + 0.5*(w[0]/ds[w[1]])
    source = ds.h.arbitrary_grid(LE, RE, res)

    field = ("deposit", "all_cic")
    num = (source[field] * source[field]).sum(axis=axis)
    num *= (RE[axis] - LE[axis])*ds['cm'] # dl
    den = (source[field]).sum(axis=axis)
    den *= (RE[axis] - LE[axis])*ds['cm'] # dl
    proj = (num/den)
    proj[proj!=proj] = 1e-100 # remove NaN's
    plt.clf()
    norm = LogNorm(colorbounds[0], colorbounds[1], clip=True)
    plt.imshow(proj.swapaxes(0,1), interpolation='nearest', origin='lower',
               norm = norm, extent = [-0.5*(w[0]/ds[w[1]]), 0.5*(w[0]/ds[w[1]]), 
                                       -0.5*(w[0]/ds[w[1]]), 0.5*(w[0]/ds[w[1]])])
    plt.xlabel(r"$%d\/ \mathrm{Mpc} / h \/(\mathrm{comoving})$" \
                   % round(ds.units["mpc"]/0.702))
                   # % round(ds.units["mpch"]))
    plt.ylabel(r"$%d\/ \mathrm{Mpc} / h \/(\mathrm{comoving})$" \
                   % round(ds.units["mpc"]/0.702))
                   # % round(ds.units["mpch"]))
    cb = plt.colorbar()
    cb.set_label(r"$\mathrm{Density}\/\/[\mathrm{g}/\mathrm{cm}^3]$")
    plt.savefig("./images/%s_%s.png" % (ds, field[1]), dpi=150, bbox_inches='tight', \
                pad_inches=0.1)
    
    #=======================
    #  [4] PROFILES
    #=======================
    sphere_radius        = 300  # kpc
    inner_radius         = 0.8  # kpc
    total_bins           = 30

    sp = ds.h.sphere(center, (sphere_radius, 'kpc'))
    prof = BinnedProfile1D(sp, total_bins, "ParticleRadiuskpc", \
                           inner_radius, sphere_radius, end_collect = True)
    prof.add_fields([("all","ParticleMassMsun")],
                    weight = None, accumulation=False)
    prof["AverageDMDensity"] = prof[("all","ParticleMassMsun")] * 6.77e-32
    for k in range(0, len(prof["AverageDMDensity"])): # g/cm^3
        if k == 0:
            prof["AverageDMDensity"][k] /= ((4.0/3.0)*np.pi*prof["ParticleRadiuskpc"][k]**3)  
        else:
            prof["AverageDMDensity"][k] /= ((4.0/3.0)*np.pi* \
                                            (prof["ParticleRadiuskpc"][k]**3 - \
                                             prof["ParticleRadiuskpc"][k-1]**3))

    plt.clf()
    plt.loglog(prof["ParticleRadiuskpc"], prof["AverageDMDensity"], '-k')
    plt.xlabel(r"$\mathrm{Radius}\/\/[\mathrm{kpc}]$")
    plt.ylabel(r"$\mathrm{Dark}\/\mathrm{Matter}\/\mathrm{Density}\/\/[\mathrm{g}/\mathrm{cm}^3]$")
    plt.ylim(1e-29, 1e-23)
    plt.savefig("./images/%s_radprof.png" % ds)

    fout = open("./images/%s_profile.dat" % ds, "w")
    fout.write("# sphere_radius:"+str(sphere_radius)+"\n")
    fout.write("# inner_radius:"+str(inner_radius)+"\n")
    fout.write("# \n")
    fout.write("# radius(kpc)   gas_DM_enclosed(Msun)   prof_DM_shell(g/cm^3)\n")
    fout.write("# \n")
    for k in range(0, total_bins):
        fout.write(str(prof["ParticleRadiuskpc"][k])+"   "+ 
                   str(prof[("all","ParticleMassMsun")][k])+"   "+ 
                   str(prof["AverageDMDensity"][k])+"\n")
    fout.close()

    #=======================
    #  [5] HOP HALOFINDER
    #=======================
    if os.path.exists("./images/%s_MergerHalos.out" % ds) == False:
        omega_matter = 0.272
        halos = HaloFinder(ds, threshold = omega_matter*360)
        field = ("deposit", "finest_density")
        halos.write_out("./images/%s_MergerHalos.out" % ds)
        
        source = ds.h.region(center, center - (w[0]/ds[w[1]])/2.0, center + (w[0]/ds[w[1]])/2.0)
        proj = ds.h.proj(field, 2, weight_field = field, data_source = source)
        pw = proj.to_pw(fields = [field], center = center, width = w)
        pw.set_zlim(field, 1e-31, 1e-24)
        pw.annotate_hop_circles(halos, annotate=True, print_halo_mass=False, min_size=1000)
        pw.save("./images/%s_Projection_z_%s_subset.png" % (ds, field[1]))

if __name__ == '__main__':
    main()
