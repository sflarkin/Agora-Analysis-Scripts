import sys, os
for spec in ["~/yt/yt-3.0", "~/yt-3.0"]:
    if os.path.isdir(os.path.expanduser(spec)):
        sys.path.insert(0, os.path.expanduser(spec))
        break

if not os.path.isdir("images"): os.makedirs("images")
from yt.config import ytcfg; ytcfg["yt","loglevel"] = "20"
from yt.mods import *
from yt.data_objects.particle_filters import \
        particle_filter, filter_registry
from yt.data_objects.particle_fields import \
        particle_deposition_functions

@particle_filter("finest", ["ParticleMassMsun"])
def finest(pfilter, data):
    return data["ParticleMassMsun"] < 340000

def process_dataset(ds, center):
    
    for ptype in ["finest", "all"]:
        t1 = time.time()
        p = ProjectionPlot(ds, "x", ("deposit", "%s_density" % ptype),
                           center = center)
        p.save("images/%s_z1_%s.png" % (ds, ptype))
        p.zoom(10)
        p.save("images/%s_z2_%s.png" % (ds, ptype))
        t2 = time.time()
        print "Took %0.3e for %s" % (t2-t1, ptype)

do_enzo = True
do_ramses = True
do_gadget = True
do_gasoline = True
do_pkdgrav = True

if do_ramses:
    ds_ramses = load("output_00101/info_00101.txt")
    ds_ramses.add_particle_filter("finest")
    center = np.array([ 0.48598457, 0.52665735, 0.48984628])
    process_dataset(ds_ramses, center)

if do_gadget:
    ds_gadget = GadgetStaticOutput("snapshot_010", unit_base = {"mpchcm": 1.0})
    ds_gadget.add_particle_filter("finest")
    center = np.array([ 29.75540543, 32.12417221, 28.28912735])
    process_dataset(ds_gadget, center)

if do_gasoline:
    cosmology_parameters = dict(current_redshift = 0.0, omega_lambda = 0.728,
        omega_matter = 0.272, hubble_constant = 0.702)
    ds_gasoline = TipsyStaticOutput("agora_1e11.00400",
        cosmology_parameters = cosmology_parameters,
        unit_base = {'mpchcm': 1.0/60.0})
    ds_gasoline.add_particle_filter("finest")
    center = np.array([-0.014738, 0.026979, -0.010535])
    process_dataset(ds_gasoline, center)

if do_pkdgrav:
    cosmology_parameters = dict(current_redshift = 0.0, omega_lambda = 0.728,
        omega_matter = 0.272, hubble_constant = 0.702)
    ds_pkdgrav = TipsyStaticOutput("halo1e11_run1.00400", endian="<",
        field_dtypes = {"Coordinates": "d"},
        cosmology_parameters = cosmology_parameters,
        unit_base = {'mpchcm': 1.0/60.0})
    ds_pkdgrav.add_particle_filter("finest")
    center = np.array([-0.01434195, 0.027505, -0.01086525])
    process_dataset(ds_pkdgrav, center)

if do_enzo:
    ds_enzo = load("DD0040/data0040")
    ds_enzo.add_particle_filter("finest")
    center = np.array([ 0.49297869, 0.50791068, 0.50727271])
    process_dataset(ds_enzo, center)

