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
from yt.data_objects.particle_filters import \
        particle_filter, filter_registry
from yt.data_objects.particle_fields import \
        particle_deposition_functions, \
        particle_scalar_functions, \
        particle_vector_functions 

@particle_filter("finest", ["ParticleMassMsun"])
def finest(pfilter, data):
    return data["ParticleMassMsun"] < 340000

@derived_field(name=("deposit", "finest_density"),
               validators = [ValidateSpatial()],
               display_name = "\\mathrm{Finest DM %s}" % name,
               units = r"\mathrm{g}/\mathrm{cm}^{3}",
               projected_units = r"\mathrm{g}/\mathrm{cm}^{2}",
               projection_conversion = 'cm')
def finest_deposit(field, data):
    pos = data["finest", "Coordinates"]
    d = data.deposit(pos, [data["finest","ParticleMassMsun"]], method="sum")
    return d / data["CellVolume"]


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

fo = filter_registry["finest"][0]

do_enzo = True
do_ramses = True
do_gadget = True

if do_enzo:
    ds_enzo = load("DD0040/data0040")
    ds_enzo.add_particle_filter(fo)
    df = particle_deposition_functions("finest",
        "Coordinates", "ParticleMass", EnzoFieldInfo)
    ds_enzo.h._derived_fields_add(df)
    center = np.array([ 0.49297869, 0.50791068, 0.50727271])
    process_dataset(ds_enzo, center)

if do_ramses:
    ds_ramses = load("output_00101/info_00101.txt")
    ds_ramses.add_particle_filter(fo)
    df = particle_deposition_functions("finest",
        "Coordinates", "ParticleMass", RAMSESFieldInfo)
    ds_ramses.h._derived_fields_add(df)
    center = np.array([ 0.48598457, 0.52665735, 0.48984628])
    process_dataset(ds_ramses, center)

if do_gadget:
    ds_gadget = GadgetStaticOutput("snapshot_010", unit_base = {"mpchcm": 1.0})
    ds_gadget.add_particle_filter(fo)
    df = particle_deposition_functions("finest",
        "Coordinates", "Mass", GadgetFieldInfo)
    ds_gadget.h._derived_fields_add(df)
    center = np.array([ 29.75540543, 32.12417221, 28.28912735])  # Gadget unit system: [0, 60]
    process_dataset(ds_gadget, center)
