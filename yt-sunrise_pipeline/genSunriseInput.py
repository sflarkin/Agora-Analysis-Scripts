'''
Generate the cameras to use in Sunrise and make projection plots
of the data for some of these cameras. Then export the data within
the fov to a FITS file in a format that Sunrise understands.

by Miguel Rocha  - miguel@scitechanalytics.com
'''
import os, sys, argparse
from glob import glob
import numpy as np
from collections import OrderedDict


if __name__ != "__main__":
    from yt.analysis_modules.sunrise_export import sunrise_octree_exporter
    from yt.analysis_modules.halo_finding.halo_objects import RockstarHaloList 
    yt.enable_parallelism()


def parse():
    '''
    Parse command line arguments
    ''' 
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description='''\
                                Generate the cameras to use in Sunrise and make projection plots
                                of the data for some of these cameras. Then export the data within
                                the fov to a FITS file in a format that Sunrise understands.
                                ''')
 
    parser.add_argument('sim_dirs', nargs='+', help='Simulation directories to be analyzed.')
    
    parser.add_argument('-s', '--snap_base', default='10MpcBox_csf512_',
                        help='Base of the snapshots file names.') 

    parser.add_argument('-d', '--distance', default=100, type=float,
                        help='Distance between cameras and the center of the galaxy (in [kpc]).')

    parser.add_argument('-f', '--fov', default=50, type=float,
                        help='Field of view of the cameras at the image plane (in [kpc]).') 

    parser.add_argument('--star_particles', default='stars',
                        help='The name given to the star particles in the yt dataset '\
                        '(as printed in ds.field_list or ds.derived_field_list).') 

    parser.add_argument('--dm_particles', default='darkmatter',
                        help='The name given to the dark matter particles in the yt dataset '\
                        '(as printed in ds.field_list or ds.derived_field_list).') 

    parser.add_argument( '--galprops_file', default='sim_dir/analysis/catalogs/*_galaxy_props.npy',
                        help='File containing the galaxy properties. A python dictionary is expected '\
                             'as generated by findGalaxyProps.py.')

    parser.add_argument('--cams_to_plot', nargs='+', default=['face','edge','45'],
                        help='Cameras for which to make slice and projection plots ')

    parser.add_argument( '--max_level', default=None, type=int,
                         help='Max level to refine when exporting the oct-tree structure.')

    parser.add_argument('--out_dir',default='sim_dir/analysis/sunrise_analysis/',
                        help='Directory where the output will be placed. A sub directory will be created '\
                            'for each snapshot') 

    parser.add_argument('--rockstar_out_dir', default='sim_dir/analysis/rockstar_output/',
                        help='Directory where to find the rockstar output, used to annotate halos '\
                            'on plots. If not found or set to None halos will not be annotated.')

    parser.add_argument('--no_plots',action='store_true',
                        help='Do not generate projection plots.') 

    parser.add_argument('--no_export',action='store_true',
                        help='Do not export data to fits for Sunrise.') 

    args = vars(parser.parse_args())
    return args


def generate_cameras(normal_vector, distance=100.0, fov=50.0):
    '''
    Set camera positions and orientations
    '''
    from yt.utilities.orientation import Orientation

    print "\nGenerating cameras"
    
    north = np.array([0.,1.,0.])
    orient = Orientation(normal_vector=normal_vector, north_vector=north)
    R=np.linalg.inv(orient.inv_mat)

    camera_set = OrderedDict([
            ['face',([0.,0.,1.],[0.,-1.,0],True)], #up is north=+y
            ['edge',([0.,1.,0.],[0.,0.,-1.],True)],#up is along z
            ['45',([0.,0.7071,0.7071],[0., 0., -1.],True)],
            ['Z-axis',([0.,0.,1.],[0.,-1.,0],False)], #up is north=+y
            ['Y-axis',([0.,1.,0.],[0.,0.,-1.],False)],#up is along z
            ])  
    segments = 10
    np.random.seed(0)
    ts = np.random.random(segments)*np.pi*2
    ps = np.random.random(segments)*np.pi-np.pi/2.0
    for i,(theta, phi) in enumerate(zip(ts,ps)):
        pos = [np.cos(theta),0.,np.sin(phi)]
        vc  = [np.cos(np.pi/2.-theta),0.,np.sin(np.pi/2.-phi)] 
        camera_set['Random_%03i'%(i)]=(pos,vc,False)

    i=0    
    cameras = OrderedDict()
    for name,(normal,north,do_rot)  in camera_set.iteritems():
        orient = Orientation(normal_vector=normal, north_vector=north)
        if do_rot:
            drot = R.copy()
        else:
            drot = np.identity(3)
        sunrise_pos = np.dot(orient.normal_vector, drot)
        sunrise_up  = L.copy()
        if np.all(np.abs(sunrise_up-sunrise_pos)<1e-3):
            sunrise_up[0] *= 0.5 
        sunrise_direction = -1.0*sunrise_pos
        sunrise_afov = 2.0*np.arctan((fov/2.0)/distance)
        norm = lambda x: x/np.sqrt(np.sum(x*x))
        if np.all(np.abs(norm(sunrise_up)-norm(sunrise_pos))<1e-3):
            sunrise_up[0]*=0.5
            sunrise_up = norm(sunrise_up)
        line = (distance*sunrise_pos, distance*sunrise_direction, sunrise_up,
                sunrise_afov, fov, distance) 
        cameras[name] = line
        i+=1

    print "Successfully generated cameras\n"
    return cameras


def write_cameras(prefix, cameras):
    print "Writing cameras to ",  prefix+'.cameras'
    fn = prefix+'.cameras'
    campos = ()
    for name,row in cameras.iteritems():
        campos += (tuple(row[1])+tuple(row[0])+tuple(row[2])+tuple([row[3]]),)
    campos = np.array(campos)
    np.savetxt(fn, campos)   
    fn = prefix+'.camnames'
    fh = open(fn,'w')
    fh.write('\n'.join([c for c in cameras.keys()]))
    fh.close()


def plot_particles(prefix, ds, center, cameras, 
                   dm_particles='darkmatter', star_particles='stars',
                   halo_file=None, cams_to_plot=['face','edge','45']):
    """
    Project the stars and DM densities on the fov of cams_to_plot. 
    Also make slice plots of the density and the LOS velocity of the stars.
    Finally project on 200 kpc and 1 Mpc scales the dark matter density 
    along the normal axis of the first camera on cams_to_plot, circling
    halos.

    """
    camnames = cameras.keys()
    cams = [cameras[n] for n in camnames]

    if halo_file:
        try:
            halo_list = RockstarHaloList(ds, halo_file) 
        except TypeError:
            halo_list = None
    else:
        halo_list = None

    for i,(name,cam) in enumerate(cameras.iteritems()):               
        if name in cams_to_plot:
            offaxisprojection(prefix+'_%s_fov'%name, ds, cam, center, 
                              particle_type=star_particles,
                              halo_list=halo_list)    
            offaxisprojection(prefix+'_%s_fov'%name, ds, cam, center,
                              particle_type=star_particles,
                              halo_list=halo_list, slice=True)
            offaxisprojection(prefix+'_%s_fov'%name, ds, cam, center,
                              particle_type=star_particles, 
                              halo_list=halo_list,
                              field=star_particles+'_Vlos', slice=True)
            offaxisprojection(prefix+'_%s_fov'%name, ds, cam, center,
                              particle_type=dm_particles,
                              halo_list=halo_list)

    offaxisprojection(prefix+'_'+camnames[0]+'_fov', ds, cams[0], center, 
                      particle_type=dm_particles,
                      halo_list=halo_list)
    offaxisprojection(prefix+'_'+camnames[0]+'_200Kpc', ds, cams[0], center,
                      particle_type=dm_particles,
                      fov=200, halo_list=halo_list, min_mass=1e9)
    offaxisprojection(prefix+'_'+camnames[0]+'_2Mpc', ds, cams[0], center,
                      particle_type=dm_particles,
                      fov=ds.arr(2.0, 'Mpc').in_units('kpc'),
                      halo_list=halo_list, min_mass=1e10)
   

def offaxisprojection(prefix, ds, camera, center, field='cic',
                      particle_type='darkmatter', fov=None,
                      slice=False, halo_list=None, min_mass=1e8):

    center = center.in_units('kpc')
    normal, distance, up, width = get_camprops(camera)
    if fov is not None:
        width = fov
    width = ds.arr(width, 'kpc')  
    distance = ds.arr(distance, 'kpc') 

    LeftEdge = center - max(distance, width)/2.0
    RightEdge =  center +  max(distance, width)/2.0 
    box = ds.box(LeftEdge, RightEdge)

    weight = ('deposit', particle_type+'_cic')
 
    if field == 'cic': 
        field = weight
        weight = None
    elif field == particle_type+'_Vlos': 
        def _Vlos(field, data):
            vx, vy, vz = (data[(particle_type+'_cic_velocity_%s'%ax)].in_units('km/s') for ax in 'xyz')
            norm = normal/np.sqrt(np.dot(normal, normal))
            vlos = vx*norm[0] + vy*norm[1] + vz*norm[2]
            return vlos 
        box.ds.add_field(('deposit', particle_type+"_Vlos"), function=_Vlos, units='km/s',
                         take_log=False, force_override=True, 
                         display_name = r'$\rm{%s\ LOS\ Velocity\ (km/s)}$'%particle_type)

    if slice:
        p=yt.OffAxisSlicePlot(box.ds, normal, field, 
                              center.in_units('code_length'),
                              width=(width.value, width.units),
                              north_vector=up)
    else:
        p=yt.OffAxisProjectionPlot(box.ds, normal, field, 
                                   center.in_units('code_length'),
                                   width=(width.value, width.units),
                                   depth=(2*distance.value, distance.units),
                                   north_vector=up, 
                                   weight_field=weight)
    if halo_list:
        p.annotate_hop_circles(halo_list, annotate=True, 
                               fixed_radius=(1.0, 'kpc'),
                               max_number=int(1e9),
                               min_mass=min_mass, min_size=1e-99,
                               width=(2*width.value, width.units))

    del(box)    
    if yt.is_root():
        p.save(prefix)


def plot_gas(prefix, ds, center, cameras, cams_to_plot=['face','edge','45']):
    """
    Make projection and slice plots of the gas density. Also make a couple of
    phase plots.
    """

    field = 'density'

    for i, (name, cam) in enumerate(cameras.iteritems()):
        if name in cams_to_plot:
            normal, distance, up, width = get_camprops(cam)
            p=yt.OffAxisProjectionPlot(ds, normal, field, 
                                       center.in_units('code_length'), 
                                       width=(width, 'kpc'),
                                       depth=(2*distance, 'kpc'),
                                       north_vector=up)
            if yt.is_root():
                p.save(prefix+'_%s_fov'%name)

            p=yt.OffAxisSlicePlot(ds, normal, 'metal_ia_density', center, 
                                  width=(width, 'kpc'), north_vector=up)
            if yt.is_root():
                p.save(prefix+'_%s_fov'%name)

    sph = ds.sphere(center, (20.0, 'kpc'))

    def _MetalMass(field, data):
        return (data['metal_ia_density']*data['cell_volume']).in_units('Msun')
    sph.ds.add_field(('gas', 'MetalMass'), function=_MetalMass, units='Msun')         
    
    def _Temperature(field, data):
        te = data['thermal_energy']
        hd = data['H_nuclei_density']
        temp = (2.0*te/(3.0*hd*yt.physical_constants.kb)).in_units('K')
        return temp
    sph.ds.add_field(('gas', 'Temperature'), function=_Temperature, units='K')

    
    p = yt.PhasePlot(sph, "density", "Temperature", "cell_mass", weight_field=None)
    p.set_unit('density', 'Msun/kpc**3')
    p.set_unit('cell_mass', 'Msun')
    if yt.is_root:
        p.save(prefix)
    
    p = yt.PhasePlot(sph, "density", "Temperature", "MetalMass", weight_field=None)
    p.set_unit('density', 'Msun/kpc**3')
    if yt.is_root:
        p.save(prefix)
    

def export_fits(ds, center, export_radius, prefix, star_particles, max_level=None):
    '''
    Convert the contents of a dataset to a FITS file format that Sunrise
    understands.
    '''

    print "\nExporting data in %s to FITS for Sunrise"%ds.parameter_filename.split('/')[-1]

    filename = prefix+'.fits'
    center = center.in_units('kpc')
    width = export_radius.in_units('kpc')
    info = {}

    fle, fre, ile, ire, nrefined, nleafs, nstars = \
        sunrise_octree_exporter.export_to_sunrise(ds, filename, star_particles, 
                                                  center, width, max_level=max_level)
    info['export_center']=center.value
    info['export_radius']=width.value
    info['export_ile']=ile
    info['export_ire']=ire
    info['export_fle']=fle
    info['export_fre']=fre
    info['export_max_level']=max_level
    info['export_nrefined']=nrefined
    info['export_nleafs']=nleafs
    info['export_nstars']=nstars

    print "Successfully generated FITS for snapshot %s"%ds.parameter_filename.split('/')[-1]
    print info,'\n'
    return info


def get_camprops(cam):
    '''
    Get the properties of the given camera 
    '''
    normal = np.array(cam[0])
    pos = np.array(cam[1])
    up = np.array(cam[2])
    afov = np.array(cam[3])
    distance=np.sqrt(np.sum(pos**2.0))
    fov= cam[4]
    return normal,distance,up,fov


def get_halo_file(rockstar_out_dir, scale):
    '''
    Find the rockstar out_?.list file with the given scale factor
    '''
    halo_file = None
    if not os.path.exists(rockstar_out_dir): pass
    else:
        halo_files = glob(rockstar_out_dir+'/*.list')
        for this_file in halo_files:
            f = file(this_file)
            lines = f.readlines()
            line = [line for line in lines if '#a =' in line]
            scale_in_file = round(float(line[0].split('=')[-1]), 4)
            if scale == scale_in_file:
                halo_file = this_file
                break
    return halo_file


if __name__ == "__main__":

    args = parse()

    import yt
    from yt.analysis_modules.sunrise_export import sunrise_octree_exporter
    from yt.analysis_modules.halo_finding.halo_objects import RockstarHaloList  
    yt.enable_parallelism()
    
    if yt.is_root():
        print '/nStarting '+ sys.argv[0]
        print 'Parsed arguments: '
        print args
        print

    # Get parsed values
    sim_dirs, snap_base = args['sim_dirs'], args['snap_base']
    print 'Analyzing ', sim_dirs

    out_dir = args['out_dir']
    modify_outdir = 0
    if  'sim_dir' in out_dir: 
        out_dir = out_dir.replace('sim_dir','')
        modify_outdir = 1 

    rockstar_out_dir = args['rockstar_out_dir']
    modify_rockstar_outdir = 0 
    if  'sim_dir' in rockstar_out_dir: 
        rockstar_out_dir =  rockstar_out_dir.replace('sim_dir','')
        modify_rockstar_outdir = 1 

    galprops_file = args['galprops_file']
    modify_galprops_file = 0
    if  'sim_dir' in galprops_file: 
        galprops_file = galprops_file.replace('sim_dir','')
        modify_galprops_file = 1

    cam_dist, cam_fov = float(args['distance']), float(args['fov'])  
    star_particles, dm_particles = args['star_particles'], args['dm_particles']
    cams_to_plot = args['cams_to_plot']
    max_level = args['max_level']
    no_plots, no_export = args['no_plots'], args['no_export']
       
    # Loop over simulation directories    
    for sim_dir in sim_dirs:
        
        # Set paths and file names
        sim_dir = os.path.expandvars(sim_dir)
        sim_dir = os.path.abspath(sim_dir)

        if modify_outdir:  out_dir = sim_dir+'/'+out_dir
        if modify_rockstar_outdir:  rockstar_out_dir = sim_dir+'/'+rockstar_out_dir
        if modify_galprops_file: galprops_file = sim_dir+'/'+galprops_file
        
        if yt.is_root():
            if not os.path.exists(out_dir): os.makedirs(out_dir)

        # Get the galaxy properties
        galprops_files = glob(galprops_file)
        if len(galprops_files) > 1:
            print 'More than one file matches %s, '\
                'the supplied file name for the galaxy properties. '\
                'Set which file you want to use with --galprops_file'\
                % (galprops_file)
            sys.exit()
        else:
            galprops_file = galprops_files[0]
            galprops = np.load(galprops_file)[()] 

        # Generate data series
        snaps = glob(sim_dir+'/'+snap_base+'*')
        ts = yt.DatasetSeries(snaps)

        # Loop over snapshot to generate cameras and projection plots, 
        # parallelization happens while generating the plots.
        for ds in reversed(ts): 

            scale = round(1.0/(ds.current_redshift+1.0),4)
            if scale not in galprops['scale']: continue

            idx = np.argwhere(galprops['scale'] == scale)[0][0]

            # Set camera positions and orientations
            stars_L = galprops['stars_L'][idx]
            gas_L = galprops['gas_L'][idx]
            try:
                L_sum = stars_L + gas_L
            except TypeError:
                L_sum = gas_L
            L = L_sum/np.sqrt(np.sum(L_sum*L_sum))
            cameras = generate_cameras(L, distance=cam_dist, fov=cam_fov)    
    
            # Write cameras to file
            scale_dir = out_dir+'a'+str(scale)+'/'
            prefix = galprops_file.replace('galaxy_props.npy', 'a'+str(scale)).split('/')[-1] 
            if yt.is_root():
                if not os.path.exists(scale_dir): os.makedirs(scale_dir)
                write_cameras(scale_dir+prefix, cameras)

            # Make plots
            if not no_plots:
                print "\nGenerating plots for snapshot %s"%ds.parameter_filename.split('/')[-1]
                plots_dir = scale_dir+'/yt_plots/'
                if yt.is_root():
                    if not os.path.exists(plots_dir): os.makedirs(plots_dir)

                gal_center = galprops['stars_hist_center'][idx]
                if np.all(gal_center):
                    gal_center = ds.arr(gal_center, 'kpc')
                else:
                    continue
    
                halo_file = get_halo_file(rockstar_out_dir, scale)
                plot_particles(plots_dir+prefix, ds, gal_center, cameras, 
                               dm_particles=dm_particles, 
                               star_particles=star_particles,
                               halo_file=halo_file,
                               cams_to_plot=cams_to_plot)
                plot_gas(plots_dir+prefix, ds, gal_center, cameras, cams_to_plot=cams_to_plot)
                print "Successfully generated plots for snapshot %s\n"%ds.parameter_filename.split('/')[-1]


        # Send one snapshots to each processor to export 
        if not no_export:
            for ds in ts.piter():
#            for ds in reversed(ts): # comment line above and uncomment this to debug on latest snaps

                scale = round(1.0/(ds.current_redshift+1.0),4)
                if scale not in galprops['scale']: continue

                idx = np.argwhere(galprops['scale'] == scale)[0][0]
                scale_dir = out_dir+'a'+str(scale)+'/'
                prefix = galprops_file.replace('galaxy_props.npy', 'a'+str(scale)).split('/')[-1] 

                # Export fits files with the data for Sunrise
                gal_center = galprops['stars_hist_center'][idx]
                if np.all(gal_center):
                    gal_center = ds.arr(gal_center, 'kpc')
                else:
                    continue
               
                export_radius = ds.arr(max(1.2*cam_dist, 1.2*cam_fov), 'kpc')
                export_info = export_fits(ds, gal_center, export_radius, 
                                          scale_dir+prefix, star_particles, 
                                          max_level=max_level)
                export_info['sim_name'] = prefix.split('_')[0]
                export_info['scale'] = scale
                export_info['halo_id'] = prefix.split('_')[1].replace('halo','')
                np.save(scale_dir+prefix+'_export_info.npy', export_info)
