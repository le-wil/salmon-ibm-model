# example run: python MOM6salmonira.py 321 10 0 1 0.01 0.1

# model sys.argv:
# 1: random seed
# 2: number of fish (>0)
# 3: kappaF - relative contribution of other fish to fish swimming direction (0-1)
# 4: kappaB - relative contribution of bathymetry field to fish swimming direction (0-1)
# 5: kappaI - relative contribution of inertia to fish swimming direction (0.01)
# 6: p - mean prey index per grid cell (0.1)

# Individual-based model of salmon migration using MOM6 velocity, bathymetry, and temperature fields
# 20 min time step, save every 12 hours, 10 days (all adjustable)

import numpy as np
from parcels import FieldSet, Field, ParticleSet, Variable, ScipyParticle, AdvectionRK4
import time as ostime
import sys
import salmonirakernels as pk
import xarray as xr

ti0 = ostime.time()

# Load MOM6 datasets
mom6_u = xr.open_dataset('/server/hpc/lol_scratch/Leah/salmonira/flow_fields/psl_flow_fields/ssu_rotate.nep.full.hcast.daily.raw.r20241015.199301-201912.nc')

mom6_v = xr.open_dataset('/server/hpc/lol_scratch/Leah/salmonira/flow_fields/psl_flow_fields/ssv_rotate.nep.full.hcast.daily.raw.r20241015.199301-201912.nc')

grid = xr.open_dataset('/server/hpc/lol_scratch/Leah/salmonira/flow_fields/psl_flow_fields/ocean_static.nc')

# DATA EXTRACTION

# all data is surface-level
u_surface = mom6_u['ssu_rotate'].transpose('time', 'yh', 'xh')
v_surface = mom6_v['ssv_rotate'].transpose('time', 'yh', 'xh')

# zonal velocities at u-points
lon_u = grid['geolon_u'].values[:, :-1]
lat_u = grid['geolat_u'].values[:, :-1]

# meridional velocities at v-points
lon_v = grid['geolon_v'].values[:-1, :]
lat_v = grid['geolat_v'].values[:-1, :]

# for prey field (centered at h-points)
lon_c = grid['geolon'].values
lat_c = grid['geolat'].values

# extract time
time = mom6_u['time'].values

# convert time from datetime64 to seconds since an epoch reference
time_origin = np.datetime64('1970-01-01')
time_seconds = (time - time_origin) / np.timedelta64(1, 's')

# extract depth
depth = grid['deptho']

# run the model
if(__name__=='__main__'):
    
    # calculate domain
    lon_min = np.nanmin(lon_u)
    lon_max = np.nanmax(lon_u)
    lat_min = np.nanmin(lat_u)
    lat_max = np.nanmax(lat_u)

    d_lon = lon_max - lon_min
    d_lat = lat_max - lat_min

    mean_lat = 0.5 * (lat_min + lat_max)
    cos_lat = np.cos(np.deg2rad(mean_lat))
    
    # habitat domain
    lx = d_lon * 111 * cos_lat
    ly = d_lat * 111
    
    print('Domain size (km): lx =', lx, ', ly =', ly)

    nfish = int(sys.argv[2]) # number of fish
    npart = nfish # total number of particles
    seed = int(sys.argv[1])
    np.random.seed(seed) # seeding of Monte Carlo simulations

    # constrain to around Columbia River mouth (add to system arguments?)
    release_lon = -124.08
    release_lat = 46.14
    delta = 0.05  # 0.05 degrees spread ~5km

    # randomly disperse particles within release area
    lon_start = np.random.uniform(release_lon - delta, release_lon + delta, npart)
    lat_start = np.random.uniform(release_lat - delta, release_lat + delta, npart)
    
    # Correct lon_start if negative
    lon_start = np.where(lon_start < 0, lon_start + 360, lon_start)

    # manually create fields (adjust if from_mom6 function is created)
    field_U = Field('U', 
                    data=u_surface.values,
                    lon=grid['geolon_u'].values[:, :-1], 
                    lat=grid['geolat_u'].values[:, :-1],
                    time=time_seconds, 
                    mesh='spherical', 
                    allow_time_extrapolation=True)

    field_V = Field('V', 
                    data=v_surface.values,
                    lon=grid['geolon_v'].values[:-1, :], 
                    lat=grid['geolat_v'].values[:-1, :],
                    time=time_seconds, 
                    mesh='spherical', 
                    allow_time_extrapolation=True)

    fieldset = FieldSet(U=field_U, V=field_V)
    
    # create a bathymetry field
    field_depth = Field(name='bathy', 
                        data=depth.values, 
                        lon=depth['xh'].values, 
                        lat=depth['yh'].values, 
                        mesh='spherical', 
                        allow_time_extrapolation=True)

    fieldset.add_field(field_depth)
    
    # ADD TEMPERATURE FIELD

    # add other fields like salinity or prey fields

    # create prey_data with spatial dimensions of lon_u and lat_u
    prey_data = np.ones((len(time_seconds), lat_u.shape[0], lat_u.shape[1]))

    field_prey = Field('prey', 
                       prey_data, 
                       lon=lon_u, 
                       lat=lat_u, 
                       time=time_seconds, 
                       mesh='spherical', 
                       allow_time_extrapolation=True)
    
    fieldset.add_field(field_prey)
    
    # check dimensions
    for field in fieldset.get_fields():
        if hasattr(field, 'data'):
            print(f"Field '{field.name}': data shape = {field.data.shape}")
    print("End of Field Info\n")
    
    # Set the parameters for the model:
    
    # Estimate mean grid spacing (in degrees) for U-grid (lon_u, lat_u)
    # Take central differences along one axis
    dy = np.nanmean(np.diff(lat_u[:, 0]))  # spacing in latitude direction
    dx = np.nanmean(np.diff(lon_u[0, :]))  # spacing in longitude direction

    # average the two to get a scalar resolution
    avg_gres = 0.5 * (dx + dy)

    fieldset.add_constant("gres", avg_gres)

    max_interaction_distance = 10 # km
    
    fieldset.add_constant("kappaF", float(sys.argv[3]))
    fieldset.add_constant("kappaB", float(sys.argv[4]))
    fieldset.add_constant("kappaI", float(sys.argv[5]))

    scale = 300
    fieldset.add_constant("epsP", 12/(24*3600)/scale) # prey depletion by fish (per second)
    fieldset.add_constant("Td", 2/(24*3600)/scale) # fish gastric evacuation rate
    fieldset.add_constant("scaleD", scale) # scale fish gastric evacuation rate

    p = float(sys.argv[6])
    fieldset.add_constant("p",  p) # p parameter of the geometric distribution
    fieldset.add_constant("nfish", nfish) # total number of fish particles
    # Set a maximum fish swimming velocity
    fieldset.add_constant("Vmax", (0.5 / 1000)) # km/s
    # Set a maximum fish-fish interaction distance
    fieldset.add_constant("RtT", 3.0)
    # the domain
    fieldset.add_constant("Lx", lx)
    fieldset.add_constant("Ly", ly)
    # Determines concentration parameter of von Mises'
    fieldset.add_constant("alpha", 3.) # swimming towards prey
    fieldset.add_constant("gamma", 2.) # swimming towards other fish

    # Create custom particle class with extra variable that indicates whether the interaction kernel should be executed on this particle.
    class TFParticle(ScipyParticle):
        mlon = Variable('mlon', dtype=np.float32, to_write=False, initial=0)
        mlat = Variable('mlat', dtype=np.float32, to_write=False, initial=0)
        # To govern the displacement of particles (used for velocity normalization):
        dla = Variable('dla', dtype=np.float32, to_write=False, initial=0.)
        dlo = Variable('dlo', dtype=np.float32, to_write=False, initial=0.)
        gradx_prey = Variable('gradx_prey', dtype=np.float32, to_write=False, initial=0.)
        grady_prey = Variable('grady_prey', dtype=np.float32, to_write=False, initial=0.)
        gradx_bathy = Variable('gradx_bathy', dtype=np.float32, to_write=False, initial=0.)
        grady_bathy = Variable('grady_bathy', dtype=np.float32, to_write=False, initial=0.)
        # Stomach fullness:
        # St = Variable('St', dtype=np.float32, to_write=False, initial=0.5)
        # Sta = Variable('Sta', dtype=np.float32, to_write=True, initial=0)
        # Stna = Variable('Stna', dtype=np.float32, to_write=True, initial=0)
        # Stac = Variable('Stac', dtype=np.float32, to_write=True, initial=0)
        # Stnac = Variable('Stnac', dtype=np.float32, to_write=True, initial=0)

    # Create ParticleSet using the surface field and surface-level particles
    pset = ParticleSet(fieldset=fieldset, pclass=TFParticle,
                    lon=lon_start, lat=lat_start, interaction_distance=max_interaction_distance)

    output_file = pset.ParticleFile(
        name=f"output/MOM6salmonira_no{seed}_nfish{nfish}_kappaF{fieldset.kappaF:.2f}_kappaB{fieldset.kappaB:.2f}_kappaI{fieldset.kappaI:.2f}.zarr",
        outputdt=4.32e4, # output twice a day
    )

    rt = 8.64e5 # 10 days of simulation
    # 1.728e5 # 2-day test run 
    
    # set up the kernels
    # kernels = pset.Kernel(pk.CaughtP) + pset.Kernel(pk.GEvacuation) 
    # ^increase fish hunger (use in future for predation?)
    kernels = AdvectionRK4 + pset.Kernel(pk.prevloc) # find fish previous location
    kernels += pset.Kernel(pk.FaugerasDiffusion) # find swimming angle (bathymetry + northward bias)
    kernels += pset.Kernel(pk.Inertia) # calculate effect of inertia
    kernels += pset.Kernel(pk.DisplaceParticle) # displace fish due to swimming
    kernels += pset.Kernel(pk.reflectiveBC) # reflective boundary conditions
    kernels += pset.Kernel(pk.BathyGrad) # calculate direction of 50m isobath
    
    # ikernels = pset.InteractionKernel(pk.Ifishfish) # calculate effect of other fish on particle
                           
    pset.execute(pyfunc=kernels,            
                 runtime=rt, 
                 dt=1.2e3, # 20 minute time step
                 output_file=output_file,
                 verbose_progress=True)
    
    # add/remove output_file=output_file, for debugging
    # add pyfunc_inter=ikernels, if add ikernels
    
print(f"Simulation complete! {nfish} salmon particles tracked for {rt/86400:.1f} days.")
print(f"Total execution time: {(ostime.time() - ti0)/60:.2f} minutes.")