# this is an old file using ROMS flow fields and attempting to adjust the release location of the particles

# this file attempts to use FieldSet.from_c_grid which is supposed to work with the new version of parcels but it doesn't seem to recognize it even though I have Parcels3.1.2

import numpy as np

from parcels import FieldSet, Field
from parcels import ParticleSet
from parcels import Variable
from parcels import ScipyParticle
# from parcels import ErrorCode
import parcels
print(parcels.__version__)

import time as ostime
import sys
ti0 = ostime.time()

import scratchtunaFADpreykernels as pk

import xarray as xr

# Load ROMS dataset
roms_data = xr.open_dataset('/server/hpc/lol_scratch/Leah/salmonira/flow_fields/stcc_his_osu_subset_0001.nc')

# print(roms_data)

# Extract surface-level velocity fields (assuming s_rho=0 is the surface)
u_surface = roms_data['u'].isel(s_rho=0).transpose('ocean_time', 'eta_u', 'xi_u')  # U component at the surface
v_surface = roms_data['v'].isel(s_rho=0).transpose('ocean_time', 'eta_v', 'xi_v')  # V component at the surface

# chops u_surface and v_surface to be same shape 
# corrected_shape = tuple(np.minimum(np.shape(u_surface), np.shape(v_surface)))
# print(corrected_shape)

# Extract the necessary latitude, longitude, and time information
lon_u = roms_data['lon_u'].values  # Longitude for U velocity
lat_u = roms_data['lat_u'].values  # Latitude for U velocity
lon_v = roms_data['lon_v'].values  # Longitude for V velocity
lat_v = roms_data['lat_v'].values  # Latitude for V velocity
time = roms_data['ocean_time'].values  # Time values from ROMS

print(u_surface)
print(np.shape(u_surface))

print(v_surface)
print(np.shape(v_surface))

# Convert time from datetime64 to seconds since the first time point
time_origin = np.datetime64('1970-01-01')  # Use an epoch reference if needed
time_seconds = (roms_data['ocean_time'].values - time_origin) / np.timedelta64(1, 's')

# def create_preyfieldRW(lx, ly, res, nprey=100000, Pavg = float(sys.argv[9])):
#     # Randomly distribute the prey over the grid
#     # Pavg is the average number of prey per grid cell
#     dataP = np.zeros(((ly//res)+2, (lx//res)+2))
#     gc = (ly/res)*(lx/res) 
#     add = 1 * gc / nprey * Pavg
#     for n in range(nprey):
#         i = np.random.randint(1,dataP.shape[0]-1)
#         j = np.random.randint(1,dataP.shape[1]-1)
#         dataP[i,j] += add
#     assert dataP.max() <= 1
#     assert dataP.min() >= 0
#     return dataP

# def create_preyfieldDG(lx, ly, res, nprey=int(1e5), Pavg = float(sys.argv[9])):
#     # Randomly distribute the prey over the grid
#     dataP = np.zeros(((ly//res)+2, (lx//res)+2))
#     gc = (ly/res)*(lx/res) 
#     add = 1 * gc / nprey * Pavg
#     for n in range(nprey):
#         bo = True
#         co = 0
#         while(bo):
#             co += 1
#             i = 1 + np.random.binomial(dataP.shape[0]-2, 0.5)
#             j = 1 + np.random.binomial(dataP.shape[1]-2, 0.3)
#             if(dataP[i,j]<=1-add):
#                 dataP[i,j] += add
#                 bo = False
#             elif(co==10):
#                 bo = False
#     # normalize the field
#     assert dataP.max() <= 1
#     assert dataP.min() >= 0
#     return dataP

# def create_preyfieldBJ(lx, ly, res, nprey=int(1e5), Pavg = float(sys.argv[9])):
#     # Randomly distribute the prey over the grid
#     dataP = np.zeros(((ly//res)+2, (lx//res)+2))
#     gc = (ly/res)*(lx/res) 
#     add = 1 * gc / nprey * Pavg
#     for n in range(nprey):
#         i = 1 + np.random.binomial(dataP.shape[0]-2, 0.5)
#         j = np.random.randint(1,dataP.shape[1]-1)
#         dataP[i,j] += add
#     assert dataP.max() <= 1
#     assert dataP.min() >= 0
#     return dataP

if(__name__=='__main__'):
    # Calculate lx and ly based on the flow field
    lx = lon_u.max() - lon_u.min()  # Length of the flow field in the longitudinal direction
    ly = lat_u.max() - lat_u.min()  # Width of the flow field in the latitudinal direction
    print(f"Calculated habitat length (lx): {lx}")
    print(f"Calculated habitat width (ly): {ly}")
    # assert lx%2==0
    # assert ly%2==0, 'if prey binomial distribution'
    nfad = int(sys.argv[2]) # number of FADs
    ntuna = int(sys.argv[3]) # number of tuna

    npart = ntuna + nfad + 1 # total number of particles
    seed = int(sys.argv[1])
    np.random.seed(seed) # seeding of Monte Carlo simulations

    min_release_lat = 20.5
    min_release_lon = 128.2
    max_release_lat = 21
    max_release_lon = 128.8

    # Set the initial locations of the particles
    X = np.random.uniform(min_release_lon, max_release_lon, npart)
    Y = np.random.uniform(min_release_lat, max_release_lat, npart)

    # Set the initial locations of the particles (at the surface, i.e., depth=0)
    lon_start = np.random.uniform(min_release_lon, max_release_lon, npart)  # Use ROMS longitude bounds
    lat_start = np.random.uniform(min_release_lat, max_release_lat, npart)  # Use ROMS latitude bounds

    # assert (X<=lx).all()
    # assert (Y<=ly).all()
    
    # Flow field configuration

    # print("Shape of u_surface:", u_surface.shape)
    # print("Shape of v_surface:", v_surface.shape)

    # # Also check the shape of the coordinate arrays
    # print("Shape of lon_u:", lon_u.shape)
    # print("Shape of lat_v:", lat_v.shape)
    # print("Shape of time:", time.shape)
    # Create the FieldSet using the consistent or common lon/lat
    # Ensure the FieldSet has the correct shape and transpose the data if needed
    # Create the FieldSet using separate grids for U and V components

    # Create separate Fields for U and V
#     field_U = Field('U', u_surface.values, lon=lon_u, lat=lat_u, time=time_seconds, mesh='spherical', allow_time_extrapolation=True)
#     field_V = Field('V', v_surface.values, lon=lon_v, lat=lat_v, time=time_seconds, mesh='spherical', allow_time_extrapolation=True)

#     # Create the FieldSet and add the U and V fields
#     fieldset = FieldSet(field_U, field_V)

    fieldset = FieldSet.from_c_grid(u_data=u_surface.values, v_data=v_surface.values,
                                lon_u=lon_u, lat_u=lat_u,
                                lon_v=lon_v, lat_v=lat_v,
                                time=time_seconds,
                                mesh='spherical')

    # [print(np.shape(field.__dict__.get("data"))) for field in fieldset.get_fields()]

    # add other fields like temperature, salinity, or prey fields

    # BJ: Bickley Jet
    # DG: Double Gyre
    # RW: Random Walk
    # ff = 'BJ'
    # assert ff in ['RW', 'DG', 'BJ']

    # define the particle types: tuna particle is 0, dFAD particle is 1
    ptype = np.zeros(npart)
    # The zeroth particle is only used in fishing strategy FS1.
    # This is article does nothing, but is only located at a random
    # tuna particle before a fishing event, where it acts as a dFAD.
    ptype[:nfad + 1] = 1

    # Define a fieldset without flow

    # Define res as the largest divisor of both lx and ly, with a maximum threshold
    def calculate_res(lx, ly, max_res=10):
        # Start from the maximum resolution (max_res) and find a divisor
        for r in range(max_res, 0, -1):
            if lx % r == 0 and ly % r == 0:
                return r
        return 1  # Fallback if no divisors are found
    
    res = calculate_res(lx, ly, max_res=10)
    print(f"Calculated resolution (res): {res}")

    # res = 10 # resolution of the field
    # if(ff=='RW'):
    #     dataP = create_preyfieldRW(lx, ly, res)
    # if(ff=='DG'):
    #     dataP = create_preyfieldDG(lx, ly, res)
    # if(ff=='BJ'):
    #     dataP = create_preyfieldBJ(lx, ly, res)
    # gridx, gridy = np.meshgrid(np.arange(-res,lx+res,res), np.arange(-res,ly+res,res))
    # gridx = np.array(gridx) + 0.5*res
    # gridy = np.array(gridy) + 0.5*res

    # sets synthetic zero-velocity flow field
    # fieldset = FieldSet.from_data({'U': np.zeros(dataP.shape), 'V': np.zeros(dataP.shape)},
    #                                {'lon': gridx, 'lat': gridy},
    #                                 mesh='flat')

    # to determine the strength of displacement due to prey field
    fieldset.add_constant('gres',res)
    # fieldset.add_constant('flowtype',ff)
    # Create the field of tuna prey
    # assert ly%res==0
    # assert lx%res==0

    # Adjust prey_data to match the spatial dimensions of lon_u and lat_u
    prey_data = np.ones(((len(time_seconds), len(lat_u), len(lon_u))), dtype=np.float32)  # Uniform prey field with correct shape
    # prey_data = np.random.rand(len(time_seconds), len(lat_u), len(lon_u)) # for a random prey field
    # print(prey_data.shape)
    # print(len(lon_u), len(lat_u), len(time_seconds))

    # Create the prey field with the corrected shape
    field_prey = Field('prey', prey_data, lon=lon_u, lat=lat_u, time=time_seconds, mesh='spherical', allow_time_extrapolation=True)
    fieldset.add_field(field_prey)

    # fieldP = Field('prey', dataP, grid=fieldset.U.grid,
    #               interp_method='nearest', mesh='flat')
    # fieldset.add_field(fieldP) # prey field added to the velocity FieldSet
    # fieldset.prey.to_write = False # enabling the writing of Field prey during execution

    if(nfad>0):
        # Add lists (which are added as an interactive field here)
        # These keep track of the FAD order from FADs with many associated 
        # tuna to FADs with little associated tuna
        # only needed when p>0 in the fishing strategy
        fieldF = Field('FADorders', np.arange(nfad), lon=np.arange(nfad), lat=np.array([0]), time=np.array([0]),
                       interp_method='nearest', mesh='spherical', allow_time_extrapolation=True)
        fieldset.add_field(fieldF) # prey field added to the velocity FieldSet
        fieldset.FADorders.to_write = False # enabling the writing of Field prey during execution
        # FAD number where fish is caught
        fieldF = Field('FADc', np.array([0]), lon=np.array([0]), lat=np.array([0]), time=np.array([0]),
                       interp_method='nearest', mesh='spherical', allow_time_extrapolation=True)
        fieldset.add_field(fieldF) # prey field added to the velocity FieldSet
        fieldset.FADc.to_write = False # enabling the writing of Field prey during execution

    # list that determines at which tuna particle to fish
    # under fishing strategy FS1
    fieldFe = Field('fe', np.array([0]), lon=np.array([0]), lat=np.array([0]), time=np.array([0]),
                       interp_method='nearest', mesh='spherical', allow_time_extrapolation=True)
    fieldset.add_field(fieldFe) # prey field added to the velocity FieldSet
    fieldset.FADc.to_write = False # enabling the writing of Field prey during execution

    # Set the parameters for the model:
    # general
    #  Taxis coefficients
    fieldset.add_constant("kappaT", float(sys.argv[4]))
    fieldset.add_constant("kappaF", float(sys.argv[5]))
    fieldset.add_constant("kappaP", float(sys.argv[6]))
    fieldset.add_constant("kappaI", float(sys.argv[7]))
    max_interaction_distance = 10
    print('realistic FAD-tuna interaction distance is around 10km (7Nm), now (km):',max_interaction_distance)
    fieldset.add_constant("RtF", 2.) # FAD association radius (km)
    fieldset.add_constant("Rtt", 3.) # tuna-tuna max interaction distance (km)
    scale = 300
    fieldset.add_constant("epsP", 12/(24*3600)/scale) # prey depletion by tuna (per second)
    fieldset.add_constant("Td", 2/(24*3600)/scale) # tuna gastric evacuation rate
    fieldset.add_constant("scaleD", scale) # scale tuna gastric evacuation rate

    fieldset.add_constant("epsT", 0.5) # Fraction of associated tuna caught
    p = float(sys.argv[8])
    fieldset.add_constant("p",  p) # p parameter of the geometric distribution
    fieldset.add_constant("nfad", nfad) # total number of FADs
    fieldset.add_constant("ntuna", ntuna) # total number of tuna particles
    # Set a maximum tuna swimming velocity
    fieldset.add_constant("Vmax", (0.4 / 1000)) # km/s
    # the domain
    fieldset.add_constant("Lx", lx)
    fieldset.add_constant("Ly", ly)
    # Determines concentration parameter of von Mises'
    fieldset.add_constant("alpha", 3.) # swimming towards prey
    fieldset.add_constant("gamma", 2.) # swimming towards other tuna

    # Random walk flow:
    # if(ff=='RW'):
    #     fieldset.add_constant_field("Kh_zonalF", 0.05/1000, mesh="flat") # in km/s
    #     fieldset.add_constant_field("Kh_meridionalF", 0.05/1000, mesh="flat") # in km/s
    #     fieldset.add_constant_field("Kh_zonalT", 0.1/1000, mesh="flat") 
    #     fieldset.add_constant_field("Kh_meridionalT", 0.05/1000, mesh="flat") 

    # Parameters of the Logistic curve, which determines
    # the dependence of FAD attraction strength on the number 
    # of associated tuna
    fieldset.add_constant("lL", 1.) # maximum of logistic curve
    fieldset.add_constant("lk", 0.35) # steepness of logistic curve
    fieldset.add_constant("lx0", 12) # value of the sigmoid midpoint
    # And for the vp
    fieldset.add_constant("pL", 2.5) # maximum of logistic curve

    # Parameter for the Bickley Jet flow
    # if(ff=='BJ'):
    #     fieldset.add_constant("Ubj", (.1 / 1000)) # maximum flow strength (km/s)

    # Parameters for the double gyre flow
    # if(ff=='DG'):
    #     fieldset.add_constant("A", (0.05 / 1000)) # flow strength
    #     fieldset.add_constant("omega", 2*np.pi/ (10*24*60*60)) # frequency of one oscillation (per second)
    #     fieldset.add_constant("epsDG", 0.2) # 

    # Create custom particle class with extra variable that indicates
    # whether the interaction kernel should be executed on this particle.
    class TFParticle(ScipyParticle):
        ptype = Variable('ptype', dtype=np.int32, to_write='once')
        caught = Variable('caught', dtype=np.float32, initial=0)
        mlon = Variable('mlon', dtype=np.float32, to_write=False, initial=0)
        mlat = Variable('mlat', dtype=np.float32, to_write=False, initial=0)
        FADkap = Variable('FADkap', dtype=np.float32, to_write=True, initial=1.)
        # To govern the displacement of particles (used for velocity normalization):
        dla = Variable('dla', dtype=np.float32, to_write=False, initial=0.)
        dlo = Variable('dlo', dtype=np.float32, to_write=False, initial=0.)
        gradx = Variable('gradx', dtype=np.float32, to_write=False, initial=0.)
        grady = Variable('grady', dtype=np.float32, to_write=False, initial=0.)
        # Stomach fullness:
        St = Variable('St', dtype=np.float32, to_write=False, initial=0.5)
        Sta = Variable('Sta', dtype=np.float32, to_write=True, initial=0)
        Stna = Variable('Stna', dtype=np.float32, to_write=True, initial=0)
        Stac = Variable('Stac', dtype=np.float32, to_write=True, initial=0)
        Stnac = Variable('Stnac', dtype=np.float32, to_write=True, initial=0)

    print('number of FADs: ',np.sum((ptype==1)))
    # pset = ParticleSet(fieldset=fieldset, pclass=TFParticle,
    #                    lon=X, lat=Y,
    #                    interaction_distance=max_interaction_distance,
    #                    ptype=ptype.astype(np.int32))

    # Create ParticleSet using the surface ROMS field and surface-level particles
    pset = ParticleSet(fieldset=fieldset, pclass=TFParticle,
                    lon=lon_start, lat=lat_start, interaction_distance=max_interaction_distance,
                    ptype=ptype.astype(np.int32))

    output_file = pset.ParticleFile(
        name="output/RRtunaFADPrey_no%d_npart%d_nfad%d_T%.2f_F%.2f_P%.2f_I%.2f_p%.1f_Pa%.1f.zarr" % (
            seed, npart, nfad, float(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6]), float(sys.argv[7]), float(sys.argv[8]), float(sys.argv[9])
        ),
        outputdt=4.32e4  # output twice a day
    )

    rt = 8.64e5 # 10 days of simulation
    print('model run time (days): ',rt/24/3600)

    # set up the kernels, which depends on the configuration used
    kernels = pset.Kernel(pk.CaughtP) + pset.Kernel(pk.GEvacuation) # increase tuna hunger
    ikernels = pset.InteractionKernel(pk.Iattraction)
    # if(ff=='DG'):
    #     kernels += pset.Kernel(pk.DoubleGyre) # Double Gyre flow 
    # elif(ff=='RW'):
    #     kernels += pset.Kernel(pk.DiffusionUniformKhP) # Random walk flow
    # if(ff=='BJ'):
    #     kernels += pset.Kernel(pk.BickleyJet) # Bickley jet flow
    #     kernels += pset.Kernel(pk.DisplaceParticle) # displace tuna due to swimming
    #     kernels += pset.Kernel(pk.zper_mrefBC) # reflective boundary conditions
    #     kernels += pset.Kernel(pk.PreyGrad_zpb) # calculate prey gradient

    #     ikernels += pset.InteractionKernel(pk.ItunaFAD_zpb)
    #     ikernels += pset.InteractionKernel(pk.Itunatuna_zpb)
    kernels += pset.Kernel(pk.DisplaceParticle) # displace tuna due to swimming
    # kernels += pset.Kernel(pk.reflectiveBC) # reflective boundary conditions
    kernels += pset.Kernel(pk.PreyGrad) # calculate prey gradient

    ikernels += pset.InteractionKernel(pk.ItunaFAD)
    ikernels += pset.InteractionKernel(pk.Itunatuna)

    kernels += pset.Kernel(pk.FaugerasDiffusion)
    kernels += pset.Kernel(pk.Inertia)
    kernels += pset.Kernel(pk.prevloc)
    kernels += pset.Kernel(pk.PreyDeplete)

    ikernels +=  pset.InteractionKernel(pk.Stcheck)
    if(p!=-2): # p==-2 means that no fish is caught
        ikernels +=  pset.InteractionKernel(pk.ItunaPredFAD)

    # def delete_particle(particle, fieldset, time):
    #     particle.delete()

    pset.execute(pyfunc=kernels,
                 pyfunc_inter=ikernels,
                              # 20 minute time step
                 runtime=rt, dt=1.2e3, output_file=output_file, recovery={"DeleteParticle": delete_particle}, verbose_progress=False)


    output_file.close()
print('total time (minutes):',(ostime.time()-ti0)/60)

print(f"FieldSet longitude range: {lon_u.min()} to {lon_u.max()}")
print(f"FieldSet latitude range: {lat_u.min()} to {lat_u.max()}")
