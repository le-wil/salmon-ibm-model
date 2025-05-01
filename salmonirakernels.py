import numpy as np
import math
import random 
import parcels.rng as ParcelsRandom
import scipy
from scipy.stats import vonmises

__all__ = [# add particle-particle interaction kernels if needed
           # add particle-field interaction kernels if needed
           'prevloc', 'Inertia', 'BathyGrad',
           'FaugerasDiffusion','DisplaceParticle', 
           'reflectiveBC' # boundary conditions
           ]

#%% Interaction kernels -------------------------------------------------------------------------------------------------------------

# def Ifishfish(particle, fieldset, time, neighbors, mutator):
#     """InterActionKernel that "pulls" all neighbor fish particles toward the fish"""
#     try:
#         if(fieldset.kappaF!=0):
#             DS = [0,0]
#             for n in neighbors:
#                 dx = np.array([n.lat - particle.lat, n.lon - particle.lon])
#                 norm = np.linalg.norm(dx)
#                 if(norm>0 and norm<fieldset.RtT): # RtT: salmon-salmon interaction radius
#                     DS[0] += (dx[0] / norm)
#                     DS[1] += (dx[1] / norm)

#             if(DS!=[0,0]):
#                 def f(particle, dlat, dlon, ddepth):
#                     particle.dla += dlat
#                     particle.dlo += dlon
#                 vm = [fieldset.gamma * np.linalg.norm(DS), np.arctan2(DS[1], DS[0])]
#                 angle = ParcelsRandom.vonmisesvariate(vm[1], vm[0])
#                 d_vec = [
#                     fieldset.kappaF * np.cos(angle),
#                     fieldset.kappaF * np.sin(angle),
#                     0
#                 ]

#                 mutator[particle.id].append((f, d_vec))
#     except Exception as e:
#         print(f"Error in Ifishfish kernel: {e}")
#         print(f"Particle state: {particle}")
#     return

# def ItunaPredFAD(particle, fieldset, time, neighbors, mutator):
#     """The predation of tuna near a FAD"""
#     try:
#         # The geometric function
#         def geometric(x, p=fieldset.p):
#             return (1-p)**(x-1)*p

#         # keep track of caught tuna (for both FAD and tuna particles)
#         def fcF(particle, dc):
#             particle.caught += dc
#         def fcN(particle, dc):
#             particle.caught += dc
#         def reset_flo(particle, plon, plat):
#             particle.lon = plon
#             particle.lat = plat

#         # if FAD, consume tuna with a probability
#         day = 86400 # if once a day
#         if(day%particle.dt!=0):
#             print('no fishing taking place!!!!! Set a different dt value.')

#         # Set the fishing location if p==-1
#         if(np.isclose(time%day,day-particle.dt) and fieldset.p==-1):
#             if(particle.id==0):
#                 fieldset.fe.data[0,0,0] = random.randint(fieldset.nfad+1,fieldset.nfad+1+fieldset.ntuna)
#             elif(particle.id==fieldset.fe.data[0,0,0]):
#                 mutator[0].append((reset_flo,[particle.lon, particle.lat]))

#         if(particle.ptype==1 and time>0 and np.isclose(time%day,0)):  # if FAD, consume tuna with a probability
#             if(particle.id==0 and fieldset.p==-1):
#                 prob = 1
#             elif(fieldset.nfad==0):
#                 prob = 0
#             elif(particle.id!=0 and fieldset.p==-1):
#                 prob = 0
#             elif(particle.id==0):
#                 prob = 0
#             else:
#                 fr = fieldset.FADorders.data[0][0].tolist()
#                 nl = [fr.index(x) for x in sorted(fr, reverse=True)[:fieldset.nfad]]
#                 ci = fieldset.FADc.data[0,0,0]-1
#                 if(particle.id==nl[int(ci)]):
#                     prob = 1
#                 else:
#                     prob = 0

#             assert prob in [0,1]
#             if(prob>0):
#                 prob *= fieldset.epsT
#                 for n in neighbors:
#                     if n.ptype==0:
#                         x_p = np.array([particle.lat, particle.lon, particle.depth]) # FAD location
#                         x_n = np.array([n.lat, n.lon, n.depth]) # tuna location
#                         dist = np.linalg.norm(x_p-x_n)
#                         if(dist<fieldset.RtF and (random.random()<prob)):
#                             mutator[particle.id].append((fcF,[1]))  # FAD catches tuna particle
#                             mutator[n.id].append((fcN, [1]))  # tuna particle is caught
#     except Exception as e:
#         print(f"Error in ItunaPredFAD kernel: {e}")
#         print(f"Particle state: {particle}")
#     return

#%% Field-particle interaction kernels ------------------------------------------------------
# def PreyDeplete(particle, fieldset, time):
#     """Depletes prey at particle location based on stomach fullness and redistributes it randomly in domain"""
#     try:
#         xi = (np.abs(np.array(fieldset.prey.lon)-particle.lon)).argmin()
#         yi = (np.abs(np.array(fieldset.prey.lat)-particle.lat)).argmin()
#         preyno = fieldset.prey.data[0,yi,xi]
#         dep = min(fieldset.epsP*particle.dt, preyno)
#         if(1-particle.St>fieldset.scaleD*dep):
#             # increase stomach fullness
#             particle.St = min(1.0, particle.St + fieldset.scaleD * dep)
#             # deplete prey from fieldset
#             fieldset.prey.data[0,yi,xi] -= dep

#             # redistribute prey elsewhere (up to 3 tries)
#             for _ in range(3):
#                 i = np.random.randint(1,fieldset.prey.data[0].shape[1]-1)
#                 j = np.random.randint(1,fieldset.prey.data[0].shape[0]-1)
#                 if fieldset.prey.data[0, j, i]<=(1-dep):
#                     fieldset.prey.data[0, j, i] += dep
#                     break
#             else:
#                 print('did not add the depletion properly')

#         assert fieldset.prey.data[0].min() >=0 
#         assert fieldset.prey.data[0].max() <=1 
#         fieldset.prey.grid.time[0] = time # updating Field prey time
#     except Exception as e:
#         print(f"Error in PreyDeplete kernel: {e}")
#         print(f"Particle state: {particle}")

# def GEvacuation(particle, fieldset, time):
#     """The fish stomach gets emptier over time"""
#     E = fieldset.Td * fieldset.scaleD * particle.dt
#     # particle.St -= min(particle.St,E)
#     particle.St = max(0, particle.St - min(particle.St, E))  # Avoid negative fullness

# def DiffusionUniformKhP(particle, fieldset, time):
#     """Same as the DiffusionUniformKh kernel,
#     but allows different Kh for different particle types"""
#     # Wiener increment with zero mean and std of sqrt(dt)
#     dWx = ParcelsRandom.normalvariate(0, math.sqrt(math.fabs(particle.dt)))
#     dWy = ParcelsRandom.normalvariate(0, math.sqrt(math.fabs(particle.dt)))

#     if(particle.ptype==0):
#         bx = math.sqrt(2 * fieldset.Kh_zonalT[particle])
#         by = math.sqrt(2 * fieldset.Kh_meridionalT[particle])
#     elif(particle.ptype==1):
#         bx = math.sqrt(2 * fieldset.Kh_zonalF[particle])
#         by = math.sqrt(2 * fieldset.Kh_meridionalF[particle])
#     else:
#         bx = 0
#         by = 0
#     assert 1e20>particle.lon, 'beginU  %.3f'%(particle.lon)
#     particle.lon += bx * dWx
#     particle.lat += by * dWy
#     assert 1e20>particle.lon, 'endU'

#%% Other kernels --------------------------------------------------------------------------
def prevloc(particle, fieldset, time):
    """Stores particle's current lon and lat before displacement for inertia-based swimming behavior"""
    particle.mlon = particle.lon
    particle.mlat = particle.lat

def Inertia(particle, fieldset, time):
    """Adds a displacement vector in the same direction as the previous movement scaled by kappaI"""
    if(time>0):
        dlon = particle.lon - particle.mlon
        dlat = particle.lat - particle.mlat
        norm = (dlon**2+dlat**2)**0.5
        if(norm!=0):
            particle.dlo += fieldset.kappaI * dlon / norm
            particle.dla += fieldset.kappaI * dlat / norm
        assert particle.dlo<1e20

# def CaughtP(particle, fieldset, time):
#     # Reposition the tuna to a random location, if it has been caught
#     if(particle.ptype!=1):
#         if(particle.caught>0):
#             particle.lon = random.uniform(0,fieldset.Lx)
#             particle.lat = random.uniform(0,fieldset.Ly)
#             particle.caught=0

# def PreyGrad(particle, fieldset, time):
#     """Calculate the prey gradient at particle location"""
#     # is there a built-in gradient method?
#     # handle edge cases differently? currently set to 0 if can't look left or right
# #     if(particle.lon>fieldset.Lx-fieldset.gres/2):
# #         gradresr = 0
# #         gradresl = fieldset.gres
# #     else:
# #         gradresr = fieldset.gres/2
# #     if(particle.lon<fieldset.gres/2):
# #         gradresl = 0
# #         gradresr = fieldset.gres
# #     else:
# #         gradresl = fieldset.gres/2    
# #     if(particle.lat>fieldset.Ly-fieldset.gres/2):
# #         gradresu = 0
# #         gradresd = fieldset.gres
# #     else:
# #         gradresu = fieldset.gres/2
# #     if(particle.lat<fieldset.gres/2):
# #         gradresd = 0
# #         gradresu = fieldset.gres
# #     else:
# #         gradresd = fieldset.gres/2
# #     particle.gradx = (fieldset.prey[0,0,particle.lat, particle.lon+gradresr] -
# #             fieldset.prey[0,0,particle.lat, particle.lon-gradresl]) 
# #     particle.grady = (fieldset.prey[0,0,particle.lat+gradresu, particle.lon] -
# #             fieldset.prey[0,0,particle.lat-gradresd, particle.lon])
#     xi = np.abs(np.array(fieldset.prey.lon[0]) - particle.lon).argmin()
#     yi = np.abs(np.array(fieldset.prey.lat[:,0]) - particle.lat).argmin()

#     dx = dy = fieldset.gres
    
#     # Apply central differences with minimal edge handling
#     if 1 <= xi < fieldset.prey.data.shape[2] - 1:
#         particle.gradx_prey = (fieldset.prey.data[0, yi, xi + 1] - fieldset.prey.data[0, yi, xi - 1]) / (2 * dx)
#     else:
#         particle.gradx_prey = 0

#     if 1 <= yi < fieldset.prey.data.shape[1] - 1:
#         particle.grady_prey = (fieldset.prey.data[0, yi + 1, xi] - fieldset.prey.data[0, yi - 1, xi]) / (2 * dy)
#     else:
#         particle.grady_prey = 0

def BathyGrad(particle, fieldset, time):
    """Calculate the local direction to the 50m isobath"""
    # Find nearest grid point indices
    xi = np.abs(fieldset.bathy.lon - particle.lon).argmin()
    yi = np.abs(fieldset.bathy.lat - particle.lat).argmin()

    dx = dy = fieldset.gres

    if 1 <= xi < fieldset.bathy.data.shape[1] - 1 and 1 <= yi < fieldset.bathy.data.shape[0] - 1:
        depth_center = fieldset.bathy.data[0, yi, xi] - 50
        depth_xplus = fieldset.bathy.data[0, yi, xi + 1] - 50
        depth_xminus = fieldset.bathy.data[0, yi, xi - 1] - 50
        depth_yplus = fieldset.bathy.data[0, yi + 1, xi] - 50
        depth_yminus = fieldset.bathy.data[0, yi - 1, xi] - 50

        particle.gradx_bathy = (depth_xplus - depth_xminus) / (2 * dx)
        particle.grady_bathy = (depth_yplus - depth_yminus) / (2 * dy)
    else:
        particle.gradx_bathy = 0
        particle.grady_bathy = 0

def FaugerasDiffusion(particle, fieldset, time):
    '''Applies swimming toward 50 m isobath with a northward bias'''
    # def LogisticCurve(x, L=fieldset.pL, k=15, x0=0.7):
    #     # x is the stomach fullness
    #     x = 1 - x  # make it stomach emptiness
    #     res = L / (1 + math.e**(-k * (x - x0)))
    #     return res
    
    # bathymetry-based gradient direction
    mu = np.arctan2(particle.grady_bathy, particle.gradx_bathy)
    
    # northward swimming preference
    north_angle = np.pi / 2
    bias_strength = 0.3  # 0 = pure bathy gradient; 1 = pure north
    mu = (1 - bias_strength) * mu + bias_strength * north_angle
    
    kappaM = fieldset.alpha * np.linalg.norm([
        particle.gradx_bathy * fieldset.gres, 
        particle.grady_bathy * fieldset.gres])  # standard deviation angle
    
    angle = ParcelsRandom.vonmisesvariate(mu, kappaM)
    # particle displacement (can add direction scaled by stomach emptiness by * LogisticCurve(particle.St))
    particle.dlo += fieldset.kappaB * np.cos(angle)
    particle.dla += fieldset.kappaB * np.sin(angle)

def DisplaceParticle(particle, fieldset, time):
    '''Displace the particle according to the summed swimming behaviour'''
    # adjust so prey modulates speed?
    # v = fieldset.Vmax*(1-fieldset.prey[0,0,particle.lat, particle.lon])
    # constant swimming speed
    v = fieldset.Vmax  # in km/s
    norm = (particle.dlo**2 + particle.dla**2)**0.5
    if norm > 0:
        v_lat = v / 111 # km -> deg lat
        v_lon = v / (111 * math.cos(math.radians(particle.lat))) # km -> deg lon
        # Apply small displacement in degrees
        particle.lon += (particle.dlo / norm * v_lon) * particle.dt
        particle.lat += (particle.dla / norm * v_lat) * particle.dt
    # set the displacement per timestep to zero again
    particle.dlo = 0
    particle.dla = 0

#%% boundary conditions ---------------------------------------------------------------------
# reflective in all directions
def reflectiveBC(particle, fieldset, time):
    if(particle.lat>fieldset.Ly):
        particle.lat = fieldset.Ly - (particle.lat-fieldset.Ly)
    elif(particle.lat<0):
        particle.lat *= -1
    if(particle.lon>fieldset.Lx):
        particle.lon = fieldset.Lx - (particle.lon-fieldset.Lx)
    elif(particle.lon<0):
        particle.lon *= -1
        
# add deleteparticle behavior?
# add preystarvationmortality kernel?