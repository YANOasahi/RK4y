import numpy as np
import variables_conditions as vc

#  set brho
j = 0.0
dp = 0.1*j
brho = vc.brho0*(1+dp/100)

# set emittance
l = 0.0
emi = l*20   # 20 steps
a_max = np.sqrt(emi/vc.betax)*0.95
a_bin = 2*a_max/15
k = 0.0
a_init = k*a_bin - a_max

# set dispersion [m/%] in the ring
dispersion = -3.3494*dp**6 + 2.2088*dp**5 + 1.5699*dp**4 - \
            0.731*dp**3 - 0.1573*dp**2 + 0.0121*dp + 7.0737

# m/q
MoQ=vc.mass/vc.z
# momentum based on brho
mom=vc.c*vc.z*brho*1E-6
# beam energy based on momentum
energy = (np.sqrt(mom**2+(vc.mass*vc.amu)**2)-(vc.mass*vc.amu))/vc.mass
# beta based on brho and q/m
beta = ((vc.c*vc.z*brho)/(vc.amu*vc.mass))/np.sqrt(((vc.c*vc.z*brho)/(vc.amu*vc.mass))*((vc.c*vc.z*brho)/(vc.amu*vc.mass))+1)
# Lorentz factor
gamma=1/np.sqrt(1-beta**2)
# velocity of the beam [mm/ns]
v=vc.c*beta*1E-6

print('*******    Initial conditions    *******')
print(f'Proton number of the beam is {vc.z}')
print(f'Mass of the beam is {vc.mass}')
print(f'Mass to charge ratio of the beam is {MoQ}')
print(f'Beam energy is {energy*1E6} MeV/u')
