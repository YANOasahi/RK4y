# modify these values to obtain enough isochronous condition

# _______about magnets_______
# current of the main magnets
main_current = 1918.1205

# _______about particles_______
# proton number
z = 32.0
# mass
mass = 77.9229
# brho of particles when dp/p=0 [Tm]
brho0 = 4.7447
# beta of the ring in the X-axis
betax = 7.817

# _______about particles_______
step_time = 0.0001   # in ns

# _______basically we don't have to change below_______
# current of the trim coils
trim_current = [592.128, 1.958, 91.045, 226.784,
                158.111, 146.1, 117.213, 280.012, 47.809, 279.008]

# trim_current = [902.398, -71.626, -91.662, 262.263,
#                 111.998, 256.673, -215.892, 710.165, -382.961, 288.332]

# calculation is normalized by this current
main_base = 1915.0

# set magnetic flux density
B0z = 1.182166

# initial position of particles
x0 = 9287.959673
y0 = 0.0

# amu [ev/c^2]
amu = 931494061

# light speed [m/s]
c = 299792458
