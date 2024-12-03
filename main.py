import variables_conditions as vc

#  set brho
j = 0.0
dp = 0.1*j
brho = vc.brho0*(1+dp/100)

# set dispersion [m/%] in the ring
dispersion = -3.3494*dp**6 + 2.2088*dp**5 + 1.5699*dp**4 - \
       0.731*dp**3 - 0.1573*dp**2 + 0.0121*dp + 7.0737