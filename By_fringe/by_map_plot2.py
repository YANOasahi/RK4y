import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

data = np.genfromtxt('./by_map.dat')

x_axis = np.array(data[:, 0])
y_axis = np.array(data[:, 1])
z_axis = np.array(data[:, 2])
plots = np.array([data[:, 0],data[:, 1],data[:, 2]])
by = np.array(data[:, 3])

# **********   3 variables fitting
# for Y vs By
# fifth order function
def fifth(x,a,b,c,d,e,f):
    return a*np.power(x,5)+b*np.power(x,4)+c*np.power(x,3)+d*np.square(x)+e*x+f

# for X vs By
# log normal finction
def lognormal(x, A, mu, sigma):
    return A / (x * sigma * np.sqrt(2 * np.pi)) * np.exp(- (np.log(x) - mu)**2 / (2 * sigma**2))

# gaussian
def gaussian(x, A, mu, sigma):
    return A * np.exp(- (x - mu)**2 / (2 * sigma**2))

# for Z vs By
# z looks like a factor
def factor(x,a):
    return a*x

# define fitting function
def combined_model(plots,
                   A_log, mu_log, sigma_log,
                   A1, mu1, sigma1,
                   A2, mu2, sigma2,
                   A3, mu3, sigma3,
                   A4, mu4, sigma4,
                   a_y, b_y, c_y, d_y, e_y, f_y,
                   a_z):
    x, y, z = plots
    lognorm = lognormal(x, A_log, mu_log, sigma_log)
    gauss1 = gaussian(x, A1, mu1, sigma1)
    gauss2 = gaussian(x, A2, mu2, sigma2)
    gauss3 = gaussian(x, A3, mu3, sigma3)
    gauss4 = gaussian(x, A4, mu4, sigma4)
    fifth_y = fifth(y, a_y, b_y, c_y, d_y, e_y, f_y)
    factor_z = factor(z,a_z)
    return factor_z * ((lognorm+gauss1+gauss2+gauss3+gauss4) * fifth_y)

# initial parameter
p0 = [0.0006958040766977115, 6.292789011773518, 0.09067786199274608,  # log normal
      -8.169333763413994e-07, 560.6826636407911, 27.58503262197233,   # gaussian
      2.54948402661953e-07, 656.3964871401387, 32.1459215691162,      # gaussian
      1.9893328877436285e-07, 703.2723248019659, 59.25730588891589,   # gaussian
      1.2306205278883204e-07, 1061.4207387188872, 405.3453929882518,  # gaussian
      1.08879695e-14, 1.92593466e-14, 4.73198961e-10, -4.69999453e-10, 3.80312597e-06, 5.43736144e-07,  # fifth_y
      29891.69217]

try:
    popt, pcov = curve_fit(combined_model, plots, by, p0=p0, maxfev=100000, xtol=1, ftol=1)
except RuntimeError:
    print("Optimal parameters not found. Returning the best found parameters.")

perr = np.sqrt(np.diag(pcov))

A_log, mu_log, sigma_log, \
A1, mu1, sigma1, \
A2, mu2, sigma2, \
A3, mu3, sigma3, \
A4, mu4, sigma4, \
a_y, b_y, c_y, d_y, e_y, f_y, \
a_z = popt

fig2D_2 = plt.figure()
y_by = fig2D_2.add_subplot()
y_by.scatter(y_axis, by, s=2, label='data')
y_by.plot(y_axis, fifth(y_axis,a_y,b_y,c_y,d_y,e_y,f_y), label='fit')

fig2D_3 = plt.figure()
z_by = fig2D_3.add_subplot()
z_by.scatter(z_axis, by, s=2, label='data')
z_by.plot(z_axis, factor(z_axis,a_z), label='fit')