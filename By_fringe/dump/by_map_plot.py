import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# data = np.genfromtxt('./by_map.dat')
data = np.genfromtxt('./by_map_z25.dat')

x_axis = np.array(data[:, 0])
y_axis = np.array(data[:, 1])
z_axis = np.array(data[:, 2])
by = np.array(data[:, 3])
plots = np.array([data[:, 0], data[:, 1], data[:, 2]])

# # **********   3D plot
# fig3D_1 = plt.figure()
# xy_by = fig3D_1.add_subplot(projection='3d')
# xy_by.scatter(x_axis, y_axis, by, s=4)

# fig3D_2 = plt.figure()
# xy_by = fig3D_2.add_subplot(projection='3d')
# xy_by.scatter(x_axis, z_axis, by, s=4)

# fig3D_3 = plt.figure()
# xy_by = fig3D_3.add_subplot(projection='3d')
# xy_by.scatter(y_axis, z_axis, by, s=4)

# # **********   2D plot
# fig2D_1 = plt.figure()
# x_by = fig2D_1.add_subplot()
# x_by.scatter(x_axis, by, s=2)

# fig2D_2 = plt.figure()
# y_by = fig2D_2.add_subplot()
# y_by.scatter(y_axis, by, s=2)

# fig2D_3 = plt.figure()
# z_by = fig2D_3.add_subplot()
# z_by.scatter(z_axis, by, s=2)

# plt.show()

# **********   3 variables fitting
# 1次関数の定義
def linear(x, a, b):
    return a*x + b

# 3次関数の定義
def third_order(x, a, b, c, d):
    return a*np.power(x,3) + b*np.square(x) + c*x + d

# ログノーマル関数の定義
def lognormal(x, A, mu, sigma):
    return A / (x * sigma * np.sqrt(2 * np.pi)) * np.exp(- (np.log(x) - mu)**2 / (2 * sigma**2))

# ガウス関数の定義
def gaussian(x, A, mu, sigma):
    return A * np.exp(- (x - mu)**2 / (2 * sigma**2))


# フィッティング関数の定義
# xに依存する項はログノーマル関数＋4つのガウシアン
# yに依存する項は3次関数
# zに依存する項は1次関数
def combined_model(plots,
                   A_log, mu_log, sigma_log,
                   A1, mu1, sigma1,
                   A2, mu2, sigma2,
                   A3, mu3, sigma3,
                   A4, mu4, sigma4,
                   a_y, b_y, c_y, d_y,
                   a_z, b_z):
    x, y, z = plots
    lognorm = lognormal(x, A_log, mu_log, sigma_log)  # ログノーマル関数
    gauss1 = gaussian(x, A1, mu1, sigma1)  # ガウシアン
    gauss2 = gaussian(x, A2, mu2, sigma2)  # ガウシアン
    gauss3 = gaussian(x, A3, mu3, sigma3)  # ガウシアン
    gauss4 = gaussian(x, A4, mu4, sigma4)  # ガウシアン
    third_y = third_order(y, a_y, b_y, c_y, d_y)  # 3次関数
    linear_z = linear(z, a_z, b_z)  # 1次関数

    return lognorm + gauss1 + gauss2 + gauss3 + gauss4 + third_y + linear_z

# def combined_model(plots,
#                    A_log, mu_log, sigma_log,
#                    A1, mu1, sigma1,
#                    A2, mu2, sigma2,
#                    A3, mu3, sigma3,
#                    A4, mu4, sigma4,
#                    a_y, b_y, c_y, d_y):
#     x, y, z = plots
#     lognorm = lognormal(x, A_log, mu_log, sigma_log)  # ログノーマル関数
#     gauss1 = gaussian(x, A1, mu1, sigma1)  # ガウシアン
#     gauss2 = gaussian(x, A2, mu2, sigma2)  # ガウシアン
#     gauss3 = gaussian(x, A3, mu3, sigma3)  # ガウシアン
#     gauss4 = gaussian(x, A4, mu4, sigma4)  # ガウシアン
#     third_y = third_order(y, a_y, b_y, c_y, d_y)  # 3次関数

#     return lognorm + gauss1 + gauss2 + gauss3 + gauss4 + third_y

# 初期値の設定
p0 = [0.0006958040766977115, 6.292789011773518, 0.09067786199274608,  # ログノーマルの初期値
      -8.169333763413994e-07, 560.6826636407911, 27.58503262197233,  # 1つ目のガウス
      2.54948402661953e-07, 656.3964871401387, 32.1459215691162,  # 2つ目のガウス
      1.9893328877436285e-07, 703.2723248019659, 59.25730588891589,  # 3つ目のガウス
      1.2306205278883204e-07, 1061.4207387188872, 405.3453929882518,  # 4つ目のガウス
      1, 10, 10, 10,  # yに依存する項
      0.002, 0  # zに依存する項
      ]

# # 初期値の設定
# p0 = [0.0006958040766977115, 6.292789011773518, 0.09067786199274608,  # ログノーマルの初期値
#       -8.169333763413994e-07, 560.6826636407911, 27.58503262197233,  # 1つ目のガウス
#       2.54948402661953e-07, 656.3964871401387, 32.1459215691162,  # 2つ目のガウス
#       1.9893328877436285e-07, 703.2723248019659, 59.25730588891589,  # 3つ目のガウス
#       1.2306205278883204e-07, 1061.4207387188872, 405.3453929882518,  # 4つ目のガウス
#       1, 10, 10, 10  # yに依存する項
#       ]

# popt, pcov = curve_fit(combined_model, plots, by, p0=p0, maxfev=50000, xtol=1e-4, ftol=1e-4)
try:
    popt, pcov = curve_fit(combined_model, plots, by, p0=p0, maxfev=100000, xtol=10, ftol=10)
except RuntimeError:
    print("Optimal parameters not found. Returning the best found parameters.")

perr = np.sqrt(np.diag(pcov))

A_log, mu_log, sigma_log, \
A1, mu1, sigma1, \
A2, mu2, sigma2, \
A3, mu3, sigma3, \
A4, mu4, sigma4, \
a_y, b_y, c_y, d_y, \
a_z, b_z = popt

# A_log, mu_log, sigma_log, \
# A1, mu1, sigma1, \
# A2, mu2, sigma2, \
# A3, mu3, sigma3, \
# A4, mu4, sigma4, \
# a_y, b_y, c_y, d_y= popt

# **********   2D plot
# fig2D_1 = plt.figure()
# x_by = fig2D_1.add_subplot()
# x_by.scatter(x_axis, by, s=2, label='data')
# x_by.scatter(x_axis, , s=2, label='data')

fig2D_2 = plt.figure()
y_by = fig2D_2.add_subplot()
y_by.scatter(y_axis, by, s=2, label='data')
y_by.plot(y_axis, third_order(y_axis,a_y,b_y,c_y,d_y), label='fit')

fig2D_3 = plt.figure()
z_by = fig2D_3.add_subplot()
z_by.scatter(z_axis, by, s=2, label='data')
z_by.plot(z_axis, linear(z_axis,a_z,b_z), label='fit')