## RK4y (under development)
<br>
===== Structure ===== <br>

&bull;FourGaussian.py;<br>
For preparing the magnetic field of the trim coils, Btrim_sum. We have 10 trim coils, and each magnetic field, Btrim, is represented by summing four gaussians. Finally, we will get the field as the sum of 10 trim coils.<br>
Variables are in variables_FourGaussians.py. There are 4 list, amp, mean, sigma, and offset. Note that actual mean values of gaussians are mean-offset.<br>
<br>

&bull;Enge.py;<br>
For preparing the fringe magnetic field, enge(), of the main coils. The function is ege function.<br>
Variables are in variables_Enge.py. This is 6 order enge function.<br>
<br>

&bull;Bz.py;<br>
For preparing the magnetic field of the main coils, Bz. <br>

<br>
&bull;position.py;<br>
For calculating a proper coordinate for each dipole magnet. Note that the coordinates of the dipoles are rotated　Euclidean space.<br>
The angle of rotation is based on each geometry of the dipole.<br>
This function is used in Diff.py.<br>

<br>
&bull;rk45.py;<br>
Run this code. We will use scipy.integrate.solve_ivt function in RK45 mode. Please note that we need to convert unit many times in this code.<br>
If you want, you can confirm the results of this code and compare to ones of an explicit Runge-Kutta code, which is in dump/main.py.<br>
