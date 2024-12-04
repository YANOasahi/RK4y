## RK4y (under development)
<br>
# Structure<br>
&bull;FourGaussian.py;<br>
For preparing the magnetic field of the trim coils, Btrim_sum. We have 10 trim coils, and each magnetic field, Btrim, is represented by summing four gaussians. Finally, we will get the field as the sum of 10 trim coils.<br>
Variables are in variables_FourGaussians.py. There are 4 list, amp, mean, sigma, and offset. Note that actual mean values of gaussians are mean-offset.<br>
<br>
&bull;Enge.py;<br>
For preparing the fringe magnetic field, enge(), of the main coils. The function is ege function.<br>
Variables are in variables_Enge.py. This is 6 order enge function.<br>
<br>
&bull;Bz.py;<br>
For preparing the magnetic field of the main coils, Bz. 