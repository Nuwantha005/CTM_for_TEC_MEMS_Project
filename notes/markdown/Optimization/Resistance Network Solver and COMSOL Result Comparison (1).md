This was done to ensure the reliability resistance network based solver before proceeding into the preliminary optimization.

# 1 Setup

Data of the model are shown in [[Radial TEC Parameter List]]'s parameter table's default column, and the material table' default column. A 3 stage conic section was modeled using SolidWorks.

![[image.png|Annular TEC Segment]]

TSV were added only to the first stage and only 3 were possible given the dimensions of the copper interconnect.

![[image-1.png|TSV]]

The chip surface and the central cylinder was modeled as a one single element and SOI layer was added on top of that.
![[image-2.png|Chip Layer+Cylinder and SOI layer (transparent)]]

# 2 Observations

> No matter the value of the heat flux (i.e. for 0.05 $W/m^2$ and 100 $100 W/m^2$) the maximum temperature of the chip layer remain around ~400K. 

This is concerning. It was also noted that the highest temperature is present at the copper outerconnects rather than the chip surface. Temperature at the copper interconnects were rather low.


# 3 Results

For all the current values, this shows the same graph. Its like there is no effect on changing the current. Like the Thermoelectric effect is not working.


![[image-4.png]]

![[image-5.png]]

Its the same even at 25A. Idk what to do. when i change the heat flux from 0.1 to 50 it stayed the same. When i removed the heat flux boundary and added constant temperature boundary, it still stayed the same. Highest temperature around 380K.

# 4 Parametric Sweep
Thermoelectric effect works fine. when i reversed the thermoelectric direction, the maximum temperature increased.

I changed the conductivity of the insulator to check whether its working as a dam. But nothing changed and that removes the possibility as well.

The problem isn't the large temperatures themselves. The problem is that nothing changes when the current or heat flux changed.

I reduced the dimensions of the copper interconnects and this increased the maximum temperature. So it suggests that joule heating works to increase the temperature, because reducing the interconnects increases their resistance.

![[image-6.png|Node 0 temperature against outerconnect width]]

That suspicion is conformed by this parametric sweep that was done between the outerconnect width and the central node temperature, which shows clear downward trend in the temperature as the width increases.

# 5 Optimization 

An optimization run was carried out to determine the dimensions of the wires that gives the minimum central node temperature.

## Run 1

![[image-7.png]]

![[image-9.png]]


![[image-8.png]]

LL_w_oc (m)	LL_w_ic (m)	LL_beat_oc (rad)	LL_beta_ic (rad)	Objective
2.0000E-4	2.0000E-4	0.43633	0.42883	314.96

>**Verdict**: We need wider copper wires. And to lower the temperature more than this, we need to change no. of stages and other parameters, thus needs to use preliminary optimization.

## Run 2

Parameter List with bounds
![[image-10.png]]


![[image-11.png]]
