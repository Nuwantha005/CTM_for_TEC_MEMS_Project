# 1 Modified Geometry with Wider Interconnects

# 2 Initial Results

### 2.1.1 Contour
Current = 0.1 A and Heat Flux = 10

![[image-2.png]]

### 2.1.2 Parametric Sweep for Flux

We can see a clear linear trend with increased flux that increases the Node 0 temperature

![[image-3.png]]

At this configuration, we can clearly sustaain 
## 2.2 Extended Flux

No we apply flux to the entire bottom surface, similar to how it could've been.

![[image-4.png]]


>[!Heat Flux Problem]
>We have been using values such as $50\;W/m^2$, but in reality the surface heat fluxes of CPU's are on the level of $200\;W/cm^2$, which when turns into ISO units, gets to around $2,000,000\;W/m^2$.
>Now I understand why they used 1e6 heal flux in mirco channel simulations in [[CFD Project]]

The heat flux on the top surface of a modern commercial desktop CPU chip, specifically the

### Google Search AI:
**Integrated Heat Spreader (IHS) surface**, is typically in the range of **30 to 70 W/cm²** under load conditions with standard air cooling. However, the heat flux at the actual silicon **die surface** is much higher and can exceed **100 W/cm²**, with modern high-performance chips designed to handle over 112 W/cm². 

Heat Flux Details

- **Integrated Heat Spreader (IHS) Surface:** This is the metal surface you see on top of the CPU when you install it. The heat flux here is generally lower than at the die because the IHS is designed to spread the heat across a larger area for the cooler to dissipate. Typical values are around 30-70 W/cm² during intensive use.

- **CPU Die Surface:** The actual silicon chip underneath the IHS has a much smaller area, leading to significantly higher heat flux (power density). Modern CPUs like the Intel i9-14900K or AMD Ryzen 9 7950X can have heat fluxes well over 100 W/cm² at the die level during heavy multi-threaded loads or stress tests, and some designs have shown capabilities of handling over 200 W/cm² with advanced liquid cooling solutions.

- **Cooling System Limits:** Standard air cooling systems often approach their practical limits around 40-50 W/cm² (at the IHS level), which is why high-end CPUs often require robust liquid cooling solutions to prevent thermal throttling.

- **Hotspots:** The heat generation on a CPU die is not uniform, with specific "hotspots" experiencing even higher localized heat fluxes than the average. 

The high heat flux presents a significant thermal management challenge, and modern CPUs are designed to run at high temperatures (up to 95-100°C) and will automatically throttle performance to prevent damage if these limits are consistently exceeded

Parametric study for current gives around 0.15 as the optimal current.

![[image-5.png]]

After increasing the interconnects radial dimensions, we can decrease the central temperature to around ~140 Celcius.

![[image-6.png]]

LL_w_ic (m)	LL_w_oc (m)	LL_w_is (m)	Objective
5.0000E-5	5.0000E-5	4.0000E-5	461.10
3.5000E-4	5.0000E-5	4.0000E-5	455.04
5.0000E-5	3.5000E-4	4.0000E-5	432.62
5.0000E-5	5.0000E-5	3.5000E-4	429.83
2.5000E-4	2.5000E-4	2.4667E-4	424.01
3.5000E-4	3.5000E-4	3.5000E-4	411.36
2.5000E-4	1.5000E-4	1.4333E-4	433.15
3.5000E-4	3.5000E-4	3.5000E-4	411.36

It  can be noted that all 3 parameters reached their upper bound. Which indicates that the system kind of neglects the Peltier effect and solely rely on increase of highly conductive solids to decrease the temperature of the central node.

> This is concerning because it says that we are better off adding a conductor in that region instead of that TEC
> 
> We have to increase the no of stages and see whether there is a effect



