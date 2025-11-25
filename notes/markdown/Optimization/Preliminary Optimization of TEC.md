Try to get approximations for following parameters in the TEC Array.

1) No of Stages
2) Basic Dimensions
	1) **Wedge angle**: Angle of a single TEC (angle * number = 360)
	2) Radial Distance or a single unit (variable distances if possible)
	3) Length of TEC legs
3) Material 
4) Configuration
	 There 2 main configurations that can be followed. Optimization should be done to these 2 differently, and their top candidates can be used to compare the approaches.
	1) Same number of peltier units in all the stages
		this approach needs the area of each stage to increase going from the center of the array to the edge. Therefore, a single TEC unit would a cone with the given Wedge angle, and this approach can be used to model the entire TEC array using one such TEC cone.
	2) different number of peltier units in each and every stage
		in this method, the area of each stage peltier unit can be adjusted as needed, thus it will become a optimization parameter. However, it will remove the ability to model a cone like structure due to non-asymmetric nature.
		`You can choose a ratio of TE number between consecitive stages instead of optimizing for every parameter.`

[[Radial TEC Parameter List]]

Note that you have to consider the fact that the bottom surface will have active heat flux throughout, compared to classic models that considers the heat flux available only at the cold side of the TEC.

### Objective Function
Current candidates for the objective function of the preliminary optimization process.
- [ ] Maximize COP
- [ ] Minimize Power
- [ ] Maximize $\Delta T$

### Parameters

#### Geometric
 1) **No of stages:** usually based on intuition, and will not be calculated mathematically. We can run optimizations for several number of stages and the compare results.
 2) Wedge Angle
 3) Leg Length
 4) Thickness of the Array
 5) Connecting Plate Thickness
#### Operating
 1) Current
 2) Voltage
 3) Setback Coefficient

#### Material
1) Setback coefficient
2) thermal conductivity
3) electrical resistance

# Optimization Problem Formulation

1) [[TEC element level modelling]]
2) [[Resistor Network of TEC modules]]
3) [[System of Equations (TEC Preliminary Optimization)]]
4) [[Objective Function (TEC Project)]]
5) [[Constraints (TEC Preliminary Optimization)]]

# Papers

1) [[A general approach in evaluating and optimizing thermoelectric coolers]]
2) [[Micro thermoelectric cooler - Planar multistage]]


