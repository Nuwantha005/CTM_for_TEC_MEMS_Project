---
Title: Boilerplate Sections - Radial TEC CTM for MERCon
date: 2026-04-11 20:58
type: publication
tags:
  - publication
  - MEMS
status: complete
---
Title: A Parametric Compact Thermal Modelling Framework for Rapid Design Space Exploration of Radial Multistage Thermoelectric Coolers
# 1. Abstract
This paper builds a Compact Thermal Model (CTM) for a multi stage thermoelectirc cooler (TEC) under the assumption of 1D heat transfer in thermally isolated wedges. The geometry of the TEC array would be parameterized using the dimensions of different components, which gets used to calculate the thermal or electrical resistances of each components. Adding all the thermal resistances together, a resistor network will be formed, and electrical resistances will be used to calculate the heat generation rate. All the calculated values will be used to form a linear system of equations with node temperatures as the unknowns. This system can be solved with easy using standard matrix inversion, and be used in the design exploration stage to identify geometric parameters, as well as to optimize for specific objectives. The CTM was verified against an high fidelity simulation in COMSOL Multiphysics.

---
# 2. Introduction
Microelectromechanical systems (MEMS) with requires ways to dissipate the heat they generated during their operations, and at low power levels, this done using passive cooling. When the power consumption of the MEMS device increases, some form of active cooling is required to maintain the healthy operating temperatures. At the scales these devices are constructed, traditional active cooling methods such as forced convection or refrigerant loop based systems becomes difficult to implement. Coolers that utilizes thermoelectric effect provides ease of miniaturization due to their solid-state nature, and be fabricated using traditional electronics fabrication techniques.
## 2.1. Thermoelectirc Effect
When 2 conductors with different material properties are connected to each other, and forced a current to flow through them, a phenomenon known as Peltier effect (Thermoelectirc Effect) occurs. The charge carriers carry heat, and arrange them in such a way so that both conductors carry heat in the same direction can lead to a solid state device that "pumps" heat from one end to another. Even though this phenomenon occurs in all type of conductors, its too weak to make any observable difference in common materials.
[image]
The governing equation for a thermoelectric (TE) junction is given by$$Q_c = (S_p - S_n) I T_c - K (T_h - T_c) - \tfrac{1}{2} I^2 R $$where $Q_c$ denotes the heat absorbed at the cold side of the junction. The three terms represent the Peltier effect, thermal conduction, and Joule heating, respectively. $S_p$ and $S_n$ are the Seebeck coefficients of the p-type and n-type materials.

In this work, the terms cold side and hot side are defined based on device operation, not absolute temperature. The cold side refers to the location where heat is extracted ($Q_c$), while the hot side refers to where heat is rejected.

It is important to note that, depending on operating conditions, the cold side may not necessarily be at a lower temperature than the hot side. In such cases, the direction of passive heat conduction is from the cold side to the hot side, and the conduction term $-K(T_h - T_c)$ becomes positive, indicating that conduction assists heat removal.

Semiconductor materials such as bismuth Telluride ($Bi_2Te_3$) has far larger setback coefficients that let their charged particles carry large amount of heat compared to regular materials. In practice, two p and n type blocks are used to make the TE junction.
## 2.2. Related Work
Mico scale Thermoelectric devices are used for many purposes such as heating, cooling, and power generation. Of these, cooling is heavily researched area due to its perceived advantages cooling electronics components, as well as MEMS devices. At this scale, Thermoelectric coolers (TEC) are mostly designed as layers, and can be designed in different arrangements based on space availability and heat transfer direction.
![[image.png|346]]
Of these, the most notable 2 design architectures are ones that transfers heat in the direction of normal to the plane layer ones that transfers heat along the layer. [1](cite: O. Owoyele, S. Ferguson, and B. T. O’Connor, “Performance analysis of a thermoelectric cooler with a corrugated architecture,” Applied Energy, vol. 147, pp. 184–191, June 2015, doi: 10.1016/j.apenergy.2015.01.132 ) Radial TEC belongs to latter and incorporates a modular and asymmetric design, which allows mathematical modelling and simulations to be run on thermal "islands" with symmetric boundary conditions rather than modelling the entire structure.

Apart from the arrangement, TEC devices can have multiple stages, that increases the temperature difference obtainable between the heat source and the sink. However, the heat pumping capacity of final stages should be larger than that of previous stages, because latter stages needs to pump added heat that was generated in previous stages. This usually requires to increase the surface area of the latter stages, which results in a pyramid like structure. (cite [[Micro thermoelectric cooler - Planar multistage]]) The radial design allows this cross section area addition naturally when stages are placed in increasing radial order, which makes it a really good candidate for multi stage design.

This paper intends to mathematically model such multi stage radial TEC and build an framework to calculate the junction temperatures, as was done in [[Micro thermoelectric cooler - Planar multistage]].

---
# 3. Exploration


# 4. References
