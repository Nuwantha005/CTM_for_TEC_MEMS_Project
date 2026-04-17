---
Title: Methodology- Geometry - Radial TEC CTM for MERCon
date: 2026-04-12 18:24
type: publication
tags:
  - publication
  - MEMS
status: complete
---
# 1. Geometry
## 1.1. overall
The TEC is configured such that a single stage occupies a concentric, with multiple stages stacked radially outward. The TEC is located on top of the  surface being cooled, means the heat sources has to be modeled underneath the TEC layer. This bottom layer is made out of high thermal conductive material such as Single Crystal Silicon, which is common in MEMS. Between the heat generation layer and the TEC layer, an insulation layer is required to prevent electricity from conducting between 2 layers but it should allow heat to pass through. Part of the bottom layer passes through the insulation layer and emerges in the middle of the radial TEC, allowing a low resistance path for heat to travel from bottom of the layer to TEC absorption side. This part is modeled as a cylinder and made out of the same material as the bottom layer.
## 1.2. TEC Element
TEC elements are sandwiched between insulators in both radial and azimuthal directions to keep the charged conductors isolated. Inside an element, there is another azimuthal insulator between p and n TE legs because the junction has to be maintained at one radial end. Radial insulators are concentric rings of constant width while azimuthal insulators expand radially.

Inside a TEC element, there is a copper interconnect that connects p and n type TE legs together, which is the high temperature side that heat is being pumped from. On the other end, there is the copper connector called "outerconnect" that connects the consecutive p and n legs between wedges. Using this repeating arrangement, for a given stage, only one wedge has to be connected to external source, significantly simplifying the wiring process.
# 2. Parameters
## 2.1. List
The full set of geometric parameters that can completely and uniquely define the geometry of the TEC is listed in the [table 01].

| Symbol                | Name                       | Physical Meaning                                                              | Figure | Units         |
| --------------------- | -------------------------- | ----------------------------------------------------------------------------- | ------ | ------------- |
| $r_{\text{base}}$     | Base Radius                | Radius of the entire TEC array                                                |        | $\mu m$       |
| $\theta$              | Wedge angle                | Azimuthal span of a single wedge                                              |        | rad           |
| $t_{\text{TEC}}$      | TEC layer Thicknes         | Thickness of the layer containing TEC                                         |        | $\mu m$       |
| $t_{\text{ins}}$      | Insulator layer thickness  | Thickness of the electrical insulation layer between TEC and generation layer |        | $\mu m$       |
| $t_{\text{gen}}$      | Generation layer thickness | Thickness of the heat generation layer.                                       |        | $\mu m$       |
| $r_{\text{cyl}}$      | Central cylinder radius    | Radius of the central                                                         |        | $\mu m$       |
| $r_{\text{in},i}$     | Inner radius               | Inner radius of a single TEC element excluding radial insulator               |        | $\mu m$       |
| $r_{\text{out},i}$    | Outer radius               | Inner radius of a single TEC element excluding radial insulator               |        | $\mu m$       |
| $L_i$                 | Length                     | Radial span of a single TEC element                                           |        | $\mu m$       |
| $W_{\text{is},i}$     | Radial insulator width     | Radial span of the radial insulator                                           |        | $\mu m$       |
| $W_{\text{az},i}$     | Azimuthal insulator width  | Azimuthal span of the azimuthal insulator in linear dimensions                |        | $\mu m$       |
| $\beta_{\text{ic},i}$ | Interconnect angle         | Azimuthal angle of the copper interconnect                                    |        | rad           |
| $t_{\text{ic},i}$     | Interconnect thickness     | Thickness of the interconnect                                                 |        | $\mu m$       |
| $W_{\text{ic},i}$     | Interconnect thickness     | Radial span of the copper interconnect                                        |        | $\mu m$       |
| $\beta_{\text{oc},i}$ | Outerconnect angle         | Azimuthal angle of the copper outerconnect                                    |        | rad           |
| $t_{\text{oc},i}$     | Outerconnect thickness     | Thickness of the outerconnect                                                 |        | $\mu m$       |
| $W_{\text{oc},i}$     | Outerconnect thickness     | Radial span of the copper outerconnect                                        |        | $\mu m$       |
| $M$                   | No. of Stages              | Total number of stages in the mutistage TEC                                   |        | $\mathbb N^+$ |

## 2.2. Constraints
To maintain the physical structure of the TEC, these geometric parameters should be constrained by rules shown in [table 02]

| Reason                                             | Formula                                                                                                 |
| -------------------------------------------------- | ------------------------------------------------------------------------------------------------------- |
| Avoiding interconnect - outerconnect contact       | $W_{\text{ic},i}+W_{\text{oc},i}<L$                                                                     |
| Radial span of all stages adding up to base radius | $r_{\text{base}} = r_{\text{cyl}} + \sum_{i=1}^{M} \left(L_i +W_{\text{is},i}\right)+W_{\text{is},M+1}$ |
| Azimuthal insulator must fit inside the element    | $2W_{\text{az},i} < r_{\text{in},i}\theta$                                                              |
| Central cylinder must be contained in the TEC      | $0<r_\text{cyl} < r_\text{base}$                                                                        |

# 3. Dimensional Reduction
Observing the [table 1], its apparent that there are 9 independent parameters for a single stage, which means the complexity of the problem is $O(M)$. Performing any kind of optimization on this many parameters will lead to performance constraints, but due to the linear nature of the system that gets formed, its not computationally intractable in modern computers. However the issue is that we will run into the "Curse of Dimensionality", which turn solution space sparse, and finding optimum solutions increasingly difficult.

Therefore, we need to perform some sort of dimensional reduction. Here we turn several parameters into ratios which depends on other parameters, and then fix these parameters across the stages. The ratio parameters, which would be the independent variables of the optimization problem are listed on [table 03].
## 3.1. Reduced Parameter Definitions

| Symbol                | Name                            | Definition                                       | Replaces          | Bounds    |
| --------------------- | ------------------------------- | ------------------------------------------------ | ----------------- | --------- |
| $N$                   | No. of wedges                   | $N=2\pi/\theta$                                  | $\theta$          |           |
| $f_L$                 | Length ratio                    | $f_L = L_{i+1}/L_i$                              | $L$               |           |
| $f_\text{ic,w}$       | Interconnect width fraction     | $f_\text{ic,w}=W_\text{ic}/L_i$                  | $W_\text{ic}$     | $(0,0.5)$ |
| $f_\text{ic,t}$       | Interconnect thickness fraction | $f_\text{ic,t}=t_\text{ic}/t_\text{TEC}$         | $t_\text{ic}$     | $(0,1)$   |
| $f_{\text{ic},\beta}$ | Interconnect angle fraction     | $f_{\text{ic},\beta}=\beta_\text{ic}/\theta$     | $\beta_\text{ic}$ | $(0,1)$   |
| $f_\text{oc,w}$       | Outerconnect width fraction     | $f_\text{oc,w}=W_\text{oc}/L_i$                  | $W_\text{oc}$     | $(0,0.5)$ |
| $f_\text{oc,t}$       | Outerconnect thickness fraction | $f_\text{oc,t}=t_\text{oc}/t_\text{TEC}$         | $t_\text{oc}$     | $(0,1)$   |
| $f_{\text{oc},\beta}$ | Outerconnect angle fraction     | $f_{\text{oc},\beta}=\beta_\text{oc}/\theta$     | $\beta_\text{oc}$ | $(0,1)$   |
| $f_\text{f}$          | Fill Factor                     | $f_\text{f}=1-W_{\text{az}}/(\text{Arc Length})$ | $W_{\text{az}}$   | $(0,1)$   |

## 3.2. Recovery Map
[Table 04] shows the relationship between original geometric parameters and optimization variables - specifically how to recover geometric parameters from optimization variables, alongside other details relevant to their relationship.

| Original Parameter    | Recovery Formula                                           | Notes                                                                            |
| --------------------- | ---------------------------------------------------------- | -------------------------------------------------------------------------------- |
| $\theta$              | $2\pi / N$                                                 | —                                                                                |
| $r_{\text{in},i}$     | $r_\text{cyl} + i\cdot W_\text{is} + \sum_{j=1}^{i-1} L_j$ | $r_{\text{in},1} = r_\text{cyl} + W_\text{is}$                                   |
| $r_{\text{out},i}$    | $r_{\text{in},i} + L_i$                                    | —                                                                                |
| $L_i$                 | $L_1 \cdot f_L^{(i-1)}$                                    | $L_1 = \frac{(r_\text{base} - r_\text{cyl} - (M+1)W_\text{is})(1-f_L)}{1-f_L^M}$ |
| $W_{\text{az}}(r)$    | $(1-f_f)\cdot r\theta$                                     | Continuous radial function, not stage-indexed                                    |
| $W_{\text{ic},i}$     | $f_{\text{ic,w}}\cdot L_i$                                 | —                                                                                |
| $t_{\text{ic},i}$     | $f_{\text{ic,t}}\cdot t_\text{TEC}$                        | —                                                                                |
| $\beta_{\text{ic},i}$ | $f_{\text{ic},\beta}\cdot \theta$                          | Stage-independent since $\theta$ is uniform                                      |
| $W_{\text{oc},i}$     | $f_{\text{oc,w}}\cdot L_i$                                 | —                                                                                |
| $t_{\text{oc},i}$     | $f_{\text{oc,t}}\cdot t_\text{TEC}$                        | —                                                                                |
| $\beta_{\text{oc},i}$ | $f_{\text{oc},\beta}\cdot \theta$                          | Stage-independent since $\theta$ is uniform                                      |

Notes:
- $W_{\text{is}}$ is fixed across all stages and its independent, so it hasn't been reduced to a fractional variable.
- $W_{\text{az},i}$ is radially expanding, and not constant across the stages. Therefore $f_\text{f}$ acts more as an angular quantity, but just expressed in arch length terms because the numerical value conveys what fraction of total radial area is occupied by TE materials.
- Generation layer thickness is only relevant for calculating heat flow resistance on that layer only and it heat generation is calculated using surface area, not volume.e


---
# 4. Exploration


# 5. References
