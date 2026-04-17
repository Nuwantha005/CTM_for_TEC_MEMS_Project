---
Title: Methodology- Resistance Calc - Radial TEC CTM for MERCon
date: 2026-04-13 14:26
type: publication
tags:
  - publication
  - MEMS
status: complete
---
To build the thermal network, the thermal resistances of each component must be calculated. Apart from that, the electrical resistances of the TEC elements must also be calculated to account for Joule heating effects.
# 1. Resistor Arrangement
The thermal resistances of components depends on the heat transfer direction which is radial. Its assumed that no azimuthal or vertical heat transfer occurs within layers other than the heat transfer between generation layer and TEC layer, which will be modeled explicitly.

Due to the existence of many components, their individual resistances must be lumped together to be used as the conduction term in governing equation. Due to the multistage nature of the system, we need to pick "nodes" that connects these resistances with each other, and the unknowns of the system would be temperatures at these nodes.

This paper uses the boundary between azimuthal insulator and radial insulator as the node in question. Therefore all the thermal resistances between these repeating interfaces would be lumped into a single resistor using their parallel and series connections as shown in figure []

Bottom generation layer also need explicit resistance modelling, and nodes in that layer are located right under the nodes in TEC layer as shown in figure [] The heat conduction from top generation layer to TEC layer will be modeled as resistances located between top and bottom nodes as shown in the same figure.
# 2. Thermal Resistances
The general formula for thermal resistance with variable cross section area is,$$ R_{\text{t}} = \int_{x=0}^{L} \frac{dx}{\kappa \cdot A(x)}$$Where $\kappa$ is the thermal conductivity. Almost all the expressions obtained in the following calculations can be reduced to following form.$$\int_{x_1}^{x_2} \frac{1}{ax}dx 
= \frac{1}{a}  \ln\left| \frac{x_2}{x_1} \right|$$
## 2.1. Special Cases
There are several calculations that needs extra attention due to their deviations from the formula mentioned above.
### 2.1.1. TE Legs
When calculating the thermal conductivity of the TE legs, we have to divide the region radially into 3 separate parts.
- Region containing interconnect: $r\in (r_{\text{in}},r_{\text{in}}+W_\text{ic})$
- Region without any connections: $(r_{\text{in}}+W_\text{ic},r_{\text{out}}-W_\text{oc})$
- Region containing outerconnect: $(r_{\text{out}}-W_\text{oc},r_{\text{out}})$
The cross section area in the first and last regions will be reduced due to the existence of connections. The full calculated table is shown below.
### 2.1.2. Vertical
We assume that the heat travels between the nodes in 2 layers passes through a resistance equivalent of the area occupied by half of each stage connected to the node. So for the $i^\text{th}$ node, the area is between following radial span.$$\left(r_{\text{in},i}-\frac{L_{i-1}}{2},r_{\text{in},i}+\frac{L_{i}}{2} \right)$$Since the area is given by $\frac{\theta}{2}(r_2^2-r_1^2)$, the resistance is$$R_{\text{ve},i} = \frac{2 t_{\text{ins}}}{\kappa_{\text{ins}}\theta  \left[\left(r_{\text{in},i}-\frac{L_{i-1}}{2}\right)^2 - \left( r_{\text{in},i}+\frac{L_{i}}{2}\right)^2\right]}$$
[insert figure if space persists]
### 2.1.3. Cylinder
In the resistor network, the central cylinder acts as a common node for both layers, denoted by \textbf{Node 0}. This is because it's a single block spanning all 3 layers as shown in \cref{fig:cylinder_iso}. Therefore, it is required to calculated the resistance between node 0 and node 1 of both TEC layer and the chip layer.

This involves calculating resistance from center $(r=0)$ to the outer surface of the cylinder $(r=r_\text{cyl})$. The standard formula [] breaks down here because there is a singularity at $r=0$. We use the standard solution for a cylinder with volumetric generation is $T_{\text{center}} - T_{\text{surf}} = \frac{Q}{4\pi k t}$ as an approximation. Adjusting for the wedge angle $\theta$ (where the full circle is $2\pi$), the factor becomes $2\theta$. From this, the thermal resistance between node 0 and node 1 in the TEC layer can be derived as$$R_{\text{cyl}} = \frac{1}{2\theta \kappa t}$$For the generation layer, this remains the same.$$ R_{\text{gen},0 \to 1} =  \frac{1}{2\theta \kappa_{\text{gen}} t_{\text{gen}}}$$But for the TEC layer, the resistance form first radial insulator should be added in series because it won't get absorbed into the TEC element. (it absorbs the insulator after the element) Hence,$$R_{\text{TEC},0 \to 1} = \frac{1}{2\theta \kappa_{\text{gen}} t_{\text{TEC}}} + R_{\text{is},0}$$
## 2.2. Summery

| **Component**                                              | **Area**                                                                                                    | **Limits**                                                                        | **Final Expression (Original)**                                                                                                                                                | **Reduced Expression (Optimization Variables)**                                                                                                                                                  |
| ---------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| **TE Leg I** ($R_{\text{t,TE,I}}$)                         | $\frac{r}{2} \left[ \theta \cdot t_{\text{TEC}} (2f_f - 1) - \beta_{\text{ic}} \cdot t_{\text{ic}} \right]$ | $(r_{\text{in}}, r_{\text{in}}+W_{\text{ic}})$                                    | $$\frac{2}{\kappa_{\text{TE}} [\theta t_{\text{TEC}} (2f_f - 1) - \beta_{\text{ic}} t_{\text{ic}}]} \ln \left( \frac{r_{\text{in}} + W_{\text{ic}}}{r_{\text{in}}} \right)$$   | $$\frac{N}{\pi \kappa_{\text{TE}} t_{\text{TEC}} [2f_f - 1 - f_{\text{ic},\beta} f_{\text{ic,t}}]} \ln \left( 1 + f_{\text{ic,w}} \frac{L_i}{r_{\text{in},i}} \right)$$                          |
| **TE Leg II** ($R_{\text{t,TE,II}}$)                       | $\frac{1}{2}r\theta \cdot t_{\text{TEC}} (2f_f - 1)$                                                        | $(r_{\text{in}}+W_{\text{ic}} , r_{\text{out}}-W_{\text{oc}})$                    | $$\frac{2}{\kappa_{\text{TE}} \theta \cdot t_{\text{TEC}} (2f_f - 1)} \ln \left( \frac{r_{\text{out}} - W_{\text{oc}}}{r_{\text{in}} + W_{\text{ic}}} \right)$$                | $$\frac{N}{\pi \kappa_{\text{TE}} t_{\text{TEC}} (2f_f - 1)} \ln \left( \frac{r_{\text{in},i} + L_i(1 - f_{\text{oc,w}})}{r_{\text{in},i} + f_{\text{ic,w}}L_i} \right)$$                        |
| **TE Leg III** ($R_{\text{t,TE,III}}$)                     | $\frac{r}{2} \left[ \theta \cdot t_{\text{TEC}} (2f_f - 1) - \beta_{\text{oc}} \cdot t_{\text{oc}} \right]$ | $(r_{\text{out}}-W_{\text{oc}}, r_{\text{out}})$                                  | $$\frac{2}{\kappa_{\text{TE}} [\theta t_{\text{TEC}} (2f_f - 1) - \beta_{\text{oc}} t_{\text{oc}}]} \ln \left( \frac{r_{\text{out}}}{r_{\text{out}} - W_{\text{oc}}} \right)$$ | $$\frac{N}{\pi \kappa_{\text{TE}} t_{\text{TEC}} [2f_f - 1 - f_{\text{oc},\beta} f_{\text{oc,t}}]} \ln \left( \frac{r_{\text{in},i} + L_i}{r_{\text{in},i} + L_i(1 - f_{\text{oc,w}})} \right)$$ |
| **Interconnect** ($R_{\text{t,ic}}$)                       | $r \cdot \beta_{\text{ic}} \cdot t_{\text{ic}}$                                                             | $(r_{\text{in}} ,r_{\text{in}}+W_{\text{ic}})$                                    | $$\frac{1}{\kappa_{\text{TE}} \cdot \beta_{\text{ic}} \cdot t_{\text{ic}}} \ln\left( \frac{r_{\text{in}} + W_{\text{ic}}}{r_{\text{in}}} \right)$$                             | $$\frac{N}{2\pi \kappa_{\text{TE}} f_{\text{ic},\beta} f_{\text{ic,t}} t_{\text{TEC}}} \ln\left( 1 + f_{\text{ic,w}}\frac{L_i}{r_{\text{in},i}} \right)$$                                        |
| **Outerconnect** ($R_{\text{t,oc}}$)                       | $r \cdot \beta_{\text{oc}} \cdot t_{\text{oc}}$                                                             | $(r_{\text{out}}-W_{\text{oc}} ,r_{\text{out}})$                                  | $$\frac{1}{\kappa_{\text{TE}} \cdot \beta_{\text{oc}} \cdot t_{\text{oc}}} \ln\left( \frac{r_{\text{out}}}{r_{\text{out}} - W_{\text{oc}}} \right)$$                           | $$\frac{N}{2\pi \kappa_{\text{TE}} f_{\text{oc},\beta} f_{\text{oc,t}} t_{\text{TEC}}} \ln\left( \frac{r_{\text{in},i} + L_i}{r_{\text{in},i} + L_i(1 - f_{\text{oc,w}})} \right)$$              |
| **Radial Insulator** ($R_{\text{t,is}}$)                   | $r \cdot \theta \cdot t_{\text{TEC}}$                                                                       | $r_{\text{out}} \to r_{\text{out}}+W_{\text{is}}$                                 | $$\frac{1}{\kappa_{\text{is}} \cdot \theta \cdot t_{\text{TEC}}} \ln\left( \frac{r_{\text{out}} + W_{\text{is}}}{r_{\text{out}}} \right)$$                                     | $$\frac{N}{2\pi \kappa_{\text{is}} t_{\text{TEC}}} \ln\left( 1 + \frac{W_{\text{is}}}{r_{\text{in},i} + L_i} \right)$$                                                                           |
| **Azimuthal Insulator** ($R_{\text{az}}$)                  | $(1-f_f) \theta \cdot t_{\text{TEC}} \cdot r$                                                               | $(r_{\text{in}}, r_{\text{out}})$                                                 | $$\frac{1}{\kappa_{\text{az}} (1-f_f) \theta \cdot t_{\text{TEC}}} \ln \left( \frac{r_{\text{out}}}{r_{\text{in}}} \right)$$                                                   | $$\frac{N}{2\pi \kappa_{\text{az}} (1-f_f) t_{\text{TEC}}} \ln \left( 1 + \frac{L_i}{r_{\text{in},i}} \right)$$                                                                                  |
| **Chip layer** ($R_{\text{gen}}$)                          | $r\cdot\theta\cdot t_{\text{gen}}$                                                                          | $(r_{\text{in}}, r_{\text{out}})$                                                 | $$\frac{1}{\kappa_{\text{gen}} \cdot\theta\cdot t_{\text{gen}}}\ln\left( \frac{r_{\text{out}}}{r_{\text{in}}}\right)$$                                                         | $$\frac{N}{2\pi \kappa_{\text{gen}} t_{\text{gen}}} \ln\left( 1 + \frac{L_i}{r_{\text{in},i}} \right)$$                                                                                          |
| **Vertical** ($R_{\text{ve},i}$) _(Special)_               | N/A                                                                                                         | $\left(r_{\text{in},i}-\frac{L_{i-1}}{2},r_{\text{in},i}+\frac{L_{i}}{2} \right)$ | $$\frac{2 t_{\text{ins}}}{\kappa_{\text{ins}}\theta \left[\left(r_{\text{in},i}-\frac{L_{i-1}}{2}\right)^2 - \left( r_{\text{in},i}+\frac{L_{i}}{2}\right)^2\right]}$$         | $$\frac{N t_{\text{ins}}}{\pi \kappa_{\text{ins}} \left[\left(r_{\text{in},i} - \frac{L_i}{2f_L}\right)^2 - \left(r_{\text{in},i} + \frac{L_i}{2}\right)^2\right]}$$                             |
| **Cyl. Generation** ($R_{\text{gen},0 \to 1}$) _(Special)_ | N/A                                                                                                         | $(0, r_{\text{cyl}})$                                                             | $$\frac{1}{2\theta \kappa_{\text{gen}} t_{\text{gen}}}$$                                                                                                                       | $$\frac{N}{4\pi \kappa_{\text{gen}} t_{\text{gen}}}$$                                                                                                                                            |
| **Cyl. TEC Layer** ($R_{\text{TEC},0 \to 1}$) _(Special)_  | N/A                                                                                                         | $(0, r_{\text{cyl}} + W_{\text{is}})$                                             | $$\frac{1}{2\theta \kappa_{\text{gen}} t_{\text{TEC}}} + R_{\text{is},0}$$                                                                                                     | $$\frac{N}{4\pi \kappa_{\text{gen}} t_{\text{TEC}}} + \frac{N \ln\left( 1 + \frac{W_{\text{is}}}{r_{\text{cyl}}} \right)}{2\pi \kappa_{\text{is}} t_{\text{TEC}}}$$                              |
| **Elec. TE Leg II** ($R_{\text{e,TE,II}}$)                 | N/A                                                                                                         | $(r_{\text{in}}+W_{\text{ic}} , r_{\text{out}}-W_{\text{oc}})$                    | $$\frac{2\rho_{\text{TE}}}{\theta \cdot t_{\text{TEC}} (2f_f - 1)} \ln \left( \frac{r_{\text{out}} - W_{\text{oc}}}{r_{\text{in}} + W_{\text{ic}}} \right)$$                   | $$\frac{N \rho_{\text{TE}}}{\pi t_{\text{TEC}} (2f_f - 1)} \ln \left( \frac{r_{\text{in},i} + L_i(1 - f_{\text{oc,w}})}{r_{\text{in},i} + f_{\text{ic,w}}L_i} \right)$$                          |
| **Elec. Interconnect** ($R_{\text{e,ic}}$)                 | N/A                                                                                                         | $(r_{\text{in}} ,r_{\text{in}}+W_{\text{ic}})$                                    | $$\frac{\rho_{\text{TE}} \beta_{\text{ic}}}{t_{\text{ic}} \ln\left( \frac{r_{\text{in}} + W_{\text{ic}}}{r_{\text{in}}} \right)}$$                                             | $$\frac{2\pi \rho_{\text{TE}} f_{\text{ic},\beta}}{N f_{\text{ic,t}} t_{\text{TEC}} \ln\left( 1 + f_{\text{ic,w}}\frac{L_i}{r_{\text{in},i}} \right)}$$                                          |
| **Elec. Outerconnect** ($R_{\text{e,oc}}$)                 | N/A                                                                                                         | $(r_{\text{out}}-W_{\text{oc}} ,r_{\text{out}})$                                  | $$\frac{\rho_{\text{TE}} \beta_{\text{oc}}}{t_{\text{oc}} \ln\left( \frac{r_{\text{out}}}{r_{\text{out}}-W_\text{oc}} \right)}$$                                               | $$\frac{2\pi \rho_{\text{TE}} f_{\text{oc},\beta}}{N f_{\text{oc,t}} t_{\text{TEC}} \ln\left( \frac{r_{\text{in},i} + L_i}{r_{\text{in},i} + L_i(1 - f_{\text{oc,w}})} \right)}$$                |

# 3. Electrical Resistances
The electrical resistance of each component will be calculated based on geometry and material properties. For a variables cross sectional area, the electrical resistance can be calculated using:$$R_{\text{e}} = \int_{x=0}^{L} \frac{\rho_{\text{TE}} \cdot dx}{A(x)}$$For the TE leg, we ignore the areas in regions I and III for simplicity and consider only area II conducts electricity and generates heat. This assumption is valid because only region II sits inside the direction connection between interconnect and outerconnect - the electricity travel path. After solving the integral, we get$$R_{\text{e,TE,II}} = \frac{2\rho_{\text{TE}}}{\theta \cdot t_{\text{TEC}} (2f_f - 1)} \ln \left( \frac{r_{\text{out}} - W_{\text{oc}}}{r_{\text{in}} + W_{\text{ic}}} \right)$$The interconnects transfer electricity mainly in azimuthal direction. Therefore we need to consider the resistance of a small radial strip with differential resistance $dR = \rho_{\text{TE}} \frac{r \beta_{\text{ic}}}{t_{\text{ic}} dr}$. Since the these resistances are in parallel, the integration should be done over the conductance. $$G_{\text{total}} = \int_{r_{\text{1}}}^{r_{\text{2}}} \frac{t_{\text{ic}}}{\rho_{\text{TE}} \beta_{\text{ic}}} \frac{1}{r} dr = \frac{t_{\text{ic}}}{\rho_{\text{TE}} \beta_{\text{ic}}} \ln\left( \frac{r_{\text{2}}}{r_{\text{1}}} \right)$$Therefore, the resistance for the interconnect is,$$R_{\text{e,ic}} = \frac{\rho_{\text{TE}} \beta_{\text{ic}}}{t_{\text{ic}} \ln\left( \frac{r_{\text{in}} + W_{\text{ic}}}{r_{\text{in}}} \right)}$$Similarly for the outerconnect,$$R_{\text{e,oc}} = \frac{\rho_{\text{TE}} \beta_{\text{oc}}}{t_{\text{oc}} \ln\left( \frac{r_{\text{out}}}{r_{\text{out}}-W_\text{oc}} \right)}$$
# 4. Lumping Resistances
Now we need to lump all the resistors between two adjacent interconnect nodes together.The resistor arrangement is shown in \cref{fig:TEC_resistor_diagram_iso} is shown in a simplified manner in \cref{fig:lumped_resistor_network}.

Now all the resistors can be lumped together to a single resistor. First, let's lump together the parallel resistors in the regions I and III together.

$$
  \frac{1}{R_{\text{t,TE-p,I,combined}}} = \frac{1}{R_{\text{t,TE-p,I}}} + \frac{1}{R_{\text{t,ic}}/2}
$$

$$
  \frac{1}{R_{\text{t,TE-p,III,combined}}} = \frac{1}{R_{\text{t,TE-p,III}}} + \frac{1}{R_{\text{t,oc}}/2}
$$

Similarly for n-TE side,

$$
  \frac{1}{R_{\text{t,TE-n,I,combined}}} = \frac{1}{R_{\text{t,TE-n,I}}} + \frac{1}{R_{\text{t,ic}}/2}
$$

$$
  \frac{1}{R_{\text{t,TE-n,III,combined}}} = \frac{1}{R_{\text{t,TE-n,III}}} + \frac{1}{R_{\text{t,oc}}/2}
$$

These combined resistors are in series with the region II resistors. Therefore, the total TE leg resistors for both p-type and n-type legs become:

$$
  \begin{aligned}
    R_{\text{t,TE-p,combined}} &= R_{\text{t,TE-p,I,combined}} + R_{\text{t,TE-p,II}}\\ &+ R_{\text{t,TE-p,III,combined}}
  \end{aligned}
$$

$$
  \begin{aligned}
    R_{\text{t,TE-n,combined}} &= R_{\text{t,TE-n,I,combined}} + R_{\text{t,TE-n,II}}\\ & + R_{\text{t,TE-n,III,combined}}
  \end{aligned}
$$

These combined resistances are in parallel with the azimuthal resistances. Therefore,$$
  \begin{aligned}
    \frac{1}{R_{\text{TE,combined}}}  &= \frac{1}{R_{\text{az}}/2} + \frac{1}{R_{\text{t,TE-p,combined}}} \\ &+ \frac{1}{R_{\text{az}}} + \frac{1}{R_{\text{t,TE-n,combined}}} +  \frac{1}{R_{\text{az}}/2}
  \end{aligned}
$$

Finally, the above resistor is in series with the radial insulator resistor, and the total lumped resistance between two adjacent interconnect nodes becomes,$$
  R_{\text{TEC}} = R_{\text{TE,combined}} + R_{\text{is}}
$$


---
# 5. Exploration


# 6. References
