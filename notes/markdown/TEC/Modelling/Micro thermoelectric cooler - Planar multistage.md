---
Title: "Micro thermoelectric cooler: Planar multistage"
tags:
  - MEMS
  - research_paper
DOi: https://dx.doi.org/10.1016/j.ijheatmasstransfer.2008.10.014
number headings: auto, first-level 1, max 6, 1.1
---
[Micro thermoelectric cooler - Planar multistage](<file:///E:\Semester 7\ME4311 - MicroNano Electro Mechanical Systems and Nanotechnology\Project\Papers\Micro thermoelectric cooler - Planar multistage.pdf>)

# 1 Overview


# 2 Problem Definition
The 3D pyramid like model in figure a) is arranged to only the perimeter in figure b) and collapsed into two dimensions to create the planar system shown in figure c)

![[Pasted image 20251011171929.png|571x939]]

A close up view is shown in the following figure. Note that the bulk Silicon solid elements are thermally isolated and are called as **thermal islands** in the paper. The $SiO_2$  substrate is used as the bridge and structural support for the TE modules which pump heat between these islands. 

![[image.png]]

Apart from that, for structural purposes, there is a **tether** running down in the bottom of the whole module acting as a bridge, similar to Silicon substrate at the top. Glass is used for this and has a serpentine structure to reduce parasitic back conduction.

>[!note]
>This essentially works as thermally isolated rectangular concentric rings which are connected to each other by TEC modules
>

This is quite similar to our radial design, albeit the vertical conduction through the Silicon would happen and we need to consider whether isolation such as this would be required.
# 3 Theory

## 3.1 Individual Node level
The extrinsic thermoelectric figure of merit is,
$$Z_{\mathrm{e,ext}}=\frac{\alpha_S^2}{R_e/R_k}$$
Where the thermal resistance $R_k$ , and the electrical resistance $R_e$. These are,
$$\begin{split} R_{k} &= \sum_{i=1}^{N_{s}} R_{k,i} \\ R_{e} &= \sum_{i=1}^{N_{s}} (R_{e,te,i} + R_{e,wire,i}) \\ \frac{1}{R_{k,i}} &= \frac{1}{R_{k,te,i}} + \frac{1}{R_{k,sub,i}} + \frac{1}{R_{k,tether,i}} \\ \frac{1}{R_{k,te,i}} &= N_{i} \left[ \frac{1}{(R_{k,te})_{p}} + \frac{1}{(R_{k,te})_{n}} \right] = N_{i} \left[ \left( \frac{A_{k,te}k}{L_{te}} \right)_{p} + \left( \frac{A_{k,te}k}{L_{te}} \right)_{n} \right] \\ R_{e,te,i} &= N_{i} \left[ \left( \frac{\rho_{e}L_{te}}{A_{k,te}} \right)_{p} + \left( \frac{\rho_{e}L_{te}}{A_{k,te}} \right)_{n} \right]. \end{split}\tag{3}$$
In Eq. (3), $N_s$ is the number of the stage, $N_i$ is the number of TE couples in the $i^{th}$ stage, $R_{k,te,i}$ is the thermal resistance of the TE material at the $i^{th}$ stage, $R_{k,sub,i}$ the thermal resistance of substrate at the ith stage, $R_{k,tether,i}$ is the thermal resistance of glass tether at the ith stage, $R_{e,te,i}$ is the electrical resistance of TE material at the ith stage, and $R_{e,wire,i}$ is the electrical resistance of inter-connecting.

>[!Modification]
>We can do the summation another ways. This way, or we can consider TEC element wise. Say that the total electrical resistance of the TEC element is $R_{e,i}$. Then we know that there are 2 n and p type TE elements inside and set of wires connected to each of them. so we can say,
$$R_{e,i}=R_{e,TE,i}+R_{e,wire}$$
Then we can divide it into 2 TE elements.
$$R_{e,TE,i}=R_{e,TE-N,i}+R_{e,TE-P,i}$$
Now divide resistance from wire,
$$R_{e,wire}=R_{e,w_i}+R_{e,ic_i}$$
Where the 2 elements represent the wire segment before current TEC element and the wire segment connecting two P and N type TE segments together (interconnect). This way we can add contact resistance as well.

Modifications are not required for the thermal resistance. ~~However, we can improve the model by adding ==heat generation term for the wire==, but it would be at thermal modelling stage and not here. We might be able to add it **alongside the vertical heat flux**, calculated based on the wire dimensions.~~

We need to include all the generated terms inside the TEC element in the fundamental equation because it should have sufficient power to pump out all that. Adding generation term at the system level not good enough because if fails to enforce this at element level. See first section of [[TEC element level modelling]].

---
## 3.2 Resistor Network Arrangement


1. **Parallel Resistances (Within a Stage):** The thermal resistance of the substrate ($R_{k,sub}$), the tether ($R_{k,tether}$), and the connecting wires ($R_{k,wire}$) are all in **parallel** with the thermal resistance of the TE couples ($R_{k,te}$). All these paths connect the **Isothermal Stage Node ($T_i$)** to the adjacent stage's platform, thus determining the overall thermal isolation between two stages.
    
2. **Series Staging (End-to-End):** The thermal and electrical network for each stage is connected **in series** with the next stage to form the entire multistage cooler, running from the cold end ($T_c$) to the hot end ($T_h$).
### 3.2.1 **Parallel Arrangement (Stage-Level)**

The components that thermally bridge two adjacent platforms (or silicon islands) are in parallel. This combination determines the total **heat leakage** for that stage.

- **Goal:** The design aims for **high thermal isolation**, meaning the total effective thermal resistance of this parallel block should be as high as possible.
    
- **Challenge:** The presence of the substrate ($R_{k,sub}$) and tether ($R_{k,tether}$) means these parallel paths significantly reduce the total resistance, leading to lower cooling performance compared to an ideal cooler where only $R_{k,te}$ exists.
    
### 3.2.2 **Series Arrangement (System-Level)**

The individual stages are connected end-to-end to achieve a large total temperature difference ($\Delta T$).

- **Thermal Series:** Heat flows from the cold end, is pumped by the first stage to the second, then to the third, and so on, until the final stage rejects the accumulated heat to the ambient sink. Each stage handles the heat load accumulated from all previous (colder) stages.
    
- **Electrical Series:** The problem implies a single current ($J_e$) is driven through all stages. The stages are electrically connected in a series path to ensure the current passes through every TE couple, cascading the cooling effect.
    
The model is effectively a **series of thermal blocks**, where each block contains **parallel resistance paths** for heat flow and specific **Joule/Peltier energy generation/absorption nodes**.

> The tether, wire and substrate are parallel resistors. but the entire thing is in series with the TE modules alternatively.


# 4 Assumptions
- Thin film TE materials have same conductivity as bulk materials.
- Convection is neglected due to vacuum packaging.

# 5 Follow Up
1) [[phonon transport]] - reduced in the thin films due to large grain boundary scattering.
2) [[co-evaporated deposition technique]] - used to deposit TEC layers and the paper says it limits the film thickness to $10\mu m$.

# 6 References and Further Reading
1) Fabrication:
	1) `L.W. da Silva, M. Kaviany, C. Uher, Thermoelectric performance of films in the bismuth–tellurium and antimony–tellurium systems, J. Appl. Phys. 97 (2005) 114903.`
	2) `B. Huang, C. Lawrence, A. Gross, G. Hwang, N. Ghafouri, S. Lee, H. Kim, C. Li, C.Uher, K. Najafi, M. Kaviany, Low-temperature characterization and micro patterning of co-evaporated Bi2 Te 3/Sb2Te 3 films, J. Appl. Phys. (2008), in press.`
	3) `A. Gross, B. Huang, G. Hwang, C. Lawrence, N. Ghafouri, S.W. Lee, H. Kim, C. Uher, M. Kaviany, K. Najafi, A mutistage in-plane micro-thermoelectric cooler, 21st IEEE International Conference on Micro Electro Mechanical Systems, Tucson, 2008, pp. 1–2`
	
2) Contact Resistance
	`L.W. da Silva, M. Kaviany, Micro-thermoelectric cooler: interfacial effects onthermal and electrical transport, Int. J. Heat Mass Transfer 47 (10–11) (2004) 2417–2435`

3) 5 Stage TEC
	`A. Gross, B. Huang, G. Hwang, C. Lawrence, N. Ghafouri, S.W. Lee, H. Kim, C. Uher, M. Kaviany, K. Najafi, A mutistage in-plane micro-thermoelectric cooler, 21st IEEE International Conference on Micro Electro Mechanical Systems, Tucson, 2008, pp. 1–2`


