# 1 Thermal Network
## 1.1 Two Layer System

We must separate the **Heat Source Plane** (The Bottom Chip) from the **Heat Pumping Plane** (The TEC).

1. **The Chip Layer (Heat Source):** This generates the heat. It will naturally have a gradient (Hot Center → Cold Edge) due to the edge water cooling.
2. **The Isolation Layer:** Deposit a layer of low-k dielectric (polymer or thick $SiO_2$​) on top of the chip. This prevents the radial temperature differences in the TEC from interacting with the radial temperature differences in the chip.
3. **Vertical Thermal Vias (The Columns):** Use dense Copper Thermal [[TSV]]s (Through Silicon Vias) to connect the specific hot spots of the chip vertically _up_ to the cold side of the TEC rings.
4. **The TEC Layer:** This sits on top of the insulation. It pumps heat radially outward.
Then the temperature profile would look something like this.

> - **Chip Profile:** Center 80∘C → Edge 20∘C.
> - **TEC Profile:** Center 75∘C → Edge 25∘C.

- The TEC pulls heat _up_ from the center (80∘C source to 75∘C sink) and dumps it at the edge.
- Because the TEC creates a gradient (75∘C→25∘C) that roughly matches the chip (80∘C→20∘C), **back conduction is minimized** because there is little potential difference between the two layers.

We should use **Thermal [[TSV]]s only in the "Evaporator Zone" (Center + Inner Rings)** and strictly isolate the "Condenser Zone" (Outer Rings).

To understand why, we have to overlay the temperature profile of your Chip with the temperature profile of your TEC.

1. **The Chip Profile:** It is naturally **Hot at the Center** and **Cold at the Edge** (because you have water cooling at the perimeter).
2. **The TEC Profile:** It is actively **Cold at the Center** and **Hot at the Edge** (because it is pumping heat outward).

if we place TSVs at **all stages**, we create a connection at the outer edge where the TEC is hot (e.g., 25∘C) and the Chip is cold (e.g., 20∘C). Heat will flow **backwards** from the TEC into your chip, fighting the cooling system.

We should divide your radial array into two distinct zones.
### 1.1.1 Zone 1: The Evaporator Zone - first $E_N$ stages

- **Action:** **Use TSVs Here.**
- **Condition:** In this region, the TEC temperature is **lower** than the Chip temperature.
- **Physics:** By connecting TSVs to stages 1, 2, and maybe 3, you increase the "intake surface area." Instead of trying to suck 400W through a tiny pinhole at the absolute center, you suck it up through a wider disc. This lowers the heat flux density (W/cm2) at the interface, making the TEC much more efficient.

### 1.1.2 Zone 2: The Pumping Zone - last $(N-E_N)$ stages

- **Action:** **NO TSVs. Use Thick Insulation.**
- **Condition:** In this region, the TEC has pumped enough heat that it is now **hotter** than the chip underneath it.
- **Physics:** These stages are just "muscles" pushing the heat to the edge. If you add TSVs here, heat will leak back down into the chip (Back Conduction), creating a loop where the TEC pumps heat out, the TSV brings it back in, and the TEC pumps it out again. This wastes power and achieves nothing.

## 1.2 TSV Connections

### 1.2.1 Connection Location
These must be connected to the cold junction of each corresponding TEC. However, we cannot connect them to TEC legs because TEC materials have small thermal conductivity - which is required to avoid back conduction.

Therefore, we connect the TSVs to the respective TEC's cold junction. Then we can model this as a resistor connected to that node vertically. Since the TECs are in parallel, a lumped connector can be used.
### 1.2.2 TSV - TEC Junction
We cannot connect the TSV directly to the TEC module, because it the TSV material can be a electrical conductor and it will remove the electrical isolation of the thermoelectric elements and cause a short circuit. 

Because of this, we have 2 choices. Either connect the TSV to the ceramic insulator layer or to the interconnect copper wire.

#### 1.2.2.1 Ceramic - TSV Connection
This would be beneficial because the copper interconnect would connect to the electrically non conductive ceramic layer, avoiding the need of any insulator layers in between. But this ceramic layers is a poor conductor compared to TEC pumping power and copper interconnect's conductivity. This means this layer will act as a thermal dam. 

Therefore something like Aluminium Nitrite $(AIN)$ with ~$170\;W/m.K$ would be better compared to Silica $(Al_2O_3)$ with ~$25\;W/m.K$. Either way ==this layer needs to be as thin (radial) as possible== to increase the conductivity. Because of this, the cross section area of the connecting TSVs would be small, increasing the thermal resistance. 

#### 1.2.2.2 Interconnect - TSV Connection
This would be easier because the radial width of the interconnect can be increased (albeit reducing TEC section area by little) to accommodate large combined TSV cross section. The issue is that both copper layers are conductive and TSVs would conduct electricity and remove the electrical isolation of the TEC module and even cause a short circuit.

As a solution, a ==Dielectric Nano-layer== of electric insulators can be deposited between TSVs and the copper interconnecting, creating a landing pad for the heat. it would be **100−500nm of Diamond-Like Carbon (DLC) or Silicon Carbide (SiC)**.

>_Why?_ These are electrically insulating but so thin that heat jumps right through them with almost zero resistance. 

We should try to decrease the [[Thermal Spreading Resistance]] by increasing TSV count and getting the combined area closer to interconnect area. 

### 1.2.3 TSV Density
We can control this via 2 methods - no of TSVs in a row and no of rows for a given TSV set.

#### 1.2.3.1 No. of TSVs in a row
We can arrange TSVs with constant pitch between them. As the Stages goes away from the center, if the pitch remains the same, it will increase the TSV count but the TSV density will remain the same. We can choose to change this pitch stage by stage using some form of ratio - as we have done with such parameters to reduce the dimensions of the solution space.

#### 1.2.3.2 No. of rows for a TEC
This number would be limited by the width of the copper interconnect, as well as the radius of the TSV cross section. Let's simplify things for now and use 1 row, let's see what we can do later.

>[!TSV Arrangement]
>In Cocentric circular arcs centering the middle of the chip - a circular pattern of TSV circles arranged with constant pitch for a given radius.
## 1.3 TSV Geometry
A cylindrical geometry would be most suitable. While the TEC is wedge-shaped, we should **not** use a solid "wedge" or "curved bar" shaped TSV.

- **Reason:** **Thermal Stress & Cracking.** Copper (TSV material) expands significantly more than Silicon (Chip material). A long, solid metal bar acts like a wedge when heated, creating massive stress concentrations at the corners that will crack the silicon die during manufacturing or operation.
    
- **Solution:** Use standard **Circular TSVs** arranged in a dense, hexagonal packed array that forms an "arc" shape. This mimics the geometry of the wedge while relieving mechanical stress.
	
- **TSV Diameter:** Typically 5μm−20μm for aggressive density, or up to 50μm for standard processes.
    
- **Packing:** Hexagonal packing allows you to achieve ∼80% equivalent metal coverage of a solid bar without the fracture risk.

>Hexagonal pack wouldn't be needed for now. (multiple rows) As we are using a single row for simplicity


# 2 Electric Network

The condition $J_{e, \text{opt}} = \frac{\alpha_S T_c}{R_e}$ (which maximizes $\Delta T_{\text{max}}$ for $Q_c=0$) is strictly valid for a **single-stage cooler**. If the goal is to extract the maximum amount of heat $Q_c$ (Active Cooling Rate) at a fixed $\Delta T$, a different current may be optimal.

See [[TEC Physics]] for more details on optimal current.
## 2.1 Wiring Method

We have several methods available for the wiring of the TEC elements.

1) All elements are series connected $\longrightarrow \mathrm{Constnt\;Current: }\;I$   
2) Single Stage is series connected $\longrightarrow$ each stage has its own current.
3) Single radial segment is series connected $\longrightarrow$ each radial segment has its own current.
4) All Elements are parallel connected - ? #todo 
5) Single Stage is parallel connected - ? #todo 
6) Each element is separately connected $\longrightarrow$ Nightmare

Let's consider the most suitable 2 cases.
### 2.1.1 Constant Current for All

In a traditional multistage cooler, the stages are connected electrically in series, meaning the current through every stage is the **same** 
$$J_{e,1} = J_{e,2} = J_{e,3} = \dots = J_e$$
### 2.1.2 Constant Current for a Stage

We can add separate currents for each stage, but this might be difficult to optimize. the best decision is to add an ratio of currents between consecutive stages, **stage ratio $\mathbf{J_{e,i}/J_{e,i+1}}$** as optimization parameters. 

 **A. Current Varies Stage-by-Stage (Non-Series Connection):** This would require **complex parallel wiring** or separate power supplies for each stage to achieve an optimal current distribution, $J_{e,1} \ne J_{e,2} \ne J_{e,3}$, etc.
 
**The Stage Ratio $\mathbf{J_{e,i}/J_{e,i+1}}$** is the ratio of the current supplied to one stage to the current supplied to the next, adjacent stage. This ratio is optimized to ensure that the heat accumulated from the colder stages is efficiently pumped out by the warmer stages.

>For a multistage cooler, the primary design challenge is that each successive stage must pump **its own Peltier heat load plus the entire heat load accumulated from all stages colder than it.**

To handle this increasing heat load efficiently, two ratios are crucial and must be optimized together:

1. **Ratio of TE Couples ($\mathbf{N_i/N_{i+1}}$):** More couples are needed at the hotter stages (larger $N_i$) to handle the higher total heat load.

2. **Ratio of Electrical Current ($\mathbf{J_{e,i}/J_{e,i+1}}$):** The current in each stage might be adjusted so that the _Peltier cooling power_ ($\propto \alpha_S T J_e N$) of stage $i$ precisely matches the total heat load arriving from stage $i+1$ (the heat $\text{load } Q_{h,i+1}$ in Equation 7) plus all the heat leakage.

Optimizing the current ratio $J_{e,i}/J_{e,i+1}$ in addition to the couple ratio $N_i/N_{i+1}$ ensures that the cooler operates with the ==highest overall system efficiency== and achieves the largest $\Delta T$ possible for the entire structure, which is generally not achieved when every stage runs at its theoretical single-stage optimal current $J_{e, \text{opt}}$.

# 3 Geometry

Following is the Basic Structure of a single TEC element without TSVs.

![[image.png]]


the image shows a top view. the circle is the middle heat spreader which is at the center of the chip. note that this is a enlarged image used for analysis

1. Red arrow - Ceramic interconnect (Dark grey)
2. blue arrows - TEC material
	1. Blue - P material
	2. Red - N material
3. Green arrow (Pink material) - Copper wires connecting this TEC to next TEC on the same stage
4. Black arrow (Yellow material) - Copper interconnect between P and N in same TEC element.

![[image-1.png]]


# 4 Boundary Condition

we have 2 types of boundary conditions to be imposed here - a constant temperature one or conjugate transfer of heat to water in the micro channel. 

## 4.1 Modelling the BC

we have 2 ways to model this. we can set the last nodes of both layers to the same condition, and add another equation at the bottom like we did for the [[TEC Thermal Network#3 Node 0 treatment|Node 0]]. So the Silicon / Chip layer is touching the water or the boundary directly and is transferring heat to it. But this is not realistic because,

- **Packaging Constraints:** In realistic 3D-IC packages, the edge of the silicon die is often sealed with under-fill or surrounded by a mold compound which has very low thermal conductivity compared to silicon or copper. The die edge is rarely exposed directly to the coolant to prevent electrical shorts and corrosion.
    
- **Design Intent:** The purpose of the radial TEC is to elevate the heat flux from the logic dies and pump it **vertically** into the pumping layer, then **radially** outward to the heat sink ring. If you model the silicon edge as touching the water (Tw​), you create a "short circuit" where heat flows through the silicon substrate instead of the TEC. While this happens in reality (leakage), relying on it for the boundary condition defeats the purpose of simulating the TEC's effectiveness.
    
- **Heat Sink Attachment:** The microchannel heat sink is typically a copper or silicon manifold attached to the **hot side** of the TEC module at the periphery, not glued to the side of the active logic die.

So we can modify this like there is no node in the Silicon layer after the node under the last TEC element's cold junction. connecting them through the either TSV or weak SOI heat transfer. 

Now we need to model this behavior and modify the TEC equations in [[System of Equations (TEC Preliminary Optimization)#2 General Nodal Equations|this section]] to reflect the boundary conditions. First consider the Silicon layer. Since the silicon layer ends here and does not touch the water, there is no "Right Neighbor" $(T_{Si, N+1})$.

- **Vertical Term:** Remains the same: $\frac{T_{Si,N} - T_{c,N}}{R_{v,N}}$.
    
- **Input from Left:** Remains the same (conduction from $T_{Si, N-1})$.
    
- **Output to Right:** This entire term **vanishes**. Mathematically, $R_{lat, N} \to \infty$ or the conduction term is set to 0.
    
- **Source:** $Q_{gen, N}$ remains.
    
Modified Silicon Eq $(i=N)$:

$$\frac{T_{Si,N-1} - T_{Si,N}}{R_{lat, N-1}} - \frac{T_{Si,N} - T_{c,N}}{R_{v,N}} + Q_{gen,N} = 0$$

(Notice the term for $T_{Si, N+1}$ is completely gone).

For the TEC layer, the equation is dependent on the boundary condition imposed and it will be discussed on the following sections.
## 4.2 Constant Temperature BC

This BC is used for simplicity and only used in the beginning stages of the preliminary optimization where microchannel performance isnt a concern.
### 4.1.2 TEC Layer Equation

This node sits at the start of the last leg. The heat travels _through_ this last leg (Stage N) and rejects to the constant temperature boundary $T_w$.

We do **not** add a new equation. Instead, we substitute the boundary condition directly into the "Output to Right" term.

- **Current Term:** $K_N (T_{c, N+1} - T_{c, N})$
    
- **Substitution:** Replace the unknown next node $T_{c, N+1}$ with the fixed boundary temperature $T_w$.

Modified TEC Eq (i=N):

$$\underbrace{\frac{T_{Si,N} - T_{c,N}}{R_{v,N}}}_{\text{From Chip}} + \underbrace{\left[ \dots \right]_{N-1}}_{\text{From Left}} - \left[ S_{N} I_{N} T_{c,N} + K_{N}(\mathbf{T_w} - T_{c,N}) - \text{Joule}_{Cold, N} \right] = 0$$


## 4.3 Conjugate HT to the Water BC



