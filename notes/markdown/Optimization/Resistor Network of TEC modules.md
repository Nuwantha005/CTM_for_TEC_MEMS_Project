Get the design details from [[Basic Radial TEC Design]]. The following modelling is done according to that criteria. 
# 1 At System Level

To include condition that the first few layers would be pumping heat to the TEC layers and others stay insulated, we need to discretize the chip layer and connect them as the resistors. With some resistors connected to the upper TEC units through TSV resistor and others staying the sane.

This way we can even model the direct connection between the chip and the micro channel, which can be used to offload some amount of heat to the micro channel directly.
## 1.1 The Two-Layer Resistance Network

Imagine two parallel rails of nodes.
- **Bottom Rail (The Chip):** Nodes $T_{Si,1}, T_{Si,2}, \dots, T_{Si,N}$. Heat is generated here.
- **Top Rail (The TEC Cold Side):** Nodes $T_{c,1}, T_{c,2}, \dots, T_{c,N}$. Heat is extracted here.
- **Vertical Connections:** Resistors connecting the two rails ($R_{v,i}$)

### Controlling Connected and Non connected Regions mathematically

We define the zones simply by changing the value of the vertical resistor $R_{v,i}$:
- **Zone 1 (Evaporator):** $R_{v,i} = R_{TSV}$ (Very low resistance). Heat flows easily from Chip $\rightarrow$ TEC.
- **Zone 2 (Pumping):** $R_{v,i} = R_{Insulation}$ (Very high resistance). Effectively $Q \approx 0$.

We apply the generated heat ($Q_{gen}$) to the **Chip Nodes**, not the TEC nodes. If Zone 2 has high vertical resistance, the heat generated underneath it is forced to conduct **laterally** through the silicon to reach either the TSVs in Zone 1 or the Microchannel edge. No heat is "lost"; the physics will force it to find a path.

---
## 1.2 The Governing Equations for 2 Layers

We need to write two balance equations for every stage $i$, for both layers and need to solve them simultaneously.
### 1.2.1 Equation A: The Chip Node ($T_{Si,i}$)

Heat enters from generation and lateral neighbors. Heat leaves vertically to the TEC.
$$\underbrace{Q_{\text{gen},i}}_{\text{Generated Here}} + \underbrace{\frac{T_{Si,i-1} - T_{Si,i}}{R_{\text{lat}}} + \frac{T_{Si,i+1} - T_{Si,i}}{R_{\text{lat}}}}_{\text{Lateral Conduction in Chip}} - \underbrace{\frac{T_{Si,i} - T_{c,i}}{R_{v,i}}}_{\text{Vertical Extraction to TEC}} = 0$$

- **$Q_{\text{gen},i}$:** The $Q_{\mathrm{chip}}$ distributed by area. (e.g., $Q_{\mathrm{chip}} \times \frac{Area_i}{TotalArea}$).
- **$R_{\text{lat}}$:** Thermal resistance of the silicon ring laterally.
- **$R_{v,i}$:** The "Switch":
    - If $i \le \text{Cutoff Stage}$: $R_{v,i} = \frac{L_{TSV}}{k_{Cu} A_{TSV}}$
    - If $i > \text{Cutoff Stage}$: $R_{v,i} = \frac{L_{Ins}}{k_{Ins} A_{Stage}}$
### 1.2.2 Equation B: The TEC Cold Node ($T_{c,i}$)
This is your original equation from [[TEC element level modelling]], but $Q_{chip}$ is replaced by the coupling term.
$$\underbrace{\frac{T_{Si,i} - T_{c,i}}{R_{v,i}}}_{\text{Input from Chip}} + \underbrace{K_{TEC}(T_{h,i} - T_{c,i})}_{\text{Back Conduction}} - \underbrace{\alpha I T_{c,i}}_{\text{Peltier}} + \frac{1}{2} I^2 R_{elec} = 0$$
---
## 1.3 Integrating the Microchannel Heatsink

The heatsink connects to the **perimeter** and this is a **Boundary Condition** for the thermal network.
### 1.3.1 Boundary 1: The TEC Edge

The hot side of the final TEC stage ($T_{h,N}$) connects to the water:
$$Q_{out, TEC} = \frac{T_{h,N} - T_{water}}{R_{\text{conv}}}$$
### 1.3.2 Boundary 2: The Chip Edge

The final Chip node ($T_{Si,N}$) connects to the water. You simply add this term to the **Equation A** for the last node $N$:
$$\dots + \underbrace{\frac{T_{Si,N} - T_{water}}{R_{\text{contact}} + R_{\text{conv}}}}_{\text{Direct Chip Cooling}} = 0$$
- If we **don't** connect the chip directly, set $R_{contact} = \infty$.
- If we **do** connect it, this term creates the "parallel path" you were worried about. The model will automatically calculate how much heat goes this way vs. up into the TEC.

>[!Integrating micro channel dynamics]
>We can model the micro channels in the preliminary stage using simple Nusselt's number relations from things learned in [[CFD Project]], and apply them to optimize things such as micro channel velocity. But this needs to be further studied.
>
>For now let's just treat micro channel as a constant temperature boundary.


# 2 At Stage Level (rings) - parallel resistors


# 3 At Element Level


# 4 On Radial Direction


## 4.1 Asymmetric Configuration



## 4.2 Non Asymmetric Configuration


