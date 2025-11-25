
# 1 TSV Thermal Resistance

## 1.1 TSV Resistance Calculations

Resistance of constant cross section TSV can be calculated using the basic formula,
$$R_{TSV}=\rho\frac{L}{A}$$
Where the area is $A=\pi R_{TSV}^2$ and the $L$ is equal to the thickness of the SOI layer.
$$R_{TSV}=\rho_{TSV}\frac{t_{SOI}}{\pi R_{TSV}^2}$$
Since all the TSV's are in parallel to each other, their combined resistance is,
$$R_{TSV,\mathrm{total}}=\frac{R_{TSV}}{N_\mathrm{TSV}}$$
We can find the TSV count by, multiplying per-row count by row count. 
$$N_\mathrm{TSV}= N_{TSV,\mathrm{rows}}\times N_{TSV,\mathrm{per-row}} $$
The TSV row count would be limited by its radius and the width of the copper interconnect.
$$N_{TSV,rows}=\left\lfloor \frac{w_{ic}}{2R_{TSV}+g_{rad}}\right\rfloor$$
Where $g_{rad}$ is the radial clearence. Per row count can be obtained by dividing the circumference from TSV Pitch and flooring to the nearest integer,
$$N_{TSV,\mathrm{per-row}}=\left\lfloor\frac{r\beta_{ic}}{P_{TSV}} \right\rfloor$$
> For now let's take $N_{TSV,\mathrm{rows}}=1$

$$R_{TSV,\mathrm{total}}=\frac{\rho_{TSV}\frac{t_{SOI}}{\pi R_{TSV}^2}}{\left\lfloor \frac{w_{ic}}{2R_{TSV}+g_{rad}}\right\rfloor\times \left\lfloor\frac{r\beta_{ic}}{P_{TSV}} \right\rfloor}$$
## 1.2 Combined Vertical Resistance

We consider the middle of the copper interconnect as the cold junction. Between it and the TSV endpoint there is the dielectric electrical insulation layer and half of the thickness of the copper interconnect layer. Then the total vertical resistance is given by,
$$R_v=R_{TSV}+R_{\mathrm{dielectric}}+R_{k,ic}$$
> For now we can assume that the thermal resistance of dielectric layer and the copper interconnect is zero. This can be attributed very thin (~500nm) dielectric layer and the higher conductivity of copper.


# 2 Insulator Layer Contributions
## 2.1 Radial Insulator Layer Resistance

Since there is a radial insulator layer between TEC modules, we need to absorb this resistance into both TECs connecting it, in order to form the resistance network. As we have done with the gaps in [[Variable Cross section area TEC]], we divide this into 2 segments and and connect it at start and the end of the TEC resistor.

The issue is that when we do that, we assume that the nodes are these weak thermal conducting ceramic layers. The nodes must in uniform temperatures, and therefore must have higher conductance. Because of this, let's consider the copper interconnect as the node. Then, the resistance of this insulator layer should be added to the previous stage's thermal resistance.

Now we calculate the resistance of this layer,
$$R_{is}= \int_{r_1}^{r_2} \frac{1}{\kappa_{is} A(r)} dr$$
The area is $A(r)_{is}=r\theta t$ and integration limits are for end of the TEC unit, since we add end insulator resistance to the TEC element.
$$R_{is} = \int_{r_1+L-w_{is}}^{r_1+L} \frac{1}{\kappa_{is} t\theta r} dr$$
$$R_{is}= \frac{1}{\kappa_{is} t \theta}ln\left(\frac{r_1+L}{r_1+L-w_{is}}\right) $$
Then the total Resistance of TEC is series addition of this resistance and the resistance from [[Variable Cross section area TEC]],
$$R_{eff,series}=R_{is}+\frac{1}{K_{total,TE}}$$
$$R_{eff,series}=\frac{1}{\kappa_{is} t \theta}ln\left(\frac{r_1+L}{r_1+L-w_{is}}\right)+\left(\frac{G}{\kappa_P+\kappa_N}\right)$$
Where G is the geometric factor for the TEC.

## 2.2 Azimuthal Insulator layer back conduction

This one can be added to the back conduction term of the TEC element. $K(T_h-T_c)$. We can calculate a Effective conductance from both TEC legs and the insulator layer.

Since this is characterized by the arc length, we can directly derive that $A(r)=w_{az}t$. Noting that this doesn't change as the radius changes, we can directly write the conductance as,

$$K_{az}=\kappa_{az}\frac{w_{az}t}{L}$$
Since the insulator and the conductor are in parallel to each other, we can directly add this term  
$$K_{total}=K_{eff,series}+K_{az}$$
## 2.3 Total Conductance of TEC

The combined conductance is as follows.
$$K_{global} = \underbrace{\left( \frac{1}{ \frac{1}{K_{total,TE}} + R_{is} } \right)}_{\text{Series Effect (Dam)}} + \underbrace{K_{az}}_{\text{Parallel Effect (Leak)}}$$
This term can be directly used inside the matrix.
# 3 Node 0 treatment

Node 0 acts as a **Boundary Condition**. The Physical Reality: At the exact center ($r=0$ to $r=r_0$), the "Chip Layer" and "TEC Layer" are mechanically and thermally merged by the heat spreader/via bundle. Therefore, they share the same temperature. The Implementation:

In the Block Matrix formulation (defined in [[System of Equations (TEC Preliminary Optimization)]]), you have two distinct vectors $\mathbf{T_{Si}}$ and $\mathbf{T_{c}}$. Node 0 effectively collapses these two layers into one point.

## 3.1 Equation

Add a specific row for $T_0$ at the top of your matrix.

- Equation for Node 0: Heat Generation = Heat leaving to Si Ring 1 + Heat leaving to TEC Ring 1.
$$Q_{gen,0} = \frac{T_0 - T_{Si,1}}{R_{Si,0\to1}} + \frac{T_0 - T_{c,1}}{R_{TEC,0\to1}}$$
$$\left(\frac{1}{R_{Si,0\to 1}}+\frac{1}{R_{TEC,0\to 1}}\right)T_0-\frac{1}{R_{Si,0\to 1}}T_{Si,1}-\frac{1}{R_{c,0\to 1}}T_{c,1}=Q_{gen,0}$$
- **Matrix Entry:**
    - Diagonal element $(0,0)$: $\left(\frac{1}{R_{Si,0\to1}} + \frac{1}{R_{TEC,0\to1}}\right)$
        
    - Off-diagonal $(0, Si_1)$: $-\frac{1}{R_{Si,0\to1}}$
        
    - Off-diagonal $(0, TEC_1)$: $-\frac{1}{R_{TEC,0\to1}}$
        
    - RHS Vector: $Q_{gen,0}$ (The portion of heat generated in the center cylinder).

## 3.2 Resistance Values

We need to find the resistance between Silicon cylinder and the TEC, as well as the other chip layer elements. We consider the node 0 as the center for the chip / TEC layer.

For a cylindrical wedge generating heat uniformly, the thermal resistance from the peak temperature (at the center) to the outer surface is a derived analytic result.
$$R_{int} = \frac{1}{2 \cdot \theta \cdot k \cdot t}$$
**Why?** The standard solution for a cylinder with volumetric generation is $T_{center} - T_{surf} = \frac{Q}{4\pi k t}$. Adjusting for the wedge angle $\theta$ (where the full circle is $2\pi$), the factor becomes $2\theta$. F

### 3.2.1 For Chip Layer
We can substitute values for the above equation.
$$R_{int,chip} = \frac{1}{2 \theta k_{Si} t_{chip}}$$
And this is directly equal to the $R_{Si,o\to 1}$.
### 3.2.2 For TEC Layer

$$R_{int,TEC} = \frac{1}{2 \theta k_{Si} t}$$
The total Resistance is the summation of resistance from the cylinder and the radial insulator layer.
$$R_{TEC, 0\to 1} = R_{int} + R_{is}$$
Since we know that,
$$R_{is}= \frac{1}{\kappa_{is} t \theta}ln\left(\frac{R_{cyl}+w_{is}}{R_{cyl}}\right) $$

We can write,
$$R_{TEC, 0\to 1} =\frac{1}{t\theta}\left[\frac{1}{2 k_{Si}}+\frac{1}{\kappa_{is}}ln\left(\frac{R_{cyl}+w_{is}}{R_{cyl}}\right)\right]$$

# Modelling Node N



# 4 Modelling a Representative Wedge

instead of combining resistors in parallel (which creates tiny resistance values and massive total currents), we should solve for **one single wedge** and scale the inputs.

1. **Resistance Network:** Calculate $R_{lat}$​,$R_{vert}$​, and $R_{TEC}$​ for a _single_ wedge using the formulas  derived in [[Variable Cross section area TEC]].
    
2. **Current:** Use the current​ flowing through a _single_ wedge (or leg pair).
    
3. **Heat Load Scaling:**
    
    - Total Chip Power: $Q_{total}​=400W$ (for example).
        
    - Total Wedges: $N_{wedges}​=360∘/\theta_{wedge​}$.
        
    - **Input to Solver:** Apply a heat load of $Q_{in}​=Q_{total}​/N_{wedges}$​ to the representative wedge.
        

**Why this is better:** It keeps your numbers relatable (e.g., 2 Amps per wedge rather than 2000 Amps total) and makes debugging easier.

# 5 Chip Layer Calculations

To calculate the physics of the "Chip Layer" (Layer 1), we cannot treat it as a 2D surface; it must have a volume. We must add the thickness of the silicon die to your parameter list. In modern 3D stacked ICs, wafers are often thinned significantly.

- **Symbol:** $t_{chip}$
- **Typical Range:** $50 \mu m$ to $300 \mu m$.
- **Material:** Silicon ($k_{Si} \approx 130-148 \text{ W}/m\cdot K$).

Without this, the lateral resistance is mathematically infinite (area = 0).
## 5.1 Lateral Resistance
This is the resistance to heat flowing **radially** through the silicon from the center of Node $i$ to the center of Node $i+1$.

**The Geometry:**
- **Flow Direction:** Radial ($r$).
- Cross-Sectional Area $A(r)$: This is the area perpendicular to the flow.
$$A_{cross}(r) = \text{Arc Length} \times \text{Thickness} = (r \theta) \times t_{chip}$$
The Integration:
Just like for the TEC legs in [[Variable Cross section area TEC]] the resistance of a wedge-shaped block is logarithmic. The resistance between two radii $r_a$ and $r_b$ in the chip is:
$$R = \int_{r_a}^{r_b} \frac{dr}{k_{Si} A_{cross}(r)} = \int_{r_a}^{r_b} \frac{dr}{k_{Si} t_{chip} \theta r}$$
$$R = \frac{1}{k_{Si} t_{chip} \theta} \ln\left(\frac{r_b}{r_a}\right)$$
## 5.2 Heat Generation Area

Inside the chip layer $Q_{gen} = \text{Flux} \times \text{Area}$.
- **The Physics:** Heat is generated on the active face of the transistor layer. This corresponds to the **Top Surface Area** (floor area) of the wedge segment.
- **The Geometry:** The Chip Layer is continuous. Unlike the TEC layer, it does **not** have gaps ($w_{is}$) or separate legs. It covers the entire wedge angle $\theta$.

The Formula For a specific Stage $i$ spanning from radius $r_{start}$ to $r_{end}$:
$$A_{top, i} = \text{Area of Sector} = \frac{\theta}{2} (r_{end}^2 - r_{start}^2)$$
$$Q_{gen, i} = q''_{flux} \times A_{top, i}$$
For the internal cylinder, the generation term comes from area of the circle cone,
$$Q_{gen,0}=q''_{flux}\times \frac{\theta}{2} R_{cyl}^2$$

