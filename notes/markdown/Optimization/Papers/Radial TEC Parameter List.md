
We have a large amount of parameters to deal with.


# 1 Parameters

| Name                         | Default     | Symbol         | Upper Bound | Lower Bound | Units          |
| ---------------------------- | ----------- | -------------- | ----------- | ----------- | -------------- |
| Wedge Angle                  | 30$\degree$ | $\theta$       |             |             | rad            |
| Interconnect azimuthal angle | 5$\degree$  | $\beta_{ic,i}$ |             | $\theta$    | rad            |
| Outerconnect azimuthal angle | 5$\degree$  | $\beta_{oc,i}$ |             | $\theta$    | rad            |
| TE Length                    |             | $L_i$          |             |             | $\mu m$        |
| Interconnect radial length   | 50          | $w_{ic,i}$     | 1           | $L_i/2$     | $\mu m$        |
| Interconnect thickness       | 20          | $t_{ic,i}$     | 1           | $t$         | $\mu m$        |
| Outerconnect radial length   | 50          | $w_{oc,i}$     | 1           | $L_i/2$     | $\mu m$        |
| Outerconnect thickness       | 20          | $t_{oc,i}$     | 1           | $t$         | $\mu m$        |
| Inter-Stage insulation width | 40          | $w_{is}$       | 1           | $L_i/2$     | $\mu m$        |
| Azimuthal insulation width   | 20          | $w_{az}$       | 1           | 100         | $\mu m$        |
| Thermal TSV radius           | 10          | $R_{TSV,i}$    |             |             | $\mu m$        |
| Thermal TSV pitch            | 20          | $P_{TSV,i}$    |             |             | $\mu m$        |
| SOI layer thickness          | 100         | $t_{SOI}$      |             |             | $\mu m$        |
| No of TSV                    | 1           | $N_{TSV,i}$    |             |             | $\mathbb N^+$  |
| Inner Cylinder Radius        | 1000        | $R_{cyl}$      |             |             | $\mu m$        |
| Chip layer thickness         | 50          | $t_{chip}$     |             |             | $\mu m$        |
| Radial expansion factor      |             | $k_r$          | 0.1         | 10          | $\mathbb{R}^+$ |
| Stage 1 TEC length           | 1000        | $L_1$          | 10          | 100         | $\mu m$        |
| Chip width and length        | 10000       | $w_{chip}$     |             |             | $\mu m$        |
| TE thickness                 | 50          | $t$            |             |             | $\mu m$        |
| No of Stages                 | 3           | $N_s$          | 1           |             | $\mathbb N^+$  |
| TEC base radius              |             | $R_{base}$     |             |             |                |

# 2 Constant Parameters

| Name                 | Description | Symbol    | Upper Bound | Lower Bound | Units   |
| -------------------- | ----------- | --------- | ----------- | ----------- | ------- |
| TSV radial clearance |             | $g_{rad}$ |             |             | $\mu m$ |
# 3 Descriptions

## 3.1 Radial Expansion Factor
This is the ratio of the lengths between consecutive stages. i.e.
$$\frac{L_{i+1}}{L_i}=k_r$$
Since we fix the length of the first segment as $L_1$, we can see that,

| $i$ | $r_{start,i}$ | $L_i$  | $r_{end,i}$      |
| --- | ------------- | ------ | ---------------- |
| 1   | $R_{cyl}$     | $L_1$  | $R_{cyl}+L_1$    |
| 2   | $r_{end,1}$   | $kL_1$ | $r_{end,1}+kL_1$ |
Therefore,
$$L_i=k^{(i-1)}L_1$$

## 3.2 TEC Base Dimensions
Since we take the chip as a $w_{chip}\times w_{chip}$ square, and the TEC needs to cover the entire thing, the base would be a circle that covers that entire squared. therefore its diameter will be,
$$D_{base}=\sqrt{2}.w_{chip}$$
Therefore its radius,
$$R_{base}=\frac{1}{\sqrt 2}w_{chip}$$
# 4 Constraints

## 4.1 Length and No. of stages Relationship

The length of all the stages and the radial insulators with the center cylinder should add up to the radius of the base.
$$R_{cyl}+\Sigma_{i=1}^N L_i=R_{base}$$
Substituting,
$$R_{Cyl}+\Sigma_{i=1}^Nk^{(i-1)}L_1=\frac{1}{\sqrt 2}w_{chip}$$
**Dimensional Reduction**  
Treat $k_r$​ as the optimization variable (Design Variable) and calculate $L_1$​ analytically inside your objective function.

**Step A: Calculate available length for TE material**  
The radial stack consists of $N_s$​ TEC stages and $(N_s+1)$ insulation gaps ($w_{is}$​).  
$$ L_{total\;active} = (R_{base} - R_{cyl}) - (N_s + 1)w_{is} $$

**Step B: Constraint Equation**  
The sum of geometric stages must equal this active length:  
$$ \sum_{i=1}^{N_s} L_i = L_1 \sum_{i=0}^{N_s-1} k_r^i = L_1 \left( \frac{1 - k_r^{N_s}}{1 - k_r} \right) = L_{total\;active} $$

**Step C: Forward Calculation**  
In your code, let the optimizer pick $k_r$​. Then immediately calculate $L_1$​:  
$$ L_1 = L_{total\;active} \cdot \left( \frac{1 - k_r}{1 - k_r^{N_s}} \right) $$  
_Note: Handle the case where $k_r=1$ separately, where $L_1 = L_{total\;active} / N_s$)_

### 3. Parametrizing Sub-Geometries (The "Ratio" Method)
Parameters like $w_{ic}$ (interconnect length) that are constrained by $L_i$ (e.g., $w_{ic} < L_i/2$). If the optimizer shrinks $L_i$, a fixed $w_{ic}$ might become invalid (negative leg length). 

**Solution:** Optimize **Ratios**, not absolute lengths. Define a new set of design variables \alpha: 
1. **Interconnect Ratio ($\alpha_{ic}$):** 
	Instead of $w_{ic}$ (microns), optimize $\alpha_{ic} \in [0.05, 0.45]$, $w_{ic, i} = \alpha_{ic} \cdot L_i$ 
2. **Outerconnect Ratio ($\alpha_{oc}$):** Optimize $\alpha_{oc} \in [0.05, 0.45].$ $w_{oc, i} = \alpha_{oc} \cdot L_i$ 
3. **Fill Factor ($\alpha_{fill}$):** Instead of optimizing leg width directly, optimize the fraction of the arc occupied by TE material vs. Azimuthal Insulation. $w_{az} \approx (1 - \alpha_{fill}) \cdot (\text{Arc Length})$ 
 
**Why this works:** No matter how small $L_i$ gets during the optimization of $k_r$, the interconnects will scale down proportionally, guaranteeing a valid geometry.

# 5 Materials

## 5.1 Materials based on Purpose

|**Name**|**Use**|**Default**|
|---|---|---|
|Thermoelectric|TE leg material|$\text{Bi}_2\text{Te}_3$|
|Electrical Connections|Interconnects, Outerconnects|$\text{Cu}$|
|Thermal TSVs|Vertical heat conduits|$\text{Cu}$|
|Radial Insulators|Inter-Stage insulation, support|$\text{AlN}$|
|**Azimuthal Insulators**|Gap between P/N legs|**$\text{SiO}_2$ (Silica)**|
|Chip Surface|Substrate/Buffering Layer|$\text{Si}$|
|Middle Cylindrical Heat Spreader|Central core material|$\text{Si}$|
|Added: Substrate Insulator|SOI Buried Oxide (BOX)|$\text{SiO}_2$|

## 5.2 Default Materials

| **Name**                            | **Notation**             | **Thermal Conductivity (κ) [W/(m·K)]** | **Electrical Resistivity (ρ) [Ω⋅m]** |
| ----------------------------------- | ------------------------ | -------------------------------------- | ------------------------------------ |
| **bismuth telluride**               | $\text{Bi}_2\text{Te}_3$ | **$\sim 1.2$** [2]                     | **$\sim 1.0 \times 10^{-5}$** [2]    |
| **Copper**                          | $\text{Cu}$              | **$\sim 400$** [3]                     | **$\sim 1.7 \times 10^{-8}$** [3]    |
| **Silicon**                         | $\text{Si}$              | **$\sim 150$** [4]                     | **$\sim 10^{-2}$ to $10^5$** [4]     |
| **Aluminium Nitride**               | $\text{AlN}$             | **$\sim 170$** [5]                     | **$> 10^{10}$** [5]                  |
| **Aluminum Oxide (Alumina)**        | $\text{Al}_2\text{O}_3$  | **$\sim 30$** [6]                      | **$> 10^{12}$** [6]                  |
| **Added: Silicon Dioxide (Silica)** | $\text{SiO}_2$           | **$\sim 1.4$** [7]                     | **$> 10^{14}$** [7]                  |

This analysis provides suitable bounds for the azimuthal insulation width and completes the material properties table for your radial thermoelectric cooler (TEC) design.

## 5.3 Sources

[1] M. H. F. Al-Ajaj and J. M. K. Al-Ani, "Parametric analysis and optimal design for two-stage thermoelectric cooler," _J. Therm. Sci. Eng. Appl._, vol. 14, no. 1, 011017, Feb. 2022.

[2] Y. Lin, C. Lee, and Y. Lin, "Optimization design and performance evaluation of thin-film thermoelectric micro-cooler," _Sensors_, vol. 18, no. 2, p. 556, Feb. 2018.

[3] G. P. Wulff, _The Structure and Properties of Materials, Vol. III: Mechanical Behavior_. New York, NY, USA: John Wiley & Sons, 1965, pp. 293-294.

[4] R. E. Jones, B. P. Linder, and P. E. Johnson, _Materials for Semiconductor Devices_. New York, NY, USA: Academic Press, 2016, pp. 110-112.

[5] G. A. Slack, R. A. Tanzilli, R. O. Pohl, and J. W. Vandersande, "The intrinsic thermal conductivity of $\text{AlN}$," _J. Phys. Chem. Solids_, vol. 48, no. 7, pp. 641-647, 1987.

[6] W. D. Kingery, H. K. Bowen, and D. R. Uhlmann, _Introduction to Ceramics_, 2nd ed. New York, NY, USA: John Wiley & Sons, 1976, pp. 949-952.

[7] M. G. V. V. N. B. Prasad, P. V. S. S. S. N. V. Prasad, and K. C. B. C. Rao, "A review on thermal conductivity of silicon dioxide thin films used in electronic devices," _J. Therm. Anal. Calorim._, vol. 145, pp. 325-334, 2021.



## 5.4 Paper 1
[Parameter analysis and optimal design for two-stage thermoelectric cooler](<file:///E:\Semester 7\ME4311 - MicroNano Electro Mechanical Systems and Nanotechnology\Project\Papers\Parameter analysis and optimal design for two-stage thermoelectric cooler.pdf>)

This paper uses multi-physics models alongside conjugate gradient related optimization to optimize the 2 stage TEC. This comes under [[Advanced Optimization]]

