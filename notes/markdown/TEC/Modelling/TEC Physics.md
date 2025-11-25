#theories #physical_modelling #MEMS

- Semiconductors gives superior thermoelectirc performance compared to traditional materials.
- Low dimensional materials gives higher performance.

## Structure

![[Pasted image 20250928153908.png]]


## Basic Relations
Source: [The on-chip thermoelectric cooler](<file:///E:\Semester 7\ME4311 - MicroNano Electro Mechanical Systems and Nanotechnology\Project\Papers\The on-chip thermoelectric cooler.pdf>)

**Assumption:** Ideal thermal interface between cold and hot sides, i.e. no [[Thermal Contact Resistance]]. The following formula is equivalent to the first formula in [[Multi-Stage TEC Modeling]]. Net cooling power,
$$Q_{\mathrm{c}}\,=\,\bigl(S_{\mathrm{p}}\,-\,S_{\mathrm{n}}\bigr)I T_{\mathrm{c}}\,\,-\,K(T_{\mathrm{h}}\,-\,T_{\mathrm{c}})\,\,-\,\,{\frac{1}{2}}I^{2}R$$
$S_p$ and $S_n$ are setback coefficient of p and n type TEC legs. $K$ is the _overall thermal conductance_ and $R$ is the _overall electrical resistance_ and they can be represented as,
$$K=\kappa_{\mathrm{p}}{\frac{L_{\mathrm{p}}}{A_{\mathrm{p}}}}+\kappa_{\mathrm{n}}{\frac{L_{\mathrm{n}}}{A_{\mathrm{n}}}}$$
$$R=\rho_{\mathrm{p}}{\frac{L_{\mathrm{p}}}{A_{\mathrm{p}}}}+\rho_{\mathrm{n}}{\frac{L_{\mathrm{n}}}{A_{\mathrm{n}}}}$$
where $\kappa,\ \rho,\ L\ \mathrm{and}\ A$ are thermal conductivity, electrical resistivity, thickness, and cross-section area. For variable cross section area TEC  modules, [[Variable Cross section area TEC]] The power consumption of the TEC unit can be calculated using,
$$P=VI=\left(S_{\mathrm{P}}~-~S_{\mathrm{n}}\right)(T_{\mathrm{h}}~-~T_{\mathrm{c}})I+I^{2}R$$
Using the above equations, we can find coefficient of performance [[COP]] of the TEC unit, which is a standard heat pump measuring parameter.
$$\mathbf{COP}=\frac{\text{net cooling power}}{\text{power consumption}}$$
$$\mathrm{COP}={\frac{Q_{\mathrm{c}}}{P}}={\frac{\left({\mathrm{S}_{\mathrm{p}}}\,-\,S_{\mathrm{n}}\right)I T_{\mathrm{c}}\,-\,K(T_{\mathrm{h}}\,-\,T_{\mathrm{c}})\,-\,{\frac{\mathrm{i}}{2}}I^{2}R}{\left(S_{\mathrm{p}}\,-\,S_{\mathrm{n}}\right)(T_{\mathrm{h}}\,-\,T_{\mathrm{c}})I+I^{2}R}}$$
Differentiating 1st equation and setting it to 0 we get the current that gives maximum cooling power,
$$\left({\frac{d Q_{\mathrm{c}}}{d I}}\right)_{\mathrm{opt}}=0\quad\rightarrow\quad I_{\mathrm{opt}}={\frac{\left(S_{\mathrm{p}}\,-\,S_{\mathrm{n}}\right)T_{\mathrm{c}}}{R}}$$
$Q=Q_{max}$ when $I=I_{opt}$
$$Q_{\mathrm{max}}=\frac{\left(S_{\mathrm{p}}\;-\;S_{\mathrm{n}}\right)^{2}T_{\mathrm{c}}^{2}}{2R}\;-\;K(T_{\mathrm{h}}\;-\;T_{\mathrm{c}})$$
You can substitute this to the first equation with $Q_c=0$ and it will provide the maximum temperature difference that can be achieved. Physically, the condition **$Q_c = 0$** means the thermoelectric cooler is **no longer absorbing heat from the cold object or MEMS device**.

In this specific state:

- The **Peltier cooling rate** ($\alpha_S T_c J_e$) is exactly equal to the sum of the heat load coming from **Joule heating** ($\frac{1}{2} J_e^2 R_e$) and the heat conducted through the elements via **thermal conduction** ($\frac{1}{R_k} (T_h - T_c)$).
    
- The cold junction temperature ($T_c$) has reached its **lowest possible steady-state temperature** for the given hot junction temperature ($T_h$) and current ($J_e$).
    
- This is the point where the **maximum temperature difference** ($\Delta T_{\text{max}}$) is achieved. The device is running, but its net cooling capacity is zero, meaning it can maintain that low temperature but cannot cool the attached device any further or reject additional heat load from it.

refer [thermoelectric coolers chapter](<file:///E:\Semester 7\ME4311 - MicroNano Electro Mechanical Systems and Nanotechnology\Project\Papers\thermoelectric coolers chapter.pdf>) for more analysis such as optimal voltage calculation.

The book [Introduction to Thermoelectricity](<file:///E:\Semester 7\ME4311 - MicroNano Electro Mechanical Systems and Nanotechnology\Project\Proposal\Introduction to Thermoelectricity.pdf>) is a great source for Physical Modelling.

---
# Thermoelectric figure of merit
For any TE material, three properties fight each other:

1. Seebeck coefficient $\alpha_S$​:  Higher is better (stronger Peltier pumping).
2. Electrical resistivity $\rho_e$​:  Lower is better (less Joule heating).
3. Thermal conductivity $k$:  Lower is better (less heat leaks from hot to cold).

Ideal TE material:  
**High Seebeck, low resistivity, low thermal conductivity.**
The classical material figure of merit is
$$Z=\frac{\alpha^2}{\rho_e k}$$
Units: 1/K.
Higher $Z$ → theoretically higher efficiency for both cooling and power generation. If you also multiply it by absolute temperature T, you get the dimensionless ZT commonly used.

For a semiconductor TE element, you combine the Z of both p side and n side.
$$Z_{e,int} \equiv \frac{\alpha_{\rm S}^2}{\left[ (k\rho_e)_p^{1/2} + (k\rho_e)_n^{1/2} \right]^2}$$
This is called as **intrinsic** value, because it solely dependent on the material.

>[!Intrinsic thermoelectric figure of merit]
>This cannot be optimized by modifying the geometry. A material with highest $Z$ has to be selected, and the **extrinsic** value of the $Z$ is dependent on the geometric parameters.

Extrinsic figure of merit can be calculated as,
$$Z_{\mathrm{e,ext}}=\frac{\alpha_S^2}{R_e/R_k}$$
where the thermal resistance $R_k$ , and the electrical resistance $R_e$. These are dependent on the material properties. For our case, check [[Variable Cross section area TEC]].