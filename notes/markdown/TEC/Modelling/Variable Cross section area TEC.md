We need to find expressions for the overall conductance and overall thermal resistance of TEC modules with variable cross section area to be used in our radial TEC modelling. There are 2 main approaches we can use.

# 1 TEC Legs
 $K$ is the _overall thermal conductance_ and $R$ is the _overall electrical resistance_ and they can be represented as,
$$K=\kappa_{\mathrm{p}}{\frac{A_{\mathrm{p}}}{L_{\mathrm{p}}}}+\kappa_{\mathrm{n}}{\frac{A_{\mathrm{n}}}{L_{\mathrm{n}}}}$$
$$R=\rho_{\mathrm{p}}{\frac{L_{\mathrm{p}}}{A_{\mathrm{p}}}}+\rho_{\mathrm{n}}{\frac{L_{\mathrm{n}}}{A_{\mathrm{n}}}}$$
where $\kappa,\ \rho,\ L\ \mathrm{and}\ A$ are thermal conductivity, electrical resistivity, thickness, and cross-section area.

---
## 1.1 General Formula

For variable geometry, we must integrate differential elements along the heat flow path:

**Thermal Conductance**: Resistances add in series
$$\frac{1}{K_{total}} = \int_{x_1}^{x_2} \frac{dx}{\kappa A(x)}$$

**Electrical Resistance**: Resistances add in series  
$$R_{total} = \int_{x_1}^{x_2} \rho \frac{dx}{A(x)}$$
For detailed derivation, check https://www.kimi.com/share/19a9b5c9-7942-8749-8000-0000b95732d4.

![[image-4.png]]
## 1.2 Specifics for Radial TEC Case
Assuming a wedge with:
- Constant angular width: $\theta = \frac{2\pi}{N}$ radians
- Constant thickness: $t$ (perpendicular to plane)
- Area varies as: $A(r) = t \cdot r \cdot \theta/2$
### 1.2.1 Thermal Conductance
General Formula,
$$\frac{1}{K} = \int_{r_1}^{r_2} \frac{dr}{\kappa A(x)}$$
Substituting for $A(x)$,
$$\frac{1}{K} = \int_{r_1}^{r_2} \frac{dr}{\kappa tr \theta/2}\tag{1}$$
$$\frac{1}{K} = \frac{2}{t\kappa\theta}\int_{r_1}^{r_2} \frac{dr}{r}$$
$$R_k=\frac{1}{K} = \frac{2}{t\kappa\theta}ln\left(\frac{r_2}{r_1}\right)\tag{2}$$
### 1.2.2 Electric Resistance
The same way,
$$R_e =\frac{2\rho}{t\theta}ln\left(\frac{r_2}{r_1}\right)\tag{3}$$
### 1.2.3 Combined for Both TE Legs
Now we do parallel resistance summation for thermals,
$$\frac{1}{R_{total,t}}=\frac{1}{R_{N,t}}+\frac{1}{R_{P,t}}$$
$$\frac{1}{R_{total,t}}=\frac{t\kappa_N\theta}{2ln\left(r_2/r_1\right)}+\frac{t\kappa_P\theta}{2ln\left(r_2/r_1\right)}$$
$$K_{total}=\frac{1}{R_{total,t}}=\frac{t\theta\left[\kappa_N+\kappa_P\right]}{2ln\left(r_2/r_1\right)}$$
Electrical resistance is positioned in series Connections,
$$R_{total,e}={R_{N,e}}+{R_{P,e}}$$
$$R_{total,e}=\frac{2\rho_N}{t\theta}ln\left(\frac{r_2}{r_1}\right)+\frac{2\rho_P}{t\theta}ln\left(\frac{r_2}{r_1}\right)$$
$$R_{total,e}=\frac{2[\rho_N+\rho_E]}{t\theta}ln\left(\frac{r_2}{r_1}\right)$$

---
## 1.3 Correction for the Wires

Divide the integral into 3 regions. From center going radially outward,
1) Region with the copper interconnect
2) Region with full occupancy of TE material
3) Region with the copper outerconnect

>[!Otuterconnect]
>is the connection wire between current TEC element in this stage and the previous and next elements in the same stage. Its shown in pink color in the figures.
### 1.3.1 Region with the copper interconnect
The multiplication of angle of the interconnect $\beta_{ic}$ and interconnect thickness $t_{ic}$ gives the reduced area. Then the new area is,
$$A(r)=\frac{t\theta-t_{ic}\beta_{ic}}{2}r$$
From eq. (1), $r_2$ becomes $r_1+w_{ic}$, where $w_{ic}$ is the width of the interconnect at $i^{th}$ stage 
$$\frac{1}{K} = \int_{r_1}^{r_1+w_{ic}} \frac{2}{\kappa(t\theta-t_{ic}\beta_{ic})r}dr$$

Since $\beta_{ic}$ is constant for a given stage, we can take it out from the integration and we get a modified version of eq. (2).
$$R_k=\frac{1}{K} = \frac{2}{\kappa(t\theta-t_{ic}\beta_{ic})}ln\left(\frac{r_1+w_{ic}}{r_1}\right)\tag{2A}$$
 In a similar manner we get the modified resistance for the region,
$$R_e =\frac{2\rho}{(t\theta-t_{ic}\beta_{ic})}ln\left(\frac{r_1+w_{ic}}{r_1}\right)\tag{3A}$$


![[image-1.png]]

![[image-2.png]]
### 1.3.2 Region with full occupancy of TE material

Regular integral for this region, we get modified versions of eq (2) and (3) with changed radius values. $r_1\rightarrow r_1+w_{ic}$ and $r_2\rightarrow r_1+L-w_{oc}$ where $L$ is the length of the current TEC element in radial direction and $w_{oc}$ is the radial width of the outerconnect copper element.
$$R_k=\frac{1}{K} = \frac{2}{t\kappa\theta}ln\left(\frac{r_1+L-w_{oc}}{r_1+w_{ic}}\right)\tag{2B}$$
$$R_e =\frac{2\rho}{t\theta}ln\left(\frac{r_1+L-w_{oc}}{r_1+w_{ic}}\right)\tag{3B}$$
### 1.3.3 Region with the copper outerconnect
The multiplication of angle of the interconnect $\beta_{oc}$ and interconnect thickness $t_{oc}$ gives the reduced area. Then the new area is,
$$A(r)=\frac{t\theta-t_{oc}\beta_{oc}}{2}r$$
From eq. (1), $r_1\rightarrow r_1+L-w_{oc}$ and $r_2\rightarrow r_1+L$ , where $w_{ic}$ is the width of the interconnect at $i^{th}$ stage 
$$\frac{1}{K} = \int_{r_1+L-w_{oc}}^{r_1+L} \frac{2}{\kappa (t\theta-t_{oc}\beta_{oc})r}dr$$
Since $\beta_{ic}$ is constant for a given stage, we can take it out from the integration and we get a modified version of eq. (2).
$$R_k=\frac{1}{K} = \frac{2}{\kappa(t\theta-t_{oc}\beta_{oc})}ln\left(\frac{r_1+L}{r_1+L-w_{oc}}\right)\tag{2C}$$
In a similar manner we get the modified resistance for the region,
$$R_e =\frac{2\rho}{(t\theta-t_{oc}\beta_{oc})}ln\left(\frac{r_1+L}{r_1+L-w_{oc}}\right)\tag{3C}$$

![[image-3.png]]


### 1.3.4 Total Bulk Properties for TE legs

We add together all the terms. For thermal resistance, $(2A+2B+2C)$ gives
$$R_k=\frac{1}{K} = \frac{2}{\kappa(t\theta-t_{ic}\beta_{ic})}ln\left(\frac{r_1+w_{ic}}{r_1}\right)+ \frac{2}{t\kappa\theta}ln\left(\frac{r_1+L-w_{oc}}{r_1+w_{ic}}\right)+\frac{2}{\kappa(t\theta-t_{oc}\beta_{oc})}ln\left(\frac{r_1+L}{r_1+L-w_{oc}}\right)$$

$$R_k=\frac{1}{K} = \frac{2}{\kappa}\left[\frac{ln\left(\frac{r_1+w_{ic}}{r_1}\right)}{(t\theta-t_{ic}\beta_{ic})}+ \frac{ln\left(\frac{r_1+L-w_{oc}}{r_1+w_{ic}}\right)}{t\theta}+\frac{ln\left(\frac{r_1+L}{r_1+L-w_{oc}}\right)}{(t\theta-t_{oc}\beta_{oc})}\right]\tag{4}$$
For electrical Resistance, we add $(3A+3B+3C)$,

$$R_e =\frac{2\rho}{(t\theta-t_{ic}\beta_{ic})}ln\left(\frac{r_1+w_{ic}}{r_1}\right)+\frac{2\rho}{t\theta}ln\left(\frac{r_1+L-w_{oc}}{r_1+w_{ic}}\right)+\frac{2\rho}{(t\theta-t_{oc}\beta_{oc})}ln\left(\frac{r_1+L}{r_1+L-w_{oc}}\right)$$
$$R_e =2\rho\left[\frac{ln\left(\frac{r_1+w_{ic}}{r_1}\right)}{(t\theta-t_{ic}\beta_{ic})}+\frac{ln\left(\frac{r_1+L-w_{oc}}{r_1+w_{ic}}\right)}{t\theta}+\frac{ln\left(\frac{r_1+L}{r_1+L-w_{oc}}\right)}{(t\theta-t_{oc}\beta_{oc})}\right]\tag{5}$$
Equations $(4)$ and $(5)$ are for a single TE leg, either P type or N type.

## 1.4 Correction for the Gaps

We need to reduce this from the overall area calculated above.

![[image-5.png]]


![[image-6.png]]

### 1.4.1 Radial Gaps - Ceramic Insulators
To isolate stages electrically but conduct heat. Need a good thermal conductor but electrical insulator. This gap needs to be as thin as manufacturable, and therefore this is a global parameter common for all the stages. Let's denote this length in radial direction as $w_{is}$ and name it as **Inter-Stage insulation width**. 

Mathematically this affects only the first part and the second part of the integration, and it only affects the limits of the integration.
$$r_1\longrightarrow r_1+w_{is}/2$$
$$r_2\longrightarrow r_2-w_{is}/2\longrightarrow (r_1+L)-w_{is}/2$$
>[!At first and last stages]
>When $i=1$ or $i=N$, the full width must be absorbed by the stage itself, and therefore 
>$$w_{is}/2\longrightarrow w_{is}$$

### 1.4.2 Azimuthal Gaps - Insulation 
To isolate elements (and P,N within one element) electrically, and to avoid back conduction. Need a both thermally and electrically insulated material for this purpose. Since its just eating away real estate and increase back conduction with area is increased, we need to make this as thin as possible. Therefore only a global parameter is needed. We have 2 options.

- So far everything in the azymuthal direction were parameterized using angles - wire dimensions and TEC dimensions. but if we parameterize the azimuthal gaps that way, and for it to constantly go out as a conical section, it would expand in size as it reaches the edge. so it will sever no real purpose but take the area that TEC elements could have occupied otherwise. not to mention the increased azimuthal dimension would make it more likely to perform back conduction degrading performance.

- Or we can parameterize it like, we take the boundary between TEC elements and offset it by some number to both sides, creating the new boundary. this will make sure the insulator remains constant and it won't eat away from the cross section area of TEC legs. 

==Second appraoch is more appropriate==. So we can simply model this using an **arc length parameter**, and lets name it _Azimuthal insulation width_ - $w_{az}$. Then the effective arc length would be,
$$S(r)=(r\theta/2-w_{az})$$
Multiplication with 2 has to come from 2 gaps - one that exist between P,N elements within the same TEC element and the one between TEC elements. However, since we consider only one leg, we have to divide it by 2. The area becomes,
$$A(r)=S(r)t=t\left(r\frac{\theta}{2}-w_{az}\right)$$
This reduces to previous expression $A(r) = t\,(\theta r/2)$ when $w_{\mathrm{az}} = 0$.

---
#### 1.4.2.1 Azimuthal Gaps - Insulation (Corrected Derivation)
We must now use this new area function to perform the integrations for electrical and thermal resistance. #### **Electrical Resistance Integration** The general formula for resistance is: $$R_e = \int_{r_1}^{r_2} \frac{\rho}{A(r)} dr$$ Substituting our new area function $A(r)$: $$R_e = \int_{r_1}^{r_2} \frac{\rho}{t(r\frac{\theta}{2} - w_{az})} dr$$ We can factor out the constants: $$R_e = \frac{\rho}{t} \int_{r_1}^{r_2} \frac{1}{r\frac{\theta}{2} - w_{az}} dr$$ This integral is of the standard form $\int \frac{1}{ax-b}dx = \frac{1}{a}\ln|ax-b|$. In our case, $a = \frac{\theta}{2}$ and $b = w_{az}$. Performing the integration gives: $$R_e = \frac{\rho}{t} \left[ \frac{1}{\theta/2} \ln\left|r\frac{\theta}{2} - w_{az}\right| \right]_{r_1}^{r_2}$$ $$R_e = \frac{2\rho}{t\theta} \left[ \ln\left|r_2\frac{\theta}{2} - w_{az}\right| - \ln\left|r_1\frac{\theta}{2} - w_{az}\right| \right]$$ Combining the logarithmic terms yields the final expression for the electrical resistance of a single leg with a constant-width azimuthal gap: $$R_e = \frac{2\rho}{t\theta} \ln\left( \frac{r_2\frac{\theta}{2} - w_{az}}{r_1\frac{\theta}{2} - w_{az}} \right)$$**Thermal Resistance Integration** The derivation for thermal resistance, $R_k = 1/K$, is identical. We replace electrical resistivity $\rho$ with the inverse of thermal conductivity, $1/\kappa$: $$R_k = \frac{1}{K} = \int_{r_1}^{r_2} \frac{1}{\kappa A(r)} dr$$ $$R_k = \frac{1}{\kappa t} \int_{r_1}^{r_2} \frac{1}{r\frac{\theta}{2} - w_{az}} dr$$ Following the exact same integration steps, we find the thermal resistance: $$R_k = \frac{2}{\kappa t\theta} \ln\left( \frac{r_2\frac{\theta}{2} - w_{az}}{r_1\frac{\theta}{2} - w_{az}} \right)$$

---
#### 1.4.2.2 Final Corrected Resistance Values for P and N Legs (Combining All Effects)

Now, we will write the final, fully corrected expressions for a single leg. We combine the three-region model for the wires (using your original, correct area formula) with the new integration results for the azimuthal gaps. The radial gaps ($w_{is}$) simply adjust the integration limits. First, let's define the integration limits and the area functions for the three regions. 

**Integration Limits:** 
* Start of leg: $r_{start} = r_1 + w_{is}/2$ 
* End of leg: $r_{end} = r_1 + L - w_{is}/2$ *(Note: As you stated, for the first and last stages, the appropriate boundary uses $w_{is}$ instead of $w_{is}/2$.)* 

**Area Functions:** 

The baseline area considering the azimuthal gap is $A_{base}(r) = t(r\frac{\theta}{2} - w_{az})$. In the wire regions, we subtract the area of the wire itself. 
1. **Inner Connect Region:** $r \in [r_{start}, r_1+w_{ic}]$ $$A_1(r) = t\left(r\frac{\theta}{2} - w_{az}\right) - t_{ic}r\frac{\beta_{ic}}{2} = r\left(t\frac{\theta}{2} - t_{ic}\frac{\beta_{ic}}{2}\right) - tw_{az}$$ 
2. **Full TE Region:** $r \in [r_1+w_{ic}, r_1+L-w_{oc}]$ $$A_2(r) = t\left(r\frac{\theta}{2} - w_{az}\right) = r\left(t\frac{\theta}{2}\right) - tw_{az}$$
3. **Outer Connect Region:** $r \in [r_1+L-w_{oc}, r_{end}]$ $$A_3(r) = t\left(r\frac{\theta}{2} - w_{az}\right) - t_{oc}r\frac{\beta_{oc}}{2} = r\left(t\frac{\theta}{2} - t_{oc}\frac{\beta_{oc}}{2}\right) - tw_{az}$$
Each resistance calculation involves summing three integrals of the form $\int \frac{dr}{Cr - D}$, which evaluates to $\frac{1}{C}\ln\left(\frac{Cr_{upper}-D}{Cr_{lower}-D}\right)$. 

#### 1.4.2.3 Total Electrical Resistance for a Single Leg ($R_e$)
The total electrical resistance is the sum of the resistances of the three regions in series, $R_e = R_{e,1} + R_{e,2} + R_{e,3}$. For a **P-type leg** (subscript `p`):
$$\begin{align} R_{e,p} &= \rho_p \left[ \frac{1}{t\frac{\theta}{2} - t_{ic}\frac{\beta_{ic}}{2}} \ln\left(\frac{(r_1+w_{ic})(t\frac{\theta}{2}-t_{ic}\frac{\beta_{ic}}{2}) - tw_{az}}{r_{start}(t\frac{\theta}{2}-t_{ic}\frac{\beta_{ic}}{2}) - tw_{az}}\right) \right. \\\\ &+ \frac{1}{t\frac{\theta}{2}} \ln\left(\frac{(r_1+L-w_{oc})(t\frac{\theta}{2}) - tw_{az}}{(r_1+w_{ic})(t\frac{\theta}{2}) - tw_{az}}\right) \\\\ & \left.+ \frac{1}{t\frac{\theta}{2} - t_{oc}\frac{\beta_{oc}}{2}} \ln\left(\frac{r_{end}(t\frac{\theta}{2}-t_{oc}\frac{\beta_{oc}}{2}) - tw_{az}}{(r_1+L-w_{oc})(t\frac{\theta}{2}-t_{oc}\frac{\beta_{oc}}{2}) - tw_{az}}\right) \right]\end{align}$$

For an **N-type leg** (subscript `n`), the expression is identical in form, using its specific electrical resistivity: 
$$R_{e,n} = \rho_n \left[ \dots (\text{same terms as for } R_{e,p}) \dots \right]$$
#### 1.4.2.4 Total Thermal Resistance for a Single Leg ($R_k$)** 
Similarly, the total thermal resistance is the sum of the series resistances, $R_k = R_{k,1} + R_{k,2} + R_{k,3}$. For a **P-type leg** (subscript `p`): 
$$\begin{align} R_{k,p} &= \frac{1}{\kappa_p} \left[ \frac{1}{t\frac{\theta}{2} - t_{ic}\frac{\beta_{ic}}{2}} \ln\left(\frac{(r_1+w_{ic})(t\frac{\theta}{2}-t_{ic}\frac{\beta_{ic}}{2}) - tw_{az}}{r_{start}(t\frac{\theta}{2}-t_{ic}\frac{\beta_{ic}}{2}) - tw_{az}}\right) \right. \\\\ &+ \frac{1}{t\frac{\theta}{2}} \ln\left(\frac{(r_1+L-w_{oc})(t\frac{\theta}{2}) - tw_{az}}{(r_1+w_{ic})(t\frac{\theta}{2}) - tw_{az}}\right) \\\\ & +\left. \frac{1}{t\frac{\theta}{2} - t_{oc}\frac{\beta_{oc}}{2}} \ln\left(\frac{r_{end}(t\frac{\theta}{2}-t_{oc}\frac{\beta_{oc}}{2}) - tw_{az}}{(r_1+L-w_{oc})(t\frac{\theta}{2}-t_{oc}\frac{\beta_{oc}}{2}) - tw_{az}}\right) \right]\end{align}$$For an **N-type leg** (subscript `n`), the expression uses its specific thermal conductivity: 
$$R_{k,n} = \frac{1}{\kappa_n} \left[ \dots (\text{same terms as for } R_{k,p}) \dots \right]$$
---
## 1.5 Final Corrected Expressions

The derivations in the previous section show that the resistance expressions for the thermoelectric legs are composed of a material property ($\rho$ or $1/\kappa$) and a complex geometric term. 

### 1.5.1 Geometric Factor
To simplify the notation, we can define this common term as a single **Geometric Factor, $G$**. The geometric factor $G$ encapsulates all the geometric complexities of a single TE leg, including the trapezoidal shape, the embedded wires, and the azimuthal insulation gaps. It is defined as: $$ \begin{align} G = & \frac{1}{t\frac{\theta}{2} - t_{ic}\frac{\beta_{ic}}{2}} \ln\left(\frac{(r_1+w_{ic})(t\frac{\theta}{2}-t_{ic}\frac{\beta_{ic}}{2}) - tw_{az}}{r_{start}(t\frac{\theta}{2}-t_{ic}\frac{\beta_{ic}}{2}) - tw_{az}}\right) \\ & + \frac{1}{t\frac{\theta}{2}} \ln\left(\frac{(r_1+L-w_{oc})(t\frac{\theta}{2}) - tw_{az}}{(r_1+w_{ic})(t\frac{\theta}{2}) - tw_{az}}\right) \\ & + \frac{1}{t\frac{\theta}{2} - t_{oc}\frac{\beta_{oc}}{2}} \ln\left(\frac{r_{end}(t\frac{\theta}{2}-t_{oc}\frac{\beta_{oc}}{2}) - tw_{az}}{(r_1+L-w_{oc})(t\frac{\theta}{2}-t_{oc}\frac{\beta_{oc}}{2}) - tw_{az}}\right) \end{align} $$
where $r_{start} = r_1 + w_{is}/2$ and $r_{end} = r_1 + L - w_{is}/2$. 

### 1.5.2 Simplification
Using this geometric factor, the electrical and thermal resistances for the P-type and N-type legs can be written concisely: **Electrical Resistance of TE Legs:** 
$$R_{e,p} = \rho_p G$$$$R_{e,n} = \rho_n G$$ **Thermal Resistance of TE Legs:** 
$$R_{k,p} = \frac{G}{\kappa_p}$$
$$R_{k,n} = \frac{G}{\kappa_n}$$
### 1.5.3 Overall Thermal Conductance of a TEC Element 

For thermal analysis, the P-type and N-type legs are in parallel. Therefore, their thermal conductances add up. The total thermal conductance, $K_{total}$, of a single P-N element is the sum of the individual conductances. $$K_{total,TE} = K_p + K_n$$ We know that conductance $K$ is the reciprocal of thermal resistance $R_k$. $$K_{total,TE} = \frac{1}{R_{k,p}} + \frac{1}{R_{k,n}}$$ Substituting the expressions using the geometric factor $G$: $$K_{total.TE} = \frac{\kappa_p}{G} + \frac{\kappa_n}{G}$$ This simplifies to a very clean final expression: $$K_{total,TE} = \frac{\kappa_p + \kappa_n}{G}$$
# 2 Electrical Wiring 

We model the interconnect as a conductive sheet. The total resistance can be found by considering it as a set of infinitesimally thin, parallel resistive strips, each at a different radius `r` and having a width `dr`. The resistance of one such infinitesimal strip to azimuthal current flow is: $$d R_{strip} = \rho \frac{\text{length}}{\text{area}} = \rho_c \frac{r \cdot \beta}{(dr) \cdot t_c}$$ Here, the length of the current path is the arc length `rβ`, and the cross-sectional area is `dr * t_c`. These strips are all in **parallel**. Therefore, we must sum their conductances (`1/dR`). The total conductance `1/R` is the integral of the individual strip conductances: $$ \frac{1}{R} = \int \frac{1}{dR_{strip}} = \int \frac{t_c \cdot dr}{\rho_c \cdot r \beta} $$ Factoring out the constants gives: $$ \frac{1}{R} = \frac{t_c}{\rho_c \beta} \int \frac{dr}{r} = \frac{t_c}{\rho_c \beta} \left[ \ln(r) \right] $$ Now we apply this general formula to the inner and outer interconnects.
## 2.1 Inner Interconnect Resistance ($R_{ic}$)
For the inner interconnect, we apply the integration limits and the correct angle: 
* Angle: $\beta = \beta_{ic}$
* Thickness: $t_c = t_{ic}$ 
* Radial Limits: from $r_1$ to $r_1 + w_{ic}$ $$ \begin{align} \frac{1}{R_{ic}} &= \frac{t_{ic}}{\rho_c \beta_{ic}} \int_{r_1}^{r_1+w_{ic}} \frac{dr}{r} \\ \frac{1}{R_{ic}} &= \frac{t_{ic}}{\rho_c \beta_{ic}} \left[ \ln(r) \right]_{r_1}^{r_1+w_{ic}} \\ \frac{1}{R_{ic}} &= \frac{t_{ic}}{\rho_c \beta_{ic}} \ln\left( \frac{r_1+w_{ic}}{r_1} \right) \end{align} $$
Inverting the expression to find the resistance $R_{ic}$: 
$$ R_{ic} = \frac{\rho_c \beta_{ic}}{t_{ic} \ln\left( \frac{r_1+w_{ic}}{r_1} \right)} $$
## 2.2 Outer Interconnect Resistance ($R_{oc}$)

For a single outer connect, the logic is identical, but we use the corresponding parameters: * Angle: $\beta = \beta_{oc}/2$ (for one piece) * Thickness: $t_c = t_{oc}$ * Radial Limits: from $r_1+L-w_{oc}$ to $r_1+L$ $$ \begin{align} \frac{1}{R_{oc}} &= \frac{t_{oc}}{\rho_c (\beta_{oc}/2)} \int_{r_1+L-w_{oc}}^{r_1+L} \frac{dr}{r} \\ \frac{1}{R_{oc}} &= \frac{2t_{oc}}{\rho_c \beta_{oc}} \left[ \ln(r) \right]_{r_1+L-w_{oc}}^{r_1+L} \\ \frac{1}{R_{oc}} &= \frac{2t_{oc}}{\rho_c \beta_{oc}} \ln\left( \frac{r_1+L}{r_1+L-w_{oc}} \right) \end{align} $$ Inverting the expression for the resistance $R_{oc}$ of a single outer piece: 
$$ R_{oc} = \frac{\rho_c \beta_{oc}}{2t_{oc} \ln\left( \frac{r_1+L}{r_1+L-w_{oc}} \right)} $$
---
The total electrical resistance of a complete TEC unit would be,
$$R_{total,e} = R_{e,p} + R_{e,n} + R_{ic} + 2R_{oc}$$

